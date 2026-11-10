//! Two-phase search.
//!
//! Phase 1: IDA* over all internal moves (16, or 24 with slices) until the position is in
//! `H_A` = <U, F, BR, BL>, with the last move outside `H_A`. Phase 2: the exact table gives the
//! optimal A-move completion. Every solution splits uniquely into (longest prefix ending in a
//! non-A move) + (A-moves), so exhausting phase-1 depths up to `best - 1` proves optimality.
//! (Slices are not A-moves, so this holds with slices too.)
//!
//! Two variants ("sides") run in one traversal: targeting class A, and targeting class B (by
//! conjugating the position with a 90° rotation that swaps the classes).

use crate::cube::*;
use crate::fast::*;
use crate::geometry::{transpose, Mat, IDENTITY};
use crate::notation::*;
use crate::prune::decode;
use crate::tables::*;
use rayon::prelude::*;
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::sync::{Mutex, PoisonError};
use std::time::{Duration, Instant};

#[derive(Clone, Debug, Default)]
pub struct Limits {
    /// stop after this many phase-2 completions (solutions of any length)
    pub max_solutions: Option<u64>,
    /// stop after this much time
    pub time: Option<Duration>,
    /// stop once a solution of at most this length is found
    pub max_length: Option<usize>,
    /// stop after the phase-1 iteration with this depth bound (used by tests)
    pub max_depth: Option<usize>,
}

impl Limits {
    pub fn is_optimal(&self) -> bool {
        self.max_solutions.is_none()
            && self.time.is_none()
            && self.max_length.is_none()
            && self.max_depth.is_none()
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum Report {
    /// print only the best solution when the search ends
    Best,
    /// print every improvement as soon as it is found
    Improvements,
    /// print nothing (library/benchmark use)
    Quiet,
}

#[derive(Clone, Debug)]
pub struct Event {
    pub time: f64,
    pub length: usize,
}

#[derive(Clone, Debug)]
pub struct Outcome {
    pub solution: Option<Vec<PhysMove>>,
    pub length: Option<usize>,
    pub solutions_found: u64,
    pub nodes: u64,
    /// phase-1 leaves reached: paths that pass every bound at their iteration's depth (on
    /// either side: a path that is a leaf on both sides counts twice)
    pub leaves: u64,
    /// the leaves of each side: `[0]` paths to `H_A`, `[1]` paths to `H_B` (searched on the
    /// position conjugated by `SWAP_AB`)
    pub leaves_by_side: [u64; 2],
    pub elapsed: f64,
    /// true if the search proved that no shorter solution exists
    pub optimal: bool,
    /// improvement history
    pub events: Vec<Event>,
}

/// One side of the search: the position in the frame where the side's target class is class A
/// (side 0: the position as given, target `H_A`; side 1: conjugated by `SWAP_AB`, target `H_B`).
struct Side {
    cube: Cube,
    sym: Mat,
}

/// Sides of a node, as a bit set: bit `s` = the node is in the search tree of side `s` (the
/// tree the separate search for that side would visit). See [`Shared::dfs`].
/// In the evaluation of a node (`expand`, `stage2`, `finish`) the bits are slots instead (see
/// [`Scratch`]): a node is evaluated as `BOTH`, or as `SIDE_A`, alive on one side only, which
/// is then held in slot 0.
const SIDE_A: u8 = 1;
const SIDE_B: u8 = 2;
const BOTH: u8 = 3;

/// Whether side `s` is in the side set `sides`.
#[inline(always)]
const fn on(sides: u8, s: usize) -> bool {
    sides >> s & 1 != 0
}

/// The state shared by the threads of one solve, for a move set of `N` internal moves.
struct Shared<'a, I: AeIndex, const N: usize> {
    t: &'a Tables<I>,
    sides: [Side; 2],
    /// the root, in both frames
    root: Task,
    frame: Frame,
    limits: Limits,
    report: Report,
    start: Instant,
    best_len: AtomicUsize,
    /// length the class-B bound prunes against: the best length, or `target + 1` in a
    /// target-length pass (see `solve`)
    prune_len: AtomicUsize,
    count: AtomicU64,
    /// a solution of at most this length is optimal (depth + 1 of the iteration in progress,
    /// 0 before the first): see `solve`
    proof_len: AtomicUsize,
    /// set when such a solution was found
    proven: AtomicBool,
    nodes: AtomicU64,
    /// nodes evaluated once for both sides (see [`Counter::shared`])
    nodes_both: AtomicU64,
    /// phase-1 leaves reached on each side (most are solutions; the tables certify the rest
    /// only approximately, so each is checked in full)
    leaves: [AtomicU64; 2],
    stop: AtomicBool,
    best: Mutex<(Option<Vec<PhysMove>>, Vec<Event>)>,
    /// the side-1 move of each move (`SWAP_AB` conjugation), indexed `& 31`
    mb: [usize; 32],
    /// the move of each side for each move: `[0]` itself, `[1]` `mb`
    mv: [[usize; 32]; 2],
    /// largest CB distance (8): with at least this many moves left CB cannot prune
    cb_max: u8,
}

/// A node of the top levels of the search tree (see [`Shared::tasks`]).
#[derive(Clone)]
struct Task {
    path: Vec<usize>,
    /// the position in the frame of each side, and its exact distances there
    x: [Fc; 2],
    h_ae: [u8; 2],
    h_cb: [u8; 2],
}

/// Per-worker counters.
#[derive(Default)]
struct Counter {
    /// evaluated nodes (children generated)
    nodes: u64,
    /// of those, the ones the two separate searches would both have evaluated (children of a
    /// node alive on both sides): the evaluations the merged traversal saves
    shared: u64,
}

/// Internal move on the SWAP_AB-conjugated position matching each move (indexed `& 31`).
fn swap_moves<const N: usize>() -> [usize; 32] {
    std::array::from_fn(|m| {
        if m < N {
            conjugate_move(m, &SWAP_AB)
        } else {
            0
        }
    })
}

const CHECK_EVERY: u64 = 1 << 14;

/// The AE table is checked first (before CB) with at least this many moves left after the
/// children: the CB table (max 8) prunes almost nothing there (+7% / +10% nodes/s at 48T / 1T
/// when added with 7+ moves, opt/ae-first; lowered to 6 in the opt/stages stack, measured
/// only as part of it: +10% / +14% nodes/s).
const AE_FIRST: usize = 6;

/// No deferred CB lookup (see `finish`).
const NO_CB: usize = usize::MAX;

/// Scratch arrays of one node's evaluation: its children and their indices and bounds in two
/// slots, one per side (`N` moves; indices are reduced `% N`, which is `& 15` for 16 moves,
/// so that the compiler sees them in bounds). On a node alive on both sides slot `s` holds
/// side `s`; on a node alive on one side, slot 0 holds that side (side `r` in the functions
/// below) and slot 1 the other. A slot's arrays are only filled where needed: for the node's
/// own sides, and with the class-B bound for the other side too.
/// They live in a per-worker stack of frames, three per depth (see [`Shared::dfs`]), rather
/// than on the call stack, so that they are not re-initialised at every node (~3 KB of stores
/// per node; +1.3% nodes/s at 1T in optimal mode, when added).
#[derive(Clone)]
struct Scratch<const N: usize> {
    ys: [[Fc; N]; 2],
    ms: [u8; N],
    /// the sides each child is still alive on
    live: [u8; N],
    has: [[u8; N]; 2],
    hcs: [[u8; N]; 2],
    idx: [[usize; N]; 2],
    cbi: [[usize; N]; 2],
    pend: [[usize; N]; 2],
    ci: [[[usize; 4]; N]; 2],
    /// the children still alive after each stage (indices into the arrays above). Stages
    /// compact this list without branching on the prune decision, instead of copying
    /// positions (neutral when added: -1.5% nodes/s at 1T, within noise).
    sel: [u8; N],
    /// number of children
    n: usize,
    /// number alive after `stage2` (in `sel`)
    k: usize,
    /// the node's own exact distances on both sides and the smallest AE distance of a child
    /// alive on each side after `stage2`, and whether `stage2` issued the joint-table lookups
    h_ae: [u8; 2],
    h_cb: [u8; 2],
    min_ha: [u8; 2],
    cbe_done: bool,
}

impl<const N: usize> Scratch<N> {
    fn new() -> Self {
        Scratch {
            ys: [[Fc::ZERO; N]; 2],
            ms: [0; N],
            live: [0; N],
            has: [[0; N]; 2],
            hcs: [[0; N]; 2],
            idx: [[0; N]; 2],
            cbi: [[0; N]; 2],
            pend: [[NO_CB; N]; 2],
            ci: [[[0; 4]; N]; 2],
            sel: [0; N],
            n: 0,
            k: 0,
            h_ae: [0; 2],
            h_cb: [0; 2],
            min_ha: [0; 2],
            cbe_done: false,
        }
    }
}

/// Moves allowed after each move (index `MAX_MOVES`: no previous move) by [`allowed_after`],
/// as bit masks over the moves of the move set: `[0]` the face metric, `[1]` with slices.
static ALLOWED: [[u32; MAX_MOVES + 1]; 2] = [allowed_masks(false), allowed_masks(true)];

const fn allowed_masks(slice: bool) -> [u32; MAX_MOVES + 1] {
    let n = if slice { MAX_MOVES } else { N_MOVES };
    let mut a = [0u32; MAX_MOVES + 1];
    let mut last = 0;
    while last <= MAX_MOVES {
        let mut m = 0;
        while m < n {
            if allowed_after(slice, last, m) {
                a[last] |= 1 << m;
            }
            m += 1;
        }
        last += 1;
    }
    a
}

/// The moves whose slot is one of `slots`, as a bit mask.
const fn slot_mask(slots: [usize; 4]) -> u32 {
    let mut a = 0u32;
    let mut m = 0;
    while m < N_MOVES {
        let s = m >> 1;
        if s == slots[0] || s == slots[1] || s == slots[2] || s == slots[3] {
            a |= 1 << m;
        }
        m += 1;
    }
    a
}

/// The class-A moves as a bit mask.
const A_MOVES: u32 = slot_mask(A_SLOTS);

/// The class-B moves (the other face-metric moves; `SWAP_AB` maps the class-A moves onto
/// them, checked by `move_masks_match_rules`): the moves of side 1's frame that are class-A
/// moves there.
const B_MOVES: u32 = slot_mask([0, 3, 5, 6]);

/// The moves that can end a phase-1 path on a side set (index: the set): on side 0 a move
/// outside class A, on side 1 a move outside class B (slices are in neither class).
const FINAL_MOVES: [u32; 4] = [0, !A_MOVES, !B_MOVES, !A_MOVES | !B_MOVES];

/// The sides on which each move can end a phase-1 path (indexed `& 31`).
static FINAL_SIDES: [u8; 32] = {
    let mut a = [0u8; 32];
    let mut m = 0;
    while m < 32 {
        a[m] = (FINAL_MOVES[1] >> m & 1) as u8 | ((FINAL_MOVES[2] >> m & 1) as u8) << 1;
        m += 1;
    }
    a
};

/// Every move on both sides (the children of a node that are not leaves).
static ALL_SIDES: [u8; 32] = [BOTH; 32];

/// Weight 8^-h of a surviving node with bound h in a task's ordering statistic.
const STAT_W: [f32; 32] = {
    let mut w = [0f32; 32];
    let mut h = 0;
    let mut x = 1.0f32;
    while h < 32 {
        w[h] = x;
        x /= 8.0;
        h += 1;
    }
    w
};

impl<I: AeIndex, const N: usize> Shared<'_, I, N> {
    fn time_up(&self) -> bool {
        if let Some(t) = self.limits.time {
            if self.start.elapsed() >= t {
                self.stop.store(true, Ordering::Relaxed);
                return true;
            }
        }
        false
    }

    /// Move `m` (of side 0's frame) in the frame of side `s`.
    #[inline(always)]
    fn side_move(&self, s: usize, m: u8) -> usize {
        self.mv[s & 1][usize::from(m) & 31]
    }

    /// A phase-1 leaf on side `s`: `path` (moves of side 0's frame) ends in the side's target
    /// group by the tables.
    fn leaf(&self, s: usize, path: &[usize]) {
        let side = &self.sides[s & 1];
        self.leaves[s & 1].fetch_add(1, Ordering::Relaxed);
        // the path in the side's frame
        let mv = |m: usize| if s == 0 { m } else { self.mb[m & 31] };
        let cube = path.iter().fold(side.cube, |c, &m| c.apply(mv(m)));
        let d1 = path.len();
        // The tables only bound the distance to H_A, so membership is checked in full here
        // (with the small tier, for example, the edge chirality is not in the tables).
        if let Some(d2) = self.t.p2.length(&cube) {
            let n = self.count.fetch_add(1, Ordering::Relaxed) + 1;
            if self.limits.max_solutions.is_some_and(|k| n >= k) {
                self.stop.store(true, Ordering::Relaxed);
            }
            let total = d1 + d2 as usize;
            if total < self.best_len.load(Ordering::Relaxed) {
                let mut g = self.best.lock().unwrap_or_else(PoisonError::into_inner);
                if total < self.best_len.load(Ordering::Relaxed) {
                    // (the position is in H_A, so phase 2 always succeeds)
                    let p2 = self.t.p2.solve(&cube);
                    debug_assert!(p2.is_some(), "phase 2 failed");
                    let Some(p2) = p2 else { return };
                    // map the side's internal moves back to the original frame
                    let inv = transpose(&side.sym);
                    let internal: Vec<usize> = path
                        .iter()
                        .map(|&m| mv(m))
                        .chain(p2.iter().copied())
                        .map(|m| conjugate_move(m, &inv))
                        .collect();
                    let phys = solution_to_physical(self.frame, &internal, self.t.moves.wide);
                    let t = self.start.elapsed().as_secs_f64();
                    if self.report == Report::Improvements {
                        println!("{:>9.3}s  {:>2}  {}", t, total, format_moves(&phys));
                    }
                    g.0 = Some(phys);
                    g.1.push(Event {
                        time: t,
                        length: total,
                    });
                    self.best_len.store(total, Ordering::Relaxed);
                    self.prune_len.fetch_min(total, Ordering::Relaxed);
                    // Every iteration below the current one is complete, so a solution of
                    // length depth + 1 is optimal at once (see `solve`): stop the iteration.
                    if self.limits.max_solutions.is_none()
                        && total <= self.proof_len.load(Ordering::Relaxed)
                    {
                        self.proven.store(true, Ordering::Relaxed);
                        self.stop.store(true, Ordering::Relaxed);
                    }
                    if self.limits.max_length.is_some_and(|l| total <= l) {
                        self.stop.store(true, Ordering::Relaxed);
                    }
                }
            }
        }
    }

    /// Depth-first search below a task node (at depth `g = path.len()`, alive on `sides`)
    /// with bound `bound`.
    ///
    /// The search runs for two targets at once, `H_A` (side 0) and `H_B` (side 1, through the
    /// position conjugated by `SWAP_AB`; `SWAP_AB` maps canonical sequences to canonical
    /// sequences, so both sides walk the same move sequences, in side 0's moves). Each node
    /// carries its position in both frames and the set of sides it is alive on: a child is
    /// alive on side `s` if the node is and the child passes side `s`'s bounds (its distance
    /// to the side's target within the moves left, and, with `B`, its distance to the other
    /// side's target against the best solution so far; see `finish`). A node alive on both
    /// sides is evaluated once, where the two separate searches each evaluated it (`Counter`
    /// counts these). The tree is the union of the two separate searches' trees, with each
    /// node's sides being exactly the searches that visit it, so the optimality argument of
    /// `solve` (both searches complete every iteration) still holds.
    /// Speedup of the merged traversal: see `solve_n`.
    ///
    /// A node is evaluated in three steps: [`Self::expand`] generates its children and issues
    /// their first table lookups, [`Self::stage2`] reads them and issues the second round,
    /// and [`Self::finish`] reads that and searches below the surviving children. Siblings
    /// are pipelined: the next sibling's first two steps run before a child's subtree is
    /// searched, so that their table reads find the data in cache. The search waited on
    /// memory about half of the time; when added: +40% nodes/s at 48T and +45% at 1T
    /// (time-limited benchmark), optimal-mode time 0.70x / 0.72x (pipelining only the first
    /// step: +22% / +18% nodes/s in a quick optimal-mode test).
    fn dfs<const B: bool>(
        &self,
        task: &Task,
        sides: u8,
        path: &mut Vec<usize>,
        bound: usize,
        ctr: &mut Counter,
        stat: &mut f32,
        frames: &mut [[Scratch<N>; 3]],
    ) {
        let g = path.len();
        if g == bound {
            for s in 0..2 {
                let fin = path.last().is_none_or(|&m| on(FINAL_SIDES[m & 31], s));
                if on(sides, s) && task.h_ae[s] == 0 && task.h_cb[s] == 0 && fin {
                    self.leaf(s, path);
                }
            }
            return;
        }
        if self.stop.load(Ordering::Relaxed) {
            return;
        }
        // (the caller provides three frames per remaining level)
        let Some(([fr, _, _], frames)) = frames.split_first_mut() else {
            debug_assert!(false, "out of search frames");
            return;
        };
        let last = path.last().copied().unwrap_or(usize::MAX);
        let budget = bound - g - 1;
        // (a node alive on side 1 only has that side in slot 0)
        let (a, b) = if sides == SIDE_B { (1, 0) } else { (0, 1) };
        let r = usize::from(sides == SIDE_B);
        let kind = if sides == BOTH { BOTH } else { SIDE_A };
        let x = [&task.x[a], &task.x[b]];
        let (ha, hc) = ([task.h_ae[a], task.h_ae[b]], [task.h_cb[a], task.h_cb[b]]);
        self.expand_any(kind, r, x, last, budget, fr, ctr);
        let lim = if B {
            self.prune_len.load(Ordering::Relaxed).saturating_sub(g + 2)
        } else {
            0
        };
        self.stage2_any::<B>(kind, r, x, ha, hc, [NO_CB; 2], budget, lim, fr);
        self.finish_any::<B>(kind, r, x, path, bound, ctr, stat, fr, frames);
    }

    /// [`Self::expand`] for a node of kind `sides` (`BOTH`, or `SIDE_A`: alive on side `r`
    /// only, held in slot 0). Not inlined: one copy per kind.
    #[inline(never)]
    fn expand_any(
        &self,
        sides: u8,
        r: usize,
        x: [&Fc; 2],
        last: usize,
        budget: usize,
        fr: &mut Scratch<N>,
        ctr: &mut Counter,
    ) {
        if sides == BOTH {
            self.expand::<BOTH>(0, x, last, budget, fr, ctr);
        } else {
            self.expand::<SIDE_A>(r, x, last, budget, fr, ctr);
        }
    }

    /// [`Self::stage2`] for a node of kind `sides` (see `expand_any`).
    #[inline(never)]
    fn stage2_any<const B: bool>(
        &self,
        sides: u8,
        r: usize,
        x: [&Fc; 2],
        h_ae: [u8; 2],
        h_cb: [u8; 2],
        cb_pend: [usize; 2],
        budget: usize,
        lim: usize,
        fr: &mut Scratch<N>,
    ) {
        if sides == BOTH {
            self.stage2::<BOTH, B>(0, x, h_ae, h_cb, cb_pend, budget, lim, fr);
        } else {
            self.stage2::<SIDE_A, B>(r, x, h_ae, h_cb, cb_pend, budget, lim, fr);
        }
    }

    /// [`Self::finish`] for a node of kind `sides` (see `expand_any`).
    #[inline(always)]
    fn finish_any<const B: bool>(
        &self,
        sides: u8,
        r: usize,
        x: [&Fc; 2],
        path: &mut Vec<usize>,
        bound: usize,
        ctr: &mut Counter,
        stat: &mut f32,
        fr: &mut Scratch<N>,
        frames: &mut [[Scratch<N>; 3]],
    ) {
        if sides == BOTH {
            self.finish::<BOTH, B>(0, x, path, bound, ctr, stat, fr, frames);
        } else {
            self.finish::<SIDE_A, B>(r, x, path, bound, ctr, stat, fr, frames);
        }
    }

    /// First step of the evaluation of a node of kind `S` (positions `x` in its slots), with
    /// `budget` moves left after its children (`last`: the move leading to it): generate the
    /// children into `fr` on the node's sides, count them, and issue their first table
    /// lookups (AE with 6+ moves left, CB otherwise). Every stage issues the lookups of all
    /// children before reading any of them, so the cache misses of different children
    /// overlap (2.0x nodes/s at 48T and 1T when added, opt/prefetch).
    #[inline(always)]
    fn expand<const S: u8>(
        &self,
        r: usize,
        x: [&Fc; 2],
        last: usize,
        budget: usize,
        fr: &mut Scratch<N>,
        ctr: &mut Counter,
    ) {
        let t = self.t;
        let f = fast();
        let mut moves = ALLOWED[usize::from(N > N_MOVES)][last.min(MAX_MOVES)];
        // the children are leaves: only moves that end a phase-1 path on one of the sides
        let fin = if budget == 0 {
            moves &= FINAL_MOVES[(S as usize) << r];
            &FINAL_SIDES
        } else {
            &ALL_SIDES
        };
        if S == BOTH {
            // (of the leaves, only those that end a path on both sides are shared)
            let both = if budget == 0 {
                moves & FINAL_MOVES[1] & FINAL_MOVES[2]
            } else {
                moves
            };
            ctr.shared += u64::from(both.count_ones());
        }
        let fw0 = if on(S, 0) {
            Self::flip_words(f, x[0])
        } else {
            (0, 0)
        };
        let fw1 = if on(S, 1) {
            Self::flip_words(f, x[1])
        } else {
            (0, 0)
        };
        let mut n = 0;
        while moves != 0 {
            let m = moves.trailing_zeros() as usize;
            moves &= moves - 1;
            if on(S, 0) {
                fr.ys[0][n % N] = f.apply_any::<N>(x[0], fw0.0, fw0.1, self.side_move(r, m as u8));
            }
            if on(S, 1) {
                fr.ys[1][n % N] = f.apply_any::<N>(x[1], fw1.0, fw1.1, self.mb[m & 31]);
            }
            fr.ms[n % N] = m as u8;
            // (a node alive on one side has it in slot 0, and the moves are its leaf moves)
            fr.live[n % N] = if S == BOTH { fin[m & 31] } else { SIDE_A };
            n += 1;
        }
        fr.n = n;
        // count every evaluated node; check the clock regularly (`finish` sees the stop)
        let before = ctr.nodes;
        ctr.nodes += n as u64;
        if before / CHECK_EVERY != ctr.nodes / CHECK_EVERY {
            self.time_up();
        }
        for s in 0..2 {
            if !on(S, s) {
                continue;
            }
            if budget >= AE_FIRST {
                t.ae_ix.index_batch_fast_pf::<N>(
                    f,
                    &fr.ys[s][..n],
                    &mut fr.idx[s][..n],
                    &t.ae.data,
                );
            } else {
                for i in 0..n {
                    fr.cbi[s][i] = CbTable::index_fast(f, &fr.ys[s][i]);
                    prefetch_at(&t.cb.tab.data, fr.cbi[s][i] >> 5);
                }
            }
        }
    }

    /// Second step (see [`Self::dfs`]): read the first lookups of the children in `fr`, and
    /// issue the lookups of the next stages for those still alive on a side. `h_ae`, `h_cb`
    /// are the node's distances on both sides; `cb_pend[s]`, if not `NO_CB`, is the index of
    /// the node's own CB entry on side `s` (prefetched by its parent), and `h_cb[s]` its
    /// parent's distance. `lim` is the children's class-B limit (see `finish`) as of now: it
    /// can only fall before `finish` applies it, so on a node alive on both sides, where both
    /// sides' first bounds are known here, a side whose class-B bound already fails is dropped
    /// at once, before its CB and joint-table lookups. In the last iteration of a full proof
    /// (always with slices) the class-B limit equals the moves left, so there almost every
    /// node is alive on both sides and every child that fails on one side fails on the other
    /// (with slices optimal 1T 0.80x, 48T 0.91x time when added; face turns neutral).
    #[inline(always)]
    fn stage2<const S: u8, const B: bool>(
        &self,
        r: usize,
        x: [&Fc; 2],
        h_ae: [u8; 2],
        h_cb: [u8; 2],
        cb_pend: [usize; 2],
        budget: usize,
        lim: usize,
        fr: &mut Scratch<N>,
    ) {
        let t = self.t;
        let f = fast();
        // the sides whose bounds the children need: their own, and with the class-B bound the
        // other one too
        let need = if B { BOTH } else { S };
        let Scratch {
            ys,
            ms,
            live,
            has,
            hcs,
            idx,
            cbi,
            ci,
            sel,
            n,
            k: kk,
            h_ae: fh_ae,
            h_cb: fh_cb,
            min_ha: fmin_ha,
            cbe_done,
            ..
        } = fr;
        let n = *n;
        let mut k = 0;
        let mut h_cb = h_cb;
        *cbe_done = false;
        if budget >= AE_FIRST {
            let mut min_ha = [u8::MAX; 2];
            for i in 0..n {
                let mut c = live[i];
                for s in 0..2 {
                    if !on(S, s) {
                        continue;
                    }
                    let ha = decode(h_ae[s], t.ae.get(idx[s][i]));
                    has[s][i] = ha;
                    c &= !(u8::from(ha as usize > budget) << s);
                }
                if B && S == BOTH {
                    // (early class-B check, see above)
                    c &= !(u8::from(has[1][i] as usize > lim)
                        | u8::from(has[0][i] as usize > lim) << 1);
                }
                for s in 0..2 {
                    if on(S, s) {
                        min_ha[s] = min_ha[s].min(if on(c, s) { has[s][i] } else { u8::MAX });
                    }
                }
                live[i] = c;
                sel[k % N] = i as u8;
                k += usize::from(c != 0);
            }
            *fmin_ha = min_ha;
            for s in 0..2 {
                if !on(S, s) {
                    continue;
                }
                for &i in &sel[..k] {
                    let i = usize::from(i) % N;
                    cbi[s][i] = CbTable::index_fast(f, &ys[s][i]);
                    prefetch_at(&t.cb.tab.data, cbi[s][i] >> 5);
                }
            }
            if k > 0 {
                // this node's own CB entries were prefetched by the parent, before its AE
                // stage
                for s in 0..2 {
                    if on(need, s) && cb_pend[s] != NO_CB {
                        h_cb[s] = decode(h_cb[s], t.cb.tab.get(cb_pend[s]));
                    }
                }
            }
            // The class-B lookups (the side the node is not alive on) are issued here too,
            // so that their latency overlaps the CB and joint-table stages (+6% / +4%
            // nodes/s at 48T / 1T, optimal 0.95x time, when added, opt/bpath).
            if B && S != BOTH && k > 0 {
                let o = usize::from(S == SIDE_A);
                let (fw, fws) = Self::flip_words(f, x[o]);
                for &i in &sel[..k] {
                    let i = usize::from(i) % N;
                    ys[o][i] = f.apply_any::<N>(x[o], fw, fws, self.side_move(o ^ r, ms[i]));
                }
                t.ae_ix
                    .index_sel_pf(f, &ys[o], &sel[..k], &mut idx[o], &t.ae.data);
                for &i in &sel[..k] {
                    let i = usize::from(i) % N;
                    cbi[o][i] = CbTable::index_fast(f, &ys[o][i]);
                    prefetch_at(&t.cb.tab.data, cbi[o][i] >> 5);
                }
            }
            // The joint table goes before reading the CB words, so their latency overlaps it
            // (part of the opt/stages stack: +10% / +14% nodes/s at 48T / 1T as a whole).
            if budget <= 8 {
                if let Some(cbe) = &t.cbe {
                    for s in 0..2 {
                        if !on(S, s) {
                            continue;
                        }
                        for &i in &sel[..k] {
                            let i = usize::from(i) % N;
                            if S != BOTH || on(live[i], s) {
                                cbe.prefetch_all(f, &ys[s][i], &mut ci[s][i]);
                            }
                        }
                    }
                    *cbe_done = true;
                }
            }
        } else {
            for i in 0..n {
                let mut c = live[i];
                for s in 0..2 {
                    if !on(S, s) {
                        continue;
                    }
                    let hc = decode(h_cb[s], t.cb.tab.get(cbi[s][i]));
                    hcs[s][i] = hc;
                    c &= !(u8::from(hc as usize > budget) << s);
                }
                if B && S == BOTH {
                    // (early class-B check, see above)
                    c &= !(u8::from(hcs[1][i] as usize > lim)
                        | u8::from(hcs[0][i] as usize > lim) << 1);
                }
                live[i] = c;
                sel[k % N] = i as u8;
                k += usize::from(c != 0);
            }
            for s in 0..2 {
                if on(S, s) {
                    t.ae_ix
                        .index_sel_pf(f, &ys[s], &sel[..k], &mut idx[s], &t.ae.data);
                }
            }
        }
        *kk = k;
        *fh_ae = h_ae;
        *fh_cb = h_cb;
    }

    /// Last step of the evaluation of a node at depth `g = path.len()` of kind `S` (positions
    /// `x` in its slots; see [`Self::dfs`]), after [`Self::expand`] and [`Self::stage2`] on
    /// `fr`, and the search below it.
    fn finish<const S: u8, const B: bool>(
        &self,
        r: usize,
        x: [&Fc; 2],
        path: &mut Vec<usize>,
        bound: usize,
        ctr: &mut Counter,
        stat: &mut f32,
        fr: &mut Scratch<N>,
        frames: &mut [[Scratch<N>; 3]],
    ) {
        if self.stop.load(Ordering::Relaxed) {
            return;
        }
        let g = path.len();
        let budget = bound - g - 1;
        let t = self.t;
        let f = fast();
        let need = if B { BOTH } else { S };
        let Scratch {
            ys,
            ms,
            live,
            has,
            hcs,
            idx,
            cbi,
            pend,
            ci,
            sel,
            n: _,
            k,
            h_ae,
            h_cb,
            min_ha,
            cbe_done,
        } = fr;
        let (mut k, h_ae, h_cb, min_ha, cbe_done) = (*k, *h_ae, *h_cb, *min_ha, *cbe_done);
        let mut n = 0;
        let mut b_done = false;
        if budget >= AE_FIRST {
            if cbe_done {
                if let Some(cbe) = &t.cbe {
                    let mut k2 = 0;
                    for j in 0..k {
                        let i = sel[j];
                        let u = usize::from(i) % N;
                        let mut c = live[u];
                        for s in 0..2 {
                            if on(S, s) && (S != BOTH || on(c, s)) {
                                let ok = cbe.bound(&ci[s][u]) as usize <= budget;
                                c &= !(u8::from(!ok) << s);
                            }
                        }
                        live[u] = c;
                        sel[k2 % N] = i;
                        k2 += usize::from(c != 0);
                    }
                    k = k2;
                }
            }
            let lim = if B {
                // a child at depth g+1 with bound h cannot give a solution shorter than g+1+h
                self.prune_len.load(Ordering::Relaxed).saturating_sub(g + 2)
            } else {
                0
            };
            // Without a joint-table stage, CB values that can neither prune nor change the
            // ordering statistic (max(ha, hc) = ha) are read by each child itself, after its
            // own AE stage, when the prefetch issued here has completed (part of the
            // opt/micro-compute stack: 1.29x / 1.28x nodes/s at 48T / 1T as a whole). On a
            // node alive on both sides each side's CB value also enters the other side's
            // class-B bound, which it cannot change either while lim >= cb_max.
            let cb_max = self.cb_max;
            let defer = |s: usize| {
                on(S, s)
                    && !cbe_done
                    && budget > AE_FIRST
                    && budget >= cb_max as usize
                    && min_ha[s] >= cb_max
                    && (!B || S != BOTH || lim >= cb_max as usize)
            };
            let defer = [defer(0), defer(1)];
            for j in 0..k {
                let i = sel[j];
                let u = usize::from(i) % N;
                let mut c = live[u];
                for s in 0..2 {
                    if !on(need, s) {
                        continue;
                    }
                    // (a deferred child carries this node's distance, at most cb_max <= budget)
                    let hc = if defer[s] {
                        h_cb[s]
                    } else {
                        decode(h_cb[s], t.cb.tab.get(cbi[s][u]))
                    };
                    pend[s][u] = if defer[s] { cbi[s][u] } else { NO_CB };
                    hcs[s][u] = hc;
                    if on(S, s) {
                        c &= !(u8::from(hc as usize > budget) << s);
                    } else {
                        has[s][u] = decode(h_ae[s], t.ae.get(idx[s][u]));
                    }
                }
                if B {
                    for s in 0..2 {
                        if on(S, s) {
                            let o = 1 - s;
                            let hb = has[o][u].max(hcs[o][u]);
                            c &= !(u8::from(hb as usize > lim) << s);
                        }
                    }
                }
                live[u] = c;
                sel[n % N] = i;
                n += usize::from(c != 0);
            }
            b_done = true;
        } else {
            for j in 0..k {
                let i = sel[j];
                let u = usize::from(i) % N;
                let mut c = live[u];
                for s in 0..2 {
                    if on(S, s) {
                        let ha = decode(h_ae[s], t.ae.get(idx[s][u]));
                        has[s][u] = ha;
                        c &= !(u8::from(ha as usize > budget) << s);
                    }
                    if on(need, s) {
                        pend[s][u] = NO_CB;
                    }
                }
                live[u] = c;
                sel[n % N] = i;
                n += usize::from(c != 0);
            }
        }
        // The joint corner/B-centre/edge-triple table only adds information near the leaves.
        if budget <= 8 && !cbe_done {
            if let Some(cbe) = &t.cbe {
                for s in 0..2 {
                    if !on(S, s) {
                        continue;
                    }
                    for &i in &sel[..n] {
                        let i = usize::from(i) % N;
                        if S != BOTH || on(live[i], s) {
                            cbe.prefetch_all(f, &ys[s][i], &mut ci[s][i]);
                        }
                    }
                }
                let mut k = 0;
                for j in 0..n {
                    let i = sel[j];
                    let u = usize::from(i) % N;
                    let mut c = live[u];
                    for s in 0..2 {
                        if on(S, s) && (S != BOTH || on(c, s)) {
                            let ok = cbe.bound(&ci[s][u]) as usize <= budget;
                            c &= !(u8::from(!ok) << s);
                        }
                    }
                    live[u] = c;
                    sel[k % N] = i;
                    k += usize::from(c != 0);
                }
                n = k;
            }
        }
        // Class-B bound against the best solution so far (optimal mode 0.29x / 0.67x time at
        // 48T / 1T when added, opt/classb): a child stays alive on side s only if its distance
        // to the other side's target cannot make the solution too long either.
        if B && n > 0 && !b_done {
            let best = self.prune_len.load(Ordering::Relaxed);
            // a child at depth g+1 with bound h cannot give a solution shorter than g+1+h
            let lim = best.saturating_sub(g + 2);
            if S != BOTH {
                // the other side's bounds (on a node alive on both sides they are known)
                let o = usize::from(S == SIDE_A);
                let (fw, fws) = Self::flip_words(f, x[o]);
                for &i in &sel[..n] {
                    let i = usize::from(i) % N;
                    ys[o][i] = f.apply_any::<N>(x[o], fw, fws, self.side_move(o ^ r, ms[i]));
                }
                t.ae_ix
                    .index_sel_pf(f, &ys[o], &sel[..n], &mut idx[o], &t.ae.data);
                for &i in &sel[..n] {
                    let i = usize::from(i) % N;
                    cbi[o][i] = CbTable::index_fast(f, &ys[o][i]);
                    prefetch_at(&t.cb.tab.data, cbi[o][i] >> 5);
                }
                for &i in &sel[..n] {
                    let i = usize::from(i) % N;
                    has[o][i] = decode(h_ae[o], t.ae.get(idx[o][i]));
                    hcs[o][i] = decode(h_cb[o], t.cb.tab.get(cbi[o][i]));
                }
            }
            let mut k = 0;
            for j in 0..n {
                let i = sel[j];
                let u = usize::from(i) % N;
                let mut c = live[u];
                for s in 0..2 {
                    if on(S, s) {
                        let o = 1 - s;
                        let hb = has[o][u].max(hcs[o][u]);
                        c &= !(u8::from(hb as usize > lim) << s);
                    }
                }
                live[u] = c;
                sel[k % N] = i;
                k += usize::from(c != 0);
            }
            n = k;
        }
        // Ordering statistic of the task: how close its subtree gets to the phase-1 targets.
        for &i in &sel[..n] {
            let i = usize::from(i) % N;
            for s in 0..2 {
                if on(S, s) && (S != BOTH || on(live[i], s)) {
                    *stat += STAT_W[has[s][i].max(hcs[s][i]) as usize];
                }
            }
        }
        if budget == 0 {
            // the children are leaves
            for &i in &sel[..n] {
                let i = usize::from(i) % N;
                let m = ms[i] as usize;
                for s in 0..2 {
                    if on(S, s) && on(live[i], s) && has[s][i] == 0 && hcs[s][i] == 0 {
                        path.push(m);
                        self.leaf(s ^ r, path);
                        path.pop();
                    }
                }
                if self.stop.load(Ordering::Relaxed) {
                    return;
                }
            }
            return;
        }
        // (the caller provides three frames per remaining level)
        let Some((trio, frames)) = frames.split_first_mut() else {
            debug_assert!(false, "out of search frames");
            return;
        };
        // Sibling pipeline: child j is finished in `f0` while child j+1 has had its first
        // lookups read (`f1`) and child j+2 its first lookups issued (`f2`).
        let [mut f0, mut f1, mut f2] = trio.each_mut();
        let cb = budget - 1;
        // the class-B limit of the children's children (at depth g + 2)
        let clim = if B {
            self.prune_len.load(Ordering::Relaxed).saturating_sub(g + 3)
        } else {
            0
        };
        let sel = *sel;
        let child = |j: usize| usize::from(sel[j % N]) % N;
        // A child alive on one side is evaluated by the one-side code with that side in slot 0
        // (`r`: the side in slot 0), so a child alive on side 1 only of a node alive on both
        // has its slots swapped: one copy of the one-side code instead of two (hot code
        // 46 -> 31 KB; optimal mode neutral on identical trees, 1T 0.99x / 1.01x time without
        // / with slices in three back-to-back runs; time-limited nodes/s 1.01x / 1.03x at
        // 48T / 1T).
        let kind = |i: usize| {
            if S == BOTH && live[i] == BOTH {
                BOTH
            } else {
                SIDE_A
            }
        };
        let slots = |i: usize| -> (usize, usize, usize) {
            if S == BOTH && live[i] == SIDE_B {
                (1, 1, 0)
            } else {
                (r, 0, 1)
            }
        };
        let xs = |i: usize| {
            let (_, a, b) = slots(i);
            [&ys[a][i], &ys[b][i]]
        };
        let hs = |i: usize| {
            let (_, a, b) = slots(i);
            (
                [has[a][i], has[b][i]],
                [hcs[a][i], hcs[b][i]],
                [pend[a][i], pend[b][i]],
            )
        };
        if n > 0 {
            let i = child(0);
            self.expand_any(kind(i), slots(i).0, xs(i), ms[i] as usize, cb, f0, ctr);
        }
        if n > 1 {
            let i = child(1);
            self.expand_any(kind(i), slots(i).0, xs(i), ms[i] as usize, cb, f1, ctr);
        }
        if n > 0 {
            let i = child(0);
            let (ha, hc, pd) = hs(i);
            self.stage2_any::<B>(kind(i), slots(i).0, xs(i), ha, hc, pd, cb, clim, f0);
        }
        for j in 0..n {
            if j + 2 < n {
                let i = child(j + 2);
                self.expand_any(kind(i), slots(i).0, xs(i), ms[i] as usize, cb, f2, ctr);
            }
            if j + 1 < n {
                let i = child(j + 1);
                let (ha, hc, pd) = hs(i);
                self.stage2_any::<B>(kind(i), slots(i).0, xs(i), ha, hc, pd, cb, clim, f1);
            }
            let i = child(j);
            path.push(ms[i] as usize);
            let r = slots(i).0;
            self.finish_any::<B>(kind(i), r, xs(i), path, bound, ctr, stat, f0, frames);
            path.pop();
            if self.stop.load(Ordering::Relaxed) {
                return;
            }
            (f0, f1, f2) = (f1, f2, f0);
        }
    }

    /// The chirality flip words of a position (the slice word only with slices; see
    /// [`FastTables::apply_any`]).
    #[inline(always)]
    fn flip_words(f: &FastTables, x: &Fc) -> (u64, u32) {
        let em = masks(&x.e);
        let fws = if N > N_MOVES {
            f.flip_word_slice(&em)
        } else {
            0
        };
        (f.flip_word(&em), fws)
    }

    /// All children of a task (unpruned), with their exact bounds on both sides. The
    /// lookups of all children are issued before any is read.
    fn expand_all(&self, t: &Task) -> Vec<Task> {
        let tb = self.t;
        let f = fast();
        let last = t.path.last().copied().unwrap_or(usize::MAX);
        let em = masks(&t.x[0].e);
        let emb = masks(&t.x[1].e);
        let mut ys = [t.x[0]; N];
        let mut yb = [t.x[1]; N];
        let mut ms = [0usize; N];
        let mut n = 0;
        for m in 0..N {
            if !allowed_after(N > N_MOVES, last, m) {
                continue;
            }
            ys[n] = f.apply_m(&t.x[0], &em, m);
            yb[n] = f.apply_m(&t.x[1], &emb, self.mb[m & 31]);
            ms[n] = m;
            n += 1;
        }
        let mut ia = [0usize; N];
        let mut iab = [0usize; N];
        let mut ic = [0usize; N];
        let mut icb = [0usize; N];
        tb.ae_ix
            .index_batch_fast_pf::<N>(f, &ys[..n], &mut ia[..n], &tb.ae.data);
        tb.ae_ix
            .index_batch_fast_pf::<N>(f, &yb[..n], &mut iab[..n], &tb.ae.data);
        for i in 0..n {
            ic[i] = CbTable::index_fast(f, &ys[i]);
            prefetch(&tb.cb.tab.data[ic[i] >> 5]);
            icb[i] = CbTable::index_fast(f, &yb[i]);
            prefetch(&tb.cb.tab.data[icb[i] >> 5]);
        }
        (0..n)
            .map(|i| {
                let mut path = Vec::with_capacity(t.path.len() + 1);
                path.extend_from_slice(&t.path);
                path.push(ms[i]);
                Task {
                    path,
                    x: [ys[i], yb[i]],
                    h_ae: [
                        decode(t.h_ae[0], tb.ae.get(ia[i])),
                        decode(t.h_ae[1], tb.ae.get(iab[i])),
                    ],
                    h_cb: [
                        decode(t.h_cb[0], tb.cb.tab.get(ic[i])),
                        decode(t.h_cb[1], tb.cb.tab.get(icb[i])),
                    ],
                }
            })
            .collect()
    }

    /// Tasks for one iteration (as node ids of `tree`, with the sides each is alive on): the
    /// root, expanded level by level until there are enough tasks for the thread pool. The
    /// depth stays small because most of an IDA* iteration's work lies in the middle levels,
    /// which must not be explored serially. The tree keeps every node evaluated by earlier
    /// iterations with all its children, so a later iteration only evaluates nodes it has not
    /// seen and re-filters. (Splitting: 48T first solution 0.50x when added, opt/split;
    /// memoized tree and adaptive task count: a further 0.51x, opt/split2 stack.)
    fn tasks(&self, tree: &mut Tree, bound: usize, threads: usize, want: usize) -> Vec<(u32, u8)> {
        if tree.nodes.is_empty() {
            tree.nodes.push(self.root.clone());
            tree.kids.push(None);
        }
        // the sides a node is alive on, given its parent's: those where its own bound fits
        // (the task levels use the distance bounds alone; the search applies the rest below)
        let alive = |t: &Task, parent: u8| -> u8 {
            let mut c = 0;
            for s in 0..2 {
                if on(parent, s) && t.h_ae[s].max(t.h_cb[s]) as usize + t.path.len() <= bound {
                    c |= 1 << s;
                }
            }
            c
        };
        let mut front: Vec<(u32, u8)> = vec![];
        let root = alive(&tree.nodes[0], BOTH);
        if root != 0 {
            front.push((0, root));
        }
        let max_k = bound.saturating_sub(4).min(MAX_SPLIT);
        let mut k = 0;
        // (single-threaded runs are split too, so per-task decisions can change mid-iteration)
        while front.len() < want && k < max_k {
            let need: Vec<u32> = front
                .iter()
                .map(|&(i, _)| i)
                .filter(|&i| tree.kids[i as usize].is_none())
                .collect();
            if !need.is_empty() {
                let nodes = &tree.nodes;
                let kids: Vec<Vec<Task>> = if threads > 1 && need.len() >= par_min() {
                    need.par_iter()
                        .map(|&i| self.expand_all(&nodes[i as usize]))
                        .collect()
                } else {
                    need.iter()
                        .map(|&i| self.expand_all(&nodes[i as usize]))
                        .collect()
                };
                for (&i, ch) in need.iter().zip(kids) {
                    let lo = tree.nodes.len() as u32;
                    tree.nodes.extend(ch);
                    tree.kids[i as usize] = Some((lo, tree.nodes.len() as u32));
                }
                tree.kids.resize(tree.nodes.len(), None);
            }
            let mut next = Vec::with_capacity(front.len() * if N > N_MOVES { 19 } else { 13 });
            for &(i, sides) in &front {
                // (every front node has been expanded above)
                if let Some((lo, hi)) = tree.kids[i as usize] {
                    next.extend((lo..hi).filter_map(|c| {
                        let a = alive(&tree.nodes[c as usize], sides);
                        (a != 0).then_some((c, a))
                    }));
                }
            }
            front = next;
            k += 1;
        }
        front
    }
}

/// The top levels of the search tree, kept across the iterations of one solve.
#[derive(Default)]
struct Tree {
    nodes: Vec<Task>,
    /// children of each node (a contiguous range of `nodes`), once evaluated
    kids: Vec<Option<(u32, u32)>>,
}

/// `FTO_FULL_PROOF=1` turns off the early end of optimality proofs (see `solve`), so that
/// `scripts/check-optimal.sh` can compare its results with the full proof.
fn full_proof() -> bool {
    static V: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *V.get_or_init(|| std::env::var("FTO_FULL_PROOF").is_ok_and(|v| v == "1"))
}

fn env_usize(name: &str, default: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(default)
}

/// Fewest nodes to expand in parallel (fewer are expanded on the calling thread).
fn par_min() -> usize {
    static V: std::sync::OnceLock<usize> = std::sync::OnceLock::new();
    *V.get_or_init(|| env_usize("FTO_PARMIN", 256))
}

/// Iterations predicted to evaluate fewer nodes than this run on the calling thread.
fn serial_nodes() -> u64 {
    static V: std::sync::OnceLock<usize> = std::sync::OnceLock::new();
    *V.get_or_init(|| env_usize("FTO_SERIAL", 5000)) as u64
}

/// Predicted nodes per parallel task (fewer tasks are made for smaller iterations).
fn task_nodes() -> u64 {
    static V: std::sync::OnceLock<usize> = std::sync::OnceLock::new();
    *V.get_or_init(|| env_usize("FTO_TASKN", 512)) as u64
}

/// Deepest level at which an iteration is split into parallel tasks.
const MAX_SPLIT: usize = 4;

/// Solve one position with the move set of the tables (`t.moves`). `frame` is the physical
/// frame of the scrambled position.
pub fn solve<I: AeIndex>(
    t: &Tables<I>,
    cube: Cube,
    frame: Frame,
    limits: &Limits,
    report: Report,
) -> Outcome {
    // The search is compiled once per move-set size, so the face metric runs the same code as
    // without slices.
    if t.moves.slice {
        solve_n::<I, MAX_MOVES>(t, cube, frame, limits, report)
    } else {
        solve_n::<I, N_MOVES>(t, cube, frame, limits, report)
    }
}

/// [`solve`] with `N` internal moves (`t.moves.n_moves()`).
fn solve_n<I: AeIndex, const N: usize>(
    t: &Tables<I>,
    cube: Cube,
    frame: Frame,
    limits: &Limits,
    report: Report,
) -> Outcome {
    let start = Instant::now();
    let Some(cube_b) = cube.conjugate(&SWAP_AB) else {
        // not a valid position in the fixed frame: nothing to search
        return Outcome {
            solution: None,
            length: None,
            solutions_found: 0,
            nodes: 0,
            leaves: 0,
            leaves_by_side: [0; 2],
            elapsed: start.elapsed().as_secs_f64(),
            optimal: false,
            events: vec![],
        };
    };
    let x = [fast().of(&cube), fast().of(&cube_b)];
    let (h0, h1) = (t.exact_fast(&x[0]), t.exact_fast(&x[1]));
    let root = Task {
        path: vec![],
        x,
        h_ae: [h0.0, h1.0],
        h_cb: [h0.1, h1.1],
    };
    // Every solution is a prefix followed by class-A moves, so one side would suffice to
    // prove optimality; both are kept because the earlier a short solution is found, the more
    // the class-B bound prunes (measured: one variant was 1.9-3.9x slower in optimal mode),
    // and because with both the proof can stop one iteration early (see the main loop).
    // The two sides are searched in one traversal (see `Shared::dfs`). Speedup when added
    // (opt/merged, against main, two interleaved rounds, large tier): optimal-mode proofs
    // 0.80x / 0.76x time at 48T / 1T (1T: 0.70x nodes; 30 more 17-move scrambles: 0.85x),
    // random positions 0.79x / 0.90x / 0.93x (48T), with slices 0.96x / 0.85x; time-limited
    // nodes/s 0.95x / 0.93x (the shared nodes count once).
    let sh = Shared::<I, N> {
        t,
        sides: [
            Side {
                cube,
                sym: IDENTITY,
            },
            Side {
                cube: cube_b,
                sym: SWAP_AB,
            },
        ],
        root,
        frame,
        limits: limits.clone(),
        report,
        start,
        best_len: AtomicUsize::new(usize::MAX),
        prune_len: AtomicUsize::new(usize::MAX),
        count: AtomicU64::new(0),
        proof_len: AtomicUsize::new(0),
        proven: AtomicBool::new(false),
        nodes: AtomicU64::new(0),
        nodes_both: AtomicU64::new(0),
        leaves: [AtomicU64::new(0), AtomicU64::new(0)],
        stop: AtomicBool::new(false),
        best: Mutex::new((None, vec![])),
        mb: swap_moves::<N>(),
        mv: [std::array::from_fn(|m| m), swap_moves::<N>()],
        cb_max: (t.cb.tab.hist.len() - 1) as u8,
    };
    let lo = (0..2)
        .map(|s| sh.root.h_ae[s].max(sh.root.h_cb[s]) as usize)
        .fold(usize::MAX, usize::min);
    let mut optimal = false;
    let threads = rayon::current_num_threads();
    let mut depth = lo;
    let mut tree = Tree::default();
    // ordering statistic of each tree node from the previous iteration (0 if not a task then)
    let mut node_stat: Vec<f32> = vec![];
    let mut prev_ids: Vec<(u32, u8)> = vec![];
    let mut prev_nodes: Option<u64> = None;
    // Target-length passes: only in optimal mode, when the best is at least this far above the
    // lower bound (0 = never). On by default for one thread only: there the best stays loose
    // for long and the passes cut the proof time (1T 17-move set 11.95 -> 5.27 s, 0.44x, when
    // added); with 48 threads a tight best is found early and the rerun iterations cost more
    // (19-move set 1.07x, random positions 1.00-1.04x). FTO_TARGET_GAP overrides. Off with
    // FTO_FULL_PROOF: the passes rest on the same argument as the early proof end.
    // The argument of the early proof end below: a solving phase-1 path of length best - 1 that
    // ends in a class-B move is found earlier by the class-B side. It fails with slices: a
    // path can end in a slice, which is in neither class, so no side finds it earlier.
    let class_b_argument = !t.moves.slice;
    let target_gap = if limits.is_optimal() && !full_proof() {
        env_usize("FTO_TARGET_GAP", if threads == 1 { 2 } else { 0 })
    } else {
        0
    };
    let mut target: Option<usize> = None;
    loop {
        // Optimality is proven one iteration before depth reaches the best length. The
        // iteration at depth best - 1 could only improve on the best with a leaf that is solved
        // (phase-2 length 0): a path of that length that ends in a class-B move and solves the
        // puzzle. Such a path is a prefix ending in a class-A move followed by class-B moves, so
        // the class-B side already found it at a smaller depth, in an iteration completed
        // before this one. (This needs both sides in every iteration: the merged traversal
        // visits the union of their trees.) When solutions are counted, that iteration still
        // yields completions (longer ones), so it is run. With slices the argument fails (see
        // `class_b_argument`), and the full proof runs (checked by scripts/check-optimal.sh and
        // the brute-force tests).
        let early = limits.max_solutions.is_none() && !full_proof() && class_b_argument;
        let proof_depth = depth + usize::from(early);
        if proof_depth >= sh.best_len.load(Ordering::Relaxed) {
            optimal = true;
            break;
        }
        if limits.max_depth.is_some_and(|d| depth > d) {
            break;
        }
        // Target-length passes (optimal mode, FTO_TARGET_GAP): pass T runs the depths lo..T-1
        // with the class-B bound set for solutions of at most T moves instead of against a
        // loose best, which prunes far more. By the argument above every solution of at most
        // T moves is found at a depth below T, so a completed pass proves that none exists.
        // With slices the pass also runs depth T (a solution of T moves may need a phase-1
        // path of T moves ending in a slice).
        if let Some(t) = target {
            if depth >= t + usize::from(!class_b_argument) {
                // pass T complete: no solution of at most T moves
                if sh.best_len.load(Ordering::Relaxed) <= t + 1 {
                    optimal = true;
                    break;
                }
                target = Some(t + 1);
                sh.prune_len.store(
                    (t + 2).min(sh.best_len.load(Ordering::Relaxed)),
                    Ordering::Relaxed,
                );
                depth = lo;
                prev_nodes = None;
                continue;
            }
            // a solution of at most T moves found in pass T is optimal at once
            sh.proof_len.store(t, Ordering::Relaxed);
        } else {
            // the iterations below `depth` are complete: no solution of at most `depth` moves
            // (of at most `depth - 1` with slices)
            let lower = depth + usize::from(class_b_argument);
            if target_gap > 0
                && depth > lo
                && sh.best_len.load(Ordering::Relaxed) >= lower + target_gap
            {
                target = Some(lower);
                sh.prune_len.store(lower + 1, Ordering::Relaxed);
                if std::env::var_os("FTO_DEBUG").is_some() {
                    eprintln!(
                        "  target-length passes from {lower} (best {:?})",
                        sh.best_len
                    );
                }
                depth = lo;
                prev_nodes = None;
                continue;
            }
            // For the same reason a solution of length depth + 1 found during this iteration
            // is optimal at once (no shorter one exists after the iterations below), and the
            // search stops there instead of finishing the iteration.
            sh.proof_len
                .store(if early { depth + 1 } else { 0 }, Ordering::Relaxed);
        }
        // Small iterations run on this thread: waking the pool costs more than they do.
        let predicted = prev_nodes.map_or(if depth < 14 { 0 } else { u64::MAX }, |n| {
            n.saturating_mul(16)
        });
        let inline = threads == 1 || predicted < serial_nodes();
        // Enough tasks to balance the load, but no more than the iteration's work warrants
        // (each extra split level costs a parallel expansion of the tree).
        let want = if threads == 1 {
            64
        } else if inline {
            1
        } else {
            ((predicted / task_nodes()) as usize).clamp(threads, 64 * threads)
        };
        let mut ids = sh.tasks(&mut tree, depth, threads, want);
        // Order the tasks by the statistics of the same prefixes in the previous iteration, so
        // that the most promising subtrees are searched first (egg-solver's job ordering;
        // 48T 5 s mean length -0.225 moves, 1T optimal 0.69x time when added, opt/order).
        node_stat.resize(tree.nodes.len(), 0.0);
        if !prev_ids.is_empty() {
            ids.sort_by(|a, b| node_stat[b.0 as usize].total_cmp(&node_stat[a.0 as usize]));
        }
        // Pruning against the best length changes how many (non-improving) solutions are
        // reached, so it is off when the user counts solutions.
        let count_mode = limits.max_solutions.is_some();
        let t_split = start.elapsed().as_secs_f64();
        let stats: Vec<std::sync::atomic::AtomicU32> = (0..ids.len())
            .map(|_| std::sync::atomic::AtomicU32::new(0))
            .collect();
        let next = std::sync::atomic::AtomicUsize::new(0);
        let nodes_before = sh.nodes.load(Ordering::Relaxed);
        let tree_ref = &tree;
        let worker = || {
            let mut path: Vec<usize> = Vec::with_capacity(64);
            // one frame per level below a task (tasks are at depth >= 0)
            let mut frames = vec![[Scratch::<N>::new(), Scratch::new(), Scratch::new()]; depth + 1];
            // the per-thread counters persist across tasks so time checks happen regularly
            let mut ctr = Counter::default();
            loop {
                let i = next.fetch_add(1, Ordering::Relaxed);
                if i >= ids.len() || sh.stop.load(Ordering::Relaxed) || sh.time_up() {
                    break;
                }
                let (id, sides) = ids[i];
                let task = &tree_ref.nodes[id as usize];
                path.clear();
                path.extend_from_slice(&task.path);
                let (before, shared_before) = (ctr.nodes, ctr.shared);
                let mut stat = 0f32;
                let use_b = !count_mode && sh.prune_len.load(Ordering::Relaxed) < usize::MAX;
                if use_b {
                    sh.dfs::<true>(
                        task,
                        sides,
                        &mut path,
                        depth,
                        &mut ctr,
                        &mut stat,
                        &mut frames,
                    );
                } else {
                    sh.dfs::<false>(
                        task,
                        sides,
                        &mut path,
                        depth,
                        &mut ctr,
                        &mut stat,
                        &mut frames,
                    );
                }
                stats[i].store(stat.to_bits(), Ordering::Relaxed);
                sh.nodes.fetch_add(ctr.nodes - before, Ordering::Relaxed);
                sh.nodes_both
                    .fetch_add(ctr.shared - shared_before, Ordering::Relaxed);
            }
        };
        if inline {
            worker();
        } else {
            // Workers take tasks strictly in order from a shared counter.
            rayon::scope(|scope| {
                for _ in 0..threads.min(ids.len()) {
                    scope.spawn(|_| worker());
                }
            });
        }
        prev_nodes = Some(sh.nodes.load(Ordering::Relaxed) - nodes_before);
        for &(i, _) in &prev_ids {
            node_stat[i as usize] = 0.0;
        }
        for (i, &(id, _)) in ids.iter().enumerate() {
            node_stat[id as usize] = f32::from_bits(stats[i].load(Ordering::Relaxed));
        }
        let n_tasks = ids.len();
        prev_ids = ids;
        if std::env::var_os("FTO_DEBUG").is_some() {
            eprintln!("  depth {depth}: split done at {t_split:.4}s");
            eprintln!(
                "  depth {depth}: {} tasks, {} nodes ({} shared by both sides), {} solutions, best {:?}, {:.3}s",
                n_tasks,
                sh.nodes.load(Ordering::Relaxed),
                sh.nodes_both.load(Ordering::Relaxed),
                sh.count.load(Ordering::Relaxed),
                sh.best_len.load(Ordering::Relaxed),
                start.elapsed().as_secs_f64()
            );
        }
        if sh.proven.load(Ordering::Relaxed) {
            optimal = true;
            break;
        }
        if sh.stop.load(Ordering::Relaxed) || sh.time_up() {
            break;
        }
        depth += 1;
    }
    let (sol, events) = sh.best.into_inner().unwrap_or_else(PoisonError::into_inner);
    let length = sol.as_ref().map(Vec::len);
    if report == Report::Best {
        if let Some(s) = &sol {
            println!("{}", format_moves(s));
        }
    }
    let leaves_by_side = sh.leaves.map(AtomicU64::into_inner);
    Outcome {
        solution: sol,
        length,
        solutions_found: sh.count.load(Ordering::Relaxed),
        nodes: sh.nodes.load(Ordering::Relaxed),
        leaves: leaves_by_side[0] + leaves_by_side[1],
        leaves_by_side,
        elapsed: start.elapsed().as_secs_f64(),
        optimal,
        events,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn move_masks_match_rules() {
        for moves in [
            MoveSet::FACE,
            MoveSet {
                wide: false,
                slice: true,
            },
        ] {
            let tab = &ALLOWED[usize::from(moves.slice)];
            let n = moves.n_moves();
            for last in (0..n).chain([MAX_MOVES]) {
                for m in 0..MAX_MOVES {
                    let l = if last == MAX_MOVES { usize::MAX } else { last };
                    assert_eq!(tab[last] >> m & 1 == 1, m < n && moves.allowed_after(l, m));
                }
            }
        }
        for m in 0..MAX_MOVES {
            assert_eq!(A_MOVES >> m & 1 == 1, is_a_move(m));
            // the class-B moves are the moves whose image under SWAP_AB is a class-A move
            assert_eq!(
                B_MOVES >> m & 1 == 1,
                is_a_move(conjugate_move(m, &SWAP_AB)),
                "move {m}"
            );
            let fin = FINAL_SIDES[m];
            assert_eq!(on(fin, 0), !is_a_move(m));
            assert_eq!(on(fin, 1), !is_a_move(conjugate_move(m, &SWAP_AB)));
        }
    }

    /// The merged traversal walks side 1 through the images of side 0's moves, so the two
    /// sides' canonical sequences must correspond: `SWAP_AB` keeps the move order.
    #[test]
    fn class_swap_keeps_the_move_order() {
        for slice in [false, true] {
            let n = if slice { MAX_MOVES } else { N_MOVES };
            let mb = |m: usize| {
                if m < n {
                    conjugate_move(m, &SWAP_AB)
                } else {
                    m
                }
            };
            for last in (0..n).chain([MAX_MOVES]) {
                for m in 0..n {
                    assert_eq!(
                        allowed_after(slice, last, m),
                        allowed_after(slice, mb(last), mb(m)),
                        "{last} {m}"
                    );
                }
            }
            // a bijection on the moves
            let mut seen = vec![false; n];
            for m in 0..n {
                assert!(mb(m) < n && !std::mem::replace(&mut seen[mb(m)], true));
            }
        }
    }
}
