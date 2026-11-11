//! IDA* straight to the solved state (optimal solver).
//!
//! Iterations with depth bound `D = h(root), h(root) + 1, ...`: a depth-first search over the
//! canonical move sequences (`crate::cube::allowed_after`) that skips a child whose lower bound
//! exceeds the moves left. Iteration `D` finds exactly the solutions of `D` moves, after the
//! iterations below it have found none, so the first solution is optimal and a completed
//! iteration proves that no solution of at most `D` moves exists.
//!
//! The top levels of each iteration are split into tasks (move prefixes of up to
//! [`MAX_SPLIT`] moves) that the threads take in order from a shared counter; the tasks are
//! ordered by how close their subtrees came to solved in the previous iteration, so that the
//! last iteration tends to meet a solution early.
//!
//! Each node generates its children, computes all their table indices and issues the
//! prefetches, then reads the values (so the cache misses of the children overlap), and
//! recurses into the children that pass. Values are exact distances mod 3, decoded with the
//! parent's exact distances.

use super::coords::*;
use super::tables::{exact, Bounds, MAX_K};
use crate::cube::*;
use crate::notation::{solution_to_physical, Frame, PhysMove};
use crate::prune::decode;
use crate::tables::prefetch_at;
use std::collections::HashMap;
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::sync::{Mutex, PoisonError};
use std::time::{Duration, Instant};

/// Stopping rules. Without any, the search returns the first optimal solution.
#[derive(Clone, Debug, Default)]
pub struct Limits {
    /// return up to this many solutions of the optimal length (all of them if fewer)
    pub max_solutions: Option<u64>,
    /// stop after this much time; the outcome then holds the lower bound proven so far
    pub time: Option<Duration>,
    /// give up when the optimum is longer than this (no iteration deeper than this runs)
    pub max_length: Option<usize>,
}

/// Why a search ended.
#[derive(Clone, Copy, Debug, PartialEq, Eq, serde::Serialize)]
pub enum Stop {
    /// solutions of the optimal length found
    Solved,
    /// the time limit ran out first
    Time,
    /// the optimum is longer than `Limits::max_length`
    MaxLength,
    /// the tables are not distance tables for the position (never with checked tables)
    BadTables,
}

/// One completed (or interrupted) iteration.
#[derive(Clone, Debug, serde::Serialize)]
pub struct Iteration {
    pub depth: usize,
    /// nodes generated in this iteration
    pub nodes: u64,
    /// table lookups in this iteration
    pub lookups: u64,
    pub tasks: usize,
    /// time since the start, at its end
    pub elapsed: f64,
}

#[derive(Clone, Debug)]
pub struct Outcome {
    /// solutions found, all of the optimal length (internal moves of the fixed frame)
    pub internal: Vec<Vec<usize>>,
    /// the same as physical moves
    pub solutions: Vec<Vec<PhysMove>>,
    /// the optimal length, if a solution was found
    pub length: Option<usize>,
    /// no solution has fewer moves (proven by the completed iterations)
    pub lower_bound: usize,
    pub stop: Stop,
    /// nodes generated (children of expanded nodes)
    pub nodes: u64,
    /// table lookups (each a read of one table entry, after computing its index)
    pub lookups: u64,
    /// children whose lookups were read, by the set of lookups (bit k: lookup k) whose value
    /// exceeds the moves left after the child: with early exits in a given order, a child
    /// costs the lookups up to the first one of its set (all `K` if the set is empty)
    pub prune_masks: Vec<u64>,
    /// children generated, by the moves left after them
    pub children_by_left: Vec<u64>,
    pub elapsed: f64,
    pub iterations: Vec<Iteration>,
    /// how many of the 48 symmetries fix the position (1: only the identity)
    pub symmetries: usize,
}

/// Moves allowed after each move (index `MAX_MOVES`: no previous move) by
/// [`allowed_after`], as bit masks over the moves of the move set: `[0]` face turns, `[1]` with
/// slices.
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

/// Weight 8^-h of a surviving child with bound h in a task's ordering statistic.
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

/// Symmetric start (see [`Shared::dedup`]); `FTO_OPT_SYMSTART=0` turns it off (to compare).
fn sym_start() -> bool {
    std::env::var("FTO_OPT_SYMSTART").map_or(true, |v| v != "0")
}

/// Deepest level of the task split.
pub const MAX_SPLIT: usize = 4;

/// The shared stop flag is read every this many nodes (and after each solution).
const CHECK_EVERY: u64 = 1 << 14;

/// Iterations predicted to generate fewer nodes than this run on the calling thread.
const SERIAL_NODES: u64 = 50_000;

/// Predicted nodes per parallel task.
const TASK_NODES: u64 = 1 << 16;

/// A node of the top levels (a task).
#[derive(Clone, Copy)]
struct Task {
    path: [u8; MAX_SPLIT],
    len: u8,
    p: Pos,
    /// exact distances of the lookups
    h: [u8; MAX_K],
}

impl Task {
    fn last(&self) -> usize {
        if self.len == 0 {
            MAX_MOVES
        } else {
            usize::from(self.path[usize::from(self.len) - 1])
        }
    }
}

/// Key of a move prefix (for the ordering statistics).
fn key(path: &[u8]) -> u32 {
    path.iter()
        .rev()
        .fold(0, |k, &m| k * (MAX_MOVES as u32 + 1) + u32::from(m) + 1)
}

/// Per-worker counters.
#[derive(Default)]
struct Counter {
    nodes: u64,
    lookups: u64,
    /// pruning statistics (see [`Outcome::prune_masks`])
    stats: PruneStats,
    /// `Shared::stop` as last read
    stop: bool,
}

/// Which lookups prune the children, for choosing a lookup order (see
/// [`Outcome::prune_masks`] and [`Outcome::children_by_left`]).
#[derive(Clone, Debug, Default)]
struct PruneStats {
    /// children by the set of lookups whose value exceeds the moves left (bit k: lookup k)
    masks: Vec<u64>,
    /// children by the moves left after them
    left: Vec<u64>,
}

impl PruneStats {
    #[inline(always)]
    fn mask(&mut self, m: usize) {
        if self.masks.len() <= m {
            self.masks.resize(1 << MAX_K, 0);
        }
        self.masks[m] += 1;
    }
    #[inline(always)]
    fn left(&mut self, left: usize, n: usize) {
        if self.left.len() <= left {
            self.left.resize(left + 1, 0);
        }
        self.left[left] += n as u64;
    }
    fn add(&mut self, o: &PruneStats) {
        for (v, w) in [(&mut self.masks, &o.masks), (&mut self.left, &o.left)] {
            if v.len() < w.len() {
                v.resize(w.len(), 0);
            }
            for (a, b) in v.iter_mut().zip(w) {
                *a += b;
            }
        }
    }
}

/// Scratch arrays of one node's children: positions, indices, values and survivors (one frame
/// per depth, reused).
#[derive(Clone)]
struct Eval<const N: usize> {
    ys: [Pos; N],
    ms: [u8; N],
    idx: [[usize; MAX_K]; N],
    hs: [[u8; MAX_K]; N],
    /// each child's bound (the largest of its values read so far)
    mx: [u8; N],
    /// each child's [`Bounds::pre`] value of the first stage
    pre: [u32; N],
    /// the order of this frame's stages ([`stage_order`])
    order: [u8; MAX_K],
    sel: [u8; N],
}

impl<const N: usize> Eval<N> {
    fn new() -> Self {
        let z = Pos::solved();
        Eval {
            ys: [z; N],
            ms: [0; N],
            idx: [[0; MAX_K]; N],
            hs: [[0; MAX_K]; N],
            mx: [0; N],
            pre: [0; N],
            order: [0; MAX_K],
            sel: [0; N],
        }
    }
}

/// The state shared by the threads of one solve, with `N` internal moves (16, or 24 with
/// slices).
struct Shared<'a, B: Bounds, const N: usize> {
    b: &'a B,
    /// `opt_moves()`, read once (no `OnceLock` access per node)
    mt: &'static OptMoves,
    solved: Pos,
    cube: Cube,
    limits: Limits,
    start: Instant,
    stop: AtomicBool,
    nodes: AtomicU64,
    lookups: AtomicU64,
    stats: Mutex<PruneStats>,
    /// read every lookup of every child and record which ones prune it (no early exit;
    /// `FTO_OPT_PRUNE_STATS=1`)
    prune_stats: bool,
    sols: Mutex<Vec<Vec<usize>>>,
    /// symmetric start: the symmetries (other than the identity) that fix the root, if the
    /// task levels are filtered (see [`Shared::dedup`]); empty otherwise
    symmetries: Vec<usize>,
}

impl<B: Bounds, const N: usize> Shared<'_, B, N> {
    fn time_up(&self) -> bool {
        if self.limits.time.is_some_and(|t| self.start.elapsed() >= t) {
            self.stop.store(true, Ordering::Relaxed);
            return true;
        }
        false
    }

    #[inline(always)]
    fn count(&self, ctr: &mut Counter, n: usize) {
        let before = ctr.nodes;
        ctr.nodes += n as u64;
        if before / CHECK_EVERY != ctr.nodes / CHECK_EVERY {
            self.time_up();
            ctr.stop = self.stop.load(Ordering::Relaxed);
        }
    }

    /// The most solutions to return.
    fn max_solutions(&self) -> usize {
        self.limits
            .max_solutions
            .map_or(1, |n| usize::try_from(n).unwrap_or(usize::MAX).max(1))
    }

    /// A solution (`path`, internal moves).
    fn found(&self, path: &[u8], ctr: &mut Counter) {
        let moves: Vec<usize> = path.iter().map(|&m| usize::from(m)).collect();
        // (the tables only bound the distance; the leaf test compares whole positions, and the
        // solution is checked on the cubie model as well)
        if self.cube.apply_seq(&moves) != Cube::solved() {
            debug_assert!(false, "a leaf passed as solved but the moves do not solve");
            return;
        }
        let mut sols = self.sols.lock().unwrap_or_else(PoisonError::into_inner);
        if sols.len() < self.max_solutions() {
            sols.push(moves);
        }
        if sols.len() >= self.max_solutions() {
            self.stop.store(true, Ordering::Relaxed);
            ctr.stop = true;
        }
    }

    /// Children of task `t` (the moves allowed after its last move) that pass the bound, with
    /// their exact distances.
    fn expand(&self, t: &Task, bound: usize, ctr: &mut Counter) -> Vec<Task> {
        let mut out = Vec::with_capacity(N);
        let mut moves = ALLOWED[usize::from(N > N_MOVES)][t.last()];
        let g = usize::from(t.len);
        while moves != 0 {
            let m = moves.trailing_zeros() as usize;
            moves &= moves - 1;
            let mut c = *t;
            c.path[g % MAX_SPLIT] = m as u8;
            c.len += 1;
            c.p = self.mt.apply_n::<N>(&t.p, m);
            out.push(c);
        }
        let n = out.len();
        ctr.nodes += n as u64;
        ctr.stats.left(bound - g - 1, n);
        let mut ev = Eval::<N>::new();
        for (i, c) in out.iter().enumerate() {
            ev.ys[i % N] = c.p;
        }
        let k = self.bound_children(&mut ev, n, t.h, bound - g - 1, ctr, false);
        let alive: Vec<Task> = ev.sel[..k]
            .iter()
            .map(|&i| {
                let i = usize::from(i) % N;
                let mut c = out[i];
                c.h = ev.hs[i];
                c
            })
            .collect();
        alive
    }

    /// The bounds of children `ev.ys[..n]` (the parent's exact values `h`, `budget` moves left
    /// after them): the tier's lookups in its order (`Bounds::ORDER`), each read only for the
    /// children that every earlier one let through (early exit); each stage computes and
    /// prefetches the indices of all those children before reading any value. Leaves the
    /// survivors (in child order, so the search visits the same tree whatever the order) in
    /// `ev.sel` and returns their number; a survivor has all its values in `ev.hs` and its
    /// bound in `ev.mx`. With `FTO_OPT_PRUNE_STATS=1` every lookup is read for every child, to
    /// record which ones prune it (see [`Outcome::prune_masks`]).
    #[inline(always)]
    fn bound_children(
        &self,
        ev: &mut Eval<N>,
        n: usize,
        h: [u8; MAX_K],
        budget: usize,
        ctr: &mut Counter,
        ready: bool,
    ) -> usize {
        let mut alive = n;
        for i in 0..n {
            ev.sel[i % N] = i as u8;
            ev.mx[i % N] = 0;
        }
        if !ready {
            ev.order = stage_order::<B>(&h, budget);
            self.first_stage(ev, n, false);
        }
        let order = ev.order;
        let order = &order[..B::K.min(MAX_K)];
        for (s, &k) in order.iter().enumerate() {
            let k = usize::from(k) % MAX_K;
            ctr.lookups += alive as u64;
            let tab = self.b.table(k);
            // (the first stage's indices are computed and prefetched already, the later
            // stages' as soon as a child passes the stage before)
            let nk = order.get(s + 1).map(|&x| usize::from(x) % MAX_K);
            let mut next = 0;
            for j in 0..alive {
                let i = ev.sel[j % N];
                let u = usize::from(i) % N;
                let v = decode(h[k], tab.get(ev.idx[u][k]));
                ev.hs[u][k] = v;
                ev.mx[u] = ev.mx[u].max(v);
                ev.sel[next % N] = i;
                let keep = self.prune_stats || usize::from(v) <= budget;
                if let (true, Some(nk)) = (keep, nk) {
                    let ix = self.b.index(nk, &ev.ys[u]);
                    ev.idx[u][nk] = ix;
                    prefetch_at(&self.b.table(nk).data, ix >> 5);
                }
                next += usize::from(keep);
            }
            alive = next;
        }
        if self.prune_stats {
            // every child read every lookup: record its pruning set, then keep the survivors
            let mut next = 0;
            for j in 0..alive {
                let i = ev.sel[j % N];
                let u = usize::from(i) % N;
                let mut mask = 0;
                for k in 0..B::K {
                    mask |= usize::from(usize::from(ev.hs[u][k]) > budget) << k;
                }
                ctr.stats.mask(mask);
                ev.sel[next % N] = i;
                next += usize::from(mask == 0);
            }
            alive = next;
        }
        alive
    }

    /// The first stage's indices of children `ev.ys[..n]` (in child order), their entries
    /// prefetched.
    /// (`have_pre`: `ev.pre` holds the children's [`Bounds::pre`] values of this lookup.)
    #[inline(always)]
    fn first_stage(&self, ev: &mut Eval<N>, n: usize, have_pre: bool) {
        if B::K == 0 {
            return;
        }
        let k = usize::from(ev.order[0]) % MAX_K;
        let tab = self.b.table(k);
        for u in 0..n {
            let y = &ev.ys[u % N];
            let pre = if have_pre {
                ev.pre[u % N]
            } else {
                self.b.pre(k, y)
            };
            let ix = self.b.index_with(k, y, pre);
            ev.idx[u % N][k] = ix;
            prefetch_at(&tab.data, ix >> 5);
        }
    }

    /// Tasks of an iteration with depth bound `bound`: the root expanded level by level (at
    /// most `MAX_SPLIT` levels, and not below the bound) until there are at least `want`.
    fn tasks(&self, root: &Task, bound: usize, want: usize, ctr: &mut Counter) -> Vec<Task> {
        let mut front = vec![*root];
        let mut k = 0;
        while front.len() < want && k < MAX_SPLIT && k < bound {
            front = front
                .iter()
                .flat_map(|t| self.expand(t, bound, ctr))
                .collect();
            if !self.symmetries.is_empty() {
                front = self.dedup(front);
            }
            k += 1;
        }
        front
    }

    /// Symmetric start: of the tasks of one level (in the order of their move prefixes) whose
    /// positions are conjugate under the root's symmetries, keep the first. Conjugate positions
    /// are equally far from solved, so only one optimal solution is needed per class (the
    /// filter is off when `-n` asks for more):
    /// let `g` be the smallest prefix of its length of a canonical optimal solution, and `q`
    /// a smaller task conjugate to it; `q` followed by an optimal solution of its position is
    /// optimal, and its canonical form starts with `q`, or with `q`'s last move replaced by a
    /// smaller one (a wide turn followed by the single turn of its axis is reordered, single
    /// first; a merge would give a shorter solution): a smaller prefix of a canonical optimal
    /// solution than `g`, a contradiction. So `g` is the first of its class, at every level.
    /// The key of a class is the least conjugate under the root's symmetry group (a group:
    /// the stabiliser of the root).
    fn dedup(&self, tasks: Vec<Task>) -> Vec<Task> {
        let sym = super::sym::sym48();
        let key = |t: &Task| {
            let c = t.p.to_cube();
            self.symmetries
                .iter()
                .map(|&k| sym.conj_by_maps(&c, k))
                .chain([c])
                .map(|x| (x.cp, x.co, x.ep, x.xa, x.xb))
                .min()
        };
        let mut seen = std::collections::HashSet::new();
        tasks.into_iter().filter(|t| seen.insert(key(t))).collect()
    }

    /// Generate the children of `p` (the moves allowed after `last`) into `fr`, and when they
    /// have moves left after them (`budget > 0`) compute and prefetch their first-stage
    /// entries. Returns the number of children. Nothing is counted here: a node's children are
    /// counted when the search enters it ([`Shared::enter`]), so the counts do not depend on
    /// how early its frame was prepared.
    #[inline(always)]
    fn prepare(
        &self,
        p: &Pos,
        h: &[u8; MAX_K],
        last: usize,
        budget: usize,
        fr: &mut Eval<N>,
    ) -> usize {
        let mut moves = ALLOWED[usize::from(N > N_MOVES)][last.min(MAX_MOVES)];
        let mut n = 0;
        // (with lookups left, the first stage's Bounds::pre in the same loop: it prefetches
        // what the first-stage index reads first)
        fr.order = stage_order::<B>(h, budget);
        let k0 = usize::from(fr.order[0]) % MAX_K;
        while moves != 0 {
            let m = moves.trailing_zeros() as usize;
            moves &= moves - 1;
            let y = self.mt.apply_n::<N>(p, m);
            if budget > 0 {
                fr.pre[n % N] = self.b.pre(k0, &y);
            }
            fr.ys[n % N] = y;
            fr.ms[n % N] = m as u8;
            n += 1;
        }
        if budget > 0 {
            self.first_stage(fr, n, true);
        }
        n
    }

    /// Depth-first search below a node at depth `g` (its distances `h`, `last` the move that
    /// led to it) with depth bound `bound`. `path` holds the moves from the root.
    fn dfs(
        &self,
        p: &Pos,
        h: [u8; MAX_K],
        last: usize,
        bound: usize,
        path: &mut Vec<u8>,
        ctr: &mut Counter,
        frames: &mut [[Eval<N>; 2]],
        stat: &mut f32,
    ) {
        let g = path.len();
        if g >= bound {
            if *p == self.solved {
                self.found(path, ctr);
            }
            return;
        }
        let Some((pair, _)) = frames.split_first_mut() else {
            debug_assert!(false, "out of search frames");
            return;
        };
        let n = self.prepare(p, &h, last, bound - g - 1, &mut pair[0]);
        self.enter(n, h, bound, path, ctr, frames, 0, stat);
    }

    /// The search below a node at depth `g` < `bound` whose `n` children are prepared in
    /// `frames[0][slot]` ([`Shared::prepare`]). Sibling pipeline: before the search descends
    /// into a surviving child, the next surviving child's children are generated and their
    /// first-stage entries prefetched (in the other frame of the level below), so those misses
    /// overlap the current child's subtree, which near the leaves is a few nodes. The same
    /// tree in the same order, so the same node counts.
    fn enter(
        &self,
        n: usize,
        h: [u8; MAX_K],
        bound: usize,
        path: &mut Vec<u8>,
        ctr: &mut Counter,
        frames: &mut [[Eval<N>; 2]],
        slot: usize,
        stat: &mut f32,
    ) {
        if ctr.stop {
            return;
        }
        let g = path.len();
        // (the caller provides two frames per remaining level)
        let Some((pair, rest)) = frames.split_first_mut() else {
            debug_assert!(false, "out of search frames");
            return;
        };
        let fr = &mut pair[slot & 1];
        let budget = bound - g - 1;
        self.count(ctr, n);
        ctr.stats.left(budget, n);
        if budget == 0 {
            // The children are leaves: every table is 0 exactly at the solved state, so only a
            // solved child can pass; it is recognised without lookups.
            for i in 0..n {
                if fr.ys[i % N] == self.solved {
                    path.push(fr.ms[i % N]);
                    self.found(path, ctr);
                    path.pop();
                }
            }
            return;
        }
        let kk = self.bound_children(fr, n, h, budget, ctr, true);
        for j in 0..kk {
            let u = usize::from(fr.sel[j % N]) % N;
            *stat += STAT_W[usize::from(fr.mx[u]) & 31];
        }
        if kk == 0 {
            return;
        }
        // the survivors, at depth g + 1 < bound: their children have budget - 1 moves left
        let child = |j: usize| usize::from(fr.sel[j % N]) % N;
        let Some((cpair, _)) = rest.split_first_mut() else {
            debug_assert!(false, "out of search frames");
            return;
        };
        let i0 = child(0);
        let mut n_cur = self.prepare(
            &fr.ys[i0],
            &fr.hs[i0],
            usize::from(fr.ms[i0]),
            budget - 1,
            &mut cpair[0],
        );
        for j in 0..kk {
            let n_next = if j + 1 < kk {
                let i1 = child(j + 1);
                let Some((cpair, _)) = rest.split_first_mut() else {
                    return;
                };
                self.prepare(
                    &fr.ys[i1],
                    &fr.hs[i1],
                    usize::from(fr.ms[i1]),
                    budget - 1,
                    &mut cpair[(j + 1) & 1],
                )
            } else {
                0
            };
            let i = child(j);
            path.push(fr.ms[i]);
            self.enter(n_cur, fr.hs[i], bound, path, ctr, rest, j & 1, stat);
            path.pop();
            if ctr.stop {
                return;
            }
            n_cur = n_next;
        }
    }
}

/// The order in which the children of a node with exact values `h` and `budget` moves left
/// after them read their lookups: first the lookups at the node's limit (`h[k] == budget + 1`,
/// the most a passing node can have: they prune every child whose value does not drop), then
/// those one below (`h[k] == budget`: they prune the children whose value rises), each group
/// in the tier's order ([`Bounds::ORDER`]); last the lookups that cannot prune any child (a
/// child's value is at most one above its parent's), whose exact values are only needed for
/// the children that pass the rest. The survivors and their values are the same in any order.
#[inline(always)]
fn stage_order<B: Bounds>(h: &[u8; MAX_K], budget: usize) -> [u8; MAX_K] {
    let mut o = [0u8; MAX_K];
    let mut n = 0;
    let order = &B::ORDER[..B::K.min(B::ORDER.len())];
    // group of lookup k: 0 at the limit, 1 one below, 2 cannot prune
    let group = |k: usize| {
        (budget + 1)
            .saturating_sub(usize::from(h[k % MAX_K]))
            .min(2)
    };
    for g in 0..3 {
        for &k in order {
            if group(k) == g {
                o[n % MAX_K] = k as u8;
                n += 1;
            }
        }
    }
    o
}

/// Solve one position optimally with the lower bounds `b`. `frame` is the physical frame of the
/// scrambled position (for the physical moves of the solutions). With `verbose`, each
/// iteration is reported on stderr.
pub fn solve<B: Bounds>(
    b: &B,
    cube: Cube,
    frame: Frame,
    limits: &Limits,
    verbose: bool,
) -> Outcome {
    // (compiled once per number of moves, so the face metric runs the same code as before)
    if b.moves().slice {
        solve_n::<B, MAX_MOVES>(b, cube, frame, limits, verbose)
    } else {
        solve_n::<B, N_MOVES>(b, cube, frame, limits, verbose)
    }
}

/// [`solve`] with `N` internal moves (`b.moves().n_moves()`).
fn solve_n<B: Bounds, const N: usize>(
    b: &B,
    cube: Cube,
    frame: Frame,
    limits: &Limits,
    verbose: bool,
) -> Outcome {
    let start = Instant::now();
    let p = Pos::of(&cube);
    let mut out = Outcome {
        internal: vec![],
        solutions: vec![],
        length: None,
        lower_bound: 0,
        stop: Stop::BadTables,
        nodes: 0,
        lookups: 0,
        prune_masks: vec![],
        children_by_left: vec![],
        elapsed: 0.0,
        iterations: vec![],
        symmetries: 1,
    };
    let Some(h0) = exact(b, &p) else {
        out.elapsed = start.elapsed().as_secs_f64();
        return out;
    };
    // the symmetries of the octahedron (other than the identity) that fix the position
    let root_syms: Vec<usize> = {
        let sym = super::sym::sym48();
        (1..super::sym::N_SYM)
            .filter(|&k| sym.conj_by_maps(&cube, k) == cube)
            .collect()
    };
    out.symmetries = root_syms.len() + 1;
    let sh = Shared::<B, N> {
        b,
        mt: opt_moves(),
        solved: Pos::solved(),
        cube,
        limits: limits.clone(),
        start,
        stop: AtomicBool::new(false),
        nodes: AtomicU64::new(0),
        lookups: AtomicU64::new(0),
        stats: Mutex::new(PruneStats::default()),
        prune_stats: std::env::var("FTO_OPT_PRUNE_STATS").is_ok_and(|v| v == "1"),
        sols: Mutex::new(vec![]),
        symmetries: if limits.max_solutions.is_none() && sym_start() {
            root_syms.clone()
        } else {
            vec![]
        },
    };
    let root = Task {
        path: [0; MAX_SPLIT],
        len: 0,
        p,
        h: h0,
    };
    let lo = usize::from(h0.iter().copied().max().unwrap_or(0));
    if verbose {
        let names: Vec<String> = (0..B::K)
            .map(|k| format!("{} {}", b.name(k), h0[k]))
            .collect();
        eprintln!(
            "root bounds: {} (search from depth {lo}); symmetries of the position: {} of 48{}",
            names.join(", "),
            root_syms.len() + 1,
            if sh.symmetries.is_empty() {
                ""
            } else {
                " (symmetric start on)"
            }
        );
    }
    let threads = rayon::current_num_threads();
    // ordering statistic of each task of the previous iteration, by its move prefix
    let mut stats: HashMap<u32, f32> = HashMap::new();
    let mut prev_nodes: Option<u64> = None;
    let mut depth = lo;
    out.lower_bound = lo;
    out.stop = loop {
        if limits.max_length.is_some_and(|l| depth > l) {
            break Stop::MaxLength;
        }
        let it_start = sh.nodes.load(Ordering::Relaxed);
        let it_lookups = sh.lookups.load(Ordering::Relaxed);
        let predicted = prev_nodes.map_or(0, |n| n.saturating_mul(13));
        let inline = threads == 1 || predicted < SERIAL_NODES;
        let want = if inline && threads > 1 {
            1
        } else {
            let per = (predicted / TASK_NODES) as usize;
            per.clamp(if threads == 1 { 64 } else { threads }, 64 * threads.max(4))
        };
        let mut ctr = Counter::default();
        let mut tasks = sh.tasks(&root, depth, want, &mut ctr);
        sh.nodes.fetch_add(ctr.nodes, Ordering::Relaxed);
        sh.lookups.fetch_add(ctr.lookups, Ordering::Relaxed);
        sh.stats
            .lock()
            .unwrap_or_else(PoisonError::into_inner)
            .add(&ctr.stats);
        // order by the previous iteration's statistic of the longest known prefix (stable:
        // ties keep the canonical order)
        if !stats.is_empty() {
            let score = |t: &Task| -> f32 {
                (0..=usize::from(t.len))
                    .rev()
                    .find_map(|l| stats.get(&key(&t.path[..l])).copied())
                    .unwrap_or(0.0)
            };
            let mut scored: Vec<(f32, Task)> = tasks.iter().map(|t| (score(t), *t)).collect();
            scored.sort_by(|a, b| b.0.total_cmp(&a.0));
            tasks = scored.into_iter().map(|(_, t)| t).collect();
        }
        let task_stats: Vec<AtomicU64> = (0..tasks.len()).map(|_| AtomicU64::new(0)).collect();
        let next = AtomicUsize::new(0);
        let worker = || {
            let mut path: Vec<u8> = Vec::with_capacity(64);
            let mut frames = vec![[Eval::<N>::new(), Eval::<N>::new()]; depth + 1];
            let mut ctr = Counter::default();
            loop {
                let i = next.fetch_add(1, Ordering::Relaxed);
                if i >= tasks.len() || sh.stop.load(Ordering::Relaxed) || sh.time_up() {
                    break;
                }
                let t = &tasks[i];
                path.clear();
                path.extend_from_slice(&t.path[..usize::from(t.len)]);
                let (before, lookups_before) = (ctr.nodes, ctr.lookups);
                ctr.stop = false;
                let mut stat = 0f32;
                sh.dfs(
                    &t.p,
                    t.h,
                    t.last(),
                    depth,
                    &mut path,
                    &mut ctr,
                    &mut frames,
                    &mut stat,
                );
                task_stats[i].store(u64::from(stat.to_bits()), Ordering::Relaxed);
                sh.nodes.fetch_add(ctr.nodes - before, Ordering::Relaxed);
                sh.lookups
                    .fetch_add(ctr.lookups - lookups_before, Ordering::Relaxed);
            }
            sh.stats
                .lock()
                .unwrap_or_else(PoisonError::into_inner)
                .add(&ctr.stats);
        };
        if inline {
            worker();
        } else {
            rayon::scope(|scope| {
                for _ in 0..threads.min(tasks.len()) {
                    scope.spawn(|_| worker());
                }
            });
        }
        let it_nodes = sh.nodes.load(Ordering::Relaxed) - it_start;
        prev_nodes = Some(it_nodes);
        stats = tasks
            .iter()
            .zip(&task_stats)
            .map(|(t, s)| {
                (
                    key(&t.path[..usize::from(t.len)]),
                    f32::from_bits(s.load(Ordering::Relaxed) as u32),
                )
            })
            .collect();
        let elapsed = start.elapsed().as_secs_f64();
        let it_lookups = sh.lookups.load(Ordering::Relaxed) - it_lookups;
        out.iterations.push(Iteration {
            depth,
            nodes: it_nodes,
            lookups: it_lookups,
            tasks: tasks.len(),
            elapsed,
        });
        let n_sols = sh.sols.lock().unwrap_or_else(PoisonError::into_inner).len();
        if verbose {
            eprintln!(
                "  depth {depth:>2}: {it_nodes:>14} nodes, {it_lookups:>14} lookups, {:>5} tasks, \
                 {elapsed:>9.3} s{}",
                tasks.len(),
                if n_sols > 0 {
                    format!(", {n_sols} solution(s)")
                } else {
                    String::new()
                }
            );
        }
        if n_sols > 0 {
            out.length = Some(depth);
            out.lower_bound = depth;
            break Stop::Solved;
        }
        if sh.stop.load(Ordering::Relaxed) || sh.time_up() {
            // (the interrupted iteration proves nothing)
            break Stop::Time;
        }
        // no solution of at most `depth` moves
        depth += 1;
        out.lower_bound = depth;
    };
    out.internal = sh.sols.into_inner().unwrap_or_else(PoisonError::into_inner);
    out.solutions = out
        .internal
        .iter()
        .map(|s| solution_to_physical(frame, s, b.moves().wide))
        .collect();
    out.nodes = sh.nodes.load(Ordering::Relaxed);
    out.lookups = sh.lookups.load(Ordering::Relaxed);
    let stats = sh
        .stats
        .into_inner()
        .unwrap_or_else(PoisonError::into_inner);
    out.prune_masks = stats.masks;
    out.prune_masks.truncate(1 << B::K);
    out.children_by_left = stats.left;
    out.elapsed = start.elapsed().as_secs_f64();
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn move_masks_match_rules() {
        for slice in [false, true] {
            let n = if slice { MAX_MOVES } else { N_MOVES };
            for last in (0..n).chain([MAX_MOVES]) {
                for m in 0..MAX_MOVES {
                    let l = if last == MAX_MOVES { usize::MAX } else { last };
                    let want = m < n && allowed_after(slice, l, m);
                    assert_eq!(ALLOWED[usize::from(slice)][last] >> m & 1 == 1, want);
                }
            }
        }
    }

    #[test]
    fn prefix_keys_are_distinct() {
        let mut seen = std::collections::HashSet::new();
        for len in 0..=MAX_SPLIT {
            let mut path = [0u8; MAX_SPLIT];
            for x in 0..MAX_MOVES.pow(len as u32) {
                let mut y = x;
                for p in path.iter_mut().take(len) {
                    *p = (y % MAX_MOVES) as u8;
                    y /= MAX_MOVES;
                }
                assert!(seen.insert(key(&path[..len])));
            }
        }
    }
}
