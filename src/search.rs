//! Two-phase search.
//!
//! Phase 1: IDA* over all 16 moves until the position is in `H_A` = <U, F, BR, BL>, with the
//! last move outside `H_A`. Phase 2: the exact table gives the optimal A-move completion.
//! Every solution splits uniquely into (longest prefix ending in a non-A move) + (A-moves),
//! so exhausting phase-1 depths up to `best - 1` proves optimality.
//!
//! Two variants run side by side: targeting class A, and targeting class B (by conjugating the
//! position with a 90° rotation that swaps the classes).

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
    /// phase-1 leaves reached: paths that pass every bound at their iteration's depth
    pub leaves: u64,
    pub elapsed: f64,
    /// true if the search proved that no shorter solution exists
    pub optimal: bool,
    /// improvement history
    pub events: Vec<Event>,
}

struct Variant {
    cube: Cube,
    sym: Mat,
    root: Fc,
    h_ae: u8,
    h_cb: u8,
    /// The root conjugated by `SWAP_AB` (its `H_A` distance is this root's `H_B` distance), and
    /// its bounds: a second lower bound on the distance to solved.
    root_b: Fc,
    hb_ae: u8,
    hb_cb: u8,
}

struct Shared<'a, I: AeIndex> {
    t: &'a Tables<I>,
    variants: Vec<Variant>,
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
    /// leaves that turned out not to be in `H_A` (rare: the tables certify the rest)
    rejected: AtomicU64,
    stop: AtomicBool,
    best: Mutex<(Option<Vec<PhysMove>>, Vec<Event>)>,
    mb: [usize; N_MOVES],
    /// largest CB distance (8): with at least this many moves left CB cannot prune
    cb_max: u8,
}

#[derive(Clone)]
struct Task {
    v: usize,
    path: Vec<usize>,
    x: Fc,
    h_ae: u8,
    h_cb: u8,
    xb: Fc,
    hb_ae: u8,
    hb_cb: u8,
}

/// Internal move on the SWAP_AB-conjugated position matching each move.
fn swap_moves() -> [usize; N_MOVES] {
    std::array::from_fn(|m| conjugate_move(m, &SWAP_AB))
}

const CHECK_EVERY: u64 = 1 << 14;

/// The AE table is checked first (before CB) with at least this many moves left after the
/// children: the CB table (max 8) prunes almost nothing there (+7% / +10% nodes/s at 48T / 1T
/// when added with 7+ moves, opt/ae-first; lowered to 6 in the opt/stages stack, measured
/// only as part of it: +10% / +14% nodes/s).
const AE_FIRST: usize = 6;

/// No deferred CB lookup (see `dfs`).
const NO_CB: usize = usize::MAX;

/// Scratch arrays of one `dfs` level: the children of a node and their indices and bounds.
/// They live in a per-worker stack of frames, three per depth (see [`Shared::dfs`]), rather
/// than on the call stack, so that they are not re-initialised at every node (~3 KB of stores
/// per node; +1.3% nodes/s at 1T in optimal mode, when added).
#[derive(Clone)]
struct Scratch {
    ys: [Fc; N_MOVES],
    /// (the children conjugated by `SWAP_AB`, only with the class-B bound)
    ybs: [Fc; N_MOVES],
    ms: [u8; N_MOVES],
    has: [u8; N_MOVES],
    hcs: [u8; N_MOVES],
    hbas: [u8; N_MOVES],
    hbcs: [u8; N_MOVES],
    idx: [usize; N_MOVES],
    cbi: [usize; N_MOVES],
    pend: [usize; N_MOVES],
    bidx: [usize; N_MOVES],
    bcbi: [usize; N_MOVES],
    ci: [[usize; 4]; N_MOVES],
    /// the children still alive after each stage (indices into the arrays above). Stages
    /// compact this list without branching on the prune decision, instead of copying
    /// positions (neutral when added: -1.5% nodes/s at 1T, within noise).
    sel: [u8; N_MOVES],
    /// number of children
    n: usize,
    /// number alive after `stage2` (in `sel`)
    k: usize,
    /// the node's own exact distances and the smallest AE distance of a child after
    /// `stage2`, and whether `stage2` issued the joint-table lookups
    h_ae: u8,
    h_cb: u8,
    min_ha: u8,
    cbe_done: bool,
}

impl Scratch {
    fn new() -> Scratch {
        Scratch {
            ys: [Fc::ZERO; N_MOVES],
            ybs: [Fc::ZERO; N_MOVES],
            ms: [0; N_MOVES],
            has: [0; N_MOVES],
            hcs: [0; N_MOVES],
            hbas: [0; N_MOVES],
            hbcs: [0; N_MOVES],
            idx: [0; N_MOVES],
            cbi: [0; N_MOVES],
            pend: [NO_CB; N_MOVES],
            bidx: [0; N_MOVES],
            bcbi: [0; N_MOVES],
            ci: [[0; 4]; N_MOVES],
            sel: [0; N_MOVES],
            n: 0,
            k: 0,
            h_ae: 0,
            h_cb: 0,
            min_ha: 0,
            cbe_done: false,
        }
    }
}

/// Moves allowed after each move (index `N_MOVES`: no previous move), as bit masks.
static ALLOWED: [u16; N_MOVES + 1] = {
    let mut a = [0u16; N_MOVES + 1];
    let mut last = 0;
    while last <= N_MOVES {
        let mut m = 0;
        while m < N_MOVES {
            // same rule as `allowed_after`
            let ok = if last == N_MOVES {
                true
            } else {
                let (ls, s) = (last >> 1, m >> 1);
                !(s == ls || ((s & 3) == (ls & 3) && s < ls))
            };
            if ok {
                a[last] |= 1 << m;
            }
            m += 1;
        }
        last += 1;
    }
    a
};

/// The class-A moves as a bit mask.
static A_MOVES: u16 = {
    let mut a = 0u16;
    let mut m = 0;
    while m < N_MOVES {
        if matches!(m >> 1, 1 | 2 | 4 | 7) {
            a |= 1 << m;
        }
        m += 1;
    }
    a
};

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

impl<I: AeIndex> Shared<'_, I> {
    fn time_up(&self) -> bool {
        if let Some(t) = self.limits.time {
            if self.start.elapsed() >= t {
                self.stop.store(true, Ordering::Relaxed);
                return true;
            }
        }
        false
    }

    fn leaf(&self, v: usize, path: &[usize]) {
        let var = &self.variants[v];
        let cube = var.cube.apply_seq(path);
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
                    // map the variant's internal moves back to the original frame
                    let inv = transpose(&var.sym);
                    let internal: Vec<usize> = path
                        .iter()
                        .chain(p2.iter())
                        .map(|&m| conjugate_move(m, &inv))
                        .collect();
                    let phys = solution_to_physical(self.frame, &internal);
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
        } else {
            self.rejected.fetch_add(1, Ordering::Relaxed);
        }
    }

    /// Depth-first search below a node at depth `g` with bound `bound`.
    /// With `B`, children are also pruned when their class-B bound shows they cannot lead to
    /// a solution shorter than the best one found so far.
    /// If `cb_pend` is not `NO_CB`, `h_cb` is the parent's CB distance and this node's CB
    /// entry (already prefetched by the parent) is at `cb_pend`.
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
        v: usize,
        path: &mut Vec<usize>,
        x: &Fc,
        h_ae: u8,
        h_cb: u8,
        cb_pend: usize,
        xb: &Fc,
        hb_ae: u8,
        hb_cb: u8,
        bound: usize,
        nodes: &mut u64,
        stat: &mut f32,
        frames: &mut [[Scratch; 3]],
    ) {
        let g = path.len();
        if g == bound {
            if h_ae == 0 && h_cb == 0 && (g == 0 || !is_a_move(path[g - 1])) {
                self.leaf(v, path);
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
        self.expand(x, last, budget, fr, nodes);
        self.stage2::<B>(h_ae, h_cb, cb_pend, xb, budget, fr);
        self.finish::<B>(v, path, xb, hb_ae, hb_cb, bound, nodes, stat, fr, frames);
    }

    /// First step of the evaluation of a node with `budget` moves left after its children
    /// (`last`: the move leading to it): generate the children into `fr`, count them, and
    /// issue their first table lookups (AE with 6+ moves left, CB otherwise). Every stage
    /// issues the lookups of all children before reading any of them, so the cache misses of
    /// different children overlap (2.0x nodes/s at 48T and 1T when added, opt/prefetch).
    #[inline(always)]
    fn expand(&self, x: &Fc, last: usize, budget: usize, fr: &mut Scratch, nodes: &mut u64) {
        let t = self.t;
        let f = fast();
        let fw = f.flip_word(&masks(&x.e));
        let mut n = 0;
        let mut moves = ALLOWED[last.min(N_MOVES)];
        if budget == 0 {
            moves &= !A_MOVES;
        }
        while moves != 0 {
            let m = moves.trailing_zeros() as usize;
            moves &= moves - 1;
            fr.ys[n & 15] = f.apply_w(x, fw, m);
            fr.ms[n & 15] = m as u8;
            n += 1;
        }
        fr.n = n;
        // count every evaluated node; check the clock regularly (`finish` sees the stop)
        let before = *nodes;
        *nodes += n as u64;
        if before / CHECK_EVERY != *nodes / CHECK_EVERY {
            self.time_up();
        }
        if budget >= AE_FIRST {
            t.ae_ix
                .index_batch_fast_pf(f, &fr.ys[..n], &mut fr.idx[..n], &t.ae.data);
        } else {
            for i in 0..n {
                fr.cbi[i] = CbTable::index_fast(f, &fr.ys[i]);
                prefetch_at(&t.cb.tab.data, fr.cbi[i] >> 5);
            }
        }
    }

    /// Second step (see [`Self::dfs`]): read the first lookups of the children in `fr`, and
    /// issue the lookups of the next stages for those still alive.
    #[inline(always)]
    fn stage2<const B: bool>(
        &self,
        h_ae: u8,
        h_cb: u8,
        cb_pend: usize,
        xb: &Fc,
        budget: usize,
        fr: &mut Scratch,
    ) {
        let t = self.t;
        let f = fast();
        let Scratch {
            ys,
            ybs,
            ms,
            has,
            hcs,
            idx,
            cbi,
            bidx,
            bcbi,
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
            let mut min_ha = u8::MAX;
            for i in 0..n {
                let ha = decode(h_ae, t.ae.get(idx[i]));
                let ok = ha as usize <= budget;
                has[i] = ha;
                min_ha = min_ha.min(if ok { ha } else { u8::MAX });
                sel[k & 15] = i as u8;
                k += usize::from(ok);
            }
            *fmin_ha = min_ha;
            for &i in &sel[..k] {
                let i = usize::from(i & 15);
                cbi[i] = CbTable::index_fast(f, &ys[i]);
                prefetch_at(&t.cb.tab.data, cbi[i] >> 5);
            }
            if k > 0 && cb_pend != NO_CB {
                // this node's own CB entry was prefetched by the parent, before its AE stage
                h_cb = decode(h_cb, t.cb.tab.get(cb_pend));
            }
            // The class-B lookups are issued here too, so that their latency overlaps the
            // CB and joint-table stages (+6% / +4% nodes/s at 48T / 1T, optimal 0.95x time,
            // when added, opt/bpath).
            if B && k > 0 {
                let fwb = f.flip_word(&masks(&xb.e));
                for &i in &sel[..k] {
                    let i = usize::from(i & 15);
                    ybs[i] = f.apply_w(xb, fwb, self.mb[usize::from(ms[i] & 15)]);
                }
                t.ae_ix.index_sel_pf(f, ybs, &sel[..k], bidx, &t.ae.data);
                for &i in &sel[..k] {
                    let i = usize::from(i & 15);
                    bcbi[i] = CbTable::index_fast(f, &ybs[i]);
                    prefetch_at(&t.cb.tab.data, bcbi[i] >> 5);
                }
            }
            // The joint table goes before reading the CB words, so their latency overlaps it
            // (part of the opt/stages stack: +10% / +14% nodes/s at 48T / 1T as a whole).
            if budget <= 8 {
                if let Some(cbe) = &t.cbe {
                    for &i in &sel[..k] {
                        let i = usize::from(i & 15);
                        cbe.prefetch_all(f, &ys[i], &mut ci[i]);
                    }
                    *cbe_done = true;
                }
            }
        } else {
            for i in 0..n {
                let hc = decode(h_cb, t.cb.tab.get(cbi[i]));
                hcs[i] = hc;
                sel[k & 15] = i as u8;
                k += usize::from(hc as usize <= budget);
            }
            t.ae_ix.index_sel_pf(f, ys, &sel[..k], idx, &t.ae.data);
        }
        *kk = k;
        *fh_ae = h_ae;
        *fh_cb = h_cb;
    }

    /// Last step of the evaluation of a node at depth `g = path.len()` (see [`Self::dfs`]),
    /// after [`Self::expand`] and [`Self::stage2`] on `fr`, and the search below it.
    fn finish<const B: bool>(
        &self,
        v: usize,
        path: &mut Vec<usize>,
        xb: &Fc,
        hb_ae: u8,
        hb_cb: u8,
        bound: usize,
        nodes: &mut u64,
        stat: &mut f32,
        fr: &mut Scratch,
        frames: &mut [[Scratch; 3]],
    ) {
        if self.stop.load(Ordering::Relaxed) {
            return;
        }
        let g = path.len();
        let budget = bound - g - 1;
        let t = self.t;
        let f = fast();
        let Scratch {
            ys,
            ybs,
            ms,
            has,
            hcs,
            hbas,
            hbcs,
            idx,
            cbi,
            pend,
            bidx,
            bcbi,
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
                        let ok = cbe.bound(&ci[usize::from(i & 15)]) as usize <= budget;
                        sel[k2 & 15] = i;
                        k2 += usize::from(ok);
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
            // opt/micro-compute stack: 1.29x / 1.28x nodes/s at 48T / 1T as a whole).
            let defer = !cbe_done && budget >= self.cb_max as usize && min_ha >= self.cb_max;
            for j in 0..k {
                let i = sel[j];
                let u = usize::from(i & 15);
                // (a deferred child carries this node's distance, at most cb_max <= budget)
                let hc = if defer {
                    h_cb
                } else {
                    decode(h_cb, t.cb.tab.get(cbi[u]))
                };
                let mut ok = hc as usize <= budget;
                pend[u] = if defer { cbi[u] } else { NO_CB };
                if B {
                    let hba = decode(hb_ae, t.ae.get(bidx[u]));
                    let hbc = decode(hb_cb, t.cb.tab.get(bcbi[u]));
                    ok &= hba.max(hbc) as usize <= lim;
                    hbas[u] = hba;
                    hbcs[u] = hbc;
                }
                hcs[u] = hc;
                sel[n & 15] = i;
                n += usize::from(ok);
            }
            b_done = true;
        } else {
            for j in 0..k {
                let i = sel[j];
                let u = usize::from(i & 15);
                let ha = decode(h_ae, t.ae.get(idx[u]));
                has[u] = ha;
                pend[u] = NO_CB;
                sel[n & 15] = i;
                n += usize::from(ha as usize <= budget);
            }
        }
        // The joint corner/B-centre/edge-triple table only adds information near the leaves.
        if budget <= 8 && !cbe_done {
            if let Some(cbe) = &t.cbe {
                for &i in &sel[..n] {
                    let i = usize::from(i & 15);
                    cbe.prefetch_all(f, &ys[i], &mut ci[i]);
                }
                let mut k = 0;
                for j in 0..n {
                    let i = sel[j];
                    let ok = cbe.bound(&ci[usize::from(i & 15)]) as usize <= budget;
                    sel[k & 15] = i;
                    k += usize::from(ok);
                }
                n = k;
            }
        }
        // Class-B bound against the best solution so far (optimal mode 0.29x / 0.67x time at
        // 48T / 1T when added, opt/classb).
        if B && n > 0 && !b_done {
            let best = self.prune_len.load(Ordering::Relaxed);
            // a child at depth g+1 with bound h cannot give a solution shorter than g+1+h
            let lim = best.saturating_sub(g + 2);
            let fwb = f.flip_word(&masks(&xb.e));
            for &i in &sel[..n] {
                let i = usize::from(i & 15);
                ybs[i] = f.apply_w(xb, fwb, self.mb[usize::from(ms[i] & 15)]);
            }
            t.ae_ix.index_sel_pf(f, ybs, &sel[..n], idx, &t.ae.data);
            for &i in &sel[..n] {
                let i = usize::from(i & 15);
                cbi[i] = CbTable::index_fast(f, &ybs[i]);
                prefetch_at(&t.cb.tab.data, cbi[i] >> 5);
            }
            let mut k = 0;
            for j in 0..n {
                let i = sel[j];
                let u = usize::from(i & 15);
                let hba = decode(hb_ae, t.ae.get(idx[u]));
                let hbc = decode(hb_cb, t.cb.tab.get(cbi[u]));
                hbas[u] = hba;
                hbcs[u] = hbc;
                sel[k & 15] = i;
                k += usize::from(hba.max(hbc) as usize <= lim);
            }
            n = k;
        }
        // Ordering statistic of the task: how close its subtree gets to the phase-1 target.
        for &i in &sel[..n] {
            let i = usize::from(i & 15);
            *stat += STAT_W[has[i].max(hcs[i]) as usize];
        }
        if budget == 0 {
            // the children are leaves
            for &i in &sel[..n] {
                let i = usize::from(i & 15);
                let m = ms[i] as usize;
                if has[i] == 0 && hcs[i] == 0 && !is_a_move(m) {
                    path.push(m);
                    self.leaf(v, path);
                    path.pop();
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
        let sel = *sel;
        let child = |j: usize| usize::from(sel[j & 15] & 15);
        let yb_of = |i: usize| if B { &ybs[i] } else { xb };
        if n > 0 {
            let i = child(0);
            self.expand(&ys[i], ms[i] as usize, cb, f0, nodes);
        }
        if n > 1 {
            let i = child(1);
            self.expand(&ys[i], ms[i] as usize, cb, f1, nodes);
        }
        if n > 0 {
            let i = child(0);
            self.stage2::<B>(has[i], hcs[i], pend[i], yb_of(i), cb, f0);
        }
        for j in 0..n {
            if j + 2 < n {
                let i = child(j + 2);
                self.expand(&ys[i], ms[i] as usize, cb, f2, nodes);
            }
            if j + 1 < n {
                let i = child(j + 1);
                self.stage2::<B>(has[i], hcs[i], pend[i], yb_of(i), cb, f1);
            }
            let i = child(j);
            path.push(ms[i] as usize);
            self.finish::<B>(
                v,
                path,
                yb_of(i),
                hbas[i],
                hbcs[i],
                bound,
                nodes,
                stat,
                f0,
                frames,
            );
            path.pop();
            if self.stop.load(Ordering::Relaxed) {
                return;
            }
            (f0, f1, f2) = (f1, f2, f0);
        }
    }

    /// All children of a task (unpruned), with their exact bounds on both classes. The
    /// lookups of all children are issued before any is read.
    fn expand_all(&self, t: &Task) -> Vec<Task> {
        let tb = self.t;
        let f = fast();
        let last = t.path.last().copied().unwrap_or(usize::MAX);
        let em = masks(&t.x.e);
        let emb = masks(&t.xb.e);
        let mut ys = [t.x; N_MOVES];
        let mut yb = [t.xb; N_MOVES];
        let mut ms = [0usize; N_MOVES];
        let mut n = 0;
        for m in 0..N_MOVES {
            if !allowed_after(last, m) {
                continue;
            }
            ys[n] = f.apply(&t.x, &em, m);
            yb[n] = f.apply(&t.xb, &emb, self.mb[m]);
            ms[n] = m;
            n += 1;
        }
        let mut ia = [0usize; N_MOVES];
        let mut iab = [0usize; N_MOVES];
        let mut ic = [0usize; N_MOVES];
        let mut icb = [0usize; N_MOVES];
        tb.ae_ix
            .index_batch_fast_pf(f, &ys[..n], &mut ia[..n], &tb.ae.data);
        tb.ae_ix
            .index_batch_fast_pf(f, &yb[..n], &mut iab[..n], &tb.ae.data);
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
                    v: t.v,
                    path,
                    x: ys[i],
                    h_ae: decode(t.h_ae, tb.ae.get(ia[i])),
                    h_cb: decode(t.h_cb, tb.cb.tab.get(ic[i])),
                    xb: yb[i],
                    hb_ae: decode(t.hb_ae, tb.ae.get(iab[i])),
                    hb_cb: decode(t.hb_cb, tb.cb.tab.get(icb[i])),
                }
            })
            .collect()
    }

    /// Tasks for one iteration (as node ids of `tree`): the roots, expanded level by level
    /// until there are enough tasks for the thread pool. The depth stays small because most
    /// of an IDA* iteration's work lies in the middle levels, which must not be explored
    /// serially. The tree keeps every node evaluated by earlier iterations with all its
    /// children, so a later iteration only evaluates nodes it has not seen and re-filters.
    /// (Splitting: 48T first solution 0.50x when added, opt/split; memoized tree and adaptive
    /// task count: a further 0.51x, opt/split2 stack.)
    fn tasks(&self, tree: &mut Tree, bound: usize, threads: usize, want: usize) -> Vec<u32> {
        if tree.nodes.is_empty() {
            for (v, var) in self.variants.iter().enumerate() {
                tree.nodes.push(Task {
                    v,
                    path: vec![],
                    x: var.root,
                    h_ae: var.h_ae,
                    h_cb: var.h_cb,
                    xb: var.root_b,
                    hb_ae: var.hb_ae,
                    hb_cb: var.hb_cb,
                });
                tree.kids.push(None);
            }
        }
        let alive = |t: &Task| t.h_ae.max(t.h_cb) as usize + t.path.len() <= bound;
        let mut front: Vec<u32> = (0..self.variants.len() as u32)
            .filter(|&i| alive(&tree.nodes[i as usize]))
            .collect();
        let max_k = bound.saturating_sub(4).min(MAX_SPLIT);
        let mut k = 0;
        // (single-threaded runs are split too, so per-task decisions can change mid-iteration)
        while front.len() < want && k < max_k {
            let need: Vec<u32> = front
                .iter()
                .copied()
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
            let mut next = Vec::with_capacity(front.len() * 13);
            for &i in &front {
                // (every front node has been expanded above)
                if let Some((lo, hi)) = tree.kids[i as usize] {
                    next.extend((lo..hi).filter(|&c| alive(&tree.nodes[c as usize])));
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

/// Solve one position. `frame` is the physical frame of the scrambled position.
pub fn solve<I: AeIndex>(
    t: &Tables<I>,
    cube: Cube,
    frame: Frame,
    limits: &Limits,
    report: Report,
) -> Outcome {
    let start = Instant::now();
    let variants: Option<Vec<Variant>> = [IDENTITY, SWAP_AB]
        .iter()
        .map(|sym| {
            let c = if *sym == IDENTITY {
                cube
            } else {
                cube.conjugate(sym)?
            };
            let root = fast().of(&c);
            let (h_ae, h_cb) = t.exact_fast(&root);
            let root_b = fast().of(&c.conjugate(&SWAP_AB)?);
            let (hb_ae, hb_cb) = t.exact_fast(&root_b);
            Some(Variant {
                cube: c,
                sym: *sym,
                root,
                h_ae,
                h_cb,
                root_b,
                hb_ae,
                hb_cb,
            })
        })
        .collect();
    let Some(variants) = variants else {
        // not a valid position in the fixed frame: nothing to search
        return Outcome {
            solution: None,
            length: None,
            solutions_found: 0,
            nodes: 0,
            leaves: 0,
            elapsed: start.elapsed().as_secs_f64(),
            optimal: false,
            events: vec![],
        };
    };
    // Every solution is a prefix followed by class-A moves, so one variant would suffice to
    // prove optimality; both are kept because the earlier a short solution is found, the more
    // the class-B bound prunes (measured: one variant was 1.9-3.9x slower in optimal mode),
    // and because with both the proof can stop one iteration early (see the main loop).
    let sh = Shared {
        t,
        variants,
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
        rejected: AtomicU64::new(0),
        stop: AtomicBool::new(false),
        best: Mutex::new((None, vec![])),
        mb: swap_moves(),
        cb_max: (t.cb.tab.hist.len() - 1) as u8,
    };
    let lo = sh
        .variants
        .iter()
        .map(|v| v.h_ae.max(v.h_cb) as usize)
        .fold(usize::MAX, usize::min);
    let mut optimal = false;
    let threads = rayon::current_num_threads();
    let mut depth = lo;
    let mut tree = Tree::default();
    // ordering statistic of each tree node from the previous iteration (0 if not a task then)
    let mut node_stat: Vec<f32> = vec![];
    let mut prev_ids: Vec<u32> = vec![];
    let mut prev_nodes: Option<u64> = None;
    // Target-length passes: only in optimal mode, when the best is at least this far above the
    // lower bound (0 = never). On by default for one thread only: there the best stays loose
    // for long and the passes cut the proof time (1T 17-move set 11.95 -> 5.27 s, 0.44x, when
    // added); with 48 threads a tight best is found early and the rerun iterations cost more
    // (19-move set 1.07x, random positions 1.00-1.04x). FTO_TARGET_GAP overrides. Off with
    // FTO_FULL_PROOF: the passes rest on the same argument as the early proof end.
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
        // the class-B variant already found it at a smaller depth, in an iteration completed
        // before this one. (This needs both variants in every iteration.) When solutions are
        // counted, that iteration still yields completions (longer ones), so it is run.
        let proof_depth = depth + usize::from(limits.max_solutions.is_none() && !full_proof());
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
        if let Some(t) = target {
            if depth >= t {
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
            let lower = depth + 1;
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
                .store(if full_proof() { 0 } else { depth + 1 }, Ordering::Relaxed);
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
            ids.sort_by(|a, b| node_stat[*b as usize].total_cmp(&node_stat[*a as usize]));
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
        let var_nodes = [AtomicU64::new(0), AtomicU64::new(0)];
        let tree_ref = &tree;
        let worker = || {
            let mut path: Vec<usize> = Vec::with_capacity(64);
            // one frame per level below a task (tasks are at depth >= 0)
            let mut frames = vec![[Scratch::new(), Scratch::new(), Scratch::new()]; depth + 1];
            // the per-thread counter persists across tasks so time checks happen regularly
            let mut nodes = 0u64;
            loop {
                let i = next.fetch_add(1, Ordering::Relaxed);
                if i >= ids.len() || sh.stop.load(Ordering::Relaxed) || sh.time_up() {
                    break;
                }
                let task = &tree_ref.nodes[ids[i] as usize];
                path.clear();
                path.extend_from_slice(&task.path);
                let before = nodes;
                let mut stat = 0f32;
                let use_b = !count_mode && sh.prune_len.load(Ordering::Relaxed) < usize::MAX;
                if use_b {
                    sh.dfs::<true>(
                        task.v,
                        &mut path,
                        &task.x,
                        task.h_ae,
                        task.h_cb,
                        NO_CB,
                        &task.xb,
                        task.hb_ae,
                        task.hb_cb,
                        depth,
                        &mut nodes,
                        &mut stat,
                        &mut frames,
                    );
                } else {
                    sh.dfs::<false>(
                        task.v,
                        &mut path,
                        &task.x,
                        task.h_ae,
                        task.h_cb,
                        NO_CB,
                        &task.xb,
                        task.hb_ae,
                        task.hb_cb,
                        depth,
                        &mut nodes,
                        &mut stat,
                        &mut frames,
                    );
                }
                stats[i].store(stat.to_bits(), Ordering::Relaxed);
                sh.nodes.fetch_add(nodes - before, Ordering::Relaxed);
                var_nodes[task.v % 2].fetch_add(nodes - before, Ordering::Relaxed);
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
        let var_counts = var_nodes.map(AtomicU64::into_inner);
        for &i in &prev_ids {
            node_stat[i as usize] = 0.0;
        }
        for (i, &id) in ids.iter().enumerate() {
            node_stat[id as usize] = f32::from_bits(stats[i].load(Ordering::Relaxed));
        }
        let n_tasks = ids.len();
        prev_ids = ids;
        if std::env::var_os("FTO_DEBUG").is_some() {
            eprintln!("  depth {depth}: split done at {t_split:.4}s");
            eprintln!(
                "  depth {depth}: {} tasks, {} nodes ({} / {} by variant), {} solutions, best {:?}, {:.3}s",
                n_tasks,
                sh.nodes.load(Ordering::Relaxed),
                var_counts[0],
                var_counts[1],
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
    Outcome {
        solution: sol,
        length,
        solutions_found: sh.count.load(Ordering::Relaxed),
        nodes: sh.nodes.load(Ordering::Relaxed),
        leaves: sh.count.load(Ordering::Relaxed) + sh.rejected.load(Ordering::Relaxed),
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
        for last in 0..=N_MOVES {
            for m in 0..N_MOVES {
                let l = if last == N_MOVES { usize::MAX } else { last };
                assert_eq!(ALLOWED[last] >> m & 1 == 1, allowed_after(l, m));
            }
        }
        for m in 0..N_MOVES {
            assert_eq!(A_MOVES >> m & 1 == 1, is_a_move(m));
        }
    }
}
