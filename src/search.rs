//! Two-phase search.
//!
//! Phase 1: IDA* over all 16 moves until the position is in H_A = <U, F, BR, BL>, with the
//! last move outside H_A. Phase 2: the exact table gives the optimal A-move completion.
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
use std::sync::Mutex;
use std::time::{Duration, Instant};

#[derive(Clone, Debug, Default)]
pub struct Limits {
    /// stop after this many phase-2 completions (solutions of any length)
    pub max_solutions: Option<u64>,
    /// stop after this much time
    pub time: Option<Duration>,
    /// stop once a solution of at most this length is found
    pub max_length: Option<usize>,
}

impl Limits {
    pub fn is_optimal(&self) -> bool {
        self.max_solutions.is_none() && self.time.is_none() && self.max_length.is_none()
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
    /// The root conjugated by SWAP_AB (its H_A distance is this root's H_B distance), and
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
    count: AtomicU64,
    nodes: AtomicU64,
    stop: AtomicBool,
    best: Mutex<(Option<Vec<PhysMove>>, Vec<Event>)>,
    mb: [usize; N_MOVES],
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

impl<'a, I: AeIndex> Shared<'a, I> {
    fn time_up(&self) -> bool {
        if let Some(t) = self.limits.time {
            if self.start.elapsed() >= t {
                self.stop.store(true, Ordering::Relaxed);
                return true;
            }
        }
        false
    }

    #[inline(always)]
    /// A child that survives the phase-1 bound, with its exact bounds on both classes.
    fn child(&self, t: &Task, m: usize, budget: usize) -> Option<Task> {
        let tb = self.t;
        let f = fast();
        let y = f.apply(&t.x, &masks(&t.x.e), m);
        let hc = decode(t.h_cb, tb.cb.tab.get(CbTable::index_fast(f, &y)));
        if hc as usize > budget {
            return None;
        }
        let ha = decode(t.h_ae, tb.ae.get(tb.ae_ix.index_fast(f, &y, f.a_rank(&y))));
        if ha as usize > budget {
            return None;
        }
        let yb = f.apply(&t.xb, &masks(&t.xb.e), self.mb[m]);
        let hbc = decode(t.hb_cb, tb.cb.tab.get(CbTable::index_fast(f, &yb)));
        let hba = decode(
            t.hb_ae,
            tb.ae.get(tb.ae_ix.index_fast(f, &yb, f.a_rank(&yb))),
        );
        let mut path = t.path.clone();
        path.push(m);
        Some(Task {
            v: t.v,
            path,
            x: y,
            h_ae: ha,
            h_cb: hc,
            xb: yb,
            hb_ae: hba,
            hb_cb: hbc,
        })
    }

    fn leaf(&self, v: usize, path: &[usize]) {
        let var = &self.variants[v];
        let cube = var.cube.apply_seq(path);
        let d1 = path.len();
        // with the small tier the edge chirality is not in the tables: verify membership
        if let Some(d2) = self.t.p2.length(&cube) {
            let n = self.count.fetch_add(1, Ordering::Relaxed) + 1;
            if self.limits.max_solutions.is_some_and(|k| n >= k) {
                self.stop.store(true, Ordering::Relaxed);
            }
            let total = d1 + d2 as usize;
            if total < self.best_len.load(Ordering::Relaxed) {
                let mut g = self.best.lock().unwrap();
                if total < self.best_len.load(Ordering::Relaxed) {
                    let p2 = self.t.p2.solve(&cube).expect("phase 2 failed");
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
                    if self.limits.max_length.is_some_and(|l| total <= l) {
                        self.stop.store(true, Ordering::Relaxed);
                    }
                }
            }
        }
    }

    /// Depth-first search below a node at depth `g` with bound `bound`.
    /// With `B`, children are also pruned when their class-B bound shows they cannot lead to
    /// a solution shorter than the best one found so far.
    fn dfs<const B: bool>(
        &self,
        v: usize,
        path: &mut Vec<usize>,
        x: &Fc,
        h_ae: u8,
        h_cb: u8,
        xb: &Fc,
        hb_ae: u8,
        hb_cb: u8,
        bound: usize,
        nodes: &mut u64,
        stat: &mut f32,
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
        let last = path.last().copied().unwrap_or(usize::MAX);
        let budget = bound - g - 1;
        let t = self.t;
        let f = fast();
        let em = masks(&x.e);
        // Evaluate all children as a batch: every stage issues prefetches for all children
        // before reading any of them, so the cache misses of different children overlap.
        let mut ys = [*x; N_MOVES];
        let mut ms = [0u8; N_MOVES];
        let mut n = 0;
        for m in 0..N_MOVES {
            if !allowed_after(last, m) || (budget == 0 && is_a_move(m)) {
                continue;
            }
            ys[n] = f.apply(x, &em, m);
            ms[n] = m as u8;
            n += 1;
        }
        // count every evaluated node; check the clock regularly
        let before = *nodes;
        *nodes += n as u64;
        if before / CHECK_EVERY != *nodes / CHECK_EVERY && self.time_up() {
            return;
        }
        let mut has = [0u8; N_MOVES];
        let mut hcs = [0u8; N_MOVES];
        let mut idx = [0usize; N_MOVES];
        let mut cbi = [0usize; N_MOVES];
        // The CB table (max 8) prunes almost nothing with 7+ moves left: check AE first there.
        let ae_first = budget >= 7;
        if ae_first {
            t.ae_ix.index_batch_fast_pf(f, &ys[..n], &mut idx[..n], &t.ae.data);
            let mut k = 0;
            for i in 0..n {
                let ha = decode(h_ae, t.ae.get(idx[i]));
                if ha as usize <= budget {
                    ys[k] = ys[i];
                    ms[k] = ms[i];
                    has[k] = ha;
                    cbi[k] = CbTable::index_fast(f, &ys[k]);
                    prefetch(&t.cb.tab.data[cbi[k] >> 5]);
                    k += 1;
                }
            }
            n = 0;
            for i in 0..k {
                let hc = decode(h_cb, t.cb.tab.get(cbi[i]));
                if hc as usize <= budget {
                    ys[n] = ys[i];
                    ms[n] = ms[i];
                    has[n] = has[i];
                    hcs[n] = hc;
                    n += 1;
                }
            }
        } else {
            for i in 0..n {
                cbi[i] = CbTable::index_fast(f, &ys[i]);
                prefetch(&t.cb.tab.data[cbi[i] >> 5]);
            }
            let mut k = 0;
            for i in 0..n {
                let hc = decode(h_cb, t.cb.tab.get(cbi[i]));
                if hc as usize <= budget {
                    ys[k] = ys[i];
                    ms[k] = ms[i];
                    hcs[k] = hc;
                    k += 1;
                }
            }
            t.ae_ix.index_batch_fast_pf(f, &ys[..k], &mut idx[..k], &t.ae.data);
            n = 0;
            for i in 0..k {
                let ha = decode(h_ae, t.ae.get(idx[i]));
                if ha as usize <= budget {
                    ys[n] = ys[i];
                    ms[n] = ms[i];
                    hcs[n] = hcs[i];
                    has[n] = ha;
                    n += 1;
                }
            }
        }
        // The joint corner/B-centre/edge-triple table only adds information near the leaves.
        if budget <= 8 {
            if let Some(cbe) = &t.cbe {
                let mut ci = [[0usize; 4]; N_MOVES];
                for i in 0..n {
                    cbe.prefetch_all(f, &ys[i], &mut ci[i]);
                }
                let mut k = 0;
                for i in 0..n {
                    if cbe.bound(&ci[i]) as usize <= budget {
                        ys[k] = ys[i];
                        ms[k] = ms[i];
                        hcs[k] = hcs[i];
                        has[k] = has[i];
                        k += 1;
                    }
                }
                n = k;
            }
        }
        // Class-B bound against the best solution so far.
        let mut ybs = [*xb; N_MOVES];
        let mut hbas = [0u8; N_MOVES];
        let mut hbcs = [0u8; N_MOVES];
        if B && n > 0 {
            let best = self.best_len.load(Ordering::Relaxed);
            // a child at depth g+1 with bound h cannot give a solution shorter than g+1+h
            let lim = best.saturating_sub(g + 2);
            let emb = masks(&xb.e);
            for i in 0..n {
                ybs[i] = f.apply(xb, &emb, self.mb[ms[i] as usize]);
            }
            t.ae_ix.index_batch_fast(f, &ybs[..n], &mut idx[..n]);
            for i in 0..n {
                prefetch(&t.ae.data[idx[i] >> 5]);
                cbi[i] = CbTable::index_fast(f, &ybs[i]);
                prefetch(&t.cb.tab.data[cbi[i] >> 5]);
            }
            let mut k = 0;
            for i in 0..n {
                let hba = decode(hb_ae, t.ae.get(idx[i]));
                let hbc = decode(hb_cb, t.cb.tab.get(cbi[i]));
                if (hba.max(hbc) as usize) <= lim {
                    ys[k] = ys[i];
                    ms[k] = ms[i];
                    hcs[k] = hcs[i];
                    has[k] = has[i];
                    ybs[k] = ybs[i];
                    hbas[k] = hba;
                    hbcs[k] = hbc;
                    k += 1;
                }
            }
            n = k;
        }
        // Ordering statistic of the task: how close its subtree gets to the phase-1 target.
        for i in 0..n {
            *stat += STAT_W[has[i].max(hcs[i]) as usize];
        }
        for i in 0..n {
            path.push(ms[i] as usize);
            self.dfs::<B>(
                v, path, &ys[i], has[i], hcs[i], &ybs[i], hbas[i], hbcs[i], bound, nodes, stat,
            );
            path.pop();
            if self.stop.load(Ordering::Relaxed) {
                return;
            }
        }
    }

    /// Children of a task that survive pruning under `bound`.
    fn expand(&self, t: &Task, bound: usize) -> Vec<Task> {
        let g = t.path.len();
        let last = t.path.last().copied().unwrap_or(usize::MAX);
        let budget = bound - g - 1;
        let mut out = vec![];
        for m in 0..N_MOVES {
            if !allowed_after(last, m) || (budget == 0 && is_a_move(m)) {
                continue;
            }
            if let Some(c) = self.child(t, m, budget) {
                out.push(c);
            }
        }
        out
    }

    /// Tasks for one iteration: the roots, expanded level by level (in parallel) until there
    /// are enough tasks for the thread pool. The depth stays small because most of an IDA*
    /// iteration's work lies in the middle levels, which must not be explored serially.
    fn tasks(&self, bound: usize, threads: usize) -> Vec<Task> {
        let mut tasks: Vec<Task> = self
            .variants
            .iter()
            .enumerate()
            .filter(|(_, var)| var.h_ae.max(var.h_cb) as usize <= bound)
            .map(|(v, var)| Task {
                v,
                path: vec![],
                x: var.root,
                h_ae: var.h_ae,
                h_cb: var.h_cb,
                xb: var.root_b,
                hb_ae: var.hb_ae,
                hb_cb: var.hb_cb,
            })
            .collect();
        let max_k = bound.saturating_sub(4).min(MAX_SPLIT);
        let mut k = 0;
        // (single-threaded runs are split too, so per-task decisions can change mid-iteration)
        while tasks.len() < 64 * threads && k < max_k {
            tasks = tasks
                .par_iter()
                .flat_map_iter(|t| self.expand(t, bound))
                .collect();
            k += 1;
        }
        tasks
    }
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
    let variants: Vec<Variant> = [IDENTITY, SWAP_AB]
        .iter()
        .map(|sym| {
            let c = if *sym == IDENTITY {
                cube
            } else {
                cube.conjugate(sym)
            };
            let root = fast().of(&c);
            let (h_ae, h_cb) = t.exact_fast(&root);
            let root_b = fast().of(&c.conjugate(&SWAP_AB));
            let (hb_ae, hb_cb) = t.exact_fast(&root_b);
            Variant {
                cube: c,
                sym: *sym,
                root,
                h_ae,
                h_cb,
                root_b,
                hb_ae,
                hb_cb,
            }
        })
        .collect();
    // Every solution is a prefix followed by class-A moves, so one variant would suffice to
    // prove optimality; both are kept because the earlier a short solution is found, the more
    // the class-B bound prunes (measured: one variant was 1.9-3.9x slower in optimal mode).
    let sh = Shared {
        t,
        variants,
        frame,
        limits: limits.clone(),
        report,
        start,
        best_len: AtomicUsize::new(usize::MAX),
        count: AtomicU64::new(0),
        nodes: AtomicU64::new(0),
        stop: AtomicBool::new(false),
        best: Mutex::new((None, vec![])),
        mb: swap_moves(),
    };
    let lo = sh
        .variants
        .iter()
        .map(|v| v.h_ae.max(v.h_cb) as usize)
        .min()
        .unwrap();
    let mut optimal = false;
    let threads = rayon::current_num_threads();
    let mut depth = lo;
    let mut prev_stats: std::collections::HashMap<(usize, [u8; 8]), f32> =
        std::collections::HashMap::new();
    loop {
        if depth >= sh.best_len.load(Ordering::Relaxed) {
            optimal = true;
            break;
        }
        let mut tasks = sh.tasks(depth, threads);
        // Order the tasks by the statistics of the same prefixes in the previous iteration, so
        // that the most promising subtrees are searched first (egg-solver's job ordering).
        let key = |t: &Task| -> (usize, [u8; 8]) {
            let mut k = [u8::MAX; 8];
            for (i, &m) in t.path.iter().take(8).enumerate() {
                k[i] = m as u8;
            }
            (t.v, k)
        };
        if !prev_stats.is_empty() {
            let mut keyed: Vec<(f32, Task)> = tasks
                .into_iter()
                .map(|t| (*prev_stats.get(&key(&t)).unwrap_or(&0.0), t))
                .collect();
            keyed.sort_by(|a, b| b.0.total_cmp(&a.0));
            tasks = keyed.into_iter().map(|(_, t)| t).collect();
        }
        // Pruning against the best length changes how many (non-improving) solutions are
        // reached, so it is off when the user counts solutions.
        let count_mode = limits.max_solutions.is_some();
        let t_split = start.elapsed().as_secs_f64();
        let stats: Vec<std::sync::atomic::AtomicU32> = (0..tasks.len())
            .map(|_| std::sync::atomic::AtomicU32::new(0))
            .collect();
        let next = std::sync::atomic::AtomicUsize::new(0);
        // Workers take tasks strictly in order from a shared counter.
        rayon::scope(|scope| {
            for _ in 0..threads {
                scope.spawn(|_| {
                    let mut path: Vec<usize> = Vec::with_capacity(64);
                    // the per-thread counter persists across tasks so time checks happen regularly
                    let mut nodes = 0u64;
                    loop {
                        let i = next.fetch_add(1, Ordering::Relaxed);
                        if i >= tasks.len() || sh.stop.load(Ordering::Relaxed) || sh.time_up() {
                            break;
                        }
                        let task = &tasks[i];
                        path.clear();
                        path.extend_from_slice(&task.path);
                        let before = nodes;
                        let mut stat = 0f32;
                        let use_b = !count_mode && sh.best_len.load(Ordering::Relaxed) < usize::MAX;
                        if use_b {
                            sh.dfs::<true>(
                                task.v, &mut path, &task.x, task.h_ae, task.h_cb, &task.xb,
                                task.hb_ae, task.hb_cb, depth, &mut nodes, &mut stat,
                            );
                        } else {
                            sh.dfs::<false>(
                                task.v, &mut path, &task.x, task.h_ae, task.h_cb, &task.xb,
                                task.hb_ae, task.hb_cb, depth, &mut nodes, &mut stat,
                            );
                        }
                        stats[i].store(stat.to_bits(), Ordering::Relaxed);
                        sh.nodes.fetch_add(nodes - before, Ordering::Relaxed);
                    }
                });
            }
        });
        prev_stats = tasks
            .iter()
            .zip(stats.iter())
            .map(|(t, s)| (key(t), f32::from_bits(s.load(Ordering::Relaxed))))
            .collect();
        if std::env::var_os("FTO_DEBUG").is_some() {
            eprintln!("  depth {depth}: split done at {t_split:.4}s");
            eprintln!(
                "  depth {depth}: {} tasks, {} nodes, {} solutions, best {:?}, {:.3}s",
                tasks.len(),
                sh.nodes.load(Ordering::Relaxed),
                sh.count.load(Ordering::Relaxed),
                sh.best_len.load(Ordering::Relaxed),
                start.elapsed().as_secs_f64()
            );
        }
        if sh.stop.load(Ordering::Relaxed) || sh.time_up() {
            break;
        }
        depth += 1;
    }
    let (sol, events) = sh.best.into_inner().unwrap();
    let length = sol.as_ref().map(|s| s.len());
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
        elapsed: start.elapsed().as_secs_f64(),
        optimal,
        events,
    }
}
