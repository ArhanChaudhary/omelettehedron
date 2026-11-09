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
}

#[derive(Clone)]
struct Task {
    v: usize,
    path: Vec<usize>,
    x: Fc,
    h_ae: u8,
    h_cb: u8,
}

const CHECK_EVERY: u64 = 1 << 14;

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
    fn child(&self, x: &Fc, h_ae: u8, h_cb: u8, m: usize, budget: usize) -> Option<(Fc, u8, u8)> {
        let t = self.t;
        let f = fast();
        let y = f.apply(x, &masks(&x.e), m);
        let hc = decode(h_cb, t.cb.tab.get(CbTable::index_fast(f, &y)));
        if hc as usize > budget {
            return None;
        }
        let ha = decode(h_ae, t.ae.get(t.ae_ix.index_fast(f, &y, f.a_rank(&y))));
        if ha as usize > budget {
            return None;
        }
        Some((y, ha, hc))
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
    fn dfs(
        &self,
        v: usize,
        path: &mut Vec<usize>,
        x: &Fc,
        h_ae: u8,
        h_cb: u8,
        bound: usize,
        nodes: &mut u64,
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
            t.ae_ix.index_batch_fast(f, &ys[..n], &mut idx[..n]);
            for &i in &idx[..n] {
                prefetch(&t.ae.data[i >> 5]);
            }
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
            t.ae_ix.index_batch_fast(f, &ys[..k], &mut idx[..k]);
            for &i in &idx[..k] {
                prefetch(&t.ae.data[i >> 5]);
            }
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
        for i in 0..n {
            path.push(ms[i] as usize);
            self.dfs(v, path, &ys[i], has[i], hcs[i], bound, nodes);
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
            if let Some((y, ha, hc)) = self.child(&t.x, t.h_ae, t.h_cb, m, budget) {
                let mut p = t.path.clone();
                p.push(m);
                out.push(Task {
                    v: t.v,
                    path: p,
                    x: y,
                    h_ae: ha,
                    h_cb: hc,
                });
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
            })
            .collect();
        let max_k = bound.saturating_sub(4).min(MAX_SPLIT);
        let mut k = 0;
        while threads > 1 && tasks.len() < 64 * threads && k < max_k {
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
            Variant {
                cube: c,
                sym: *sym,
                root,
                h_ae,
                h_cb,
            }
        })
        .collect();
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
    loop {
        if depth >= sh.best_len.load(Ordering::Relaxed) {
            optimal = true;
            break;
        }
        let tasks = sh.tasks(depth, threads);
        let t_split = start.elapsed().as_secs_f64();
        tasks.par_iter().for_each_init(
            || (Vec::with_capacity(64), 0u64),
            |(path, nodes), task| {
                if sh.stop.load(Ordering::Relaxed) || sh.time_up() {
                    return;
                }
                path.clear();
                path.extend_from_slice(&task.path);
                let before = *nodes;
                // the per-thread counter persists across tasks so time checks happen regularly
                sh.dfs(task.v, path, &task.x, task.h_ae, task.h_cb, depth, nodes);
                sh.nodes.fetch_add(*nodes - before, Ordering::Relaxed);
            },
        );
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
