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
        // The CB table (max 8) prunes almost nothing with 6+ moves left: check AE first there.
        let ae_first = budget >= 6;
        let mut cbe_done = false;
        let mut b_done = false;
        let mut ybs = [*xb; N_MOVES];
        let mut hbas = [0u8; N_MOVES];
        let mut hbcs = [0u8; N_MOVES];
        let mut bidx = [0usize; N_MOVES];
        let mut bcbi = [0usize; N_MOVES];
        if ae_first {
            t.ae_ix.index_batch_fast_pf(f, &ys[..n], &mut idx[..n], &t.ae.data);
            // Without a joint-table stage to hide the CB latency, fetch the CB words of all
            // children while the AE words are in flight.
            let early = budget > 8;
            if early {
                for i in 0..n {
                    cbi[i] = CbTable::index_fast(f, &ys[i]);
                    prefetch(&t.cb.tab.data[cbi[i] >> 5]);
                }
            }
            let mut k = 0;
            for i in 0..n {
                let ha = decode(h_ae, t.ae.get(idx[i]));
                if ha as usize <= budget {
                    ys[k] = ys[i];
                    ms[k] = ms[i];
                    has[k] = ha;
                    if early {
                        cbi[k] = cbi[i];
                    } else {
                        cbi[k] = CbTable::index_fast(f, &ys[k]);
                        prefetch(&t.cb.tab.data[cbi[k] >> 5]);
                    }
                    k += 1;
                }
            }
            // The class-B lookups are issued here too, so that their latency overlaps the
            // CB and joint-table stages.
            let lim = if B {
                // a child at depth g+1 with bound h cannot give a solution shorter than g+1+h
                self.best_len.load(Ordering::Relaxed).saturating_sub(g + 2)
            } else {
                0
            };
            if B && k > 0 {
                let emb = masks(&xb.e);
                for i in 0..k {
                    ybs[i] = f.apply(xb, &emb, self.mb[ms[i] as usize]);
                }
                t.ae_ix.index_batch_fast_pf(f, &ybs[..k], &mut bidx[..k], &t.ae.data);
                for i in 0..k {
                    bcbi[i] = CbTable::index_fast(f, &ybs[i]);
                    prefetch(&t.cb.tab.data[bcbi[i] >> 5]);
                }
            }
            // The joint table goes before reading the CB words, so their latency overlaps it.
            if budget <= 8 {
                if let Some(cbe) = &t.cbe {
                    let mut ci = [[0usize; 4]; N_MOVES];
                    for i in 0..k {
                        cbe.prefetch_all(f, &ys[i], &mut ci[i]);
                    }
                    let mut k2 = 0;
                    for i in 0..k {
                        if cbe.bound(&ci[i]) as usize <= budget {
                            ys[k2] = ys[i];
                            ms[k2] = ms[i];
                            has[k2] = has[i];
                            cbi[k2] = cbi[i];
                            if B {
                                ybs[k2] = ybs[i];
                                bidx[k2] = bidx[i];
                                bcbi[k2] = bcbi[i];
                            }
                            k2 += 1;
                        }
                    }
                    k = k2;
                    cbe_done = true;
                }
            }
            n = 0;
            for i in 0..k {
                let hc = decode(h_cb, t.cb.tab.get(cbi[i]));
                if hc as usize > budget {
                    continue;
                }
                if B {
                    let hba = decode(hb_ae, t.ae.get(bidx[i]));
                    let hbc = decode(hb_cb, t.cb.tab.get(bcbi[i]));
                    if hba.max(hbc) as usize > lim {
                        continue;
                    }
                    ybs[n] = ybs[i];
                    hbas[n] = hba;
                    hbcs[n] = hbc;
                }
                ys[n] = ys[i];
                ms[n] = ms[i];
                has[n] = has[i];
                hcs[n] = hc;
                n += 1;
            }
            b_done = true;
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
        if budget <= 8 && !cbe_done {
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
        if B && n > 0 && !b_done {
            let best = self.best_len.load(Ordering::Relaxed);
            // a child at depth g+1 with bound h cannot give a solution shorter than g+1+h
            let lim = best.saturating_sub(g + 2);
            let emb = masks(&xb.e);
            for i in 0..n {
                ybs[i] = f.apply(xb, &emb, self.mb[ms[i] as usize]);
            }
            t.ae_ix.index_batch_fast_pf(f, &ybs[..n], &mut idx[..n], &t.ae.data);
            for i in 0..n {
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
        tb.ae_ix.index_batch_fast_pf(f, &ys[..n], &mut ia[..n], &tb.ae.data);
        tb.ae_ix.index_batch_fast_pf(f, &yb[..n], &mut iab[..n], &tb.ae.data);
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
                let (lo, hi) = tree.kids[i as usize].unwrap();
                next.extend((lo..hi).filter(|&c| alive(&tree.nodes[c as usize])));
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
    let mut tree = Tree::default();
    // ordering statistic of each tree node from the previous iteration (0 if not a task then)
    let mut node_stat: Vec<f32> = vec![];
    let mut prev_ids: Vec<u32> = vec![];
    let mut prev_nodes: Option<u64> = None;
    loop {
        if depth >= sh.best_len.load(Ordering::Relaxed) {
            optimal = true;
            break;
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
        // that the most promising subtrees are searched first (egg-solver's job ordering).
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
        let tree_ref = &tree;
        let worker = || {
            let mut path: Vec<usize> = Vec::with_capacity(64);
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
                "  depth {depth}: {} tasks, {} nodes, {} solutions, best {:?}, {:.3}s",
                n_tasks,
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
