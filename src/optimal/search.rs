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
    pub elapsed: f64,
    pub iterations: Vec<Iteration>,
}

/// Moves allowed after each move (index `N_MOVES`: no previous move), as bit masks.
static ALLOWED: [u32; N_MOVES + 1] = {
    let mut a = [0u32; N_MOVES + 1];
    let mut last = 0;
    while last <= N_MOVES {
        let l = if last == N_MOVES { MAX_MOVES } else { last };
        let mut m = 0;
        while m < N_MOVES {
            if allowed_after(false, l, m) {
                a[last] |= 1 << m;
            }
            m += 1;
        }
        last += 1;
    }
    a
};

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
            N_MOVES
        } else {
            usize::from(self.path[usize::from(self.len) - 1])
        }
    }
}

/// Key of a move prefix (for the ordering statistics).
fn key(path: &[u8]) -> u32 {
    path.iter().rev().fold(0, |k, &m| k * 17 + u32::from(m) + 1)
}

/// Per-worker counters.
#[derive(Default)]
struct Counter {
    nodes: u64,
    /// `Shared::stop` as last read
    stop: bool,
}

/// Scratch arrays of one node's children (one frame per depth, reused).
#[derive(Clone)]
struct Scratch {
    ys: [Pos; N_MOVES],
    ms: [u8; N_MOVES],
    idx: [[usize; MAX_K]; N_MOVES],
    hs: [[u8; MAX_K]; N_MOVES],
    sel: [u8; N_MOVES],
}

impl Scratch {
    fn new() -> Self {
        let z = Pos::solved();
        Scratch {
            ys: [z; N_MOVES],
            ms: [0; N_MOVES],
            idx: [[0; MAX_K]; N_MOVES],
            hs: [[0; MAX_K]; N_MOVES],
            sel: [0; N_MOVES],
        }
    }
}

/// The state shared by the threads of one solve.
struct Shared<'a, B: Bounds> {
    b: &'a B,
    /// `opt_moves()`, read once (no `OnceLock` access per node)
    mt: &'static OptMoves,
    solved: Pos,
    cube: Cube,
    limits: Limits,
    start: Instant,
    stop: AtomicBool,
    nodes: AtomicU64,
    sols: Mutex<Vec<Vec<usize>>>,
}

impl<B: Bounds> Shared<'_, B> {
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
        let mut out = Vec::with_capacity(N_MOVES);
        let mut moves = ALLOWED[t.last()];
        let g = usize::from(t.len);
        while moves != 0 {
            let m = moves.trailing_zeros() as usize;
            moves &= moves - 1;
            let mut c = *t;
            c.path[g % MAX_SPLIT] = m as u8;
            c.len += 1;
            c.p = self.mt.apply(&t.p, m);
            out.push(c);
        }
        ctr.nodes += out.len() as u64;
        let mut idx = [[0usize; MAX_K]; N_MOVES];
        for (i, c) in out.iter().enumerate() {
            for k in 0..B::K {
                idx[i][k] = self.b.index(k, &c.p);
                prefetch_at(&self.b.table(k).data, idx[i][k] >> 5);
            }
        }
        for (i, c) in out.iter_mut().enumerate() {
            for k in 0..B::K {
                c.h[k] = decode(t.h[k], self.b.table(k).get(idx[i][k]));
            }
        }
        out.retain(|c| usize::from(c.h.iter().copied().max().unwrap_or(0)) + g < bound);
        out
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
            k += 1;
        }
        front
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
        frames: &mut [Scratch],
        stat: &mut f32,
    ) {
        let g = path.len();
        if g >= bound {
            if *p == self.solved {
                self.found(path, ctr);
            }
            return;
        }
        if ctr.stop {
            return;
        }
        // (the caller provides one frame per remaining level)
        let Some((fr, rest)) = frames.split_first_mut() else {
            debug_assert!(false, "out of search frames");
            return;
        };
        let budget = bound - g - 1;
        let mut moves = ALLOWED[last.min(N_MOVES)];
        let mut n = 0;
        while moves != 0 {
            let m = moves.trailing_zeros() as usize;
            moves &= moves - 1;
            fr.ys[n % N_MOVES] = self.mt.apply(p, m);
            fr.ms[n % N_MOVES] = m as u8;
            n += 1;
        }
        self.count(ctr, n);
        if budget == 0 {
            // The children are leaves: every table is 0 exactly at the solved state, so only a
            // solved child can pass; it is recognised without lookups.
            for i in 0..n {
                if fr.ys[i % N_MOVES] == self.solved {
                    path.push(fr.ms[i % N_MOVES]);
                    self.found(path, ctr);
                    path.pop();
                }
            }
            return;
        }
        // all indices and prefetches first, then the values
        for i in 0..n {
            for k in 0..B::K {
                let ix = self.b.index(k, &fr.ys[i % N_MOVES]);
                fr.idx[i % N_MOVES][k] = ix;
                prefetch_at(&self.b.table(k).data, ix >> 5);
            }
        }
        let mut kk = 0;
        for i in 0..n {
            let u = i % N_MOVES;
            let mut mx = 0;
            for k in 0..B::K {
                let v = decode(h[k], self.b.table(k).get(fr.idx[u][k]));
                fr.hs[u][k] = v;
                mx = mx.max(v);
            }
            let ok = usize::from(mx) <= budget;
            if ok {
                *stat += STAT_W[usize::from(mx) & 31];
            }
            fr.sel[kk % N_MOVES] = i as u8;
            kk += usize::from(ok);
        }
        for j in 0..kk {
            let i = usize::from(fr.sel[j % N_MOVES]) % N_MOVES;
            let m = fr.ms[i];
            path.push(m);
            let (y, hy) = (fr.ys[i], fr.hs[i]);
            self.dfs(&y, hy, usize::from(m), bound, path, ctr, rest, stat);
            path.pop();
            if ctr.stop {
                return;
            }
        }
    }
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
    let start = Instant::now();
    let p = Pos::of(&cube);
    let mut out = Outcome {
        internal: vec![],
        solutions: vec![],
        length: None,
        lower_bound: 0,
        stop: Stop::BadTables,
        nodes: 0,
        elapsed: 0.0,
        iterations: vec![],
    };
    let Some(h0) = exact(b, &p) else {
        out.elapsed = start.elapsed().as_secs_f64();
        return out;
    };
    let sh = Shared {
        b,
        mt: opt_moves(),
        solved: Pos::solved(),
        cube,
        limits: limits.clone(),
        start,
        stop: AtomicBool::new(false),
        nodes: AtomicU64::new(0),
        sols: Mutex::new(vec![]),
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
        // (for a future symmetric-start filter: a symmetric position has equivalent subtrees)
        let sym = super::sym::sym48();
        let n_sym = (0..super::sym::N_SYM)
            .filter(|&k| sym.conj_by_maps(&cube, k) == cube)
            .count();
        eprintln!(
            "root bounds: {} (search from depth {lo}); symmetries of the position: {n_sym} of 48",
            names.join(", ")
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
            let mut frames = vec![Scratch::new(); depth + 1];
            let mut ctr = Counter::default();
            loop {
                let i = next.fetch_add(1, Ordering::Relaxed);
                if i >= tasks.len() || sh.stop.load(Ordering::Relaxed) || sh.time_up() {
                    break;
                }
                let t = &tasks[i];
                path.clear();
                path.extend_from_slice(&t.path[..usize::from(t.len)]);
                let before = ctr.nodes;
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
            }
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
        out.iterations.push(Iteration {
            depth,
            nodes: it_nodes,
            tasks: tasks.len(),
            elapsed,
        });
        let n_sols = sh.sols.lock().unwrap_or_else(PoisonError::into_inner).len();
        if verbose {
            eprintln!(
                "  depth {depth:>2}: {it_nodes:>14} nodes, {:>5} tasks, {elapsed:>9.3} s{}",
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
        .map(|s| solution_to_physical(frame, s, false))
        .collect();
    out.nodes = sh.nodes.load(Ordering::Relaxed);
    out.elapsed = start.elapsed().as_secs_f64();
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn move_masks_match_rules() {
        for last in 0..=N_MOVES {
            for m in 0..N_MOVES {
                let l = if last == N_MOVES { usize::MAX } else { last };
                assert_eq!(ALLOWED[last] >> m & 1 == 1, allowed_after(false, l, m));
            }
        }
    }

    #[test]
    fn prefix_keys_are_distinct() {
        let mut seen = std::collections::HashSet::new();
        for len in 0..=MAX_SPLIT {
            let mut path = [0u8; MAX_SPLIT];
            for x in 0..16usize.pow(len as u32) {
                let mut y = x;
                for p in path.iter_mut().take(len) {
                    *p = (y % 16) as u8;
                    y /= 16;
                }
                assert!(seen.insert(key(&path[..len])));
            }
        }
    }
}
