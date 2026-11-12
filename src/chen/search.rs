//! The block-split search: Kociemba's two-phase loop over the split of [`super::space`].
//!
//! * **Phase 1**: IDA* over all 16 internal moves to `H'` (A-centres and edges solved), with
//!   the two-phase solver's A-centre/edge table as the bound. The tables of the medium tier
//!   and above give the exact distance; the small and tiny tiers give lower bounds, so every
//!   leaf is checked in full. A phase-1 path must end in a class-B move (a class-A move keeps a
//!   position in or out of `H'`, so a path that ends in one is a shorter phase-1 path plus a
//!   phase-2 move).
//! * **Phase 2**, from each phase-1 leaf: IDA* with the allowed moves only, bounded by the
//!   AC, SB and CC tables ([`super::tables`]), up to `best - phase-1 length - 1` moves. A
//!   phase-2 path ends when it reaches `H_A` (with a class-B move, or at once): the two-phase
//!   solver's exact `H_A` table then gives the rest. Every allowed path to solved ends in a run
//!   of class-A moves inside `H_A`, so this finds the shortest allowed completion of the leaf.
//!
//! Phase 1 runs for two targets, like the two-phase solver: `H'` itself (side 0), and the
//! class-B version, through the position conjugated by `SWAP_AB` (side 1). The phase-1 depth
//! grows by one per iteration; each iteration's top levels are split into tasks that the
//! threads take in turn, and all threads share the best length.
//!
//! The split restricts phase 2, so a solution it finds need not be optimal, and running out of
//! phase-1 depths below the best proves nothing: [`Outcome::optimal`] is always false.

use super::space::*;
use super::tables::*;
use crate::cube::*;
use crate::fast::*;
use crate::geometry::{transpose, Mat, IDENTITY};
use crate::notation::*;
use crate::progress::{rate, Bound, Watch};
use crate::prune::decode;
use crate::search::{Event, Limits, Outcome, Report};
use crate::tables::{prefetch_at, AeIndex, Tables};
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::sync::{Mutex, PoisonError};
use std::time::Instant;

/// Longest phase 2 tried before the first solution is known (the restricted distance is far
/// below this; it only keeps a search finite).
const P2_CAP: usize = 40;

/// Phase-1 levels expanded into tasks (at most).
const SPLIT: usize = 3;

/// Nodes between two checks of the stop flag and the clock.
const CHECK_EVERY: u32 = 1 << 12;

/// One target of phase 1 (side 0: the position as given, target `H'`; side 1: conjugated by
/// `SWAP_AB`, target the class-B version of `H'`).
struct Side {
    cube: Cube,
    sym: Mat,
}

/// A phase-1 node at the split depth.
struct Task {
    side: usize,
    path: Vec<usize>,
    x: Fc,
    h: u8,
}

/// The tasks of one phase-1 depth.
struct Level {
    depth: usize,
    tasks: Vec<Task>,
    /// tasks handed out and finished
    next: AtomicUsize,
    done: AtomicUsize,
    /// when it was made, and the node counts (phase 1, phase 2) then
    t0: f64,
    nodes: [u64; 2],
}

/// A phase-1 leaf in `H'`, ready for its phase-2 search.
struct Leaf {
    side: usize,
    path: Vec<usize>,
    p: P2,
    /// its exact phase-2 table distances
    h: [u8; N_TABLES],
}

impl Leaf {
    /// Its phase-2 lower bound.
    fn h2(&self) -> u8 {
        self.h.iter().copied().max().unwrap_or(0)
    }
}

/// Per-worker counters and caches.
#[derive(Default)]
struct Ctr {
    nodes1: u64,
    nodes2: u64,
    leaves: [u64; 2],
    /// nodes until the next check of the stop flag and the clock
    countdown: u32,
    /// the stop flag as last read
    stop: bool,
    /// phase-2 node count at which the current leaf is given up
    limit2: u64,
}

/// State shared by the threads of one solve.
struct Shared<'a, I: AeIndex> {
    t: &'a Tables<I>,
    ct: &'a ChenTables,
    s: &'static Space,
    f: &'static FastTables,
    sides: [Side; 2],
    frame: Frame,
    limits: Limits,
    report: Report,
    start: Instant,
    best_len: AtomicUsize,
    count: AtomicU64,
    nodes1: AtomicU64,
    nodes2: AtomicU64,
    leaves: [AtomicU64; 2],
    stop: &'a AtomicBool,
    watch: Watch<'a>,
    best: Mutex<(Option<Vec<PhysMove>>, Vec<Event>)>,
    /// phase-2 nodes per leaf before it is given up (see `leaf_budget`)
    budget: u64,
    /// `FTO_DEBUG`: a line per phase-2 iteration of each leaf on stderr
    debug: bool,
}

impl<I: AeIndex> Shared<'_, I> {
    fn time_up(&self) -> bool {
        let el = self.start.elapsed();
        let over = |t: Option<std::time::Duration>| t.is_some_and(|t| el >= t);
        if over(self.limits.time) || over(self.limits.cancel_after) {
            self.stop.store(true, Ordering::Relaxed);
            return true;
        }
        false
    }

    /// Count a node; every [`CHECK_EVERY`] nodes read the stop flag and the clock.
    #[inline]
    fn tick(&self, ctr: &mut Ctr) -> bool {
        if ctr.countdown == 0 {
            ctr.countdown = CHECK_EVERY;
            ctr.stop = self.stop.load(Ordering::Relaxed) || self.time_up();
        }
        ctr.countdown -= 1;
        ctr.stop
    }

    /// Phase-1 bound of a child `y` whose parent's bound is `parent`.
    #[inline]
    fn ae_bound(&self, y: &Fc, parent: u8) -> u8 {
        let t = self.t;
        match t.ae_ix.bound(self.f, y, &t.ae) {
            Some(h) => h,
            None => decode(
                parent,
                t.ae.get(t.ae_ix.index_fast(self.f, y, self.f.a_rank(y))),
            ),
        }
    }

    /// The children of a phase-1 node within `depth` (`g` moves made): (move, position,
    /// bound).
    fn children(&self, x: &Fc, h: u8, last: usize, g: usize, depth: usize) -> Vec<(usize, Fc, u8)> {
        let fw = self.f.flip_word(&masks(&x.e));
        (0..N_MOVES)
            .filter(|&m| allowed_after(false, last, m))
            .filter_map(|m| {
                let y = self.f.apply_w(x, fw, m);
                let hy = self.ae_bound(&y, h);
                (g + 1 + hy as usize <= depth).then_some((m, y, hy))
            })
            .collect()
    }

    /// Phase-1 IDA* below a node at depth `path.len()`, with bound `depth`.
    fn phase1(
        &self,
        side: usize,
        x: &Fc,
        h: u8,
        path: &mut Vec<usize>,
        depth: usize,
        ctr: &mut Ctr,
    ) {
        let g = path.len();
        let last = path.last().copied().unwrap_or(usize::MAX);
        if g == depth {
            if h == 0 && path.last().is_none_or(|&m| !is_a_move(m)) {
                if let Some(leaf) = self.leaf_check(side, path, ctr) {
                    self.leaf(&leaf, ctr);
                }
            }
            return;
        }
        if self.tick(ctr) {
            return;
        }
        ctr.nodes1 += 1;
        let fw = self.f.flip_word(&masks(&x.e));
        for m in 0..N_MOVES {
            if !allowed_after(false, last, m) {
                continue;
            }
            let y = self.f.apply_w(x, fw, m);
            let hy = self.ae_bound(&y, h);
            if g + 1 + hy as usize > depth {
                continue;
            }
            path.push(m);
            self.phase1(side, &y, hy, path, depth, ctr);
            path.pop();
            if ctr.stop {
                return;
            }
        }
    }

    /// A phase-1 leaf candidate: `path` (moves of the side's frame) reaches `H'` by the tables.
    /// Check it in full; the leaf with its phase-2 table distances if it is in `H'` and could
    /// still improve on the best.
    fn leaf_check(&self, side: usize, path: &[usize], ctr: &mut Ctr) -> Option<Leaf> {
        let cube = self.sides[side].cube.apply_seq(path);
        let p = P2::of(&cube)?;
        if !self.s.in_hprime(p.cfg) {
            return None;
        }
        ctr.leaves[side] += 1;
        let h = self.ct.exact(&p)?;
        let leaf = Leaf {
            side,
            path: path.to_vec(),
            p,
            h,
        };
        let best = self.best_len.load(Ordering::Relaxed);
        (best == usize::MAX || path.len() + usize::from(leaf.h2()) < best).then_some(leaf)
    }

    /// Phase 2 from a leaf: IDA* from its bound up to the longest phase 2 that improves on the
    /// best, within the leaf's node budget.
    fn leaf(&self, leaf: &Leaf, ctr: &mut Ctr) {
        let (side, path, p, h) = (leaf.side, &leaf.path[..], &leaf.p, leaf.h);
        let h2 = usize::from(leaf.h2());
        let d1 = path.len();
        let last = path.last().copied().unwrap_or(usize::MAX);
        let mut p2path = Vec::with_capacity(32);
        let mut b2 = h2;
        // the leaf's phase-2 node budget (see `leaf_budget`)
        ctr.limit2 = ctr.nodes2.saturating_add(self.budget);
        loop {
            // the longest phase 2 that would still improve on the best (re-read: other threads
            // may have found better)
            let best = self.best_len.load(Ordering::Relaxed);
            let lim = if best == usize::MAX {
                P2_CAP
            } else {
                best.saturating_sub(d1 + 1)
            };
            if b2 > lim || best <= d1 {
                return;
            }
            p2path.clear();
            let before = ctr.nodes2;
            let found = self.phase2(side, p, h, b2, last, path, &mut p2path, ctr);
            if self.debug {
                eprintln!(
                    "  leaf side {side} d1 {d1} h {h:?} bound {b2}: {} nodes{}",
                    ctr.nodes2 - before,
                    if found { ", found" } else { "" }
                );
            }
            if found || ctr.stop || ctr.nodes2 >= ctr.limit2 {
                return;
            }
            b2 += 1;
        }
    }

    /// Phase-2 IDA* below `p` (`p2path.len()` moves made, exact table distances `h`) with
    /// bound `bound`. Returns whether it found a solution (then recorded).
    fn phase2(
        &self,
        side: usize,
        p: &P2,
        h: [u8; N_TABLES],
        bound: usize,
        last: usize,
        p1path: &[usize],
        p2path: &mut Vec<usize>,
        ctr: &mut Ctr,
    ) -> bool {
        let g = p2path.len();
        // in H_A: the H_A table completes the solution (only after a class-B move: after a
        // class-A move the parent was in H_A too, with an H_A distance at most one more)
        if p2path.last().is_none_or(|&m| !is_a_move(m)) {
            if let Some(i) = self.ct.ha_index(self.s, p) {
                let d = self.t.p2.dist[i] as usize;
                if g + d <= bound {
                    self.record(side, p1path, p2path);
                    return true;
                }
            }
        }
        if g >= bound || ctr.nodes2 >= ctr.limit2 || self.tick(ctr) {
            return false;
        }
        ctr.nodes2 += 1;
        let (s, f, ct) = (self.s, self.f, self.ct);
        let left = bound - g - 1;
        // The children and their table indices first, each table word prefetched as soon as
        // its index is known, then the reads: 0.72x the time to the same first solutions of 4
        // random positions on one thread (4.91 -> 3.55 s, two rounds each, identical node
        // counts), against making and reading each child in turn.
        let mut kids = [*p; N_MOVES];
        let mut ms = [0u8; N_MOVES];
        let mut ix = [[0usize; N_TABLES]; N_MOVES];
        let mut n = 0;
        for m in 0..N_MOVES {
            if !allowed_after(false, last, m) {
                continue;
            }
            let Some(q) = p.apply(s, m) else { continue };
            let i = [
                ChenTables::ac_index(&q),
                ct.sb
                    .as_ref()
                    .map_or(0, |_| ChenTables::sb_index_with(s, f, &q)),
                ct.px
                    .as_ref()
                    .map_or(0, |_| ChenTables::px_index_with(s, f, &q)),
                ct.cc_index(&q),
            ];
            prefetch_at(&ct.ac.data, i[0] >> 5);
            if let Some(t) = &ct.sb {
                prefetch_at(&t.data, i[1] >> 5);
            }
            if let Some(t) = &ct.px {
                prefetch_at(&t.data, i[2] >> 5);
            }
            prefetch_at(&ct.cc.data, i[3] >> 5);
            kids[n] = q;
            ms[n] = m as u8;
            ix[n] = i;
            n += 1;
        }
        // (Searching the children in order of their bound instead made no difference: first
        // solutions of 4 random positions after 14.1M phase-2 nodes in all, against 14.4M.)
        for k in 0..n {
            let i = ix[k];
            let ac = decode(h[0], ct.ac.get(i[0]));
            if ac as usize > left {
                continue;
            }
            let px = ct.px.as_ref().map_or(0, |t| decode(h[2], t.get(i[2])));
            if px as usize > left {
                continue;
            }
            let sb = ct.sb.as_ref().map_or(0, |t| decode(h[1], t.get(i[1])));
            if sb as usize > left {
                continue;
            }
            let cc = decode(h[3], ct.cc.get(i[3]));
            if cc as usize > left {
                continue;
            }
            let m = usize::from(ms[k]);
            p2path.push(m);
            if self.phase2(
                side,
                &kids[k],
                [ac, sb, px, cc],
                bound,
                m,
                p1path,
                p2path,
                ctr,
            ) {
                return true;
            }
            p2path.pop();
            if ctr.stop {
                return false;
            }
        }
        false
    }

    /// Record a solution: phase 1 `p1path`, phase 2 `p2path` (both in the side's frame), then
    /// the `H_A` table's moves.
    fn record(&self, side: usize, p1path: &[usize], p2path: &[usize]) {
        let n = self.count.fetch_add(1, Ordering::Relaxed) + 1;
        if self.limits.max_solutions.is_some_and(|k| n >= k) {
            self.stop.store(true, Ordering::Relaxed);
        }
        let sd = &self.sides[side];
        let cube = sd.cube.apply_seq(p1path).apply_seq(p2path);
        let Some(p3) = self.t.p2.solve(&cube) else {
            debug_assert!(false, "H_A completion failed");
            return;
        };
        let total = p1path.len() + p2path.len() + p3.len();
        let mut g = self.best.lock().unwrap_or_else(PoisonError::into_inner);
        if total >= self.best_len.load(Ordering::Relaxed) {
            return;
        }
        let inv = transpose(&sd.sym);
        let internal: Vec<usize> = p1path
            .iter()
            .chain(p2path)
            .chain(&p3)
            .map(|&m| conjugate_move(m, &inv))
            .collect();
        debug_assert_eq!(
            self.sides[0].cube.apply_seq(&internal),
            Cube::solved(),
            "not a solution"
        );
        let phys = solution_to_physical(self.frame, &internal, self.t.moves.wide);
        let t = self.start.elapsed().as_secs_f64();
        if self.report == Report::Improvements {
            println!("{:>9.3}s  {:>2}  {}", t, total, format_moves(&phys));
        }
        self.watch.bound(Bound::Upper(&phys));
        self.watch.line(format_args!(
            "best {total} moves (phase 1 {}, phase 2 {} + {})",
            p1path.len(),
            p2path.len(),
            p3.len()
        ));
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

    /// A progress line for a level whose tasks are all done (the node counts include those of
    /// the levels running at the same time).
    fn level_done(&self, l: &Level, watch: &Watch) {
        let Some(p) = watch.progress else { return };
        let dt = self.start.elapsed().as_secs_f64() - l.t0;
        let d1 = self.nodes1.load(Ordering::Relaxed) - l.nodes[0];
        let d2 = self.nodes2.load(Ordering::Relaxed) - l.nodes[1];
        let best = self.best_len.load(Ordering::Relaxed);
        p.line(format_args!(
            "phase-1 depth {} {} in {dt:.3} s: {} tasks, {d1:.3e} phase-1 nodes, {d2:.3e} \
             phase-2 nodes since it began ({}); best {}",
            l.depth,
            if self.stop.load(Ordering::Relaxed) {
                "stopped"
            } else {
                "done"
            },
            l.tasks.len(),
            rate(d1 + d2, dt),
            if best == usize::MAX {
                "-".to_string()
            } else {
                best.to_string()
            }
        ));
    }

    /// The tasks of the iteration with bound `depth`: the nodes `min(depth, SPLIT)` moves
    /// deep that fit the bound, the two sides interleaved.
    fn tasks(&self, roots: &[(Fc, u8); 2], depth: usize) -> Vec<Task> {
        let k = depth.min(SPLIT);
        let mut per_side: [Vec<Task>; 2] = [vec![], vec![]];
        for side in 0..2 {
            let (x, h) = roots[side];
            if h as usize > depth {
                continue;
            }
            let mut front = vec![Task {
                side,
                path: vec![],
                x,
                h,
            }];
            for g in 0..k {
                let mut next = vec![];
                for t in &front {
                    let last = t.path.last().copied().unwrap_or(usize::MAX);
                    for (m, y, hy) in self.children(&t.x, t.h, last, g, depth) {
                        let mut path = t.path.clone();
                        path.push(m);
                        next.push(Task {
                            side,
                            path,
                            x: y,
                            h: hy,
                        });
                    }
                }
                front = next;
            }
            per_side[side] = front;
        }
        let [a, b] = per_side;
        let mut out = Vec::with_capacity(a.len() + b.len());
        let (mut ia, mut ib) = (a.into_iter(), b.into_iter());
        loop {
            match (ia.next(), ib.next()) {
                (None, None) => break,
                (x, y) => out.extend(x.into_iter().chain(y)),
            }
        }
        out
    }
}

/// Phase-2 nodes spent on one phase-1 leaf before it is given up for the next.
/// `FTO_CHEN_BUDGET` overrides (0: no limit). Without a limit the first leaves keep the search
/// for long: first solutions of 4 random positions (small tier, one thread) after 1.1e7, more
/// than 1.6e7 (none in 100 s), 1.7e7 and 1.2e8 phase-2 nodes (29, -, 29 and 31 moves); with
/// 2^16, 2^18 and 2^20 nodes per leaf after 1.6e7, 1.4e7 and 1.7e7 in all (27 or 28 moves
/// each).
fn leaf_budget() -> u64 {
    match std::env::var("FTO_CHEN_BUDGET")
        .ok()
        .and_then(|v| v.parse::<u64>().ok())
    {
        Some(0) => u64::MAX,
        Some(n) => n,
        None => 1 << 18,
    }
}

/// Solve one position with the block split: the two-phase solver's tables `t` (phase 1 and
/// the `H_A` table) and the phase-2 tables `ct`. `frame` is the physical frame of the
/// scrambled position. Only the face metric is supported.
pub fn solve<I: AeIndex>(
    t: &Tables<I>,
    ct: &ChenTables,
    cube: Cube,
    frame: Frame,
    limits: &Limits,
    report: Report,
    watch: &Watch,
) -> Outcome {
    let start = Instant::now();
    let own_stop = AtomicBool::new(false);
    let empty = |start: Instant| Outcome {
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
    let Some(cube_b) = cube.conjugate(&SWAP_AB) else {
        return empty(start);
    };
    let f = fast();
    let x = [f.of(&cube), f.of(&cube_b)];
    let roots = [(x[0], t.exact_fast(&x[0]).0), (x[1], t.exact_fast(&x[1]).0)];
    let sh = Shared {
        t,
        ct,
        s: space(),
        f,
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
        frame,
        limits: limits.clone(),
        report,
        start,
        best_len: AtomicUsize::new(usize::MAX),
        count: AtomicU64::new(0),
        nodes1: AtomicU64::new(0),
        nodes2: AtomicU64::new(0),
        leaves: [AtomicU64::new(0), AtomicU64::new(0)],
        stop: watch.stop.unwrap_or(&own_stop),
        watch: *watch,
        best: Mutex::new((None, vec![])),
        budget: leaf_budget(),
        debug: std::env::var_os("FTO_DEBUG").is_some(),
    };
    let threads = rayon::current_num_threads();
    let first = roots[0].1.min(roots[1].1) as usize;
    // Phase-1 depths as levels of tasks. A thread that finds no task left at any level makes
    // the next level instead of waiting for the others to finish theirs: with few phase-1
    // solutions per depth (9 at the first depth of a random position), waiting left most
    // threads idle. 48 threads, 10 random positions, small tier (two rounds each): first
    // solution after 0.0065 s instead of 0.078 s, best after 0.01 s 25.2 moves (9 of 10 solved)
    // instead of 26.0 (3 of 10), after 0.1 s 23.7 (10 of 10) instead of 24.6 (7 of 10); the
    // same after 1 and 2 s.
    let levels: Mutex<Vec<std::sync::Arc<Level>>> = Mutex::new(vec![]);
    let make_level = |depth: usize| -> Option<std::sync::Arc<Level>> {
        // a phase-1 path of `depth` moves gives at least `depth` moves
        if depth >= sh.best_len.load(Ordering::Relaxed)
            || limits.max_depth.is_some_and(|d| depth > d)
        {
            return None;
        }
        Some(std::sync::Arc::new(Level {
            depth,
            tasks: sh.tasks(&roots, depth),
            next: AtomicUsize::new(0),
            done: AtomicUsize::new(0),
            t0: start.elapsed().as_secs_f64(),
            nodes: [
                sh.nodes1.load(Ordering::Relaxed),
                sh.nodes2.load(Ordering::Relaxed),
            ],
        }))
    };
    // the next task: the first one not handed out, at the lowest level that has one
    let claim = || -> Option<(std::sync::Arc<Level>, usize)> {
        let mut lv = levels.lock().unwrap_or_else(PoisonError::into_inner);
        loop {
            // (a level as deep as the best can no longer improve on it)
            let best = sh.best_len.load(Ordering::Relaxed);
            if let Some(l) = lv
                .iter()
                .find(|l| l.depth < best && l.next.load(Ordering::Relaxed) < l.tasks.len())
            {
                let i = l.next.fetch_add(1, Ordering::Relaxed);
                return Some((l.clone(), i));
            }
            let depth = lv.last().map_or(first, |l| l.depth + 1);
            let l = make_level(depth)?;
            let empty = l.tasks.is_empty();
            lv.push(l.clone());
            if empty {
                sh.level_done(&l, watch);
            }
        }
    };
    let worker = || {
        let mut ctr = Ctr::default();
        let mut path = Vec::with_capacity(32);
        while !sh.stop.load(Ordering::Relaxed) && !sh.time_up() {
            let Some((l, i)) = claim() else { break };
            if let Some(task) = l.tasks.get(i) {
                path.clear();
                path.extend_from_slice(&task.path);
                sh.phase1(task.side, &task.x, task.h, &mut path, l.depth, &mut ctr);
            }
            // (the counters go to the shared totals per task, for the progress lines)
            sh.nodes1
                .fetch_add(std::mem::take(&mut ctr.nodes1), Ordering::Relaxed);
            sh.nodes2
                .fetch_add(std::mem::take(&mut ctr.nodes2), Ordering::Relaxed);
            for s in 0..2 {
                sh.leaves[s].fetch_add(std::mem::take(&mut ctr.leaves[s]), Ordering::Relaxed);
            }
            if l.done.fetch_add(1, Ordering::Relaxed) + 1 == l.tasks.len() {
                sh.level_done(&l, watch);
            }
        }
    };
    // (Collecting an iteration's leaves first and searching them in order of their phase-2
    // bound made no difference: first solutions of 4 random positions after 14.8M phase-2
    // nodes in all, against 14.4M.)
    if threads == 1 {
        worker();
    } else {
        rayon::scope(|scope| {
            for _ in 0..threads {
                scope.spawn(|_| worker());
            }
        });
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
        nodes: sh.nodes1.load(Ordering::Relaxed) + sh.nodes2.load(Ordering::Relaxed),
        leaves: leaves_by_side[0] + leaves_by_side[1],
        leaves_by_side,
        elapsed: start.elapsed().as_secs_f64(),
        optimal: false,
        events,
    }
}
