//! The radius-N phase-2 distance: the shortest solution with every move allowed in which every
//! position along the way has relaxed phase-1 distance at most `N` (its A-centres and edges at
//! most `N` moves from a configuration of the phase-2 space). `N = 0` is the block rule's
//! distance only in the sense of the positions allowed (a move inside the configurations is a
//! block move, since a cut hexagon leaves the configurations); larger `N` lets a solution break
//! the blocks for a while.
//!
//! IDA* over the fixed-frame moves, bounded below by the optimal solver's tables (any lower
//! bound of the distance with every move allowed is one of the radius-N distance), each child
//! dropped when its relaxed phase-1 distance exceeds `N` ([`NearSet`]: the A-centre/edge
//! entries within `N` of a target). The first iterations can be skipped with a known lower bound
//! (the distance with every move allowed). On a time limit the result is the bound proven.

use super::joint::AeDepths;
use crate::coord::Coords;
use crate::cube::{Cube, N_MOVES};
use crate::optimal::coords::{opt_moves, Pos};
use crate::optimal::tables::{exact, Bounds, MAX_K};
use crate::prune::decode;
use crate::tables::{AeIndex, AeSym};
use rayon::prelude::*;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::time::Instant;

/// The relaxed phase-1 distance of the A-centre/edge entries within `n` of a target.
pub struct NearSet {
    /// the largest distance held
    pub n: u8,
    /// (entry, distance), sorted by entry
    dist: Vec<(usize, u8)>,
}

impl NearSet {
    /// The entries at distance at most `n` (from the exact distances).
    pub fn new(dep: &AeDepths, n: u8) -> NearSet {
        let mut dist: Vec<(usize, u8)> = dep
            .depth
            .par_iter()
            .enumerate()
            .filter_map(|(i, d)| {
                let d = d.load(Ordering::Relaxed);
                (d <= n).then_some((i, d))
            })
            .collect();
        dist.sort_unstable();
        NearSet { n, dist }
    }

    /// Entries in the set.
    pub fn len(&self) -> usize {
        self.dist.len()
    }

    /// Whether the set is empty.
    pub fn is_empty(&self) -> bool {
        self.dist.is_empty()
    }

    /// The relaxed phase-1 distance of entry `i` if it is at most [`NearSet::n`].
    #[inline]
    pub fn of_entry(&self, i: usize) -> Option<u8> {
        self.dist
            .binary_search_by_key(&i, |e| e.0)
            .ok()
            .map(|k| self.dist[k].1)
    }

    /// The relaxed phase-1 distance of a position if it is at most [`NearSet::n`].
    pub fn d1(&self, ix: &AeSym, p: &Pos) -> Option<u8> {
        let mut ep = [0u8; 12];
        ep.copy_from_slice(&p.ep[..12]);
        let mut xa = [0u8; 12];
        xa.copy_from_slice(&p.xa[..12]);
        let c = Cube {
            ep,
            xa,
            ..Cube::solved()
        };
        self.of_entry(ix.index(&Coords::of(&c)))
    }
}

/// What a search found.
#[derive(Clone, Debug)]
pub struct RadiusResult {
    /// the largest length proven impossible plus one (a lower bound)
    pub lower: usize,
    /// the shortest solution found (internal moves), if any
    pub solution: Option<Vec<usize>>,
    /// nodes searched
    pub nodes: u64,
    /// whether the search stopped at the time limit
    pub timed_out: bool,
}

/// The search: tables, index and radius.
pub struct Radius<'a, B: Bounds> {
    /// lower bounds with every move allowed
    pub b: &'a B,
    /// the A-centre/edge index of the relaxed table
    pub ix: &'a AeSym,
    /// the relaxed phase-1 distances near the targets (to at least `n`)
    pub near: &'a NearSet,
    /// the radius
    pub n: u8,
}

struct Ctx<'a> {
    stop: &'a AtomicBool,
    found: &'a AtomicBool,
    deadline: Instant,
    nodes: &'a AtomicU64,
}

impl<B: Bounds> Radius<'_, B> {
    /// The children of `p` (exact table values `h`) that may lie on a solution of `left` more
    /// moves after them: within the radius and no table value above `left`. Lookup 0 of the
    /// bounds must be the relaxed table's A-centre/edge index (the large tier's first lookup
    /// is): its index tells the radius. As in the optimal solver, each stage computes and
    /// prefetches the indices of all the children still alive before reading any value.
    /// Returns how many of `out` are filled.
    #[inline]
    fn children(
        &self,
        p: &Pos,
        h: &[u8; MAX_K],
        left: u8,
        last: usize,
        out: &mut [(usize, Pos, [u8; MAX_K]); N_MOVES],
    ) -> usize {
        let mt = opt_moves();
        let mut pre = [0u32; N_MOVES];
        let mut idx = [[0usize; MAX_K]; N_MOVES];
        let mut n = 0;
        for m in 0..N_MOVES {
            // (a face turned twice in a row is one turn of that face, with no position between)
            if last < N_MOVES && m / 2 == last / 2 {
                continue;
            }
            let q = mt.apply(p, m);
            pre[n] = self.b.pre(0, &q);
            out[n] = (m, q, [0; MAX_K]);
            n += 1;
        }
        // the first lookup: its index is the relaxed index (the radius), then its value
        let t0 = self.b.table(0);
        let mut alive = 0;
        let mut sel = [0usize; N_MOVES];
        for u in 0..n {
            let i0 = self.b.index_with(0, &out[u].1, pre[u]);
            if self.near.of_entry(i0).is_none_or(|d| d > self.n) {
                continue;
            }
            idx[u][0] = i0;
            crate::tables::prefetch_at(&t0.data, i0 >> 5);
            sel[alive] = u;
            alive += 1;
        }
        let order = &B::ORDER[..B::K.min(B::ORDER.len())];
        for (st, &k) in order.iter().enumerate() {
            let tab = self.b.table(k);
            let nk = order.get(st + 1).copied();
            let mut next = 0;
            for j in 0..alive {
                let u = sel[j];
                let v = decode(h[k], tab.get(idx[u][k]));
                out[u].2[k] = v;
                if v <= left {
                    if let Some(nk) = nk {
                        let ix = self.b.index(nk, &out[u].1);
                        idx[u][nk] = ix;
                        crate::tables::prefetch_at(&self.b.table(nk).data, ix >> 5);
                    }
                    sel[next] = u;
                    next += 1;
                }
            }
            alive = next;
        }
        // the survivors, in move order
        for j in 0..alive {
            out[j] = out[sel[j]];
        }
        alive
    }

    fn dfs(
        &self,
        p: &Pos,
        h: &[u8; MAX_K],
        left: u8,
        last: usize,
        path: &mut Vec<usize>,
        cx: &Ctx,
        nodes: &mut u64,
    ) -> bool {
        *nodes += 1;
        if nodes.trailing_zeros() >= 16 {
            cx.nodes.fetch_add(0x1_0000, Ordering::Relaxed);
            if Instant::now() >= cx.deadline {
                cx.stop.store(true, Ordering::Relaxed);
            }
        }
        if cx.stop.load(Ordering::Relaxed) || cx.found.load(Ordering::Relaxed) {
            return false;
        }
        if left == 0 {
            return h.iter().take(B::K).all(|&v| v == 0) && *p == Pos::of(&Cube::solved());
        }
        let mut kids = [(0usize, *p, [0u8; MAX_K]); N_MOVES];
        let n = self.children(p, h, left - 1, last, &mut kids);
        for &(m, q, hq) in &kids[..n] {
            path.push(m);
            if self.dfs(&q, &hq, left - 1, m, path, cx, nodes) {
                return true;
            }
            path.pop();
        }
        false
    }

    /// The radius-N distance of `cube` (fixed frame), searching lengths `from..=max` (`from`: a
    /// known lower bound, e.g. the distance with every move allowed), until `seconds`.
    pub fn solve(
        &self,
        cube: &Cube,
        from: usize,
        max: usize,
        seconds: f64,
    ) -> Option<RadiusResult> {
        let root = Pos::of(cube);
        let h0 = exact(self.b, &root)?;
        if self.near.of_entry(self.b.index(0, &root)) != self.near.d1(self.ix, &root)
            || self.near.d1(self.ix, &root).is_none_or(|d| d > self.n)
        {
            // (the first lookup is not the relaxed index, or the root is outside the radius)
            return None;
        }
        let solved = Pos::of(&Cube::solved());
        let start = Instant::now();
        let deadline = start + std::time::Duration::from_secs_f64(seconds);
        let stop = AtomicBool::new(false);
        let nodes = AtomicU64::new(0);
        let hmax = h0.iter().take(B::K).copied().max().unwrap_or(0) as usize;
        if root == solved {
            return Some(RadiusResult {
                lower: 0,
                solution: Some(vec![]),
                nodes: 0,
                timed_out: false,
            });
        }
        let mut bound = from.max(hmax).max(1);
        while bound <= max {
            // tasks: the nodes 3 moves down (or fewer for a short bound)
            let split = 3.min(bound - 1);
            let mut level: Vec<(Vec<usize>, Pos, [u8; MAX_K])> = vec![(vec![], root, h0)];
            for _ in 0..split {
                let mut next = vec![];
                for (path, p, h) in &level {
                    let left = (bound - path.len() - 1) as u8;
                    let last = path.last().copied().unwrap_or(N_MOVES);
                    let mut kids = [(0usize, *p, [0u8; MAX_K]); N_MOVES];
                    let n = self.children(p, h, left, last, &mut kids);
                    for &(m, q, hq) in &kids[..n] {
                        let mut np = path.clone();
                        np.push(m);
                        next.push((np, q, hq));
                    }
                }
                level = next;
            }
            let found = AtomicBool::new(false);
            let cx = Ctx {
                stop: &stop,
                found: &found,
                deadline,
                nodes: &nodes,
            };
            let sol = level.par_iter().find_map_any(|(path, p, h)| {
                let mut path = path.clone();
                let last = path.last().copied().unwrap_or(N_MOVES);
                let left = (bound - path.len()) as u8;
                let mut n = 0u64;
                let ok = self.dfs(p, h, left, last, &mut path, &cx, &mut n);
                cx.nodes.fetch_add(n & 0xffff, Ordering::Relaxed);
                if ok {
                    found.store(true, Ordering::Relaxed);
                    Some(path)
                } else {
                    None
                }
            });
            if let Some(path) = sol {
                return Some(RadiusResult {
                    lower: bound,
                    solution: Some(path),
                    nodes: nodes.load(Ordering::Relaxed),
                    timed_out: false,
                });
            }
            if stop.load(Ordering::Relaxed) {
                return Some(RadiusResult {
                    lower: bound,
                    solution: None,
                    nodes: nodes.load(Ordering::Relaxed),
                    timed_out: true,
                });
            }
            bound += 1;
        }
        Some(RadiusResult {
            lower: bound,
            solution: None,
            nodes: nodes.load(Ordering::Relaxed),
            timed_out: false,
        })
    }
}
