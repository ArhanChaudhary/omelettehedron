//! Independent searches to check the disk search against.
//!
//! * [`slot_bfs`]: the same slots in memory, one slot at a time, with the general table
//!   builder ([`crate::prune::bfs`]: 2 bits per slot, forward and backward steps, equivalent
//!   slots through [`crate::prune::Graph::equivalents`]); for the instances that fit.
//! * [`position_bfs`]: whole positions in the fixed-corner frame without symmetry (the
//!   configuration move table, the corner and B-centre moves of the cube model), a hash set per
//!   layer; positions per depth, for the first depths.

use super::instance::{Instance, NO_EDGE, N_SYM};
use crate::chen::space::{space, NONE};
use crate::chen::tables::P2;
use crate::coord::mset;
use crate::cube::{Cube, N_MOVES};
use crate::prune::{bfs, Graph, Table2};
use std::collections::HashSet;

/// The slot graph of an instance (slot `((c * n_x) + x) * n_p + p`).
pub struct SlotGraph<'a> {
    /// the instance
    pub inst: &'a Instance,
}

impl SlotGraph<'_> {
    /// The index of slot `(c, x, p)`.
    #[inline]
    pub fn index(&self, c: usize, x: usize, p: usize) -> usize {
        (c * self.inst.n_x + x) * self.inst.n_p + p
    }
}

impl Graph for SlotGraph<'_> {
    const MOVES: usize = N_MOVES;
    type Node = (usize, usize, usize, usize);
    type Cache = ();
    fn node(&self, i: usize, (): &mut ()) -> (usize, usize, usize, usize) {
        let n_p = self.inst.n_p;
        let n_x = self.inst.n_x;
        (i, i / n_p / n_x, i / n_p % n_x, i % n_p)
    }
    #[inline]
    fn step(&self, &(i, c, x, p): &(usize, usize, usize, usize), m: usize) -> usize {
        let inst = self.inst;
        let (c2, k) = inst.edges[c][m];
        if c2 == NO_EDGE {
            return i;
        }
        let k = k as usize;
        let x2 = inst.xs[k * inst.n_x + inst.xm[x][m] as usize] as usize;
        let p2 = inst.fiber_perm(m, k)[p] as usize;
        self.index(c2 as usize, x2, p2)
    }
    fn equivalents(&self, j: usize, out: &mut Vec<usize>) {
        let inst = self.inst;
        let (_, c, x, p) = self.node(j, &mut ());
        for &k in &inst.stab[c] {
            let k = k as usize;
            out.push(self.index(
                c,
                inst.xs[k * inst.n_x + x] as usize,
                inst.ps[p * N_SYM + k] as usize,
            ));
        }
    }
}

/// The slot table of an instance, in memory: its histogram is the slots per depth.
pub fn slot_bfs(inst: &Instance, verbose: bool) -> Table2 {
    let g = SlotGraph { inst };
    let (c0, x0, p0) = inst.solved;
    let t = g.index(c0 as usize, x0 as usize, p0 as usize);
    let mut targets = vec![t];
    g.equivalents(t, &mut targets);
    targets.sort_unstable();
    targets.dedup();
    bfs(
        &g,
        inst.n_slots() as usize,
        &targets,
        "p2bfs-slots",
        verbose,
    )
}

/// Positions at each depth up to `max_depth`, by a plain search in the fixed-corner frame (no
/// symmetry): configuration, corner id and B-centre colouring rank.
pub fn position_bfs(max_depth: usize) -> Vec<u64> {
    let s = space();
    let ms = mset();
    let key = |p: &P2| -> u64 {
        let mut l = [0u8; 12];
        l.copy_from_slice(&p.xb[..12]);
        (u64::from(p.cfg) * 960 + u64::from(p.c)) * 369_600 + u64::from(ms.rank(&l))
    };
    let Some(start) = P2::of(&Cube::solved()) else {
        return vec![];
    };
    let mut prev: HashSet<u64> = HashSet::new();
    let mut cur: HashSet<u64> = HashSet::from([key(&start)]);
    let mut front = vec![start];
    let mut counts = vec![1u64];
    for _ in 0..max_depth {
        let mut next_set: HashSet<u64> = HashSet::new();
        let mut next = vec![];
        for p in &front {
            for m in 0..N_MOVES {
                if s.next[p.cfg as usize][m] == NONE {
                    continue;
                }
                let Some(q) = p.apply(s, m) else { continue };
                let k = key(&q);
                if !cur.contains(&k) && !prev.contains(&k) && next_set.insert(k) {
                    next.push(q);
                }
            }
        }
        counts.push(next.len() as u64);
        prev = cur;
        cur = next_set;
        front = next;
    }
    counts
}
