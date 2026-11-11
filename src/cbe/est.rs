//! Estimate (experiment, feature `est`): what a table with the CBE index but distances to the
//! SOLVED projection would prune in the two-phase optimal mode, without building it.
//!
//! A breadth-first ball of radius `R` around the solved projection (in the reduced CBE index
//! space, with equivalent indices), and above that radius a depth-first search over canonical
//! move sequences decides whether a projection's distance is at most a threshold.

#![allow(clippy::unwrap_used, clippy::expect_used, clippy::indexing_slicing)]

use super::*;
use crate::prune::Graph;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::OnceLock;

pub struct Ball {
    /// open addressing: 0 = empty, else (index + 1) << 4 | distance
    slots: Vec<u64>,
    shift: u32,
    pub r: u8,
    pub sizes: Vec<usize>,
    /// optional bit set of the entries within `r + 1` (empty if not built)
    pub bits: Vec<u64>,
}

impl Ball {
    #[inline]
    fn h(&self, i: usize) -> usize {
        ((i as u64 + 1).wrapping_mul(0x9E37_79B9_7F4A_7C15) >> self.shift) as usize
    }

    #[inline]
    pub fn get(&self, i: usize) -> Option<u8> {
        let mask = self.slots.len() - 1;
        let key = (i as u64 + 1) << 4;
        let mut p = self.h(i);
        loop {
            let v = self.slots[p];
            if v == 0 {
                return None;
            }
            if v & !15 == key {
                return Some((v & 15) as u8);
            }
            p = (p + 1) & mask;
        }
    }

    fn insert(&mut self, i: usize, d: u8) {
        let mask = self.slots.len() - 1;
        let key = (i as u64 + 1) << 4;
        let mut p = self.h(i);
        loop {
            let v = self.slots[p];
            if v == 0 {
                self.slots[p] = key | u64::from(d);
                return;
            }
            if v & !15 == key {
                return;
            }
            p = (p + 1) & mask;
        }
    }

    /// The ball of radius `r` around the solved projection.
    pub fn build(ix: &CbeIndex, r: u8) -> Ball {
        let f = fast();
        let g = CbeGraph::<N_MOVES>::new(ix);
        let solved = f.of(&Cube::solved());
        let i0 = ix.index(f, &solved, &masks(&solved.e), 0);
        let mut l0 = vec![i0];
        ix.equivalents(f, i0, &mut l0);
        l0.sort_unstable();
        l0.dedup();
        let mut levels: Vec<Vec<usize>> = vec![l0];
        for d in 1..=usize::from(r) {
            let t0 = std::time::Instant::now();
            let cur = &levels[d - 1];
            let prev = if d >= 2 { Some(&levels[d - 2]) } else { None };
            let mut next: Vec<usize> = cur
                .par_chunks(1 << 12)
                .flat_map_iter(|ch| {
                    let mut cache = None;
                    let mut out = vec![];
                    for &i in ch {
                        let nd = g.node(i, &mut cache);
                        for m in 0..N_MOVES {
                            let j = g.step(&nd, m);
                            if cur.binary_search(&j).is_err()
                                && prev.is_none_or(|p| p.binary_search(&j).is_err())
                            {
                                out.push(j);
                            }
                        }
                    }
                    out
                })
                .collect();
            next.par_sort_unstable();
            next.dedup();
            let mut eq: Vec<usize> = next
                .par_iter()
                .flat_map_iter(|&j| {
                    let mut v = vec![];
                    ix.equivalents(f, j, &mut v);
                    v
                })
                .collect();
            next.append(&mut eq);
            next.par_sort_unstable();
            next.dedup();
            eprintln!(
                "est ball: depth {d}: {} entries ({:.1} s)",
                next.len(),
                t0.elapsed().as_secs_f64()
            );
            levels.push(next);
        }
        let total: usize = levels.iter().map(Vec::len).sum();
        let cap = (total * 2).next_power_of_two().max(1024);
        let mut b = Ball {
            slots: vec![0; cap],
            shift: 64 - cap.trailing_zeros(),
            r,
            sizes: levels.iter().map(Vec::len).collect(),
            bits: vec![],
        };
        for (d, l) in levels.iter().enumerate() {
            for &i in l {
                b.insert(i, d as u8);
            }
        }
        if std::env::var_os("FTO_EST_BITS").is_some() {
            let t0 = std::time::Instant::now();
            let words = ix.size().div_ceil(64);
            let bits: Vec<AtomicU64> = (0..words).map(|_| AtomicU64::new(0)).collect();
            let set = |j: usize| -> bool {
                let old = bits[j >> 6].fetch_or(1 << (j & 63), Ordering::Relaxed);
                old >> (j & 63) & 1 == 0
            };
            levels.par_iter().for_each(|l| {
                l.par_iter().for_each(|&i| {
                    set(i);
                });
            });
            let last = &levels[levels.len() - 1];
            let new: u64 = last
                .par_chunks(1 << 12)
                .map(|ch| {
                    let mut cache = None;
                    let mut cnt = 0u64;
                    let mut eq = vec![];
                    for &i in ch {
                        let nd = g.node(i, &mut cache);
                        for m in 0..N_MOVES {
                            let j = g.step(&nd, m);
                            if set(j) {
                                cnt += 1;
                                eq.clear();
                                ix.equivalents(f, j, &mut eq);
                                for &e in &eq {
                                    cnt += u64::from(set(e));
                                }
                            }
                        }
                    }
                    cnt
                })
                .sum();
            eprintln!(
                "est ball: depth {} as a bit set: {new} entries ({:.1} s)",
                r + 1,
                t0.elapsed().as_secs_f64()
            );
            b.bits = bits.into_iter().map(AtomicU64::into_inner).collect();
        }
        b
    }

    /// Whether entry `i` is within `r + 1` (needs `bits`).
    #[inline]
    pub fn within_r1(&self, i: usize) -> bool {
        self.bits.get(i >> 6).is_some_and(|w| w >> (i & 63) & 1 == 1)
    }
}

pub struct Est {
    pub ix: CbeIndex,
    pub ball: Ball,
    /// thresholds above this are not decided (the child is kept)
    pub tmax: u8,
    /// the saturating codes' bases of the candidates
    pub bases: Vec<u8>,
    pub dfs_nodes: AtomicU64,
}

/// Candidates: bit 0 none (sanity), 1 exact (8 lookups), 2 exact (side 0's 4 lookups),
/// 3.. saturating with base `bases[c - 3]` (8 lookups).
pub const MAX_C: usize = 8;

#[derive(Default, Clone)]
pub struct Totals {
    pub nodes: [u64; MAX_C],
    /// per candidate and moves left L: (checked, pruned)
    pub by_l: [[(u64, u64); 24]; MAX_C],
    /// children whose `sat9` decision was compared with the to-solved table's, and mismatches
    pub compared: u64,
    pub mismatch: u64,
}

impl Totals {
    pub fn add(&mut self, o: &Totals) {
        for c in 0..MAX_C {
            self.nodes[c] += o.nodes[c];
            for l in 0..24 {
                self.by_l[c][l].0 += o.by_l[c][l].0;
                self.by_l[c][l].1 += o.by_l[c][l].1;
            }
        }
        self.compared += o.compared;
        self.mismatch += o.mismatch;
    }
}

/// Per-worker state: candidate masks per depth, and counts.
#[derive(Default)]
pub struct EstCtr {
    pub mask: Vec<u8>,
    pub t: Totals,
}

impl EstCtr {
    pub fn set(&mut self, g: usize, m: u8) {
        if self.mask.len() <= g {
            self.mask.resize(g + 1, 0);
        }
        self.mask[g] = m;
    }
    pub fn get(&self, g: usize) -> u8 {
        self.mask.get(g).copied().unwrap_or(0)
    }
}

pub fn est() -> Option<&'static Est> {
    static E: OnceLock<Option<Est>> = OnceLock::new();
    E.get_or_init(|| {
        if std::env::var_os("FTO_EST").is_none() {
            return None;
        }
        let env = |n: &str, d: u8| {
            std::env::var(n)
                .ok()
                .and_then(|s| s.parse().ok())
                .unwrap_or(d)
        };
        let r = env("FTO_EST_R", 8);
        let ix = CbeIndex::new();
        let t0 = std::time::Instant::now();
        let threads = std::env::var("FTO_EST_THREADS")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(16);
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("pool");
        let ball = pool.install(|| Ball::build(&ix, r));
        eprintln!(
            "est ball: radius {r}, {} entries, {:.1} s",
            ball.sizes.iter().sum::<usize>(),
            t0.elapsed().as_secs_f64()
        );
        Some(Est {
            ix,
            ball,
            tmax: env("FTO_EST_TMAX", 11),
            bases: vec![8, 9, 10, 11],
            dfs_nodes: AtomicU64::new(0),
        })
    })
    .as_ref()
}

impl Est {
    pub fn all_mask(&self) -> u8 {
        ((1u32 << (3 + self.bases.len())) - 1) as u8
    }

    /// Whether the projection of triple `k` of `y` is within `t` moves of solved.
    pub fn within(&self, y: &Fc, k: usize, t: u8, last: usize, nodes: &mut u64) -> bool {
        *nodes += 1;
        let f = fast();
        let em = masks(&y.e);
        let i = self.ix.index(f, y, &em, k);
        if let Some(d) = self.ball.get(i) {
            return d <= t;
        }
        if t <= self.ball.r {
            return false;
        }
        if !self.ball.bits.is_empty() {
            if self.ball.within_r1(i) {
                // (exactly r + 1)
                return true;
            }
            if t <= self.ball.r + 1 {
                return false;
            }
        }
        for m in 0..N_MOVES {
            if !allowed_after(false, last, m) {
                continue;
            }
            let z = f.apply(y, &em, m);
            if self.within(&z, k, t - 1, m, nodes) {
                return true;
            }
        }
        false
    }

    /// Whether some of the eight projections (four triples of each position) is more than `t`
    /// moves from solved.
    pub fn any_beyond(&self, side0: &Fc, side1: &Fc, t: u8) -> bool {
        let mut nodes = 0u64;
        (0..4).any(|k| !self.within(side0, k, t, usize::MAX, &mut nodes))
            || (0..4).any(|k| !self.within(side1, k, t, usize::MAX, &mut nodes))
    }

    /// The candidates (of `mask`) that prune a child with positions `side0`, `side1` and
    /// `l` moves left for the whole solution.
    pub fn prune_mask(&self, mask: u8, side0: &Fc, side1: &Fc, l: usize, ctr: &mut EstCtr) -> u8 {
        let mut nodes = 0u64;
        // per threshold: whether some projection exceeds it (4: side 0 only, 8: both)
        let mut c4: [Option<bool>; 16] = [None; 16];
        let mut c8: [Option<bool>; 16] = [None; 16];
        let lc = l.min(23);
        // the threshold of each candidate (prune iff some projection's distance exceeds it)
        let mut ts = [usize::MAX; MAX_C];
        for c in 1..3 + self.bases.len() {
            if mask >> c & 1 == 0 {
                continue;
            }
            ctr.t.by_l[c][lc].0 += 1;
            ts[c] = if c <= 2 {
                l
            } else {
                let b = usize::from(self.bases[c - 3]);
                if l > b + 2 {
                    usize::MAX
                } else {
                    l.max(b)
                }
            };
        }
        // decide the largest thresholds first: exceeding t implies exceeding every t' < t
        let mut order: Vec<usize> = (1..3 + self.bases.len()).filter(|&c| ts[c] <= usize::from(self.tmax)).collect();
        order.sort_by_key(|&c| std::cmp::Reverse(ts[c]));
        let mut out = 0u8;
        for c in order {
            let t = ts[c] as u8;
            let four = c == 2;
            let known = |cache: &[Option<bool>; 16]| -> Option<bool> {
                if let Some(v) = cache[usize::from(t)] {
                    return Some(v);
                }
                // a larger threshold exceeded
                if (usize::from(t) + 1..16).any(|u| cache[u] == Some(true)) {
                    return Some(true);
                }
                // a smaller threshold not exceeded
                if (0..usize::from(t)).any(|u| cache[u] == Some(false)) {
                    return Some(false);
                }
                None
            };
            let p4 = |nodes: &mut u64, c4: &mut [Option<bool>; 16]| -> bool {
                if let Some(v) = known(c4) {
                    return v;
                }
                let v = (0..4).any(|k| !self.within(side0, k, t, usize::MAX, nodes));
                c4[usize::from(t)] = Some(v);
                v
            };
            let p = if four {
                p4(&mut nodes, &mut c4)
            } else if let Some(v) = known(&c8) {
                v
            } else {
                let v = p4(&mut nodes, &mut c4)
                    || (0..4).any(|k| !self.within(side1, k, t, usize::MAX, &mut nodes));
                c8[usize::from(t)] = Some(v);
                v
            };
            if p {
                out |= 1 << c;
                ctr.t.by_l[c][lc].1 += 1;
            }
        }
        self.dfs_nodes.fetch_add(nodes, Ordering::Relaxed);
        out
    }

    pub fn report(&self, t: &Totals) {
        let names: Vec<String> = ["none".to_string(), "exact8".into(), "exact4".into()]
            .into_iter()
            .chain(self.bases.iter().map(|b| format!("sat{b}")))
            .collect();
        let base = t.nodes[0].max(1) as f64;
        let mut s = String::from("est:");
        for (c, n) in names.iter().enumerate() {
            s += &format!(" {n}={} ({:.4})", t.nodes[c], t.nodes[c] as f64 / base);
        }
        s += &format!(" compared={} mismatch={}", t.compared, t.mismatch);
        s += "\n";
        for (c, n) in names.iter().enumerate().skip(1) {
            s += &format!("est-L {n}:");
            for l in 0..24 {
                let (a, b) = t.by_l[c][l];
                if a > 0 {
                    s += &format!(" {l}:{a}/{b}");
                }
            }
            s += "\n";
        }
        eprint!("{s}");

    }
}
