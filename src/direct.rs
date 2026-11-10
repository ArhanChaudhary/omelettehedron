//! Direct lower bounds on the distance to the *solved* state.
//!
//! The phase-1 tables bound the distance to H_A (or, conjugated, to H_B). Near H_A those
//! bounds are small even when the phase-2 part still needs many moves, so a proof of
//! optimality expands many nodes close to H_A that cannot lead to a shorter solution. The
//! tables here bound the distance to solved directly from the edges and from the corners.

use crate::cube::*;
use crate::fast::*;
use crate::tables::prefetch;
use rayon::prelude::*;
use std::sync::atomic::{AtomicU64, AtomicU8, Ordering};

/// Even permutations of the 12 edges (the reachable ones).
pub const N_EPERM: usize = 239_500_800;

/// Rank of an even permutation of 0..12 (lanes 0..12 of `p`) in 0..12!/2: the first ten
/// Lehmer digits in mixed radix 12, 11, ..., 3 (the last two digits follow from the parity).
#[inline(always)]
pub fn eperm_rank(p: &[u8; 16]) -> u32 {
    let mut used = 0u32;
    let mut r = 0u32;
    for (i, &v) in p.iter().take(10).enumerate() {
        let v = v as u32;
        let d = v - (used & ((1 << v) - 1)).count_ones();
        r = r * (12 - i as u32) + d;
        used |= 1 << v;
    }
    r
}

/// Inverse of [`eperm_rank`].
pub fn eperm_unrank(mut r: u32) -> [u8; 16] {
    let mut d = [0u32; 10];
    for i in (0..10).rev() {
        d[i] = r % (12 - i as u32);
        r /= 12 - i as u32;
    }
    let mut avail: u32 = 0xfff;
    let mut p = [0u8; 16];
    let mut par = 0;
    for i in 0..10 {
        // the d-th smallest unused value
        let mut a = avail;
        for _ in 0..d[i] {
            a &= a - 1;
        }
        let v = a.trailing_zeros();
        p[i] = v as u8;
        avail &= !(1 << v);
        par ^= d[i] & 1;
    }
    let lo = avail.trailing_zeros() as u8;
    let hi = 31 - avail.leading_zeros() as u8;
    // last two Lehmer digits (par, 0) make the permutation even
    if par == 0 {
        p[10] = lo;
        p[11] = hi;
    } else {
        p[10] = hi;
        p[11] = lo;
    }
    p
}

pub struct DirectTables {
    /// exact distance to solved of the edge permutation, one byte per even permutation
    pub edges: Vec<u8>,
    /// exact distance to solved of the corner state (padded to 1024 ids)
    pub corners: Box<[u8; 1024]>,
    pub edge_hist: Vec<u64>,
}

impl DirectTables {
    pub fn build(verbose: bool) -> DirectTables {
        let t0 = std::time::Instant::now();
        let d = defs();
        let gathers: [[u8; 16]; N_MOVES] = std::array::from_fn(|m| {
            let mut g = [0x80u8; 16];
            for p in 0..12 {
                g[p] = d.moves[m].ep[p];
            }
            g
        });
        let dist: Vec<AtomicU8> = (0..N_EPERM).into_par_iter().map(|_| AtomicU8::new(u8::MAX)).collect();
        let mut id = [0u8; 16];
        for p in 0..12 {
            id[p] = p as u8;
        }
        dist[eperm_rank(&id) as usize].store(0, Ordering::Relaxed);
        let mut hist = vec![1u64];
        let mut depth = 0u8;
        loop {
            let newc = AtomicU64::new(0);
            (0..N_EPERM).into_par_iter().chunks(1 << 16).for_each(|chunk| {
                let mut c = 0u64;
                for i in chunk {
                    if dist[i].load(Ordering::Relaxed) != depth {
                        continue;
                    }
                    let p = eperm_unrank(i as u32);
                    for g in &gathers {
                        let q = shuffle(&p, g);
                        let j = eperm_rank(&q) as usize;
                        if dist[j].load(Ordering::Relaxed) == u8::MAX
                            && dist[j]
                                .compare_exchange(u8::MAX, depth + 1, Ordering::Relaxed, Ordering::Relaxed)
                                .is_ok()
                        {
                            c += 1;
                        }
                    }
                }
                newc.fetch_add(c, Ordering::Relaxed);
            });
            let c = newc.into_inner();
            if c == 0 {
                break;
            }
            hist.push(c);
            depth += 1;
        }
        let total: u64 = hist.iter().sum();
        assert_eq!(total as usize, N_EPERM, "edge permutations not all reached");
        let edges: Vec<u8> = dist.into_par_iter().map(|a| a.into_inner()).collect();
        // corners: BFS over the 960 corner ids
        let f = fast();
        let mut corners = Box::new([u8::MAX; 1024]);
        let c0 = crate::coord::corners().id(&Cube::solved()) as usize;
        corners[c0] = 0;
        let mut frontier = vec![c0];
        let mut cd = 0u8;
        while !frontier.is_empty() {
            let mut next = vec![];
            for &c in &frontier {
                for m in 0..N_MOVES {
                    let y = f.corner_move(c as u16, m) as usize;
                    if corners[y] == u8::MAX {
                        corners[y] = cd + 1;
                        next.push(y);
                    }
                }
            }
            frontier = next;
            cd += 1;
        }
        if verbose {
            let mean = hist.iter().enumerate().map(|(i, &c)| i as f64 * c as f64).sum::<f64>() / total as f64;
            eprintln!(
                "direct: edges max {} mean {mean:.3} hist {:?}; corners max {}; {:.1}s",
                hist.len() - 1,
                hist,
                cd - 1,
                t0.elapsed().as_secs_f64()
            );
        }
        DirectTables {
            edges,
            corners,
            edge_hist: hist,
        }
    }

    #[inline(always)]
    pub fn index(x: &Fc) -> usize {
        eperm_rank(&x.ep) as usize
    }

    #[inline(always)]
    pub fn prefetch(&self, i: usize) {
        prefetch(&self.edges[i]);
    }

    /// Lower bound on the distance to solved (`i` = [`DirectTables::index`] of `x`).
    #[inline(always)]
    pub fn bound(&self, x: &Fc, i: usize) -> u8 {
        self.edges[i].max(self.corners[x.c as usize & 1023])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rank_roundtrip() {
        use rand::{Rng, SeedableRng};
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(5);
        for _ in 0..10000 {
            let r = rng.gen_range(0..N_EPERM as u32);
            let p = eperm_unrank(r);
            assert_eq!(perm_parity(&p[..12]), 0);
            assert_eq!(eperm_rank(&p), r);
        }
        assert_eq!(eperm_rank(&eperm_unrank(N_EPERM as u32 - 1)), N_EPERM as u32 - 1);
    }
}
