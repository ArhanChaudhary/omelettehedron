//! Index of the corner/B-centre (CB) table reduced by the 24 class-preserving symmetries.
//!
//! The CB coordinate (corner state x B-centre colouring, 3.5e8 entries, 89 MB at 2 bits) is
//! conjugated by the symmetry that maps its corner state to the smallest one of its class, and
//! the front corner is re-fixed with a whole-puzzle rotation (the distance to `H_A` is invariant
//! under both), as the joint table does for its corners. That leaves 51 corner classes: 1.9e7
//! entries, 4.7 MB, which stays in the cache. (+1.3% / +1.0% nodes/s at 48T / 1T measured by the
//! round-3 analysis, archive/analysis-round3.)

use crate::cbe::conj_renorm_corner_with;
use crate::coord::*;
use crate::cube::*;
use crate::fast::*;
use crate::prune::Graph;
use crate::sym::*;
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::OnceLock;

pub struct CbSymIndex {
    /// per corner id (padded to 1024): (class id, symmetry M, reposition R*M)
    kt: Box<[(u16, u8, u8); 1024]>,
    /// class id -> representative corner id
    rep: Vec<u16>,
    /// per class: (M, R*M) pairs other than the identity that fix the representative
    stab: Vec<Vec<(u8, u8)>>,
    /// per symmetry (padded to 32): B colour relabelling and B position gather
    bmap: [[u8; 16]; 32],
    bgather: [[u8; 16]; 32],
}

/// The index, built once (about 0.1 s).
pub fn cbsym() -> &'static CbSymIndex {
    static IX: OnceLock<CbSymIndex> = OnceLock::new();
    IX.get_or_init(CbSymIndex::new)
}

impl CbSymIndex {
    fn new() -> CbSymIndex {
        let s = syms();
        let wps: Vec<Vec<usize>> = s.mats.iter().map(|r| defs().sm.whole(r)).collect();
        // conj[m][c]: corner id and re-fixing rotation of corner state c conjugated by m (a
        // derivation that cannot fail for the fixed geometry; identity as the fallback)
        let conj: Vec<Vec<(u16, usize)>> = (0..N_GEO)
            .into_par_iter()
            .map(|m| {
                (0..N_CORNER as u16)
                    .map(|c| {
                        let r = conj_renorm_corner_with(c, m, &wps);
                        debug_assert!(r.is_some(), "no rotation re-fixes the front corner");
                        r.unwrap_or((c, 0))
                    })
                    .collect()
            })
            .collect();
        let min_image = |c: usize| {
            (0..N_GEO)
                .map(|m| (conj[m][c].0, m))
                .min()
                .unwrap_or((c as u16, 0))
        };
        let mut rep: Vec<u16> = (0..N_CORNER).map(|c| min_image(c).0).collect();
        rep.sort_unstable();
        rep.dedup();
        let cls_id: HashMap<u16, u16> = rep
            .iter()
            .enumerate()
            .map(|(i, &r)| (r, i as u16))
            .collect();
        let mut kt = Box::new([(0u16, 0u8, 0u8); 1024]);
        for (c, k) in kt.iter_mut().enumerate().take(N_CORNER) {
            let (c2, m) = min_image(c);
            let r = conj[m][c].1;
            *k = (cls_id.get(&c2).copied().unwrap_or(0), m as u8, s.mul[r][m]);
        }
        let stab = rep
            .iter()
            .map(|&r| {
                (1..N_GEO)
                    .filter(|&m| conj[m][r as usize].0 == r)
                    .map(|m| (m as u8, s.mul[conj[m][r as usize].1][m]))
                    .collect()
            })
            .collect();
        let bmap = std::array::from_fn(|m| {
            let mut v = [0u8; 16];
            if m < N_GEO {
                v[..4].copy_from_slice(&s.b_col[m]);
            }
            v
        });
        let bgather = std::array::from_fn(|t| {
            let mut g = [0x80u8; 16];
            if t < N_GEO {
                for p in 0..12 {
                    g[s.b_pos[t][p] as usize] = p as u8;
                }
            }
            g
        });
        CbSymIndex {
            kt,
            rep,
            stab,
            bmap,
            bgather,
        }
    }

    /// Number of entries.
    pub fn size(&self) -> usize {
        self.rep.len() * N_MSET
    }

    /// Number of corner classes.
    pub fn classes(&self) -> usize {
        self.rep.len()
    }

    /// Rank of the B colouring `b` after symmetry `m` and reposition `t`.
    #[inline(always)]
    fn row(&self, f: &FastTables, b: &[u8; 16], m: usize, t: usize) -> usize {
        let b = shuffle(&shuffle(&self.bmap[m & 31], b), &self.bgather[t & 31]);
        f.rank(&masks(&b)) as usize
    }

    /// Index of a position (only its corners and B-centres matter).
    #[inline(always)]
    pub fn index(&self, f: &FastTables, x: &Fc) -> usize {
        let (cls, m, t) = self.kt[x.c as usize & 1023];
        cls as usize * N_MSET + self.row(f, &x.b, m as usize, t as usize)
    }

    /// A position with index `i` (corners and B-centres only).
    pub fn decode(&self, i: usize) -> Fc {
        let mut x = Fc {
            a: [0; 16],
            b: [0; 16],
            e: [0; 16],
            c: self.rep.get(i / N_MSET).copied().unwrap_or(0),
            ebits: 0,
        };
        x.b[..12].copy_from_slice(&mset().unrank((i % N_MSET) as u32));
        x
    }

    /// The other indices of the position with index `i` (its class representative has a
    /// non-trivial stabiliser).
    pub fn equivalents(&self, i: usize, out: &mut Vec<usize>) {
        let cls = i / N_MSET;
        let Some(stab) = self.stab.get(cls) else {
            return;
        };
        if stab.is_empty() {
            return;
        }
        let f = fast();
        let x = self.decode(i);
        for &(m, t) in stab {
            out.push(cls * N_MSET + self.row(f, &x.b, m as usize, t as usize));
        }
    }
}

/// Moves between reduced CB entries (for table generation), `N` moves (16, or 24 with the
/// slices: the reduction stays valid, as every symmetry maps slices to slices).
pub struct CbSymGraph<const N: usize>;

impl<const N: usize> Graph for CbSymGraph<N> {
    const MOVES: usize = N;
    type Node = Fc;
    type Cache = ();
    fn node(&self, i: usize, (): &mut ()) -> Fc {
        cbsym().decode(i)
    }
    fn step(&self, x: &Fc, m: usize) -> usize {
        let f = fast();
        // (the edge labels are all 0, so no chirality flips are needed)
        cbsym().index(f, &f.apply_m(x, &[0xfff, 0, 0, 0], m))
    }
    fn equivalents(&self, j: usize, out: &mut Vec<usize>) {
        cbsym().equivalents(j, out);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    /// The reduced index depends only on the corners and B-centres, and a decoded entry has
    /// the index it was decoded from.
    #[test]
    fn index_depends_on_cb_only() {
        let ix = cbsym();
        assert_eq!(ix.classes(), 51);
        let f = fast();
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(9);
        for _ in 0..2000 {
            let c = Cube::random(&mut rng);
            let x = f.of(&c);
            let i = ix.index(f, &x);
            // same corners and B-centres, other A-centres and edges: same index
            let mut c2 = Cube::random(&mut rng);
            c2.cp = c.cp;
            c2.co = c.co;
            c2.xb = c.xb;
            assert_eq!(ix.index(f, &f.of(&c2)), i);
            // the decoded representative has the same index
            assert_eq!(ix.index(f, &ix.decode(i)), i);
        }
    }
}
