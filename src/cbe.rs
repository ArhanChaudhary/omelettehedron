//! Joint corner x B-centre x edge-triple table ("CBE").
//!
//! The two main phase-1 tables see the A-centres/edges (AE) and the corners/B-centres (CB)
//! separately. This table sees the corners, the B-centre colouring and the positions and
//! chirality of the three edges of one A-face together, which bounds the middle of the
//! search much better than max(AE, CB).
//!
//! Triple k is mapped onto triple 0 by conjugating with a class-preserving symmetry M with
//! M(k) = 0 and then re-fixing the front corner with a whole-puzzle rotation R (the distance to
//! `H_A` is invariant under both). The corner coordinate is reduced by the 6 symmetries fixing
//! triple 0. Values are stored as a saturating code: 0 = at most 6 moves, 1 = 7, 2 = 8,
//! 3 = at least 9; only the 7..9 range adds information over the other tables.
//!
//! (2.05x fewer nodes; first solution 0.75x / 0.65x at 48T / 1T when added, opt/cbe.)

use crate::coord::*;
use crate::cube::*;
use crate::fast::*;
use crate::geometry::*;
use crate::mem::HugeBuf;
use crate::prune::*;
use crate::sym::*;
use crate::tables::{prefetch, CbTable, Phase2};
use rayon::prelude::*;
use std::collections::HashMap;
use std::path::Path;

/// 220 position sets x 1 chirality bit
pub const N_E1: usize = 440;
pub const N_ROW: usize = N_MSET * N_E1;

/// Lower bound encoded by each saturating code.
pub const CODE_LB: [u8; 4] = [0, 7, 8, 9];

pub struct CbeIndex {
    /// per corner id (padded to 1024) and triple k: (class id, symmetry M, reposition R*M)
    kt: Box<[[(u16, u8, u8); 4]; 1024]>,
    /// class id -> representative corner id
    pub rep: Vec<u16>,
    /// per class: (M, R*M) pairs other than the identity that fix the representative
    stab: Vec<Vec<(u8, u8)>>,
    /// (per symmetry, padded to 32 so that `& 31` replaces bounds checks)
    bmap: [[u8; 16]; 32],
    bgather: [[u8; 16]; 32],
    egather: [[u8; 16]; 32],
    relab_flip: [[u8; 4]; 32],
    c12_unrank: Vec<u16>,
}

/// Conjugate the corner state `c` by symmetry `m` and re-fix the front corner; returns the new
/// corner id and the index of the rotation used.
#[cfg(test)]
fn conj_renorm_corner(c: u16, m: usize) -> Option<(u16, usize)> {
    let wps: Vec<Vec<usize>> = syms().mats.iter().map(|r| defs().sm.whole(r)).collect();
    conj_renorm_corner_with(c, m, &wps)
}

/// [`conj_renorm_corner`] with the sticker permutations of the symmetries precomputed. `None`
/// if no rotation re-fixes the front corner (never: the rotations act transitively on the
/// vertices).
fn conj_renorm_corner_with(c: u16, m: usize, wps: &[Vec<usize>]) -> Option<(u16, usize)> {
    let s = syms();
    let (cp, co) = corners().id_to_state[c as usize];
    let st = Cube {
        cp,
        co,
        ..Cube::solved()
    }
    .to_stickers();
    let mm = &s.mats[m];
    let wp = &wps[m];
    let mut n = [0u8; 72];
    for i in 0..72 {
        n[wp[i]] = face_of_normal(mat_vec(mm, FACES[st[i] as usize])) as u8;
    }
    for (ri, r) in s.mats.iter().enumerate() {
        if det(r) != 1 {
            continue;
        }
        let rp = &wps[ri];
        let mut n2 = [0u8; 72];
        for i in 0..72 {
            n2[rp[i]] = n[i];
        }
        if let Some(cube) = Cube::from_stickers(&n2) {
            if cube.cp[0] == 0 && cube.co[0] == 0 {
                return Some((corners().id(&cube), ri));
            }
        }
    }
    None
}

impl CbeIndex {
    pub fn new() -> CbeIndex {
        let s = syms();
        // symmetries mapping A-face k to A-face 0
        let to0 = |k: usize| -> Vec<usize> { (0..N_GEO).filter(|&m| s.a_col[m][k] == 0).collect() };
        let s0 = to0(0);
        let wps: Vec<Vec<usize>> = s.mats.iter().map(|r| defs().sm.whole(r)).collect();
        let conj: Vec<Vec<(u16, usize)>> = (0..N_GEO)
            .into_par_iter()
            .map(|m| {
                (0..N_CORNER as u16)
                    .map(|c| {
                        let r = conj_renorm_corner_with(c, m, &wps);
                        debug_assert!(r.is_some(), "corner {c}, symmetry {m}");
                        r.unwrap_or((c, 0))
                    })
                    .collect()
            })
            .collect();
        // class representatives under the stabiliser of triple 0 (s0 holds the identity)
        let mut rep_of = vec![0u16; N_CORNER];
        for c in 0..N_CORNER {
            rep_of[c] = s0.iter().map(|&m| conj[m][c].0).fold(u16::MAX, u16::min);
        }
        let mut rep: Vec<u16> = rep_of.clone();
        rep.sort_unstable();
        rep.dedup();
        let cls_id: HashMap<u16, u16> = rep
            .iter()
            .enumerate()
            .map(|(i, &r)| (r, i as u16))
            .collect();
        let kt4: [Vec<(u16, u8, u8)>; 4] = std::array::from_fn(|k| {
            (0..N_CORNER)
                .map(|c| {
                    // (some symmetry maps A-face k to A-face 0)
                    let best = to0(k)
                        .into_iter()
                        .map(|m| (m, conj[m][c]))
                        .min_by_key(|&(_, (c2, _))| c2);
                    debug_assert!(best.is_some());
                    let (m, (c2, r)) = best.unwrap_or((0, conj[0][c]));
                    (cls_id[&c2], m as u8, s.mul[r][m])
                })
                .collect()
        });
        let stab: Vec<Vec<(u8, u8)>> = rep
            .iter()
            .map(|&r| {
                s0.iter()
                    .filter(|&&m| m != 0 && conj[m][r as usize].0 == r)
                    .map(|&m| (m as u8, s.mul[conj[m][r as usize].1][m]))
                    .collect()
            })
            .collect();
        let mut kt = Box::new([[(0u16, 0u8, 0u8); 4]; 1024]);
        for c in 0..N_CORNER {
            for k in 0..4 {
                kt[c][k] = kt4[k][c];
            }
        }
        let pad_map = |m: usize| {
            let mut v = [0u8; 16];
            if m < N_GEO {
                v[..4].copy_from_slice(&s.b_col[m]);
            }
            v
        };
        let gather = |pos: &[[u8; 12]; N_GEO], t: usize| {
            let mut g = [0x80u8; 16];
            if t < N_GEO {
                for p in 0..12 {
                    g[pos[t][p] as usize] = p as u8;
                }
            }
            g
        };
        let mut c12_unrank = vec![0u16; 220];
        let f = fast();
        for m in 0..4096usize {
            if m.count_ones() == 3 {
                c12_unrank[f.c12[m] as usize] = m as u16;
            }
        }
        CbeIndex {
            kt,
            rep,
            stab,
            bmap: std::array::from_fn(pad_map),
            bgather: std::array::from_fn(|t| gather(&s.b_pos, t)),
            egather: std::array::from_fn(|t| gather(&s.e_pos, t)),
            relab_flip: std::array::from_fn(|m| {
                if m < N_GEO {
                    std::array::from_fn(|k| (s.relab_bits[m][0] >> s.a_col[m][k]) & 1)
                } else {
                    [0; 4]
                }
            }),
            c12_unrank,
        }
    }

    pub fn size(&self) -> usize {
        self.rep.len() * N_ROW
    }

    /// Row offset within a class for the (B, triple k) part after symmetry (m, t), where m
    /// maps triple k to triple 0.
    #[inline(always)]
    fn row(&self, f: &FastTables, x: &Fc, em: [u16; 4], k: usize, m: usize, t: usize) -> usize {
        let (m, t) = (m & 31, t & 31);
        let b = shuffle(&shuffle(&self.bmap[m], &x.b), &self.bgather[t]);
        // positions of triple k after repositioning (it is relabelled to 0)
        let mask = label_mask(&shuffle(&x.e, &self.egather[t]), k as u8) as usize;
        let ft = &f.flip_t[t];
        let q = em[k] as usize & 0xfff;
        let bit =
            ((x.ebits >> k) & 1) ^ self.relab_flip[m][k] ^ ((ft[q >> 6] >> (q & 63)) as u8 & 1);
        f.rank(&masks(&b)) as usize * N_E1 + f.c12[mask & 0xfff] as usize * 2 + bit as usize
    }

    /// Index of triple `k` of position `x` (`em` = edge masks of `x`).
    #[inline(always)]
    pub fn index(&self, f: &FastTables, x: &Fc, em: &[u16; 4], k: usize) -> usize {
        let (cls, m, t) = self.kt[x.c as usize & 1023][k];
        cls as usize * N_ROW + self.row(f, x, *em, k, m as usize, t as usize)
    }

    /// A position with the given index (triple 0 only; other edges labelled 1).
    fn decode(&self, i: usize) -> Fc {
        let cls = i / N_ROW;
        let r = i % N_ROW;
        let (b, e1) = (r / N_E1, r % N_E1);
        let mask = self.c12_unrank[e1 / 2];
        let mut x = Fc {
            a: [0; 16],
            b: [0; 16],
            e: [0; 16],
            c: self.rep[cls],
            ebits: (e1 % 2) as u8,
        };
        x.b[..12].copy_from_slice(&mset().unrank(b as u32));
        for p in 0..12 {
            x.e[p] = u8::from(mask >> p & 1 == 0);
        }
        x
    }

    fn equivalents(&self, f: &FastTables, i: usize, out: &mut Vec<usize>) {
        let cls = i / N_ROW;
        if self.stab[cls].is_empty() {
            return;
        }
        let x = self.decode(i);
        let em = masks(&x.e);
        for &(m, t) in &self.stab[cls] {
            out.push(cls * N_ROW + self.row(f, &x, em, 0, m as usize, t as usize));
        }
    }
}

/// Moves between CBE entries (all of triple 0).
struct CbeGraph<'a> {
    ix: &'a CbeIndex,
    f: &'static FastTables,
}

impl Graph for CbeGraph<'_> {
    const MOVES: usize = N_MOVES;
    /// position and its flip word
    type Node = (Fc, u64);
    #[inline(always)]
    fn node(&self, i: usize) -> (Fc, u64) {
        let x = self.ix.decode(i);
        (x, self.f.flip_word(&masks(&x.e)))
    }
    #[inline(always)]
    fn step(&self, (x, fw): &(Fc, u64), m: usize) -> usize {
        let y = self.f.apply_w(x, *fw, m);
        self.ix.index(self.f, &y, &masks(&y.e), 0)
    }
    fn equivalents(&self, j: usize, out: &mut Vec<usize>) {
        self.ix.equivalents(self.f, j, out);
    }
}

pub struct CbeTable {
    pub ix: CbeIndex,
    pub codes: HugeBuf<u64>,
    pub hist: Vec<u64>,
}

impl CbeTable {
    #[inline(always)]
    pub fn code(&self, i: usize) -> u8 {
        ((self.codes[i >> 5] >> ((i & 31) << 1)) & 3) as u8
    }

    #[inline(always)]
    pub fn word(&self, i: usize) -> &u64 {
        &self.codes[i >> 5]
    }

    pub fn build(verbose: bool) -> CbeTable {
        let ix = CbeIndex::new();
        let f = fast();
        let n = ix.size();
        if verbose {
            eprintln!(
                "cbe: {} corner classes, {} entries ({:.1} GB per 2-bit array)",
                ix.rep.len(),
                n,
                n as f64 / 4e9
            );
        }
        // targets: all of H_A (corners free, B attached, edges in H_A arrangements)
        let p2 = Phase2::new();
        let attach = CbTable::attach_targets();
        let edges: Vec<[u8; 12]> = p2.edge_id.keys().copied().collect();
        let mut targets: Vec<usize> = (0..N_CORNER)
            .into_par_iter()
            .flat_map_iter(|c| {
                let (cp, co) = corners().id_to_state[c];
                let xb = mset().unrank(attach[c]);
                let ix = &ix;
                edges.iter().flat_map(move |&ep| {
                    let cube = Cube {
                        cp,
                        co,
                        ep,
                        xb,
                        ..Cube::solved()
                    };
                    let x = f.of(&cube);
                    let mut v = vec![ix.index(f, &x, &masks(&x.e), 0)];
                    ix.equivalents(f, v[0], &mut v);
                    v
                })
            })
            .collect();
        targets.sort_unstable();
        targets.dedup();
        // The codes are made directly (see `saturating`).
        let t = saturating(&CbeGraph { ix: &ix, f }, n, &targets, 6, "cbe", verbose);
        CbeTable {
            ix,
            codes: t.data,
            hist: t.hist,
        }
    }

    /// # Errors
    ///
    /// If the file cannot be written.
    pub fn save(&self, path: &Path) -> std::io::Result<()> {
        save_words(path, self.ix.size(), &self.hist, &self.codes)
    }

    /// # Errors
    ///
    /// If the file cannot be read or does not hold a table of the expected size.
    pub fn load(path: &Path) -> std::io::Result<CbeTable> {
        let ix = CbeIndex::new();
        let t = Table2::load(path, ix.size())?;
        Ok(CbeTable {
            ix,
            codes: t.data,
            hist: t.hist,
        })
    }

    /// Lower bound from all four triples of `x`, prefetching first.
    #[inline(always)]
    pub fn prefetch_all(&self, f: &FastTables, x: &Fc, out: &mut [usize; 4]) {
        let em = masks(&x.e);
        for k in 0..4 {
            out[k] = self.ix.index(f, x, &em, k);
            prefetch(self.word(out[k]));
        }
    }

    #[inline(always)]
    pub fn bound(&self, idx: &[usize; 4]) -> u8 {
        idx.iter()
            .fold(0, |b, &i| b.max(CODE_LB[self.code(i) as usize]))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    /// The fast symmetry transform equals conjugation + re-fixing at sticker level.
    #[test]
    fn transform_matches_stickers() {
        let s = syms();
        let d = defs();
        let f = fast();
        let ix = CbeIndex::new();
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(41);
        for _ in 0..50 {
            let cube = Cube::random(&mut rng);
            let x = f.of(&cube);
            let em = masks(&x.e);
            for m in 0..N_GEO {
                let k = (0..4).find(|&k| s.a_col[m][k] == 0).unwrap();
                // sticker-level conjugation by m, then rotate so the front corner is fixed
                let st = cube.to_stickers();
                let mm = &s.mats[m];
                let wp = d.sm.whole(mm);
                let mut n = [0u8; 72];
                for i in 0..72 {
                    n[wp[i]] = face_of_normal(mat_vec(mm, FACES[st[i] as usize])) as u8;
                }
                let (_, r) = conj_renorm_corner(x.c, m).unwrap();
                let rp = d.sm.whole(&s.mats[r]);
                let mut n2 = [0u8; 72];
                for i in 0..72 {
                    n2[rp[i]] = n[i];
                }
                let y = Cube::from_stickers(&n2).expect("valid");
                let fy = f.of(&y);
                let t = s.mul[r][m] as usize;
                let row = ix.row(f, &x, em, k, m, t);
                let fe = masks(&fy.e)[0] as usize;
                let want = f.rank(&masks(&fy.b)) as usize * N_E1
                    + f.c12[fe] as usize * 2
                    + (fy.ebits & 1) as usize;
                assert_eq!(row, want, "m={m} k={k}");
                assert_eq!(fy.c, conj_renorm_corner(x.c, m).unwrap().0);
            }
        }
    }
}
