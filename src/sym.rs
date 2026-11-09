//! Symmetry reduction of the A-centre/edge ("AE") phase-1 coordinate.
//!
//! The distance to the phase-2 group H_A = <U, F, BR, BL> is unchanged by
//! * relabelling colours/pieces with a class-preserving symmetry `g`, and
//! * moving positions with a class-preserving symmetry `t` (whole-puzzle rotations and
//!   reflections; reflections mirror the move set),
//!
//! provided `g` and `t` are both proper or both improper. That gives 288 symmetries acting on
//! (A-centre colouring, edge triple labelling, chirality).

use crate::coord::*;
use crate::cube::*;
use crate::geometry::*;
use crate::mem::HugeBuf;
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::OnceLock;

pub const N_GEO: usize = 24;

pub struct Syms {
    pub mats: Vec<Mat>,
    /// A-colour (= A-triple label) relabelling
    pub a_col: [[u8; 4]; N_GEO],
    /// A-centre position map p -> t(p)
    pub a_pos: [[u8; 12]; N_GEO],
    /// edge position map
    pub e_pos: [[u8; 12]; N_GEO],
    /// edge piece relabelling
    pub e_piece: [[u8; 12]; N_GEO],
    /// chirality bits after relabelling with g: relab_bits[g][bits]
    pub relab_bits: [[u8; 16]; N_GEO],
    /// the 288 (g, t) pairs
    pub pairs: Vec<(u8, u8)>,
    /// relabelled triple labelling rank: relab_e[g * N_MSET + rank]
    pub relab_e: HugeBuf<u32>,
    /// repositioned triple labelling: repos_e[t * N_MSET + rank] = new rank | flip mask << 20
    pub repos_e: HugeBuf<u32>,
}

fn build() -> Syms {
    let d = defs();
    let sm = &d.sm;
    let mats = class_preserving_group();
    assert_eq!(mats.len(), N_GEO);
    let mut a_col = [[0u8; 4]; N_GEO];
    let mut a_pos = [[0u8; 12]; N_GEO];
    let mut e_pos = [[0u8; 12]; N_GEO];
    let mut e_piece = [[0u8; 12]; N_GEO];
    let edir = |e: usize| -> V3 {
        let p1 = sm.st[d.edge_slots[e][0]].p;
        let p2 = sm.st[d.edge_slots[e][1]].p;
        let f = |x: i32, y: i32| {
            if x.abs() == 4 && y.abs() == 4 {
                x.signum()
            } else {
                0
            }
        };
        [f(p1[0], p2[0]), f(p1[1], p2[1]), f(p1[2], p2[2])]
    };
    for (gi, m) in mats.iter().enumerate() {
        for i in 0..4 {
            let f = face_of_normal(mat_vec(m, FACES[A_FACES[i]]));
            a_col[gi][i] = A_FACES.iter().position(|&x| x == f).unwrap() as u8;
        }
        for p in 0..12 {
            let q = sm.index_of(mat_vec(m, sm.st[d.a_pos[p]].p));
            a_pos[gi][p] = d.a_pos.iter().position(|&x| x == q).unwrap() as u8;
            let dir = mat_vec(m, edir(p));
            e_pos[gi][p] = (0..12).find(|&e| edir(e) == dir).unwrap() as u8;
        }
        // edge pieces are named by their home position
        e_piece[gi] = e_pos[gi];
    }
    let tp = triple_pieces();
    let mut relab_bits = [[0u8; 16]; N_GEO];
    for g in 0..N_GEO {
        for bits in 0..16u8 {
            let mut nb = 0;
            for t in 0..4 {
                let [k0, k1, k2] = tp[t];
                let ep = &e_piece[g];
                let flip = parity3(ep[k0 as usize], ep[k1 as usize], ep[k2 as usize]);
                nb |= ((bits >> t & 1) ^ flip) << a_col[g][t];
            }
            relab_bits[g][bits as usize] = nb;
        }
    }
    let mut pairs = vec![];
    for g in 0..N_GEO {
        for t in 0..N_GEO {
            if det(&mats[g]) == det(&mats[t]) {
                pairs.push((g as u8, t as u8));
            }
        }
    }
    let ms = mset();
    let relab_e: Vec<Vec<u32>> = (0..N_GEO)
        .map(|g| {
            (0..N_MSET as u32)
                .into_par_iter()
                .map(|r| {
                    let l = ms.unrank(r);
                    ms.rank(&l.map(|x| a_col[g][x as usize]))
                })
                .collect()
        })
        .collect();
    let repos_e: Vec<Vec<u32>> = (0..N_GEO)
        .map(|t| {
            (0..N_MSET as u32)
                .into_par_iter()
                .map(|r| {
                    let l = ms.unrank(r);
                    let mut n = [0u8; 12];
                    for p in 0..12 {
                        n[e_pos[t][p] as usize] = l[p];
                    }
                    ms.rank(&n) | (triple_flip_mask(&l, &e_pos[t]) as u32) << 20
                })
                .collect()
        })
        .collect();
    let relab_e = HugeBuf::from_slice(&relab_e.concat());
    let repos_e = HugeBuf::from_slice(&repos_e.concat());
    Syms {
        mats,
        a_col,
        a_pos,
        e_pos,
        e_piece,
        relab_bits,
        pairs,
        relab_e,
        repos_e,
    }
}

pub fn syms() -> &'static Syms {
    static S: OnceLock<Syms> = OnceLock::new();
    S.get_or_init(build)
}

impl Syms {
    /// Apply symmetry pair (g, t) to an A-centre colouring.
    pub fn apply_a(&self, g: usize, t: usize, l: &[u8; 12]) -> [u8; 12] {
        let mut n = [0u8; 12];
        for p in 0..12 {
            n[self.a_pos[t][p] as usize] = self.a_col[g][l[p] as usize];
        }
        n
    }
    /// Apply (g, t) to the edge coordinate (triple rank, 4 chirality bits).
    #[inline(always)]
    pub fn apply_e(&self, g: usize, t: usize, e: u32, bits: u8) -> (u32, u8) {
        let e1 = self.relab_e[g * N_MSET + e as usize];
        let b1 = self.relab_bits[g][bits as usize];
        let p = self.repos_e[t * N_MSET + e1 as usize];
        (p & 0xfffff, b1 ^ (p >> 20) as u8)
    }
}

/// Symmetry classes of A-centre colourings.
pub struct AClasses {
    /// rank -> class id (low 16 bits) | symmetry pair index << 16
    pub sym: HugeBuf<u32>,
    /// class id -> representative rank
    pub rep: Vec<u32>,
}

impl AClasses {
    /// Classes under the given list of (g, t) pairs (first pair must be the identity).
    pub fn new(pairs: &[(u8, u8)]) -> AClasses {
        let s = syms();
        let ms = mset();
        let best: Vec<(u32, u16)> = (0..N_MSET as u32)
            .into_par_iter()
            .map(|r| {
                let l = ms.unrank(r);
                let mut bst = (u32::MAX, 0u16);
                for (i, &(g, t)) in pairs.iter().enumerate() {
                    let x = ms.rank(&s.apply_a(g as usize, t as usize, &l));
                    if x < bst.0 {
                        bst = (x, i as u16);
                    }
                }
                bst
            })
            .collect();
        let mut id: HashMap<u32, u32> = HashMap::new();
        let mut rep = vec![];
        let mut sym = vec![0u32; N_MSET];
        for r in 0..N_MSET {
            let (x, si) = best[r];
            let next = rep.len() as u32;
            let c = *id.entry(x).or_insert(next);
            if c == next {
                rep.push(x);
            }
            sym[r] = c | (si as u32) << 16;
        }
        AClasses {
            sym: HugeBuf::from_slice(&sym),
            rep,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    /// The symmetries must map the move set to itself (as permutations of the coordinates):
    /// applying (g,t) then a move equals applying the conjugated move then (g,t).
    #[test]
    fn symmetries_commute_with_physical_moves() {
        let s = syms();
        let d = defs();
        let sm = &d.sm;
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(9);
        let ms = mset();
        for _ in 0..20 {
            let c = Cube::random(&mut rng);
            for t in 0..N_GEO {
                let m = &s.mats[t];
                // physical face turn f maps to face t(f), mirrored if improper
                for f in 0..8 {
                    for cw in [true, false] {
                        let perm = sm.face_turn(f, cw);
                        let f2 = face_of_normal(mat_vec(m, FACES[f]));
                        let cw2 = cw == (det(m) == 1);
                        let perm2 = sm.face_turn(f2, cw2);
                        let st = c.to_stickers();
                        // A colouring at sticker level, then repositioned
                        let a1 = permute_stickers(&st, &perm);
                        let xa1: [u8; 12] = std::array::from_fn(|p| {
                            A_FACES
                                .iter()
                                .position(|&x| x as u8 == a1[d.a_pos[p]])
                                .unwrap() as u8
                        });
                        let lhs = s.apply_a(0, t, &xa1);
                        let mut st2 = [0u8; 72];
                        let wp = sm.whole(m);
                        for i in 0..72 {
                            st2[wp[i]] = st[i];
                        }
                        let b1 = permute_stickers(&st2, &perm2);
                        let rhs: [u8; 12] = std::array::from_fn(|p| {
                            A_FACES
                                .iter()
                                .position(|&x| x as u8 == b1[d.a_pos[p]])
                                .unwrap() as u8
                        });
                        assert_eq!(ms.rank(&lhs), ms.rank(&rhs));
                    }
                }
            }
        }
    }

    #[test]
    fn class_count() {
        let a = AClasses::new(&syms().pairs);
        eprintln!("A classes under 288 symmetries: {}", a.rep.len());
        assert!(a.rep.len() < 2000);
    }
}
