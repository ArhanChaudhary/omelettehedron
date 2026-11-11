//! The 48 symmetries of the octahedron, acting on positions of the fixed-corner frame (optimal
//! solver).
//!
//! Conjugating a position by a symmetry `M` of the octahedron (stickers moved by `M` and
//! recoloured by `M`) maps the solved state to the solved state and face turns to face turns (a
//! reflection mirrors the turns), and a whole-puzzle rotation only changes the frame. So the
//! distance to solved is the same for a position and for `R M P M^-1`, where the rotation `R`
//! brings the front corner piece home again (the fixed frame of [`crate::cube`]); `R` depends
//! on the corner state. On the coordinates, symmetry `s` with re-fixing rotation `R`:
//! * corners: `c -> conj[s][c]` (the new corner id and `R`);
//! * centres and edges: relabelled by `s` (colours, edge pieces) and repositioned by `t = R s`.
//!
//! The class-swapping symmetries (24..48) exchange the A- and B-centre orbits.
//!
//! These fix the solved state, unlike the two-phase solver's symmetries (`src/sym.rs`), which
//! may relabel and reposition independently because they only need to preserve `H_A`.

use crate::coord::{corners, N_CORNER};
use crate::cube::*;
use crate::geometry::*;
use rayon::prelude::*;
use std::sync::OnceLock;

/// Number of symmetries.
pub const N_SYM: usize = 48;
/// The class-preserving symmetries are `0..N_SYM_CP` (in the order of
/// [`class_preserving_group`], identity first, so the indices agree with `crate::sym::syms()`);
/// `N_SYM_CP + i` is `SWAP_AB` times symmetry `i`.
pub const N_SYM_CP: usize = 24;

pub struct Sym48 {
    pub mats: Vec<Mat>,
    /// `mul[i][j]`: index of `mats[i] * mats[j]` (apply `j` first)
    pub mul: Vec<[u8; N_SYM]>,
    /// index of the inverse
    pub inv: [u8; N_SYM],
    /// index of `SWAP_AB` (the 90° rotation about the fixed corner's axis)
    pub swap: usize,
    /// `conj[s][c]` = (corner id of `R M c M^-1`, index of `R`), for corner id `c`
    pub conj: Vec<Vec<(u16, u8)>>,
    /// centre colours: `ccol[s][o][k]` is the colour index (in the image orbit) of colour `k` of
    /// orbit `o` (0 = A, 1 = B); the image orbit is `o ^ 1` for a class-swapping symmetry
    pub ccol: Vec<[[u8; 4]; 2]>,
    /// centre positions: `cpos[s][o][p]` is the image (in the image orbit) of position `p`
    pub cpos: Vec<[[u8; 12]; 2]>,
    /// edge positions: `epos[s][p]` is the image of edge position `p`; edge pieces are named by
    /// their home position, so it is also the relabelling of the pieces
    pub epos: Vec<[u8; 12]>,
}

/// Whether symmetry `s` exchanges the face classes.
#[inline]
pub const fn swaps_classes(s: usize) -> bool {
    s >= N_SYM_CP
}

/// The symmetries, built once (about 0.1 s, in parallel).
pub fn sym48() -> &'static Sym48 {
    static S: OnceLock<Sym48> = OnceLock::new();
    S.get_or_init(build)
}

/// Unit direction of edge position `e` (the signs of the two coordinates of magnitude 4 that
/// its two stickers share).
fn edge_dir(d: &Defs, e: usize) -> V3 {
    let p1 = d.sm.st[d.edge_slots[e][0]].p;
    let p2 = d.sm.st[d.edge_slots[e][1]].p;
    let f = |x: i32, y: i32| {
        if x.abs() == 4 && y.abs() == 4 {
            x.signum()
        } else {
            0
        }
    };
    [f(p1[0], p2[0]), f(p1[1], p2[1]), f(p1[2], p2[2])]
}

/// The corner id and rotation index of `R M c M^-1` (see the module docs); `rots` are the
/// rotations (index, sticker permutation) to try. `None` only if no rotation re-fixes the
/// front corner, which cannot happen for the fixed geometry (the class-preserving rotations act
/// transitively on the vertices, and the two fixing a vertex differ by a half turn of its
/// corner).
fn conj_corner(
    m: &Mat,
    wp: &[usize],
    rots: &[(usize, &Vec<usize>)],
    c: usize,
) -> Option<(u16, u8)> {
    let (cp, co) = corners().id_to_state[c];
    let st = Cube {
        cp,
        co,
        ..Cube::solved()
    }
    .to_stickers();
    let mut n = [0u8; 72];
    for i in 0..72 {
        n[wp[i]] = face_of_normal(mat_vec(m, FACES[st[i] as usize])) as u8;
    }
    for &(ri, rp) in rots {
        let mut n2 = [0u8; 72];
        for i in 0..72 {
            n2[rp[i]] = n[i];
        }
        if let Some(cube) = Cube::from_stickers(&n2) {
            if cube.cp[0] == 0 && cube.co[0] == 0 {
                return Some((corners().id(&cube), ri as u8));
            }
        }
    }
    None
}

fn build() -> Sym48 {
    let d = defs();
    let sm = &d.sm;
    let cp = class_preserving_group();
    debug_assert_eq!(cp.len(), N_SYM_CP);
    let mut mats = cp.clone();
    mats.extend(cp.iter().map(|g| mat_mul(&SWAP_AB, g)));
    let find = |m: &Mat| mats.iter().position(|x| x == m);
    // (the 48 matrices form a group: products and inverses are found; 0 is never used)
    let mul: Vec<[u8; N_SYM]> = (0..N_SYM)
        .map(|i| std::array::from_fn(|j| find(&mat_mul(&mats[i], &mats[j])).unwrap_or(0) as u8))
        .collect();
    let inv: [u8; N_SYM] =
        std::array::from_fn(|i| (0..N_SYM).find(|&j| mul[i][j] == 0).unwrap_or(0) as u8);
    let swap = find(&SWAP_AB).unwrap_or(N_SYM_CP);
    debug_assert_eq!(swap, N_SYM_CP);
    let wps: Vec<Vec<usize>> = mats.iter().map(|m| sm.whole(m)).collect();
    // the class-preserving rotations re-fix the front corner
    let rots: Vec<(usize, &Vec<usize>)> = (0..N_SYM_CP)
        .filter(|&r| det(&mats[r]) == 1)
        .map(|r| (r, &wps[r]))
        .collect();
    let conj: Vec<Vec<(u16, u8)>> = (0..N_SYM)
        .into_par_iter()
        .map(|s| {
            (0..N_CORNER)
                .map(|c| {
                    let r = conj_corner(&mats[s], &wps[s], &rots, c);
                    debug_assert!(r.is_some(), "no rotation re-fixes corner {c} after {s}");
                    r.unwrap_or((c as u16, 0))
                })
                .collect()
        })
        .collect();
    let mut ccol = vec![[[0u8; 4]; 2]; N_SYM];
    let mut cpos = vec![[[0u8; 12]; 2]; N_SYM];
    let mut epos = vec![[0u8; 12]; N_SYM];
    let dirs: Vec<V3> = (0..12).map(|e| edge_dir(d, e)).collect();
    for (s, m) in mats.iter().enumerate() {
        for o in 0..2 {
            let (faces, pos) = if o == 0 {
                (&A_FACES, &d.a_pos)
            } else {
                (&B_FACES, &d.b_pos)
            };
            for k in 0..4 {
                ccol[s][o][k] = class_index(face_of_normal(mat_vec(m, FACES[faces[k]]))) as u8;
            }
            for p in 0..12 {
                let q = sm.index_of(mat_vec(m, sm.st[pos[p]].p));
                debug_assert_eq!(is_a_face(q / 9), (o == 0) != swaps_classes(s));
                cpos[s][o][p] = centre_index(q) as u8;
            }
        }
        for p in 0..12 {
            let dir = mat_vec(m, dirs[p]);
            let e = dirs.iter().position(|&x| x == dir);
            debug_assert!(e.is_some(), "edge direction {dir:?}");
            epos[s][p] = e.unwrap_or(0) as u8;
        }
    }
    Sym48 {
        mats,
        mul,
        inv,
        swap,
        conj,
        ccol,
        cpos,
        epos,
    }
}

impl Sym48 {
    /// The reposition symmetry `R s` of corner id `c` under `s` (see the module docs).
    #[inline]
    pub fn repos(&self, s: usize, c: usize) -> usize {
        let r = self.conj[s][c].1 as usize;
        self.mul[r][s] as usize
    }

    /// The conjugate `R M P M^-1` of a whole position through the sticker model, and `R`'s index
    /// (the reference for the coordinate maps; `None` for an invalid position).
    pub fn conj_cube(&self, cube: &Cube, s: usize) -> Option<(Cube, usize)> {
        let sm = &defs().sm;
        let m = &self.mats[s];
        let st = cube.to_stickers();
        let wp = sm.whole(m);
        let mut n = [0u8; 72];
        for i in 0..72 {
            n[wp[i]] = face_of_normal(mat_vec(m, FACES[st[i] as usize])) as u8;
        }
        for r in (0..N_SYM_CP).filter(|&r| det(&self.mats[r]) == 1) {
            let rp = sm.whole(&self.mats[r]);
            let mut n2 = [0u8; 72];
            for i in 0..72 {
                n2[rp[i]] = n[i];
            }
            if let Some(c) = Cube::from_stickers(&n2) {
                if c.cp[0] == 0 && c.co[0] == 0 {
                    return Some((c, r));
                }
            }
        }
        None
    }

    /// The conjugate of a whole position by `s`, from the coordinate maps (corners from
    /// `conj`, centres and edges relabelled by `s` and repositioned by `R s`).
    pub fn conj_by_maps(&self, cube: &Cube, s: usize) -> Cube {
        let cs = corners();
        let c = cs.id(cube) as usize % N_CORNER;
        let (c2, _) = self.conj[s][c];
        let t = self.repos(s, c);
        let (cp, co) = cs.id_to_state[c2 as usize % N_CORNER];
        let mut out = Cube { cp, co, ..*cube };
        let sw = usize::from(swaps_classes(s));
        // orbit o's centres land in orbit o ^ sw
        for (o, src) in [(0, &cube.xa), (1, &cube.xb)] {
            let dst = if o ^ sw == 0 {
                &mut out.xa
            } else {
                &mut out.xb
            };
            for p in 0..12 {
                dst[self.cpos[t][o][p] as usize] = self.ccol[s][o][src[p] as usize];
            }
        }
        for p in 0..12 {
            out.ep[self.epos[t][p] as usize] = self.epos[s][cube.ep[p] as usize];
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    #[test]
    fn group_structure() {
        let s = sym48();
        assert_eq!(s.mats.len(), N_SYM);
        let mut sorted = s.mats.clone();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(sorted.len(), N_SYM);
        for i in 0..N_SYM {
            assert_eq!(s.mul[i][s.inv[i] as usize], 0);
            assert_eq!(swaps_classes(i), {
                let f = face_of_normal(mat_vec(&s.mats[i], FACES[A_FACES[0]]));
                !is_a_face(f)
            });
        }
        // the first 24 are the two-phase solver's symmetries, in the same order
        assert_eq!(&s.mats[..N_SYM_CP], &crate::sym::syms().mats[..]);
    }

    /// The coordinate maps give the same conjugate as the sticker model, for every symmetry,
    /// and the conjugate has the same moves (conjugated) as the position.
    #[test]
    fn maps_match_sticker_conjugation() {
        let s = sym48();
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(48);
        for _ in 0..60 {
            let c = Cube::random(&mut rng);
            for k in 0..N_SYM {
                let (want, r) = s.conj_cube(&c, k).unwrap();
                assert_eq!(s.conj[k][corners().id(&c) as usize].1 as usize, r);
                assert_eq!(s.conj_by_maps(&c, k), want, "symmetry {k}");
            }
        }
        // the solved state is fixed by every symmetry
        for k in 0..N_SYM {
            assert_eq!(s.conj_by_maps(&Cube::solved(), k), Cube::solved());
        }
    }

    /// `SWAP_AB` fixes the front corner, so its conjugation is `Cube::conjugate`.
    #[test]
    fn swap_is_cube_conjugate() {
        let s = sym48();
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(49);
        for _ in 0..100 {
            let c = Cube::random(&mut rng);
            assert_eq!(
                Some(s.conj_by_maps(&c, s.swap)),
                c.conjugate(&SWAP_AB),
                "{c:?}"
            );
            assert_eq!(s.repos(s.swap, corners().id(&c) as usize), s.swap);
        }
    }
}
