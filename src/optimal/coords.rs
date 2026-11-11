//! Positions of the optimal solver and the coordinates its tables index.
//!
//! A position ([`Pos`]) is kept as byte vectors (edge permutation, A- and B-centre colourings),
//! so that a move is one byte shuffle each, plus the corner id. The tables index:
//! * the edge permutation: 12 edges, even permutations, 239,500,800 ([`eperm_rank`]);
//! * corners x B-centre colouring, reduced by the 24 class-preserving symmetries
//!   ([`CcIndex`]: 51 corner classes x 369,600 colourings). The corners x A-centres are read
//!   from the same table: conjugating by `SWAP_AB` (which fixes the solved state and the front
//!   corner) turns the A-centres into B-centres, and the index composes that conjugation with
//!   the reduction, so no second position is carried.

use super::sym::*;
use crate::coord::{corners, mset, N_CORNER, N_MSET};
use crate::cube::*;
use crate::fast::{fast, masks, mset_labels, shuffle, FastTables};
use std::collections::HashMap;
use std::sync::OnceLock;

/// Even permutations of the 12 edges.
pub const N_EPERM: usize = 239_500_800;

/// A position: lanes 0..12 of each vector hold the data, lanes 12..16 are zero.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[repr(C, align(16))]
pub struct Pos {
    /// edge piece at each edge position
    pub ep: [u8; 16],
    /// A-colour index at each A-centre position
    pub xa: [u8; 16],
    /// B-colour index at each B-centre position
    pub xb: [u8; 16],
    /// corner id ([`crate::coord::corners`]; 0 is solved)
    pub c: u16,
}

fn pad(v: &[u8; 12]) -> [u8; 16] {
    let mut o = [0u8; 16];
    o[..12].copy_from_slice(v);
    o
}

impl Pos {
    pub fn of(cube: &Cube) -> Pos {
        Pos {
            ep: pad(&cube.ep),
            xa: pad(&cube.xa),
            xb: pad(&cube.xb),
            c: corners().id(cube),
        }
    }

    pub fn solved() -> Pos {
        Pos::of(&Cube::solved())
    }

    pub fn to_cube(&self) -> Cube {
        let (cp, co) = corners().id_to_state[self.c as usize % N_CORNER];
        let mut c = Cube {
            cp,
            co,
            ..Cube::solved()
        };
        c.ep.copy_from_slice(&self.ep[..12]);
        c.xa.copy_from_slice(&self.xa[..12]);
        c.xb.copy_from_slice(&self.xb[..12]);
        c
    }
}

/// Move tables of the 16 face-metric moves.
pub struct OptMoves {
    /// gather shuffles: `new[p] = old[src[p]]` (lanes 12..16 select zero)
    shuf_e: [[u8; 16]; N_MOVES],
    shuf_a: [[u8; 16]; N_MOVES],
    shuf_b: [[u8; 16]; N_MOVES],
    /// corner id after each move, padded to 1024 ids so that `& 1023` replaces bounds checks
    corner: Box<[[u16; N_MOVES]; 1024]>,
}

fn gather(src: &[u8; 12]) -> [u8; 16] {
    let mut s = [0x80u8; 16];
    s[..12].copy_from_slice(src);
    s
}

impl OptMoves {
    fn new() -> OptMoves {
        let d = defs();
        let cs = corners();
        #[allow(clippy::large_stack_arrays)] // 32 KB, once at startup
        let mut corner = Box::new([[0u16; N_MOVES]; 1024]);
        for (id, &(cp, co)) in cs.id_to_state.iter().enumerate() {
            let cube = Cube {
                cp,
                co,
                ..Cube::solved()
            };
            corner[id] = std::array::from_fn(|m| cs.id(&cube.apply(m)));
        }
        OptMoves {
            shuf_e: std::array::from_fn(|m| gather(&d.moves[m].ep)),
            shuf_a: std::array::from_fn(|m| gather(&d.moves[m].xa)),
            shuf_b: std::array::from_fn(|m| gather(&d.moves[m].xb)),
            corner,
        }
    }

    /// Apply move `m` (`m < N_MOVES`).
    #[inline(always)]
    pub fn apply(&self, p: &Pos, m: usize) -> Pos {
        let m = m & 15;
        Pos {
            ep: shuffle(&p.ep, &self.shuf_e[m]),
            xa: shuffle(&p.xa, &self.shuf_a[m]),
            xb: shuffle(&p.xb, &self.shuf_b[m]),
            c: self.corner[p.c as usize & 1023][m],
        }
    }

    /// Apply move `m` to the edge permutation only.
    #[inline(always)]
    pub fn apply_ep(&self, ep: &[u8; 16], m: usize) -> [u8; 16] {
        shuffle(ep, &self.shuf_e[m & 15])
    }

    /// Corner id after move `m`.
    #[inline(always)]
    pub fn corner(&self, c: u16, m: usize) -> u16 {
        self.corner[c as usize & 1023][m & 15]
    }

    /// Apply move `m` to the corners and B-centres only.
    #[inline(always)]
    pub fn apply_cb(&self, c: u16, xb: &[u8; 16], m: usize) -> (u16, [u8; 16]) {
        let m = m & 15;
        (
            self.corner[c as usize & 1023][m],
            shuffle(xb, &self.shuf_b[m]),
        )
    }
}

/// The move tables, built once.
pub fn opt_moves() -> &'static OptMoves {
    static M: OnceLock<OptMoves> = OnceLock::new();
    M.get_or_init(OptMoves::new)
}

/// Rank of an even permutation of 12 (lanes 0..12; values 0..12): the mixed-radix Lehmer code of
/// its first 10 entries (the last two follow from the parity), `0..N_EPERM`. The identity is 0.
#[inline(always)]
pub fn eperm_rank(ep: &[u8; 16]) -> u32 {
    let mut seen = 0u32;
    let mut r = 0u32;
    for i in 0..10 {
        let p = u32::from(ep[i] & 15);
        // (entries below 12: the shift stays in range)
        let c = p - (seen & ((1 << p) - 1)).count_ones();
        r = r * (12 - i as u32) + c;
        seen |= 1 << p;
    }
    r
}

/// The even permutation of rank `r` ([`eperm_rank`]).
pub fn eperm_unrank(r: u32) -> [u8; 16] {
    let mut digits = [0u32; 10];
    let mut x = r;
    for i in (0..10).rev() {
        let k = 12 - i as u32;
        digits[i] = x % k;
        x /= k;
    }
    let mut avail = 0xfffu32;
    let mut ep = [0u8; 16];
    let mut inversions = 0;
    for i in 0..10 {
        // the digits[i]-th remaining value
        let mut a = avail;
        for _ in 0..digits[i] {
            a &= a - 1;
        }
        let p = a.trailing_zeros();
        ep[i] = p as u8;
        avail &= !(1 << p);
        inversions += digits[i];
    }
    let lo = avail.trailing_zeros() as u8;
    let hi = avail.ilog2() as u8;
    // the last two in the order that makes the permutation even
    if inversions % 2 == 0 {
        (ep[10], ep[11]) = (lo, hi);
    } else {
        (ep[10], ep[11]) = (hi, lo);
    }
    ep
}

/// Index of the corner x centre-colouring table (distance to solved of the corners and the
/// B-centres), reduced by the 24 class-preserving symmetries; the A-centres are looked up
/// through `SWAP_AB` (see the module docs).
///
/// Entry `class * N_MSET + row`: `class` is the symmetry class of the corner state (51), and
/// `row` the rank of the B-colouring after the symmetry that maps the corner state to its class
/// representative. For a lookup through the A-centres, the symmetry is that one composed with
/// `SWAP_AB` (class swapping), which maps the A-centres onto the B-centre positions and
/// colours.
pub struct CcIndex {
    /// per corner id (padded to 1024) and orbit (0: the A-centres through `SWAP_AB`, 1: the
    /// B-centres): (class, relabelling symmetry s, repositioning symmetry t)
    kt: Box<[[(u16, u8, u8); 2]; 1024]>,
    /// class -> representative corner id
    pub rep: Vec<u16>,
    /// per class: the (s, t) pairs other than the identity that fix the representative
    stab: Vec<Vec<(u8, u8)>>,
    /// per symmetry s (padded to 64): the colour map from the orbit it maps onto the B-centres
    /// (B for a class-preserving s, A for a class-swapping one)
    relab: Box<[[u8; 16]; 64]>,
    /// per symmetry t (padded to 64): gather of the B positions from that orbit's positions
    gath: Box<[[u8; 16]; 64]>,
    /// `fast()` (for the multiset rank)
    f: &'static FastTables,
}

/// Entries of the corner x centre table: corner classes x colourings.
pub fn cc_size(ix: &CcIndex) -> usize {
    ix.rep.len() * N_MSET
}

impl CcIndex {
    pub fn new() -> CcIndex {
        let s = sym48();
        // class representative: the smallest corner id among the images under the 24
        // class-preserving symmetries, and the first symmetry that gives it
        let best = |c: usize| -> (u16, usize) {
            (0..N_SYM_CP)
                .map(|k| (s.conj[k][c].0, k))
                .min()
                .unwrap_or((c as u16, 0))
        };
        let mut rep: Vec<u16> = (0..N_CORNER).map(|c| best(c).0).collect();
        rep.sort_unstable();
        rep.dedup();
        let cls: HashMap<u16, u16> = rep
            .iter()
            .enumerate()
            .map(|(i, &r)| (r, i as u16))
            .collect();
        #[allow(clippy::large_stack_arrays)] // 10 KB, once
        let mut kt = Box::new([[(0u16, 0u8, 0u8); 2]; 1024]);
        for c in 0..N_CORNER {
            let (r, k) = best(c);
            kt[c][1] = (
                cls.get(&r).copied().unwrap_or(0),
                k as u8,
                s.repos(k, c) as u8,
            );
        }
        // the A-centres: conjugate by SWAP_AB (into B-centres), then reduce as above
        let w = s.swap;
        for c in 0..N_CORNER {
            let (cw, _) = s.conj[w][c];
            let tw = s.repos(w, c);
            let (cl, k, t) = kt[cw as usize % N_CORNER][1];
            kt[c][0] = (cl, s.mul[k as usize][w], s.mul[t as usize][tw]);
        }
        let stab = rep
            .iter()
            .map(|&r| {
                (1..N_SYM_CP)
                    .filter(|&k| s.conj[k][r as usize].0 == r)
                    .map(|k| (k as u8, s.repos(k, r as usize) as u8))
                    .collect()
            })
            .collect();
        // a symmetry maps orbit src onto the B-centres
        let src = |k: usize| usize::from(!swaps_classes(k));
        let mut relab = Box::new([[0u8; 16]; 64]);
        let mut gath = Box::new([[0x80u8; 16]; 64]);
        for k in 0..N_SYM {
            relab[k][..4].copy_from_slice(&s.ccol[k][src(k)]);
            for p in 0..12 {
                gath[k][s.cpos[k][src(k)][p] as usize] = p as u8;
            }
        }
        CcIndex {
            kt,
            rep,
            stab,
            relab,
            gath,
            f: fast(),
        }
    }

    /// Rank of the colouring `x` after symmetry pair (s, t).
    #[inline(always)]
    fn row(&self, x: &[u8; 16], s: usize, t: usize) -> usize {
        let y = shuffle(&shuffle(&self.relab[s & 63], x), &self.gath[t & 63]);
        self.f.rank(&masks(&y)) as usize
    }

    /// Index of the corners and B-centres `xb`.
    #[inline(always)]
    pub fn index_b(&self, c: u16, xb: &[u8; 16]) -> usize {
        let (cl, s, t) = self.kt[c as usize & 1023][1];
        cl as usize * N_MSET + self.row(xb, s as usize, t as usize)
    }

    /// Index of the corners and A-centres `xa` (the B-centre entry of the position conjugated
    /// by `SWAP_AB`).
    #[inline(always)]
    pub fn index_a(&self, c: u16, xa: &[u8; 16]) -> usize {
        let (cl, s, t) = self.kt[c as usize & 1023][0];
        cl as usize * N_MSET + self.row(xa, s as usize, t as usize)
    }

    /// A position (corner id, B-colouring) of index `i`.
    pub fn decode(&self, i: usize) -> (u16, [u8; 16]) {
        let c = self.rep.get(i / N_MSET).copied().unwrap_or(0);
        (c, mset_labels()[i % N_MSET])
    }

    /// The other indices of the position of index `i` (its corner representative has a
    /// non-trivial stabiliser).
    pub fn equivalents(&self, i: usize, out: &mut Vec<usize>) {
        let cl = i / N_MSET;
        let Some(stab) = self.stab.get(cl) else {
            return;
        };
        if stab.is_empty() {
            return;
        }
        let (_, x) = self.decode(i);
        for &(s, t) in stab {
            out.push(cl * N_MSET + self.row(&x, s as usize, t as usize));
        }
    }

    /// The index of the solved state and its equivalents.
    pub fn targets(&self) -> Vec<usize> {
        let p = Pos::solved();
        let i = self.index_b(p.c, &p.xb);
        let mut v = vec![i];
        self.equivalents(i, &mut v);
        v.sort_unstable();
        v.dedup();
        v
    }
}

/// Index of the corner x edge-permutation table (distance to solved of the corners and edges),
/// reduced by all 48 symmetries of the octahedron (the class-swapping ones too: an edge has one
/// face of each class, so they map edges to edges).
///
/// Entry `class * N_EPERM + row`: `class` is the symmetry class of the corner state, and `row`
/// the rank of the edge permutation after the symmetry that maps the corner state to its class
/// representative (edge pieces relabelled by `s`, positions moved by `t = R s`).
pub struct CeIndex {
    /// per corner id (padded to 1024): (class, relabelling symmetry s, repositioning
    /// symmetry t)
    kt: Box<[(u16, u8, u8); 1024]>,
    /// class -> representative corner id
    pub rep: Vec<u16>,
    /// per class: the (s, t) pairs other than the identity that fix the representative
    stab: Vec<Vec<(u8, u8)>>,
    /// per symmetry (padded to 64): edge piece relabelling and edge position gather
    relab: Box<[[u8; 16]; 64]>,
    gath: Box<[[u8; 16]; 64]>,
}

impl CeIndex {
    pub fn new() -> CeIndex {
        let s = sym48();
        let best = |c: usize| -> (u16, usize) {
            (0..N_SYM)
                .map(|k| (s.conj[k][c].0, k))
                .min()
                .unwrap_or((c as u16, 0))
        };
        let mut rep: Vec<u16> = (0..N_CORNER).map(|c| best(c).0).collect();
        rep.sort_unstable();
        rep.dedup();
        let cls: HashMap<u16, u16> = rep
            .iter()
            .enumerate()
            .map(|(i, &r)| (r, i as u16))
            .collect();
        let mut kt = Box::new([(0u16, 0u8, 0u8); 1024]);
        for c in 0..N_CORNER {
            let (r, k) = best(c);
            kt[c] = (
                cls.get(&r).copied().unwrap_or(0),
                k as u8,
                s.repos(k, c) as u8,
            );
        }
        let stab = rep
            .iter()
            .map(|&r| {
                (1..N_SYM)
                    .filter(|&k| s.conj[k][r as usize].0 == r)
                    .map(|k| (k as u8, s.repos(k, r as usize) as u8))
                    .collect()
            })
            .collect();
        let mut relab = Box::new([[0u8; 16]; 64]);
        let mut gath = Box::new([[0x80u8; 16]; 64]);
        for k in 0..N_SYM {
            relab[k][..12].copy_from_slice(&s.epos[k]);
            for p in 0..12 {
                gath[k][s.epos[k][p] as usize] = p as u8;
            }
        }
        CeIndex {
            kt,
            rep,
            stab,
            relab,
            gath,
        }
    }

    /// Entries: corner classes x even edge permutations.
    /// (More than 2^32: 64-bit targets only; saturates elsewhere.)
    pub fn size(&self) -> usize {
        self.rep.len().saturating_mul(N_EPERM)
    }

    /// Rank of the edge permutation `ep` after symmetry pair (s, t).
    #[inline(always)]
    fn row(&self, ep: &[u8; 16], s: usize, t: usize) -> usize {
        eperm_rank(&shuffle(
            &shuffle(&self.relab[s & 63], ep),
            &self.gath[t & 63],
        )) as usize
    }

    /// Index of the corners `c` and edges `ep`.
    #[inline(always)]
    pub fn index(&self, c: u16, ep: &[u8; 16]) -> usize {
        let (cl, s, t) = self.kt[c as usize & 1023];
        cl as usize * N_EPERM + self.row(ep, s as usize, t as usize)
    }

    /// A position (corner id, edge permutation) of index `i`.
    pub fn decode(&self, i: usize) -> (u16, [u8; 16]) {
        let c = self.rep.get(i / N_EPERM).copied().unwrap_or(0);
        (c, eperm_unrank((i % N_EPERM) as u32))
    }

    /// The other indices of the position of index `i`.
    pub fn equivalents(&self, i: usize, out: &mut Vec<usize>) {
        let cl = i / N_EPERM;
        let Some(stab) = self.stab.get(cl) else {
            return;
        };
        if stab.is_empty() {
            return;
        }
        let (_, ep) = self.decode(i);
        for &(s, t) in stab {
            out.push(cl * N_EPERM + self.row(&ep, s as usize, t as usize));
        }
    }

    /// The index of the solved state (its equivalents are itself).
    pub fn target(&self) -> usize {
        let p = Pos::solved();
        self.index(p.c, &p.ep)
    }
}

/// Multiset rank of a labelling (lanes 0..12), for tests and code off the hot path.
pub fn mset_rank(x: &[u8; 16]) -> u32 {
    let mut l = [0u8; 12];
    l.copy_from_slice(&x[..12]);
    mset().rank(&l)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{Rng, SeedableRng};

    #[test]
    fn eperm_rank_roundtrip() {
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(5);
        assert_eq!(eperm_rank(&Pos::solved().ep), 0);
        for _ in 0..20_000 {
            let r = rng.gen_range(0..N_EPERM as u32);
            let ep = eperm_unrank(r);
            assert_eq!(perm_parity(&ep[..12]), 0);
            let mut sorted = ep;
            sorted[..12].sort_unstable();
            assert_eq!(
                &sorted[..12],
                &std::array::from_fn::<u8, 12, _>(|i| i as u8)
            );
            assert_eq!(eperm_rank(&ep), r);
        }
        // the extremes
        for r in [0, 1, N_EPERM as u32 - 1] {
            assert_eq!(eperm_rank(&eperm_unrank(r)), r);
        }
        // every even permutation reached by moves has a rank in range
        let mut c = Cube::solved();
        for _ in 0..1000 {
            c = c.apply(rng.gen_range(0..N_MOVES));
            let ep = Pos::of(&c).ep;
            assert!((eperm_rank(&ep) as usize) < N_EPERM);
            assert_eq!(eperm_unrank(eperm_rank(&ep)), ep);
        }
    }

    #[test]
    fn moves_match_cube() {
        let mt = opt_moves();
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(6);
        for _ in 0..300 {
            let c = Cube::random(&mut rng);
            let p = Pos::of(&c);
            assert_eq!(p.to_cube(), c);
            for m in 0..N_MOVES {
                assert_eq!(mt.apply(&p, m), Pos::of(&c.apply(m)));
            }
        }
    }

    /// The corner x edge index is invariant (up to equivalent indices) under all 48 symmetries,
    /// and depends only on the corners and edges.
    #[test]
    fn ce_index_is_symmetric() {
        let ix = CeIndex::new();
        let s = sym48();
        eprintln!("corner classes under 48 symmetries: {}", ix.rep.len());
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(8);
        for _ in 0..300 {
            let c = Cube::random(&mut rng);
            let p = Pos::of(&c);
            let i = ix.index(p.c, &p.ep);
            let mut eq = vec![i];
            ix.equivalents(i, &mut eq);
            for k in 0..N_SYM {
                let q = Pos::of(&s.conj_by_maps(&c, k));
                assert!(eq.contains(&ix.index(q.c, &q.ep)), "symmetry {k}");
            }
            let (dc, de) = ix.decode(i);
            assert_eq!(ix.index(dc, &de), i);
        }
        let mut eq = vec![];
        ix.equivalents(ix.target(), &mut eq);
        assert!(eq.iter().all(|&j| j == ix.target()));
    }

    /// The A lookup equals the B lookup of the `SWAP_AB` conjugate; the B index depends only on
    /// the corners and B-centres and is invariant under the 24 class-preserving symmetries.
    #[test]
    fn cc_index_matches_conjugation() {
        let ix = CcIndex::new();
        assert_eq!(ix.rep.len(), 51);
        let s = sym48();
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(7);
        for _ in 0..500 {
            let c = Cube::random(&mut rng);
            let p = Pos::of(&c);
            let w = Pos::of(&c.conjugate(&SWAP_AB).unwrap());
            assert_eq!(ix.index_a(p.c, &p.xa), ix.index_b(w.c, &w.xb));
            let i = ix.index_b(p.c, &p.xb);
            // a conjugate has the same index, or an equivalent one (when the corner class
            // representative has a non-trivial stabiliser)
            let mut eq = vec![i];
            ix.equivalents(i, &mut eq);
            for k in 0..N_SYM_CP {
                let q = Pos::of(&s.conj_by_maps(&c, k));
                assert!(eq.contains(&ix.index_b(q.c, &q.xb)), "symmetry {k}");
            }
            // a decoded entry has its own index
            let (dc, dx) = ix.decode(i);
            assert_eq!(ix.index_b(dc, &dx), i);
        }
        assert_eq!(ix.targets().len(), 1);
    }
}
