//! Coordinates and move tables.
//!
//! * corner coordinate: 960 states (5 movable corners, even permutations, 16 orientations)
//! * A-centre / B-centre colouring: multiset rank of 12 positions with 4 colours x 3 (369600)
//! * edges: "triple" labelling (which A-face each edge belongs to, 369600) plus 4 chirality
//!   bits (parity of the arrangement of each triple's three edges)

use crate::cube::*;
use crate::mem::HugeBuf;
use rayon::prelude::*;
use std::sync::OnceLock;

pub const N_CORNER: usize = 960;
pub const N_MSET: usize = 369_600;
/// edge coset (triple labelling x 3 free chirality bits)
pub const N_ECOSET: usize = N_MSET * 8;

const fn binom(n: usize, k: usize) -> usize {
    if k > n {
        return 0;
    }
    let mut r = 1;
    let mut i = 0;
    while i < k {
        r = r * (n - i) / (i + 1);
        i += 1;
    }
    r
}

/// Ranking of 3-subsets of n positions given as bitmasks.
struct Comb {
    rank: Vec<u16>,   // mask -> rank
    unrank: Vec<u16>, // rank -> mask
}

impl Comb {
    fn new(n: usize) -> Comb {
        let mut rank = vec![u16::MAX; 1 << n];
        let mut unrank = vec![];
        for m in 0..(1usize << n) {
            if m.count_ones() == 3 {
                // colex rank
                let mut r = 0;
                let mut k = 0;
                for i in 0..n {
                    if m >> i & 1 == 1 {
                        k += 1;
                        r += binom(i, k);
                    }
                }
                rank[m] = r as u16;
            }
        }
        unrank.resize(binom(n, 3), 0);
        for m in 0..(1usize << n) {
            if rank[m] != u16::MAX {
                unrank[rank[m] as usize] = m as u16;
            }
        }
        Comb { rank, unrank }
    }
}

/// The bits of `x` at the positions set in `mask`, packed to the bottom (`pext`), or the
/// inverse (`pdep`), bit by bit. Only used to build the 6-bit tables below at compile time.
const fn compress_slow(x: u32, mask: u32) -> u32 {
    let (mut r, mut k, mut m) = (0, 0, mask);
    while m != 0 {
        let b = m.trailing_zeros();
        r |= (x >> b & 1) << k;
        k += 1;
        m &= m - 1;
    }
    r
}

const fn expand_slow(x: u32, mask: u32) -> u32 {
    let (mut r, mut k, mut m) = (0, 0, mask);
    while m != 0 {
        let b = m.trailing_zeros();
        r |= (x >> k & 1) << b;
        k += 1;
        m &= m - 1;
    }
    r
}

/// `[mask][x]` for 6-bit `mask` and `x`, so that a 12-bit compress or expand is two lookups.
/// Tables rather than BMI2 `pext`/`pdep`, which are microcoded and slow on AMD before Zen 3.
/// Speedup: with `MaskPerm`, the start-up tables on one thread: move tables 0.52 -> 0.13 s,
/// symmetry tables 1.06 -> 0.21 s (the tables alone: 0.37 s and 0.58 s).
static PEXT6: [[u8; 64]; 64] = bit_table(false);
static PDEP6: [[u8; 64]; 64] = bit_table(true);

const fn bit_table(dep: bool) -> [[u8; 64]; 64] {
    let mut t = [[0u8; 64]; 64];
    let mut m = 0;
    while m < 64 {
        let mut x = 0;
        while x < 64 {
            let v = if dep {
                expand_slow(x as u32, m as u32)
            } else {
                compress_slow(x as u32, m as u32)
            };
            t[m][x] = v as u8;
            x += 1;
        }
        m += 1;
    }
    t
}

/// `compress_slow(x, mask)` for a `mask` of at most 12 bits.
#[inline]
fn compress(x: u32, mask: u32) -> u32 {
    let lo = (mask & 63) as usize;
    let hi = (mask >> 6 & 63) as usize;
    u32::from(PEXT6[lo][(x & 63) as usize])
        | u32::from(PEXT6[hi][(x >> 6 & 63) as usize]) << lo.count_ones()
}

/// `expand_slow(x, mask)` for a `mask` of at most 12 bits.
#[inline]
fn expand(x: u32, mask: u32) -> u32 {
    let lo = (mask & 63) as usize;
    let hi = (mask >> 6 & 63) as usize;
    let k = lo.count_ones();
    u32::from(PDEP6[lo][(x & ((1 << k) - 1)) as usize])
        | u32::from(PDEP6[hi][(x >> k & 63) as usize]) << 6
}

pub struct Mset {
    c12: Comb,
    c9: Comb,
    c6: Comb,
}

impl Mset {
    fn new() -> Mset {
        Mset {
            c12: Comb::new(12),
            c9: Comb::new(9),
            c6: Comb::new(6),
        }
    }
    /// rank of a labelling of 12 positions with labels 0..4, three each
    pub fn rank(&self, l: &[u8; 12]) -> u32 {
        let mut m = [0u32; 4];
        for (p, &x) in l.iter().enumerate() {
            m[x as usize] |= 1 << p;
        }
        self.rank_masks(m[0], m[1], m[2])
    }
    /// The rank from the position masks of labels 0, 1 and 2 (disjoint, three bits each).
    #[inline]
    pub fn rank_masks(&self, m0: u32, m1: u32, m2: u32) -> u32 {
        let r0 = u32::from(self.c12.rank[m0 as usize]);
        let rest = 0xfff & !m0;
        let r1 = u32::from(self.c9.rank[compress(m1, rest) as usize]);
        let rest2 = rest & !m1;
        let r2 = u32::from(self.c6.rank[compress(m2, rest2) as usize]);
        r0 + 220 * (r1 + 84 * r2)
    }
    /// The position masks of labels 0, 1 and 2 in the labelling of rank `r`.
    #[inline]
    pub fn unrank_masks(&self, r: u32) -> [u32; 3] {
        let (r0, r1, r2) = (r % 220, (r / 220) % 84, r / (220 * 84));
        let m0 = u32::from(self.c12.unrank[r0 as usize]);
        let rest = 0xfff & !m0;
        let m1 = expand(u32::from(self.c9.unrank[r1 as usize]), rest);
        let rest2 = rest & !m1;
        let m2 = expand(u32::from(self.c6.unrank[r2 as usize]), rest2);
        [m0, m1, m2]
    }
    pub fn unrank(&self, r: u32) -> [u8; 12] {
        let [m0, m1, m2] = self.unrank_masks(r);
        // the masks are disjoint: label 3 - 3*[in m0] - 2*[in m1] - [in m2]
        std::array::from_fn(|p| (3 - 3 * (m0 >> p & 1) - 2 * (m1 >> p & 1) - (m2 >> p & 1)) as u8)
    }
}

pub fn mset() -> &'static Mset {
    static M: OnceLock<Mset> = OnceLock::new();
    M.get_or_init(Mset::new)
}

/// Parity of a 3-element sequence (0 = increasing order is an even permutation of it).
#[inline]
pub fn parity3(a: u8, b: u8, c: u8) -> u8 {
    ((a > b) as u8) ^ ((a > c) as u8) ^ ((b > c) as u8)
}

/// Pieces of each A-triple, ascending.
pub fn triple_pieces() -> [[u8; 3]; 4] {
    let d = defs();
    let mut t = [[0u8; 3]; 4];
    let mut n = [0usize; 4];
    for k in 0..12 {
        let a = d.edge_atri[k] as usize;
        t[a][n[a]] = k as u8;
        n[a] += 1;
    }
    t
}

/// Chirality bits of the edge arrangement.
pub fn edge_chirality(ep: &[u8; 12]) -> u8 {
    let mut pos = [0u8; 12];
    for (p, &k) in ep.iter().enumerate() {
        pos[k as usize] = p as u8;
    }
    let tp = triple_pieces();
    let mut bits = 0;
    for t in 0..4 {
        let [a, b, c] = tp[t];
        bits |= parity3(pos[a as usize], pos[b as usize], pos[c as usize]) << t;
    }
    bits
}

/// Per-triple parity change when the positions are mapped by `dest` (labels `l` before).
pub fn triple_flip_mask(l: &[u8; 12], dest: &[u8; 12]) -> u8 {
    // destinations of each triple's positions, in ascending position order (no allocation)
    let mut q = [[0u8; 3]; 4];
    let mut n = [0usize; 4];
    for p in 0..12 {
        let t = l[p] as usize;
        q[t][n[t]] = dest[p];
        n[t] += 1;
    }
    let mut mask = 0;
    for t in 0..4 {
        mask |= parity3(q[t][0], q[t][1], q[t][2]) << t;
    }
    mask
}

/// Parity of the reference permutation of a triple labelling (pieces of each triple placed
/// in ascending order onto its positions in ascending order).
pub fn ref_parity(l: &[u8; 12]) -> u8 {
    ref_parity_with(&triple_pieces(), l)
}

/// `ref_parity` with `triple_pieces()` computed once by the caller.
fn ref_parity_with(tp: &[[u8; 3]; 4], l: &[u8; 12]) -> u8 {
    let mut perm = [0u8; 12];
    let mut n = [0usize; 4];
    for p in 0..12 {
        let t = l[p] as usize;
        perm[p] = tp[t][n[t]];
        n[t] += 1;
    }
    perm_parity(&perm)
}

/// Corner coordinate helpers.
pub struct Corners {
    pub raw_to_id: Vec<u16>,
    pub id_to_state: Vec<([u8; 6], [u8; 6])>,
}

fn corner_raw(cp: [u8; 6], co: [u8; 6]) -> usize {
    // Lehmer code of cp[1..6] (values 1..5) and 5 orientation bits
    let mut r = 0;
    for i in 1..6 {
        let smaller = (i + 1..6).filter(|&j| cp[j] < cp[i]).count();
        r = r * (6 - i) + smaller;
    }
    let mut o = 0;
    for i in 1..6 {
        o |= (co[i] as usize) << (i - 1);
    }
    r * 32 + o
}

impl Corners {
    fn new() -> Corners {
        let mut raw_to_id = vec![u16::MAX; 120 * 32];
        let mut id_to_state = vec![];
        let s0 = Cube::solved();
        let mut q = vec![(s0.cp, s0.co)];
        raw_to_id[corner_raw(s0.cp, s0.co)] = 0;
        let mut i = 0;
        while i < q.len() {
            let (cp, co) = q[i];
            i += 1;
            id_to_state.push((cp, co));
            let c = Cube { cp, co, ..s0 };
            for m in 0..N_MOVES {
                let n = c.apply(m);
                let r = corner_raw(n.cp, n.co);
                if raw_to_id[r] == u16::MAX {
                    raw_to_id[r] = q.len() as u16;
                    q.push((n.cp, n.co));
                }
            }
        }
        debug_assert_eq!(id_to_state.len(), N_CORNER);
        Corners {
            raw_to_id,
            id_to_state,
        }
    }
    #[inline]
    pub fn id(&self, c: &Cube) -> u16 {
        self.raw_to_id[corner_raw(c.cp, c.co)]
    }
}

pub fn corners() -> &'static Corners {
    static C: OnceLock<Corners> = OnceLock::new();
    C.get_or_init(Corners::new)
}

/// All coordinates of a position used by phase 1.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Coords {
    pub c: u16,
    pub a: u32,
    pub b: u32,
    pub e: u32,
    pub ebits: u8,
}

impl Coords {
    pub fn of(cube: &Cube) -> Coords {
        let ms = mset();
        Coords {
            c: corners().id(cube),
            a: ms.rank(&cube.xa),
            b: ms.rank(&cube.xb),
            e: ms.rank(&cube.etri()),
            ebits: edge_chirality(&cube.ep),
        }
    }
}

/// Move tables. `e` entries pack the new triple rank (low 20 bits) and the chirality flip mask.
pub struct MoveTables {
    pub c: Vec<[u16; N_MOVES]>,
    pub a: HugeBuf<[u32; N_MOVES]>,
    pub b: HugeBuf<[u32; N_MOVES]>,
    pub e: HugeBuf<[u32; N_MOVES]>,
    /// reference parity of each triple labelling (bit 3 of the chirality is implied)
    pub eref: Vec<u8>,
}

/// A permutation of the 12 positions, applied to position masks with two lookups.
/// Speedup: the start-up move and symmetry tables permute a labelling's three masks with this
/// instead of unranking to labels, permuting them and ranking again.
pub struct MaskPerm {
    lo: [u16; 64],
    hi: [u16; 64],
    /// `flip[m]` for a mask of three positions p0 < p1 < p2: the parity of the order of their
    /// destinations (`triple_flip_mask` for one triple)
    flip: Vec<u8>,
}

impl MaskPerm {
    /// `dest[p]` is where position `p` goes.
    pub fn new(dest: &[u8; 12]) -> MaskPerm {
        let image = |m: usize, off: usize| -> u16 {
            (0..6)
                .filter(|&i| m >> i & 1 == 1)
                .fold(0, |a, i| a | 1 << dest[i + off])
        };
        let flip = (0..1usize << 12)
            .map(|m| {
                let mut q = (0..12).filter(|&p| m >> p & 1 == 1).map(|p| dest[p]);
                match (q.next(), q.next(), q.next()) {
                    (Some(a), Some(b), Some(c)) => parity3(a, b, c),
                    _ => 0,
                }
            })
            .collect();
        MaskPerm {
            lo: std::array::from_fn(|m| image(m, 0)),
            hi: std::array::from_fn(|m| image(m, 6)),
            flip,
        }
    }
    /// The permutation that makes `n[p] = l[src[p]]` from a labelling `l`.
    pub fn gather(src: &[u8; 12]) -> MaskPerm {
        let mut dest = [0u8; 12];
        for p in 0..12 {
            dest[src[p] as usize] = p as u8;
        }
        MaskPerm::new(&dest)
    }
    #[inline]
    pub fn apply(&self, m: u32) -> u32 {
        u32::from(self.lo[(m & 63) as usize] | self.hi[(m >> 6 & 63) as usize])
    }
    /// The rank of the permuted labelling whose label masks are `m`.
    #[inline]
    pub fn rank(&self, ms: &Mset, m: &[u32; 3]) -> u32 {
        ms.rank_masks(self.apply(m[0]), self.apply(m[1]), self.apply(m[2]))
    }
    /// `triple_flip_mask` of the labelling with label masks `m` under this permutation.
    #[inline]
    pub fn flips(&self, m: &[u32; 3]) -> u8 {
        let m3 = 0xfff ^ m[0] ^ m[1] ^ m[2];
        self.flip[m[0] as usize]
            | self.flip[m[1] as usize] << 1
            | self.flip[m[2] as usize] << 2
            | self.flip[m3 as usize] << 3
    }
}

pub fn label_perm_table(src: &[[u8; 12]; N_MOVES]) -> Vec<[u32; N_MOVES]> {
    let ms = mset();
    let perms: Vec<MaskPerm> = src.iter().map(MaskPerm::gather).collect();
    (0..N_MSET as u32)
        .into_par_iter()
        .map(|r| {
            let m = ms.unrank_masks(r);
            std::array::from_fn(|k| perms[k].rank(ms, &m))
        })
        .collect()
}

impl MoveTables {
    fn new() -> MoveTables {
        let d = defs();
        let cs = corners();
        let c: Vec<[u16; N_MOVES]> = cs
            .id_to_state
            .iter()
            .map(|&(cp, co)| {
                let cube = Cube {
                    cp,
                    co,
                    ..Cube::solved()
                };
                std::array::from_fn(|m| cs.id(&cube.apply(m)))
            })
            .collect();
        let xa_src: [[u8; 12]; N_MOVES] = std::array::from_fn(|m| d.moves[m].xa);
        let xb_src: [[u8; 12]; N_MOVES] = std::array::from_fn(|m| d.moves[m].xb);
        let ep_src: [[u8; 12]; N_MOVES] = std::array::from_fn(|m| d.moves[m].ep);
        let a = label_perm_table(&xa_src);
        let b = label_perm_table(&xb_src);
        let ms = mset();
        let perms: Vec<MaskPerm> = ep_src.iter().map(MaskPerm::gather).collect();
        let e: Vec<[u32; N_MOVES]> = (0..N_MSET as u32)
            .into_par_iter()
            .map(|r| {
                let l = ms.unrank_masks(r);
                std::array::from_fn(|m| perms[m].rank(ms, &l) | u32::from(perms[m].flips(&l)) << 20)
            })
            .collect();
        let tp = triple_pieces();
        let eref: Vec<u8> = (0..N_MSET as u32)
            .into_par_iter()
            .map(|r| ref_parity_with(&tp, &ms.unrank(r)))
            .collect();
        MoveTables {
            c,
            a: HugeBuf::from_slice(&a),
            b: HugeBuf::from_slice(&b),
            e: HugeBuf::from_slice(&e),
            eref,
        }
    }

    #[inline(always)]
    pub fn apply(&self, x: &Coords, m: usize) -> Coords {
        let pe = self.e[x.e as usize][m];
        Coords {
            c: self.c[x.c as usize][m],
            a: self.a[x.a as usize][m],
            b: self.b[x.b as usize][m],
            e: pe & 0xfffff,
            ebits: x.ebits ^ (pe >> 20) as u8,
        }
    }
}

pub fn move_tables() -> &'static MoveTables {
    static T: OnceLock<MoveTables> = OnceLock::new();
    T.get_or_init(MoveTables::new)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    #[test]
    fn bit_tables_match_bit_loops() {
        for mask in 0..1u32 << 12 {
            for x in 0..1u32 << 12 {
                assert_eq!(compress(x, mask), compress_slow(x, mask));
                if x < 1 << mask.count_ones() {
                    assert_eq!(expand(x, mask), expand_slow(x, mask));
                }
            }
        }
    }

    #[test]
    fn move_tables_match_label_construction() {
        let (d, ms, mt) = (defs(), mset(), move_tables());
        let tp = triple_pieces();
        for r in 0..N_MSET as u32 {
            let l = ms.unrank(r);
            assert_eq!(mt.eref[r as usize], ref_parity(&l));
            assert_eq!(ref_parity_with(&tp, &l), ref_parity(&l));
            for m in 0..N_MOVES {
                let gather = |src: &[u8; 12]| -> u32 {
                    ms.rank(&std::array::from_fn(|p| l[src[p] as usize]))
                };
                assert_eq!(mt.a[r as usize][m], gather(&d.moves[m].xa));
                assert_eq!(mt.b[r as usize][m], gather(&d.moves[m].xb));
                let ep = &d.moves[m].ep;
                let mut dest = [0u8; 12];
                for p in 0..12 {
                    dest[ep[p] as usize] = p as u8;
                }
                let want = gather(ep) | u32::from(triple_flip_mask(&l, &dest)) << 20;
                assert_eq!(mt.e[r as usize][m], want);
            }
        }
    }

    #[test]
    fn mset_roundtrip() {
        let ms = mset();
        for r in (0..N_MSET as u32).step_by(997) {
            assert_eq!(ms.rank(&ms.unrank(r)), r);
        }
        assert_eq!(
            ms.rank(&Cube::solved().xa),
            ms.rank(&[0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3])
        );
    }

    #[test]
    fn coords_follow_moves() {
        let mt = move_tables();
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(3);
        for _ in 0..300 {
            let c = Cube::random(&mut rng);
            let x = Coords::of(&c);
            // bit 3 of the chirality is implied by the others and the labelling
            let b = x.ebits;
            assert_eq!((b ^ b >> 1 ^ b >> 2 ^ b >> 3) & 1, mt.eref[x.e as usize]);
            for m in 0..N_MOVES {
                assert_eq!(mt.apply(&x, m), Coords::of(&c.apply(m)));
            }
        }
    }

    #[test]
    fn corner_orientation_sum_is_even() {
        assert_eq!(corners().id_to_state.len(), N_CORNER);
        for &(_, co) in &corners().id_to_state {
            assert_eq!(co.iter().fold(0, |a, &b| a ^ b), 0);
        }
    }
}
