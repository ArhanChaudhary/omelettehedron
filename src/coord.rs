//! Coordinates and move tables.
//!
//! * corner coordinate: 960 states (5 movable corners, even permutations, 16 orientations)
//! * A-centre / B-centre colouring: multiset rank of 12 positions with 4 colours x 3 (369600)
//! * edges: "triple" labelling (which A-face each edge belongs to, 369600) plus 4 chirality
//!   bits (parity of the arrangement of each triple's three edges)

use crate::cube::*;
use rayon::prelude::*;
use std::sync::OnceLock;

pub const N_CORNER: usize = 960;
pub const N_MSET: usize = 369600;
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

#[inline]
fn compress(x: u32, mask: u32) -> u32 {
    let mut r = 0;
    let mut k = 0;
    let mut m = mask;
    while m != 0 {
        let b = m.trailing_zeros();
        if x >> b & 1 == 1 {
            r |= 1 << k;
        }
        k += 1;
        m &= m - 1;
    }
    r
}

#[inline]
fn expand(x: u32, mask: u32) -> u32 {
    let mut r = 0;
    let mut k = 0;
    let mut m = mask;
    while m != 0 {
        let b = m.trailing_zeros();
        if x >> k & 1 == 1 {
            r |= 1 << b;
        }
        k += 1;
        m &= m - 1;
    }
    r
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
        let r0 = self.c12.rank[m[0] as usize] as u32;
        let rest = 0xfff & !m[0];
        let r1 = self.c9.rank[compress(m[1], rest) as usize] as u32;
        let rest2 = rest & !m[1];
        let r2 = self.c6.rank[compress(m[2], rest2) as usize] as u32;
        r0 + 220 * (r1 + 84 * r2)
    }
    pub fn unrank(&self, r: u32) -> [u8; 12] {
        let (r0, r1, r2) = (r % 220, (r / 220) % 84, r / (220 * 84));
        let m0 = self.c12.unrank[r0 as usize] as u32;
        let rest = 0xfff & !m0;
        let m1 = expand(self.c9.unrank[r1 as usize] as u32, rest);
        let rest2 = rest & !m1;
        let m2 = expand(self.c6.unrank[r2 as usize] as u32, rest2);
        let mut l = [3u8; 12];
        for p in 0..12 {
            if m0 >> p & 1 == 1 {
                l[p] = 0;
            } else if m1 >> p & 1 == 1 {
                l[p] = 1;
            } else if m2 >> p & 1 == 1 {
                l[p] = 2;
            }
        }
        l
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
    let mut mask = 0;
    for t in 0..4u8 {
        let q: Vec<u8> = (0..12u8).filter(|&p| l[p as usize] == t).collect();
        mask |= parity3(
            dest[q[0] as usize],
            dest[q[1] as usize],
            dest[q[2] as usize],
        ) << t;
    }
    mask
}

/// Parity of the reference permutation of a triple labelling (pieces of each triple placed
/// in ascending order onto its positions in ascending order).
pub fn ref_parity(l: &[u8; 12]) -> u8 {
    let tp = triple_pieces();
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

fn corner_raw(cp: &[u8; 6], co: &[u8; 6]) -> usize {
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
        raw_to_id[corner_raw(&s0.cp, &s0.co)] = 0;
        let mut i = 0;
        while i < q.len() {
            let (cp, co) = q[i];
            i += 1;
            id_to_state.push((cp, co));
            let c = Cube { cp, co, ..s0 };
            for m in 0..N_MOVES {
                let n = c.apply(m);
                let r = corner_raw(&n.cp, &n.co);
                if raw_to_id[r] == u16::MAX {
                    raw_to_id[r] = q.len() as u16;
                    q.push((n.cp, n.co));
                }
            }
        }
        assert_eq!(id_to_state.len(), N_CORNER);
        Corners {
            raw_to_id,
            id_to_state,
        }
    }
    #[inline]
    pub fn id(&self, c: &Cube) -> u16 {
        self.raw_to_id[corner_raw(&c.cp, &c.co)]
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
    pub a: Vec<[u32; N_MOVES]>,
    pub b: Vec<[u32; N_MOVES]>,
    pub e: Vec<[u32; N_MOVES]>,
    /// reference parity of each triple labelling (bit 3 of the chirality is implied)
    pub eref: Vec<u8>,
}

pub fn label_perm_table(src: &[[u8; 12]; N_MOVES]) -> Vec<[u32; N_MOVES]> {
    let ms = mset();
    (0..N_MSET as u32)
        .into_par_iter()
        .map(|r| {
            let l = ms.unrank(r);
            std::array::from_fn(|m| {
                let n: [u8; 12] = std::array::from_fn(|p| l[src[m][p] as usize]);
                ms.rank(&n)
            })
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
        let dests: [[u8; 12]; N_MOVES] = std::array::from_fn(|m| {
            let mut dst = [0u8; 12];
            for p in 0..12 {
                dst[ep_src[m][p] as usize] = p as u8;
            }
            dst
        });
        let e: Vec<[u32; N_MOVES]> = (0..N_MSET as u32)
            .into_par_iter()
            .map(|r| {
                let l = ms.unrank(r);
                std::array::from_fn(|m| {
                    let n: [u8; 12] = std::array::from_fn(|p| l[ep_src[m][p] as usize]);
                    ms.rank(&n) | (triple_flip_mask(&l, &dests[m]) as u32) << 20
                })
            })
            .collect();
        let eref: Vec<u8> = (0..N_MSET as u32)
            .into_par_iter()
            .map(|r| ref_parity(&ms.unrank(r)))
            .collect();
        MoveTables { c, a, b, e, eref }
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
        for &(_, co) in &corners().id_to_state {
            assert_eq!(co.iter().fold(0, |a, &b| a ^ b), 0);
        }
    }
}
