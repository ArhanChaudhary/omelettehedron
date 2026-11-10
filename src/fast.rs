//! Cache-resident coordinates for the phase-1 search.
//!
//! Centre colourings and the edge triple labelling are kept as 12-byte label vectors, so a
//! move or a symmetry is a single byte shuffle (`pshufb`). Ranks are computed from per-label
//! bit masks with small tables, and edge chirality is updated with bitset lookups. All of these
//! tables fit in L1/L2, so the only random DRAM accesses left per node are the pruning tables
//! and the 1.5 MB symmetry-class table.

use crate::coord::*;
use crate::cube::*;
use crate::sym::*;
use std::sync::OnceLock;

/// Phase-1 position: labels in lanes 0..12, lanes 12..16 are zero padding.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(C, align(16))]
pub struct Fc {
    pub a: [u8; 16],
    pub b: [u8; 16],
    pub e: [u8; 16],
    pub c: u16,
    pub ebits: u8,
}

pub struct FastTables {
    /// gather shuffles of the 16 moves: new[p] = old[src[p]]
    shuf_a: [[u8; 16]; N_MOVES],
    shuf_b: [[u8; 16]; N_MOVES],
    shuf_e: [[u8; 16]; N_MOVES],
    /// chirality flip of a triple occupying a 12-bit position mask, per move (bitsets)
    flip_m: Box<[[u64; 64]; N_MOVES]>,
    corner: Vec<[u16; N_MOVES]>,
    pub(crate) c12: Box<[u16; 4096]>,
    /// rank tables pre-scaled by their weight in the multiset rank
    r12: Box<[u32; 4096]>,
    r9: Box<[u32; 512]>,
    r6: [u32; 64],
    cmp6: Box<[[u8; 64]; 64]>,
    /// per geometric symmetry: label map (relabel), position gather (reposition), flips
    /// (padded to 32 symmetries so that `& 31` replaces bounds checks)
    relab_map: [[u8; 16]; 32],
    repos_gather: [[u8; 16]; 32],
    relab_bits: [[u8; 16]; 32],
    pub(crate) flip_t: Box<[[u64; 64]; 32]>,
}

fn boxed<T: Copy, const N: usize>(v: &[T]) -> Box<[T; N]> {
    v.to_vec()
        .into_boxed_slice()
        .try_into()
        .ok()
        .expect("table size")
}

fn comb_rank(n: usize) -> Vec<u16> {
    // colex rank of 3-subsets of n positions, indexed by bitmask
    let mut r = vec![0u16; 1 << n];
    for m in 0..(1usize << n) {
        if m.count_ones() == 3 {
            let mut x = 0usize;
            let mut k = 0;
            for i in 0..n {
                if m >> i & 1 == 1 {
                    k += 1;
                    let mut b = 1usize;
                    for j in 0..k {
                        b = b * (i - j) / (j + 1);
                    }
                    x += b;
                }
            }
            r[m] = x as u16;
        }
    }
    r
}

fn bitset_flip(dest: &[u8; 12]) -> [u64; 64] {
    let mut bs = [0u64; 64];
    for m in 0..4096usize {
        if m.count_ones() == 3 {
            let q: Vec<usize> = (0..12).filter(|&p| m >> p & 1 == 1).collect();
            if parity3(dest[q[0]], dest[q[1]], dest[q[2]]) == 1 {
                bs[m >> 6] |= 1 << (m & 63);
            }
        }
    }
    bs
}

fn pad(src: &[u8; 12]) -> [u8; 16] {
    let mut s = [0x80u8; 16];
    s[..12].copy_from_slice(src);
    s
}

impl FastTables {
    fn new() -> FastTables {
        let d = defs();
        let mt_c: Vec<[u16; N_MOVES]> = {
            let cs = corners();
            cs.id_to_state
                .iter()
                .map(|&(cp, co)| {
                    let cube = Cube {
                        cp,
                        co,
                        ..Cube::solved()
                    };
                    std::array::from_fn(|m| cs.id(&cube.apply(m)))
                })
                .collect()
        };
        let shuf_a = std::array::from_fn(|m| pad(&d.moves[m].xa));
        let shuf_b = std::array::from_fn(|m| pad(&d.moves[m].xb));
        let shuf_e = std::array::from_fn(|m| pad(&d.moves[m].ep));
        let flip_m: Vec<[u64; 64]> = (0..N_MOVES)
            .map(|m| {
                let mut dest = [0u8; 12];
                for p in 0..12 {
                    dest[d.moves[m].ep[p] as usize] = p as u8;
                }
                bitset_flip(&dest)
            })
            .collect();
        let mut cmp6 = vec![[0u8; 64]; 64];
        for mask in 0..64usize {
            for x in 0..64usize {
                let mut r = 0;
                let mut k = 0;
                for b in 0..6 {
                    if mask >> b & 1 == 1 {
                        if x >> b & 1 == 1 {
                            r |= 1 << k;
                        }
                        k += 1;
                    }
                }
                cmp6[mask][x] = r;
            }
        }
        let s = syms();
        let relab_map = std::array::from_fn(|g| {
            let mut m = [0u8; 16];
            if g < N_GEO {
                m[..4].copy_from_slice(&s.a_col[g]);
            }
            m
        });
        let repos_gather = std::array::from_fn(|t| {
            let mut gth = [0x80u8; 16];
            if t < N_GEO {
                for p in 0..12 {
                    gth[s.e_pos[t][p] as usize] = p as u8;
                }
            }
            gth
        });
        let relab_bits = std::array::from_fn(|g| if g < N_GEO { s.relab_bits[g] } else { [0; 16] });
        let flip_t: Vec<[u64; 64]> = (0..32)
            .map(|t| {
                if t < N_GEO {
                    bitset_flip(&s.e_pos[t])
                } else {
                    [0; 64]
                }
            })
            .collect();
        let c12 = comb_rank(12);
        let c9 = comb_rank(9);
        let c6 = comb_rank(6);
        let r12: Vec<u32> = c12.iter().map(|&x| x as u32).collect();
        let r9: Vec<u32> = c9.iter().map(|&x| 220 * x as u32).collect();
        let r6: Vec<u32> = c6.iter().map(|&x| 220 * 84 * x as u32).collect();
        FastTables {
            shuf_a,
            shuf_b,
            shuf_e,
            flip_m: boxed(&flip_m),
            corner: mt_c,
            c12: boxed(&c12),
            r12: boxed(&r12),
            r9: boxed(&r9),
            r6: r6.try_into().unwrap(),
            cmp6: boxed(&cmp6),
            relab_map,
            repos_gather,
            relab_bits,
            flip_t: boxed(&flip_t),
        }
    }

    pub fn of(&self, cube: &Cube) -> Fc {
        let mut e = [0u8; 12];
        let d = defs();
        for p in 0..12 {
            e[p] = d.edge_atri[cube.ep[p] as usize];
        }
        let mut f = Fc {
            a: [0; 16],
            b: [0; 16],
            e: [0; 16],
            c: corners().id(cube),
            ebits: edge_chirality(&cube.ep),
        };
        f.a[..12].copy_from_slice(&cube.xa);
        f.b[..12].copy_from_slice(&cube.xb);
        f.e[..12].copy_from_slice(&e);
        f
    }

    /// Apply move `m`. `em` are the edge-label masks of `x` (from [`masks`]).
    #[inline(always)]
    pub fn apply(&self, x: &Fc, em: &[u16; 4], m: usize) -> Fc {
        let fm = &self.flip_m[m & 15];
        let mut flips = 0u8;
        for t in 0..4 {
            let k = em[t] as usize & 0xfff;
            flips |= ((fm[k >> 6] >> (k & 63)) as u8 & 1) << t;
        }
        Fc {
            a: shuffle(&x.a, &self.shuf_a[m]),
            b: shuffle(&x.b, &self.shuf_b[m]),
            e: shuffle(&x.e, &self.shuf_e[m]),
            c: self.corner[x.c as usize][m],
            ebits: x.ebits ^ flips,
        }
    }

    /// Bits of `x` at the set positions of the 12-bit `mask`, packed (x must be a subset).
    #[inline(always)]
    fn compress(&self, x: u32, mask: u32) -> u32 {
        let lo = (mask & 63) as usize;
        let hi = ((mask >> 6) & 63) as usize;
        let a = self.cmp6[lo][(x & 63) as usize] as u32;
        let b = self.cmp6[hi][((x >> 6) & 63) as usize] as u32;
        a | b << (lo as u32).count_ones()
    }

    /// Multiset rank (same numbering as [`Mset::rank`]) from the label masks.
    #[inline(always)]
    pub fn rank(&self, m: &[u16; 4]) -> u32 {
        let (m0, m1, m2) = (m[0] as u32 & 0xfff, m[1] as u32, m[2] as u32);
        let free = 0xfff & !m0;
        let free2 = free & !m1;
        self.r12[m0 as usize]
            + self.r9[self.compress(m1, free) as usize & 511]
            + self.r6[self.compress(m2, free2) as usize & 63]
    }

    #[inline(always)]
    pub fn a_rank(&self, x: &Fc) -> u32 {
        self.rank(&masks(&x.a))
    }
    #[inline(always)]
    pub fn b_rank(&self, x: &Fc) -> u32 {
        self.rank(&masks(&x.b))
    }
    #[inline(always)]
    pub fn e_rank(&self, x: &Fc) -> u32 {
        self.rank(&masks(&x.e))
    }

    /// Apply the geometric symmetry pair (g, t) to the edge labelling and chirality, returning
    /// (rank, chirality bits).
    #[inline(always)]
    pub fn sym_e(&self, x: &Fc, em: &[u16; 4], g: usize, t: usize) -> (u32, u8) {
        // Reposition flips depend only on the occupied positions, and relabelling is a bit
        // permutation applied after xor-ing constants, so the flips can be applied before
        // relabelling: relab_bits[g][bits ^ flips(t, positions of each original triple)].
        let ft = &self.flip_t[t & 31];
        let mut flips = 0u8;
        for k in 0..4 {
            let q = em[k] as usize & 0xfff;
            flips |= ((ft[q >> 6] >> (q & 63)) as u8 & 1) << k;
        }
        let b = self.relab_bits[g & 31][((x.ebits ^ flips) & 15) as usize];
        // label map: e1[p] = map[e[p]]; then gather positions
        let e2 = shuffle(
            &shuffle(&self.relab_map[g & 31], &x.e),
            &self.repos_gather[t & 31],
        );
        (self.rank(&masks(&e2)), b)
    }

    /// Relabel only (large tier).
    #[inline(always)]
    pub fn relabel_e(&self, x: &Fc, g: usize) -> (u32, u8) {
        let e2 = shuffle(&self.relab_map[g & 31], &x.e);
        (
            self.rank(&masks(&e2)),
            self.relab_bits[g & 31][(x.ebits & 15) as usize],
        )
    }
}

pub fn fast() -> &'static FastTables {
    static F: OnceLock<FastTables> = OnceLock::new();
    F.get_or_init(FastTables::new)
}

/// Byte gather: out[i] = v[s[i]] (0 where s[i] has the high bit set).
#[inline(always)]
pub fn shuffle(v: &[u8; 16], s: &[u8; 16]) -> [u8; 16] {
    #[cfg(all(target_arch = "x86_64", target_feature = "ssse3"))]
    // SAFETY: SSSE3 is enabled for this build; loads/stores are unaligned-safe
    unsafe {
        use core::arch::x86_64::*;
        let r = _mm_shuffle_epi8(
            _mm_loadu_si128(v.as_ptr() as *const __m128i),
            _mm_loadu_si128(s.as_ptr() as *const __m128i),
        );
        let mut o = [0u8; 16];
        _mm_storeu_si128(o.as_mut_ptr() as *mut __m128i, r);
        o
    }
    #[cfg(not(all(target_arch = "x86_64", target_feature = "ssse3")))]
    {
        std::array::from_fn(|i| {
            if s[i] & 0x80 != 0 {
                0
            } else {
                v[(s[i] & 15) as usize]
            }
        })
    }
}

/// Bit masks of the positions (0..12) holding labels 0, 1, 2, 3.
#[inline(always)]
pub fn masks(v: &[u8; 16]) -> [u16; 4] {
    #[cfg(all(target_arch = "x86_64", target_feature = "sse2"))]
    // SAFETY: SSE2 is baseline on x86_64
    unsafe {
        use core::arch::x86_64::*;
        let x = _mm_loadu_si128(v.as_ptr() as *const __m128i);
        std::array::from_fn(|l| {
            (_mm_movemask_epi8(_mm_cmpeq_epi8(x, _mm_set1_epi8(l as i8))) & 0xfff) as u16
        })
    }
    #[cfg(not(all(target_arch = "x86_64", target_feature = "sse2")))]
    {
        let mut m = [0u16; 4];
        for p in 0..12 {
            m[v[p] as usize] |= 1 << p;
        }
        m
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    #[test]
    fn fast_coords_match_rank_coords() {
        let f = fast();
        let s = syms();
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(31);
        for _ in 0..300 {
            let cube = Cube::random(&mut rng);
            let x = f.of(&cube);
            let r = Coords::of(&cube);
            assert_eq!(
                (f.a_rank(&x), f.b_rank(&x), f.e_rank(&x), x.c, x.ebits),
                (r.a, r.b, r.e, r.c, r.ebits)
            );
            let em = masks(&x.e);
            for m in 0..N_MOVES {
                let y = f.apply(&x, &em, m);
                assert_eq!(y, f.of(&cube.apply(m)));
            }
            for &(g, t) in &s.pairs {
                let (e, b) = f.sym_e(&x, &em, g as usize, t as usize);
                assert_eq!((e, b), s.apply_e(g as usize, t as usize, r.e, r.ebits));
            }
        }
    }
}
