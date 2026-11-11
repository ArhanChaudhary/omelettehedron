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
use crate::tables::{prefetch_at, CbTable, Phase2};
use rayon::prelude::*;
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, OnceLock};

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
pub(crate) fn conj_renorm_corner_with(
    c: u16,
    m: usize,
    wps: &[Vec<usize>],
) -> Option<(u16, usize)> {
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
            b: mset_labels()[b],
            e: [0; 16],
            c: self.rep[cls],
            ebits: (e1 % 2) as u8,
        };
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

/// Moves between CBE entries (all of triple 0), `N` moves (the face metric, or with slices).
///
/// Consecutive entries share their corner class and B colouring (groups of `N_E1`), so the
/// corner/B part of each neighbour's index is made once per group and move. A step then only
/// moves the triple: one byte shuffle (the move and the symmetry composed), a mask and two
/// small lookups. (CBE generation 38.1 -> 21.2 s, large tier 46.2 -> 29.5 s at 48T when added,
/// opt/gen.)
struct CbeGraph<'a, const N: usize> {
    ix: &'a CbeIndex,
    f: &'static FastTables,
    /// labelling of each triple-0 position set (0 = triple 0, 1 = other edges)
    e_labels: Vec<[u8; 16]>,
    /// per move and repositioning symmetry t: the move's gather followed by the gather of t
    comp: Box<[[[u8; 16]; 32]; N]>,
    /// per t: whether repositioning flips the chirality of triple 0, as a bitset over the
    /// triple's positions after repositioning
    flip_post: Box<[[u64; 64]; 32]>,
}

/// The corner/B part of the neighbours of a group of entries.
#[derive(Clone, Copy)]
struct CbeGroup<const N: usize> {
    /// per move: class * `N_ROW` + B rank * `N_E1` of the neighbour
    base: [u64; N],
    /// per move: repositioning symmetry of the neighbour
    t: [u8; N],
    /// bit m: relabelling flips the chirality bit of the neighbour by move m
    flip: u32,
}

#[derive(Clone, Copy)]
struct CbeNode<const N: usize> {
    g: CbeGroup<N>,
    e: [u8; 16],
    /// bit m: chirality bit of the neighbour by move m before repositioning
    flip: u32,
}

impl<'a, const N: usize> CbeGraph<'a, N> {
    fn new(ix: &'a CbeIndex) -> Self {
        let f = fast();
        let e_labels = ix
            .c12_unrank
            .iter()
            .map(|&q| std::array::from_fn(|p| u8::from(p < 12 && q >> p & 1 == 0)))
            .collect();
        let d = defs();
        let comp = Box::new(std::array::from_fn(|m| {
            std::array::from_fn(|t| {
                let mut g = [0x80u8; 16];
                if t < N_GEO {
                    for p in 0..12 {
                        g[p] = d.moves[m].ep[ix.egather[t][p] as usize];
                    }
                }
                g
            })
        }));
        let mut flip_post = Box::new([[0u64; 64]; 32]);
        for t in 0..N_GEO {
            for q in (0..4096usize).filter(|q| q.count_ones() == 3) {
                let l: [u8; 16] = std::array::from_fn(|p| u8::from(p >= 12 || q >> p & 1 == 0));
                let post = label_mask(&shuffle(&l, &ix.egather[t]), 0) as usize & 0xfff;
                if f.flip_t[t][q >> 6] >> (q & 63) & 1 == 1 {
                    flip_post[t][post >> 6] |= 1 << (post & 63);
                }
            }
        }
        CbeGraph {
            ix,
            f,
            e_labels,
            comp,
            flip_post,
        }
    }

    /// Group `key` = class * `N_MSET` + B rank.
    fn group(&self, key: usize) -> CbeGroup<N> {
        let (cls, b) = (key / N_MSET, key % N_MSET);
        let x = Fc {
            a: [0; 16],
            b: mset_labels()[b],
            e: [0; 16],
            c: self.ix.rep[cls],
            ebits: 0,
        };
        let mut g = CbeGroup {
            base: [0; N],
            t: [0; N],
            flip: 0,
        };
        for m in 0..N {
            let y = self.f.apply_any::<N>(&x, 0, 0, m);
            let (cls2, m2, t2) = self.ix.kt[y.c as usize & 1023][0];
            let (m2, t2) = (m2 as usize & 31, t2 as usize & 31);
            let bb = shuffle(&shuffle(&self.ix.bmap[m2], &y.b), &self.ix.bgather[t2]);
            g.base[m] = (cls2 as usize * N_ROW + self.f.rank(&masks(&bb)) as usize * N_E1) as u64;
            g.t[m] = t2 as u8;
            g.flip |= u32::from(self.ix.relab_flip[m2][0]) << m;
        }
        g
    }
}

impl<const N: usize> Graph for CbeGraph<'_, N> {
    const MOVES: usize = N;
    type Node = CbeNode<N>;
    /// the last group and its key
    type Cache = Option<(usize, CbeGroup<N>)>;
    #[inline(always)]
    fn node(&self, i: usize, cache: &mut Option<(usize, CbeGroup<N>)>) -> CbeNode<N> {
        let key = i / N_E1;
        let g = match cache {
            Some((k, g)) if *k == key => *g,
            _ => {
                let g = self.group(key);
                *cache = Some((key, g));
                g
            }
        };
        let e1 = i % N_E1;
        let q = self.ix.c12_unrank[e1 / 2] as usize & 0xfff;
        let bit = if e1 & 1 == 1 { u32::MAX } else { 0 };
        let mut moved = u32::from(self.f.flip_moves[q]);
        if N > N_MOVES {
            moved |= u32::from(slice_tables().flip[q]) << N_MOVES;
        }
        CbeNode {
            g,
            e: self.e_labels[e1 / 2],
            flip: g.flip ^ moved ^ bit,
        }
    }
    #[inline(always)]
    fn step(&self, x: &CbeNode<N>, m: usize) -> usize {
        let t = x.g.t[m] as usize & 31;
        let mask = label_mask(&shuffle(&x.e, &self.comp[m % N][t]), 0) as usize & 0xfff;
        let post = (self.flip_post[t][mask >> 6] >> (mask & 63)) as usize & 1;
        x.g.base[m] as usize + self.f.c12[mask] as usize * 2 + ((x.flip >> m) as usize & 1 ^ post)
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

    /// The table in the metric of `moves`.
    pub fn build(moves: MoveSet, verbose: bool) -> CbeTable {
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
        let name = format!("cbe{}", moves.table_suffix());
        let t = if moves.slice {
            let g = CbeGraph::<MAX_MOVES>::new(&ix);
            saturating(&g, n, &targets, 6, &name, verbose)
        } else {
            let g = CbeGraph::<N_MOVES>::new(&ix);
            saturating(&g, n, &targets, 6, &name, verbose)
        };
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
            prefetch_at(&self.codes, out[k] >> 5);
        }
    }

    /// (`CODE_LB` is increasing, so this is the bound of the largest code)
    #[inline(always)]
    pub fn bound(&self, idx: &[usize; 4]) -> u8 {
        CODE_LB[idx.iter().fold(0, |b, &i| b.max(self.code(i))) as usize & 3]
    }
}

/// Distance to the SOLVED projection on the CBE index (corners, B-centre colouring, and the
/// position set and chirality of one A-face's edge triple; the same index and symmetries as
/// [`CbeTable`], read through [`CbeTable::ix`]), one bit per entry: whether the distance is at
/// least [`CBS_MIN`]. The distance to solved is invariant under the class-preserving symmetries,
/// the re-fixing rotations and `SWAP_AB`, so one table serves all four A-triples of both sides,
/// eight lookups per position: a lower bound for optimal mode (see `search.rs`).
///
/// One bit is all the search uses: its existing bounds leave almost no node with fewer than 10
/// moves left for the whole solution, and reading the table where 11 or more are left cost more
/// than it saved (measured with 2-bit codes of distances 9 to 12: reading it with at most 11
/// moves left 0.900x / 0.933x time at 1T / 48T, at most 10 left 0.843x / 0.874x).
pub struct CbsTable {
    pub bits: HugeBuf<u64>,
    /// numbers of entries at each distance
    pub hist: Vec<u64>,
}

/// The bit of an entry of [`CbsTable`] is set if its distance to solved is at least this.
pub const CBS_MIN: u8 = 11;

impl CbsTable {
    /// Whether entry `i` is at least [`CBS_MIN`] moves from solved.
    #[inline(always)]
    pub fn far(&self, i: usize) -> bool {
        (self.bits[i >> 6] >> (i & 63)) & 1 == 1
    }

    /// The indices of the solved projection (with its equivalents).
    pub fn targets(ix: &CbeIndex) -> Vec<usize> {
        let f = fast();
        let x = f.of(&Cube::solved());
        let mut v = vec![ix.index(f, &x, &masks(&x.e), 0)];
        ix.equivalents(f, v[0], &mut v);
        v.sort_unstable();
        v.dedup();
        v
    }

    /// The table (face turns only): saturating 2-bit codes from a breadth-first search (code
    /// 2 or 3 = distance at least `CBS_MIN`), packed into bits. Peak memory: 1.5 times the
    /// table's 2-bit size (7.0 + 3.5 GB).
    pub fn build(ix: &CbeIndex, verbose: bool) -> CbsTable {
        let g = CbeGraph::<N_MOVES>::new(ix);
        let t = saturating(
            &g,
            ix.size(),
            &Self::targets(ix),
            CBS_MIN - 2,
            "cbe-solved",
            verbose,
        );
        let mut bits = HugeBuf::<u64>::for_overwrite(ix.size().div_ceil(64));
        let codes = &t.data;
        bits.par_iter_mut().enumerate().for_each(|(w, b)| {
            // the high bits of the 2-bit codes of entries 64w .. 64w + 64
            let half = |j: usize| -> u64 {
                let x = codes.get(2 * w + j).copied().unwrap_or(0) >> 1 & 0x5555_5555_5555_5555;
                unspread(x)
            };
            *b = half(0) | half(1) << 32;
        });
        // (entries past the end hold code 3 in the BFS array: clear their bits)
        let n = ix.size();
        if n % 64 != 0 {
            if let Some(last) = bits.last_mut() {
                *last &= (1u64 << (n % 64)) - 1;
            }
        }
        CbsTable { bits, hist: t.hist }
    }

    /// Number of 2-bit slots the file header records (the bits are stored as a 2-bit table of
    /// half as many entries, so that [`Table2`]'s file format and checksum apply).
    fn slots(ix: &CbeIndex) -> usize {
        ix.size().div_ceil(2)
    }

    /// # Errors
    ///
    /// If the file cannot be written.
    pub fn save(&self, ix: &CbeIndex, path: &Path) -> std::io::Result<()> {
        save_words(path, Self::slots(ix), &self.hist, &self.bits)
    }

    /// # Errors
    ///
    /// If the file cannot be read or does not hold a table of the expected size.
    pub fn load(ix: &CbeIndex, path: &Path) -> std::io::Result<CbsTable> {
        let t = Table2::load(path, Self::slots(ix))?;
        Ok(CbsTable {
            bits: t.data,
            hist: t.hist,
        })
    }

    /// The four triple indices of `x` (into `out`), each word prefetched.
    #[inline(always)]
    pub fn prefetch_all(&self, ix: &CbeIndex, f: &FastTables, x: &Fc, out: &mut [usize; 4]) {
        let em = masks(&x.e);
        for k in 0..4 {
            out[k] = ix.index(f, x, &em, k);
            prefetch_at(&self.bits, out[k] >> 6);
        }
    }

    /// Whether some of the four triples (indices `idx`) is at least [`CBS_MIN`] from solved.
    #[inline(always)]
    pub fn any_far(&self, idx: &[usize; 4]) -> bool {
        idx.iter().any(|&i| self.far(i))
    }

    /// Check that the bits fit this layout: positions up to `CBS_MIN - 2` random moves from
    /// solved and all their neighbours are within `CBS_MIN - 1` (bit clear), and about the
    /// fraction of random positions given by the histogram is at least `CBS_MIN` away. `false`
    /// if a check fails (a file made for another layout).
    pub fn check_walks(&self, ix: &CbeIndex) -> bool {
        use rand::{Rng, SeedableRng};
        let f = fast();
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(0x50_1ed);
        let any_far = |x: &Fc| -> bool {
            let em = masks(&x.e);
            (0..4).any(|k| self.far(ix.index(f, x, &em, k)))
        };
        for _ in 0..16 {
            let mut x = f.of(&Cube::solved());
            for _ in 0..rng.gen_range(0..=CBS_MIN - 2) {
                x = f.apply(&x, &masks(&x.e), rng.gen_range(0..N_MOVES));
            }
            let em = masks(&x.e);
            if any_far(&x) || (0..N_MOVES).any(|m| any_far(&f.apply(&x, &em, m))) {
                return false;
            }
        }
        let total: u64 = self.hist.iter().sum();
        let far: u64 = self.hist.iter().skip(usize::from(CBS_MIN)).sum();
        if total == 0 {
            return false;
        }
        let want = far as f64 / total as f64;
        let samples = 4000;
        let got = (0..samples)
            .filter(|_| {
                let x = f.of(&Cube::random(&mut rng));
                self.far(ix.index(f, &x, &masks(&x.e), rng.gen_range(0..4)))
            })
            .count() as f64
            / f64::from(samples);
        (got - want).abs() < 0.05
    }
}

/// A tier's [`CbsTable`], read from its file only when a search can use it: optimal mode, once
/// some child may have fewer than `CBS_MIN` moves left (see `search.rs`). Time-limited solves
/// and proofs that end before that never read its 3.5 GB (start-up of the large tier 1.71 ->
/// 2.13 s at 1T, 0.27 -> 0.34 s at 48T when it was read with the other tables). The first such
/// search starts reading it on a background thread and goes on without the bound; the
/// workers use it from the moment it is ready, and later searches from their start.
pub struct LazyCbs {
    path: PathBuf,
    cell: OnceLock<Option<CbsTable>>,
    started: AtomicBool,
}

impl LazyCbs {
    /// A table that is read later from `path`.
    pub fn new(path: PathBuf) -> LazyCbs {
        LazyCbs {
            path,
            cell: OnceLock::new(),
            started: AtomicBool::new(false),
        }
    }

    /// A table already in memory (just generated).
    pub fn with(path: PathBuf, t: CbsTable) -> LazyCbs {
        let l = LazyCbs::new(path);
        let _ = l.cell.set(Some(t));
        l.started.store(true, Ordering::Relaxed);
        l
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    /// The table if it has been read (`None` while it is being read, or if it failed its
    /// checks).
    #[inline]
    pub fn loaded(&self) -> Option<&CbsTable> {
        self.cell.get().and_then(Option::as_ref)
    }

    /// The table, read now if no one has yet (another caller waits for the first). `None` if the
    /// file cannot be read or fails its checks; a message on stderr says why, and the bound is
    /// then off (the searches stay correct without it).
    pub fn get(&self, ix: &CbeIndex) -> Option<&CbsTable> {
        self.started.store(true, Ordering::Relaxed);
        self.cell
            .get_or_init(|| {
                let fail = |why: String| {
                    eprintln!(
                        "{}: {why}; the distance-to-solved bound is off (delete the file to \
                         regenerate it)",
                        self.path.display()
                    );
                    None
                };
                match CbsTable::load(ix, &self.path) {
                    Ok(t) if t.check_walks(ix) => Some(t),
                    Ok(_) => {
                        fail("not a distance-to-solved table for this version's layout".into())
                    }
                    Err(e) => fail(e.to_string()),
                }
            })
            .as_ref()
    }

    /// Start reading the table on a background thread (one thread, so that it takes one core
    /// from the search for the 0.4 s or so it needs), unless that has started already.
    pub fn start(self: &Arc<Self>) {
        if self.started.swap(true, Ordering::Relaxed) {
            return;
        }
        let me = Arc::clone(self);
        let spawned = std::thread::Builder::new()
            .name("cbe-solved".into())
            .spawn(move || {
                let read = || {
                    let ix = CbeIndex::new();
                    me.get(&ix).is_some()
                };
                // (its parallel loops on one thread of their own, not on the search's pool)
                match rayon::ThreadPoolBuilder::new().num_threads(1).build() {
                    Ok(p) => p.install(read),
                    Err(_) => read(),
                }
            });
        if spawned.is_err() {
            // no threads (e.g. WebAssembly without them): the bound stays off
            self.started.store(true, Ordering::Relaxed);
        }
    }
}

/// The even bits of `x` (bits 0, 2, .., 62; the odd bits must be clear) packed into bits 0..32.
fn unspread(mut x: u64) -> u64 {
    x = (x | x >> 1) & 0x3333_3333_3333_3333;
    x = (x | x >> 2) & 0x0f0f_0f0f_0f0f_0f0f;
    x = (x | x >> 4) & 0x00ff_00ff_00ff_00ff;
    x = (x | x >> 8) & 0x0000_ffff_0000_ffff;
    (x | x >> 16) & 0x0000_0000_ffff_ffff
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    /// The grouped neighbour indices of the generator equal the index of each moved position.
    #[test]
    fn graph_steps_match_index() {
        check_graph::<N_MOVES>();
        check_graph::<MAX_MOVES>();
    }

    fn check_graph<const N: usize>() {
        use rand::Rng;
        let f = fast();
        let ix = CbeIndex::new();
        let g = CbeGraph::<N>::new(&ix);
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(43);
        let mut cache = None;
        for _ in 0..300 {
            let i0 = rng.gen_range(0..ix.size() - 40);
            for i in i0..i0 + 40 {
                let x = ix.decode(i);
                let em = masks(&x.e);
                let node = g.node(i, &mut cache);
                for m in 0..N {
                    let y = f.apply_m(&x, &em, m);
                    assert_eq!(
                        g.step(&node, m),
                        ix.index(f, &y, &masks(&y.e), 0),
                        "{i} {m}"
                    );
                }
            }
        }
    }

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

    /// A missing or broken distance-to-solved file turns the bound off (`None`) instead of
    /// failing, read at once or on the background thread.
    #[test]
    fn lazy_table_failures_turn_the_bound_off() {
        let ix = CbeIndex::new();
        let dir = std::env::temp_dir().join(format!("fto-lazy-cbs-{:016x}", rand::random::<u64>()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("cbe-solved11.prun");
        let l = LazyCbs::new(p.clone());
        assert!(l.loaded().is_none());
        assert!(l.get(&ix).is_none());
        std::fs::write(&p, b"FTOPRUN2 not a table").unwrap();
        assert!(LazyCbs::new(p.clone()).get(&ix).is_none());
        // (WebAssembly test runs have no threads: `start` then leaves the bound off)
        if cfg!(not(target_family = "wasm")) {
            let l = Arc::new(LazyCbs::new(p));
            l.start();
            l.start();
            let t0 = std::time::Instant::now();
            while l.cell.get().is_none() && t0.elapsed().as_secs() < 120 {
                std::thread::sleep(std::time::Duration::from_millis(20));
            }
            assert!(l.cell.get().is_some(), "background read did not finish");
            assert!(l.loaded().is_none());
        }
        std::fs::remove_dir_all(&dir).unwrap();
    }

    /// `unspread` packs the even bits.
    #[test]
    fn unspread_packs_even_bits() {
        use rand::Rng;
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(5);
        for _ in 0..1000 {
            let x: u64 = rng.gen::<u64>() & 0x5555_5555_5555_5555;
            let want = (0..32).fold(0u64, |a, j| a | (x >> (2 * j) & 1) << j);
            assert_eq!(unspread(x), want);
        }
    }

    /// Exact distances to the solved projection within `r` moves (breadth-first, on the
    /// table's indices, equivalents included).
    fn solved_ball(ix: &CbeIndex, r: u8) -> HashMap<usize, u8> {
        let g = CbeGraph::<N_MOVES>::new(ix);
        let f = fast();
        let mut d: HashMap<usize, u8> = CbsTable::targets(ix).into_iter().map(|i| (i, 0)).collect();
        let mut front: Vec<usize> = d.keys().copied().collect();
        for k in 1..=r {
            let mut next = vec![];
            let mut cache = None;
            for &i in &front {
                let nd = g.node(i, &mut cache);
                for m in 0..N_MOVES {
                    let j = g.step(&nd, m);
                    if !d.contains_key(&j) {
                        let mut v = vec![j];
                        ix.equivalents(f, j, &mut v);
                        for e in v {
                            if d.insert(e, k).is_none() {
                                next.push(e);
                            }
                        }
                    }
                }
            }
            front = next;
        }
        d
    }

    /// Whether the projection of triple `k` of `x` is within `t` moves of solved, from a ball of
    /// radius `r` and a search over canonical sequences above it.
    fn within(
        ix: &CbeIndex,
        ball: &HashMap<usize, u8>,
        r: u8,
        x: &Fc,
        k: usize,
        t: u8,
        last: usize,
    ) -> bool {
        let f = fast();
        let em = masks(&x.e);
        if let Some(&d) = ball.get(&ix.index(f, x, &em, k)) {
            return d <= t;
        }
        t > r
            && (0..N_MOVES).any(|m| {
                allowed_after(false, last, m)
                    && within(ix, ball, r, &f.apply(x, &em, m), k, t - 1, m)
            })
    }

    /// The distance-to-solved bits equal exact distances found by search (positions a random
    /// walk of 8 to 14 moves from solved, all four triples). Needs the table
    /// (`FTO_TABLES` holding `cbe-solved11.prun`, 3.5 GB).
    #[test]
    #[ignore = "needs the 3.5 GB distance-to-solved table in FTO_TABLES"]
    fn to_solved_codes_match_search() {
        use rand::Rng;
        let ix = CbeIndex::new();
        let dir = crate::tables::table_dir(None);
        let t = CbsTable::load(&ix, &dir.join(format!("cbe-solved{CBS_MIN}.prun"))).unwrap();
        assert!(t.check_walks(&ix));
        let r = 7;
        let ball = solved_ball(&ix, r);
        let f = fast();
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(47);
        let mut seen = [0u32; 2];
        for _ in 0..60 {
            let mut x = f.of(&Cube::solved());
            for _ in 0..rng.gen_range(8..=14) {
                x = f.apply(&x, &masks(&x.e), rng.gen_range(0..N_MOVES));
            }
            let em = masks(&x.e);
            for k in 0..4 {
                let far = t.far(ix.index(f, &x, &em, k));
                let want = !within(&ix, &ball, r, &x, k, CBS_MIN - 1, usize::MAX);
                assert_eq!(far, want, "triple {k}");
                seen[usize::from(far)] += 1;
            }
        }
        // (both values occur)
        assert!(seen.iter().all(|&n| n > 0), "{seen:?}");
    }
}
