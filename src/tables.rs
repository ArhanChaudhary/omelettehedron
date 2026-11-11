//! Pruning tables for phase 1 and the exact distance table for phase 2.
//!
//! Phase 1 must reach `H_A` = <U, F, BR, BL>. Membership means, up to a whole-puzzle rotation:
//! (AE) A-centres solved and every edge on its own A-face in the right cyclic order, and
//! (CB) every B-centre has the colour of the adjacent corner sticker. The two parts get
//! separate tables; the heuristic is their maximum.

use crate::coord::*;
use crate::cube::*;
use crate::fast::*;
use crate::prune::*;
use crate::sym::*;
use rayon::prelude::*;
use std::collections::HashMap;
use std::path::{Path, PathBuf};

#[derive(
    Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum, serde::Serialize, serde::Deserialize,
)]
pub enum Tier {
    /// ~135 MB: edge triples without chirality, 288-fold symmetry
    Small,
    /// ~1.0 GB: full edge coset, 288-fold symmetry
    Medium,
    /// medium + the joint corner/B-centre/edge-triple table (~8 GB)
    Large,
    /// relabel-only edge-coset table (~23 GB) + the joint table (~30 GB in total)
    Huge,
}

/// The tables of a tier (see [`Tables::plan`]), with their sizes in memory.
pub fn tier_plan(dir: &Path, tier: Tier, moves: MoveSet) -> Vec<crate::preflight::Planned> {
    match tier {
        Tier::Small => Tables::plan(dir, &AeSym::new(false), false, moves),
        Tier::Medium => Tables::plan(dir, &AeSym::new(true), false, moves),
        Tier::Large => Tables::plan(dir, &AeSym::new(true), true, moves),
        Tier::Huge => Tables::plan(dir, &AeRelabel::new(), true, moves),
    }
}

/// Index function of an AE table.
pub trait AeIndex: Send + Sync + 'static {
    fn size(&self) -> usize;
    fn index(&self, x: &Coords) -> usize;
    /// A position (A rank, triple rank, chirality bits) with the given index.
    fn decode(&self, i: usize) -> (u32, u32, u8);
    fn name(&self) -> &'static str;
    /// Append every other index that denotes a position equivalent to index `i`.
    fn equivalents(&self, _i: usize, _out: &mut Vec<usize>) {}
    /// `index` over a batch, software-pipelined so that the dependent table lookups of
    /// different positions overlap.
    #[inline(always)]
    fn index_batch(&self, xs: &[Coords], out: &mut [usize]) {
        for (x, o) in xs.iter().zip(out.iter_mut()) {
            *o = self.index(x);
        }
    }
    /// Index of a fast-coordinate position, given its A rank.
    #[inline(always)]
    fn index_fast(&self, f: &FastTables, x: &Fc, a_rank: u32) -> usize {
        self.index_entry(f, x, *self.class_entry(a_rank))
    }
    /// Index of a fast-coordinate position, given the symmetry-class table entry of its A
    /// rank (only the edges of `x` are read).
    fn index_entry(&self, f: &FastTables, x: &Fc, v: u32) -> usize;
    /// Entries per A-centre class (the index is class * `class_size()` + edge part).
    fn class_size(&self) -> usize;
    /// [`AeIndex::index_entry`] in parts, for the generator: `(base, g, flips)` such that the
    /// index is `base + (f.relab(g, x.ebits ^ flips) & 7)` whatever the chirality bits of `x`
    /// are; `None` if the index holds no chirality bits.
    fn entry_parts(&self, f: &FastTables, x: &Fc, v: u32) -> Option<(usize, usize, u8)>;
    /// The symmetry-class table read first by `index_fast` (for prefetching).
    fn class_entry(&self, a_rank: u32) -> &u32;
    /// Prefetch the symmetry-class entry read first by `index_fast`.
    fn prefetch_class(&self, a_rank: u32);
    /// Pipelined `index_fast` over a batch of at most `N` positions.
    #[inline(always)]
    fn index_batch_fast<const N: usize>(&self, f: &FastTables, xs: &[Fc], out: &mut [usize]) {
        let n = xs.len();
        let mut ar = [0u32; N];
        for i in 0..n {
            ar[i] = f.a_rank(&xs[i]);
            self.prefetch_class(ar[i]);
        }
        for i in 0..n {
            out[i] = self.index_fast(f, &xs[i], ar[i]);
        }
    }
    /// `index_batch_fast` that prefetches each table word as soon as its index is known, so
    /// the first lookups' latency overlaps the index computation of the later positions.
    #[inline(always)]
    fn index_batch_fast_pf<const N: usize>(
        &self,
        f: &FastTables,
        xs: &[Fc],
        out: &mut [usize],
        tab: &[u64],
    ) {
        let n = xs.len();
        let mut ar = [0u32; N];
        for i in 0..n {
            ar[i] = f.a_rank(&xs[i]);
            self.prefetch_class(ar[i]);
        }
        for i in 0..n {
            out[i] = self.index_fast(f, &xs[i], ar[i]);
            prefetch_at(tab, out[i] >> 5);
        }
    }
    /// [`AeIndex::index_batch_fast_pf`] over the positions `xs[i]` for `i` in `sel`, writing
    /// `out[i]` (`N` moves; the `% N` keep the indices in bounds without checks, and are
    /// `& 15` for 16 moves).
    #[inline(always)]
    fn index_sel_pf<const N: usize>(
        &self,
        f: &FastTables,
        xs: &[Fc; N],
        sel: &[u8],
        out: &mut [usize; N],
        tab: &[u64],
    ) {
        let mut ar = [0u32; N];
        for (j, &i) in sel.iter().enumerate() {
            ar[j % N] = f.a_rank(&xs[usize::from(i) % N]);
            self.prefetch_class(ar[j % N]);
        }
        for (j, &i) in sel.iter().enumerate() {
            let i = usize::from(i) % N;
            out[i] = self.index_fast(f, &xs[i], ar[j % N]);
            prefetch_at(tab, out[i] >> 5);
        }
    }
}

/// Hint that `p` will be read soon.
#[inline(always)]
pub fn prefetch<T>(p: &T) {
    // WebAssembly has no prefetch, but a load whose value is dropped starts the cache miss
    // just the same while execution runs on. `p` is a reference, so the byte is always
    // readable.
    #[cfg(target_arch = "wasm32")]
    // SAFETY: `p` is a valid reference, so its first byte can be read
    unsafe {
        std::ptr::read_volatile(std::ptr::from_ref(p).cast::<u8>());
    }
    #[cfg(not(target_arch = "wasm32"))]
    prefetch_ptr(std::ptr::from_ref(p));
}

/// Hint that `s[i]` will be read soon. Where prefetching is a hint there is no bounds check:
/// prefetching outside `s` has no effect on the program (with the single-lookup CBE bound:
/// +0.8% nodes/s at 1T in optimal mode, when added). On WebAssembly the prefetch is a real
/// load, so an index outside `s` is skipped.
#[inline(always)]
pub fn prefetch_at<T>(s: &[T], i: usize) {
    #[cfg(target_arch = "wasm32")]
    if let Some(x) = s.get(i) {
        prefetch(x);
    }
    #[cfg(not(target_arch = "wasm32"))]
    prefetch_ptr(s.as_ptr().wrapping_add(i));
}

/// A prefetch instruction where there is one (never faults, whatever the address).
#[cfg(not(target_arch = "wasm32"))]
#[inline(always)]
fn prefetch_ptr<T>(p: *const T) {
    #[cfg(target_arch = "x86_64")]
    // SAFETY: prefetch is a hint and never faults, whatever the address
    unsafe {
        core::arch::x86_64::_mm_prefetch(p.cast::<i8>(), core::arch::x86_64::_MM_HINT_T0);
    }
    #[cfg(target_arch = "aarch64")]
    // SAFETY: PRFM is a hint and never faults, whatever the address
    unsafe {
        core::arch::asm!(
            "prfm pldl1keep, [{0}]",
            in(reg) p,
            options(nostack, preserves_flags, readonly)
        );
    }
    #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
    let _ = p;
}

pub struct AeSym {
    /// `cls.sym` holds class id | g << 16 | t << 24 (the symmetry pair itself, not its index)
    pub cls: AClasses,
    pub chiral: bool,
    /// symmetry pairs (other than the identity) fixing each class representative
    pub stab: Vec<Vec<(u8, u8)>>,
    /// the symmetry and move tables, held here so that the per-entry methods read no
    /// `OnceLock` (see [`AeGraph`])
    s: &'static Syms,
    mt: &'static MoveTables,
}

impl AeSym {
    pub fn new(chiral: bool) -> AeSym {
        let s = syms();
        let mut cls = AClasses::new(&s.pairs);
        for v in cls.sym.iter_mut() {
            let (g, t) = s.pairs[(*v >> 16) as usize];
            *v = (*v & 0xffff) | (g as u32) << 16 | (t as u32) << 24;
        }
        let ms = mset();
        let stab = cls
            .rep
            .iter()
            .map(|&r| {
                let l = ms.unrank(r);
                s.pairs[1..]
                    .iter()
                    .copied()
                    .filter(|&(g, t)| s.apply_a(g as usize, t as usize, &l) == l)
                    .collect()
            })
            .collect();
        AeSym {
            cls,
            chiral,
            stab,
            s,
            mt: move_tables(),
        }
    }
}

impl AeIndex for AeSym {
    fn size(&self) -> usize {
        self.cls.rep.len() * if self.chiral { N_ECOSET } else { N_MSET }
    }
    #[inline(always)]
    fn index(&self, x: &Coords) -> usize {
        let s = self.s;
        let v = self.cls.sym[x.a as usize];
        let (g, t) = (v >> 16 & 0xff, v >> 24);
        let (e, b) = s.apply_e(g as usize, t as usize, x.e, x.ebits);
        let c = (v & 0xffff) as usize;
        if self.chiral {
            c * N_ECOSET + e as usize * 8 + (b & 7) as usize
        } else {
            c * N_MSET + e as usize
        }
    }
    #[inline(always)]
    fn index_batch(&self, xs: &[Coords], out: &mut [usize]) {
        let s = self.s;
        let n = xs.len();
        let mut v = [0u32; MAX_MOVES];
        let mut e1 = [0u32; MAX_MOVES];
        for i in 0..n {
            prefetch(&self.cls.sym[xs[i].a as usize]);
        }
        for i in 0..n {
            v[i] = self.cls.sym[xs[i].a as usize];
            let g = v[i] >> 16 & 0xff;
            prefetch(&s.relab_e[g as usize * N_MSET + xs[i].e as usize]);
        }
        for i in 0..n {
            let (g, t) = (v[i] >> 16 & 0xff, v[i] >> 24);
            e1[i] = s.relab_e[g as usize * N_MSET + xs[i].e as usize];
            prefetch(&s.repos_e[t as usize * N_MSET + e1[i] as usize]);
        }
        for i in 0..n {
            let (g, t) = (v[i] >> 16 & 0xff, v[i] >> 24);
            let b1 = s.relab_bits[g as usize][xs[i].ebits as usize];
            let p = s.repos_e[t as usize * N_MSET + e1[i] as usize];
            let (e, b) = (p & 0xfffff, b1 ^ (p >> 20) as u8);
            let c = (v[i] & 0xffff) as usize;
            out[i] = if self.chiral {
                c * N_ECOSET + e as usize * 8 + (b & 7) as usize
            } else {
                c * N_MSET + e as usize
            };
        }
    }
    fn decode(&self, i: usize) -> (u32, u32, u8) {
        let mt = self.mt;
        if self.chiral {
            let c = i / N_ECOSET;
            let r = i % N_ECOSET;
            let e = (r / 8) as u32;
            let b3 = (r % 8) as u8;
            let bit3 = mt.eref[e as usize] ^ (b3 ^ b3 >> 1 ^ b3 >> 2) & 1;
            (self.cls.rep[c], e, b3 | bit3 << 3)
        } else {
            let e = (i % N_MSET) as u32;
            (self.cls.rep[i / N_MSET], e, mt.eref[e as usize] << 3)
        }
    }
    #[inline(always)]
    fn index_entry(&self, f: &FastTables, x: &Fc, v: u32) -> usize {
        let (g, t) = (v >> 16 & 0xff, v >> 24);
        let (e, b) = f.sym_e(x, &masks(&x.e), g as usize, t as usize);
        let c = (v & 0xffff) as usize;
        if self.chiral {
            c * N_ECOSET + e as usize * 8 + (b & 7) as usize
        } else {
            c * N_MSET + e as usize
        }
    }
    fn class_size(&self) -> usize {
        if self.chiral {
            N_ECOSET
        } else {
            N_MSET
        }
    }
    #[inline(always)]
    fn entry_parts(&self, f: &FastTables, x: &Fc, v: u32) -> Option<(usize, usize, u8)> {
        if !self.chiral {
            return None;
        }
        let (g, t) = ((v >> 16 & 0xff) as usize, (v >> 24) as usize);
        let (e, flips) = f.sym_e_parts(x, &masks(&x.e), g, t);
        Some(((v & 0xffff) as usize * N_ECOSET + e as usize * 8, g, flips))
    }
    fn class_entry(&self, a_rank: u32) -> &u32 {
        &self.cls.sym[a_rank as usize]
    }
    #[inline(always)]
    fn prefetch_class(&self, a_rank: u32) {
        prefetch_at(&self.cls.sym, a_rank as usize);
    }
    fn name(&self) -> &'static str {
        if self.chiral {
            "ae-sym288-chiral"
        } else {
            "ae-sym288-triples"
        }
    }
    fn equivalents(&self, i: usize, out: &mut Vec<usize>) {
        let stride = if self.chiral { N_ECOSET } else { N_MSET };
        let c = i / stride;
        if self.stab[c].is_empty() {
            return;
        }
        let s = self.s;
        let (_, e, bits) = self.decode(i);
        for &(g, t) in &self.stab[c] {
            let (e2, b2) = s.apply_e(g as usize, t as usize, e, bits);
            out.push(if self.chiral {
                c * N_ECOSET + e2 as usize * 8 + (b2 & 7) as usize
            } else {
                c * N_MSET + e2 as usize
            });
        }
    }
}

/// Full edge coset, reduced only by colour relabelling with the 12 rotations.
pub struct AeRelabel {
    pub cls: AClasses,
    /// as in [`AeSym`]
    s: &'static Syms,
    mt: &'static MoveTables,
}

impl AeRelabel {
    pub fn new() -> AeRelabel {
        let pairs: Vec<(u8, u8)> = (0..N_GEO as u8)
            .filter(|&g| crate::geometry::det(&syms().mats[g as usize]) == 1)
            .map(|g| (g, 0))
            .collect();
        AeRelabel {
            cls: AClasses::new(&pairs),
            s: syms(),
            mt: move_tables(),
        }
    }
}

impl AeIndex for AeRelabel {
    fn size(&self) -> usize {
        self.cls.rep.len() * N_ECOSET
    }
    #[inline(always)]
    fn index(&self, x: &Coords) -> usize {
        let s = self.s;
        let v = self.cls.sym[x.a as usize];
        // pair index == rotation index in the identity-first list of rotations
        let g = ROT_IDX[(v >> 16) as usize] as usize;
        let e = s.relab_e[g * N_MSET + x.e as usize];
        let b = s.relab_bits[g][x.ebits as usize];
        (v & 0xffff) as usize * N_ECOSET + e as usize * 8 + (b & 7) as usize
    }
    #[inline(always)]
    fn index_batch(&self, xs: &[Coords], out: &mut [usize]) {
        let s = self.s;
        let n = xs.len();
        let mut v = [0u32; MAX_MOVES];
        for i in 0..n {
            prefetch(&self.cls.sym[xs[i].a as usize]);
        }
        for i in 0..n {
            v[i] = self.cls.sym[xs[i].a as usize];
            let g = ROT_IDX[(v[i] >> 16) as usize] as usize;
            prefetch(&s.relab_e[g * N_MSET + xs[i].e as usize]);
        }
        for i in 0..n {
            let g = ROT_IDX[(v[i] >> 16) as usize] as usize;
            let e = s.relab_e[g * N_MSET + xs[i].e as usize];
            let b = s.relab_bits[g][xs[i].ebits as usize];
            out[i] = (v[i] & 0xffff) as usize * N_ECOSET + e as usize * 8 + (b & 7) as usize;
        }
    }
    fn decode(&self, i: usize) -> (u32, u32, u8) {
        let mt = self.mt;
        let c = i / N_ECOSET;
        let r = i % N_ECOSET;
        let e = (r / 8) as u32;
        let b3 = (r % 8) as u8;
        let bit3 = mt.eref[e as usize] ^ (b3 ^ b3 >> 1 ^ b3 >> 2) & 1;
        (self.cls.rep[c], e, b3 | bit3 << 3)
    }
    #[inline(always)]
    fn index_entry(&self, f: &FastTables, x: &Fc, v: u32) -> usize {
        let g = ROT_IDX[(v >> 16) as usize] as usize;
        let (e, b) = f.relabel_e(x, g);
        (v & 0xffff) as usize * N_ECOSET + e as usize * 8 + (b & 7) as usize
    }
    fn class_size(&self) -> usize {
        N_ECOSET
    }
    #[inline(always)]
    fn entry_parts(&self, f: &FastTables, x: &Fc, v: u32) -> Option<(usize, usize, u8)> {
        let g = ROT_IDX[(v >> 16) as usize] as usize;
        let base = (v & 0xffff) as usize * N_ECOSET + f.relabel_rank(x, g) as usize * 8;
        Some((base, g, 0))
    }
    fn class_entry(&self, a_rank: u32) -> &u32 {
        &self.cls.sym[a_rank as usize]
    }
    #[inline(always)]
    fn prefetch_class(&self, a_rank: u32) {
        prefetch_at(&self.cls.sym, a_rank as usize);
    }
    fn name(&self) -> &'static str {
        "ae-relabel12-chiral"
    }
}

/// `class_preserving_group()` lists the identity first and the 12 rotations before reflections.
const ROT_IDX: [u8; 12] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11];

/// Exact phase-2 table over `H_A` (933120 positions = 960 corner states x 972 edge states).
pub struct Phase2 {
    pub edge_id: HashMap<[u8; 12], u16>,
    /// In `H_A` the A-centres follow from the edge arrangement and the B-centres from the corner
    /// state (checked over all 933,120 states): the centres of each edge id and corner id.
    edge_xa: Vec<[u8; 12]>,
    corner_xb: Vec<[u8; 12]>,
    pub dist: Vec<u8>,
    pub emoves: Vec<[u16; N_MOVES]>,
}

pub const N_P2_EDGE: usize = 972;

impl Phase2 {
    pub fn new() -> Phase2 {
        let s0 = Cube::solved();
        let amoves: Vec<usize> = (0..N_MOVES).filter(|&m| is_a_move(m)).collect();
        let mut edge_id: HashMap<[u8; 12], u16> = HashMap::new();
        // breadth-first over whole positions, keeping the first one seen per edge arrangement
        let mut ecubes = vec![s0];
        edge_id.insert(s0.ep, 0);
        let mut i = 0;
        while i < ecubes.len() {
            let c = ecubes[i];
            i += 1;
            for &m in &amoves {
                let n = c.apply(m);
                if let std::collections::hash_map::Entry::Vacant(e) = edge_id.entry(n.ep) {
                    e.insert(ecubes.len() as u16);
                    ecubes.push(n);
                }
            }
        }
        let elist: Vec<[u8; 12]> = ecubes.iter().map(|c| c.ep).collect();
        let edge_xa: Vec<[u8; 12]> = ecubes.iter().map(|c| c.xa).collect();
        // the same per corner state
        let mut corner_xb = vec![s0.xb; N_CORNER];
        let mut seen = vec![false; N_CORNER];
        seen[corners().id(&s0) as usize] = true;
        let mut ccubes = vec![s0];
        let mut i = 0;
        while i < ccubes.len() {
            let c = ccubes[i];
            i += 1;
            for &m in &amoves {
                let n = c.apply(m);
                let id = corners().id(&n) as usize;
                if !seen[id] {
                    seen[id] = true;
                    corner_xb[id] = n.xb;
                    ccubes.push(n);
                }
            }
        }
        debug_assert!(seen.iter().all(|&s| s));
        debug_assert_eq!(elist.len(), N_P2_EDGE);
        let emoves: Vec<[u16; N_MOVES]> = elist
            .iter()
            .map(|&ep| {
                let c = Cube { ep, ..s0 };
                std::array::from_fn(|m| {
                    if is_a_move(m) {
                        edge_id[&c.apply(m).ep]
                    } else {
                        u16::MAX
                    }
                })
            })
            .collect();
        let mt = move_tables();
        let n = N_CORNER * N_P2_EDGE;
        let mut dist = vec![u8::MAX; n];
        dist[0] = 0;
        let mut fr = vec![0usize];
        let mut d = 0u8;
        while !fr.is_empty() {
            let mut nf = vec![];
            for &x in &fr {
                let (c, e) = (x / N_P2_EDGE, x % N_P2_EDGE);
                for &m in &amoves {
                    let y = mt.c[c][m] as usize * N_P2_EDGE + emoves[e][m] as usize;
                    if dist[y] == u8::MAX {
                        dist[y] = d + 1;
                        nf.push(y);
                    }
                }
            }
            fr = nf;
            d += 1;
        }
        debug_assert!(dist.iter().all(|&x| x != u8::MAX));
        Phase2 {
            edge_id,
            edge_xa,
            corner_xb,
            dist,
            emoves,
        }
    }

    /// Index of a position in the phase-2 table, or None if it is not in `H_A`.
    #[inline]
    pub fn index(&self, cube: &Cube) -> Option<usize> {
        let e = *self.edge_id.get(&cube.ep)? as usize;
        let c = corners().id(cube) as usize;
        (cube.xa == self.edge_xa[e] && cube.xb == self.corner_xb[c]).then_some(c * N_P2_EDGE + e)
    }

    /// Phase-2 length of a position (None if it is not in `H_A`).
    #[inline]
    pub fn length(&self, cube: &Cube) -> Option<u8> {
        Some(self.dist[self.index(cube)?])
    }

    /// Optimal phase-2 move sequence (A moves only) for a position in `H_A`.
    pub fn solve(&self, cube: &Cube) -> Option<Vec<usize>> {
        let mut x = *cube;
        let id = |y: &Cube| -> Option<usize> {
            Some(corners().id(y) as usize * N_P2_EDGE + *self.edge_id.get(&y.ep)? as usize)
        };
        let mut d = self.dist[id(&x)?];
        let mut out = vec![];
        while d > 0 {
            let m = (0..N_MOVES)
                .filter(|&m| is_a_move(m))
                .find(|&m| id(&x.apply(m)).is_some_and(|i| self.dist[i] == d - 1))?;
            out.push(m);
            x = x.apply(m);
            d -= 1;
        }
        if x == Cube::solved() {
            Some(out)
        } else {
            None
        }
    }

    pub fn histogram(&self) -> Vec<u64> {
        let mut h = vec![0u64; 16];
        for &d in &self.dist {
            h[d as usize] += 1;
        }
        while h.last() == Some(&0) {
            h.pop();
        }
        h
    }
}

/// CB part: corner state x B-centre colouring; targets are the 960 states of `H_A`.
pub struct CbTable {
    pub tab: Table2,
    /// B colouring attached to each corner state inside `H_A`
    pub attach: Vec<u32>,
}

/// Entries of the CB table: 51 corner classes under the 24 class-preserving symmetries x the
/// B-centre colourings (see [`crate::cbsym`]).
pub const N_CB: usize = 51 * N_MSET;

impl CbTable {
    pub fn attach_targets() -> Vec<u32> {
        let mt = move_tables();
        let s0 = Cube::solved();
        let b0 = mset().rank(&s0.xb);
        let mut attach = vec![u32::MAX; N_CORNER];
        attach[0] = b0;
        let mut fr = vec![0usize];
        while !fr.is_empty() {
            let mut nf = vec![];
            for &c in &fr {
                for m in 0..N_MOVES {
                    if !is_a_move(m) {
                        continue;
                    }
                    let c2 = mt.c[c][m] as usize;
                    let b2 = mt.b[attach[c] as usize][m];
                    if attach[c2] == u32::MAX {
                        attach[c2] = b2;
                        nf.push(c2);
                    } else {
                        debug_assert_eq!(
                            attach[c2], b2,
                            "B colouring is not a function of the corners in H_A"
                        );
                    }
                }
            }
            fr = nf;
        }
        attach
    }

    /// The table in the metric of `moves`.
    pub fn build(moves: MoveSet, verbose: bool) -> Table2 {
        let attach = Self::attach_targets();
        let ix = crate::cbsym::cbsym();
        let f = fast();
        let mut targets: Vec<usize> = vec![];
        for (c, &b) in attach.iter().enumerate().take(N_CORNER) {
            let mut x = ix.decode(0);
            x.c = c as u16;
            x.b[..12].copy_from_slice(&mset().unrank(b));
            let i = ix.index(f, &x);
            targets.push(i);
            ix.equivalents(i, &mut targets);
        }
        targets.sort_unstable();
        targets.dedup();
        let name = format!("cb{}", moves.table_suffix());
        if moves.slice {
            let g = crate::cbsym::CbSymGraph::<MAX_MOVES>::default();
            bfs(&g, N_CB, &targets, &name, verbose)
        } else {
            let g = crate::cbsym::CbSymGraph::<N_MOVES>::default();
            bfs(&g, N_CB, &targets, &name, verbose)
        }
    }

    /// Index of a position given by its coordinates.
    pub fn index(x: &Coords) -> usize {
        let ix = crate::cbsym::cbsym();
        let mut y = ix.decode(0);
        y.c = x.c;
        y.b[..12].copy_from_slice(&mset().unrank(x.b));
        ix.index(fast(), &y)
    }

    /// Index of a position (reduced by the 24 class-preserving symmetries, so that the table
    /// is 4.7 MB instead of 89 MB and stays in the cache).
    #[inline(always)]
    pub fn index_fast(f: &FastTables, x: &Fc) -> usize {
        f.cbsym.index(f, x)
    }
}

/// All indices of the solved AE class.
pub fn ae_targets<I: AeIndex>(ix: &I) -> Vec<usize> {
    let t = ix.index(&Coords::of(&Cube::solved()));
    let mut v = vec![t];
    ix.equivalents(t, &mut v);
    v.sort_unstable();
    v.dedup();
    v
}

/// Moves between AE entries, with the cache-resident coordinates of the search: a move is a
/// byte shuffle and the index needs only small tables (the move and symmetry tables of
/// [`Coords`] are megabytes each). An entry's A-centres are its class representative, so the
/// A part of each neighbour (its class-table entry) is tabulated per class and move. (Fast
/// coordinates: small tier 2.2 -> 1.5 s, AE medium generation 9.4 -> 8.1 s at 48T; the
/// per-class table: AE medium 7.6 -> 6.1 s, medium tier 8.1 -> 6.7 s; when added, opt/gen.)
struct AeGraph<'a, I: AeIndex, const N: usize> {
    ix: &'a I,
    f: &'static FastTables,
    /// [`mset_labels`], held here: a `OnceLock` read per entry is an atomic load, cheap
    /// natively but slow with threads in WebAssembly on `JavaScriptCore`, where threads made
    /// generation slower than one thread. With this, [`AeSym`]'s and [`crate::cbsym`]'s held
    /// tables: small tier in Bun, 6 threads 48.6 -> 7.35 s; one thread and native unchanged
    /// (opt/wasm-gen).
    labels: &'static [[u8; 16]],
    /// class-table entry of each class representative after each move
    next: Vec<[u32; N]>,
    /// whether the index holds chirality bits (see [`AeIndex::entry_parts`])
    chiral: bool,
}

/// The neighbours of the 8 entries that differ only in their chirality bits (same A class and
/// edge labelling): the neighbour by move m of the entry with bits b is
/// `base[m] + (relab(g[m], b ^ flips[m]) & 7)`.
#[derive(Clone, Copy)]
struct AeGroup<const N: usize> {
    base: [usize; N],
    g: [u8; N],
    flips: [u8; N],
}

impl<'a, I: AeIndex, const N: usize> AeGraph<'a, I, N> {
    fn new(ix: &'a I) -> Self {
        let f = fast();
        // `decode` reads the move tables; their lazy initialisation runs parallel work, so it
        // must not start inside a rayon job (the initialising thread could steal a job that
        // waits for it).
        move_tables();
        let labels = mset_labels();
        let cs = ix.class_size();
        let next = (0..ix.size() / cs)
            .into_par_iter()
            .map(|c| {
                let x = Fc {
                    a: labels[ix.decode(c * cs).0 as usize],
                    b: [0; 16],
                    e: [0; 16],
                    c: 0,
                    ebits: 0,
                };
                std::array::from_fn(|m| *ix.class_entry(f.a_rank(&f.apply_any::<N>(&x, 0, 0, m))))
            })
            .collect();
        let zero = Fc {
            a: [0; 16],
            b: [0; 16],
            e: [0; 16],
            c: 0,
            ebits: 0,
        };
        let chiral = ix.entry_parts(f, &zero, *ix.class_entry(0)).is_some();
        AeGraph {
            ix,
            f,
            labels,
            next,
            chiral,
        }
    }

    /// The group of entries `8 * key .. 8 * key + 8` (chiral indices only).
    fn group(&self, key: usize) -> AeGroup<N> {
        let i0 = key * 8;
        let (_, e, _) = self.ix.decode(i0);
        let c = i0 / self.ix.class_size();
        let x = Fc {
            a: [0; 16],
            b: [0; 16],
            e: self.labels[e as usize],
            c: 0,
            ebits: 0,
        };
        let em = masks(&x.e);
        let fw = self.f.flip_word(&em);
        let fws = if N > N_MOVES {
            self.f.flip_word_slice(&em)
        } else {
            0
        };
        let mut g = AeGroup {
            base: [0; N],
            g: [0; N],
            flips: [0; N],
        };
        for m in 0..N {
            // (x has no chirality bits, so y has the flips of the move)
            let y = self.f.apply_any::<N>(&x, fw, fws, m);
            if let Some((base, gm, fl)) = self.ix.entry_parts(self.f, &y, self.next[c][m]) {
                g.base[m] = base;
                g.g[m] = gm as u8;
                g.flips[m] = y.ebits ^ fl;
            }
        }
        g
    }
}

impl<I: AeIndex, const N: usize> Graph for AeGraph<'_, I, N> {
    const MOVES: usize = N;
    /// edges (the A-centres are not needed), flip words (face moves, slices) and A class
    type Node = (Fc, u64, u32, usize);
    /// the last group of a forward step and its key
    type Cache = Option<(usize, AeGroup<N>)>;
    #[inline(always)]
    fn node(&self, i: usize, _: &mut Option<(usize, AeGroup<N>)>) -> (Fc, u64, u32, usize) {
        let (_, e, bits) = self.ix.decode(i);
        let x = Fc {
            a: [0; 16],
            b: [0; 16],
            e: self.labels[e as usize],
            c: 0,
            ebits: bits,
        };
        let em = masks(&x.e);
        let fws = if N > N_MOVES {
            self.f.flip_word_slice(&em)
        } else {
            0
        };
        (x, self.f.flip_word(&em), fws, i / self.ix.class_size())
    }
    #[inline(always)]
    fn step(&self, (x, fw, fws, c): &(Fc, u64, u32, usize), m: usize) -> usize {
        let y = self.f.apply_any::<N>(x, *fw, *fws, m);
        self.ix.index_entry(self.f, &y, self.next[*c][m])
    }
    fn equivalents(&self, j: usize, out: &mut Vec<usize>) {
        self.ix.equivalents(j, out);
    }
    /// For chiral indices the neighbours are made per group of 8 entries (most of the work
    /// is shared: only the chirality part differs). (Medium tier 6.9 -> 6.6 s at 48T, 185.7 ->
    /// 167.0 s at 1T when added, opt/gen.)
    #[inline(always)]
    fn neighbours(&self, i: usize, cache: &mut Self::Cache, mut out: impl FnMut(usize)) {
        if !self.chiral {
            let x = self.node(i, cache);
            for m in 0..N {
                out(self.step(&x, m));
            }
            return;
        }
        let key = i / 8;
        if !matches!(cache, Some((k, _)) if *k == key) {
            *cache = Some((key, self.group(key)));
        }
        let Some((_, g)) = cache else {
            return;
        };
        let bits = self.ix.decode(i).2;
        for m in 0..N {
            out(g.base[m] + (self.f.relab(g.g[m] as usize, bits ^ g.flips[m]) & 7) as usize);
        }
    }
}

/// The AE table of index `ix` in the metric of `moves`.
pub fn build_ae<I: AeIndex>(ix: &I, moves: MoveSet, verbose: bool) -> Table2 {
    let name = format!("{}{}", ix.name(), moves.table_suffix());
    let targets = ae_targets(ix);
    if moves.slice {
        let g = AeGraph::<I, MAX_MOVES>::new(ix);
        bfs(&g, ix.size(), &targets, &name, verbose)
    } else {
        let g = AeGraph::<I, N_MOVES>::new(ix);
        bfs(&g, ix.size(), &targets, &name, verbose)
    }
}

/// All tables needed by the solver, for one AE index type and move set.
pub struct Tables<I: AeIndex> {
    /// the move set: the phase-1 tables hold distances in its metric, and the search uses it
    pub moves: MoveSet,
    pub ae_ix: I,
    pub ae: Table2,
    pub ae_targets: Vec<usize>,
    pub cb: CbTable,
    pub p2: Phase2,
    pub cbe: Option<crate::cbe::CbeTable>,
    /// distances to solved on the joint table's index (large and huge tiers, face turns): a
    /// bound for optimal mode (see `search.rs`), read when first needed
    pub cbs: Option<std::sync::Arc<crate::cbe::LazyCbs>>,
}

pub fn table_dir(dir: Option<&Path>) -> PathBuf {
    dir.map(Path::to_path_buf)
        .or_else(|| std::env::var_os("FTO_TABLES").map(PathBuf::from))
        .unwrap_or_else(|| PathBuf::from("tables"))
}

fn load_or_build(
    path: &Path,
    n: usize,
    verbose: bool,
    build: impl FnOnce() -> Table2,
) -> std::io::Result<Table2> {
    match Table2::load(path, n) {
        Ok(t) => {
            if verbose {
                eprintln!(
                    "loaded {} ({:.1} MB)",
                    path.display(),
                    t.bytes() as f64 / 1e6
                );
            }
            Ok(t)
        }
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
            if let Some(dir) = path.parent() {
                std::fs::create_dir_all(dir)?;
            }
            if verbose {
                eprintln!("generating {} ...", path.display());
            }
            let t = build();
            t.save(path)?;
            Ok(t)
        }
        Err(e) => Err(e),
    }
}

impl<I: AeIndex> Tables<I> {
    /// Build all tables in memory without touching the disk.
    pub fn generate(ae_ix: I, moves: MoveSet, verbose: bool) -> Tables<I> {
        let ae = build_ae(&ae_ix, moves, verbose);
        let cb = CbTable {
            tab: CbTable::build(moves, verbose),
            attach: CbTable::attach_targets(),
        };
        let ae_targets = ae_targets(&ae_ix);
        Tables {
            moves,
            ae_ix,
            ae,
            ae_targets,
            cb,
            p2: Phase2::new(),
            cbe: None,
            cbs: None,
        }
    }

    /// [`Self::load_with`] without the joint table.
    ///
    /// # Errors
    ///
    /// As [`Self::load_with`].
    pub fn load(dir: &Path, ae_ix: I, moves: MoveSet, verbose: bool) -> std::io::Result<Tables<I>> {
        Self::load_with(dir, ae_ix, false, moves, verbose)
    }

    /// The tables [`Self::load_with`] reads or generates, with their sizes in memory (for the
    /// checks of [`crate::preflight`]). 2 bits per entry.
    /// (generation needs no memory beyond the tables themselves, except the 7.0 GB search
    /// array of the 3.5 GB distance-to-solved table: the large tier peaks at 18.6 GB for
    /// 11.6 GB of tables)
    pub fn plan(
        dir: &Path,
        ae_ix: &I,
        with_cbe: bool,
        moves: MoveSet,
    ) -> Vec<crate::preflight::Planned> {
        let sfx = moves.table_suffix();
        let ae_path = dir.join(format!("{}{sfx}.prun", ae_ix.name()));
        let cb_path = dir.join(format!("cb-sym24{sfx}.prun"));
        let cbe_path = dir.join(format!("cbe{sfx}.prun"));
        let with_cbs = with_cbe && !moves.slice;
        let cbs_path = dir.join(format!("cbe-solved{}.prun", crate::cbe::CBS_MIN));
        let planned = |path: PathBuf, size: &dyn Fn() -> usize| {
            let bytes =
                std::fs::metadata(&path).map_or_else(|_| size().div_ceil(4) as u64, |m| m.len());
            crate::preflight::Planned {
                path,
                bytes,
                working_set: 0,
            }
        };
        let mut plan = vec![
            planned(ae_path, &|| ae_ix.size()),
            planned(cb_path, &|| N_CB),
        ];
        if with_cbe {
            plan.push(planned(cbe_path, &|| crate::cbe::CbeIndex::new().size()));
        }
        if with_cbs {
            // (one bit per entry, made from a 2-bit breadth-first search array)
            let n = || crate::cbe::CbeIndex::new().size();
            let exists = std::fs::metadata(&cbs_path).is_ok();
            let mut p = planned(cbs_path, &|| n().div_ceil(2));
            if !exists {
                p.working_set = n().div_ceil(4) as u64;
            }
            plan.push(p);
        }
        plan
    }

    /// Load (generating if needed) the tables, optionally with the joint CBE table, in the
    /// metric of `moves` (with slices, the tables have their own files, `*-slice.prun`).
    ///
    /// # Errors
    ///
    /// If a table file cannot be read or written, or holds a table of the wrong format, size
    /// or checksum.
    pub fn load_with(
        dir: &Path,
        ae_ix: I,
        with_cbe: bool,
        moves: MoveSet,
        verbose: bool,
    ) -> std::io::Result<Tables<I>> {
        let sfx = moves.table_suffix();
        let ae_path = dir.join(format!("{}{sfx}.prun", ae_ix.name()));
        let cb_path = dir.join(format!("cb-sym24{sfx}.prun"));
        let cbe_path = dir.join(format!("cbe{sfx}.prun"));
        // (face turns only: the slice metric would need its own table, not made)
        let with_cbs = with_cbe && !moves.slice;
        let cbs_path = dir.join(format!("cbe-solved{}.prun", crate::cbe::CBS_MIN));
        // Check all tables of the tier before generating any (a table that exists is loaded:
        // it only needs memory).
        crate::preflight::check_plan(&Self::plan(dir, &ae_ix, with_cbe, moves))?;
        let ae = load_or_build(&ae_path, ae_ix.size(), verbose, || {
            build_ae(&ae_ix, moves, verbose)
        })?;
        let cbt = load_or_build(&cb_path, N_CB, verbose, || CbTable::build(moves, verbose))?;
        let cb = CbTable {
            tab: cbt,
            attach: CbTable::attach_targets(),
        };
        let ae_targets = ae_targets(&ae_ix);
        let cbe = if with_cbe {
            let path = cbe_path;
            Some(match crate::cbe::CbeTable::load(&path) {
                Ok(t) => {
                    if verbose {
                        eprintln!(
                            "loaded {} ({:.1} MB)",
                            path.display(),
                            t.codes.len() as f64 * 8e-6
                        );
                    }
                    t
                }
                Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
                    if verbose {
                        eprintln!("generating {} ...", path.display());
                    }
                    let t = crate::cbe::CbeTable::build(moves, verbose);
                    t.save(&path)?;
                    t
                }
                Err(e) => return Err(e),
            })
        } else {
            None
        };
        // The distance-to-solved table is read only when a search can use it (see `LazyCbs`);
        // a missing one is generated now.
        let cbs = match (&cbe, with_cbs) {
            (Some(cbe), true) => Some(std::sync::Arc::new(
                if std::fs::metadata(&cbs_path).is_ok() {
                    if verbose {
                        eprintln!(
                            "found {} (read when optimal mode can use it)",
                            cbs_path.display()
                        );
                    }
                    let lazy = crate::cbe::LazyCbs::new(cbs_path);
                    // FTO_CBS_EAGER=1: read it now, so that every search uses it from its
                    // start (deterministic node counts on one thread, for benchmarks)
                    if std::env::var_os("FTO_CBS_EAGER").is_some_and(|v| v == "1") {
                        let _ = lazy.get(&cbe.ix);
                    }
                    lazy
                } else {
                    if verbose {
                        eprintln!("generating {} ...", cbs_path.display());
                    }
                    let t = crate::cbe::CbsTable::build(&cbe.ix, verbose);
                    t.save(&cbe.ix, &cbs_path)?;
                    if !t.check_walks(&cbe.ix) {
                        return Err(std::io::Error::new(
                            std::io::ErrorKind::InvalidData,
                            format!("{}: generated table fails its checks", cbs_path.display()),
                        ));
                    }
                    crate::cbe::LazyCbs::with(cbs_path, t)
                },
            )),
            _ => None,
        };
        let t = Tables {
            moves,
            ae_ix,
            ae,
            ae_targets,
            cb,
            p2: Phase2::new(),
            cbe,
            cbs,
        };
        t.check_walks(&ae_path, &cb_path)?;
        Ok(t)
    }

    /// A table file that passes its checksum may still have been made for another layout (an
    /// older version's index): its values are then meaningless, and walking down it would not
    /// end. Walk down both mod-3 tables from a few random positions: each walk must reach a
    /// target within the table's largest distance. (A few hundred lookups; the joint table's
    /// saturating codes cannot be walked, so only its size and checksum are checked.)
    ///
    /// # Errors
    ///
    /// If a walk does not reach a target.
    pub fn check_walks(&self, ae_path: &Path, cb_path: &Path) -> std::io::Result<()> {
        use rand::SeedableRng;
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(0x7ab1e);
        for _ in 0..8 {
            let x = fast().of(&Cube::random(&mut rng));
            let ((_, ae_ok), (_, cb_ok)) = self.walks(&x);
            for (ok, path) in [(ae_ok, ae_path), (cb_ok, cb_path)] {
                if !ok {
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!(
                            "{}: not a distance table for this version's layout (made by a \
                             different version?); delete it to regenerate",
                            path.display()
                        ),
                    ));
                }
            }
        }
        Ok(())
    }

    /// Exact AE distance of a position in the metric of the tables (walks down the mod-3
    /// table). If the table is not a distance table (no neighbour is one step closer), the
    /// steps walked so far.
    pub fn ae_exact(&self, x: &Coords) -> u8 {
        self.exact_fast(&fast().of_coords(x)).0
    }

    /// Exact (AE, CB) distances of a fast-coordinate position in the metric of the tables.
    /// If a table is not a distance table (no neighbour is one step closer), the steps walked
    /// so far (a lower bound).
    pub fn exact_fast(&self, x: &Fc) -> (u8, u8) {
        let ((ha, _), (hc, _)) = self.walks(x);
        (ha, hc)
    }

    /// The walks of [`Self::exact_fast`]: (steps, whether a target was reached) for AE and CB.
    /// A walk stops after the table's largest distance, so a table of another layout cannot
    /// make it run forever.
    fn walks(&self, x: &Fc) -> ((u8, bool), (u8, bool)) {
        let f = fast();
        let n_moves = self.moves.n_moves();
        let ae_ix = |y: &Fc| self.ae_ix.index_fast(f, y, f.a_rank(y));
        let cb_ix = |y: &Fc| CbTable::index_fast(f, y);
        // Each step reads the values of all neighbours, with every lookup issued first.
        let walk =
            |ix: &dyn Fn(&Fc) -> usize, tab: &Table2, done: &dyn Fn(&Fc) -> bool| -> (u8, bool) {
                let mut y = *x;
                let mut d = 0;
                let max = tab.hist.len().saturating_sub(1);
                loop {
                    if done(&y) {
                        return (d, true);
                    }
                    if usize::from(d) >= max {
                        return (d, false);
                    }
                    let want = (tab.get(ix(&y)) + 2) % 3;
                    let em = masks(&y.e);
                    let zs: [Fc; MAX_MOVES] = std::array::from_fn(|m| {
                        if m < n_moves {
                            f.apply_m(&y, &em, m)
                        } else {
                            y
                        }
                    });
                    let mut is = [0usize; MAX_MOVES];
                    for m in 0..n_moves {
                        is[m] = ix(&zs[m]);
                        prefetch(&tab.data[is[m] >> 5]);
                    }
                    let Some(m) = (0..n_moves).find(|&m| tab.get(is[m]) == want) else {
                        return (d, false);
                    };
                    y = zs[m];
                    d += 1;
                }
            };
        let ha = walk(&ae_ix, &self.ae, &|y| self.ae_targets.contains(&ae_ix(y)));
        let hc = walk(&cb_ix, &self.cb.tab, &|y| {
            self.cb.attach[y.c as usize] == f.b_rank(y)
        });
        (ha, hc)
    }

    /// Exact CB distance of a position in the metric of the tables. If the table is not a
    /// distance table (no neighbour is one step closer), the steps walked so far.
    pub fn cb_exact(&self, x: &Coords) -> u8 {
        self.exact_fast(&fast().of_coords(x)).1
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{Rng, SeedableRng};

    /// Moves between entries of the unreduced CB coordinate (corner state x B colouring), `N`
    /// moves (16, or 24 with the slices).
    struct RawCbGraph<const N: usize> {
        mt: &'static MoveTables,
        /// B-centre rank after each slice move (slices move no corner)
        slice_b: Vec<[u32; 8]>,
    }

    impl<const N: usize> RawCbGraph<N> {
        fn new() -> Self {
            let d = defs();
            RawCbGraph {
                mt: move_tables(),
                slice_b: label_perm_table(&std::array::from_fn(|j| d.moves[N_MOVES + j].xb)),
            }
        }
    }

    impl<const N: usize> Graph for RawCbGraph<N> {
        const MOVES: usize = N;
        type Node = (usize, usize);
        type Cache = ();
        fn node(&self, i: usize, (): &mut ()) -> (usize, usize) {
            (i / N_MSET, i % N_MSET)
        }
        fn step(&self, &(c, b): &(usize, usize), m: usize) -> usize {
            if m < N_MOVES {
                self.mt.c[c][m] as usize * N_MSET + self.mt.b[b][m] as usize
            } else {
                c * N_MSET + self.slice_b[b][m - N_MOVES] as usize
            }
        }
    }

    /// The unreduced CB table in the metric with `N` moves.
    fn raw_cb<const N: usize>() -> Table2 {
        let attach = CbTable::attach_targets();
        let targets: Vec<usize> = (0..N_CORNER)
            .map(|c| c * N_MSET + attach[c] as usize)
            .collect();
        bfs(
            &RawCbGraph::<N>::new(),
            N_CORNER * N_MSET,
            &targets,
            "cb raw",
            false,
        )
    }

    /// The symmetry-reduced table equals the unreduced one at the positions and all their
    /// neighbours (by every move of the metric).
    fn check_reduced_cb(raw: &Table2, red: &Table2, moves: usize) {
        assert_eq!(red.hist.len(), raw.hist.len());
        let f = fast();
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(31);
        for _ in 0..20_000 {
            let x = f.of(&Cube::random(&mut rng));
            let em = masks(&x.e);
            for m in 0..moves {
                let y = f.apply_m(&x, &em, m);
                let raw_i = y.c as usize * N_MSET + f.b_rank(&y) as usize;
                assert_eq!(raw.get(raw_i), red.get(CbTable::index_fast(f, &y)));
            }
        }
    }

    /// The symmetry-reduced CB table gives the same distances as the unreduced one, whose
    /// depth histogram was computed independently with a Python/GAP sticker model.
    #[test]
    fn reduced_cb_table_matches_the_unreduced_one() {
        let raw = raw_cb::<N_MOVES>();
        assert_eq!(
            raw.hist,
            vec![
                960,
                7_488,
                83_400,
                841_890,
                7_856_382,
                53_951_003,
                178_902_609,
                111_497_121,
                1_675_147
            ]
        );
        check_reduced_cb(&raw, &CbTable::build(MoveSet::FACE, false), N_MOVES);
    }

    /// The same with slices: the reduction by symmetry holds in the slice metric too (the
    /// slices are mapped to slices by every symmetry).
    #[test]
    fn reduced_cb_table_matches_the_unreduced_one_with_slices() {
        let slice = MoveSet {
            wide: false,
            slice: true,
        };
        let raw = raw_cb::<MAX_MOVES>();
        // (as generated by the first slice-metric build, before the table was reduced)
        assert_eq!(
            raw.hist,
            vec![
                960,
                7_488,
                102_480,
                1_302_226,
                14_205_048,
                100_058_880,
                209_949_042,
                29_187_640,
                2_236
            ]
        );
        check_reduced_cb(&raw, &CbTable::build(slice, false), MAX_MOVES);
    }

    /// The generator's neighbour indices equal the index of each moved position.
    fn check_graph<I: AeIndex, const N: usize>(ix: &I) {
        let g = AeGraph::<I, N>::new(ix);
        let f = fast();
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(44);
        for _ in 0..2000 {
            let i = rng.gen_range(0..ix.size());
            let (a, e, ebits) = ix.decode(i);
            let x = f.of_coords(&Coords {
                c: 0,
                a,
                b: 0,
                e,
                ebits,
            });
            let node = g.node(i, &mut None);
            let em = masks(&x.e);
            for m in 0..N {
                let y = f.apply_m(&x, &em, m);
                assert_eq!(
                    g.step(&node, m),
                    ix.index_fast(f, &y, f.a_rank(&y)),
                    "{i} {m}"
                );
            }
        }
        // grouped neighbours of consecutive entries
        let mut cache = None;
        for _ in 0..200 {
            let i0 = rng.gen_range(0..ix.size() - 20);
            for i in i0..i0 + 20 {
                let node = g.node(i, &mut None);
                let mut v = vec![];
                g.neighbours(i, &mut cache, |j| v.push(j));
                let want: Vec<usize> = (0..N).map(|m| g.step(&node, m)).collect();
                assert_eq!(v, want, "{i}");
            }
        }
    }

    /// The medium index has 4,127,692,800 entries, below 2^32, so 32-bit targets (WebAssembly)
    /// can use it. Every index is `class * N_ECOSET + edge rank * 8 + chirality bits`, a sum of
    /// products that are each at most the index itself, so it never overflows a 32-bit
    /// `usize`. Checked here on the host: recomputed in u32 with every step checked, the
    /// indices of random positions and of all their neighbours agree with the 64-bit ones, by
    /// every path the solver and the generator take (`index`, `index_batch`, `index_fast`,
    /// `entry_parts`, `decode`, `equivalents`).
    #[test]
    fn medium_index_fits_32_bits() {
        let ix = AeSym::new(true);
        let classes = u32::try_from(ix.cls.rep.len()).unwrap();
        let stride = u32::try_from(N_ECOSET).unwrap();
        assert_eq!(classes.checked_mul(stride), Some(4_127_692_800));
        assert_eq!(ix.size(), 4_127_692_800);
        // u32 arithmetic, every step checked
        let index32 = |c: u32, e: u32, b: u8| -> usize {
            let i = c
                .checked_mul(stride)
                .and_then(|p| p.checked_add(e.checked_mul(8)?))
                .and_then(|p| p.checked_add(u32::from(b & 7)))
                .filter(|&i| (i as usize) < ix.size())
                .unwrap();
            i as usize
        };
        let f = fast();
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(0x3232);
        let (mut eq, mut top) = (vec![], 0);
        for _ in 0..20_000 {
            let x = f.of(&Cube::random(&mut rng));
            let em = masks(&x.e);
            let ys: Vec<Fc> = (0..N_MOVES).map(|m| f.apply_m(&x, &em, m)).collect();
            let coords: Vec<Coords> = ys.iter().map(|y| fast_to_coords(f, y)).collect();
            let mut batch = vec![0; N_MOVES];
            ix.index_batch(&coords, &mut batch);
            for ((y, cy), &bi) in ys.iter().zip(&coords).zip(&batch) {
                let v = *ix.class_entry(f.a_rank(y));
                let (g, t) = ((v >> 16 & 0xff) as usize, (v >> 24) as usize);
                let (e, b) = f.sym_e(y, &masks(&y.e), g, t);
                let i = index32(v & 0xffff, e, b);
                assert_eq!(ix.index_fast(f, y, f.a_rank(y)), i);
                assert_eq!(ix.index(cy), i);
                assert_eq!(bi, i);
                let (base, gp, flips) = ix.entry_parts(f, y, v).unwrap();
                assert_eq!(base + (f.relab(gp, y.ebits ^ flips) & 7) as usize, i);
                // decode, and the indices of the same position by the class's stabiliser
                let (a, e2, b2) = ix.decode(i);
                assert_eq!(a, ix.cls.rep[(v & 0xffff) as usize]);
                assert_eq!(index32(v & 0xffff, e2, b2), i);
                eq.clear();
                ix.equivalents(i, &mut eq);
                for &j in &eq {
                    let (a3, e3, b3) = ix.decode(j);
                    assert_eq!(a3, a);
                    assert_eq!(index32(v & 0xffff, e3, b3), j);
                }
                top = top.max(i);
            }
        }
        // (the random positions reach the top of the index range)
        assert!(top > ix.size() - ix.size() / 100, "{top}");
    }

    /// The coordinates of a fast-coordinate position (through its cube-independent ranks).
    fn fast_to_coords(f: &FastTables, y: &Fc) -> Coords {
        Coords {
            c: y.c,
            a: f.a_rank(y),
            b: f.b_rank(y),
            e: f.e_rank(y),
            ebits: y.ebits,
        }
    }

    /// The medium index of positions reached by seeded random move sequences (with only
    /// target-independent arithmetic, so the same on every target), condensed into a hash
    /// recorded on x86-64: on a 32-bit target (`scripts/check.sh` runs the library tests on
    /// wasm32) the same value shows that its 32-bit `usize` gives the same indices.
    #[test]
    fn medium_index_fingerprint() {
        let ix = AeSym::new(true);
        let f = fast();
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(0x0e66);
        let mut h = 0xcbf2_9ce4_8422_2325_u64;
        let mut hash = |i: usize| h = (h ^ i as u64).wrapping_mul(0x0100_0000_01b3);
        for _ in 0..4096 {
            let mut cube = Cube::solved();
            for _ in 0..40 {
                cube = cube.apply((rand::RngCore::next_u64(&mut rng) % N_MOVES as u64) as usize);
            }
            let x = f.of(&cube);
            let em = masks(&x.e);
            for m in 0..N_MOVES {
                let y = f.apply_m(&x, &em, m);
                hash(ix.index_fast(f, &y, f.a_rank(&y)));
            }
            hash(ix.index(&Coords::of(&cube)));
        }
        assert_eq!(h, 0x869f_8b6a_6288_3942);
    }

    #[test]
    fn ae_graph_steps_match_index() {
        check_graph::<_, N_MOVES>(&AeSym::new(false));
        check_graph::<_, MAX_MOVES>(&AeSym::new(false));
        // (the medium index, 4.13e9 entries, fits 32-bit targets too; the huge one, 9.1e10
        // entries, needs 64-bit sizes)
        check_graph::<_, N_MOVES>(&AeSym::new(true));
        check_graph::<_, MAX_MOVES>(&AeSym::new(true));
        #[cfg(target_pointer_width = "64")]
        {
            check_graph::<_, N_MOVES>(&AeRelabel::new());
            check_graph::<_, MAX_MOVES>(&AeRelabel::new());
        }
    }
}
