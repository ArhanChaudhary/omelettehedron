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
use std::collections::HashMap;
use std::path::{Path, PathBuf};

#[derive(
    Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum, serde::Serialize, serde::Deserialize,
)]
pub enum Tier {
    /// ~260 MB: edge triples without chirality, 288-fold symmetry
    Small,
    /// ~1.2 GB: full edge coset, 288-fold symmetry
    Medium,
    /// medium + the joint corner/B-centre/edge-triple table (~8 GB)
    Large,
    /// relabel-only edge-coset table (~23 GB) + the joint table (~31 GB in total)
    Huge,
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
    fn index_fast(&self, f: &FastTables, x: &Fc, a_rank: u32) -> usize;
    /// The symmetry-class table read first by `index_fast` (for prefetching).
    fn class_entry(&self, a_rank: u32) -> &u32;
    /// Pipelined `index_fast` over a batch.
    #[inline(always)]
    fn index_batch_fast(&self, f: &FastTables, xs: &[Fc], out: &mut [usize]) {
        let n = xs.len();
        let mut ar = [0u32; N_MOVES];
        for i in 0..n {
            ar[i] = f.a_rank(&xs[i]);
            prefetch(self.class_entry(ar[i]));
        }
        for i in 0..n {
            out[i] = self.index_fast(f, &xs[i], ar[i]);
        }
    }
    /// `index_batch_fast` that prefetches each table word as soon as its index is known, so
    /// the first lookups' latency overlaps the index computation of the later positions.
    #[inline(always)]
    fn index_batch_fast_pf(&self, f: &FastTables, xs: &[Fc], out: &mut [usize], tab: &[u64]) {
        let n = xs.len();
        let mut ar = [0u32; N_MOVES];
        for i in 0..n {
            ar[i] = f.a_rank(&xs[i]);
            prefetch(self.class_entry(ar[i]));
        }
        for i in 0..n {
            out[i] = self.index_fast(f, &xs[i], ar[i]);
            prefetch(&tab[out[i] >> 5]);
        }
    }
}

/// Hint that `p` will be read soon.
#[inline(always)]
pub fn prefetch<T>(p: &T) {
    #[cfg(target_arch = "x86_64")]
    // SAFETY: prefetch is a hint and never faults
    unsafe {
        core::arch::x86_64::_mm_prefetch(
            std::ptr::from_ref(p).cast::<i8>(),
            core::arch::x86_64::_MM_HINT_T0,
        );
    }
}

pub struct AeSym {
    /// `cls.sym` holds class id | g << 16 | t << 24 (the symmetry pair itself, not its index)
    pub cls: AClasses,
    pub chiral: bool,
    /// symmetry pairs (other than the identity) fixing each class representative
    pub stab: Vec<Vec<(u8, u8)>>,
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
        AeSym { cls, chiral, stab }
    }
}

impl AeIndex for AeSym {
    fn size(&self) -> usize {
        self.cls.rep.len() * if self.chiral { N_ECOSET } else { N_MSET }
    }
    #[inline(always)]
    fn index(&self, x: &Coords) -> usize {
        let s = syms();
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
        let s = syms();
        let n = xs.len();
        let mut v = [0u32; N_MOVES];
        let mut e1 = [0u32; N_MOVES];
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
        let mt = move_tables();
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
    fn index_fast(&self, f: &FastTables, x: &Fc, a_rank: u32) -> usize {
        let v = self.cls.sym[a_rank as usize];
        let (g, t) = (v >> 16 & 0xff, v >> 24);
        let (e, b) = f.sym_e(x, &masks(&x.e), g as usize, t as usize);
        let c = (v & 0xffff) as usize;
        if self.chiral {
            c * N_ECOSET + e as usize * 8 + (b & 7) as usize
        } else {
            c * N_MSET + e as usize
        }
    }
    #[inline(always)]
    fn class_entry(&self, a_rank: u32) -> &u32 {
        &self.cls.sym[a_rank as usize]
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
        let s = syms();
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
}

impl AeRelabel {
    pub fn new() -> AeRelabel {
        let pairs: Vec<(u8, u8)> = (0..N_GEO as u8)
            .filter(|&g| crate::geometry::det(&syms().mats[g as usize]) == 1)
            .map(|g| (g, 0))
            .collect();
        AeRelabel {
            cls: AClasses::new(&pairs),
        }
    }
}

impl AeIndex for AeRelabel {
    fn size(&self) -> usize {
        self.cls.rep.len() * N_ECOSET
    }
    #[inline(always)]
    fn index(&self, x: &Coords) -> usize {
        let s = syms();
        let v = self.cls.sym[x.a as usize];
        // pair index == rotation index in the identity-first list of rotations
        let g = ROT_IDX[(v >> 16) as usize] as usize;
        let e = s.relab_e[g * N_MSET + x.e as usize];
        let b = s.relab_bits[g][x.ebits as usize];
        (v & 0xffff) as usize * N_ECOSET + e as usize * 8 + (b & 7) as usize
    }
    #[inline(always)]
    fn index_batch(&self, xs: &[Coords], out: &mut [usize]) {
        let s = syms();
        let n = xs.len();
        let mut v = [0u32; N_MOVES];
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
        let mt = move_tables();
        let c = i / N_ECOSET;
        let r = i % N_ECOSET;
        let e = (r / 8) as u32;
        let b3 = (r % 8) as u8;
        let bit3 = mt.eref[e as usize] ^ (b3 ^ b3 >> 1 ^ b3 >> 2) & 1;
        (self.cls.rep[c], e, b3 | bit3 << 3)
    }
    #[inline(always)]
    fn index_fast(&self, f: &FastTables, x: &Fc, a_rank: u32) -> usize {
        let v = self.cls.sym[a_rank as usize];
        let g = ROT_IDX[(v >> 16) as usize] as usize;
        let (e, b) = f.relabel_e(x, g);
        (v & 0xffff) as usize * N_ECOSET + e as usize * 8 + (b & 7) as usize
    }
    #[inline(always)]
    fn class_entry(&self, a_rank: u32) -> &u32 {
        &self.cls.sym[a_rank as usize]
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
    pub dist: Vec<u8>,
    pub emoves: Vec<[u16; N_MOVES]>,
}

pub const N_P2_EDGE: usize = 972;

impl Phase2 {
    pub fn new() -> Phase2 {
        let s0 = Cube::solved();
        let amoves: Vec<usize> = (0..N_MOVES).filter(|&m| is_a_move(m)).collect();
        let mut edge_id: HashMap<[u8; 12], u16> = HashMap::new();
        let mut elist = vec![s0.ep];
        edge_id.insert(s0.ep, 0);
        let mut i = 0;
        while i < elist.len() {
            let c = Cube { ep: elist[i], ..s0 };
            i += 1;
            for &m in &amoves {
                let n = c.apply(m).ep;
                if let std::collections::hash_map::Entry::Vacant(e) = edge_id.entry(n) {
                    e.insert(elist.len() as u16);
                    elist.push(n);
                }
            }
        }
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
            dist,
            emoves,
        }
    }

    /// Phase-2 length of a position known to be in `H_A` up to its edge permutation
    /// (None if the edges are not an `H_A` arrangement).
    #[inline]
    pub fn length(&self, cube: &Cube) -> Option<u8> {
        let e = *self.edge_id.get(&cube.ep)?;
        Some(self.dist[corners().id(cube) as usize * N_P2_EDGE + e as usize])
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
    /// experiment (FTO_CBSYM=1): `tab` is the symmetry-reduced table with this index
    pub sym: Option<crate::cbsym::CbSymIndex>,
}

pub const N_CB: usize = N_CORNER * N_MSET;

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

    pub fn build(verbose: bool) -> Table2 {
        let mt = move_tables();
        let attach = Self::attach_targets();
        let targets: Vec<usize> = (0..N_CORNER)
            .map(|c| c * N_MSET + attach[c] as usize)
            .collect();
        bfs(
            N_CB,
            &targets,
            |i, out| {
                let (c, b) = (i / N_MSET, i % N_MSET);
                let mc = &mt.c[c];
                let mb = &mt.b[b];
                for m in 0..N_MOVES {
                    out.push(mc[m] as usize * N_MSET + mb[m] as usize);
                }
            },
            "cb",
            verbose,
        )
    }

    #[inline(always)]
    pub fn index(x: &Coords) -> usize {
        x.c as usize * N_MSET + x.b as usize
    }

    #[inline(always)]
    pub fn index_fast(f: &FastTables, x: &Fc) -> usize {
        x.c as usize * N_MSET + f.b_rank(x) as usize
    }

    /// Index into `tab` (plain or symmetry-reduced).
    #[inline(always)]
    pub fn ix(&self, f: &FastTables, x: &Fc) -> usize {
        match &self.sym {
            Some(s) => s.index(f, x),
            None => Self::index_fast(f, x),
        }
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

pub fn build_ae<I: AeIndex>(ix: &I, verbose: bool) -> Table2 {
    let mt = move_tables();
    bfs(
        ix.size(),
        &ae_targets(ix),
        |i, out| {
            let (a, e, bits) = ix.decode(i);
            let x = Coords {
                c: 0,
                a,
                b: 0,
                e,
                ebits: bits,
            };
            for m in 0..N_MOVES {
                let j = ix.index(&mt.apply(&x, m));
                out.push(j);
                ix.equivalents(j, out);
            }
        },
        ix.name(),
        verbose,
    )
}

/// All tables needed by the solver, for one AE index type.
pub struct Tables<I: AeIndex> {
    pub ae_ix: I,
    pub ae: Table2,
    pub ae_targets: Vec<usize>,
    pub cb: CbTable,
    pub p2: Phase2,
    pub cbe: Option<crate::cbe::CbeTable>,
    /// experiment (FTO_CBA=1): joint corner x B-centre x A-colour table
    pub cba: Option<crate::cba::CbaTable>,
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
            if verbose {
                eprintln!("generating {} ...", path.display());
            }
            let t = build();
            std::fs::create_dir_all(path.parent().unwrap_or(Path::new(".")))?;
            t.save(path)?;
            Ok(t)
        }
        Err(e) => Err(e),
    }
}

impl<I: AeIndex> Tables<I> {
    /// Build all tables in memory without touching the disk.
    pub fn generate(ae_ix: I, verbose: bool) -> Tables<I> {
        let ae = build_ae(&ae_ix, verbose);
        let cb = CbTable {
            tab: CbTable::build(verbose),
            attach: CbTable::attach_targets(),
            sym: None,
        };
        let ae_targets = ae_targets(&ae_ix);
        Tables {
            ae_ix,
            ae,
            ae_targets,
            cb,
            p2: Phase2::new(),
            cbe: None,
            cba: None,
        }
    }

    /// [`Self::load_with`] without the joint table.
    ///
    /// # Errors
    ///
    /// As [`Self::load_with`].
    pub fn load(dir: &Path, ae_ix: I, verbose: bool) -> std::io::Result<Tables<I>> {
        Self::load_with(dir, ae_ix, false, verbose)
    }

    /// Load (generating if needed) the tables, optionally with the joint CBE table.
    ///
    /// # Errors
    ///
    /// If a table file cannot be read or written, or holds a table of the wrong format, size
    /// or checksum.
    pub fn load_with(
        dir: &Path,
        ae_ix: I,
        with_cbe: bool,
        verbose: bool,
    ) -> std::io::Result<Tables<I>> {
        let ae_path = dir.join(format!("{}.prun", ae_ix.name()));
        let ae = load_or_build(&ae_path, ae_ix.size(), verbose, || {
            build_ae(&ae_ix, verbose)
        })?;
        let cbt = load_or_build(&dir.join("cb.prun"), N_CB, verbose, || {
            CbTable::build(verbose)
        })?;
        let mut cb = CbTable {
            tab: cbt,
            attach: CbTable::attach_targets(),
            sym: None,
        };
        if std::env::var("FTO_CBSYM").is_ok_and(|v| v == "1") {
            let t0 = std::time::Instant::now();
            let ix = crate::cbsym::CbSymIndex::new();
            let tab = ix.build(verbose);
            // check against the plain table on random positions and their neighbours
            let f = fast();
            let mut rng = <rand_xoshiro::Xoshiro256PlusPlus as rand::SeedableRng>::seed_from_u64(5);
            for _ in 0..2000 {
                let x = f.of(&Cube::random(&mut rng));
                let em = masks(&x.e);
                for m in 0..N_MOVES {
                    let y = f.apply(&x, &em, m);
                    assert_eq!(
                        tab.get(ix.index(f, &y)),
                        cb.tab.get(CbTable::index_fast(f, &y)),
                        "cb-sym mismatch"
                    );
                }
            }
            eprintln!(
                "cb-sym: {} classes, {:.1} MB, hist {:?}, {:.2}s",
                ix.rep.len(),
                tab.bytes() as f64 / 1e6,
                tab.hist,
                t0.elapsed().as_secs_f64()
            );
            cb.tab = tab;
            cb.sym = Some(ix);
        }
        let ae_targets = ae_targets(&ae_ix);
        let cbe = if with_cbe {
            let path = dir.join("cbe.prun");
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
                    let t = crate::cbe::CbeTable::build(verbose);
                    t.save(&path)?;
                    t
                }
                Err(e) => return Err(e),
            })
        } else {
            None
        };
        Ok(Tables {
            ae_ix,
            ae,
            ae_targets,
            cb,
            p2: Phase2::new(),
            cbe,
            cba: if std::env::var("FTO_CBA").is_ok_and(|v| v == "1") {
                Some(crate::cba::CbaTable::build(verbose))
            } else {
                None
            },
        })
    }

    /// Exact AE distance of a position (walks down the mod-3 table). If the table is not a
    /// distance table (no neighbour is one step closer), the steps walked so far.
    pub fn ae_exact(&self, x: &Coords) -> u8 {
        let mt = move_tables();
        let mut x = *x;
        let mut d = 0;
        loop {
            let i = self.ae_ix.index(&x);
            if self.ae_targets.contains(&i) {
                return d;
            }
            let v = self.ae.get(i);
            let want = (v + 2) % 3;
            let next = (0..N_MOVES)
                .map(|m| mt.apply(&x, m))
                .find(|y| self.ae.get(self.ae_ix.index(y)) == want);
            let Some(next) = next else {
                return d;
            };
            x = next;
            d += 1;
        }
    }

    /// Exact (AE, CB) distances of a fast-coordinate position. If a table is not a distance
    /// table (no neighbour is one step closer), the steps walked so far (a lower bound).
    pub fn exact_fast(&self, x: &Fc) -> (u8, u8) {
        let f = fast();
        let ae_ix = |y: &Fc| self.ae_ix.index_fast(f, y, f.a_rank(y));
        let cb_ix = |y: &Fc| self.cb.ix(f, y);
        // Each step reads the values of all neighbours, with every lookup issued first.
        let walk = |ix: &dyn Fn(&Fc) -> usize, tab: &Table2, done: &dyn Fn(&Fc) -> bool| -> u8 {
            let mut y = *x;
            let mut d = 0;
            while !done(&y) {
                let want = (tab.get(ix(&y)) + 2) % 3;
                let em = masks(&y.e);
                let zs: [Fc; N_MOVES] = std::array::from_fn(|m| f.apply(&y, &em, m));
                let mut is = [0usize; N_MOVES];
                for m in 0..N_MOVES {
                    is[m] = ix(&zs[m]);
                    prefetch(&tab.data[is[m] >> 5]);
                }
                let Some(m) = (0..N_MOVES).find(|&m| tab.get(is[m]) == want) else {
                    break;
                };
                y = zs[m];
                d += 1;
            }
            d
        };
        let ha = walk(&ae_ix, &self.ae, &|y| self.ae_targets.contains(&ae_ix(y)));
        let hc = walk(&cb_ix, &self.cb.tab, &|y| {
            self.cb.attach[y.c as usize] == f.b_rank(y)
        });
        (ha, hc)
    }

    /// Exact CB distance of a position. If the table is not a distance table (no neighbour is
    /// one step closer), the steps walked so far.
    pub fn cb_exact(&self, x: &Coords) -> u8 {
        let mt = move_tables();
        let mut x = *x;
        let mut d = 0;
        loop {
            if self.cb.attach[x.c as usize] == x.b {
                return d;
            }
            let v = self.cb.tab.get(CbTable::index(&x));
            let want = (v + 2) % 3;
            let next = (0..N_MOVES)
                .map(|m| mt.apply(&x, m))
                .find(|y| self.cb.tab.get(CbTable::index(y)) == want);
            let Some(next) = next else {
                return d;
            };
            x = next;
            d += 1;
        }
    }
}
