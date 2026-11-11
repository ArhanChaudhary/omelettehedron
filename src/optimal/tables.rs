//! Pattern databases of the optimal solver: distances to the solved state.
//!
//! Every table is the exact distance (mod 3, two bits per entry, see [`crate::prune`]) of a
//! projection of the position to the projection of the solved state, computed by breadth-first
//! search from that single state. A projection's distance is a lower bound on the distance of
//! the position itself (every solution projects to a path of the same length), and it changes
//! by at most one per move, so the search carries exact values down the tree and decodes each
//! child's value from its parent's. The bound of a position is the maximum over the tables.
//!
//! Tables of a tier are listed by [`Bounds`]: a tier is `K` lookups, each an index function
//! of the position and a table. New tiers (stronger tables) add an implementation.
//!
//! Besides its own tables (distances to solved), every tier also reads the two-phase solver's
//! A-centre/edge table ([`AeBound`]), through that solver's public index, on the position and
//! on its `SWAP_AB` conjugate: it holds distances to `H_A`, the positions solvable with the
//! class-A faces, and solved is in `H_A` (and in `H_B`, which the conjugate's distance to
//! `H_A` measures), so these are lower bounds on the distance to solved as well. They are not
//! distances to a single state, but each is still an exact BFS distance (mod 3) that changes by
//! at most one per move, so they are decoded like the others.

use super::coords::*;
use crate::cbe::CbeTable;
use crate::cube::*;
use crate::fast::{fast, FastTables};
use crate::preflight::{check_plan, Planned};
use crate::prune::{bfs, Graph, Table2};
use crate::tables::{ae_targets, build_ae, prefetch_at, AeIndex, AeSym};
use std::path::{Path, PathBuf};

/// Most lookups per position of any tier (the search's scratch arrays have this many slots).
pub const MAX_K: usize = 8;

/// Table size of the optimal solver (`omelettehedron optimal ... --tier`).
#[derive(
    Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum, serde::Serialize, serde::Deserialize,
)]
pub enum OptTier {
    /// ~195 MB: corners x centres (one 4.7 MB table for both centre orbits), the edge
    /// permutation (60 MB) and the two-phase A-centre/edge-triple table (129 MB)
    Small,
    /// ~3.4 GB: corners x centres, corners x edges (48 symmetries, 2.4 GB) and the two-phase
    /// A-centre/edge table with chirality (1.03 GB); 64-bit targets only
    Medium,
    /// ~10.4 GB: medium and the two-phase corner/B-centre/edge-triple table (7.0 GB) on both
    /// sides; 64-bit targets only
    Large,
}

impl OptTier {
    pub fn name(self) -> &'static str {
        match self {
            OptTier::Small => "small",
            OptTier::Medium => "medium",
            OptTier::Large => "large",
        }
    }
}

/// The lower bounds of a tier: `K` table lookups per position.
pub trait Bounds: Sync {
    /// Lookups per position (at most [`MAX_K`]).
    const K: usize;
    /// Index of lookup `k` (`k < K`) of position `p`.
    fn index(&self, k: usize, p: &Pos) -> usize;
    /// The table of lookup `k`.
    fn table(&self, k: usize) -> &Table2;
    /// Whether index `i` of lookup `k` is the solved state's (distance 0).
    fn is_target(&self, k: usize, i: usize) -> bool;
    /// Short name of lookup `k` (for messages).
    fn name(&self, k: usize) -> &'static str;
    /// Saturating lower bounds (codes, not exact distances, so read without the parent's
    /// values): `S` groups of 4 lookups (at most [`MAX_S`]), each giving a bound of at most
    /// `SAT_MAX`. A child is checked against them only when they can prune it: when fewer than
    /// `SAT_MAX` moves are left after it.
    const S: usize = 0;
    const SAT_MAX: u8 = 0;
    /// Compute the indices of saturating group `s` of `p` into `out` and prefetch them.
    #[inline(always)]
    fn sat_index(&self, _s: usize, _p: &Pos, _out: &mut [usize; 4]) {}
    /// The lower bound of saturating group `s` from its indices.
    #[inline(always)]
    fn sat_bound(&self, _s: usize, _idx: &[usize; 4]) -> u8 {
        0
    }
}

/// Most saturating groups of any tier.
pub const MAX_S: usize = 2;

/// Moves between entries of the corner x B-centre table.
struct CcGraph<'a> {
    ix: &'a CcIndex,
    mt: &'static OptMoves,
}

impl Graph for CcGraph<'_> {
    const MOVES: usize = N_MOVES;
    type Node = (u16, [u8; 16]);
    type Cache = ();
    fn node(&self, i: usize, (): &mut ()) -> (u16, [u8; 16]) {
        self.ix.decode(i)
    }
    #[inline(always)]
    fn step(&self, &(c, x): &(u16, [u8; 16]), m: usize) -> usize {
        let (c2, x2) = self.mt.apply_cb(c, &x, m);
        self.ix.index_b(c2, &x2)
    }
    fn equivalents(&self, j: usize, out: &mut Vec<usize>) {
        self.ix.equivalents(j, out);
    }
}

/// Moves between edge permutations.
struct EdgeGraph {
    mt: &'static OptMoves,
}

impl Graph for EdgeGraph {
    const MOVES: usize = N_MOVES;
    type Node = [u8; 16];
    type Cache = ();
    fn node(&self, i: usize, (): &mut ()) -> [u8; 16] {
        eperm_unrank(i as u32)
    }
    #[inline(always)]
    fn step(&self, x: &[u8; 16], m: usize) -> usize {
        eperm_rank(&self.mt.apply_ep(x, m)) as usize
    }
}

/// The corner x centre table (distance to solved of the corners and the B-centres).
pub fn build_cc(ix: &CcIndex, verbose: bool) -> Table2 {
    let g = CcGraph {
        ix,
        mt: opt_moves(),
    };
    bfs(&g, cc_size(ix), &ix.targets(), "opt-cc", verbose)
}

/// The edge permutation table.
pub fn build_edges(verbose: bool) -> Table2 {
    bfs(
        &EdgeGraph { mt: opt_moves() },
        N_EPERM,
        &[0],
        "opt-edges",
        verbose,
    )
}

/// Moves between entries of the corner x edge table.
struct CeGraph<'a> {
    ix: &'a CeIndex,
    mt: &'static OptMoves,
}

impl Graph for CeGraph<'_> {
    const MOVES: usize = N_MOVES;
    type Node = (u16, [u8; 16]);
    type Cache = ();
    fn node(&self, i: usize, (): &mut ()) -> (u16, [u8; 16]) {
        self.ix.decode(i)
    }
    #[inline(always)]
    fn step(&self, &(c, ep): &(u16, [u8; 16]), m: usize) -> usize {
        self.ix
            .index(self.mt.corner(c, m), &self.mt.apply_ep(&ep, m))
    }
    fn equivalents(&self, j: usize, out: &mut Vec<usize>) {
        self.ix.equivalents(j, out);
    }
}

/// The corner x edge table.
pub fn build_ce(ix: &CeIndex, verbose: bool) -> Table2 {
    let g = CeGraph {
        ix,
        mt: opt_moves(),
    };
    bfs(&g, ix.size(), &[ix.target()], "opt-ce", verbose)
}

/// The two-phase solver's A-centre/edge table (distance to `H_A`) and its index, read-only:
/// the triples-only table for the small tier, the one with edge chirality for the medium tier
/// (the two-phase small and medium tiers' files, shared with that solver).
pub struct AeBound {
    pub ix: AeSym,
    pub tab: Table2,
    /// the indices of `H_A` (the solved class and its equivalents)
    targets: Vec<usize>,
    /// the two-phase coordinates of a position and of its conjugate
    maps: FcMaps,
    f: &'static FastTables,
}

impl AeBound {
    /// The index for a tier: with edge chirality from the medium tier on.
    fn index_for(tier: OptTier) -> AeSym {
        AeSym::new(tier != OptTier::Small)
    }

    fn new(ix: AeSym, tab: Table2) -> AeBound {
        AeBound {
            targets: ae_targets(&ix),
            ix,
            tab,
            maps: FcMaps::new(),
            f: fast(),
        }
    }

    /// The two-phase coordinates of the position, or (`swap`) of its `SWAP_AB` conjugate.
    #[inline(always)]
    fn fc(&self, p: &Pos, swap: bool) -> crate::fast::Fc {
        if swap {
            self.maps.fc_swap(p)
        } else {
            self.maps.fc(p)
        }
    }

    /// Index of the position (`swap == false`: distance to `H_A`) or of its `SWAP_AB`
    /// conjugate (`swap == true`: the position's distance to `H_B`).
    #[inline(always)]
    pub fn index(&self, p: &Pos, swap: bool) -> usize {
        let x = self.fc(p, swap);
        self.ix.index_fast(self.f, &x, self.f.a_rank(&x))
    }

    fn is_target(&self, i: usize) -> bool {
        self.targets.contains(&i)
    }
}

/// All tables of a tier.
pub struct OptTables {
    pub tier: OptTier,
    pub cc_ix: CcIndex,
    /// corners x B-centres (and, through `SWAP_AB`, corners x A-centres)
    pub cc: Table2,
    cc_target: usize,
    /// the edge permutation (small tier)
    pub edges: Option<Table2>,
    /// corners x edge permutation (medium tier)
    pub ce: Option<(CeIndex, Table2)>,
    /// the two-phase A-centre/edge table (every tier)
    pub ae: AeBound,
    /// the two-phase corner/B-centre/edge-triple table (large tier)
    pub cbe: Option<CbeTable>,
}

/// File names of the tables (all start with `opt-`, apart from the two-phase solver's).
const CC_FILE: &str = "opt-cc-sym24.prun";
const EDGES_FILE: &str = "opt-edges.prun";
const CE_FILE: &str = "opt-ce-sym48.prun";
/// The two-phase solver's joint table (its large tier's file, shared).
const CBE_FILE: &str = "cbe.prun";

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

/// The tier needs 64-bit indices (more than 2^32 entries).
fn too_large(tier: OptTier) -> bool {
    cfg!(target_pointer_width = "32") && tier != OptTier::Small
}

impl OptTables {
    /// Build the tables of a tier in memory, without touching the disk.
    ///
    /// # Errors
    ///
    /// If the tier needs 64-bit indices on a 32-bit target.
    pub fn generate(tier: OptTier, verbose: bool) -> std::io::Result<OptTables> {
        if too_large(tier) {
            return Err(std::io::Error::other(
                "only --tier small is supported on 32-bit targets",
            ));
        }
        let cc_ix = CcIndex::new();
        let cc = build_cc(&cc_ix, verbose);
        let (edges, ce) = match tier {
            OptTier::Small => (Some(build_edges(verbose)), None),
            OptTier::Medium | OptTier::Large => {
                let ix = CeIndex::new();
                let t = build_ce(&ix, verbose);
                (None, Some((ix, t)))
            }
        };
        let ae_ix = AeBound::index_for(tier);
        let ae_tab = build_ae(&ae_ix, MoveSet::FACE, verbose);
        let ae = AeBound::new(ae_ix, ae_tab);
        let cbe = (tier == OptTier::Large).then(|| CbeTable::build(MoveSet::FACE, verbose));
        Ok(Self::assemble(tier, cc_ix, cc, edges, ce, ae, cbe))
    }

    fn assemble(
        tier: OptTier,
        cc_ix: CcIndex,
        cc: Table2,
        edges: Option<Table2>,
        ce: Option<(CeIndex, Table2)>,
        ae: AeBound,
        cbe: Option<CbeTable>,
    ) -> OptTables {
        let cc_target = cc_ix.targets().first().copied().unwrap_or(0);
        OptTables {
            tier,
            cc_ix,
            cc,
            cc_target,
            edges,
            ce,
            ae,
            cbe,
        }
    }

    /// Load the tables of a tier from `dir`, generating (and saving) the missing ones. All
    /// tables are checked against the free disk space and memory before any is generated.
    ///
    /// # Errors
    ///
    /// If a table file cannot be read or written, holds a table of the wrong format, size,
    /// checksum or layout, or the tables do not fit (see [`crate::preflight`]), or the tier
    /// needs 64-bit indices on a 32-bit target.
    pub fn load(dir: &Path, tier: OptTier, verbose: bool) -> std::io::Result<OptTables> {
        if too_large(tier) {
            return Err(std::io::Error::other(
                "only --tier small is supported on 32-bit targets",
            ));
        }
        let cc_ix = CcIndex::new();
        let cc_path = dir.join(CC_FILE);
        let edges_path = dir.join(EDGES_FILE);
        let ce_path = dir.join(CE_FILE);
        let ce_ix = (tier != OptTier::Small).then(CeIndex::new);
        let cbe_path = dir.join(CBE_FILE);
        let ae_ix = AeBound::index_for(tier);
        // (the two-phase solver's file of the same table: generated by either solver, shared)
        let ae_path = dir.join(format!("{}.prun", ae_ix.name()));
        let planned = |path: &PathBuf, entries: usize| {
            let bytes =
                std::fs::metadata(path).map_or_else(|_| entries.div_ceil(4) as u64, |m| m.len());
            Planned {
                path: path.clone(),
                bytes,
                working_set: 0,
            }
        };
        let mut plan = vec![
            planned(&cc_path, cc_size(&cc_ix)),
            planned(&ae_path, ae_ix.size()),
        ];
        match &ce_ix {
            None => plan.push(planned(&edges_path, N_EPERM)),
            Some(ix) => plan.push(planned(&ce_path, ix.size())),
        }
        if tier == OptTier::Large {
            // (the index is built only for its size when the file is missing)
            let bytes = std::fs::metadata(&cbe_path).map_or_else(
                |_| crate::cbe::CbeIndex::new().size().div_ceil(4) as u64,
                |m| m.len(),
            );
            plan.push(Planned {
                path: cbe_path.clone(),
                bytes,
                working_set: 0,
            });
        }
        check_plan(&plan)?;
        let cc = load_or_build(&cc_path, cc_size(&cc_ix), verbose, || {
            build_cc(&cc_ix, verbose)
        })?;
        let (edges, ce, second) = match ce_ix {
            None => (
                Some(load_or_build(&edges_path, N_EPERM, verbose, || {
                    build_edges(verbose)
                })?),
                None,
                edges_path,
            ),
            Some(ix) => {
                let t = load_or_build(&ce_path, ix.size(), verbose, || build_ce(&ix, verbose))?;
                (None, Some((ix, t)), ce_path)
            }
        };
        let ae_tab = load_or_build(&ae_path, ae_ix.size(), verbose, || {
            build_ae(&ae_ix, MoveSet::FACE, verbose)
        })?;
        let ae = AeBound::new(ae_ix, ae_tab);
        let cbe = if tier == OptTier::Large {
            Some(match CbeTable::load(&cbe_path) {
                Ok(t) => {
                    if verbose {
                        eprintln!("loaded {}", cbe_path.display());
                    }
                    t
                }
                Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
                    if verbose {
                        eprintln!("generating {} ...", cbe_path.display());
                    }
                    let t = CbeTable::build(MoveSet::FACE, verbose);
                    t.save(&cbe_path)?;
                    t
                }
                Err(e) => return Err(e),
            })
        } else {
            None
        };
        let t = Self::assemble(tier, cc_ix, cc, edges, ce, ae, cbe);
        t.check_walks(&[&second, &cc_path, &ae_path])?;
        Ok(t)
    }

    /// The lower bounds of the small tier (`None` if the tables are another tier's).
    pub fn small(&self) -> Option<SmallBounds<'_>> {
        Some(SmallBounds {
            cc_ix: &self.cc_ix,
            cc: &self.cc,
            cc_target: self.cc_target,
            edges: self.edges.as_ref()?,
            ae: &self.ae,
        })
    }

    /// The lower bounds of the large tier (`None` if the tables are another tier's).
    pub fn large(&self) -> Option<LargeBounds<'_>> {
        Some(LargeBounds {
            m: self.medium()?,
            cbe: self.cbe.as_ref()?,
        })
    }

    /// The lower bounds of the medium tier (`None` if the tables are another tier's).
    pub fn medium(&self) -> Option<MediumBounds<'_>> {
        let (ce_ix, ce) = self.ce.as_ref()?;
        Some(MediumBounds {
            cc_ix: &self.cc_ix,
            cc: &self.cc,
            cc_target: self.cc_target,
            ce_ix,
            ce,
            ce_target: ce_ix.target(),
            ae: &self.ae,
        })
    }

    /// A table file that passes its checksum may still have been made for another layout: walk
    /// down every table from a few random positions; each walk must reach the solved entry
    /// within the table's largest distance.
    ///
    /// # Errors
    ///
    /// If a walk does not reach the solved entry (`paths`: the tier's own file, the corner x
    /// centre file and the A-centre/edge file, for the message).
    fn check_walks(&self, paths: &[&Path]) -> std::io::Result<()> {
        use rand::SeedableRng;
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(0x0071_ab1e);
        let positions: Vec<Pos> = (0..8).map(|_| Pos::of(&Cube::random(&mut rng))).collect();
        // (lookup 0 is the tier's own table, 1 and 2 share the corner x centre table, 3 and 4
        // the A-centre/edge table)
        let bad = |k: usize| -> std::io::Error {
            let file = paths[match k {
                0 => 0,
                1 | 2 => 1,
                _ => 2,
            }
            .min(paths.len().saturating_sub(1))];
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "{}: not a distance table for this version's layout (made by a different \
                     version?); delete it to regenerate",
                    file.display()
                ),
            )
        };
        // (the joint table's saturating codes cannot be walked: only its size and checksum are
        // checked, as in the two-phase solver)
        let r = if let Some(b) = self.small() {
            check_walks_of(&b, &positions)
        } else if let Some(b) = self.medium() {
            check_walks_of(&b, &positions)
        } else {
            Err(0)
        };
        r.map_err(bad)
    }

    /// Bytes of the tables in memory.
    pub fn bytes(&self) -> usize {
        self.cc.bytes()
            + self.edges.as_ref().map_or(0, Table2::bytes)
            + self.ce.as_ref().map_or(0, |(_, t)| t.bytes())
            + self.ae.tab.bytes()
            + self.cbe.as_ref().map_or(0, |t| t.codes.len() * 8)
    }
}

/// `Err(k)` if a walk down lookup `k` from one of `ps` does not reach the solved entry.
fn check_walks_of<B: Bounds>(b: &B, ps: &[Pos]) -> Result<(), usize> {
    for p in ps {
        for k in 0..B::K {
            if walk(b, k, p).is_none() {
                return Err(k);
            }
        }
    }
    Ok(())
}

/// The small tier: edges, corners x B-centres, corners x A-centres, and the A-centre/edge
/// table on the position and on its `SWAP_AB` conjugate.
pub struct SmallBounds<'a> {
    cc_ix: &'a CcIndex,
    cc: &'a Table2,
    cc_target: usize,
    edges: &'a Table2,
    ae: &'a AeBound,
}

impl Bounds for SmallBounds<'_> {
    const K: usize = 5;
    #[inline(always)]
    fn index(&self, k: usize, p: &Pos) -> usize {
        match k {
            0 => eperm_rank(&p.ep) as usize,
            1 => self.cc_ix.index_b(p.c, &p.xb),
            2 => self.cc_ix.index_a(p.c, &p.xa),
            3 => self.ae.index(p, false),
            _ => self.ae.index(p, true),
        }
    }
    #[inline(always)]
    fn table(&self, k: usize) -> &Table2 {
        match k {
            0 => self.edges,
            1 | 2 => self.cc,
            _ => &self.ae.tab,
        }
    }
    fn is_target(&self, k: usize, i: usize) -> bool {
        match k {
            0 => i == 0,
            1 | 2 => i == self.cc_target,
            _ => self.ae.is_target(i),
        }
    }
    fn name(&self, k: usize) -> &'static str {
        [
            "edges",
            "corners+B",
            "corners+A",
            "AE (H_A)",
            "AE of SWAP_AB (H_B)",
        ][k.min(4)]
    }
}

/// The medium tier: corners x edges, corners x B-centres, corners x A-centres, and the
/// A-centre/edge table on the position and on its `SWAP_AB` conjugate.
pub struct MediumBounds<'a> {
    cc_ix: &'a CcIndex,
    cc: &'a Table2,
    cc_target: usize,
    ce_ix: &'a CeIndex,
    ce: &'a Table2,
    ce_target: usize,
    ae: &'a AeBound,
}

impl Bounds for MediumBounds<'_> {
    const K: usize = 5;
    #[inline(always)]
    fn index(&self, k: usize, p: &Pos) -> usize {
        match k {
            0 => self.ce_ix.index(p.c, &p.ep),
            1 => self.cc_ix.index_b(p.c, &p.xb),
            2 => self.cc_ix.index_a(p.c, &p.xa),
            3 => self.ae.index(p, false),
            _ => self.ae.index(p, true),
        }
    }
    #[inline(always)]
    fn table(&self, k: usize) -> &Table2 {
        match k {
            0 => self.ce,
            1 | 2 => self.cc,
            _ => &self.ae.tab,
        }
    }
    fn is_target(&self, k: usize, i: usize) -> bool {
        match k {
            0 => i == self.ce_target,
            1 | 2 => i == self.cc_target,
            _ => self.ae.is_target(i),
        }
    }
    fn name(&self, k: usize) -> &'static str {
        [
            "corners+edges",
            "corners+B",
            "corners+A",
            "AE (H_A)",
            "AE of SWAP_AB (H_B)",
        ][k.min(4)]
    }
}

/// The large tier: the medium tier's lookups, and the two-phase corner/B-centre/edge-triple
/// table (`cbe`, distance to `H_A` as a saturating code: at most 6, 7, 8, at least 9) on the
/// position and on its `SWAP_AB` conjugate, each read for all four A-faces' edge triples.
pub struct LargeBounds<'a> {
    m: MediumBounds<'a>,
    cbe: &'a CbeTable,
}

impl Bounds for LargeBounds<'_> {
    const K: usize = MediumBounds::K;
    #[inline(always)]
    fn index(&self, k: usize, p: &Pos) -> usize {
        self.m.index(k, p)
    }
    #[inline(always)]
    fn table(&self, k: usize) -> &Table2 {
        self.m.table(k)
    }
    fn is_target(&self, k: usize, i: usize) -> bool {
        self.m.is_target(k, i)
    }
    fn name(&self, k: usize) -> &'static str {
        self.m.name(k)
    }
    const S: usize = 2;
    const SAT_MAX: u8 = crate::cbe::CODE_LB[3];
    #[inline(always)]
    fn sat_index(&self, s: usize, p: &Pos, out: &mut [usize; 4]) {
        let x = self.m.ae.fc(p, s == 1);
        self.cbe.prefetch_all(self.m.ae.f, &x, out);
    }
    #[inline(always)]
    fn sat_bound(&self, _s: usize, idx: &[usize; 4]) -> u8 {
        self.cbe.bound(idx)
    }
}

/// Exact distance of lookup `k` of `p` (its table walked down to the solved entry: at each step a
/// neighbour one closer, which the mod-3 value identifies). `None` if the table is not a
/// distance table (no neighbour is one step closer, or no target within its largest distance).
pub fn walk<B: Bounds>(b: &B, k: usize, p: &Pos) -> Option<u8> {
    let mt = opt_moves();
    let tab = b.table(k);
    let max = tab.hist.len().saturating_sub(1);
    let mut y = *p;
    let mut d = 0u8;
    loop {
        let i = b.index(k, &y);
        if b.is_target(k, i) {
            return Some(d);
        }
        if usize::from(d) >= max {
            return None;
        }
        let want = (tab.get(i) + 2) % 3;
        let zs: [Pos; N_MOVES] = std::array::from_fn(|m| mt.apply(&y, m));
        let mut is = [0usize; N_MOVES];
        for m in 0..N_MOVES {
            is[m] = b.index(k, &zs[m]);
            prefetch_at(&tab.data, is[m] >> 5);
        }
        let m = (0..N_MOVES).find(|&m| tab.get(is[m]) == want)?;
        y = zs[m];
        d += 1;
    }
}

/// Exact distances of all lookups of `p` (`None` if a table is not a distance table).
pub fn exact<B: Bounds>(b: &B, p: &Pos) -> Option<[u8; MAX_K]> {
    let mut h = [0u8; MAX_K];
    for k in 0..B::K.min(MAX_K) {
        h[k] = walk(b, k, p)?;
    }
    Some(h)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coord::N_MSET;
    use rand::SeedableRng;

    /// The unreduced corner x B-centre graph, for checking the reduced table.
    struct RawCc(&'static OptMoves);

    impl Graph for RawCc {
        const MOVES: usize = N_MOVES;
        type Node = (u16, [u8; 16]);
        type Cache = ();
        fn node(&self, i: usize, (): &mut ()) -> (u16, [u8; 16]) {
            ((i / N_MSET) as u16, crate::fast::mset_labels()[i % N_MSET])
        }
        fn step(&self, &(c, x): &(u16, [u8; 16]), m: usize) -> usize {
            let (c2, x2) = self.0.apply_cb(c, &x, m);
            c2 as usize * N_MSET + mset_rank(&x2) as usize
        }
    }

    /// The symmetry-reduced corner x centre table agrees with the unreduced BFS on random
    /// positions and all their neighbours, through the B-centres and through the A-centres
    /// (the A lookup against the unreduced table of the `SWAP_AB` conjugate).
    #[test]
    #[cfg_attr(
        target_family = "wasm",
        ignore = "an unreduced BFS of 3.5e8 entries: minutes single-threaded under WASI"
    )]
    fn reduced_cc_matches_unreduced() {
        let mt = opt_moves();
        let ix = CcIndex::new();
        let red = build_cc(&ix, false);
        let s = Pos::solved();
        let raw = bfs(
            &RawCc(mt),
            crate::coord::N_CORNER * N_MSET,
            &[s.c as usize * N_MSET + mset_rank(&s.xb) as usize],
            "raw cc",
            false,
        );
        // (the reduced table counts indices, not positions: only the depths agree)
        assert_eq!(raw.hist.len(), red.hist.len());
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(77);
        for _ in 0..3000 {
            let c = Cube::random(&mut rng);
            for m in 0..=N_MOVES {
                let y = if m < N_MOVES { c.apply(m) } else { c };
                let p = Pos::of(&y);
                let raw_b = raw.get(p.c as usize * N_MSET + mset_rank(&p.xb) as usize);
                assert_eq!(red.get(ix.index_b(p.c, &p.xb)), raw_b);
                let w = Pos::of(&y.conjugate(&SWAP_AB).unwrap());
                let raw_a = raw.get(w.c as usize * N_MSET + mset_rank(&w.xb) as usize);
                assert_eq!(red.get(ix.index_a(p.c, &p.xa)), raw_a);
            }
        }
    }
}
