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

use super::cbe::{build_cbe, CbeIndex, N_LOOKUPS};
use super::coords::*;
use crate::cube::*;
use crate::fast::{fast, masks, FastTables};
use crate::preflight::{check_plan, MemGuard, Planned};
use crate::prune::{bfs, Graph, Table2, TableFile};
use crate::tables::{
    ae_targets, build_ae, prefetch_at, AeIndex, AeSym, AE_CHIRAL_ENTRIES, AE_CHIRAL_FILE,
    AE_TRIPLES_ENTRIES, AE_TRIPLES_FILE,
};
use std::path::{Path, PathBuf};

/// Most lookups per position of any tier (the search's scratch arrays have this many slots).
pub const MAX_K: usize = 10;

/// Table size of the optimal solver (`omelettehedron solve optimal --tier`).
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
    /// ~22.0 GB: the A-centre/edge table on both sides and a corner x B-centre x edge-triple
    /// table to solved (21.0 GB, eight lookups); 64-bit targets only
    Large,
}

impl OptTier {
    /// The tiers, smallest first.
    pub const ALL: [OptTier; 3] = [OptTier::Small, OptTier::Medium, OptTier::Large];

    pub fn name(self) -> &'static str {
        match self {
            OptTier::Small => "small",
            OptTier::Medium => "medium",
            OptTier::Large => "large",
        }
    }

    /// The table files of the tier in the metric of `moves` (the names [`OptTables::load`]
    /// reads), with their entries; `None` if the tier cannot be used with `moves` on this
    /// target (the large tier with slices, or a tier beyond small on a 32-bit target).
    pub fn files(self, moves: MoveSet) -> Option<Vec<TableFile>> {
        if unsupported(self, moves).is_some() {
            return None;
        }
        let sfx = moves.table_suffix();
        let file = |name: &str, entries: u64| TableFile {
            name: format!("{name}{sfx}.prun"),
            entries,
        };
        let (ae, own) = match self {
            OptTier::Small => (
                file(AE_TRIPLES_FILE, AE_TRIPLES_ENTRIES),
                file(EDGES_FILE, N_EPERM as u64),
            ),
            OptTier::Medium => (
                file(AE_CHIRAL_FILE, AE_CHIRAL_ENTRIES),
                file(CE_FILE, CE_ENTRIES),
            ),
            OptTier::Large => (
                file(AE_CHIRAL_FILE, AE_CHIRAL_ENTRIES),
                file(CBE_FILE, CBE_ENTRIES),
            ),
        };
        Some(vec![file(CC_FILE, CC_ENTRIES), ae, own])
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
    /// First step of [`Bounds::index`] for lookups whose index reads a table that misses the
    /// cache before the final entry is known (the A-centre/edge lookups' symmetry-class
    /// table): computes what that read needs and prefetches it. The search calls it for all
    /// children of a stage before any of their indices, and passes the result to
    /// [`Bounds::index_with`].
    #[inline(always)]
    fn pre(&self, _k: usize, _p: &Pos) -> u32 {
        0
    }
    /// [`Bounds::index`] of `p` given `pre(k, p)`.
    #[inline(always)]
    fn index_with(&self, k: usize, p: &Pos, _pre: u32) -> usize {
        self.index(k, p)
    }
    /// The order in which the search reads the lookups (the first `K` entries: a permutation of
    /// `0..K`); a child is dropped at the first value that exceeds its moves left, so the most
    /// selective lookups go first (measured with `FTO_OPT_PRUNE_STATS=1`).
    const ORDER: &'static [usize] = &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
    /// The move set (the metric of the tables, and the moves the search uses).
    fn moves(&self) -> MoveSet;
}

/// Moves between entries of the corner x B-centre table (`N` moves: 16, or 24 with slices).
struct CcGraph<'a, const N: usize> {
    ix: &'a CcIndex,
    mt: &'static OptMoves,
}

impl<const N: usize> Graph for CcGraph<'_, N> {
    const MOVES: usize = N;
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
struct EdgeGraph<const N: usize> {
    mt: &'static OptMoves,
}

impl<const N: usize> Graph for EdgeGraph<N> {
    const MOVES: usize = N;
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

/// The table's name in messages: `name`, with `-slice` for the metric with slices.
fn table_name(name: &str, moves: MoveSet) -> String {
    format!("{name}{}", moves.table_suffix())
}

/// The corner x centre table (distance to solved of the corners and the B-centres) in the
/// metric of `moves`. (The symmetry reduction holds with slices: every symmetry maps slices to
/// slices.)
pub fn build_cc(ix: &CcIndex, moves: MoveSet, verbose: bool) -> Table2 {
    let (mt, n) = (opt_moves(), cc_size(ix));
    let name = table_name("opt-cc", moves);
    if moves.slice {
        bfs(
            &CcGraph::<MAX_MOVES> { ix, mt },
            n,
            &ix.targets(),
            &name,
            verbose,
        )
    } else {
        bfs(
            &CcGraph::<N_MOVES> { ix, mt },
            n,
            &ix.targets(),
            &name,
            verbose,
        )
    }
}

/// The edge permutation table in the metric of `moves`.
pub fn build_edges(moves: MoveSet, verbose: bool) -> Table2 {
    let mt = opt_moves();
    let name = table_name("opt-edges", moves);
    if moves.slice {
        bfs(
            &EdgeGraph::<MAX_MOVES> { mt },
            N_EPERM,
            &[0],
            &name,
            verbose,
        )
    } else {
        bfs(&EdgeGraph::<N_MOVES> { mt }, N_EPERM, &[0], &name, verbose)
    }
}

/// Moves between entries of the corner x edge table.
struct CeGraph<'a, const N: usize> {
    ix: &'a CeIndex,
    mt: &'static OptMoves,
}

impl<const N: usize> Graph for CeGraph<'_, N> {
    const MOVES: usize = N;
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

/// The corner x edge table in the metric of `moves`.
pub fn build_ce(ix: &CeIndex, moves: MoveSet, verbose: bool) -> Table2 {
    let mt = opt_moves();
    let name = table_name("opt-ce", moves);
    if moves.slice {
        bfs(
            &CeGraph::<MAX_MOVES> { ix, mt },
            ix.size(),
            &[ix.target()],
            &name,
            verbose,
        )
    } else {
        bfs(
            &CeGraph::<N_MOVES> { ix, mt },
            ix.size(),
            &[ix.target()],
            &name,
            verbose,
        )
    }
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

    /// Index of the position (`swap == false`: distance to `H_A`) or of its `SWAP_AB`
    /// conjugate (`swap == true`: the position's distance to `H_B`).
    #[inline(always)]
    pub fn index(&self, p: &Pos, swap: bool) -> usize {
        let x = if swap {
            self.maps.fc_swap(p)
        } else {
            self.maps.fc(p)
        };
        self.ix.index_fast(self.f, &x, self.f.a_rank(&x))
    }

    fn is_target(&self, i: usize) -> bool {
        self.targets.contains(&i)
    }

    /// The A-centre rank of the position (or of its conjugate), its symmetry-class entry
    /// prefetched (see [`Bounds::pre`]).
    #[inline(always)]
    pub fn pre(&self, p: &Pos, swap: bool) -> u32 {
        let r = self.f.rank(&masks(&self.maps.a_of(p, swap)));
        self.ix.prefetch_class(r);
        r
    }

    /// [`AeBound::index`] given the A-centre rank from [`AeBound::pre`].
    #[inline(always)]
    pub fn index_with(&self, p: &Pos, swap: bool, a_rank: u32) -> usize {
        let x = if swap {
            self.maps.fc_swap(p)
        } else {
            self.maps.fc(p)
        };
        self.ix.index_fast(self.f, &x, a_rank)
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
    /// corners x edge permutation (medium and large tiers)
    pub ce: Option<(CeIndex, Table2)>,
    /// corners x B-centres x edge triple, and the indices of the solved state (large tier)
    pub cbe: Option<(CbeIndex, Table2, Vec<usize>)>,
    /// the two-phase A-centre/edge table (every tier)
    pub ae: AeBound,
    /// the move set: every table holds distances in its metric (with slices, the tables have
    /// their own files, `*-slice.prun`)
    pub moves: MoveSet,
}

/// File names of the tables (all start with `opt-`, apart from the two-phase solver's), before
/// the move set's suffix and `.prun`.
const CC_FILE: &str = "opt-cc-sym24";
const EDGES_FILE: &str = "opt-edges";
const CE_FILE: &str = "opt-ce-sym48";
const CBE_FILE: &str = "opt-cbe-sym6";

/// Entries of the tables (as the file headers give them), so that table files can be judged
/// without building the indices (`tier_entries_match_the_indices` checks them).
const CC_ENTRIES: u64 = 18_849_600;
const CE_ENTRIES: u64 = 9_580_032_000;
const CBE_ENTRIES: u64 = 83_913_984_000;

fn load_or_build(
    path: &Path,
    n: usize,
    verbose: bool,
    build: impl FnOnce() -> std::io::Result<Table2>,
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
            let t = build()?;
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

/// Why a tier cannot be used with a move set on this target, if it cannot: tiers beyond small
/// need 64-bit indices, and the large tier's triple table exists in the face metric only.
fn unsupported(tier: OptTier, moves: MoveSet) -> Option<std::io::Error> {
    if too_large(tier) {
        return Some(std::io::Error::other(
            "only --tier small is supported on 32-bit targets",
        ));
    }
    (tier == OptTier::Large && moves.slice).then(|| {
        std::io::Error::other(
            "the large tier has no tables for slice turns: use --tier medium with --moves slice",
        )
    })
}

impl OptTables {
    /// Build the tables of a tier in the metric of `moves` in memory, without touching the disk.
    ///
    /// # Errors
    ///
    /// If the tier needs 64-bit indices on a 32-bit target, or is the large tier with slices.
    pub fn generate(tier: OptTier, moves: MoveSet, verbose: bool) -> std::io::Result<OptTables> {
        if let Some(e) = unsupported(tier, moves) {
            return Err(e);
        }
        let cc_ix = CcIndex::new();
        let cc = build_cc(&cc_ix, moves, verbose);
        let (edges, ce) = match tier {
            OptTier::Small => (Some(build_edges(moves, verbose)), None),
            OptTier::Medium => {
                let ix = CeIndex::new();
                let t = build_ce(&ix, moves, verbose);
                (None, Some((ix, t)))
            }
            OptTier::Large => (None, None),
        };
        let cbe = match tier {
            OptTier::Large => {
                let ix = CbeIndex::new();
                let t = build_cbe_guarded(&ix, verbose)?;
                Some((ix, t))
            }
            _ => None,
        };
        let ae_ix = AeBound::index_for(tier);
        let ae_tab = build_ae(&ae_ix, moves, verbose);
        let ae = AeBound::new(ae_ix, ae_tab);
        Ok(Self::assemble(tier, cc_ix, cc, edges, ce, cbe, ae, moves))
    }

    fn assemble(
        tier: OptTier,
        cc_ix: CcIndex,
        cc: Table2,
        edges: Option<Table2>,
        ce: Option<(CeIndex, Table2)>,
        cbe: Option<(CbeIndex, Table2)>,
        ae: AeBound,
        moves: MoveSet,
    ) -> OptTables {
        let cc_target = cc_ix.targets().first().copied().unwrap_or(0);
        OptTables {
            tier,
            cc_ix,
            cc,
            cc_target,
            edges,
            ce,
            cbe: cbe.map(|(ix, t)| {
                let targets = ix.targets();
                (ix, t, targets)
            }),
            ae,
            moves,
        }
    }

    /// Load the tables of a tier in the metric of `moves` from `dir`, generating (and saving)
    /// the missing ones. All tables are checked against the free disk space and memory before
    /// any is generated.
    ///
    /// # Errors
    ///
    /// If a table file cannot be read or written, holds a table of the wrong format, size,
    /// checksum or layout, or the tables do not fit (see [`crate::preflight`]), or the tier
    /// needs 64-bit indices on a 32-bit target, or is the large tier with slices.
    pub fn load(
        dir: &Path,
        tier: OptTier,
        moves: MoveSet,
        verbose: bool,
    ) -> std::io::Result<OptTables> {
        if let Some(e) = unsupported(tier, moves) {
            return Err(e);
        }
        let cc_ix = CcIndex::new();
        let sfx = moves.table_suffix();
        let cc_path = dir.join(format!("{CC_FILE}{sfx}.prun"));
        let edges_path = dir.join(format!("{EDGES_FILE}{sfx}.prun"));
        let ce_path = dir.join(format!("{CE_FILE}{sfx}.prun"));
        let cbe_path = dir.join(format!("{CBE_FILE}{sfx}.prun"));
        let ce_ix = (tier == OptTier::Medium).then(CeIndex::new);
        let cbe_ix = (tier == OptTier::Large).then(CbeIndex::new);
        let ae_ix = AeBound::index_for(tier);
        // (the two-phase solver's file of the same table: generated by either solver, shared)
        let ae_path = dir.join(format!("{}{sfx}.prun", ae_ix.name()));
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
        match (&ce_ix, tier) {
            (Some(ix), _) => plan.push(planned(&ce_path, ix.size())),
            (None, OptTier::Small) => plan.push(planned(&edges_path, N_EPERM)),
            (None, _) => {}
        }
        if let Some(ix) = &cbe_ix {
            plan.push(planned(&cbe_path, ix.size()));
        }
        check_plan(&plan)?;
        let cc = load_or_build(&cc_path, cc_size(&cc_ix), verbose, || {
            Ok(build_cc(&cc_ix, moves, verbose))
        })?;
        // (the large tier reads neither the edge nor the corner x edge table: its own file is
        // the triple table)
        let (edges, ce, second) = match (ce_ix, tier) {
            (Some(ix), _) => {
                let t = load_or_build(&ce_path, ix.size(), verbose, || {
                    Ok(build_ce(&ix, moves, verbose))
                })?;
                (None, Some((ix, t)), ce_path)
            }
            (None, OptTier::Small) => (
                Some(load_or_build(&edges_path, N_EPERM, verbose, || {
                    Ok(build_edges(moves, verbose))
                })?),
                None,
                edges_path,
            ),
            (None, _) => (None, None, cbe_path.clone()),
        };
        let cbe = match cbe_ix {
            None => None,
            Some(ix) => {
                let t = load_or_build(&cbe_path, ix.size(), verbose, || {
                    build_cbe_guarded(&ix, verbose)
                })?;
                Some((ix, t))
            }
        };
        let ae_tab = load_or_build(&ae_path, ae_ix.size(), verbose, || {
            Ok(build_ae(&ae_ix, moves, verbose))
        })?;
        let ae = AeBound::new(ae_ix, ae_tab);
        let t = Self::assemble(tier, cc_ix, cc, edges, ce, cbe, ae, moves);
        t.check_walks(&[&second, &cc_path, &ae_path, &cbe_path])?;
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
            moves: self.moves,
        })
    }

    /// The lower bounds of the large tier (`None` if the tables are another tier's).
    pub fn large(&self) -> Option<LargeBounds<'_>> {
        let (cbe_ix, cbe, cbe_targets) = self.cbe.as_ref()?;
        Some(LargeBounds {
            ae: &self.ae,
            cbe_ix,
            cbe,
            cbe_targets,
            moves: self.moves,
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
            moves: self.moves,
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
        // (the file of a failing lookup: 0 the tier's own table, 1 the corner x centre table, 2
        // the A-centre/edge table, 3 the large tier's triple table)
        let bad = |f: usize| -> std::io::Error {
            let file = paths[f.min(paths.len().saturating_sub(1))];
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "{}: not a distance table for this version's layout (made by a different \
                     version?); delete it to regenerate",
                    file.display()
                ),
            )
        };
        // small and medium: lookup 0 the tier's own table, 1 and 2 the corner x centre table,
        // 3 and 4 the A-centre/edge table
        let medium_file = |k: usize| match k {
            0 => 0,
            1 | 2 => 1,
            _ => 2,
        };
        // large: lookups 0 and 1 the A-centre/edge table, 2.. the triples
        let large_file = |k: usize| if k < 2 { 2 } else { 3 };
        let r = match self.tier {
            OptTier::Small => self
                .small()
                .map(|b| check_walks_of(&b, &positions).map_err(medium_file)),
            OptTier::Medium => self
                .medium()
                .map(|b| check_walks_of(&b, &positions).map_err(medium_file)),
            OptTier::Large => self
                .large()
                .map(|b| check_walks_of(&b, &positions).map_err(large_file)),
        }
        .unwrap_or(Err(0));
        r.map_err(bad)
    }

    /// Bytes of the tables in memory.
    pub fn bytes(&self) -> usize {
        self.cc.bytes()
            + self.edges.as_ref().map_or(0, Table2::bytes)
            + self.ce.as_ref().map_or(0, |(_, t)| t.bytes())
            + self.cbe.as_ref().map_or(0, |(_, t, _)| t.bytes())
            + self.ae.tab.bytes()
    }
}

/// The corner x B-centre x edge-triple table, stopped cleanly (an error, nothing written) if
/// the free memory falls below the floor while it is generated ([`MemGuard`]).
fn build_cbe_guarded(ix: &CbeIndex, verbose: bool) -> std::io::Result<Table2> {
    let guard = MemGuard::start();
    build_cbe(ix, verbose, guard.flag()).ok_or_else(|| MemGuard::error(CBE_FILE))
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
    moves: MoveSet,
}

impl Bounds for SmallBounds<'_> {
    const K: usize = 5;
    fn moves(&self) -> MoveSet {
        self.moves
    }
    /// AE, its `SWAP_AB` twin, the corner x centre table (cache-resident), the edge table: on
    /// the 15-move set (seed 1717) AE alone prunes 78% of the children, the twin 76%, the edge
    /// and corner x centre tables 16-18%; 1.46 lookups per node instead of 5.00.
    const ORDER: &'static [usize] = &[3, 4, 1, 2, 0];
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
    fn pre(&self, k: usize, p: &Pos) -> u32 {
        match k {
            3 => self.ae.pre(p, false),
            4 => self.ae.pre(p, true),
            _ => 0,
        }
    }
    #[inline(always)]
    fn index_with(&self, k: usize, p: &Pos, pre: u32) -> usize {
        match k {
            3 => self.ae.index_with(p, false, pre),
            4 => self.ae.index_with(p, true, pre),
            _ => self.index(k, p),
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
    moves: MoveSet,
}

impl Bounds for MediumBounds<'_> {
    const K: usize = 5;
    fn moves(&self) -> MoveSet {
        self.moves
    }
    /// AE, its `SWAP_AB` twin, corners x edges, the corner x centre table: on the 15-move set
    /// (seed 1717) AE alone prunes 74% of the children, the twin 72%, corners x edges 61%, the
    /// corner x centre lookups 1.4% each; 1.52 lookups per node instead of 5.00.
    const ORDER: &'static [usize] = &[3, 4, 0, 1, 2];
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
    fn pre(&self, k: usize, p: &Pos) -> u32 {
        match k {
            3 => self.ae.pre(p, false),
            4 => self.ae.pre(p, true),
            _ => 0,
        }
    }
    #[inline(always)]
    fn index_with(&self, k: usize, p: &Pos, pre: u32) -> usize {
        match k {
            3 => self.ae.index_with(p, false, pre),
            4 => self.ae.index_with(p, true, pre),
            _ => self.index(k, p),
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

/// The large tier: the A-centre/edge table on the position and on its `SWAP_AB` conjugate,
/// then the corner x B-centre x edge-triple table on the four A-triples and, through
/// `SWAP_AB`, the four B-triples. The medium tier's corner x centre lookups are left out: an
/// A-triple lookup projects onto the corners and B-centres and more, so it is never below the
/// corner x B-centre value (and a B-triple lookup never below the corner x A-centre value).
/// So is corners x edges (2.4 GB): after the others it prunes 0.3% of the children, which costs
/// more time than it saves (without it, against the tier with it: 0.979x time on one thread,
/// 0.954x on 48, 1.001-1.028x the nodes), and the tier needs 22.0 GB instead of 24.4.
pub struct LargeBounds<'a> {
    ae: &'a AeBound,
    cbe_ix: &'a CbeIndex,
    cbe: &'a Table2,
    cbe_targets: &'a [usize],
    moves: MoveSet,
}

/// The large tier's first corner x B-centre x edge-triple lookup (0 and 1 are the A-centre/edge
/// table on the position and on its conjugate).
const LARGE_CBE: usize = 2;

impl Bounds for LargeBounds<'_> {
    const K: usize = LARGE_CBE + N_LOOKUPS;
    /// AE, the B-triple 3 lookup, the AE twin, then the other triples: greedy by the children
    /// each prunes among those the earlier ones keep (17-move set, seed 1717: alone AE prunes
    /// 65%, its twin 58%, the triple lookups 48-59%; 7.8% survive). 2.27 lookups per node on the
    /// 17-move set instead of 10. (The triples differ because of the fixed front corner.)
    const ORDER: &'static [usize] = &[0, 9, 1, 6, 4, 7, 8, 5, 3, 2];
    fn moves(&self) -> MoveSet {
        self.moves
    }
    #[inline(always)]
    fn index(&self, k: usize, p: &Pos) -> usize {
        if k < LARGE_CBE {
            self.ae.index(p, k == 1)
        } else {
            self.cbe_ix.index(k - LARGE_CBE, p)
        }
    }
    #[inline(always)]
    fn pre(&self, k: usize, p: &Pos) -> u32 {
        if k < LARGE_CBE {
            self.ae.pre(p, k == 1)
        } else {
            0
        }
    }
    #[inline(always)]
    fn index_with(&self, k: usize, p: &Pos, pre: u32) -> usize {
        if k < LARGE_CBE {
            self.ae.index_with(p, k == 1, pre)
        } else {
            self.cbe_ix.index(k - LARGE_CBE, p)
        }
    }
    #[inline(always)]
    fn table(&self, k: usize) -> &Table2 {
        if k < LARGE_CBE {
            &self.ae.tab
        } else {
            self.cbe
        }
    }
    fn is_target(&self, k: usize, i: usize) -> bool {
        if k < LARGE_CBE {
            self.ae.is_target(i)
        } else {
            self.cbe_targets.contains(&i)
        }
    }
    fn name(&self, k: usize) -> &'static str {
        [
            "AE (H_A)",
            "AE of SWAP_AB (H_B)",
            "CBE A-triple 0",
            "CBE A-triple 1",
            "CBE A-triple 2",
            "CBE A-triple 3",
            "CBE B-triple 0",
            "CBE B-triple 1",
            "CBE B-triple 2",
            "CBE B-triple 3",
        ][k.min(9)]
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
        let n = b.moves().n_moves();
        let zs: [Pos; MAX_MOVES] = std::array::from_fn(|m| mt.apply(&y, m));
        let mut is = [0usize; MAX_MOVES];
        for m in 0..n {
            is[m] = b.index(k, &zs[m]);
            prefetch_at(&tab.data, is[m] >> 5);
        }
        let m = (0..n).find(|&m| tab.get(is[m]) == want)?;
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

    /// The entries [`OptTier::files`] expects are the indices' sizes, and the names are the
    /// ones [`OptTables::load`] reads. (64-bit only: the larger indices exceed 2^32 entries.)
    #[cfg(target_pointer_width = "64")]
    #[test]
    fn tier_entries_match_the_indices() {
        assert_eq!(cc_size(&CcIndex::new()) as u64, CC_ENTRIES);
        assert_eq!(CeIndex::new().size() as u64, CE_ENTRIES);
        assert_eq!(CbeIndex::new().size() as u64, CBE_ENTRIES);
        let names = |t: OptTier, slice: bool| -> Option<Vec<String>> {
            let moves = MoveSet {
                slice,
                ..MoveSet::FACE
            };
            Some(t.files(moves)?.into_iter().map(|f| f.name).collect())
        };
        assert_eq!(
            names(OptTier::Small, false).unwrap(),
            [
                "opt-cc-sym24.prun",
                "ae-sym288-triples.prun",
                "opt-edges.prun"
            ]
        );
        assert_eq!(
            names(OptTier::Medium, true).unwrap(),
            [
                "opt-cc-sym24-slice.prun",
                "ae-sym288-chiral-slice.prun",
                "opt-ce-sym48-slice.prun"
            ]
        );
        assert_eq!(
            names(OptTier::Large, false).unwrap(),
            [
                "opt-cc-sym24.prun",
                "ae-sym288-chiral.prun",
                "opt-cbe-sym6.prun"
            ]
        );
        // the large tier has no slice tables
        assert_eq!(names(OptTier::Large, true), None);
    }

    /// Each tier's lookup order is a permutation of its lookups.
    #[test]
    fn lookup_orders_are_permutations() {
        fn check<B: Bounds>() {
            let mut o: Vec<usize> = B::ORDER.iter().take(B::K).copied().collect();
            o.sort_unstable();
            assert_eq!(o, (0..B::K).collect::<Vec<_>>());
        }
        check::<SmallBounds>();
        check::<MediumBounds>();
        check::<LargeBounds>();
    }

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
        let red = build_cc(&ix, MoveSet::FACE, false);
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
