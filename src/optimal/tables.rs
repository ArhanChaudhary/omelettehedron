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

use super::coords::*;
use crate::cube::*;
use crate::preflight::{check_plan, Planned};
use crate::prune::{bfs, Graph, Table2};
use crate::tables::prefetch_at;
use std::path::{Path, PathBuf};

/// Most lookups per position of any tier (the search's scratch arrays have this many slots).
pub const MAX_K: usize = 4;

/// Table size of the optimal solver (`omelettehedron optimal ... --tier`).
#[derive(
    Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum, serde::Serialize, serde::Deserialize,
)]
pub enum OptTier {
    /// ~65 MB: corners x centres (one 4.7 MB table for both centre orbits) and the edge
    /// permutation (60 MB)
    Small,
    /// ~1.7 GB: corners x centres and corners x edges (48 symmetries; 64-bit targets only)
    Medium,
}

impl OptTier {
    pub fn name(self) -> &'static str {
        match self {
            OptTier::Small => "small",
            OptTier::Medium => "medium",
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
}

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
}

/// File names of the tables (all start with `opt-`, apart from the two-phase solver's).
const CC_FILE: &str = "opt-cc-sym24.prun";
const EDGES_FILE: &str = "opt-edges.prun";
const CE_FILE: &str = "opt-ce-sym48.prun";

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
            OptTier::Medium => {
                let ix = CeIndex::new();
                let t = build_ce(&ix, verbose);
                (None, Some((ix, t)))
            }
        };
        Ok(Self::assemble(tier, cc_ix, cc, edges, ce))
    }

    fn assemble(
        tier: OptTier,
        cc_ix: CcIndex,
        cc: Table2,
        edges: Option<Table2>,
        ce: Option<(CeIndex, Table2)>,
    ) -> OptTables {
        let cc_target = cc_ix.targets().first().copied().unwrap_or(0);
        OptTables {
            tier,
            cc_ix,
            cc,
            cc_target,
            edges,
            ce,
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
        let ce_ix = (tier == OptTier::Medium).then(CeIndex::new);
        let planned = |path: &PathBuf, entries: usize| {
            let bytes =
                std::fs::metadata(path).map_or_else(|_| entries.div_ceil(4) as u64, |m| m.len());
            Planned {
                path: path.clone(),
                bytes,
                working_set: 0,
            }
        };
        let mut plan = vec![planned(&cc_path, cc_size(&cc_ix))];
        match &ce_ix {
            None => plan.push(planned(&edges_path, N_EPERM)),
            Some(ix) => plan.push(planned(&ce_path, ix.size())),
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
        let t = Self::assemble(tier, cc_ix, cc, edges, ce);
        t.check_walks(&[&cc_path, &second])?;
        Ok(t)
    }

    /// The lower bounds of the small tier (`None` if the tables are another tier's).
    pub fn small(&self) -> Option<SmallBounds<'_>> {
        Some(SmallBounds {
            cc_ix: &self.cc_ix,
            cc: &self.cc,
            cc_target: self.cc_target,
            edges: self.edges.as_ref()?,
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
        })
    }

    /// A table file that passes its checksum may still have been made for another layout: walk
    /// down every table from a few random positions; each walk must reach the solved entry
    /// within the table's largest distance.
    ///
    /// # Errors
    ///
    /// If a walk does not reach the solved entry (`paths`: the corner x centre file and the
    /// tier's other file, for the message).
    fn check_walks(&self, paths: &[&Path]) -> std::io::Result<()> {
        use rand::SeedableRng;
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(0x0071_ab1e);
        let positions: Vec<Pos> = (0..8).map(|_| Pos::of(&Cube::random(&mut rng))).collect();
        // (lookup 0 is the tier's own table, 1 and 2 share the corner x centre table)
        let bad = |k: usize| -> std::io::Error {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "{}: not a distance table for this version's layout (made by a different \
                     version?); delete it to regenerate",
                    paths[usize::from(k == 0)].display()
                ),
            )
        };
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

/// The small tier: edges, corners x B-centres, corners x A-centres.
pub struct SmallBounds<'a> {
    cc_ix: &'a CcIndex,
    cc: &'a Table2,
    cc_target: usize,
    edges: &'a Table2,
}

impl Bounds for SmallBounds<'_> {
    const K: usize = 3;
    #[inline(always)]
    fn index(&self, k: usize, p: &Pos) -> usize {
        match k {
            0 => eperm_rank(&p.ep) as usize,
            1 => self.cc_ix.index_b(p.c, &p.xb),
            _ => self.cc_ix.index_a(p.c, &p.xa),
        }
    }
    #[inline(always)]
    fn table(&self, k: usize) -> &Table2 {
        if k == 0 {
            self.edges
        } else {
            self.cc
        }
    }
    fn is_target(&self, k: usize, i: usize) -> bool {
        if k == 0 {
            i == 0
        } else {
            i == self.cc_target
        }
    }
    fn name(&self, k: usize) -> &'static str {
        ["edges", "corners+B", "corners+A"][k.min(2)]
    }
}

/// The medium tier: corners x edges, corners x B-centres, corners x A-centres.
pub struct MediumBounds<'a> {
    cc_ix: &'a CcIndex,
    cc: &'a Table2,
    cc_target: usize,
    ce_ix: &'a CeIndex,
    ce: &'a Table2,
    ce_target: usize,
}

impl Bounds for MediumBounds<'_> {
    const K: usize = 3;
    #[inline(always)]
    fn index(&self, k: usize, p: &Pos) -> usize {
        match k {
            0 => self.ce_ix.index(p.c, &p.ep),
            1 => self.cc_ix.index_b(p.c, &p.xb),
            _ => self.cc_ix.index_a(p.c, &p.xa),
        }
    }
    #[inline(always)]
    fn table(&self, k: usize) -> &Table2 {
        if k == 0 {
            self.ce
        } else {
            self.cc
        }
    }
    fn is_target(&self, k: usize, i: usize) -> bool {
        if k == 0 {
            i == self.ce_target
        } else {
            i == self.cc_target
        }
    }
    fn name(&self, k: usize) -> &'static str {
        ["corners+edges", "corners+B", "corners+A"][k.min(2)]
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
