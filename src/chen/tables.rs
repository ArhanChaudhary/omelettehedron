//! Phase-2 pruning tables of the block split.
//!
//! A phase-2 position is its configuration ([`super::space`]), its corner state and its
//! B-centre colouring. Three tables bound its distance to solved with the allowed moves:
//!
//! * **AC** (`chen-acs.prun`): configuration x corner state, in the pattern view under the 24
//!   symmetries (1,611 classes x 11,520 corner states = 18.6M entries, 4.6 MB): the exact
//!   distance of that projection with the allowed moves;
//! * **SB** (`chen-sb.prun`): pattern-view shape x B-centre colouring (460 x 369,600 = 170M
//!   entries, 42.5 MB), the same in the pattern view, read through the frame of the
//!   configuration (whether a move is allowed depends only on the shape, so this is a
//!   projection of the restricted move graph too; in the fixed-corner frame the shapes are 12
//!   times as many);
//! * **CC**: the optimal solver's corner x B-centre table (`opt-cc-sym24.prun`, 4.7 MB),
//!   distances to solved with every move allowed, which are lower bounds of the restricted
//!   ones.
//!
//! Each holds the distance mod 3 (two bits, see [`crate::prune`]): every allowed move can be
//! undone by an allowed move, so one move changes each distance by at most one, and the search
//! decodes a child's distance from its parent's.

use super::corners::{PCorners, N_PCORNER, N_SYM};
use super::space::*;
use crate::coord::{label_perm_table, move_tables, N_CORNER, N_MSET};
use crate::cube::*;
use crate::fast::{fast, masks, mset_labels, shuffle, FastTables};
use crate::optimal::coords::{cc_size, opt_moves, CcIndex};
use crate::optimal::tables::{build_cc, CC_FILE};
use crate::prune::{bfs, Graph, Table2};
use crate::tables::load_or_build;
use std::path::Path;

/// File names of the tables (with `.prun`).
pub const AC_FILE: &str = "chen-acs.prun";
pub const SB_FILE: &str = "chen-sb.prun";
pub const PX_FILE: &str = "chen-px.prun";

/// Entries of the tables.
/// (1,611 symmetry classes of pattern-view configurations x 11,520 pattern-view corner states)
pub const N_AC: usize = 1611 * N_PCORNER;
pub const N_SB: usize = N_PSHAPE * N_MSET;
/// (1,611 symmetry classes of pattern-view configurations, checked by the tests)
pub const N_PX: usize = 1611 * N_MSET;

/// A phase-2 position: configuration id, corner id and B-centre colouring (labels in lanes
/// 0..12).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct P2 {
    /// configuration id ([`Space`])
    pub cfg: u32,
    /// corner id ([`crate::coord::Corners`])
    pub c: u16,
    /// B-centre colouring (labels 0..4 in lanes 0..12, zero padding)
    pub xb: [u8; 16],
}

impl P2 {
    /// The phase-2 position of `cube`, if its configuration is reachable.
    pub fn of(cube: &Cube) -> Option<P2> {
        let mut xb = [0u8; 16];
        xb[..12].copy_from_slice(&cube.xb);
        Some(P2 {
            cfg: space().id(cube)?,
            c: crate::coord::corners().id(cube),
            xb,
        })
    }

    /// The position after move `m`, if the move is allowed.
    #[inline]
    pub fn apply(&self, s: &Space, m: usize) -> Option<P2> {
        let cfg = *s.next.get(self.cfg as usize)?.get(m & 15)?;
        if cfg == NONE {
            return None;
        }
        let (c, xb) = opt_moves().apply_cb(self.c, &self.xb, m);
        Some(P2 { cfg, c, xb })
    }
}

/// Moves between AC entries (index `class * N_PCORNER + corner state`, the corner state
/// conjugated like the configuration to its class representative), in the pattern view. It
/// holds the distances of the configuration x corner id table of the fixed frame (447,120 x 960
/// entries; the tests compare the two) in 4.6 MB instead of 107.3 MB, made in 1.4 s instead
/// of 25.7 s on one thread (0.1 s and 0.7 s on 48); its index costs more: 1.03x the time to
/// the same first solutions of 4 random positions on one thread (3.58 -> 3.70 s, two rounds).
struct AcGraph<'a> {
    pat: &'a Patterns,
    pc: &'a PCorners,
}

impl AcGraph<'_> {
    #[inline]
    fn index(&self, q: usize, pc: usize) -> usize {
        let (cls, k) = self.pat.class[q];
        cls as usize * N_PCORNER + self.pc.sym[pc][k as usize % N_SYM] as usize
    }
}

impl Graph for AcGraph<'_> {
    const MOVES: usize = N_MOVES;
    type Node = (usize, usize, usize);
    type Cache = ();
    fn node(&self, i: usize, (): &mut ()) -> (usize, usize, usize) {
        (i, self.pat.rep[i / N_PCORNER] as usize, i % N_PCORNER)
    }
    #[inline]
    fn step(&self, &(i, p, pc): &(usize, usize, usize), m: usize) -> usize {
        let q = self.pat.next[p][m];
        if q == u16::MAX {
            return i;
        }
        self.index(q as usize, self.pc.next[pc][m] as usize)
    }
    fn equivalents(&self, j: usize, out: &mut Vec<usize>) {
        let cls = j / N_PCORNER;
        let pc = j % N_PCORNER;
        for &k in &self.pat.stab[cls] {
            out.push(cls * N_PCORNER + self.pc.sym[pc][k as usize % N_SYM] as usize);
        }
    }
}

/// The AC graph of the space.
fn ac_graph() -> AcGraph<'static> {
    AcGraph {
        pat: &space().pat,
        pc: &space().corners,
    }
}

/// Moves between SB entries (index `pattern shape * N_MSET + B-colouring rank`), the
/// pattern-view face turns.
struct SbGraph<'a> {
    s: &'a Space,
    b: &'a [[u32; N_PMOVES]],
}

impl Graph for SbGraph<'_> {
    const MOVES: usize = N_PMOVES;
    type Node = (usize, usize);
    type Cache = ();
    fn node(&self, i: usize, (): &mut ()) -> (usize, usize) {
        (i / N_MSET, i % N_MSET)
    }
    #[inline]
    fn step(&self, &(sh, b): &(usize, usize), m: usize) -> usize {
        let n = self.s.pshape_next[sh][m];
        if n == u16::MAX {
            sh * N_MSET + b
        } else {
            n as usize * N_MSET + self.b[b][m] as usize
        }
    }
}

/// Moves between PX entries (index `class * N_MSET + B-colouring rank`, the colouring
/// conjugated like the configuration to its class representative), in the pattern view.
struct PxGraph<'a> {
    pat: &'a Patterns,
    f: &'a FastTables,
    labels: &'a [[u8; 16]],
    xb_src: [[u8; 16]; N_MOVES],
}

impl PxGraph<'_> {
    #[inline]
    fn index(&self, q: usize, xb: &[u8; 16]) -> usize {
        let (cls, k) = self.pat.class[q];
        cls as usize * N_MSET + self.f.rank(&masks(&self.pat.sym_xb(k as usize, xb))) as usize
    }
}

impl Graph for PxGraph<'_> {
    const MOVES: usize = N_MOVES;
    type Node = (usize, usize, [u8; 16]);
    type Cache = ();
    fn node(&self, i: usize, (): &mut ()) -> (usize, usize, [u8; 16]) {
        let p = self.pat.rep[i / N_MSET] as usize;
        (i, p, self.labels[i % N_MSET])
    }
    #[inline]
    fn step(&self, &(i, p, xb): &(usize, usize, [u8; 16]), m: usize) -> usize {
        let q = self.pat.next[p][m];
        if q == u16::MAX {
            return i;
        }
        self.index(q as usize, &shuffle(&xb, &self.xb_src[m]))
    }
    fn equivalents(&self, j: usize, out: &mut Vec<usize>) {
        let cls = j / N_MSET;
        let xb = self.labels[j % N_MSET];
        for &k in &self.pat.stab[cls] {
            out.push(
                cls * N_MSET + self.f.rank(&masks(&self.pat.sym_xb(k as usize, &xb))) as usize,
            );
        }
    }
}

/// The PX graph of the space.
fn px_graph() -> PxGraph<'static> {
    PxGraph {
        pat: &space().pat,
        f: fast(),
        labels: mset_labels(),
        xb_src: std::array::from_fn(|m| {
            let mut g = [0x80u8; 16];
            g[..12].copy_from_slice(&space().pat.xb_src[m]);
            g
        }),
    }
}

/// The PX table, by breadth-first search from solved in the pattern view.
fn build_px(verbose: bool) -> Table2 {
    let g = px_graph();
    let s0 = Cube::solved();
    let mut xb = [0u8; 16];
    xb[..12].copy_from_slice(&s0.xb);
    let t = g.index(0, &xb);
    let mut targets = vec![t];
    g.equivalents(t, &mut targets);
    bfs(&g, N_PX, &targets, "chen-px", verbose)
}

/// Why the phase-2 tables cannot be loaded. Errors of a file carry its path.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum ChenTableError {
    /// Reading, writing or checking the file failed (wrong size, checksum or format included).
    #[error("{path}")]
    Io {
        /// the file
        path: String,
        /// what went wrong
        source: std::io::Error,
    },
    /// The files pass their checksums but are not distance tables of this layout (made by
    /// another version): walking down them does not reach solved.
    #[error(
        "{dir}: a chen-*.prun table is not a table of this version's layout (delete it to \
         regenerate)"
    )]
    WrongLayout {
        /// the table directory
        dir: String,
    },
}

/// Which optional phase-2 tables to use (AC and CC are always used).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TableSet {
    /// the SB table (pattern-view shape x B-centres)
    pub sb: bool,
    /// the PX table (pattern-view configuration x B-centres)
    pub px: bool,
}

/// Number of phase-2 tables ([`ChenTables::exact`] gives their distances in this order: AC,
/// SB, PX, CC; a table not loaded counts 0).
pub const N_TABLES: usize = 4;

/// The phase-2 tables, with what the search needs to use them.
pub struct ChenTables {
    /// configuration x corners
    pub ac: Table2,
    /// pattern-view shape x B-centres, if loaded
    pub sb: Option<Table2>,
    /// pattern-view configuration x B-centres (24 symmetries), if loaded
    pub px: Option<Table2>,
    /// corners x B-centres with every move allowed (the optimal solver's table)
    pub cc: Table2,
    /// the index of `cc`
    pub cc_ix: CcIndex,
    /// index of the solved state in each table (AC, PX and CC: and its equivalents)
    ac_targets: Vec<usize>,
    sb_target: usize,
    px_targets: Vec<usize>,
    cc_targets: Vec<usize>,
    /// the B-centre colouring of `H_A` for each corner id (in `H_A` it follows from the
    /// corners)
    pub corner_xb: Vec<[u8; 16]>,
}

/// The B-centre colouring that goes with each corner state in `H_A` (by breadth-first search
/// over the class-A moves).
fn corner_xb() -> Vec<[u8; 16]> {
    let cs = crate::coord::corners();
    let s0 = Cube::solved();
    let mut out = vec![[0u8; 16]; N_CORNER];
    let mut seen = vec![false; N_CORNER];
    let mut q = vec![s0];
    seen[cs.id(&s0) as usize] = true;
    out[cs.id(&s0) as usize][..12].copy_from_slice(&s0.xb);
    let mut i = 0;
    while i < q.len() {
        let c = q[i];
        i += 1;
        for m in (0..N_MOVES).filter(|&m| is_a_move(m)) {
            let n = c.apply(m);
            let id = cs.id(&n) as usize;
            if !seen[id] {
                seen[id] = true;
                out[id][..12].copy_from_slice(&n.xb);
                q.push(n);
            }
        }
    }
    out
}

/// The solved state as a phase-2 position.
fn solved_p2() -> P2 {
    let s0 = Cube::solved();
    let mut xb = [0u8; 16];
    xb[..12].copy_from_slice(&s0.xb);
    P2 { cfg: 0, c: 0, xb }
}

/// The AC table, by breadth-first search from solved in the pattern view.
fn build_ac(verbose: bool) -> Table2 {
    let g = ac_graph();
    let t = g.index(0, usize::from(space().corners.solved));
    let mut targets = vec![t];
    g.equivalents(t, &mut targets);
    bfs(&g, N_AC, &targets, "chen-acs", verbose)
}

/// The SB table, by breadth-first search from solved in the pattern view.
fn build_sb(verbose: bool) -> Table2 {
    let b = label_perm_table(&pattern_xb_moves());
    let g = SbGraph { s: space(), b: &b };
    let t = ChenTables::sb_index_with(space(), fast(), &solved_p2());
    bfs(&g, N_SB, &[t], "chen-sb", verbose)
}

impl ChenTables {
    /// Build the tables in memory without touching the disk.
    pub fn generate(set: TableSet, verbose: bool) -> ChenTables {
        let cc_ix = CcIndex::new();
        let cc = build_cc(&cc_ix, MoveSet::FACE, verbose);
        Self::assemble(
            build_ac(verbose),
            set.sb.then(|| build_sb(verbose)),
            set.px.then(|| build_px(verbose)),
            cc,
            cc_ix,
        )
    }

    /// Load the tables of `set` from `dir`, generating (and saving) any that is missing.
    ///
    /// # Errors
    ///
    /// [`ChenTableError::Io`] if a table file cannot be read or written, or holds a table of
    /// the wrong size or checksum; [`ChenTableError::WrongLayout`] if the tables are of
    /// another layout.
    pub fn load(dir: &Path, set: TableSet, verbose: bool) -> Result<ChenTables, ChenTableError> {
        // (the configuration space and the move tables are built before any table: their
        // initialisation runs parallel work, which must not start inside the generator's)
        let _ = (space(), move_tables());
        let load = |name: &str, n: usize, build: &dyn Fn() -> Table2| {
            let path = dir.join(name);
            load_or_build(&path, n, verbose, build).map_err(|source| ChenTableError::Io {
                path: path.display().to_string(),
                source,
            })
        };
        let ac = load(AC_FILE, N_AC, &|| build_ac(verbose))?;
        let sb = if set.sb {
            Some(load(SB_FILE, N_SB, &|| build_sb(verbose))?)
        } else {
            None
        };
        let px = if set.px {
            Some(load(PX_FILE, N_PX, &|| build_px(verbose))?)
        } else {
            None
        };
        let cc_ix = CcIndex::new();
        let cc = load(&format!("{CC_FILE}.prun"), cc_size(&cc_ix), &|| {
            build_cc(&cc_ix, MoveSet::FACE, verbose)
        })?;
        let t = Self::assemble(ac, sb, px, cc, cc_ix);
        // a file of another layout passes its checksum but is no distance table: walk down
        // from a few positions
        if !t.check_walks() {
            return Err(ChenTableError::WrongLayout {
                dir: dir.display().to_string(),
            });
        }
        Ok(t)
    }

    fn assemble(
        ac: Table2,
        sb: Option<Table2>,
        px: Option<Table2>,
        cc: Table2,
        cc_ix: CcIndex,
    ) -> ChenTables {
        let (s, f) = (space(), fast());
        let solved = solved_p2();
        let g = px_graph();
        let t = Self::px_index_with(s, f, &solved);
        let mut px_targets = vec![t];
        g.equivalents(t, &mut px_targets);
        let ga = ac_graph();
        let t = Self::ac_index_with(s, &solved);
        let mut ac_targets = vec![t];
        ga.equivalents(t, &mut ac_targets);
        let cc_targets = cc_ix.targets();
        ChenTables {
            ac,
            sb,
            px,
            cc,
            cc_ix,
            ac_targets,
            sb_target: Self::sb_index_with(s, f, &solved),
            px_targets,
            cc_targets,
            corner_xb: corner_xb(),
        }
    }

    /// Whether walking down the tables reaches solved from a few random phase-2 positions.
    fn check_walks(&self) -> bool {
        use rand::{Rng, SeedableRng};
        let s = space();
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(0xc4e7);
        let mut p = solved_p2();
        (0..8).all(|_| {
            for _ in 0..40 {
                let m = rng.gen_range(0..N_MOVES);
                p = p.apply(s, m).unwrap_or(p);
            }
            self.exact(&p).is_some()
        })
    }

    /// AC index of a position: the class of its pattern-view configuration and its corners
    /// seen through its frame, conjugated to the class representative.
    #[inline]
    pub fn ac_index_with(s: &Space, p: &P2) -> usize {
        let pat = &s.pat;
        let cfg = p.cfg as usize;
        let q = pat.of_cfg.get(cfg).copied().unwrap_or(0) as usize;
        let (cls, k) = pat.class.get(q).copied().unwrap_or((0, 0));
        let r = s.rot.get(cfg).copied().unwrap_or(0) as usize;
        let pc = s
            .corners
            .of_fixed
            .get(r * N_CORNER + p.c as usize)
            .copied()
            .unwrap_or(0) as usize;
        let pk = s
            .corners
            .sym
            .get(pc)
            .map_or(0, |row| row[k as usize % N_SYM]);
        cls as usize * N_PCORNER + pk as usize
    }

    /// SB index of a position: its pattern-view shape and B-centres seen through its frame.
    #[inline]
    pub fn sb_index_with(s: &Space, f: &FastTables, p: &P2) -> usize {
        let ps = s.pshape.get(p.cfg as usize).copied().unwrap_or(0);
        ps as usize * N_MSET + f.rank(&masks(&s.pattern_xb(p.cfg, &p.xb))) as usize
    }

    /// PX index of a position: the class of its pattern-view configuration and its B-centres
    /// seen through its frame, conjugated to the class representative.
    #[inline]
    pub fn px_index_with(s: &Space, f: &FastTables, p: &P2) -> usize {
        let pat = &s.pat;
        let q = pat.of_cfg.get(p.cfg as usize).copied().unwrap_or(0) as usize;
        let (cls, k) = pat.class.get(q).copied().unwrap_or((0, 0));
        let x = pat.sym_xb(k as usize, &s.pattern_xb(p.cfg, &p.xb));
        cls as usize * N_MSET + f.rank(&masks(&x)) as usize
    }

    /// CC index of a position.
    #[inline]
    pub fn cc_index(&self, p: &P2) -> usize {
        self.cc_ix.index_b(p.c, &p.xb)
    }

    /// Whether a position is in `H_A` (configuration in `H'`, B-centres attached to the
    /// corners), and then its `H_A` table index.
    #[inline]
    pub fn ha_index(&self, s: &Space, p: &P2) -> Option<usize> {
        let e = *s.p2e.get(p.cfg as usize)?;
        (e != NONE && self.corner_xb.get(p.c as usize)? == &p.xb)
            .then(|| p.c as usize * crate::tables::N_P2_EDGE + e as usize)
    }

    /// Exact distances (AC, SB, PX, CC; 0 for a table not loaded) of a position, by walking
    /// down each table. `None` if a table is not a distance table of this layout (a walk does
    /// not reach its target).
    pub fn exact(&self, p: &P2) -> Option<[u8; N_TABLES]> {
        let (s, f) = (space(), fast());
        let ac = walk(
            &self.ac,
            p,
            |q| Self::ac_index_with(s, q),
            |i| self.ac_targets.contains(&i),
            true,
        )?;
        let sb = match &self.sb {
            Some(t) => walk(
                t,
                p,
                |q| Self::sb_index_with(s, f, q),
                |i| i == self.sb_target,
                true,
            )?,
            None => 0,
        };
        let px = match &self.px {
            Some(t) => walk(
                t,
                p,
                |q| Self::px_index_with(s, f, q),
                |i| self.px_targets.contains(&i),
                true,
            )?,
            None => 0,
        };
        let cc = walk(
            &self.cc,
            p,
            |q| self.cc_index(q),
            |i| self.cc_targets.contains(&i),
            false,
        )?;
        Some([ac, sb, px, cc])
    }

    /// Bytes of the tables.
    pub fn bytes(&self) -> usize {
        self.ac.bytes()
            + self.sb.as_ref().map_or(0, Table2::bytes)
            + self.px.as_ref().map_or(0, Table2::bytes)
            + self.cc.bytes()
    }
}

/// Distance of `p` in `tab` by walking down it: each step goes to a neighbour one closer
/// (value one less mod 3) until `done`. With `restricted`, only the allowed moves are walked
/// (AC and SB are distances with them); CC is a distance with every move, so its walk moves
/// the corners and B-centres only.
fn walk(
    tab: &Table2,
    p: &P2,
    index: impl Fn(&P2) -> usize,
    done: impl Fn(usize) -> bool,
    restricted: bool,
) -> Option<u8> {
    let s = space();
    let om = opt_moves();
    let mut q = *p;
    let max = tab.hist.len();
    for d in 0..=max {
        let i = index(&q);
        if done(i) {
            return Some(d as u8);
        }
        let want = (tab.get(i) + 2) % 3;
        let step = (0..N_MOVES).find_map(|m| {
            let n = if restricted {
                q.apply(s, m)?
            } else {
                let (c, xb) = om.apply_cb(q.c, &q.xb, m);
                P2 { c, xb, ..q }
            };
            (tab.get(index(&n)) == want).then_some(n)
        })?;
        q = step;
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{Rng, SeedableRng};

    /// The fixed-frame configuration x corner id graph (447,120 x 960 entries).
    struct FixedAc<'a> {
        s: &'a Space,
        corner: &'a [[u16; N_MOVES]],
    }

    impl Graph for FixedAc<'_> {
        const MOVES: usize = N_MOVES;
        type Node = (usize, usize);
        type Cache = ();
        fn node(&self, i: usize, (): &mut ()) -> (usize, usize) {
            (i / N_CORNER, i % N_CORNER)
        }
        fn step(&self, &(cfg, c): &(usize, usize), m: usize) -> usize {
            let n = self.s.next[cfg][m];
            if n == NONE {
                cfg * N_CORNER + c
            } else {
                n as usize * N_CORNER + self.corner[c][m] as usize
            }
        }
    }

    /// The AC table under symmetry gives the distances of the fixed-frame configuration x
    /// corner table (both walked down, from positions spread over the phase-2 space).
    #[test]
    fn symmetric_ac_matches_the_fixed_frame_table() {
        let s = space();
        let fixed = bfs(
            &FixedAc {
                s,
                corner: &move_tables().c,
            },
            N_CFG * N_CORNER,
            &[0],
            "fixed ac",
            false,
        );
        let ac = build_ac(false);
        let ga = ac_graph();
        let t = ChenTables::ac_index_with(s, &solved_p2());
        let mut targets = vec![t];
        ga.equivalents(t, &mut targets);
        assert_eq!(fixed.hist.iter().sum::<u64>(), (N_CFG * N_CORNER) as u64);
        assert_eq!(fixed.hist.len(), ac.hist.len());
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(7);
        for _ in 0..500 {
            let mut p = solved_p2();
            for _ in 0..rng.gen_range(0..60) {
                p = p.apply(s, rng.gen_range(0..N_MOVES)).unwrap_or(p);
            }
            let d_fixed = walk(
                &fixed,
                &p,
                |q| q.cfg as usize * N_CORNER + q.c as usize,
                |i| i == 0,
                true,
            );
            let d_sym = walk(
                &ac,
                &p,
                |q| ChenTables::ac_index_with(s, q),
                |i| targets.contains(&i),
                true,
            );
            assert!(d_fixed.is_some());
            assert_eq!(d_fixed, d_sym, "{p:?}");
        }
    }
}
