//! The relaxed phase-1 target: any configuration of the phase-2 space, not only `H'`.
//!
//! Phase 1 of the block split ends in `H'`, the 972 configurations with every hexagon intact.
//! Phase 2 can start from any of the 447,120 configurations ([`super::space`]), so phase 1 may
//! stop as soon as the A-centres and edges form one of them (a shape of colour groups that the
//! allowed moves reach, with its edges in one of their 81 block turns): the target of 460 x 3^4
//! configurations per frame suggested by Chen Shuang. Phase-1 paths get shorter; phase 2 starts
//! further from solved.
//!
//! Whether a position's A-centres and edges form a configuration depends only on its
//! A-centre/edge (AE) coordinate (the 81 block turns of a shape share one), so the target is a
//! set of AE indices, and the table is [`crate::tables::build_ae_with`]'s breadth-first search
//! from all of them, on the tier's AE index:
//!
//! * with the edge chirality (medium tier and up) the set is closed under the index's 288
//!   symmetries (colour relabelling and conjugation keep the move rule and map the frames onto
//!   each other), so the table holds exact distances;
//! * without it (small tier) an index also holds the other chirality, which is no
//!   configuration, so the table is a lower bound.
//!
//! Either way the search checks a leaf in full ([`super::space::Space::id`]). The tiny tier's
//! index has no table of distances to build from, so the relaxed target needs the small tier
//! or larger.
//!
//! The 447,120 configurations fall into 27 target indices of the 288-fold indices. Phase-1
//! distances (over the table entries): medium max 11, mean 8.45 (`H'`: 14, 10.94); small, a
//! lower bound, max 10, mean 7.61 (`H'`: 12, 10.10). Generation: 4.7 s at 48 threads and
//! 147 s on one (medium, 1.03 GB), 0.7 s and 19.7 s (small, 129 MB).
//!
//! Measured on 100 random positions (seed 1, 10 s each; mean best length, paired difference
//! with a 95% interval): against the target `H'`, shorter at every checkpoint, by 0.70 +- 0.22
//! (small) and 0.58 +- 0.23 (medium) moves at 10 s on one thread, by 0.6-1.3 moves at 0.1 s
//! on 16 and 48 threads and by 0.02-0.16 (within the intervals) at 10 s there; first solutions
//! about 1 move shorter and 0.70-0.78x the time. Against the two-phase solver at 10 s: on par on
//! the small tier on one thread (-0.17 +- 0.23), behind otherwise (+0.21 to +0.45). It reads the
//! tier's own AE table too, so the medium tier needs 2.5 GB instead of 1.4 GB.

use super::space::space;
use super::tables::ChenTableError;
use crate::coord::Coords;
use crate::cube::{Cube, MoveSet, N_MOVES};
use crate::fast::{fast, masks, Fc};
use crate::prune::Table2;
use crate::tables::{build_ae_with, load_or_build, AeIndex};
use std::path::Path;

/// The relaxed phase-1 table and its targets.
pub struct RelaxedAe {
    /// distances (mod 3) to the nearest target, on the tier's AE index
    pub tab: Table2,
    /// the target indices, sorted
    pub targets: Vec<usize>,
}

/// The AE indices of every configuration (and the indices equivalent to them), sorted.
pub fn relaxed_targets<I: AeIndex>(ix: &I) -> Vec<usize> {
    let s0 = Cube::solved();
    let mut v: Vec<usize> = space()
        .cfgs
        .iter()
        .map(|&(xa, ep)| ix.index(&Coords::of(&Cube { ep, xa, ..s0 })))
        .collect();
    v.sort_unstable();
    v.dedup();
    let mut eq = vec![];
    for &t in &v {
        ix.equivalents(t, &mut eq);
    }
    v.extend(eq);
    v.sort_unstable();
    v.dedup();
    v
}

/// The file name of the relaxed table of an index.
pub fn file_name<I: AeIndex>(ix: &I) -> String {
    format!("{}-relaxed.prun", ix.name())
}

impl RelaxedAe {
    /// Build the table in memory.
    pub fn generate<I: AeIndex>(ix: &I, verbose: bool) -> RelaxedAe {
        let targets = relaxed_targets(ix);
        let name = format!("{}-relaxed", ix.name());
        let tab = build_ae_with(ix, MoveSet::FACE, &targets, &name, verbose);
        RelaxedAe { tab, targets }
    }

    /// Load the table of `ix` from `dir`, generating (and saving) it if missing.
    ///
    /// # Errors
    ///
    /// [`ChenTableError::NoRelaxedTable`] for the tiny tier's index; [`ChenTableError::Io`] if
    /// the file cannot be read or written, or holds a table of the wrong size or checksum;
    /// [`ChenTableError::WrongLayout`] if it is not a distance table of this layout.
    pub fn load<I: AeIndex>(
        dir: &Path,
        ix: &I,
        verbose: bool,
    ) -> Result<RelaxedAe, ChenTableError> {
        if I::DIRECT {
            return Err(ChenTableError::NoRelaxedTable {
                index: ix.name().to_string(),
            });
        }
        let targets = relaxed_targets(ix);
        let path = dir.join(file_name(ix));
        let name = format!("{}-relaxed", ix.name());
        let tab = load_or_build(&path, ix.size(), verbose, || {
            build_ae_with(ix, MoveSet::FACE, &targets, &name, verbose)
        })
        .map_err(|source| ChenTableError::Io {
            path: path.display().to_string(),
            source,
        })?;
        let t = RelaxedAe { tab, targets };
        // a file of another layout passes its checksum but is no distance table
        let mut x = fast().of(&Cube::solved());
        let f = fast();
        for m in [0, 5, 9, 2, 12, 7, 1, 14] {
            x = f.apply_m(&x, &masks(&x.e), m);
        }
        if t.exact(ix, &x).is_none() {
            return Err(ChenTableError::WrongLayout {
                dir: dir.display().to_string(),
            });
        }
        Ok(t)
    }

    /// Whether index `i` is a target.
    #[inline]
    pub fn is_target(&self, i: usize) -> bool {
        self.targets.binary_search(&i).is_ok()
    }

    /// The table's distance of a position, by walking down the table (each step to a
    /// neighbour one closer) until a target. `None` if the walk does not reach one within the
    /// table's largest distance (a table of another layout).
    pub fn exact<I: AeIndex>(&self, ix: &I, x: &Fc) -> Option<u8> {
        let f = fast();
        let index = |y: &Fc| ix.index_fast(f, y, f.a_rank(y));
        let mut y = *x;
        for d in 0..self.tab.hist.len() {
            let i = index(&y);
            if self.is_target(i) {
                return Some(d as u8);
            }
            let want = (self.tab.get(i) + 2) % 3;
            let em = masks(&y.e);
            y = (0..N_MOVES)
                .map(|m| f.apply_m(&y, &em, m))
                .find(|z| self.tab.get(index(z)) == want)?;
        }
        None
    }
}

/// How many target indices have a configuration as their decoded representative, of how many
/// (all of them when the target set is closed under the index's symmetries).
pub fn target_check<I: AeIndex>(ix: &I, targets: &[usize]) -> (usize, usize) {
    let s0 = Cube::solved();
    let mut coords: Vec<(u32, u32, u8)> = space()
        .cfgs
        .iter()
        .map(|&(xa, ep)| {
            let c = Coords::of(&Cube { ep, xa, ..s0 });
            (c.a, c.e, c.ebits)
        })
        .collect();
    coords.sort_unstable();
    coords.dedup();
    let hits = targets
        .iter()
        .filter(|&&i| coords.binary_search(&ix.decode(i)).is_ok())
        .count();
    (hits, targets.len())
}
