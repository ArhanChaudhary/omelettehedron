//! Coset solver for `H_A` = <U, F, BR, BL>, after the method that proved the Rubik's Cube's
//! diameter is 20 (Rokicki, Kociemba, Davidson and Dethridge; their `hcoset` program).
//!
//! * [`space`]: the coset's elements indexed by states of `H_A`, their depth marks and the
//!   prepass.
//! * [`enumerate`]: every phase-1 sequence of a given length that ends in a class-B move.
//! * this module: the driver ([`solve`]), and the elements and solutions behind the marks.
//!
//! # Cosets
//!
//! A *coset* is the set `Q` of the 933,120 positions that one phase-1 sequence takes into `H_A`
//! together: `Q = T u` for a word `u`, where `T` is the set of positions of `H_A` (fixed-corner
//! frame, [`crate::cube`]) and `Q s = T` for every phase-1 sequence `s` of `Q` (all elements of
//! a coset share their phase-1 sequences). There are about 3.4 x 10^16 of them. The solver
//! takes any element `y` of `Q` and bounds the distance of every element from above, exactly
//! up to the length to which phase 1 is enumerated.
//!
//! # Indexing, and which side moves act on
//!
//! Positions are patterns (same-coloured centres are interchangeable) and moves act on them on
//! the right ([`Cube::apply`]). Corners and edges are distinct pieces, so their part of a
//! position (`x_CE`, slot -> piece) is a permutation, and a word `w` composes with its
//! permutation: `(x w)_CE = x_CE . w_CE`. In `H_A` the centres follow from the corners and
//! edges ([`Phase2::index`] checks it), so an element of `T` is its CE part, and these parts
//! form a group of 933,120 (corners and edges independent). So does every coset: within `Q` the
//! CE part determines the position.
//!
//! The element `z` of `Q` has the index `I(z)` = the element of `T` whose CE part is
//! `y_CE . z_CE^-1`. If the word `w` solves `z` (`z_CE = w_CE^-1`) then `I(z) = y w`, so:
//!
//! * a phase-1 sequence `s` of `n` moves marks the state `y s` of `T` it ends in, at depth `n`:
//!   the element solved by `s` alone (`z = s^-1` applied to solved);
//! * the prepass at depth `n` applies each class-A move to the states marked at `n - 1`: a
//!   solution `s p` (`p` class-A moves) of `z` marks `I(z) = y s p`.
//!
//! Both steps are the ordinary right action of moves on positions of `T`, as in `hcoset`, whose
//! bitmap is indexed by the phase-2 state reached by applying the phase-1 sequence to the coset
//! representative. A left action ("premoves" on the element itself) would not be well defined
//! here: two words that make the same pattern differ by a permutation of identical centres,
//! which a premove does not commute with.
//!
//! The element with index `p` is `z = q r` applied to solved, where `q` solves `p` in `H_A`
//! (`p q` = solved) and `r` makes `y` from solved: `(q r)_CE = p_CE^-1 . y_CE`. Here `r` is the
//! inverse of the first two-phase solution of `y` found (a phase-1 sequence and its phase-2
//! completion), so a pattern is enough as input. The tests check these sides directly
//! (`tests/coset.rs`): the element of each index, applied to every sequence found, ends in the
//! state the marks say.
//!
//! # Exactness
//!
//! Every optimal solution can be taken canonical ([`allowed_after`]: reordering commuting
//! same-axis turns does not change its length) and split into `s`, up to its last class-B move,
//! and `p`, class-A moves only. `s` ends in a class-B move and is found when phase 1 is
//! enumerated at its length. So if every length up to `n` has been enumerated, the elements
//! marked at depths up to `n` are exactly those at distance up to `n`, at their distances. Above
//! the longest enumerated length the marks (and their projection by the remaining prepasses)
//! are upper bounds.

pub mod enumerate;
pub mod space;

use crate::cube::*;
use crate::tables::*;
use enumerate::{Endpoints, Phase1, Stats};
use space::{Marks, Space, N_COSET, UNSET};
use std::sync::atomic::AtomicBool;
use std::time::Instant;

/// What to compute.
#[derive(Clone, Copy, Debug)]
pub struct Options {
    /// Fill the coset to this depth: every element is bounded by the end, or listed as not
    /// covered.
    pub target: usize,
    /// Longest phase-1 sequences to enumerate. Lengths beyond it get only the prepass.
    pub max_search: usize,
    /// Stop enumerating once every element is within `target` (the remaining depths then
    /// get only the prepass, as above `max_search`). Exact distances then reach only to the
    /// last complete length.
    pub stop_when_covered: bool,
    /// Keep the phase-1 sequences found (for solutions of chosen elements).
    pub keep_paths: bool,
    /// Groups the tasks of a length are searched in (more: the coverage is checked more
    /// often, so `stop_when_covered` ends a length sooner).
    pub chunks: usize,
    /// Mark only the sequences that end in `H_A` with the B-centres told apart (the default,
    /// and the only sound choice: see the module documentation). `false` marks every sequence
    /// that ends in `H_A` as a pattern, which is only for showing why not.
    pub labelled: bool,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            target: 23,
            max_search: 23,
            stop_when_covered: true,
            keep_paths: false,
            chunks: 1,
            labelled: true,
        }
    }
}

/// How far a length was enumerated.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Searched {
    /// every phase-1 sequence of this length was found
    Full,
    /// the enumeration stopped once every element was within the target
    Partial,
    /// not enumerated (below the phase-1 lower bound, or above the last length searched)
    No,
}

/// One depth of the fill.
#[derive(Clone, Debug)]
pub struct Level {
    pub depth: usize,
    pub searched: Searched,
    /// phase-1 sequences of this length found (ending in a class-B move)
    pub sequences: u64,
    /// elements first marked at this depth by the prepass and by the sequences
    pub from_prepass: u64,
    pub from_sequences: u64,
    /// elements marked at this depth or below
    pub covered: u64,
    /// the depth at which every element would be covered with no more sequences
    pub projected_max: Option<u8>,
    pub stats: Stats,
    pub seconds: f64,
}

/// A coset filled to a target depth.
pub struct Coset {
    /// the class-A moves on the index
    pub space: Space,
    /// depth of the first mark of each element (`UNSET`: none)
    pub marks: Marks,
    /// the final upper bound of each element (the marks extended by every later prepass)
    pub bounds: Vec<u8>,
    pub levels: Vec<Level>,
    /// the coset's phase-1 depth, if it was reached
    pub phase1_depth: Option<usize>,
    /// marks at depths up to this are exact distances (every length up to it enumerated)
    pub exact_through: Option<usize>,
    /// a word that makes the given position from solved (internal moves): the inverse of its
    /// first two-phase solution
    pub rep: Option<Vec<usize>>,
    /// the phase-1 sequences found, by length (if kept)
    pub sequences: Vec<Endpoints>,
    /// the sequences that end in `H_A` only as patterns (if kept; see the module
    /// documentation)
    pub others: Vec<Endpoints>,
}

/// Fill the coset of `cube` to `opt.target` on the current rayon pool, calling `progress` after
/// each depth.
///
/// # Errors
///
/// If the tables are of a kind the enumeration does not read (see [`Phase1::new`]).
pub fn solve<I: AeIndex>(
    t: &Tables<I>,
    cube: &Cube,
    opt: &Options,
    mut progress: impl FnMut(&Level),
) -> Result<Coset, String> {
    // The lift: the one that the first two-phase solution found solves with its B-centres told
    // apart (so the coset's phase-1 depth is the position's).
    let (lift, rep) = if opt.labelled {
        let p0 = Phase1::new(t, *cube, None)?;
        let (first, _) = p0.depth(enumerate::MAX_LEN - 1);
        let (_, s) = first.ok_or("no phase-1 sequence found")?;
        let q = t.p2.solve(&cube.apply_seq(&s)).ok_or("phase 2 failed")?;
        let sol: Vec<usize> = s.into_iter().chain(q).collect();
        let r = invert(&sol);
        let lab = Cube {
            xb: std::array::from_fn(|i| i as u8),
            ..Cube::solved()
        };
        (Some(lab.apply_seq(&r)), Some(r))
    } else {
        (None, None)
    };
    let p1 = Phase1::new(t, *cube, lift)?;
    let space = Space::new(&t.p2);
    let mut marks = Marks::new();
    let mut levels = vec![];
    let mut sequences = vec![];
    let mut others = vec![];
    let mut phase1_depth = None;
    let mut exact_through = None;
    let mut rep: Option<Vec<usize>> = rep;
    let lb = p1.lower_bound();
    let stop = AtomicBool::new(false);
    // set once every element is within the target: only prepasses after that
    let mut covered_at: Option<usize> = None;
    for n in 0..=opt.target {
        let start = Instant::now();
        let from_prepass = marks.prepass(&space, n);
        let mut level = Level {
            depth: n,
            searched: Searched::No,
            sequences: 0,
            from_prepass,
            from_sequences: 0,
            covered: 0,
            projected_max: None,
            stats: Stats::default(),
            seconds: 0.0,
        };
        let search = n >= lb && n <= opt.max_search && covered_at.is_none();
        if search {
            let mut found = Endpoints {
                len: n,
                ..Endpoints::default()
            };
            let mut stats = Stats::default();
            let mut new = 0u64;
            let mut done_early = false;
            let mut other = found.clone();
            let complete = p1.level(n, opt.keep_paths, opt.chunks, &stop, |e, o, s| {
                for &i in &e.index {
                    new += u64::from(marks.mark(i, n));
                }
                stats.add(&s);
                found.index.extend_from_slice(&e.index);
                found.paths.extend_from_slice(&e.paths);
                other.index.extend_from_slice(&o.index);
                other.paths.extend_from_slice(&o.paths);
                if opt.stop_when_covered && !e.index.is_empty() {
                    let max = marks.project(&space).into_iter().max().unwrap_or(UNSET);
                    if usize::from(max) <= opt.target {
                        done_early = true;
                        return false;
                    }
                }
                true
            });
            level.searched = if complete {
                exact_through = Some(n);
                Searched::Full
            } else {
                Searched::Partial
            };
            if done_early {
                covered_at = Some(n);
            }
            level.sequences = found.index.len() as u64;
            level.from_sequences = new;
            level.stats = stats;
            if phase1_depth.is_none() && !found.index.is_empty() {
                phase1_depth = Some(n);
            }
            if rep.is_none() {
                rep = representative(t, cube, &found, &p1);
            }
            if opt.keep_paths {
                sequences.push(found);
                others.push(other);
            }
        }
        level.covered = marks.count();
        if level.covered > 0 {
            let p = marks.project(&space);
            level.projected_max = p.iter().copied().max().filter(|&m| m != UNSET);
        }
        level.seconds = start.elapsed().as_secs_f64();
        progress(&level);
        levels.push(level);
        if marks.count() == N_COSET as u64 {
            break;
        }
    }
    let bounds = marks.project(&space);
    Ok(Coset {
        space,
        marks,
        bounds,
        levels,
        phase1_depth,
        exact_through,
        rep,
        sequences,
        others,
    })
}

/// The inverse of a two-phase solution of `cube` through one of the sequences `found`: a word
/// that makes `cube` from solved.
fn representative<I: AeIndex>(
    t: &Tables<I>,
    cube: &Cube,
    found: &Endpoints,
    p1: &Phase1<I>,
) -> Option<Vec<usize>> {
    // the sequences themselves may not be kept: find one again
    let s: Vec<usize> = if let Some(p) = found.path(0) {
        p.iter().map(|&m| usize::from(m)).collect()
    } else {
        if found.index.is_empty() {
            return None;
        }
        p1.first(found.len).0?
    };
    let end = cube.apply_seq(&s);
    let q = t.p2.solve(&end)?;
    let sol: Vec<usize> = s.into_iter().chain(q).collect();
    if cube.apply_seq(&sol) != Cube::solved() {
        return None;
    }
    // with the B-centres told apart the solution must solve the lift itself, so that the
    // inverse makes the lift (the coset's elements are then `q r` for this `r`)
    let labelled_solved: [u8; 12] = std::array::from_fn(|i| i as u8);
    if p1.labelled() && p1.lift().apply_seq(&sol).xb != labelled_solved {
        return None;
    }
    Some(invert(&sol))
}

/// The inverse of a sequence of internal moves.
pub fn invert(seq: &[usize]) -> Vec<usize> {
    seq.iter().rev().map(|&m| inverse_move(m)).collect()
}

impl Coset {
    /// A word (internal moves) that makes the element with index `p` from solved, and the
    /// element itself (see the module documentation: `q r`).
    pub fn element(&self, p2: &Phase2, p: u32) -> Option<(Vec<usize>, Cube)> {
        let rep = self.rep.as_ref()?;
        let q = self.space.solve(p2, p)?;
        let word: Vec<usize> = q.into_iter().chain(rep.iter().copied()).collect();
        let z = Cube::solved().apply_seq(&word);
        Some((word, z))
    }

    /// The shortest solution of element `z` through the kept phase-1 sequences: each sequence
    /// followed by its phase-2 completion (the per-element computation that the marks replace).
    pub fn solution(&self, p2: &Phase2, z: &Cube) -> Option<Vec<usize>> {
        let mut best: Option<(usize, Vec<usize>)> = None;
        for e in &self.sequences {
            for k in 0..e.index.len() {
                let Some(s) = e.path(k) else { continue };
                let end = s.iter().fold(*z, |c, &m| c.apply(usize::from(m)));
                let Some(d2) = p2.length(&end) else { continue };
                let total = s.len() + usize::from(d2);
                if best.as_ref().is_none_or(|b| total < b.0) {
                    let q = p2.solve(&end)?;
                    let sol = s.iter().map(|&m| usize::from(m)).chain(q).collect();
                    best = Some((total, sol));
                }
            }
        }
        best.map(|b| b.1)
    }

    /// The elements not within `depth`.
    pub fn uncovered(&self, depth: usize) -> Vec<u32> {
        (0..N_COSET as u32)
            .filter(|&i| usize::from(self.bounds[i as usize]) > depth)
            .collect()
    }
}
