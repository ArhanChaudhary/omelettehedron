//! Corner states in the pattern view, for the AC table reduced by symmetry.
//!
//! In the pattern view (see [`super::space`]) all six corners move: 360 even permutations x 32
//! orientations = [`N_PCORNER`] states. A fixed-frame corner state seen through one of the 12
//! frames is one of them, and each of them is one fixed-frame state seen through one frame.
//! The 24 symmetries that keep the face classes act on them by conjugation, as on the
//! pattern-view configurations, so the AC table (configuration x corners) is stored once per
//! symmetry class of configurations: 1,611 x 11,520 entries instead of 447,120 x 960.
//!
//! Every map here (a face turn, a whole-puzzle turn, a conjugation) moves each corner as a
//! whole, so it is tabulated per corner place, piece and orientation from the sticker model.

use super::space::{pattern_moves, QuickMap, N_PMOVES};
use crate::coord::N_CORNER;
use crate::cube::*;
use crate::geometry::{class_preserving_group, face_of_normal, mat_vec, Mat, FACES};

/// Corner states in the pattern view (checked when built).
pub const N_PCORNER: usize = 11_520;
/// Symmetries that keep the face classes (in the order of [`class_preserving_group`]).
pub const N_SYM: usize = 24;

/// A map of the corners: where the corner at each place goes, and what piece and orientation
/// it becomes there: `[place][piece][orientation] = (place, piece, orientation)`.
type CornerMap = [[[(u8, u8, u8); 2]; 6]; 6];

/// The corner map of a sticker map: sticker `i` goes to `perm[i]` with its colour (a face)
/// mapped by `colour`.
fn corner_map(perm: &[usize], colour: impl Fn(u8) -> u8) -> CornerMap {
    let d = defs();
    // the corner place and slot of each corner sticker
    let mut at = [(0usize, 0usize); 72];
    for v in 0..6 {
        for j in 0..4 {
            at[d.corner_slots[v][j]] = (v, j);
        }
    }
    let mut out = [[[(0u8, 0u8, 0u8); 2]; 6]; 6];
    for v in 0..6 {
        for k in 0..6 {
            for o in 0..2 {
                // the piece's four stickers after the map, by slot of their new place
                let mut cols = [0u8; 4];
                let mut to = 0;
                for j in 0..4 {
                    let (v2, j2) = at[perm[d.corner_slots[v][j]] % 72];
                    to = v2;
                    cols[j2] = colour(d.corner_cols[k][(j + 2 * o) % 4]);
                }
                let found = (0..6).find_map(|k2| {
                    (0..2)
                        .find(|&o2| (0..4).all(|j| d.corner_cols[k2][(j + 2 * o2) % 4] == cols[j]))
                        .map(|o2| (k2 as u8, o2 as u8))
                });
                // (a reachable map always gives a corner piece in one of its orientations)
                debug_assert!(found.is_some(), "corner map");
                let (k2, o2) = found.unwrap_or((0, 0));
                out[v][k][o] = (to as u8, k2, o2);
            }
        }
    }
    out
}

/// A corner state (piece and orientation at each place) after a corner map.
fn apply(map: &CornerMap, s: &([u8; 6], [u8; 6])) -> ([u8; 6], [u8; 6]) {
    let mut cp = [0u8; 6];
    let mut co = [0u8; 6];
    for v in 0..6 {
        let (v2, k2, o2) = map[v][s.0[v] as usize % 6][s.1[v] as usize & 1];
        cp[v2 as usize % 6] = k2;
        co[v2 as usize % 6] = o2;
    }
    (cp, co)
}

/// The pattern-view corner states and their maps.
pub struct PCorners {
    /// the pattern-view state of each fixed-frame corner id seen through each frame
    /// (index `frame * 960 + corner id`)
    pub of_fixed: Vec<u16>,
    /// the state after each internal move made in frame 0 (the pattern-view move
    /// [`super::space::Patterns::pmove`])
    pub next: Vec<[u16; N_MOVES]>,
    /// the state conjugated by each symmetry
    pub sym: Vec<[u16; N_SYM]>,
    /// the solved state
    pub solved: u16,
}

impl PCorners {
    /// Build the tables; `rots` are the frames' rotations and `pmove` the pattern-view move of
    /// each internal move (as in [`super::space::Space`]).
    pub(crate) fn new(rots: &[Mat], pmove: &[u8; N_MOVES]) -> PCorners {
        let d = defs();
        let sm = &d.sm;
        let moves: Vec<CornerMap> = pattern_moves()
            .iter()
            .map(|perm| corner_map(perm, |c| c))
            .collect();
        let syms: Vec<CornerMap> = class_preserving_group()
            .iter()
            .map(|m| {
                corner_map(&sm.whole(m), |c| {
                    face_of_normal(mat_vec(m, FACES[c as usize])) as u8
                })
            })
            .collect();
        let turns: Vec<CornerMap> = rots
            .iter()
            .map(|m| corner_map(&sm.whole(m), |c| c))
            .collect();
        // the states, breadth-first from solved over the pattern-view moves
        let s0 = Cube::solved();
        let mut states = vec![(s0.cp, s0.co)];
        let mut ids: QuickMap<([u8; 6], [u8; 6]), u16> = QuickMap::default();
        ids.insert((s0.cp, s0.co), 0);
        let mut pnext: Vec<[u16; N_PMOVES]> = vec![];
        let mut i = 0;
        while i < states.len() {
            let s = states[i];
            i += 1;
            let row = std::array::from_fn(|pm| {
                let n = apply(&moves[pm], &s);
                *ids.entry(n).or_insert_with(|| {
                    states.push(n);
                    (states.len() - 1) as u16
                })
            });
            pnext.push(row);
        }
        let id = |s: &([u8; 6], [u8; 6])| ids.get(s).copied().unwrap_or(0);
        let fixed = &crate::coord::corners().id_to_state;
        let of_fixed = turns
            .iter()
            .flat_map(|t| fixed.iter().map(move |s| (t, s)))
            .map(|(t, s)| id(&apply(t, s)))
            .collect();
        let next = pnext
            .iter()
            .map(|row| std::array::from_fn(|m| row[pmove[m] as usize % N_PMOVES]))
            .collect();
        let sym = states
            .iter()
            .map(|s| std::array::from_fn(|k| id(&apply(&syms[k], s))))
            .collect();
        debug_assert_eq!(states.len(), N_PCORNER);
        debug_assert_eq!(fixed.len(), N_CORNER);
        PCorners {
            of_fixed,
            next,
            sym,
            solved: 0,
        }
    }
}
