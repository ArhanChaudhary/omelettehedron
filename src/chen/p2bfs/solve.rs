//! Positions of the search's slots, and their shortest phase-2 solutions.
//!
//! A slot of the full instance is a pattern-view position (the class's representative
//! configuration, the corner state, the B-centre colouring); [`PState`] holds it, makes its
//! cube (all six corners where they are: no corner held) and pattern, and turns it into the
//! fixed-corner frame of the phase-2 tables ([`super::super::tables`]). [`Solver`] finds a
//! shortest solution with the allowed moves by IDA*, bounded by the AC, PX and CC tables
//! (exact distances of projections of the restricted move graph, or of all moves: lower
//! bounds), in the pattern view, so the moves are the 16 face turns.

use super::super::space::{space, Space};
use super::super::tables::{ChenTables, P2};
use super::instance::{Instance, Kind};
use crate::coord::{mset, N_MSET};
use crate::cube::{Cube, N_MOVES};
use crate::fast::{fast, masks, shuffle, FastTables};
use crate::notation::{invert_moves, PhysMove};
use crate::pattern;
use crate::prune::decode;

/// A pattern-view position.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PState {
    /// pattern-view configuration
    pub q: u16,
    /// pattern-view corner state
    pub pc: u16,
    /// B-centre colouring (labels in lanes 0..12)
    pub xb: [u8; 16],
}

/// Conversions between the pattern view and the fixed-corner frame.
pub struct Views {
    /// per pattern-view corner state: the frame and the fixed-frame corner id
    frame_of: Vec<(u8, u16)>,
    /// per (pattern configuration, frame): the fixed-frame configuration
    cfg_of: Vec<u32>,
    /// per frame: the gather from the pattern view back to the fixed frame
    back: Vec<[u8; 16]>,
    /// the B-centre gathers of the internal moves in the pattern view
    gathers: [[u8; 16]; N_MOVES],
}

impl Views {
    /// Build the conversions.
    pub fn new() -> Views {
        let s = space();
        let n_frames = s.rot_xb.len();
        let mut frame_of = vec![(0u8, 0u16); s.corners.next.len()];
        for r in 0..n_frames {
            for c in 0..crate::coord::N_CORNER {
                let pc = s.corners.of_fixed[r * crate::coord::N_CORNER + c] as usize;
                frame_of[pc] = (r as u8, c as u16);
            }
        }
        let mut cfg_of = vec![u32::MAX; s.pat.cfg.len() * n_frames];
        for cfg in 0..s.cfgs.len() {
            let q = s.pat.of_cfg[cfg] as usize;
            cfg_of[q * n_frames + s.rot[cfg] as usize] = cfg as u32;
        }
        let back = s
            .rot_xb
            .iter()
            .map(|g| {
                let mut b = [0x80u8; 16];
                for (i, &j) in g.iter().enumerate().take(12) {
                    if let Some(slot) = b.get_mut(j as usize) {
                        *slot = i as u8;
                    }
                }
                b
            })
            .collect();
        let gathers = std::array::from_fn(|m| {
            let mut g = [0x80u8; 16];
            g[..12].copy_from_slice(&s.pat.xb_src[m]);
            g
        });
        Views {
            frame_of,
            cfg_of,
            back,
            gathers,
        }
    }

    /// The fixed-frame phase-2 position of a pattern-view position.
    pub fn to_p2(&self, st: &PState) -> P2 {
        let s = space();
        let (r, c) = self.frame_of[st.pc as usize];
        let cfg = self.cfg_of[st.q as usize * s.rot_xb.len() + r as usize];
        let mut xb = shuffle(&st.xb, &self.back[r as usize]);
        for v in &mut xb[12..] {
            *v = 0;
        }
        P2 { cfg, c, xb }
    }

    /// The fixed-frame corner id and B-centre colouring of a pattern-view position's corners
    /// and B-centres (for the CC table).
    #[inline]
    fn cb_fixed(&self, pc: u16, xb: &[u8; 16]) -> (u16, [u8; 16]) {
        let (r, c) = self.frame_of[pc as usize];
        let mut x = shuffle(xb, &self.back[r as usize]);
        for v in &mut x[12..] {
            *v = 0;
        }
        (c, x)
    }

    /// The position after internal move `m` (the pattern-view face turn
    /// [`super::super::space::Patterns::pmove`]`[m]`), if allowed.
    #[inline]
    pub fn apply(&self, st: &PState, m: usize) -> Option<PState> {
        let s = space();
        let q = s.pat.next[st.q as usize][m];
        if q == u16::MAX {
            return None;
        }
        let mut xb = shuffle(&st.xb, &self.gathers[m]);
        for v in &mut xb[12..] {
            *v = 0;
        }
        Some(PState {
            q,
            pc: s.corners.next[st.pc as usize][m],
            xb,
        })
    }
}

impl PState {
    /// The solved position.
    pub fn solved() -> PState {
        let mut xb = [0u8; 16];
        xb[..12].copy_from_slice(&Cube::solved().xb);
        PState {
            q: 0,
            pc: space().corners.solved,
            xb,
        }
    }

    /// The position of slot `(c, x, p)` of the full instance (the class's representative).
    pub fn of_slot(inst: &Instance, c: usize, x: usize, p: usize) -> Option<PState> {
        if inst.kind != Kind::FULL || x >= N_MSET || c >= inst.n_class {
            return None;
        }
        let mut xb = [0u8; 16];
        xb[..12].copy_from_slice(&mset().unrank(x as u32));
        Some(PState {
            q: inst.rep[c],
            pc: *inst.p_to_pc.get(p)?,
            xb,
        })
    }

    /// Its cube: the pieces where they are in the pattern view.
    pub fn cube(&self) -> Cube {
        let s = space();
        let (xa, ep) = s.cfgs[s.pat.cfg[self.q as usize] as usize];
        let (cp, co) = s.corners.states[self.pc as usize];
        let mut xb = [0u8; 12];
        xb.copy_from_slice(&self.xb[..12]);
        Cube { cp, co, ep, xa, xb }
    }

    /// Its pattern (72 stickers, standard colours, as `omelettehedron show` prints).
    pub fn pattern(&self) -> String {
        pattern::format(&pattern::of_cube(&self.cube()))
    }
}

/// The face turn of internal move `m`.
pub fn phys_move(m: usize) -> PhysMove {
    let pm = space().pat.pmove[m] as usize;
    PhysMove::face(pm / 2, pm.is_multiple_of(2))
}

/// IDA* with the allowed moves, bounded by the AC, PX and CC tables.
pub struct Solver<'a> {
    s: &'a Space,
    f: &'a FastTables,
    t: &'a ChenTables,
    v: Views,
    /// the face of each internal move
    face: [usize; N_MOVES],
}

impl<'a> Solver<'a> {
    /// A solver with loaded tables.
    pub fn new(t: &'a ChenTables) -> Solver<'a> {
        let s = space();
        let pm = s.pat.pmove;
        Solver {
            s,
            f: fast(),
            t,
            v: Views::new(),
            face: std::array::from_fn(|m| pm[m] as usize / 2),
        }
    }

    /// The conversions.
    pub fn views(&self) -> &Views {
        &self.v
    }

    /// AC, PX and CC indices of a position.
    #[inline]
    fn indices(&self, st: &PState) -> [usize; 3] {
        let pat = &self.s.pat;
        let (cls, k) = pat.class[st.q as usize];
        let k = k as usize;
        let ac = cls as usize * super::super::corners::N_PCORNER
            + self.s.corners.sym[st.pc as usize][k] as usize;
        let px = cls as usize * N_MSET + self.f.rank(&masks(&pat.sym_xb(k, &st.xb))) as usize;
        let (c, xb) = self.v.cb_fixed(st.pc, &st.xb);
        let cc = self.t.cc_ix.index_b(c, &xb);
        [ac, px, cc]
    }

    /// The exact table distances (AC, PX, CC) of a position, `None` if a table does not walk
    /// down.
    pub fn exact(&self, st: &PState) -> Option<[u8; 3]> {
        let e = self.t.exact(&self.v.to_p2(st))?;
        Some([e[0], e[2], e[3]])
    }

    /// Whether the table indices of the pattern view agree with the fixed frame's.
    pub fn indices_agree(&self, st: &PState) -> bool {
        let p2 = self.v.to_p2(st);
        let [ac, px, cc] = self.indices(st);
        ac == ChenTables::ac_index_with(self.s, &p2)
            && px == ChenTables::px_index_with(self.s, self.f, &p2)
            && cc == self.t.cc_index(&p2)
    }

    /// A shortest solution (internal moves) with the allowed moves, at most `max` moves.
    pub fn solve(&self, st: &PState, max: usize) -> Option<Vec<usize>> {
        let h = self.exact(st)?;
        let solved = PState::solved();
        let mut bound = *h.iter().max().unwrap_or(&0) as usize;
        let mut path = Vec::new();
        while bound <= max {
            if self.dfs(st, &solved, h, bound, N_MOVES, &mut path) {
                return Some(path);
            }
            bound += 1;
        }
        None
    }

    /// [`Solver::solve`] on all threads: each bound's tree split at depth 3 into tasks.
    pub fn solve_par(&self, st: &PState, max: usize) -> Option<Vec<usize>> {
        use rayon::prelude::*;
        use std::sync::atomic::{AtomicBool, Ordering};
        let h = self.exact(st)?;
        let solved = PState::solved();
        // the nodes at depth 0..=3 that the tables do not cut, with their paths
        let mut bound = *h.iter().max().unwrap_or(&0) as usize;
        while bound <= max {
            let mut level: Vec<(PState, [u8; 3], Vec<usize>)> = vec![(*st, h, vec![])];
            for _ in 0..3 {
                let mut next = vec![];
                for (n, hn, path) in &level {
                    if n == &solved {
                        return Some(path.clone());
                    }
                    let g = path.len();
                    if g >= bound {
                        continue;
                    }
                    let left = (bound - g - 1) as u8;
                    for m in 0..N_MOVES {
                        if path.last().is_some_and(|&l| self.face[m] == self.face[l]) {
                            continue;
                        }
                        let Some(c) = self.v.apply(n, m) else {
                            continue;
                        };
                        let Some(hc) = self.child_h(&c, *hn, left) else {
                            continue;
                        };
                        let mut p = path.clone();
                        p.push(m);
                        next.push((c, hc, p));
                    }
                }
                level = next;
            }
            let found = AtomicBool::new(false);
            let sol = level.par_iter().find_map_any(|(n, hn, path)| {
                if found.load(Ordering::Relaxed) {
                    return None;
                }
                let mut p = path.clone();
                let last = p.last().copied().unwrap_or(N_MOVES);
                if self.dfs(n, &solved, *hn, bound, last, &mut p) {
                    found.store(true, Ordering::Relaxed);
                    Some(p)
                } else {
                    None
                }
            });
            if sol.is_some() {
                return sol;
            }
            bound += 1;
        }
        None
    }

    /// The exact table distances of a child (from its parent's `h`), `None` if one is above
    /// `left`.
    #[inline]
    fn child_h(&self, n: &PState, h: [u8; 3], left: u8) -> Option<[u8; 3]> {
        let i = self.indices(n);
        let ac = decode(h[0], self.t.ac.get(i[0]));
        if ac > left {
            return None;
        }
        let px = decode(h[1], self.t.px.as_ref()?.get(i[1]));
        if px > left {
            return None;
        }
        let cc = decode(h[2], self.t.cc.get(i[2]));
        (cc <= left).then_some([ac, px, cc])
    }

    fn dfs(
        &self,
        st: &PState,
        solved: &PState,
        h: [u8; 3],
        bound: usize,
        last: usize,
        path: &mut Vec<usize>,
    ) -> bool {
        if st == solved {
            return true;
        }
        let g = path.len();
        if g >= bound {
            return false;
        }
        let left = (bound - g - 1) as u8;
        for m in 0..N_MOVES {
            // (a face turned twice in a row is one turn of that face)
            if last < N_MOVES && self.face[m] == self.face[last] {
                continue;
            }
            let Some(n) = self.v.apply(st, m) else {
                continue;
            };
            let Some(hn) = self.child_h(&n, h, left) else {
                continue;
            };
            path.push(m);
            if self.dfs(&n, solved, hn, bound, m, path) {
                return true;
            }
            path.pop();
        }
        false
    }
}

/// A position's solution and scramble in face turns, and whether the scramble's pattern is the
/// position's pattern (a check of the conversions).
pub struct Written {
    /// the solution
    pub solution: Vec<PhysMove>,
    /// its inverse: a scramble that gives the position
    pub scramble: Vec<PhysMove>,
    /// the pattern
    pub pattern: String,
    /// whether the scramble's pattern equals the position's
    pub consistent: bool,
}

/// Write a solution (internal moves) of a position in face turns.
pub fn write(st: &PState, sol: &[usize]) -> Written {
    let solution: Vec<PhysMove> = sol.iter().map(|&m| phys_move(m)).collect();
    let scramble = invert_moves(&solution);
    let pat = st.pattern();
    let from_scramble = pattern::format(&pattern::of_moves(&scramble));
    Written {
        consistent: from_scramble == pat,
        solution,
        scramble,
        pattern: pat,
    }
}
