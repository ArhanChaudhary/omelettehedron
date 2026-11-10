//! Cubie-level representation in the fixed-corner frame.
//!
//! Face turns can rotate the whole puzzle, so the front corner (vertex 0, `+x`, shared by
//! U, R, L and F) is held fixed. A turn of a front face is replaced by the equivalent
//! two-layer turn of the opposite back face followed by a whole-puzzle rotation, giving 16
//! internal moves: single turns of B, BR, BL, D and wide turns of B, BR, BL, D. With slice
//! moves (see [`MoveSet`]) there are 8 more: the middle layer of each of these four axes. A
//! slice holds no corner, so it never moves the front corner either.
//!
//! Internal move index `m = 2 * slot + dir`, `dir` 0 = clockwise, 1 = anticlockwise (seen from
//! the back face). Slots 0..4 are single turns of faces 4..8 (B, BR, BL, D); slots 4..8 are
//! wide turns of the same faces (equivalent to F, L, R, U plus a rotation); slots 8..12 are
//! slice turns of the same axes, turned like the back face.

use crate::geometry::*;
use std::sync::OnceLock;

/// Internal moves of the face-turn metric (moves `0..N_MOVES`).
pub const N_MOVES: usize = 16;
/// Internal moves with the slice moves (`N_MOVES..MAX_MOVES` are the slices).
pub const MAX_MOVES: usize = 24;

#[inline]
pub fn slot(m: usize) -> usize {
    m >> 1
}
#[inline]
pub fn slot_face(s: usize) -> usize {
    4 + (s & 3)
}
#[inline]
pub fn slot_is_wide(s: usize) -> bool {
    s >> 2 == 1
}
#[inline]
pub fn slot_is_slice(s: usize) -> bool {
    s >> 2 == 2
}
#[inline]
pub fn move_cw(m: usize) -> bool {
    m & 1 == 0
}
#[inline]
pub fn inverse_move(m: usize) -> usize {
    m ^ 1
}

/// Class-A internal moves: BR, BL, Bw (= F), Dw (= U).
pub const A_SLOTS: [usize; 4] = [1, 2, 4, 7];

/// Whether `m` is a class-A move (a generator of `H_A`; slices are not).
#[inline]
pub fn is_a_move(m: usize) -> bool {
    matches!(slot(m), 1 | 2 | 4 | 7)
}

/// Canonical move order (`last >= MAX_MOVES`: no previous move). Never the same slot twice in
/// a row. The three layers of an axis commute, and modulo the free rotation about the axis
/// they generate 9 elements; of the same-axis pairs only one representative of each element
/// that needs two moves is allowed:
/// * face metric (`slice == false`): single before wide, either direction (the four
///   products are distinct and none is a single move);
/// * with slices: single then wide of the same direction only (X Xw, X' Xw'). The products
///   X Xw' and X' Xw are slices, and every pair with a slice is a single move or equals one
///   of the two allowed pairs.
///
/// (Checked against the axis group by `canonical_pairs_match_axis_group`, and against
/// whole-puzzle breadth-first search by `canonical_sequences_are_distinct_positions`.)
pub const fn allowed_after(slice: bool, last: usize, m: usize) -> bool {
    if last >= MAX_MOVES {
        return true;
    }
    let (ls, s) = (last >> 1, m >> 1);
    if s == ls {
        return false;
    }
    if s & 3 != ls & 3 {
        return true;
    }
    if slice {
        ls < 4 && s >> 2 == 1 && last & 1 == m & 1
    } else {
        s > ls
    }
}

/// The moves a solution may use. Face turns are always allowed; every allowed move counts as
/// one move.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct MoveSet {
    /// Wide (two-layer) turns. A wide turn is a turn of the opposite face plus a whole-puzzle
    /// rotation, and any orientation counts as solved, so they never shorten a solution: the
    /// search is unchanged and only the output differs (the front corner stays still, so
    /// front-face turns are printed as wide turns of the opposite back face).
    pub wide: bool,
    /// Slice (middle-layer) turns: 8 more internal moves, with their own pruning tables.
    pub slice: bool,
}

impl MoveSet {
    /// Face turns only.
    pub const FACE: MoveSet = MoveSet {
        wide: false,
        slice: false,
    };

    /// Number of internal moves (`0..n_moves()`).
    #[inline]
    pub const fn n_moves(self) -> usize {
        if self.slice {
            MAX_MOVES
        } else {
            N_MOVES
        }
    }

    /// [`allowed_after`] for this move set.
    #[inline]
    pub const fn allowed_after(self, last: usize, m: usize) -> bool {
        allowed_after(self.slice, last, m)
    }

    /// Suffix of the pruning-table file names: the tables hold distances in the move metric,
    /// so slices need their own (wide turns do not change the metric).
    pub fn table_suffix(self) -> &'static str {
        if self.slice {
            "-slice"
        } else {
            ""
        }
    }

    /// Short name, e.g. `face+slice`.
    pub fn name(self) -> String {
        let mut s = String::from("face");
        if self.wide {
            s.push_str("+wide");
        }
        if self.slice {
            s.push_str("+slice");
        }
        s
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub struct Cube {
    /// corner piece at each vertex position (position 0 is fixed)
    pub cp: [u8; 6],
    /// corner orientation (0/1) at each position
    pub co: [u8; 6],
    /// edge piece at each edge position
    pub ep: [u8; 12],
    /// A-colour index (0..4) of the centre at each A-centre position
    pub xa: [u8; 12],
    /// B-colour index of the centre at each B-centre position
    pub xb: [u8; 12],
}

pub struct Defs {
    pub sm: StickerModel,
    /// sticker indices of the 4 corner slots of each vertex, anticlockwise from the reference A slot
    pub corner_slots: [[usize; 4]; 6],
    /// home colours of each corner piece in slot order
    pub corner_cols: [[u8; 4]; 6],
    /// (A-slot, B-slot) sticker indices of each edge position
    pub edge_slots: [[usize; 2]; 12],
    /// (A face, B face) of each edge piece
    pub edge_cols: [[u8; 2]; 12],
    /// A-colour index (0..4) of the A face of each edge piece (its "A-triple")
    pub edge_atri: [u8; 12],
    /// B-colour index of the B face of each edge piece
    pub edge_btri: [u8; 12],
    pub a_pos: [usize; 12],
    pub b_pos: [usize; 12],
    pub solved: Cube,
    /// move tables as "source" arrays: new[p] = old[src[p]] (corner orientation added)
    pub moves: [Cube; MAX_MOVES],
    /// sticker permutations of the internal moves
    pub move_perms: Vec<Vec<usize>>,
}

pub fn defs() -> &'static Defs {
    static D: OnceLock<Defs> = OnceLock::new();
    D.get_or_init(build_defs)
}

fn rot90(u: V3, w: V3) -> V3 {
    // anticlockwise (right-handed) quarter turn about unit axis u
    let d = dot(u, w);
    let c = cross(u, w);
    [u[0] * d + c[0], u[1] * d + c[1], u[2] * d + c[2]]
}

fn build_defs() -> Defs {
    let sm = StickerModel::new();
    let st = &sm.st;
    // corners
    let mut corner_slots = [[0usize; 4]; 6];
    let mut corner_cols = [[0u8; 4]; 6];
    for (v, &vv) in VERTS.iter().enumerate() {
        // corner stickers at v (vv is a unit axis vector)
        let at_v: Vec<usize> = (0..72)
            .filter(|&i| st[i].kind == Kind::Corner && dot(st[i].p, vv) == 7)
            .collect();
        debug_assert_eq!(at_v.len(), 4);
        // the reference slot is on the lowest-numbered A face (every vertex touches two)
        let r0 = at_v
            .iter()
            .copied()
            .filter(|&i| is_a_face(st[i].face))
            .min_by_key(|&i| st[i].face);
        debug_assert!(r0.is_some());
        let mut cur = st[r0.unwrap_or_default()].p;
        for j in 0..4 {
            corner_slots[v][j] = sm.index_of(cur);
            corner_cols[v][j] = st[corner_slots[v][j]].face as u8;
            cur = rot90(vv, cur);
        }
    }
    // edges
    let edir = |p: V3| -> V3 {
        [
            if p[0].abs() == 4 { p[0].signum() } else { 0 },
            if p[1].abs() == 4 { p[1].signum() } else { 0 },
            if p[2].abs() == 4 { p[2].signum() } else { 0 },
        ]
    };
    let mut dirs: Vec<V3> = (0..72)
        .filter(|&i| st[i].kind == Kind::Edge)
        .map(|i| edir(st[i].p))
        .collect();
    dirs.sort_unstable();
    dirs.dedup();
    debug_assert_eq!(dirs.len(), 12);
    let mut edge_slots = [[0usize; 2]; 12];
    let mut edge_cols = [[0u8; 2]; 12];
    let mut edge_atri = [0u8; 12];
    let mut edge_btri = [0u8; 12];
    for (e, d) in dirs.iter().enumerate() {
        for i in 0..72 {
            if st[i].kind == Kind::Edge && edir(st[i].p) == *d {
                let k = usize::from(!is_a_face(st[i].face));
                edge_slots[e][k] = i;
                edge_cols[e][k] = st[i].face as u8;
            }
        }
        edge_atri[e] = class_index(edge_cols[e][0] as usize) as u8;
        edge_btri[e] = class_index(edge_cols[e][1] as usize) as u8;
    }
    // centre stickers of the A (B) faces in sticker order: stickers 6..9 of each face
    let a_pos: [usize; 12] = std::array::from_fn(|i| 9 * A_FACES[i / 3] + 6 + i % 3);
    let b_pos: [usize; 12] = std::array::from_fn(|i| 9 * B_FACES[i / 3] + 6 + i % 3);
    debug_assert!(a_pos
        .iter()
        .chain(&b_pos)
        .all(|&i| st[i].kind == Kind::Centre));
    let mut d = Defs {
        sm,
        corner_slots,
        corner_cols,
        edge_slots,
        edge_cols,
        edge_atri,
        edge_btri,
        a_pos,
        b_pos,
        solved: Cube {
            cp: [0; 6],
            co: [0; 6],
            ep: [0; 12],
            xa: [0; 12],
            xb: [0; 12],
        },
        moves: [Cube {
            cp: [0; 6],
            co: [0; 6],
            ep: [0; 12],
            xa: [0; 12],
            xb: [0; 12],
        }; MAX_MOVES],
        move_perms: vec![],
    };
    // every piece at home; three centres per face, faces in class order
    d.solved = Cube {
        cp: std::array::from_fn(|v| v as u8),
        co: [0; 6],
        ep: std::array::from_fn(|e| e as u8),
        xa: std::array::from_fn(|i| (i / 3) as u8),
        xb: std::array::from_fn(|i| (i / 3) as u8),
    };
    let solved_st = solved_stickers(&d);
    debug_assert_eq!(stickers_to_cube_with(&d, &solved_st), Some(d.solved));
    for m in 0..MAX_MOVES {
        let s = slot(m);
        let f = slot_face(s);
        let perm = if slot_is_slice(s) {
            d.sm.slice_turn(f, move_cw(m))
        } else if slot_is_wide(s) {
            d.sm.wide_turn(f, move_cw(m))
        } else {
            d.sm.face_turn(f, move_cw(m))
        };
        let after = permute_stickers(&solved_st, &perm);
        // a turn of the solved puzzle is a valid position (the fallback is never taken)
        let mc = stickers_to_cube_with(&d, &after);
        debug_assert!(mc.is_some(), "move {m} breaks the cube");
        let mut mc = mc.unwrap_or(d.solved);
        // centre sources from the sticker permutation (colours repeat, so they cannot be
        // inferred): the centre landing on position i came from the sticker src[a_pos[i]]
        let mut src = [0usize; 72];
        for (i, &p) in perm.iter().enumerate() {
            src[p] = i;
        }
        for i in 0..12 {
            mc.xa[i] = centre_index(src[d.a_pos[i]]) as u8;
            mc.xb[i] = centre_index(src[d.b_pos[i]]) as u8;
        }
        d.moves[m] = mc;
        d.move_perms.push(perm);
    }
    d
}

/// Index of a centre sticker within [`Defs::a_pos`] or [`Defs::b_pos`] (whichever holds it).
pub fn centre_index(sticker: usize) -> usize {
    3 * class_index(sticker / 9) + sticker % 3
}

pub type Stickers = [u8; 72];

pub fn solved_stickers(d: &Defs) -> Stickers {
    let mut s = [0u8; 72];
    for (i, x) in d.sm.st.iter().enumerate() {
        s[i] = x.face as u8;
    }
    s
}

/// `perm[i]` is the destination of sticker i.
pub fn permute_stickers(s: &Stickers, perm: &[usize]) -> Stickers {
    let mut n = [0u8; 72];
    for i in 0..72 {
        n[perm[i]] = s[i];
    }
    n
}

fn stickers_to_cube_with(d: &Defs, s: &Stickers) -> Option<Cube> {
    let mut c = Cube {
        cp: [0; 6],
        co: [0; 6],
        ep: [0; 12],
        xa: [0; 12],
        xb: [0; 12],
    };
    let mut seen = [false; 6];
    for v in 0..6 {
        let cols: [u8; 4] = std::array::from_fn(|j| s[d.corner_slots[v][j]]);
        let mut found = false;
        for k in 0..6 {
            let h = d.corner_cols[k];
            for o in 0..2 {
                if (0..4).all(|j| cols[j] == h[(j + 2 * o) % 4]) {
                    c.cp[v] = k as u8;
                    c.co[v] = o as u8;
                    found = true;
                }
            }
        }
        if !found || seen[c.cp[v] as usize] {
            return None;
        }
        seen[c.cp[v] as usize] = true;
    }
    let mut eseen = [false; 12];
    for e in 0..12 {
        let cols = [s[d.edge_slots[e][0]], s[d.edge_slots[e][1]]];
        let k = (0..12).find(|&k| d.edge_cols[k] == cols)?;
        if eseen[k] {
            return None;
        }
        eseen[k] = true;
        c.ep[e] = k as u8;
    }
    for i in 0..12 {
        c.xa[i] = A_FACES.iter().position(|&f| f == s[d.a_pos[i]] as usize)? as u8;
        c.xb[i] = B_FACES.iter().position(|&f| f == s[d.b_pos[i]] as usize)? as u8;
    }
    Some(c)
}

impl Cube {
    pub fn solved() -> Cube {
        defs().solved
    }

    #[inline]
    #[must_use]
    pub fn apply(&self, m: usize) -> Cube {
        let mv = &defs().moves[m];
        let mut n = *self;
        for d in 0..6 {
            let s = mv.cp[d] as usize;
            n.cp[d] = self.cp[s];
            n.co[d] = self.co[s] ^ mv.co[d];
        }
        for d in 0..12 {
            n.ep[d] = self.ep[mv.ep[d] as usize];
            n.xa[d] = self.xa[mv.xa[d] as usize];
            n.xb[d] = self.xb[mv.xb[d] as usize];
        }
        n
    }

    #[must_use]
    pub fn apply_seq(&self, seq: &[usize]) -> Cube {
        seq.iter().fold(*self, |c, &m| c.apply(m))
    }

    pub fn to_stickers(&self) -> Stickers {
        let d = defs();
        let mut s = [0u8; 72];
        for v in 0..6 {
            let k = self.cp[v] as usize;
            let o = self.co[v] as usize;
            for j in 0..4 {
                s[d.corner_slots[v][j]] = d.corner_cols[k][(j + 2 * o) % 4];
            }
        }
        for e in 0..12 {
            let k = self.ep[e] as usize;
            s[d.edge_slots[e][0]] = d.edge_cols[k][0];
            s[d.edge_slots[e][1]] = d.edge_cols[k][1];
        }
        for i in 0..12 {
            s[d.a_pos[i]] = A_FACES[self.xa[i] as usize] as u8;
            s[d.b_pos[i]] = B_FACES[self.xb[i] as usize] as u8;
        }
        s
    }

    pub fn from_stickers(s: &Stickers) -> Option<Cube> {
        stickers_to_cube_with(defs(), s)
    }

    /// Edge A-triple label at each position.
    #[inline]
    pub fn etri(&self) -> [u8; 12] {
        let d = defs();
        std::array::from_fn(|p| d.edge_atri[self.ep[p] as usize])
    }

    /// Conjugate by a whole-puzzle symmetry `m` that fixes the front corner position
    /// (vertex 0): `m * self * m^-1`, i.e. the same position viewed through `m`. `None` if the
    /// result is not a valid position in the fixed frame (only if `self` is not, or `m` moves
    /// the front corner).
    pub fn conjugate(&self, m: &Mat) -> Option<Cube> {
        let d = defs();
        let st = self.to_stickers();
        let perm = d.sm.whole(m);
        let mut n = [0u8; 72];
        for i in 0..72 {
            let f = face_of_normal(mat_vec(m, FACES[st[i] as usize]));
            n[perm[i]] = f as u8;
        }
        Cube::from_stickers(&n)
    }

    /// Uniformly random reachable position (all four parts are independent).
    pub fn random(rng: &mut impl rand::Rng) -> Cube {
        use rand::seq::SliceRandom;
        let mut c = Cube::solved();
        // corners: random even permutation of positions 1..6, random orientations with even sum
        loop {
            let mut p = [1u8, 2, 3, 4, 5];
            p.shuffle(rng);
            if perm_parity(&p) == 0 {
                c.cp[1..].copy_from_slice(&p);
                break;
            }
        }
        let mut sum = 0;
        for i in 1..5 {
            c.co[i] = rng.gen_range(0..2);
            sum ^= c.co[i];
        }
        c.co[5] = sum;
        loop {
            let mut e: [u8; 12] = std::array::from_fn(|i| i as u8);
            e.shuffle(rng);
            if perm_parity(&e) == 0 {
                c.ep = e;
                break;
            }
        }
        c.xa.shuffle(rng);
        c.xb.shuffle(rng);
        c
    }
}

pub fn perm_parity(p: &[u8]) -> u8 {
    let mut inv = 0u32;
    for i in 0..p.len() {
        for j in 0..i {
            if p[j] > p[i] {
                inv += 1;
            }
        }
    }
    (inv & 1) as u8
}

/// The 90° rotation about the axis through the fixed corner; it swaps the A and B face classes.
pub const SWAP_AB: Mat = [[1, 0, 0], [0, 0, -1], [0, 1, 0]];

/// Internal move index after conjugating by a symmetry that fixes vertex 0 (it maps back faces
/// to back faces, so single, wide and slice turns to turns of the same kind).
pub fn conjugate_move(m: usize, sym: &Mat) -> usize {
    let s = slot(m);
    let f = face_of_normal(mat_vec(sym, FACES[slot_face(s)]));
    debug_assert!(f >= 4);
    let ns = (f - 4) + 4 * (s >> 2);
    let cw = move_cw(m) == (det(sym) == 1);
    2 * ns + usize::from(!cw)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    #[test]
    fn solved_matches_stickers() {
        let d = defs();
        assert_eq!(
            Cube::from_stickers(&solved_stickers(d)),
            Some(Cube::solved())
        );
        for v in 0..6 {
            let at_v = (0..72)
                .filter(|&i| d.sm.st[i].kind == Kind::Corner && dot(d.sm.st[i].p, VERTS[v]) == 7);
            assert_eq!(at_v.count(), 4);
        }
    }

    #[test]
    fn moves_fix_corner_and_have_order_three() {
        let d = defs();
        for m in 0..MAX_MOVES {
            assert_eq!(d.moves[m].cp[0], 0);
            assert_eq!(d.moves[m].co[0], 0);
            let c = Cube::solved().apply(m).apply(m).apply(m);
            assert_eq!(c, Cube::solved());
            assert_eq!(
                Cube::solved().apply(m).apply(inverse_move(m)),
                Cube::solved()
            );
        }
        // slices move no corner at all
        for m in N_MOVES..MAX_MOVES {
            assert_eq!(d.moves[m].cp, Cube::solved().cp);
            assert_eq!(d.moves[m].co, [0; 6]);
            assert_ne!(Cube::solved().apply(m), Cube::solved());
        }
    }

    /// Every internal move, slices included, agrees with the sticker model.
    #[test]
    fn cubie_apply_matches_stickers() {
        let d = defs();
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(1);
        for _ in 0..200 {
            let c = Cube::random(&mut rng);
            let st = c.to_stickers();
            assert_eq!(Cube::from_stickers(&st), Some(c));
            for m in 0..MAX_MOVES {
                let a = c.apply(m);
                let b = Cube::from_stickers(&permute_stickers(&st, &d.move_perms[m])).unwrap();
                assert_eq!(a, b);
            }
        }
    }

    /// A slice is the wide turn followed by the inverse single turn of the same face.
    #[test]
    fn slice_is_wide_then_inverse_single() {
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(4);
        for _ in 0..20 {
            let c = Cube::random(&mut rng);
            for k in 0..4 {
                for dir in 0..2 {
                    let (single, wide, slice) = (2 * k + dir, 2 * (k + 4) + dir, 2 * (k + 8) + dir);
                    assert_eq!(c.apply(slice), c.apply(wide).apply(inverse_move(single)));
                }
            }
        }
    }

    /// The six moves of an axis (single, wide and slice turns of one back face).
    fn axis_moves(k: usize) -> [usize; 6] {
        [
            2 * k,
            2 * k + 1,
            2 * (k + 4),
            2 * (k + 4) + 1,
            2 * (k + 8),
            2 * (k + 8) + 1,
        ]
    }

    #[test]
    fn axis_layers_commute() {
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(3);
        let c = Cube::random(&mut rng);
        for k in 0..4 {
            for a in axis_moves(k) {
                for b in axis_moves(k) {
                    assert_eq!(c.apply(a).apply(b), c.apply(b).apply(a));
                }
            }
        }
    }

    /// The move-order rules derived from the axis group: by breadth-first search over the
    /// moves of one axis, the group has 9 elements (6 single moves and 2 two-move elements
    /// with slices; 4 and 4 without). An allowed same-axis pair must give a two-move element,
    /// and every two-move element must have exactly one allowed pair.
    #[test]
    fn canonical_pairs_match_axis_group() {
        use std::collections::HashMap;
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(5);
        let c = Cube::random(&mut rng);
        for slice in [false, true] {
            for k in 0..4 {
                let moves = &axis_moves(k)[..if slice { 6 } else { 4 }];
                let mut dist = HashMap::from([(c, 0usize)]);
                let mut front = vec![c];
                let mut d = 0;
                while !front.is_empty() {
                    d += 1;
                    let mut next = vec![];
                    for x in &front {
                        for &m in moves {
                            let y = x.apply(m);
                            if let std::collections::hash_map::Entry::Vacant(e) = dist.entry(y) {
                                e.insert(d);
                                next.push(y);
                            }
                        }
                    }
                    front = next;
                }
                let mut hist = [0; 3];
                for &d in dist.values() {
                    hist[d] += 1;
                }
                assert_eq!(hist, if slice { [1, 6, 2] } else { [1, 4, 4] });
                let mut reps: HashMap<Cube, usize> = HashMap::new();
                for &u in moves {
                    for &v in moves {
                        let p = c.apply(u).apply(v);
                        if allowed_after(slice, u, v) {
                            assert_eq!(dist[&p], 2, "pair {u} {v} (slice {slice})");
                            *reps.entry(p).or_default() += 1;
                        }
                    }
                }
                assert_eq!(reps.len(), hist[2]);
                assert!(reps.values().all(|&n| n == 1));
                // moves of different axes are never restricted
                for m in 0..MAX_MOVES {
                    if (m >> 1) & 3 != k {
                        assert!(allowed_after(slice, moves[0], m));
                    }
                }
            }
        }
    }

    #[test]
    fn swap_conjugation_maps_moves() {
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(2);
        let c = Cube::random(&mut rng);
        for m in 0..MAX_MOVES {
            let n = conjugate_move(m, &SWAP_AB);
            let lhs = c.apply(m).conjugate(&SWAP_AB).unwrap();
            let rhs = c.conjugate(&SWAP_AB).unwrap().apply(n);
            assert_eq!(lhs, rhs);
            if m < N_MOVES {
                // the swap exchanges the classes of single and wide turns
                assert!(n < N_MOVES);
                assert_ne!(is_a_move(m), is_a_move(n));
            } else {
                // slices stay slices (and are in neither class)
                assert!(n >= N_MOVES);
                assert!(!is_a_move(m) && !is_a_move(n));
            }
        }
    }
}
