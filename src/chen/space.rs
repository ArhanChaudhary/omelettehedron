//! The phase-2 space of the block split: hexagon configurations and the move rule.
//!
//! Each class-A face (U, F, BR, BL) has a *hexagon*: its 3 centres and the A-face stickers of
//! its 3 edges (everything on the face except the 3 corner tips). A sticker's *colour* here is
//! its A-colour: a centre's colour, or the A-face of the edge piece it belongs to.
//!
//! * Phase 1 reaches `H'`: every hexagon holds its own face's 3 centres and 3 edges, turned by
//!   any multiple of 120° (up to a whole-puzzle turn: in the fixed-corner frame the blocks may
//!   sit in any of the 12 rotated frames, as in `H_A`). That is exactly "A-centres and edges
//!   solved" for the two-phase solver's A-centre/edge (AE) tables.
//! * Phase 2 goes from `H'` to solved with a restricted move set. A class-A move (an internal
//!   move in [`A_SLOTS`]) turns whole hexagons and is always allowed. A class-B move cuts each
//!   of the three hexagons next to its turning layer into the 3 stickers inside the layer and
//!   the 3 outside; it is allowed only if, on each of them, the 3 inside have one colour and the
//!   3 outside have one colour, so that it never splits a colour group further. (A wide move's
//!   cut is that of the opposite face's layer, which stays still.)
//!
//! The *configuration* of a position is its A-centre colouring with its edge arrangement (the
//! edges are pieces, so this also fixes their B-stickers); its *shape* is the colours of the 24
//! hexagon stickers alone. Whether a move is allowed depends only on the shape.
//!
//! Seen without a fixed corner (the *pattern view*: the 8 face turns, nothing ever rotates the
//! puzzle), the allowed moves reach 37,260 configurations from solved, in [`N_PSHAPE`] shapes
//! of 81 each, and never move a block to another face. In the fixed-corner frame a wide turn
//! also turns the whole puzzle, so the blocks are seen in any of 12 frames: [`N_CFG`] = 12 x
//! 37,260 configurations, none of them a turn of another. So the phase-2 space has 447,120 x
//! 960 corner states x 369,600 B-centre colourings = 1.59 x 10^14 positions; the pattern view
//! counts the same (37,260 x 11,520 corner states x 369,600), with no two positions a turn of
//! each other. The SB table ([`super::tables`]) works in the pattern view, through the frame
//! of each configuration ([`Space::rot`]).

use crate::cube::*;
use crate::geometry::{
    class_index, face_of_normal, mat_vec, rotation_group, transpose, Mat, A_FACES, B_FACES, FACES,
};
use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};
use std::sync::OnceLock;

/// A multiply-rotate hasher with the splitmix64 finaliser, for the small keys of the space's
/// maps (MEASURE: against the standard `SipHash`).
#[derive(Default, Clone, Copy)]
pub struct QuickHash(u64);

impl Hasher for QuickHash {
    fn finish(&self) -> u64 {
        let mut x = self.0;
        x = (x ^ x >> 30).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        x = (x ^ x >> 27).wrapping_mul(0x94d0_49bb_1331_11eb);
        x ^ x >> 31
    }
    fn write(&mut self, bytes: &[u8]) {
        for c in bytes.chunks(8) {
            let mut w = [0u8; 8];
            w[..c.len()].copy_from_slice(c);
            self.write_u64(u64::from_le_bytes(w));
        }
    }
    fn write_u64(&mut self, x: u64) {
        self.0 = (self.0.rotate_left(5) ^ x).wrapping_mul(0x517c_c1b7_2722_0a95);
    }
    fn write_u128(&mut self, x: u128) {
        self.write_u64(x as u64);
        self.write_u64((x >> 64) as u64);
    }
}

/// A hash map with [`QuickHash`].
type QuickMap<K, V> = HashMap<K, V, BuildHasherDefault<QuickHash>>;

/// Configurations reachable from solved with the allowed moves in the fixed-corner frame
/// (checked when built).
pub const N_CFG: usize = 447_120;
/// Their colour shapes in the pattern view.
pub const N_PSHAPE: usize = 460;
/// The configurations of `H'` (12 frames x 81 block turns; the same as `H_A`'s edge states).
pub const N_HPRIME: usize = 972;
/// Marks a move that is not allowed, or a configuration outside `H'`.
pub const NONE: u32 = u32::MAX;
/// Pattern-view moves: the face turns of the 8 faces, `2 * face + dir` (0 clockwise).
pub const N_PMOVES: usize = 16;

/// The configuration space with its move tables.
pub struct Space {
    /// configuration after each move, or [`NONE`] where the move is not allowed
    pub next: Vec<[u32; N_MOVES]>,
    /// the A-centre colouring and edge arrangement of each configuration
    pub cfgs: Vec<([u8; 12], [u8; 12])>,
    /// configuration id of each (A-centre colouring, edge arrangement) key ([`key`])
    ids: QuickMap<u128, u32>,
    /// edge id of the two-phase solver's `H_A` table ([`crate::tables::Phase2`]) of each
    /// configuration in `H'`, [`NONE`] for the others
    pub p2e: Vec<u32>,
    /// the frame of each configuration: the rotation (index into [`Space::rot_xb`]) that turns
    /// the position into the pattern view with the blocks at home
    pub rot: Vec<u8>,
    /// pattern-view shape of each configuration (seen through its frame)
    pub pshape: Vec<u16>,
    /// pattern-view shape after each pattern-view move, or `u16::MAX` where it is not allowed
    pub pshape_next: Vec<[u16; N_PMOVES]>,
    /// per frame: the gather shuffle that turns the B-centre colouring into the pattern view
    /// (see [`Space::pattern_xb`])
    pub rot_xb: Vec<[u8; 16]>,
    /// the pattern-view configurations and their symmetry classes
    pub pat: Patterns,
}

/// The configurations in the pattern view (the 37,260 reachable from solved there, which are
/// the configurations of frame 0), with their moves and the classes of the 24 symmetries that
/// keep the face classes (`T_d`; conjugation, so solved and the move rule stay the same).
pub struct Patterns {
    /// pattern-view configuration of each configuration (seen through its frame)
    pub of_cfg: Vec<u16>,
    /// the configuration id (frame 0) of each pattern-view configuration
    pub cfg: Vec<u32>,
    /// pattern-view configuration after each internal move (the move made in frame 0 and seen
    /// through the frame it leads to: the 16 internal moves give the 16 face turns), or
    /// `u16::MAX` where the move is not allowed
    pub next: Vec<[u16; N_MOVES]>,
    /// the B-centre gathers (`new[p] = old[src[p]]`) of the same moves
    pub xb_src: [[u8; 12]; N_MOVES],
    /// symmetry class of each pattern-view configuration, and the symmetry (index into
    /// [`Patterns::sym_xb`]) that maps it to its class representative
    pub class: Vec<(u16, u8)>,
    /// the representative of each class
    pub rep: Vec<u16>,
    /// per class: the symmetries other than the identity that fix the representative
    pub stab: Vec<Vec<u8>>,
    /// per symmetry: the `(relabel, gather)` shuffles of a B-centre colouring
    pub sym_xb: Vec<([u8; 16], [u8; 16])>,
}

impl Patterns {
    /// A B-centre colouring conjugated by symmetry `k`.
    #[inline]
    pub fn sym_xb(&self, k: usize, xb: &[u8; 16]) -> [u8; 16] {
        let (relabel, gather) = self.sym_xb.get(k).copied().unwrap_or(([0; 16], [0x80; 16]));
        crate::fast::shuffle(&crate::fast::shuffle(&relabel, xb), &gather)
    }
}

/// The places and colours of a whole-puzzle symmetry `m` that keeps the face classes: where
/// each A-centre, edge and B-centre place goes, and which colour (A, B) or edge piece each
/// becomes when the colours are mapped too.
struct SymMaps {
    a: [usize; 12],
    e: [usize; 12],
    b: [usize; 12],
    acol: [u8; 4],
    bcol: [u8; 4],
    piece: [u8; 12],
}

fn sym_maps(m: &Mat) -> SymMaps {
    let d = defs();
    let sm = &d.sm;
    let to = |s: usize| sm.index_of(mat_vec(m, sm.st[s].p));
    let face = |f: usize| face_of_normal(mat_vec(m, FACES[f]));
    let edge_at = |s: usize| (0..12).find(|&e| d.edge_slots[e][0] == s).unwrap_or(0);
    let piece_of = |cols: [u8; 2]| (0..12).find(|&k| d.edge_cols[k] == cols).unwrap_or(0) as u8;
    SymMaps {
        a: std::array::from_fn(|p| centre_index(to(d.a_pos[p]))),
        e: std::array::from_fn(|e| edge_at(to(d.edge_slots[e][0]))),
        b: std::array::from_fn(|i| centre_index(to(d.b_pos[i]))),
        acol: std::array::from_fn(|c| class_index(face(A_FACES[c])) as u8),
        bcol: std::array::from_fn(|c| class_index(face(B_FACES[c])) as u8),
        piece: std::array::from_fn(|k| {
            let [a, b] = d.edge_cols[k];
            piece_of([face(a as usize) as u8, face(b as usize) as u8])
        }),
    }
}

impl SymMaps {
    /// A configuration with the whole puzzle turned (`recolour` false: the stickers keep their
    /// colours) or conjugated (`recolour` true: the colours and pieces are mapped too).
    fn apply(&self, xa: &[u8; 12], ep: &[u8; 12], recolour: bool) -> ([u8; 12], [u8; 12]) {
        let mut a = [0u8; 12];
        let mut e = [0u8; 12];
        for p in 0..12 {
            a[self.a[p]] = if recolour {
                self.acol[xa[p] as usize & 3]
            } else {
                xa[p]
            };
            e[self.e[p]] = if recolour {
                self.piece[ep[p] as usize % 12]
            } else {
                ep[p]
            };
        }
        (a, e)
    }

    /// The `(relabel, gather)` shuffles of a B-centre colouring conjugated by the symmetry.
    fn xb_shuffles(&self) -> ([u8; 16], [u8; 16]) {
        let mut relabel = [0u8; 16];
        relabel[..4].copy_from_slice(&self.bcol);
        let mut gather = [0x80u8; 16];
        for i in 0..12 {
            gather[self.b[i]] = i as u8;
        }
        (relabel, gather)
    }
}

/// The pattern-view configurations of a space (see [`Patterns`]).
fn patterns(
    cfgs: &[([u8; 12], [u8; 12])],
    ids: &QuickMap<u128, u32>,
    next: &[[u32; N_MOVES]],
    rot: &[u8],
    rots: &[Mat],
    rot_xb: &[[u8; 16]],
) -> Patterns {
    let d = defs();
    let id = |(xa, ep): ([u8; 12], [u8; 12])| ids.get(&key(&xa, &ep)).copied().unwrap_or(0);
    // the configurations of frame 0, numbered in order
    let cfg: Vec<u32> = (0..cfgs.len() as u32)
        .filter(|&c| rot[c as usize] == 0)
        .collect();
    let mut dense = vec![u16::MAX; cfgs.len()];
    for (i, &c) in cfg.iter().enumerate() {
        dense[c as usize] = i as u16;
    }
    let turn: Vec<SymMaps> = rots.iter().map(sym_maps).collect();
    let of_cfg: Vec<u16> = cfgs
        .iter()
        .zip(rot)
        .map(|(&(xa, ep), &r)| {
            let t = &turn[r as usize % turn.len()];
            dense[id(t.apply(&xa, &ep, false)) as usize]
        })
        .collect();
    let pnext: Vec<[u16; N_MOVES]> = cfg
        .iter()
        .map(|&c| {
            std::array::from_fn(|m| {
                let n = next[c as usize][m];
                if n == NONE {
                    u16::MAX
                } else {
                    of_cfg[n as usize]
                }
            })
        })
        .collect();
    // each move's frame change is the same for every configuration of frame 0: take it from
    // solved, where every move is allowed
    let xb_src: [[u8; 12]; N_MOVES] = std::array::from_fn(|m| {
        let g = rot_xb[rot[next[0][m] as usize] as usize % rot_xb.len()];
        std::array::from_fn(|j| d.moves[m].xb[g[j] as usize % 12])
    });
    // symmetry classes under T_d
    let syms: Vec<SymMaps> = crate::geometry::class_preserving_group()
        .iter()
        .map(sym_maps)
        .collect();
    let image = |p: usize, k: usize| -> u16 {
        let (xa, ep) = cfgs[cfg[p] as usize];
        dense[id(syms[k].apply(&xa, &ep, true)) as usize]
    };
    let mut class = vec![(u16::MAX, 0u8); cfg.len()];
    let mut rep = vec![];
    let mut stab = vec![];
    for p in 0..cfg.len() {
        if class[p].0 != u16::MAX {
            continue;
        }
        // p is the smallest of its class: the representative
        let c = rep.len() as u16;
        rep.push(p as u16);
        let mut st = vec![];
        for k in 0..syms.len() {
            let q = image(p, k) as usize;
            if q == p && k > 0 {
                st.push(k as u8);
            }
            // the symmetry mapping q to p is the inverse of k
            if let Some(slot) = class.get_mut(q) {
                if slot.0 == u16::MAX {
                    *slot = (c, inverse_sym(&syms, k) as u8);
                }
            }
        }
        stab.push(st);
    }
    Patterns {
        of_cfg,
        cfg,
        next: pnext,
        xb_src,
        class,
        rep,
        stab,
        sym_xb: syms.iter().map(SymMaps::xb_shuffles).collect(),
    }
}

/// The index of the inverse of symmetry `k` (as place maps).
fn inverse_sym(syms: &[SymMaps], k: usize) -> usize {
    (0..syms.len())
        .find(|&j| (0..12).all(|p| syms[j].a[syms[k].a[p]] == p && syms[j].e[syms[k].e[p]] == p))
        .unwrap_or(0)
}

/// Key of a configuration: the A-centre colouring (2 bits per position) and the edge
/// arrangement (4 bits per position).
pub fn key(xa: &[u8; 12], ep: &[u8; 12]) -> u128 {
    let mut k = 0u128;
    for p in 0..12 {
        k |= u128::from(xa[p] & 3) << (2 * p);
        k |= u128::from(ep[p] & 15) << (24 + 4 * p);
    }
    k
}

/// The 24 hexagon positions: `0..12` the A-centre positions (indices into `Cube::xa`),
/// `12..24` the edge positions (`12 + e`, the edge's A-face sticker).
const N_HEX_POS: usize = 24;

/// The hexagons: the 6 positions (see [`N_HEX_POS`]) of each A-face, in class order.
fn hexagons() -> [[usize; 6]; 4] {
    let d = defs();
    std::array::from_fn(|i| {
        let mut h = [0usize; 6];
        for j in 0..3 {
            h[j] = 3 * i + j;
        }
        let mut k = 3;
        for e in 0..12 {
            if d.edge_slots[e][0] / 9 == A_FACES[i] && k < 6 {
                h[k] = 12 + e;
                k += 1;
            }
        }
        h
    })
}

/// The sticker index of each hexagon position.
fn hex_sticker(p: usize) -> usize {
    let d = defs();
    if p < 12 {
        d.a_pos[p]
    } else {
        d.edge_slots[p - 12][0]
    }
}

/// The hexagon position of each sticker (`N_HEX_POS` for the others).
fn hex_positions() -> [usize; 72] {
    let mut out = [N_HEX_POS; 72];
    for p in 0..N_HEX_POS {
        out[hex_sticker(p)] = p;
    }
    out
}

/// The cut of internal move `m` (see [`cuts_of`]).
fn cuts(m: usize) -> Vec<(u32, u32)> {
    cuts_of(&defs().move_perms[m])
}

/// The cut of a move with sticker permutation `perm` (`perm[i]` = destination of sticker i):
/// for each hexagon it splits, the positions inside the turning layer (the stickers it moves)
/// and the positions outside, as bit masks over the 24 hexagon positions.
fn cuts_of(perm: &[usize]) -> Vec<(u32, u32)> {
    hexagons()
        .iter()
        .filter_map(|h| {
            let (mut moved, mut still) = (0u32, 0u32);
            for &p in h {
                let s = hex_sticker(p);
                if perm[s] == s {
                    still |= 1 << p;
                } else {
                    moved |= 1 << p;
                }
            }
            (moved != 0 && still != 0).then_some((moved, still))
        })
        .collect()
}

/// The colours of the 24 hexagon positions of a configuration.
pub fn colours(xa: &[u8; 12], ep: &[u8; 12]) -> [u8; N_HEX_POS] {
    let d = defs();
    std::array::from_fn(|p| {
        if p < 12 {
            xa[p]
        } else {
            d.edge_atri[ep[p - 12] as usize]
        }
    })
}

/// Whether every position of `mask` has the same colour.
fn one_colour(col: &[u8; N_HEX_POS], mask: u32) -> bool {
    let mut c = None;
    (0..N_HEX_POS).filter(|&p| mask >> p & 1 == 1).all(|p| {
        let first = *c.get_or_insert(col[p]);
        first == col[p]
    })
}

/// Whether the move with cut `cut` ([`cuts`]) is allowed on a configuration with colours `col`.
fn allowed(col: &[u8; N_HEX_POS], cut: &[(u32, u32)]) -> bool {
    cut.iter()
        .all(|&(a, b)| one_colour(col, a) && one_colour(col, b))
}

/// The sticker permutations of the pattern-view moves (`2 * face + dir`).
fn pattern_moves() -> Vec<Vec<usize>> {
    let sm = &defs().sm;
    (0..8)
        .flat_map(|f| [true, false].map(|cw| sm.face_turn(f, cw)))
        .collect()
}

/// A shape after the sticker permutation `perm` (which maps hexagon stickers to hexagon
/// stickers; `pos` from [`hex_positions`]).
fn permute_shape(col: &[u8; N_HEX_POS], perm: &[usize], pos: &[usize; 72]) -> [u8; N_HEX_POS] {
    let mut out = [0u8; N_HEX_POS];
    for p in 0..N_HEX_POS {
        out[pos[perm[hex_sticker(p)]] % N_HEX_POS] = col[p];
    }
    out
}

/// A shape seen with the whole puzzle turned by the rotation `m`: the stickers move to the
/// turned places and keep their colours.
fn rotate_shape(col: &[u8; N_HEX_POS], m: &Mat, pos: &[usize; 72]) -> [u8; N_HEX_POS] {
    let sm = &defs().sm;
    let mut out = [0u8; N_HEX_POS];
    for p in 0..N_HEX_POS {
        let to = sm.index_of(mat_vec(m, sm.st[hex_sticker(p)].p));
        out[pos[to] % N_HEX_POS] = col[p];
    }
    out
}

/// The gather shuffle (`new[p] = old[g[p]]`, lanes 12..16 cleared) that turns a B-centre
/// colouring with the whole puzzle turned by the rotation `m`.
fn rotate_xb(m: &Mat) -> [u8; 16] {
    let d = defs();
    let mut gather = [0x80u8; 16];
    for i in 0..12 {
        let to = d.sm.index_of(mat_vec(m, d.sm.st[d.b_pos[i]].p));
        gather[centre_index(to)] = i as u8;
    }
    gather
}

impl Space {
    fn new() -> Space {
        let s0 = Cube::solved();
        let cut: Vec<Vec<(u32, u32)>> = (0..N_MOVES).map(cuts).collect();
        // the allowed moves of each shape, as a bit mask (computed once per shape)
        let mut legal: QuickMap<[u8; N_HEX_POS], u16> = QuickMap::default();
        // breadth-first from solved over (A-centres, edges)
        let mut cfgs = vec![(s0.xa, s0.ep)];
        let mut ids: QuickMap<u128, u32> = QuickMap::default();
        ids.insert(key(&s0.xa, &s0.ep), 0);
        let mut next: Vec<[u32; N_MOVES]> = vec![];
        let mut i = 0;
        while i < cfgs.len() {
            let (xa, ep) = cfgs[i];
            i += 1;
            let col = colours(&xa, &ep);
            let ok = *legal.entry(col).or_insert_with(|| {
                (0..N_MOVES)
                    .filter(|&m| allowed(&col, &cut[m]))
                    .fold(0, |a, m| a | 1 << m)
            });
            let c = Cube { ep, xa, ..s0 };
            let mut row = [NONE; N_MOVES];
            for m in 0..N_MOVES {
                if ok >> m & 1 == 0 {
                    continue;
                }
                let n = c.apply(m);
                let k = key(&n.xa, &n.ep);
                let id = if let Some(&id) = ids.get(&k) {
                    id
                } else {
                    let id = cfgs.len() as u32;
                    ids.insert(k, id);
                    cfgs.push((n.xa, n.ep));
                    id
                };
                row[m] = id;
            }
            next.push(row);
        }
        // pattern-view shapes, breadth-first from solved (solved first)
        let pos = hex_positions();
        let pm = pattern_moves();
        let pcut: Vec<Vec<(u32, u32)>> = pm.iter().map(|p| cuts_of(p)).collect();
        let solved_col = colours(&s0.xa, &s0.ep);
        let mut pshapes = vec![solved_col];
        let mut pids: HashMap<[u8; N_HEX_POS], u16> = HashMap::from([(solved_col, 0)]);
        let mut pshape_next: Vec<[u16; N_PMOVES]> = vec![];
        let mut i = 0;
        while i < pshapes.len() {
            let col = pshapes[i];
            i += 1;
            let mut row = [u16::MAX; N_PMOVES];
            for m in 0..N_PMOVES {
                if allowed(&col, &pcut[m]) {
                    let n = permute_shape(&col, &pm[m], &pos);
                    let id = *pids.entry(n).or_insert_with(|| {
                        pshapes.push(n);
                        (pshapes.len() - 1) as u16
                    });
                    row[m] = id;
                }
            }
            pshape_next.push(row);
        }
        // the frame of each configuration: the rotation that turns its shape into a
        // pattern-view shape (one per shape; the 12 frames hold disjoint shapes)
        let rots: Vec<Mat> = rotation_group().iter().map(transpose).collect();
        let mut frame_of: QuickMap<[u8; N_HEX_POS], (u8, u16)> = QuickMap::default();
        let mut rot = Vec::with_capacity(cfgs.len());
        let mut pshape = Vec::with_capacity(cfgs.len());
        for (xa, ep) in &cfgs {
            let col = colours(xa, ep);
            let (r, ps) = *frame_of.entry(col).or_insert_with(|| {
                rots.iter()
                    .enumerate()
                    .find_map(|(r, m)| {
                        pids.get(&rotate_shape(&col, m, &pos))
                            .map(|&ps| (r as u8, ps))
                    })
                    .unwrap_or((0, 0))
            });
            rot.push(r);
            pshape.push(ps);
        }
        let rot_xb: Vec<[u8; 16]> = rots.iter().map(rotate_xb).collect();
        let pat = patterns(&cfgs, &ids, &next, &rot, &rots, &rot_xb);
        // the configurations in H': in the two-phase solver's H_A table with corners and
        // B-centres solved
        let p2 = crate::tables::Phase2::new();
        let p2e = cfgs
            .iter()
            .map(|&(xa, ep)| {
                p2.index(&Cube { ep, xa, ..s0 })
                    .map_or(NONE, |i| (i % crate::tables::N_P2_EDGE) as u32)
            })
            .collect();
        Space {
            next,
            cfgs,
            ids,
            p2e,
            rot,
            pshape,
            pshape_next,
            rot_xb,
            pat,
        }
    }

    /// Configuration id of a position, if its configuration is reachable.
    pub fn id(&self, cube: &Cube) -> Option<u32> {
        self.ids.get(&key(&cube.xa, &cube.ep)).copied()
    }

    /// Whether configuration `cfg` is in `H'`.
    #[inline]
    pub fn in_hprime(&self, cfg: u32) -> bool {
        self.p2e.get(cfg as usize).is_some_and(|&e| e != NONE)
    }

    /// The B-centre colouring `xb` (labels in lanes 0..12) of a position with configuration
    /// `cfg`, seen in the pattern view through the configuration's frame.
    #[inline]
    pub fn pattern_xb(&self, cfg: u32, xb: &[u8; 16]) -> [u8; 16] {
        let r = self.rot.get(cfg as usize).copied().unwrap_or(0) as usize;
        let gather = self.rot_xb.get(r).copied().unwrap_or([0x80; 16]);
        crate::fast::shuffle(xb, &gather)
    }
}

/// The gather vectors (`new[p] = old[src[p]]`) of the pattern-view moves on the B-centre
/// positions.
pub fn pattern_xb_moves() -> [[u8; 12]; N_PMOVES] {
    let d = defs();
    let pm = pattern_moves();
    std::array::from_fn(|m| {
        let mut src = [0u8; 12];
        for i in 0..12 {
            src[centre_index(pm[m][d.b_pos[i]])] = i as u8;
        }
        src
    })
}

/// The configuration space, built once (about a second).
pub fn space() -> &'static Space {
    static S: OnceLock<Space> = OnceLock::new();
    S.get_or_init(Space::new)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coord::{corners, mset};
    use std::collections::HashSet;

    /// The hexagons partition the 24 positions, every class-B move cuts three of them 3 + 3
    /// (one edge and two centres inside a single turn's layer; a wide turn moves the other
    /// three), and class-A moves cut none.
    #[test]
    fn cuts_follow_the_rule() {
        let hx = hexagons();
        let all: HashSet<usize> = hx.iter().flatten().copied().collect();
        assert_eq!(all.len(), N_HEX_POS);
        for m in 0..N_MOVES {
            let c = cuts(m);
            if is_a_move(m) {
                assert!(c.is_empty(), "move {m}");
            } else {
                assert_eq!(c.len(), 3, "move {m}");
                for &(inside, outside) in &c {
                    assert_eq!(inside.count_ones(), 3);
                    assert_eq!(outside.count_ones(), 3);
                    let edges_in = (inside >> 12).count_ones();
                    assert_eq!(edges_in, if slot_is_wide(slot(m)) { 2 } else { 1 });
                }
            }
        }
    }

    /// 447,120 configurations; 460 pattern-view shapes and 37,260 pattern-view configurations;
    /// the configurations of each frame and pattern shape are 81; `H'` is 972 of them; every allowed move can be undone (the tables'
    /// breadth-first search needs an undirected graph).
    #[test]
    fn configuration_counts() {
        let s = space();
        assert_eq!(s.cfgs.len(), N_CFG);
        assert_eq!(s.pshape_next.len(), N_PSHAPE);
        // the pattern view: 37,260 configurations (460 shapes x 81), 1,611 symmetry classes
        assert_eq!(s.pat.cfg.len(), 37_260);
        assert_eq!(
            s.pat.rep.len() * crate::coord::N_MSET,
            crate::chen::tables::N_PX
        );
        let mut per = HashMap::new();
        for c in 0..N_CFG {
            *per.entry((s.rot[c], s.pshape[c])).or_insert(0) += 1;
        }
        assert_eq!(per.len(), 12 * N_PSHAPE);
        assert!(per.values().all(|&n| n == 81));
        assert_eq!(s.p2e.iter().filter(|&&e| e != NONE).count(), N_HPRIME);
        for (id, row) in s.next.iter().enumerate() {
            for m in 0..N_MOVES {
                if row[m] != NONE {
                    assert_eq!(s.next[row[m] as usize][inverse_move(m)], id as u32);
                }
                if is_a_move(m) {
                    assert_ne!(row[m], NONE);
                }
            }
        }
        for (ps, row) in s.pshape_next.iter().enumerate() {
            for m in 0..N_PMOVES {
                if row[m] != u16::MAX {
                    assert_eq!(s.pshape_next[row[m] as usize][m ^ 1], ps as u16);
                }
            }
        }
    }

    /// The configurations of H' are exactly those of `H_A`, and every move is allowed there.
    #[test]
    fn hprime_is_the_ha_frame_set() {
        let s = space();
        let mut rng = <rand_xoshiro::Xoshiro256PlusPlus as rand::SeedableRng>::seed_from_u64(9);
        let a_moves: Vec<usize> = (0..N_MOVES).filter(|&m| is_a_move(m)).collect();
        for _ in 0..2000 {
            let mut c = Cube::solved();
            for _ in 0..30 {
                c = c.apply(a_moves[rand::Rng::gen_range(&mut rng, 0..a_moves.len())]);
            }
            let id = s.id(&c).unwrap();
            assert!(s.in_hprime(id));
            assert!(s.next[id as usize].iter().all(|&n| n != NONE));
            // in H' every hexagon has one colour: the pattern shape is the solved one
            assert_eq!(s.pshape[id as usize], 0);
        }
    }

    /// A fixed-frame move is a pattern-view move seen through the frame: the pattern shape and
    /// the B-centres seen through the frame follow a pattern-view move.
    #[test]
    fn frames_follow_the_moves() {
        let s = space();
        let pxm = pattern_xb_moves();
        let ms = mset();
        let mut rng = <rand_xoshiro::Xoshiro256PlusPlus as rand::SeedableRng>::seed_from_u64(3);
        for _ in 0..3000 {
            let cfg = rand::Rng::gen_range(&mut rng, 0..N_CFG as u32);
            let mut xb = [0u8; 16];
            let mut l = Cube::solved().xb;
            rand::seq::SliceRandom::shuffle(&mut l[..], &mut rng);
            xb[..12].copy_from_slice(&l);
            let px = s.pattern_xb(cfg, &xb);
            for m in 0..N_MOVES {
                let n = s.next[cfg as usize][m];
                if n == NONE {
                    continue;
                }
                let c = Cube {
                    xb: l,
                    ..Cube::solved()
                }
                .apply(m);
                let mut xb2 = [0u8; 16];
                xb2[..12].copy_from_slice(&c.xb);
                let px2 = s.pattern_xb(n, &xb2);
                // some pattern-view move takes (shape, B-centres) to the child's
                let ok = (0..N_PMOVES).any(|pm| {
                    let gathered: [u8; 12] = std::array::from_fn(|p| px[pxm[pm][p] as usize]);
                    s.pshape_next[s.pshape[cfg as usize] as usize][pm] == s.pshape[n as usize]
                        && ms.rank(&gathered) == ms.rank(&std::array::from_fn(|p| px2[p]))
                });
                assert!(ok, "cfg {cfg} move {m}");
            }
        }
    }

    /// Breadth-first search over whole positions in the pattern view (all 8 faces turned, both
    /// directions, as sticker permutations; centres of a colour alike), with the rule of this
    /// module applied to the hexagon stickers: 1, 16, 160, 1344, 10464 positions at distances
    /// 0..4, as an independent sticker model gave.
    #[test]
    fn pattern_bfs_counts() {
        type Move = (Vec<(u32, u32)>, Vec<usize>);
        let moves: Vec<Move> = pattern_moves()
            .into_iter()
            .map(|perm| (cuts_of(&perm), perm))
            .collect();
        let sm = &defs().sm;
        let solved: [u8; 72] = std::array::from_fn(|i| sm.st[i].face as u8);
        let mut seen: HashSet<[u8; 72]> = HashSet::from([solved]);
        let mut front = vec![solved];
        let mut counts = vec![1usize];
        for _ in 0..4 {
            let mut next = vec![];
            for st in &front {
                let col: [u8; N_HEX_POS] = std::array::from_fn(|p| st[hex_sticker(p)]);
                for (cut, perm) in &moves {
                    if !allowed(&col, cut) {
                        continue;
                    }
                    let mut t = [0u8; 72];
                    for i in 0..72 {
                        t[perm[i]] = st[i];
                    }
                    if seen.insert(t) {
                        next.push(t);
                    }
                }
            }
            counts.push(next.len());
            front = next;
        }
        assert_eq!(counts, [1, 16, 160, 1344, 10464]);
    }

    /// The same counts in the fixed-corner frame, from the configuration move table, the corner
    /// table and the B-centre colouring (the solver's phase-2 state): the 16 internal moves are
    /// the pattern view's 16 face turns up to rotation.
    #[test]
    fn fixed_frame_bfs_counts() {
        type State = (u32, u16, u32);
        let s = space();
        let s0 = Cube::solved();
        let ms = mset();
        let cs = corners();
        let state = |cfg: u32, c: &Cube| -> State { (cfg, cs.id(c), ms.rank(&c.xb)) };
        let mut seen: HashSet<State> = HashSet::from([state(0, &s0)]);
        let mut front = vec![(0u32, s0)];
        let mut counts = vec![1usize];
        for _ in 0..4 {
            let mut next = vec![];
            for (cfg, c) in &front {
                for m in 0..N_MOVES {
                    let n = s.next[*cfg as usize][m];
                    if n == NONE {
                        continue;
                    }
                    let c2 = c.apply(m);
                    assert_eq!(s.id(&c2), Some(n));
                    if seen.insert(state(n, &c2)) {
                        next.push((n, c2));
                    }
                }
            }
            counts.push(next.len());
            front = next;
        }
        assert_eq!(counts, [1, 16, 160, 1344, 10464]);
    }
}
