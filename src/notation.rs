//! Physical move notation and conversion to/from internal (fixed-corner) moves.
//!
//! `X` is a turn of face X, `Xw` a wide (two-layer) turn of X and `Xs` a slice turn: the
//! middle layer between X and its opposite face, turned like X (so `Us` = `Ds'`). Solutions
//! name slices by their front face (U, R, F, L).

use crate::cube::*;
use crate::geometry::*;
use std::fmt;

/// Which layers of a face's axis a physical move turns.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Layer {
    /// the face's own layer
    #[default]
    Face,
    /// the face's layer and the middle layer
    Wide,
    /// the middle layer only
    Slice,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PhysMove {
    pub face: usize,
    pub layer: Layer,
    pub cw: bool,
}

impl PhysMove {
    /// A turn of the single layer of `face`.
    pub fn face(face: usize, cw: bool) -> PhysMove {
        PhysMove {
            face,
            layer: Layer::Face,
            cw,
        }
    }
}

impl fmt::Display for PhysMove {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}{}{}",
            FACE_NAMES[self.face],
            match self.layer {
                Layer::Face => "",
                Layer::Wide => "w",
                Layer::Slice => "s",
            },
            if self.cw { "" } else { "'" }
        )
    }
}

pub fn format_moves(moves: &[PhysMove]) -> String {
    moves
        .iter()
        .map(ToString::to_string)
        .collect::<Vec<_>>()
        .join(" ")
}

/// Parse a scramble such as `U R' BL F2 D Rw Us'`. `X2` is a double turn (= two moves),
/// `X'`/`X3` variants of the prime are accepted. A `w` after the face name makes a wide turn
/// and an `s` a slice turn; both are always accepted, whatever move set the solver uses.
///
/// # Errors
///
/// If a token is not a face name with an optional `w` or `s` and an optional `'`, `2`, `2'`
/// or `3` suffix.
pub fn parse_moves(s: &str) -> Result<Vec<PhysMove>, String> {
    let mut out = vec![];
    for tok in s
        .split(|c: char| c.is_whitespace() || c == ',')
        .filter(|t| !t.is_empty())
    {
        let t = tok.replace('’', "'");
        let (name, suffix) = match t.find(|c: char| c == '\'' || c.is_ascii_digit()) {
            Some(i) => t.split_at(i),
            None => (t.as_str(), ""),
        };
        // (no face name ends in w, s or v)
        let (name, layer, rotation) = match name.char_indices().last() {
            Some((i, 'w' | 'W')) => (&name[..i], Layer::Wide, false),
            Some((i, 's' | 'S')) => (&name[..i], Layer::Slice, false),
            Some((i, 'v' | 'V')) => (&name[..i], Layer::Wide, true),
            _ => (name, Layer::Face, false),
        };
        let face = FACE_NAMES
            .iter()
            .position(|&n| n.eq_ignore_ascii_case(name))
            .ok_or_else(|| format!("unknown move `{tok}`"))?;
        let (n, cw) = match suffix {
            "" | "1" => (1, true),
            "'" | "1'" | "3" => (1, false),
            "2" => (2, true),
            "2'" => (2, false),
            _ => return Err(format!("unknown move suffix in `{tok}`")),
        };
        for _ in 0..n {
            out.push(PhysMove { face, layer, cw });
            if rotation {
                // `Xv` turns the whole puzzle like X: all three layers of X's axis, i.e. the
                // wide turn Xw and the opposite face turned the other way (as seen from outside
                // it), which commute
                out.push(PhysMove {
                    face: opposite(face),
                    layer: Layer::Face,
                    cw: !cw,
                });
            }
        }
    }
    Ok(out)
}

/// Orientation of the internal frame: maps internal coordinates to physical coordinates.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Frame(pub Mat);

impl Frame {
    pub fn identity() -> Frame {
        Frame(IDENTITY)
    }

    /// Physical face at internal face `f`.
    fn physical_face(&self, f: usize) -> usize {
        face_of_normal(mat_vec(&self.0, FACES[f]))
    }

    /// Internal move equivalent to a physical move, and the new frame.
    ///
    /// Count the layers of an axis from its back face X (back layer, slice, front layer) and
    /// turns in the sense of X. A move turns them by (a, b, c) thirds; the whole-puzzle
    /// rotation by c is free, so it is the internal move (a - c, b - c, 0) plus a rotation:
    /// X = (1, 0, 0) is the single turn, Xw = (1, 1, 0) the wide turn and Xs = (0, 1, 0) the
    /// slice turn of X. A turn of the front face Y (opposite X) is (0, 0, -1) = Xw plus a
    /// rotation, Yw = (0, -1, -1) is X plus a rotation, and Ys = (0, -1, 0) is Xs'.
    pub fn to_internal(&self, pm: PhysMove) -> (usize, Frame) {
        let f = face_of_normal(mat_vec(&transpose(&self.0), FACES[pm.face]));
        let dir = usize::from(!pm.cw);
        if f >= 4 {
            let kind = match pm.layer {
                Layer::Face => 0,
                Layer::Wide => 1,
                Layer::Slice => 2,
            };
            (2 * (f - 4 + 4 * kind) + dir, *self)
        } else {
            let g = opposite(f) - 4;
            let rotated = Frame(mat_mul(&self.0, &face_rotation(f, pm.cw)));
            match pm.layer {
                Layer::Face => (2 * (g + 4) + dir, rotated),
                Layer::Wide => (2 * g + dir, rotated),
                Layer::Slice => (2 * (g + 8) + (1 - dir), *self),
            }
        }
    }

    /// Physical move equivalent to an internal move, and the new frame. Wide internal moves
    /// become front-face turns plus a rotation, or with `wide` the back-face wide turns
    /// themselves (then the frame never changes). Slices are named by their front face.
    pub fn to_physical(&self, m: usize, wide: bool) -> (PhysMove, Frame) {
        let s = slot(m);
        let cw = move_cw(m);
        let x = slot_face(s);
        if slot_is_slice(s) {
            let p = self.physical_face(x);
            let pm = if p >= 4 {
                PhysMove {
                    face: opposite(p),
                    layer: Layer::Slice,
                    cw: !cw,
                }
            } else {
                PhysMove {
                    face: p,
                    layer: Layer::Slice,
                    cw,
                }
            };
            (pm, *self)
        } else if slot_is_wide(s) && wide {
            let pm = PhysMove {
                face: self.physical_face(x),
                layer: Layer::Wide,
                cw,
            };
            (pm, *self)
        } else if slot_is_wide(s) {
            let f = opposite(x);
            let pm = PhysMove::face(self.physical_face(f), cw);
            (pm, Frame(mat_mul(&self.0, &face_rotation(f, cw))))
        } else {
            (PhysMove::face(self.physical_face(x), cw), *self)
        }
    }
}

/// Internal position reached by a physical scramble, and the frame afterwards.
pub fn scramble_to_cube(moves: &[PhysMove]) -> (Cube, Frame) {
    let mut fr = Frame::identity();
    let mut c = Cube::solved();
    for &pm in moves {
        let (m, nf) = fr.to_internal(pm);
        c = c.apply(m);
        fr = nf;
    }
    (c, fr)
}

/// Convert an internal solution, starting in frame `fr`, to physical moves (with `wide`, wide
/// internal moves are printed as wide turns: see [`Frame::to_physical`]).
pub fn solution_to_physical(fr: Frame, sol: &[usize], wide: bool) -> Vec<PhysMove> {
    let mut fr = fr;
    let mut out = vec![];
    for &m in sol {
        let (pm, nf) = fr.to_physical(m, wide);
        out.push(pm);
        fr = nf;
    }
    out
}

/// Apply physical moves to a physical sticker state (used for verification).
pub fn apply_physical(st: &Stickers, moves: &[PhysMove]) -> Stickers {
    let sm = &defs().sm;
    let mut s = *st;
    for pm in moves {
        let perm = match pm.layer {
            Layer::Face => sm.face_turn(pm.face, pm.cw),
            Layer::Wide => sm.wide_turn(pm.face, pm.cw),
            Layer::Slice => sm.slice_turn(pm.face, pm.cw),
        };
        s = permute_stickers(&s, &perm);
    }
    s
}

/// True if the sticker state is solved up to a whole-puzzle rotation.
pub fn is_solved_physical(st: &Stickers) -> bool {
    // every face monochrome <=> solved up to a whole-puzzle rotation
    let d = defs();
    (0..8).all(|f| {
        let on_face: Vec<usize> = (0..72).filter(|&i| d.sm.st[i].face == f).collect();
        on_face.iter().all(|&i| st[i] == st[on_face[0]])
    })
}

/// The inverse of a move sequence: the moves reversed, each turned the other way. The position
/// that a sequence solves is its inverse applied to the solved puzzle.
pub fn invert_moves(moves: &[PhysMove]) -> Vec<PhysMove> {
    moves
        .iter()
        .rev()
        .map(|m| PhysMove { cw: !m.cw, ..*m })
        .collect()
}

/// The labelled sticker states of the whole-puzzle rotations a solution can end in (sticker i
/// of the solved puzzle labelled i), found by closing `Uv`, `Rv` and `Fv`: the 12 rotations
/// with an even corner permutation (face turns 3-cycle corners).
fn rotation_labels() -> &'static Vec<Stickers> {
    static R: std::sync::OnceLock<Vec<Stickers>> = std::sync::OnceLock::new();
    R.get_or_init(|| {
        let id: Stickers = std::array::from_fn(|i| i as u8);
        let gens: Vec<Vec<PhysMove>> = ["Uv", "Rv", "Fv"]
            .iter()
            .filter_map(|g| parse_moves(g).ok())
            .collect();
        let mut seen = vec![id];
        let mut i = 0;
        while i < seen.len() {
            for g in &gens {
                let n = apply_physical(&seen[i], g);
                if !seen.contains(&n) {
                    seen.push(n);
                }
            }
            i += 1;
        }
        seen
    })
}

/// True if `solution` solves the position of `scramble` as a pure FTO, where every piece and
/// every centre triangle is unique (not only up to colour), up to a whole-puzzle rotation.
pub fn solves_pure(scramble: &[PhysMove], solution: &[PhysMove]) -> bool {
    let end = apply_physical(&pure_stickers(scramble), solution);
    rotation_labels().contains(&end)
}

/// The pure FTO after `moves` from solved: which sticker of the solved puzzle is in each place.
pub fn pure_stickers(moves: &[PhysMove]) -> Stickers {
    let id: Stickers = std::array::from_fn(|i| i as u8);
    apply_physical(&id, moves)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{Rng, SeedableRng};

    #[test]
    fn inverse_undoes() {
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(9);
        let solved = solved_stickers(defs());
        for _ in 0..20 {
            let m = random_moves(&mut rng, 25);
            let st = apply_physical(&apply_physical(&solved, &m), &invert_moves(&m));
            assert_eq!(st, solved);
        }
        // a rotation token inverts too
        let m = parse_moves("U F' Lv").unwrap();
        let st = apply_physical(&apply_physical(&solved, &m), &invert_moves(&m));
        assert_eq!(st, solved);
    }

    #[test]
    fn pure_check() {
        // the 12 rotations about face axes; the other 12 of the octahedron (quarter turns about
        // a vertex) permute the corners oddly, which no sequence of face turns can do
        assert_eq!(rotation_labels().len(), 12);
        let scr = parse_moves("U R' BL F D").unwrap();
        let inv: Vec<PhysMove> = scr
            .iter()
            .rev()
            .map(|m| PhysMove { cw: !m.cw, ..*m })
            .collect();
        assert!(solves_pure(&scr, &inv));
        let mut rot = inv.clone();
        rot.extend(parse_moves("Lv").unwrap());
        assert!(solves_pure(&scr, &rot));
        // one move short: not solved at all
        assert!(!solves_pure(&scr, &inv[..inv.len() - 1]));
    }

    /// Random physical moves of every layer kind.
    fn random_moves(rng: &mut impl Rng, n: usize) -> Vec<PhysMove> {
        (0..n)
            .map(|_| PhysMove {
                face: rng.gen_range(0..8),
                layer: [Layer::Face, Layer::Wide, Layer::Slice][rng.gen_range(0..3)],
                cw: rng.gen(),
            })
            .collect()
    }

    #[test]
    fn parse_and_print() {
        let m = parse_moves("U R' BL F2 d br'").unwrap();
        assert_eq!(format_moves(&m), "U R' BL F F D BR'");
        assert!(parse_moves("X").is_err());
        assert!(parse_moves("Xs").is_err());
        assert!(parse_moves("Uq").is_err());
        let m = parse_moves("Uw Rw' BRw2 Us' bls Ds2' FS").unwrap();
        assert_eq!(format_moves(&m), "Uw Rw' BRw BRw Us' BLs Ds' Ds' Fs");
        // every move prints as text that parses back to it
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(7);
        let sc = random_moves(&mut rng, 200);
        assert_eq!(parse_moves(&format_moves(&sc)).unwrap(), sc);
    }

    /// `Us` = `Ds'` etc.: the two names of a slice are the same move.
    #[test]
    fn slice_names() {
        let st = solved_stickers(defs());
        for f in 0..8 {
            for cw in [true, false] {
                let a = PhysMove {
                    face: f,
                    layer: Layer::Slice,
                    cw,
                };
                let b = PhysMove {
                    face: opposite(f),
                    layer: Layer::Slice,
                    cw: !cw,
                };
                assert_eq!(apply_physical(&st, &[a]), apply_physical(&st, &[b]));
                assert_eq!(scramble_to_cube(&[a]), scramble_to_cube(&[b]));
            }
        }
    }

    #[test]
    fn rotations_turn_the_whole_puzzle() {
        let solved = solved_stickers(defs());
        for name in FACE_NAMES {
            for dir in ["", "'"] {
                let one = parse_moves(&format!("{name}v{dir}")).unwrap();
                let st = apply_physical(&solved, &one);
                // a rotation: still solved up to a rotation, but not the identity
                assert!(is_solved_physical(&st), "{name}v{dir}");
                assert_ne!(st, solved, "{name}v{dir}");
                // three of them are the identity, sticker for sticker
                let three =
                    parse_moves(&format!("{name}v{dir} {name}v{dir} {name}v{dir}")).unwrap();
                assert_eq!(apply_physical(&solved, &three), solved, "{name}v{dir} x3");
            }
        }
        // a rotation then a turn is the turn of the face that came to that place
        let a = apply_physical(&solved, &parse_moves("Uv R").unwrap());
        let (c, _) = scramble_to_cube(&parse_moves("Uv R").unwrap());
        let (d, _) = scramble_to_cube(&parse_moves("Uv R Uv'").unwrap());
        assert!(!is_solved_physical(&a));
        assert_eq!(
            c, d,
            "a trailing rotation changes nothing up to orientation"
        );
    }

    #[test]
    fn scramble_matches_physical_simulation() {
        let d = defs();
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(5);
        for _ in 0..100 {
            let sc = random_moves(&mut rng, 30);
            let phys = apply_physical(&solved_stickers(d), &sc);
            let (c, fr) = scramble_to_cube(&sc);
            // rotate the internal sticker state into the physical frame
            let st = c.to_stickers();
            let perm = d.sm.whole(&fr.0);
            let mut rotated = [0u8; 72];
            for i in 0..72 {
                rotated[perm[i]] = st[i]; // rotating the puzzle moves stickers, colours stay
            }
            assert_eq!(rotated, phys);
        }
    }

    /// The user's pure corner 6-flip position: the internal state must match the stickers.
    #[test]
    fn corner_six_flip_matches_physical_simulation() {
        let d = defs();
        let sc = parse_moves(
            "U F' U' R' L' R L U F U' F' L' R' L R BL L F' L U F' R' U R' BR' U B BR' B BL BR' D' BL D'",
        )
        .unwrap();
        let phys = apply_physical(&solved_stickers(d), &sc);
        let (c, fr) = scramble_to_cube(&sc);
        let st = c.to_stickers();
        let perm = d.sm.whole(&fr.0);
        let mut rotated = [0u8; 72];
        for i in 0..72 {
            rotated[perm[i]] = st[i];
        }
        assert_eq!(rotated, phys);
    }

    /// Internal -> physical -> internal gives the same move and frame, in every frame.
    #[test]
    fn physical_round_trip() {
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(8);
        for _ in 0..30 {
            let (_, fr) = scramble_to_cube(&random_moves(&mut rng, 20));
            for m in 0..MAX_MOVES {
                for wide in [false, true] {
                    let (pm, nf) = fr.to_physical(m, wide);
                    assert_eq!(fr.to_internal(pm), (m, nf));
                    match pm.layer {
                        Layer::Slice => assert!(pm.face < 4, "slices are named by a front face"),
                        Layer::Wide => assert!(wide && nf == fr),
                        Layer::Face => {}
                    }
                }
            }
        }
    }

    #[test]
    fn inverse_scramble_solves() {
        let d = defs();
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(6);
        for _ in 0..50 {
            let sc = random_moves(&mut rng, 25);
            let (c, fr) = scramble_to_cube(&sc);
            // internal inverse of the internal scramble
            let mut f2 = Frame::identity();
            let mut internal = vec![];
            for &pm in &sc {
                let (m, nf) = f2.to_internal(pm);
                internal.push(m);
                f2 = nf;
            }
            let inv: Vec<usize> = internal.iter().rev().map(|&m| inverse_move(m)).collect();
            assert_eq!(c.apply_seq(&inv), Cube::solved());
            for wide in [false, true] {
                let phys = solution_to_physical(fr, &inv, wide);
                let st = apply_physical(&apply_physical(&solved_stickers(d), &sc), &phys);
                assert!(is_solved_physical(&st));
            }
        }
    }
}
