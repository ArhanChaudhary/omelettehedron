//! Physical move notation and conversion to/from internal (fixed-corner) moves.

use crate::cube::*;
use crate::geometry::*;
use std::fmt;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PhysMove {
    pub face: usize,
    pub cw: bool,
}

impl fmt::Display for PhysMove {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}{}",
            FACE_NAMES[self.face],
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

/// Parse a scramble such as `U R' BL F2 D`. `X2` is a double turn (= two moves), `X'`/`X3`
/// variants of the prime are accepted.
///
/// # Errors
///
/// If a token is not a face name with an optional `'`, `2`, `2'` or `3` suffix.
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
            out.push(PhysMove { face, cw });
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

    /// Internal move equivalent to a physical move, and the new frame.
    pub fn to_internal(&self, pm: PhysMove) -> (usize, Frame) {
        let f = face_of_normal(mat_vec(&transpose(&self.0), FACES[pm.face]));
        let dir = usize::from(!pm.cw);
        if f >= 4 {
            (2 * (f - 4) + dir, *self)
        } else {
            let g = opposite(f);
            let m = 2 * (g - 4 + 4) + dir;
            (m, Frame(mat_mul(&self.0, &face_rotation(f, pm.cw))))
        }
    }

    /// Physical move equivalent to an internal move, and the new frame.
    pub fn to_physical(&self, m: usize) -> (PhysMove, Frame) {
        let s = slot(m);
        let cw = move_cw(m);
        if slot_is_wide(s) {
            let f = opposite(slot_face(s));
            let pm = PhysMove {
                face: face_of_normal(mat_vec(&self.0, FACES[f])),
                cw,
            };
            (pm, Frame(mat_mul(&self.0, &face_rotation(f, cw))))
        } else {
            let f = slot_face(s);
            (
                PhysMove {
                    face: face_of_normal(mat_vec(&self.0, FACES[f])),
                    cw,
                },
                *self,
            )
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

/// Convert an internal solution, starting in frame `fr`, to physical moves.
pub fn solution_to_physical(fr: Frame, sol: &[usize]) -> Vec<PhysMove> {
    let mut fr = fr;
    let mut out = vec![];
    for &m in sol {
        let (pm, nf) = fr.to_physical(m);
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
        s = permute_stickers(&s, &sm.face_turn(pm.face, pm.cw));
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

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{Rng, SeedableRng};

    fn random_moves(rng: &mut impl Rng, n: usize) -> Vec<PhysMove> {
        (0..n)
            .map(|_| PhysMove {
                face: rng.gen_range(0..8),
                cw: rng.gen(),
            })
            .collect()
    }

    #[test]
    fn parse_and_print() {
        let m = parse_moves("U R' BL F2 d br'").unwrap();
        assert_eq!(format_moves(&m), "U R' BL F F D BR'");
        assert!(parse_moves("X").is_err());
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
            let phys = solution_to_physical(fr, &inv);
            let st = apply_physical(&apply_physical(&solved_stickers(d), &sc), &phys);
            assert!(is_solved_physical(&st));
        }
    }
}
