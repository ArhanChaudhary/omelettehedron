//! What can be said about a position without pruning tables (`omelettehedron show`): whether
//! it is solved, its symmetries, how many pieces are solved, whether it is in the two-phase
//! solver's subgroups `H_A` and `H_B`, and a lower bound from small distance tables built in
//! memory (corners, the centres of each orbit and the edges' A-faces taken one at a time; tens
//! of milliseconds).

use crate::coord::{move_tables, Coords, N_CORNER, N_MSET};
use crate::cube::{defs, Cube, Stickers, N_MOVES, SWAP_AB};
use crate::geometry::{det, rotation_group};
use crate::optimal::sym::{swaps_classes, sym48, N_SYM};
use crate::tables::Phase2;
use serde::Serialize;
use std::sync::OnceLock;

/// The symmetries of the octahedron that fix a position: `M P M^-1 = P` up to a whole-puzzle
/// rotation (turned or mirrored by `M`, each colour renamed after the face it moves to: the
/// same position). The optimal solver searches such positions with a symmetric start.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub struct Symmetries {
    /// how many of the 48 (at least 1: the identity)
    pub count: usize,
    /// of those, rotations (the identity included) and reflections
    pub rotations: usize,
    pub reflections: usize,
    /// of those, how many exchange the face classes (`U F BR BL` and `R L B D`)
    pub class_swapping: usize,
}

/// The pieces, with the puzzle held the way (of the 12 that turns can reach) that puts the most
/// stickers on a face of their colour.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub struct Pieces {
    /// of 6
    pub corners_solved: usize,
    /// in their place but turned over (not counted as solved)
    pub corners_turned_over: usize,
    /// of 12 (an edge in its place is solved: edges have no orientation)
    pub edges_solved: usize,
    /// centres on a face of their colour, of 12 in each orbit (faces `U F BR BL`, and `R L B D`)
    pub a_centres_solved: usize,
    pub b_centres_solved: usize,
}

/// The centres as a pure FTO, where every centre is a piece of its own (known only for a
/// position given by moves): of the centres on a face of their colour, in the same hold as
/// [`Pieces`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub struct Pure {
    /// in their own place
    pub own_place: usize,
    /// in the place of another centre of the same colour
    pub exchanged: usize,
}

/// Lower bounds (face turns) from the distances to solved of parts of the position, each part
/// alone (exact distances of small projections, built in memory).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub struct Bounds {
    /// the largest of the others
    pub moves: u8,
    pub corners: u8,
    pub a_centres: u8,
    pub b_centres: u8,
    /// which A-face's edge sits in each edge place
    pub edges: u8,
}

/// What `show` reports about a position.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub struct Info {
    pub solved: bool,
    pub symmetries: Symmetries,
    pub pieces: Pieces,
    pub pure: Option<Pure>,
    /// moves to solve the position with `U F BR BL` alone, if it is in `H_A` (phase 2 of the
    /// two-phase solver: the exact table of `H_A`)
    pub h_a: Option<u8>,
    /// the same with `R L B D` alone, if it is in `H_B`
    pub h_b: Option<u8>,
    pub lower_bound: Bounds,
}

/// The facts about `cube` (a position in the frame of [`crate::cube`]); `labels` is the sticker
/// state of the pure puzzle that made it ([`crate::notation::pure_stickers`]), if known.
pub fn info(cube: &Cube, labels: Option<&Stickers>) -> Info {
    let p2 = phase2();
    let (pieces, pure) = pieces(cube, labels);
    Info {
        solved: *cube == Cube::solved(),
        symmetries: symmetries(cube),
        pieces,
        pure,
        h_a: p2.length(cube),
        h_b: cube.conjugate(&SWAP_AB).and_then(|c| p2.length(&c)),
        lower_bound: bounds(cube),
    }
}

fn phase2() -> &'static Phase2 {
    static P: OnceLock<Phase2> = OnceLock::new();
    P.get_or_init(Phase2::new)
}

/// The symmetries fixing `cube`, as the optimal solver finds them.
pub fn symmetries(cube: &Cube) -> Symmetries {
    let s = sym48();
    let fixing: Vec<usize> = (0..N_SYM)
        .filter(|&k| k == 0 || s.conj_by_maps(cube, k) == *cube)
        .collect();
    let rotations = fixing.iter().filter(|&&k| det(&s.mats[k]) == 1).count();
    Symmetries {
        count: fixing.len(),
        rotations,
        reflections: fixing.len() - rotations,
        class_swapping: fixing.iter().filter(|&&k| swaps_classes(k)).count(),
    }
}

/// The pieces in the hold that puts the most stickers on a face of their colour, and the pure
/// centres in that hold (if `labels` is given). `labels` also gives the colours if given (a
/// position made by moves may be held differently from the frame of `cube`; the counts do not
/// depend on the frame).
fn pieces(cube: &Cube, labels: Option<&Stickers>) -> (Pieces, Option<Pure>) {
    let d = defs();
    let colour: Stickers = match labels {
        Some(l) => std::array::from_fn(|i| d.sm.st[usize::from(l[i])].face as u8),
        None => cube.to_stickers(),
    };
    let mut best: Option<((usize, usize), Pieces, Option<Pure>)> = None;
    for r in rotation_group() {
        // the solved puzzle held turned by `r`: sticker i of the solved puzzle is in place
        // `whole[i]`
        let whole = d.sm.whole(&r);
        let mut target = [0u8; 72];
        let mut home = [0u8; 72];
        for (i, &w) in whole.iter().enumerate() {
            target[w] = d.sm.st[i].face as u8;
            home[w] = i as u8;
        }
        let ok = |i: usize| colour[i] == target[i];
        let corners_solved = d
            .corner_slots
            .iter()
            .filter(|s| s.iter().all(|&i| ok(i)))
            .count();
        let corners_turned_over = d
            .corner_slots
            .iter()
            .filter(|s| {
                !s.iter().all(|&i| ok(i)) && (0..4).all(|j| colour[s[j]] == target[s[(j + 2) % 4]])
            })
            .count();
        let edges_solved = d
            .edge_slots
            .iter()
            .filter(|s| s.iter().all(|&i| ok(i)))
            .count();
        let a = d.a_pos.iter().filter(|&&i| ok(i)).count();
        let b = d.b_pos.iter().filter(|&&i| ok(i)).count();
        let pure = labels.map(|l| {
            let own = d
                .a_pos
                .iter()
                .chain(&d.b_pos)
                .filter(|&&i| ok(i) && l[i] == home[i])
                .count();
            Pure {
                own_place: own,
                exchanged: a + b - own,
            }
        });
        let right = (0..72).filter(|&i| ok(i)).count();
        let key = (right, pure.map_or(0, |p| p.own_place));
        if best.as_ref().is_none_or(|b| key > b.0) {
            best = Some((
                key,
                Pieces {
                    corners_solved,
                    corners_turned_over,
                    edges_solved,
                    a_centres_solved: a,
                    b_centres_solved: b,
                },
                pure,
            ));
        }
    }
    best.map_or(
        (
            Pieces {
                corners_solved: 0,
                corners_turned_over: 0,
                edges_solved: 0,
                a_centres_solved: 0,
                b_centres_solved: 0,
            },
            None,
        ),
        |(_, p, q)| (p, q),
    )
}

/// Distances to solved of the corner state, the A- and B-centre colourings and the edges'
/// A-faces (each a projection of the position: a lower bound on its distance).
struct Distances {
    c: Vec<u8>,
    a: Vec<u8>,
    b: Vec<u8>,
    e: Vec<u8>,
}

/// Breadth-first distances from `start` over `n` states (`step(x, m)`: the state after move
/// `m`; the moves are closed under inverses, so these are the distances to `start` too).
fn bfs(n: usize, start: usize, step: impl Fn(usize, usize) -> usize) -> Vec<u8> {
    let mut dist = vec![u8::MAX; n];
    dist[start] = 0;
    let mut front = vec![start];
    let mut d = 0u8;
    while !front.is_empty() {
        let mut next = vec![];
        for &x in &front {
            for m in 0..N_MOVES {
                let y = step(x, m);
                if dist[y] == u8::MAX {
                    dist[y] = d + 1;
                    next.push(y);
                }
            }
        }
        front = next;
        d += 1;
    }
    dist
}

fn distances() -> &'static Distances {
    static D: OnceLock<Distances> = OnceLock::new();
    D.get_or_init(|| {
        let mt = move_tables();
        let s = Coords::of(&Cube::solved());
        let ((c, a), (b, e)) = rayon::join(
            || {
                rayon::join(
                    || bfs(N_CORNER, s.c.into(), |x, m| mt.c[x][m].into()),
                    || bfs(N_MSET, s.a as usize, |x, m| mt.a[x][m] as usize),
                )
            },
            || {
                rayon::join(
                    || bfs(N_MSET, s.b as usize, |x, m| mt.b[x][m] as usize),
                    || bfs(N_MSET, s.e as usize, |x, m| (mt.e[x][m] & 0xfffff) as usize),
                )
            },
        );
        Distances { c, a, b, e }
    })
}

/// The lower bounds of `cube`.
pub fn bounds(cube: &Cube) -> Bounds {
    let d = distances();
    let x = Coords::of(cube);
    let get = |t: &[u8], i: usize| t.get(i).copied().unwrap_or(0);
    let corners = get(&d.c, x.c.into());
    let a_centres = get(&d.a, x.a as usize);
    let b_centres = get(&d.b, x.b as usize);
    let edges = get(&d.e, x.e as usize);
    Bounds {
        moves: corners.max(a_centres).max(b_centres).max(edges),
        corners,
        a_centres,
        b_centres,
        edges,
    }
}

impl Info {
    /// The facts as lines of plain text.
    pub fn lines(&self) -> Vec<String> {
        let yes = |b: bool| if b { "yes" } else { "no" };
        let s = &self.symmetries;
        let sym = if s.count == 1 {
            "none (only the identity)".to_string()
        } else {
            let swap = if s.class_swapping > 0 {
                format!(
                    "; {} exchange the faces U F BR BL with R L B D",
                    s.class_swapping
                )
            } else {
                String::new()
            };
            format!(
                "{} of 48 ({} rotations, {} reflections{swap})",
                s.count, s.rotations, s.reflections
            )
        };
        let p = &self.pieces;
        let turned = match p.corners_turned_over {
            0 => String::new(),
            n => format!(", {n} more in place but turned over"),
        };
        let centres = p.a_centres_solved + p.b_centres_solved;
        let pure = match self.pure {
            None => "not known (only a scramble tells which centre is which)".to_string(),
            Some(_) if centres == 0 => "no centre is on a face of its colour".to_string(),
            Some(q) if q.exchanged == 0 => {
                format!("each of the {centres} is its own piece (none exchanged)")
            }
            Some(q) => format!(
                "{} of the {centres} are their own piece, {} exchanged with a centre of the same \
                 colour",
                q.own_place, q.exchanged
            ),
        };
        let moves = |n: u8| format!("{n} move{}", if n == 1 { "" } else { "s" });
        let sub = |name: &str, faces: &str, d: Option<u8>| {
            let answer = d.map_or("no".to_string(), |n| format!("yes, in {}", moves(n)));
            format!("in {name} (solvable with {faces} alone): {answer}")
        };
        let b = &self.lower_bound;
        vec![
            format!("solved: {}", yes(self.solved)),
            format!("symmetries: {sym}"),
            "pieces (with the puzzle held so that most stickers are right):".to_string(),
            format!("  corners: {} of 6 solved{turned}", p.corners_solved),
            format!("  edges: {} of 12 solved", p.edges_solved),
            format!(
                "  centres: {centres} of 24 on a face of their colour ({} of 12 on U F BR BL, \
                 {} of 12 on R L B D)",
                p.a_centres_solved, p.b_centres_solved
            ),
            format!("  pure FTO (every centre unique): {pure}"),
            sub("H_A", "U F BR BL", self.h_a),
            sub("H_B", "R L B D", self.h_b),
            format!(
                "lower bound: {} (the corners alone need {}, the centres {} and {}, the edges {})",
                moves(b.moves),
                b.corners,
                b.a_centres,
                b.b_centres,
                b.edges
            ),
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::notation::{parse_moves, pure_stickers, scramble_to_cube};
    use rand::SeedableRng;

    fn labels_of(text: &str) -> (Cube, Stickers) {
        let m = parse_moves(text).unwrap();
        (scramble_to_cube(&m).0, pure_stickers(&m))
    }

    #[test]
    fn solved_position() {
        let (c, l) = labels_of("");
        let i = info(&c, Some(&l));
        assert!(i.solved);
        assert_eq!(i.symmetries.count, 48);
        assert_eq!(i.symmetries.rotations, 24);
        assert_eq!(i.symmetries.class_swapping, 24);
        assert_eq!(
            i.pieces,
            Pieces {
                corners_solved: 6,
                corners_turned_over: 0,
                edges_solved: 12,
                a_centres_solved: 12,
                b_centres_solved: 12
            }
        );
        assert_eq!(
            i.pure,
            Some(Pure {
                own_place: 24,
                exchanged: 0
            })
        );
        assert_eq!((i.h_a, i.h_b), (Some(0), Some(0)));
        assert_eq!(i.lower_bound.moves, 0);
    }

    #[test]
    fn one_turn() {
        // one turn, of a class-A face (U, BL) or a class-B one (R, D): in H_A or in H_B, one
        // move from solved. It moves 3 corners, 3 edges and 9 centres, 3 of them (the face's
        // own) to places of their colour, whatever the frame.
        for text in ["U", "R'", "BL", "D"] {
            let (c, l) = labels_of(text);
            let i = info(&c, Some(&l));
            assert!(!i.solved);
            let a = matches!(text, "U" | "BL");
            assert_eq!(i.h_a, a.then_some(1), "{text}");
            assert_eq!(i.h_b, (!a).then_some(1), "{text}");
            assert_eq!(i.pieces.corners_solved, 3, "{text}");
            assert_eq!(i.pieces.corners_turned_over, 0, "{text}");
            assert_eq!(i.pieces.edges_solved, 9, "{text}");
            assert_eq!(
                i.pieces.a_centres_solved + i.pieces.b_centres_solved,
                18,
                "{text}"
            );
            assert_eq!(
                i.pure,
                Some(Pure {
                    own_place: 15,
                    exchanged: 3
                }),
                "{text}"
            );
            assert_eq!(i.lower_bound.moves, 1, "{text}");
            // the rotations about the turned face's axis
            assert!(i.symmetries.rotations >= 3, "{text}: {:?}", i.symmetries);
        }
    }

    /// Three centres of one colour cycled: solved as a normal FTO, three centres exchanged as a
    /// pure one.
    #[test]
    fn exchanged_centres() {
        let d = defs();
        let on_u: Vec<usize> = d
            .a_pos
            .iter()
            .copied()
            .filter(|&i| d.sm.st[i].face == 0)
            .collect();
        assert_eq!(on_u.len(), 3);
        let mut l: Stickers = std::array::from_fn(|i| i as u8);
        l.swap(on_u[0], on_u[1]);
        l.swap(on_u[1], on_u[2]);
        let i = info(&Cube::solved(), Some(&l));
        assert!(i.solved);
        assert_eq!(
            i.pure,
            Some(Pure {
                own_place: 21,
                exchanged: 3
            })
        );
        assert!(i.lines().iter().any(|x| x.contains("3 exchanged")));
    }

    #[test]
    fn symmetries_match_the_optimal_solver_and_colour_renaming() {
        // the count with each colour renamed after the face it moves to equals the count with
        // the colours renamed in any way (the pattern looks the same up to colour names)
        let sm = &defs().sm;
        let mats = &sym48().mats;
        let any_renaming = |c: &Cube| -> usize {
            let st = c.to_stickers();
            mats.iter()
                .filter(|m| {
                    let w = sm.whole(m);
                    let mut map = [u8::MAX; 8];
                    (0..72).all(|i| {
                        let (from, to) = (usize::from(st[i]), st[w[i]]);
                        if map[from] == u8::MAX {
                            map[from] = to;
                        }
                        map[from] == to
                    })
                })
                .count()
        };
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(3);
        let mut cubes: Vec<Cube> = (0..20).map(|_| Cube::random(&mut rng)).collect();
        for text in ["U", "U R", "U BL' D", "R L", "U F BR BL", "U2 B2"] {
            cubes.push(labels_of(text).0);
        }
        for line in include_str!("../bench/hard.tsv").lines() {
            let cols: Vec<&str> = line.split('\t').collect();
            if line.starts_with('#') || cols.len() < 4 {
                continue;
            }
            cubes.push(labels_of(cols[3]).0);
        }
        let mut symmetric = 0;
        for c in &cubes {
            let s = symmetries(c);
            assert_eq!(s.count, any_renaming(c));
            assert_eq!(s.rotations + s.reflections, s.count);
            if s.count > 1 {
                symmetric += 1;
            }
        }
        assert!(symmetric >= 6, "{symmetric}");
        // a random position has no symmetry
        assert_eq!(symmetries(&cubes[0]).count, 1);
    }

    #[test]
    fn bounds_are_lower_bounds() {
        // scrambles of n moves are at most n moves from solved
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(4);
        for n in 0..12 {
            let m: Vec<_> = (0..n)
                .map(|_| {
                    use rand::Rng;
                    crate::notation::PhysMove::face(rng.gen_range(0..8), rng.gen_bool(0.5))
                })
                .collect();
            let (c, _) = scramble_to_cube(&m);
            let b = bounds(&c);
            assert!(usize::from(b.moves) <= n, "{n}: {b:?}");
        }
        // a random position is several moves out
        let c = Cube::random(&mut rng);
        assert!(bounds(&c).moves >= 6, "{:?}", bounds(&c));
    }
}
