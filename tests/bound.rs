//! Checks behind docs/BOUND.md: the relaxed set `H'` (A-centres and edges solved up to block
//! rotation and whole-puzzle rotation), the phase-2 group `H_A` inside it, and why a coset
//! ("block") method that parametrises 933,120-element blocks through `H_A` is not valid for
//! patterns.
//!
//! Parts:
//! * structure of `H_A` and `H'` (no tables): `H_A` acts faithfully on corners and edges, `H'`
//!   is the 972 A-centre/edge configurations of `H_A` and equals "A-blocks intact", only the
//!   solved corners give uniform B faces;
//! * the counterexample (no tables): pure 3-cycles of the B-centres of one face, and the
//!   solved block, where a walk does not act as a right translation;
//! * with the medium-tier tables (`#[ignore]`, `FTO_TABLES`): the A-centre/edge table against
//!   an independent brute force, the corner/B-centre table against an independent breadth-first
//!   search, the walk search against an unpruned one, real phase-1 walks of random positions of
//!   `H'`, a sample of per-block bounds (a measurement for the cost estimate), and a position of
//!   `H'` 14 moves from `H_A` (also checked without either table).

use omelettehedron::coord::*;
use omelettehedron::cube::*;
use omelettehedron::fast::*;
use omelettehedron::geometry::*;
use omelettehedron::notation::*;
use omelettehedron::prune::decode;
use omelettehedron::sym::*;
use omelettehedron::tables::*;
use rand::{Rng, SeedableRng};
use std::collections::{HashMap, HashSet};
use std::sync::OnceLock;

/// A raw A-centre/edge configuration: A-colouring, edge triple labelling, chirality bits.
type AeKey = ([u8; 12], [u8; 12], u8);

/// The solved position with every centre labelled by its position (`xa`, `xb` = 0..12), so that
/// identical centres can be told apart. Moves act on labels exactly as on colours.
fn labelled_solved() -> Cube {
    let id: [u8; 12] = std::array::from_fn(|i| u8::try_from(i).unwrap());
    Cube {
        xa: id,
        xb: id,
        ..Cube::solved()
    }
}

/// The pattern (colours only) of a labelled position: label `i` has the colour of position `i`.
fn pattern(c: &Cube) -> Cube {
    Cube {
        xa: c.xa.map(|l| l / 3),
        xb: c.xb.map(|l| l / 3),
        ..*c
    }
}

/// The raw A-centre/edge configuration of a position, without ranks or symmetries.
fn ae_key(c: &Cube) -> AeKey {
    (c.xa, c.etri(), edge_chirality(&c.ep))
}

/// All elements of `H_A` with labelled centres, by breadth-first search with the A-moves.
fn labelled_h_a() -> &'static [Cube] {
    static H: OnceLock<Vec<Cube>> = OnceLock::new();
    H.get_or_init(|| {
        let amoves: Vec<usize> = (0..N_MOVES).filter(|&m| is_a_move(m)).collect();
        let mut seen = HashSet::from([labelled_solved()]);
        let mut all = vec![labelled_solved()];
        let mut i = 0;
        while i < all.len() {
            let c = all[i];
            i += 1;
            for &m in &amoves {
                let n = c.apply(m);
                if seen.insert(n) {
                    all.push(n);
                }
            }
        }
        all
    })
}

/// The A-centre/edge configurations (edge arrangement, A-colouring) of `H_A`.
fn h_a_configs() -> HashSet<([u8; 12], [u8; 12])> {
    labelled_h_a()
        .iter()
        .map(|c| (c.ep, pattern(c).xa))
        .collect()
}

fn inverse_word(w: &[usize]) -> Vec<usize> {
    w.iter().rev().map(|&m| inverse_move(m)).collect()
}

/// `H_A` acts faithfully on corners and edges: with labelled centres it still has 933,120
/// elements, one per pattern, and its centres follow (A-centres from the edges, B-centres from
/// the corners, as `Phase2` and `CbTable::attach_targets` assume).
#[test]
fn h_a_is_faithful_on_corners_and_edges() {
    let h = labelled_h_a();
    assert_eq!(h.len(), 933_120);
    let p2 = Phase2::new();
    let patterns: HashSet<Cube> = h.iter().map(pattern).collect();
    assert_eq!(patterns.len(), 933_120);
    assert!(patterns.iter().all(|c| p2.index(c).is_some()));
    let corners_edges: HashSet<_> = h.iter().map(|c| (c.cp, c.co, c.ep)).collect();
    assert_eq!(corners_edges.len(), 933_120);
    assert_eq!(h_a_configs().len(), 972);
    let attach = CbTable::attach_targets();
    for c in &patterns {
        let id = usize::from(corners().id(c));
        assert_eq!(mset().rank(&c.xb), attach[id]);
    }
}

/// The relaxed set `H'` (AE distance 0) is exactly "A-blocks intact": the 972 A-centre/edge
/// configurations of `H_A` are the 12 rotations of the solved puzzle times the 81 rotations of
/// the four A-face hexagons (three A-centres and three edges each), they form 12 raw AE values
/// of 81 edge arrangements each, the 288 symmetries of the AE table keep the set, and the table
/// index maps all of them to its target.
#[test]
fn relaxed_set_is_a_blocks_intact() {
    let configs = h_a_configs();
    let d = defs();
    let sm = &d.sm;
    let solved = solved_stickers(d);
    let mut built = HashSet::new();
    for r in rotation_group() {
        let rotated = permute_stickers(&solved, &sm.whole(&r));
        for k in 0..81usize {
            let mut st = rotated;
            let mut rest = k;
            for &f in &A_FACES {
                let turns = rest % 3;
                rest /= 3;
                // the hexagon of face f: its centres and both stickers of its edges
                let in_hexagon = |p: V3| -> bool {
                    let i = sm.index_of(p);
                    match sm.st[i].kind {
                        Kind::Centre => sm.st[i].face == f,
                        Kind::Edge => d
                            .edge_slots
                            .iter()
                            .any(|e| e.contains(&i) && sm.st[e[0]].face == f),
                        Kind::Corner => false,
                    }
                };
                let turn = sm.perm_where(&face_rotation(f, true), in_hexagon);
                for _ in 0..turns {
                    st = permute_stickers(&st, &turn);
                }
            }
            let c = Cube::from_stickers(&st).unwrap();
            built.insert((c.ep, c.xa));
        }
    }
    assert_eq!(built, configs);
    let mut raw: HashMap<(u32, u32, u8), usize> = HashMap::new();
    for (ep, xa) in &configs {
        let x = Coords::of(&Cube {
            ep: *ep,
            xa: *xa,
            ..Cube::solved()
        });
        *raw.entry((x.a, x.e, x.ebits)).or_default() += 1;
    }
    assert_eq!(raw.len(), 12);
    assert!(raw.values().all(|&n| n == 81));
    let (s, es, ms) = (syms(), edge_syms(), mset());
    for &(a, e, bits) in raw.keys() {
        let l = ms.unrank(a);
        for &(g, t) in &s.pairs {
            let na = ms.rank(&s.apply_a(g.into(), t.into(), &l));
            let (ne, nb) = es.apply_e(g.into(), t.into(), e, bits);
            assert!(raw.contains_key(&(na, ne, nb & 15)));
        }
    }
    let ix = AeSym::new(true);
    let targets = ae_targets(&ix);
    for &(a, e, ebits) in raw.keys() {
        let x = Coords {
            c: 0,
            a,
            b: 0,
            e,
            ebits,
        };
        assert!(targets.contains(&ix.index(&x)));
    }
}

/// Only the solved corner state of `H_A` shows one colour on the three B-centres of every B
/// face (the B-centres of `H_A` copy the adjacent corner stickers).
#[test]
fn only_solved_corners_give_uniform_b_faces() {
    let ms = mset();
    let uniform: Vec<usize> = CbTable::attach_targets()
        .iter()
        .enumerate()
        .filter(|&(_, &b)| {
            let l = ms.unrank(b);
            (0..4).all(|f| l[3 * f] == l[3 * f + 1] && l[3 * f] == l[3 * f + 2])
        })
        .map(|(c, _)| c)
        .collect();
    assert_eq!(uniform, [usize::from(corners().id(&Cube::solved()))]);
}

/// Pure 3-cycles of the three centres of each B face (internal moves; slices written as a wide
/// turn and the inverse single turn): a setup, the commutator `[BR, BL BRs BL']`, the setup
/// undone. Found by search over setups of at most four moves; each starts with a move outside
/// `H_A`, so each word ends with one.
const FACE_CYCLES: [&[usize]; 4] = [
    &[0, 9, 3, 2, 4, 10, 3, 5, 3, 4, 11, 2, 5, 2, 8, 1],
    &[12, 10, 9, 3, 2, 4, 10, 3, 5, 3, 4, 11, 2, 5, 2, 8, 11, 13],
    &[11, 9, 3, 2, 4, 10, 3, 5, 3, 4, 11, 2, 5, 2, 8, 10],
    &[10, 9, 3, 2, 4, 10, 3, 5, 3, 4, 11, 2, 5, 2, 8, 11],
];

/// `FACE_CYCLES[0]` in face turns: a pure 3-cycle of the three centres of face R.
const W_R: &str = "B F' BR D U' BL' U' BL D' U BR' U F B'";

/// Each word of [`FACE_CYCLES`] leaves the pattern solved and moves exactly the three
/// B-centres of its face: the solved position has walks that permute identical B-centres.
#[test]
fn b_faces_have_pure_centre_cycles() {
    let s = labelled_solved();
    for (f, w) in FACE_CYCLES.iter().enumerate() {
        let c = s.apply_seq(w);
        assert_eq!(pattern(&c), Cube::solved());
        let moved: Vec<usize> = (0..12).filter(|&p| usize::from(c.xb[p]) != p).collect();
        assert_eq!(moved, [3 * f, 3 * f + 1, 3 * f + 2]);
        assert!(w.last().is_some_and(|&m| !is_a_move(m)));
        // the same on the sticker model: the pattern is solved, the pure puzzle is not
        let phys = solution_to_physical(Frame::identity(), w, false);
        println!("face {}: {}", FACE_NAMES[B_FACES[f]], format_moves(&phys));
        let st = apply_physical(&solved_stickers(defs()), &phys);
        assert!(is_solved_physical(&st));
        assert!(!solves_pure(&[], &phys));
    }
}

/// The counterexample to the block method on the block of the solved position (`s0` empty, the
/// block is the `H_A` set itself). The walk `w` = `FACE_CYCLES[0]` keeps the solved pattern, so
/// the method reads it as the right translation by the identity of `C x E_A` (`w` moves no
/// corner or edge) and credits every element `z` with the solution `w` + (phase 2 of `z`).
/// For `z` = U that is "w U'", which does not solve U: `U w` keeps the corners and edges of
/// `U` but is not in `H_A`.
#[test]
fn walk_is_not_a_translation_on_the_solved_block() {
    let p2 = Phase2::new();
    let w = FACE_CYCLES[0];
    let mut outside = 0;
    let mut total = 0;
    for z in labelled_h_a().iter().map(pattern).step_by(97) {
        let zw = z.apply_seq(w);
        assert_eq!((zw.cp, zw.co, zw.ep), (z.cp, z.co, z.ep));
        total += 1;
        if p2.index(&zw).is_none() {
            outside += 1;
        }
    }
    println!("{outside} of {total} sampled elements z of H_A have z * w outside H_A");
    assert!(outside > 0);
    // the same word in face turns, with `U' U` cancelled (14 moves)
    let wp = parse_moves(W_R).unwrap();
    assert_eq!(
        format_moves(&solution_to_physical(Frame::identity(), w, false)).replace("U' U ", ""),
        W_R
    );
    let solved = solved_stickers(defs());
    assert!(is_solved_physical(&apply_physical(&solved, &wp)));
    assert!(!solves_pure(&[], &wp));
    let u = parse_moves("U").unwrap();
    let claimed = parse_moves(&format!("{W_R} U'")).unwrap();
    assert!(!is_solved_physical(&apply_physical(
        &solved,
        &[u.clone(), claimed].concat()
    )));
    let (cu, _) = scramble_to_cube(&u);
    assert!(p2.index(&cu).is_some());
    assert!(p2.index(&cu.apply_seq(w)).is_none());
}

/// The medium-tier tables from `FTO_TABLES`.
fn medium() -> &'static Tables<AeSym> {
    static T: OnceLock<Tables<AeSym>> = OnceLock::new();
    T.get_or_init(|| {
        Tables::load(&table_dir(None), AeSym::new(true), MoveSet::FACE, false)
            .expect("medium tier tables in FTO_TABLES")
    })
}

/// Canonical sequences of `left` more moves after `last`, with the minimum of (moves so far +
/// ball distance) over the positions they reach.
fn forward(c: &Cube, last: usize, left: usize, done: usize, ball: &HashMap<AeKey, u8>) -> usize {
    let mut best = ball
        .get(&ae_key(c))
        .map_or(usize::MAX, |&d| done + usize::from(d));
    if left > 0 {
        for m in 0..N_MOVES {
            if allowed_after(false, last, m) {
                best = best.min(forward(&c.apply(m), m, left - 1, done + 1, ball));
            }
        }
    }
    best
}

/// The AE table gives the distance to `H'` as a set of patterns: on 60 positions (random ones
/// and random walks from `H'`) its exact distance equals a brute force that knows only `H'`'s
/// 972 configurations (all sequences of up to 7 moves from the position against a ball of
/// radius 5 around `H'`, exact up to 12).
#[test]
#[ignore = "needs the medium tier tables in FTO_TABLES"]
fn ae_table_matches_brute_force() {
    use rayon::prelude::*;
    const BALL: u8 = 5;
    const FWD: usize = 7;
    let ball = ae_ball(BALL);
    assert_eq!(ball.values().filter(|&&d| d == 0).count(), 12);
    let t = medium();
    println!("AE table histogram (entries per distance): {:?}", t.ae.hist);
    let h = labelled_h_a();
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(11);
    for i in 0..60 {
        let c = if i % 2 == 0 {
            Cube::random(&mut rng)
        } else {
            let a = pattern(&h[rng.gen_range(0..h.len())]);
            let mut y = Cube {
                ep: a.ep,
                xa: a.xa,
                ..Cube::random(&mut rng)
            };
            for _ in 0..rng.gen_range(1..=9) {
                y = y.apply(rng.gen_range(0..N_MOVES));
            }
            y
        };
        let table = usize::from(t.ae_exact(&Coords::of(&c)));
        let first: Vec<usize> = (0..N_MOVES).collect();
        let brute = first
            .par_iter()
            .map(|&m| forward(&c.apply(m), m, FWD - 1, 1, &ball))
            .min()
            .unwrap_or(usize::MAX)
            .min(forward(&c, MAX_MOVES, 0, 0, &ball));
        println!("position {i}: table {table}, brute force {brute}");
        if table <= FWD + usize::from(BALL) {
            assert_eq!(brute, table);
        } else {
            assert_eq!(brute, usize::MAX);
        }
    }
}

/// All phase-1 walks of exactly `depth` moves from `root`: canonical sequences that end in
/// `H_A` with a move outside `H_A` (they may pass through `H_A` on the way). Exhaustive: a
/// child is dropped only when its exact AE or CB distance exceeds the moves left. Returns the
/// walks and the number of nodes visited.
fn phase1_walks(t: &Tables<AeSym>, root: &Cube, depth: usize) -> (Vec<Vec<usize>>, u64) {
    phase1_walks_with(t, root, depth, None)
}

/// [`phase1_walks`], with the CB distances read from `raw_cb` ([`raw_cb_distances`]) instead
/// of the CB table if given.
fn phase1_walks_with(
    tables: &Tables<AeSym>,
    root: &Cube,
    depth: usize,
    raw_cb: Option<&[u8]>,
) -> (Vec<Vec<usize>>, u64) {
    use rayon::prelude::*;
    struct Walker<'a> {
        tables: &'a Tables<AeSym>,
        fast: &'static FastTables,
        raw_cb: Option<&'a [u8]>,
        root: Cube,
        depth: usize,
    }
    impl Walker<'_> {
        /// The exact CB distance of `pos`, given that of its parent.
        fn cb(&self, pos: &Fc, parent: u8) -> u8 {
            match self.raw_cb {
                Some(raw) => raw[usize::from(pos.c) * N_MSET + self.fast.b_rank(pos) as usize],
                None => decode(
                    parent,
                    self.tables.cb.tab.get(CbTable::index_fast(self.fast, pos)),
                ),
            }
        }
        /// The child of `pos` (distances `ae`, `cb`) by move `m`, with its exact distances.
        fn child(&self, pos: &Fc, masks: [u16; 4], ae: u8, cb: u8, m: usize) -> (Fc, u8, u8) {
            let (fast, tables) = (self.fast, self.tables);
            let next = fast.apply(pos, &masks, m);
            let index = tables.ae_ix.index_fast(fast, &next, fast.a_rank(&next));
            let next_ae = decode(ae, tables.ae.get(index));
            let next_cb = self.cb(&next, cb);
            (next, next_ae, next_cb)
        }
        fn dfs(
            &self,
            pos: &Fc,
            ae: u8,
            cb: u8,
            path: &mut Vec<usize>,
            out: &mut Vec<Vec<usize>>,
        ) -> u64 {
            let left = self.depth - path.len();
            if left == 0 {
                if ae == 0 && cb == 0 && path.last().is_some_and(|&m| !is_a_move(m)) {
                    // the tables' targets are exactly H_A
                    assert!(self.tables.p2.index(&self.root.apply_seq(path)).is_some());
                    out.push(path.clone());
                }
                return 1;
            }
            let em = masks(&pos.e);
            let last = path.last().copied().unwrap_or(MAX_MOVES);
            let mut nodes = 1;
            for m in 0..N_MOVES {
                if !allowed_after(false, last, m) {
                    continue;
                }
                let (next, next_ae, next_cb) = self.child(pos, em, ae, cb, m);
                if usize::from(next_ae.max(next_cb)) < left {
                    path.push(m);
                    nodes += self.dfs(&next, next_ae, next_cb, path, out);
                    path.pop();
                }
            }
            nodes
        }
    }
    let walker = Walker {
        tables,
        fast: fast(),
        raw_cb,
        root: *root,
        depth,
    };
    let x0 = walker.fast.of(root);
    let (ae0, table_cb0) = tables.exact_fast(&x0);
    let cb0 = raw_cb.map_or(table_cb0, |_| walker.cb(&x0, 0));
    if depth < 2 {
        let mut out = vec![];
        let n = walker.dfs(&x0, ae0, cb0, &mut vec![], &mut out);
        return (out, n);
    }
    // split at the first two moves
    let mut prefixes = vec![];
    for m1 in 0..N_MOVES {
        for m2 in 0..N_MOVES {
            if allowed_after(false, m1, m2) {
                prefixes.push([m1, m2]);
            }
        }
    }
    let res: Vec<(Vec<Vec<usize>>, u64)> = prefixes
        .par_iter()
        .map(|&[m1, m2]| {
            let mut out = vec![];
            let (y1, a1, c1) = walker.child(&x0, masks(&x0.e), ae0, cb0, m1);
            if usize::from(a1.max(c1)) >= depth {
                return (out, 0);
            }
            let (y2, a2, c2) = walker.child(&y1, masks(&y1.e), a1, c1, m2);
            if usize::from(a2.max(c2)) > depth - 2 {
                return (out, 0);
            }
            let n = walker.dfs(&y2, a2, c2, &mut vec![m1, m2], &mut out);
            (out, n)
        })
        .collect();
    let nodes = res.iter().map(|r| r.1).sum();
    (res.into_iter().flat_map(|r| r.0).collect(), nodes)
}

/// The exact phase-1 distance (to `H_A`) of `y`, if at most `max`.
fn phase1_distance(t: &Tables<AeSym>, y: &Cube, max: usize) -> Option<usize> {
    if t.p2.index(y).is_some() {
        return Some(0);
    }
    (1..=max).find(|&d| !phase1_walks(t, y, d).0.is_empty())
}

/// A random word of `n` A-moves.
fn random_a_word(rng: &mut impl Rng, n: usize) -> Vec<usize> {
    let a: Vec<usize> = (0..N_MOVES).filter(|&m| is_a_move(m)).collect();
    (0..n).map(|_| a[rng.gen_range(0..a.len())]).collect()
}

/// A random position of `H'`: random corners and B-centres, the edges and A-centres of a
/// random element of `H_A`.
fn random_h_prime(rng: &mut impl Rng) -> Cube {
    let h = Cube::solved().apply_seq(&random_a_word(rng, 40));
    Cube {
        ep: h.ep,
        xa: h.xa,
        ..Cube::random(rng)
    }
}

/// All phase-1 walks of `y0` of lengths d1 and d1 + 1 (d1 its phase-1 distance), at most
/// `max` moves.
fn walks_to_next_depth(
    t: &Tables<AeSym>,
    y0: &Cube,
    max: usize,
) -> (Vec<Vec<usize>>, Option<usize>) {
    let mut ws = vec![];
    let mut d1 = None;
    for depth in 1..=max {
        let (w, _) = phase1_walks(t, y0, depth);
        if !w.is_empty() && d1.is_none() {
            d1 = Some(depth);
        }
        ws.extend(w);
        if d1.is_some_and(|d| depth > d) {
            break;
        }
    }
    (ws, d1)
}

/// The `H_A` pattern with the corners and edges of `c` (centres from the corners and edges).
fn h_a_pattern_with(c: &Cube) -> Cube {
    let xa: HashMap<[u8; 12], [u8; 12]> = h_a_configs().into_iter().collect();
    let attach = CbTable::attach_targets();
    Cube {
        xa: xa[&c.ep],
        xb: mset().unrank(attach[usize::from(corners().id(c))]),
        ..*c
    }
}

/// Real phase-1 walks of random positions `y0` of `H'` (all walks of lengths d1 and d1 + 1):
/// * on the planned block of `y0` (`z * s0^-1` for `z` in `H_A`, `s0` the first walk), most
///   other walks `s` take elements outside `H_A`;
/// * the phase-1 distance is not constant on that block;
/// * the element that walk `s` would mark as solved by `s` alone is not solved by it, and has
///   no solution of `|s|` moves at all (exhaustive walks of that element);
/// * every walk still works for the positions with the same corners and B-centres and any of
///   the 972 A-centre/edge configurations of `H'` (the walks depend only on those).
#[test]
#[ignore = "needs the medium tier tables in FTO_TABLES"]
fn real_walks_do_not_act_on_the_planned_blocks() {
    let t = medium();
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(7);
    let configs: Vec<_> = h_a_configs().into_iter().collect();
    let mut d1_varies = false;
    for k in 0..5 {
        let y0 = random_h_prime(&mut rng);
        let (ws, d1) = walks_to_next_depth(t, &y0, 13);
        if ws.len() < 2 {
            println!(
                "position {k}: d1 {d1:?}, {} walks of at most 13 moves: {y0:?}",
                ws.len()
            );
            continue;
        }
        let s0 = &ws[0];
        let (mut outside, mut total) = (0, 0);
        for s in &ws[1..] {
            for _ in 0..20 {
                let z = Cube::solved().apply_seq(&random_a_word(&mut rng, 40));
                let y = z.apply_seq(&inverse_word(s0));
                assert!(t.p2.index(&y.apply_seq(s0)).is_some());
                total += 1;
                if t.p2.index(&y.apply_seq(s)).is_none() {
                    outside += 1;
                }
            }
        }
        println!("position {k}: d1 {d1:?}, {} walks; planned block x other walk: {outside} of {total} outside H_A", ws.len());
        assert!(outside > 0);
        let ds: Vec<Option<usize>> = (0..4)
            .map(|_| {
                let z = Cube::solved().apply_seq(&random_a_word(&mut rng, 40));
                phase1_distance(t, &z.apply_seq(&inverse_word(s0)), s0.len())
            })
            .collect();
        println!("  d1 of four elements of the planned block: {ds:?}");
        d1_varies |= ds.iter().any(|&d| d != d1);
        // the element that the last walk would mark: corners and edges of solved * s^-1
        let s = &ws[ws.len() - 1];
        let r = Cube::solved().apply_seq(&inverse_word(s)).apply_seq(s0);
        let marked = h_a_pattern_with(&r).apply_seq(&inverse_word(s0));
        let end = marked.apply_seq(s);
        assert_eq!(
            (end.cp, end.co, end.ep),
            (Cube::solved().cp, Cube::solved().co, Cube::solved().ep)
        );
        if end != Cube::solved() {
            // no solution of |s| moves: every solution is a walk plus A-moves
            // (the empty walk when the element is in H_A itself)
            let best = (1..=s.len())
                .flat_map(|d| phase1_walks(t, &marked, d).0)
                .map(|w| w.len() + usize::from(t.p2.length(&marked.apply_seq(&w)).unwrap()))
                .chain(t.p2.length(&marked).map(usize::from))
                .min();
            println!(
                "  marked element: the method credits it with the {}-move walk alone; through \
                 walks of at most {} moves its shortest solution has {best:?} moves",
                s.len(),
                s.len()
            );
            assert!(best.is_none_or(|b| b > s.len()));
        }
        for (ep, xa) in configs.iter().step_by(37) {
            let y = Cube {
                ep: *ep,
                xa: *xa,
                ..y0
            };
            assert!(ws.iter().all(|s| t.p2.index(&y.apply_seq(s)).is_some()));
        }
    }
    assert!(d1_varies);
}

/// Measurement for the cost estimate in docs/BOUND.md (prints, asserts only consistency): for
/// 32 random blocks of the valid kind (one corner state and B-colouring, all 972 A-centre/edge
/// configurations of `H'`), the walks of lengths d1 and d1 + 1 (at most 14 moves, stopping
/// after a minute) and the resulting upper bound of every element (walk length + phase-2
/// length), with the nodes visited.
#[test]
#[ignore = "needs the medium tier tables in FTO_TABLES"]
fn block_bounds_sample() {
    let t = medium();
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(1234);
    let configs: Vec<_> = h_a_configs().into_iter().collect();
    let mut d1_hist = [0usize; 16];
    for k in 0..32 {
        let st = std::time::Instant::now();
        let y0 = random_h_prime(&mut rng);
        let mut best = vec![usize::MAX; configs.len()];
        let mut d1 = None;
        let mut nodes = 0;
        let mut per_depth = vec![];
        for depth in 1..=14usize {
            let (ws, n) = phase1_walks(t, &y0, depth);
            nodes += n;
            per_depth.push(ws.len());
            if !ws.is_empty() && d1.is_none() {
                d1 = Some(depth);
            }
            for s in &ws {
                for (b, (ep, xa)) in best.iter_mut().zip(&configs) {
                    let y = Cube {
                        ep: *ep,
                        xa: *xa,
                        ..y0
                    };
                    let l = t.p2.length(&y.apply_seq(s)).unwrap();
                    *b = (*b).min(depth + usize::from(l));
                }
            }
            if d1.is_some_and(|d| depth > d) || st.elapsed().as_secs() > 60 {
                break;
            }
        }
        let mut h: Vec<(usize, usize)> = vec![];
        for &b in &best {
            match h.iter_mut().find(|e| e.0 == b) {
                Some(e) => e.1 += 1,
                None => h.push((b, 1)),
            }
        }
        h.sort_unstable();
        if let Some(d) = d1 {
            d1_hist[d] += 1;
        }
        println!(
            "block {k}: d1 {d1:?}, walks per depth {per_depth:?}, {nodes} nodes, {:.1} s, \
             bounds over the 972 elements (bound, count) {h:?}",
            st.elapsed().as_secs_f64()
        );
    }
    println!("d1 histogram over the sampled blocks {d1_hist:?}");
}

/// The pruned walk enumeration finds exactly the walks of an unpruned search (all canonical
/// sequences), on positions whose phase-1 distance is small enough for the unpruned search.
#[test]
#[ignore = "needs the medium tier tables in FTO_TABLES"]
fn phase1_walks_match_unpruned_search() {
    fn unpruned(
        c: &Cube,
        last: usize,
        path: &mut Vec<usize>,
        depth: usize,
        p2: &Phase2,
        out: &mut Vec<Vec<usize>>,
    ) {
        if path.len() == depth {
            if path.last().is_some_and(|&m| !is_a_move(m)) && p2.index(c).is_some() {
                out.push(path.clone());
            }
            return;
        }
        for m in 0..N_MOVES {
            if allowed_after(false, last, m) {
                path.push(m);
                unpruned(&c.apply(m), m, path, depth, p2, out);
                path.pop();
            }
        }
    }
    let t = medium();
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(5);
    let mut found = 0;
    for _ in 0..24 {
        // an element of H_A followed by a few random moves
        let mut y = Cube::solved().apply_seq(&random_a_word(&mut rng, 40));
        for _ in 0..rng.gen_range(1..=5) {
            y = y.apply(rng.gen_range(0..N_MOVES));
        }
        for depth in 1..=6 {
            let (mut a, _) = phase1_walks(t, &y, depth);
            let mut b = vec![];
            unpruned(&y, MAX_MOVES, &mut vec![], depth, &t.p2, &mut b);
            a.sort_unstable();
            b.sort_unstable();
            assert_eq!(a, b);
            found += a.len();
        }
    }
    println!("{found} walks, all equal");
    assert!(found > 0);
}

/// Exact distances of every (corner state, B-colouring) pair to the 960 pairs of `H_A`, by a
/// breadth-first search of its own (no symmetry, one byte per pair; 355 MB), indexed
/// `corner id * N_MSET + B-colouring rank`.
fn raw_cb_distances() -> &'static [u8] {
    static D: OnceLock<Vec<u8>> = OnceLock::new();
    D.get_or_init(|| {
        use rayon::prelude::*;
        use std::sync::atomic::{AtomicU8, Ordering};
        let tables = move_tables();
        let size = N_CORNER * N_MSET;
        let dist: Vec<AtomicU8> = (0..size).map(|_| AtomicU8::new(u8::MAX)).collect();
        for (corner, &colouring) in CbTable::attach_targets().iter().enumerate() {
            dist[corner * N_MSET + colouring as usize].store(0, Ordering::Relaxed);
        }
        for level in 0..u8::MAX - 1 {
            let found: usize = (0..size)
                .into_par_iter()
                .with_min_len(1 << 16)
                .filter(|&i| dist[i].load(Ordering::Relaxed) == level)
                .map(|i| {
                    let (corner, colouring) = (i / N_MSET, i % N_MSET);
                    let mut new = 0;
                    for m in 0..N_MOVES {
                        let next = usize::from(tables.c[corner][m]) * N_MSET
                            + tables.b[colouring][m] as usize;
                        if dist[next]
                            .compare_exchange(
                                u8::MAX,
                                level + 1,
                                Ordering::Relaxed,
                                Ordering::Relaxed,
                            )
                            .is_ok()
                        {
                            new += 1;
                        }
                    }
                    new
                })
                .sum();
            if found == 0 {
                break;
            }
        }
        dist.into_iter().map(AtomicU8::into_inner).collect()
    })
}

/// The corner/B-centre table against [`raw_cb_distances`] on random pairs.
#[test]
#[ignore = "needs the medium tier tables in FTO_TABLES"]
fn cb_table_matches_raw_search() {
    let t = medium();
    let raw = raw_cb_distances();
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(3);
    let mut hist = [0usize; 16];
    for _ in 0..20_000 {
        let c = Cube::random(&mut rng);
        let x = Coords::of(&c);
        let want = raw[usize::from(x.c) * N_MSET + x.b as usize];
        assert_eq!(t.exact_fast(&fast().of(&c)).1, want);
        hist[usize::from(want)] += 1;
    }
    let mut all = [0u64; 16];
    for &d in raw {
        all[usize::from(d.min(15))] += 1;
    }
    println!(
        "raw CB distances over all {} pairs: {all:?}; sample {hist:?}",
        raw.len()
    );
}

/// A position of `H'` with no phase-1 walk of at most 13 moves, found by the real-walk test
/// (its fifth position), and its pattern (the `--pattern` format).
fn deep_position() -> Cube {
    Cube {
        cp: [0, 5, 3, 2, 4, 1],
        co: [0, 0, 0, 1, 1, 0],
        ep: [7, 8, 11, 4, 1, 10, 9, 0, 6, 3, 5, 2],
        xa: [3, 3, 3, 2, 2, 2, 1, 1, 1, 0, 0, 0],
        xb: [2, 0, 3, 0, 3, 3, 2, 1, 1, 1, 2, 0],
    }
}

const DEEP_PATTERN: &str =
    "OOOOGOOOW PRRYBBBBR AAAAAAAAG BYPYBRRPP / YPYPPYBRR WWWWWWWWA RBYRYBPPY OGGGOGGGG";

/// One of the 14-move walks of [`deep_position`].
const DEEP_WALK: &str = "B U' R L BR' F' BL' F BR BL L' R' U B'";

/// [`deep_position`] needs 14 moves to reach `H_A`: exhaustive walk searches with the AE table
/// and the corner/B-centre distances of [`raw_cb_distances`] find none of at most 13 moves and
/// 48 of 14. Its pattern and one walk are checked on their own (parsed pattern, physical moves).
/// Also prints its distance to `H_B` (the distance of the position conjugated by `SWAP_AB` to
/// `H_A`), if at most 14.
#[test]
#[ignore = "needs the medium tier tables in FTO_TABLES"]
fn deep_phase1_position() {
    use omelettehedron::pattern;
    let t = medium();
    let y0 = deep_position();
    assert_eq!(pattern::format(&pattern::of_cube(&y0)), DEEP_PATTERN);
    let parsed = pattern::position(&pattern::parse(DEEP_PATTERN).unwrap()).unwrap();
    assert_eq!(parsed, y0);
    let walk = parse_moves(DEEP_WALK).unwrap();
    let after = pattern::position(&pattern::after(&pattern::of_cube(&y0), &walk)).unwrap();
    assert!(t.p2.index(&after).is_some());
    assert_eq!(t.exact_fast(&fast().of(&y0)).0, 0);
    let swapped = y0.conjugate(&SWAP_AB).unwrap();
    let (ae_b, cb_b) = t.exact_fast(&fast().of(&swapped));
    let to_h_b =
        (usize::from(ae_b.max(cb_b))..=14).find(|&d| !phase1_walks(t, &swapped, d).0.is_empty());
    println!("distance to H_B: {to_h_b:?} (exact AE {ae_b}, CB {cb_b}; searched up to 14)");
    for depth in 1..=14 {
        let st = std::time::Instant::now();
        let (ws, nodes) = phase1_walks_with(t, &y0, depth, Some(raw_cb_distances()));
        println!(
            "depth {depth}: {} walks, {nodes} nodes, {:.1} s",
            ws.len(),
            st.elapsed().as_secs_f64()
        );
        assert_eq!(ws.is_empty(), depth < 14);
        if depth == 14 {
            assert_eq!(ws.len(), 48);
        }
    }
}

/// The ball of radius `radius` around `H'` in raw A-centre/edge configurations, with exact
/// distances, by breadth-first search from the 972 configurations of `H_A` (independent of
/// the AE table).
fn ae_ball(radius: u8) -> HashMap<AeKey, u8> {
    let mut ball: HashMap<AeKey, u8> = HashMap::new();
    let mut front = vec![];
    for (ep, xa) in h_a_configs() {
        let c = Cube {
            ep,
            xa,
            ..Cube::solved()
        };
        if ball.insert(ae_key(&c), 0).is_none() {
            front.push(c);
        }
    }
    for d in 1..=radius {
        let mut next = vec![];
        for c in &front {
            for m in 0..N_MOVES {
                let n = c.apply(m);
                if let std::collections::hash_map::Entry::Vacant(e) = ball.entry(ae_key(&n)) {
                    e.insert(d);
                    next.push(n);
                }
            }
        }
        front = next;
    }
    ball
}

/// [`deep_phase1_position`] without the AE table: every canonical sequence of at most 13
/// moves, pruned only by the corner/B-centre distances of [`raw_cb_distances`] and, with at
/// most 5 moves left, by the exact AE ball of radius 5. None ends in `H_A`.
#[test]
#[ignore = "slow: about 10^9 nodes (a few minutes at 48 threads)"]
fn deep_phase1_position_without_ae_table() {
    use rayon::prelude::*;
    struct Search<'a> {
        ball: &'a HashMap<AeKey, u8>,
        raw_cb: &'a [u8],
        p2: &'a Phase2,
        depth: usize,
    }
    impl Search<'_> {
        /// Walks of exactly `self.depth` moves below `pos` (`done` moves made, last `last`).
        fn count(&self, pos: &Cube, last: usize, done: usize) -> (u64, u64) {
            let left = self.depth - done;
            if left <= 5 {
                match self.ball.get(&ae_key(pos)) {
                    Some(&d) if usize::from(d) <= left => {}
                    _ => return (1, 0),
                }
            }
            // (corner/B-centre distances are at most 8)
            if left < 8 {
                let x = Coords::of(pos);
                if usize::from(self.raw_cb[usize::from(x.c) * N_MSET + x.b as usize]) > left {
                    return (1, 0);
                }
            }
            if left == 0 {
                let walk = !is_a_move(last) && self.p2.index(pos).is_some();
                return (1, u64::from(walk));
            }
            let mut total = (1, 0);
            for m in 0..N_MOVES {
                if allowed_after(false, last, m) {
                    let (n, w) = self.count(&pos.apply(m), m, done + 1);
                    total = (total.0 + n, total.1 + w);
                }
            }
            total
        }
    }
    let y0 = deep_position();
    let ball = ae_ball(5);
    let p2 = Phase2::new();
    for depth in 1..=13 {
        let search = Search {
            ball: &ball,
            raw_cb: raw_cb_distances(),
            p2: &p2,
            depth,
        };
        let mut prefixes = vec![];
        for m1 in 0..N_MOVES {
            for m2 in 0..N_MOVES {
                if allowed_after(false, m1, m2) {
                    prefixes.push([m1, m2]);
                }
            }
        }
        let st = std::time::Instant::now();
        let (nodes, walks) = if depth < 2 {
            search.count(&y0, MAX_MOVES, 0)
        } else {
            prefixes
                .par_iter()
                .map(|&[m1, m2]| search.count(&y0.apply(m1).apply(m2), m2, 2))
                .reduce(|| (0, 0), |a, b| (a.0 + b.0, a.1 + b.1))
        };
        println!(
            "depth {depth}: {walks} walks, {nodes} nodes, {:.1} s",
            st.elapsed().as_secs_f64()
        );
        assert_eq!(walks, 0);
    }
}
