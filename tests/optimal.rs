//! End-to-end tests of the optimal solver (single-phase IDA* to the solved state; its small
//! tier is built in memory once). The two-phase solver has its own tests (tests/solver.rs).

use omelettehedron::cube::*;
use omelettehedron::notation::*;
use omelettehedron::optimal::coords::*;
use omelettehedron::optimal::search::{solve, Limits, Outcome, Stop};
use omelettehedron::optimal::tables::*;
use omelettehedron::prune::decode;
use rand::{Rng, SeedableRng};
use std::collections::{HashMap, HashSet};
use std::sync::OnceLock;
use std::time::Duration;

fn tables() -> &'static OptTables {
    static T: OnceLock<OptTables> = OnceLock::new();
    T.get_or_init(|| OptTables::generate(OptTier::Small, MoveSet::FACE, false).unwrap())
}

fn pool(threads: usize) -> &'static rayon::ThreadPool {
    static P1: OnceLock<rayon::ThreadPool> = OnceLock::new();
    static P8: OnceLock<rayon::ThreadPool> = OnceLock::new();
    let p = if threads == 1 { &P1 } else { &P8 };
    p.get_or_init(|| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
    })
}

/// Solve on 1 or 8 threads.
fn solve_on(threads: usize, c: Cube, fr: Frame, limits: &Limits) -> Outcome {
    let b = tables().small().unwrap();
    pool(threads).install(|| solve(&b, c, fr, limits, false))
}

fn random_scramble(rng: &mut impl Rng, n: usize) -> Vec<PhysMove> {
    (0..n)
        .map(|_| PhysMove::face(rng.gen_range(0..8), rng.gen()))
        .collect()
}

/// The solution solves the scramble on the sticker model (up to a whole-puzzle rotation).
fn check_solution(scr: &[PhysMove], sol: &[PhysMove]) {
    let st = apply_physical(&apply_physical(&solved_stickers(defs()), scr), sol);
    assert!(
        is_solved_physical(&st),
        "solution does not solve: {} / {}",
        format_moves(scr),
        format_moves(sol)
    );
}

/// Exact distances up to 2 R by meet in the middle (fixed frame = solved up to rotation),
/// independent of the solver's tables and search (as in tests/solver.rs).
struct Mitm {
    ball: HashMap<Cube, u8>,
    moves: MoveSet,
}
const R: usize = 5;
impl Mitm {
    fn new() -> Mitm {
        Self::with(MoveSet::FACE)
    }
    /// The brute force in the metric of `moves`.
    fn with(moves: MoveSet) -> Mitm {
        let mut ball = HashMap::new();
        ball.insert(Cube::solved(), 0u8);
        let mut fr = vec![Cube::solved()];
        for d in 1..=R as u8 {
            let mut nf = vec![];
            for c in &fr {
                for m in 0..moves.n_moves() {
                    let n = c.apply(m);
                    ball.entry(n).or_insert_with(|| {
                        nf.push(n);
                        d
                    });
                }
            }
            fr = nf;
        }
        Mitm { ball, moves }
    }
    fn distance(&self, c: &Cube) -> Option<usize> {
        let mut best = None::<usize>;
        for k in 0..=R {
            let mut stack = vec![(*c, 0usize, usize::MAX)];
            while let Some((x, g, last)) = stack.pop() {
                if g == k {
                    if let Some(&b) = self.ball.get(&x) {
                        let t = k + b as usize;
                        if best.is_none_or(|bb| t < bb) {
                            best = Some(t);
                        }
                    }
                    continue;
                }
                for m in 0..self.moves.n_moves() {
                    if self.moves.allowed_after(last, m) {
                        stack.push((x.apply(m), g + 1, m));
                    }
                }
            }
            if let Some(b) = best {
                if b <= k + R {
                    return Some(b);
                }
            }
        }
        best
    }
}

fn mitm() -> &'static Mitm {
    static M: OnceLock<Mitm> = OnceLock::new();
    M.get_or_init(Mitm::new)
}

/// Optimal lengths equal the brute force on short scrambles (up to 10 moves), on 1 and 8
/// threads, and every solution solves on the sticker model.
#[test]
fn optimal_matches_brute_force() {
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(101);
    let mut checked = 0;
    for len in [0usize, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] {
        for _ in 0..5 {
            let scr = random_scramble(&mut rng, len);
            let (c, fr) = scramble_to_cube(&scr);
            let Some(truth) = mitm().distance(&c) else {
                continue;
            };
            for threads in [1, 8] {
                let o = solve_on(threads, c, fr, &Limits::default());
                assert_eq!(o.stop, Stop::Solved);
                assert_eq!(o.length, Some(truth), "scramble {}", format_moves(&scr));
                assert_eq!(o.lower_bound, truth);
                assert_eq!(o.solutions.len(), 1);
                assert_eq!(o.solutions[0].len(), truth);
                check_solution(&scr, &o.solutions[0]);
            }
            checked += 1;
        }
    }
    assert!(checked >= 55, "{checked}");
}

/// Many more short scrambles (slow; `cargo test --release -- --ignored`).
#[test]
#[ignore = "slow: 300 positions"]
fn optimal_matches_brute_force_many() {
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(102);
    for i in 0..300 {
        let scr = random_scramble(&mut rng, 6 + i % 7);
        let (c, fr) = scramble_to_cube(&scr);
        let Some(truth) = mitm().distance(&c) else {
            continue;
        };
        let o = solve_on(8, c, fr, &Limits::default());
        assert_eq!(o.length, Some(truth), "scramble {}", format_moves(&scr));
        check_solution(&scr, &o.solutions[0]);
    }
}

/// The two-phase solver's optimal mode is exact too: the lengths must agree on medium
/// scrambles (beyond the brute force's reach).
#[test]
fn optimal_matches_two_phase_optimal_mode() {
    use omelettehedron::search as tp;
    use omelettehedron::tables::{AeSym, Tables};
    let two_phase = Tables::generate(AeSym::new(false), MoveSet::FACE, false);
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(103);
    for len in [11usize, 12, 13, 14, 14, 15] {
        let scr = random_scramble(&mut rng, len);
        let (c, fr) = scramble_to_cube(&scr);
        let want = pool(8)
            .install(|| tp::solve(&two_phase, c, fr, &tp::Limits::default(), tp::Report::Quiet));
        assert!(want.optimal);
        let o = solve_on(8, c, fr, &Limits::default());
        assert_eq!(o.length, want.length, "scramble {}", format_moves(&scr));
        check_solution(&scr, &o.solutions[0]);
    }
}

/// All canonical move sequences of length `len` that solve `c` (brute force).
fn all_solutions(c: &Cube, len: usize) -> HashSet<Vec<usize>> {
    fn rec(x: &Cube, len: usize, path: &mut Vec<usize>, out: &mut HashSet<Vec<usize>>) {
        if path.len() == len {
            if *x == Cube::solved() {
                out.insert(path.clone());
            }
            return;
        }
        let last = path.last().copied().unwrap_or(usize::MAX);
        for m in 0..N_MOVES {
            if allowed_after(false, last, m) {
                path.push(m);
                rec(&x.apply(m), len, path, out);
                path.pop();
            }
        }
    }
    let mut out = HashSet::new();
    rec(c, len, &mut vec![], &mut out);
    out
}

#[test]
fn stopping_rules() {
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(104);
    // -n: up to N optimal solutions, all of them when there are fewer, the same set on 1 and
    // on 8 threads (checked against all canonical solutions of that length). Short scrambles
    // almost always have a single optimal canonical solution; the position solved by both
    // sides of the relation `BLw' B' BLw' BRw' BLw' = BRw' BLw' BRw' B' BRw'` (internal moves,
    // docs/PHASES.md section 8) has at least two.
    let rel_a = [13usize, 1, 13, 11, 13];
    let rel_b = [11usize, 13, 11, 1, 11];
    let inv_a: Vec<usize> = rel_a.iter().rev().map(|&m| inverse_move(m)).collect();
    let mut scrambles = vec![solution_to_physical(Frame::identity(), &inv_a, false)];
    for len in [3, 4, 5, 6] {
        scrambles.push(random_scramble(&mut rng, len));
    }
    let mut multi = 0;
    for scr in &scrambles {
        let (c, fr) = scramble_to_cube(scr);
        let truth = mitm().distance(&c).unwrap();
        let every = all_solutions(&c, truth);
        let lim = |n| Limits {
            max_solutions: Some(n),
            ..Limits::default()
        };
        for threads in [1, 8] {
            let o = solve_on(threads, c, fr, &lim(1_000_000));
            assert_eq!(o.length, Some(truth));
            let got: HashSet<Vec<usize>> = o.internal.iter().cloned().collect();
            assert_eq!(got.len(), o.internal.len(), "duplicate solutions");
            assert_eq!(got, every, "scramble {}", format_moves(scr));
            for s in &o.solutions {
                check_solution(scr, s);
            }
            let o = solve_on(threads, c, fr, &lim(1));
            assert_eq!(o.solutions.len(), 1);
            let o = solve_on(threads, c, fr, &lim(2));
            assert_eq!(o.solutions.len(), every.len().min(2));
            for s in &o.internal {
                assert!(every.contains(s));
            }
        }
        multi += usize::from(every.len() > 1);
    }
    // (the relation's position: both sides are optimal solutions)
    let (c, _) = scramble_to_cube(&scrambles[0]);
    assert_eq!(c.apply_seq(&rel_a), Cube::solved());
    assert_eq!(c.apply_seq(&rel_b), Cube::solved());
    assert!(multi > 0, "no position with several optimal solutions");
    // -l: give up when the optimum is longer
    let scr = random_scramble(&mut rng, 9);
    let (c, fr) = scramble_to_cube(&scr);
    let truth = mitm().distance(&c).unwrap();
    let lim = |l| Limits {
        max_length: Some(l),
        ..Limits::default()
    };
    let o = solve_on(8, c, fr, &lim(truth - 1));
    assert_eq!(o.stop, Stop::MaxLength);
    assert!(o.solutions.is_empty() && o.length.is_none());
    assert_eq!(o.lower_bound, truth);
    for l in [truth, truth + 3] {
        let o = solve_on(8, c, fr, &lim(l));
        assert_eq!(o.length, Some(truth));
        check_solution(&scr, &o.solutions[0]);
    }
    // -t: a random position cannot be solved in 0.2 s with the small tier; the lower bound is
    // at least the root bound, and the search stops promptly
    let c = Cube::random(&mut rng);
    let h0 = exact(&tables().small().unwrap(), &Pos::of(&c)).unwrap();
    for threads in [1, 8] {
        let o = solve_on(
            threads,
            c,
            Frame::identity(),
            &Limits {
                time: Some(Duration::from_millis(200)),
                ..Limits::default()
            },
        );
        assert_eq!(o.stop, Stop::Time);
        assert!(o.solutions.is_empty());
        assert!(o.lower_bound >= usize::from(*h0.iter().max().unwrap()));
        assert!(o.elapsed < 2.0, "{}", o.elapsed);
    }
}

/// One thread searches the same tree every time: identical node counts; and the per-iteration
/// counts add up.
#[test]
fn node_counts_are_deterministic_on_one_thread() {
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(105);
    for _ in 0..3 {
        let scr = random_scramble(&mut rng, 12);
        let (c, fr) = scramble_to_cube(&scr);
        let a = solve_on(1, c, fr, &Limits::default());
        let b = solve_on(1, c, fr, &Limits::default());
        assert_eq!(a.nodes, b.nodes);
        assert_eq!(a.internal, b.internal);
        assert_eq!(
            a.iterations.iter().map(|i| i.nodes).sum::<u64>(),
            a.nodes,
            "iterations"
        );
        let m = solve_on(8, c, fr, &Limits::default());
        assert_eq!(m.length, a.length);
    }
}

/// Exact table values: a child's decoded value equals its own walked distance, values change
/// by at most one per move, the solved state is 0 everywhere and every other position has a
/// positive maximum, and the bound never exceeds the true distance.
#[test]
fn tables_are_consistent() {
    let b = tables().small().unwrap();
    let mt = opt_moves();
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(106);
    let solved = Pos::solved();
    assert_eq!(exact(&b, &solved), Some([0; MAX_K]));
    for i in 0..200 {
        let c = if i % 2 == 0 {
            Cube::random(&mut rng)
        } else {
            Cube::solved().apply_seq(
                &(0..rng.gen_range(1..8))
                    .map(|_| rng.gen_range(0..N_MOVES))
                    .collect::<Vec<_>>(),
            )
        };
        let p = Pos::of(&c);
        let h = exact(&b, &p).unwrap();
        if p != solved {
            assert!(h.iter().any(|&x| x > 0));
        }
        if i % 2 == 1 {
            if let Some(d) = mitm().distance(&c) {
                assert!(h.iter().all(|&x| usize::from(x) <= d));
            }
        }
        for m in 0..N_MOVES {
            let q = mt.apply(&p, m);
            let hq = exact(&b, &q).unwrap();
            for k in 0..SmallBounds::K {
                assert!(h[k].abs_diff(hq[k]) <= 1, "lookup {k}");
                assert_eq!(decode(h[k], b.table(k).get(b.index(k, &q))), hq[k]);
            }
        }
    }
}

/// The edge table's first levels equal an independent breadth-first search over edge
/// permutations (the cubie model's moves), and its histogram covers every even permutation.
#[test]
fn edge_table_matches_independent_bfs() {
    let edges = tables().edges.as_ref().unwrap();
    assert_eq!(edges.hist.iter().sum::<u64>(), N_EPERM as u64);
    let mut seen: HashSet<[u8; 12]> = HashSet::from([Cube::solved().ep]);
    let mut front = vec![Cube::solved()];
    for d in 1..=5 {
        let mut next = vec![];
        for c in &front {
            for m in 0..N_MOVES {
                let y = c.apply(m);
                if seen.insert(y.ep) {
                    next.push(y);
                }
            }
        }
        assert_eq!(next.len() as u64, edges.hist[d], "depth {d}");
        front = next;
    }
}

/// The medium tier (corners x edges, 2.4 GB: generated or loaded from `$FTO_TABLES`): the same
/// optimal lengths as the brute force, consistent tables, and never a weaker root bound than the
/// small tier's edge table.
#[test]
#[ignore = "needs the medium tier's tables (2.4 GB; `omelettehedron optimal gen --tier medium`)"]
fn medium_tier_matches_brute_force() {
    let dir = omelettehedron::tables::table_dir(None);
    let t = OptTables::load(&dir, OptTier::Medium, MoveSet::FACE, false).unwrap();
    let b = t.medium().unwrap();
    let small = tables().small().unwrap();
    let mt = opt_moves();
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(107);
    for len in [0usize, 2, 4, 6, 8, 9, 10, 11, 12] {
        for _ in 0..4 {
            let scr = random_scramble(&mut rng, len);
            let (c, fr) = scramble_to_cube(&scr);
            let Some(truth) = mitm().distance(&c) else {
                continue;
            };
            let o = pool(8).install(|| solve(&b, c, fr, &Limits::default(), false));
            assert_eq!(o.length, Some(truth), "scramble {}", format_moves(&scr));
            check_solution(&scr, &o.solutions[0]);
        }
    }
    for _ in 0..100 {
        let p = Pos::of(&Cube::random(&mut rng));
        let h = exact(&b, &p).unwrap();
        // corners x edges >= edges
        assert!(h[0] >= exact(&small, &p).unwrap()[0]);
        for m in 0..N_MOVES {
            let q = mt.apply(&p, m);
            let hq = exact(&b, &q).unwrap();
            for k in 0..MediumBounds::K {
                assert!(h[k].abs_diff(hq[k]) <= 1);
                assert_eq!(decode(h[k], b.table(k).get(b.index(k, &q))), hq[k]);
            }
        }
    }
}

/// The large tier exists in the face metric only: with slices it is refused before any table is
/// read or generated.
#[test]
fn large_tier_refuses_slices() {
    let slice = MoveSet {
        slice: true,
        ..MoveSet::FACE
    };
    let dir = std::env::temp_dir().join("omelettehedron-no-tables");
    let e = OptTables::load(&dir, OptTier::Large, slice, false)
        .err()
        .unwrap();
    assert!(e.to_string().contains("slice"), "{e}");
    assert!(!dir.join("opt-cbe-sym6-slice.prun").exists());
}

/// The large tier (the medium tier and the 21 GB corner x B-centre x edge-triple table: loaded
/// from `$FTO_TABLES`, generated if missing): the same optimal lengths as the brute force,
/// bounds never above the brute-force distance (admissible), consistent tables, the medium
/// tier's corner x edge and A-centre/edge values, and triple lookups never below the corner x
/// centre values they replace.
#[test]
#[ignore = "needs the large tier's tables (24.4 GB; `omelettehedron optimal gen --tier large`)"]
fn large_tier_matches_brute_force() {
    let dir = omelettehedron::tables::table_dir(None);
    let t = OptTables::load(&dir, OptTier::Large, MoveSet::FACE, false).unwrap();
    let b = t.large().unwrap();
    let medium = t.medium().unwrap();
    let mt = opt_moves();
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(108);
    assert_eq!(exact(&b, &Pos::solved()), Some([0; MAX_K]));
    for len in [0usize, 2, 4, 6, 8, 9, 10, 11, 12] {
        for _ in 0..4 {
            let scr = random_scramble(&mut rng, len);
            let (c, fr) = scramble_to_cube(&scr);
            let Some(truth) = mitm().distance(&c) else {
                continue;
            };
            let h = exact(&b, &Pos::of(&c)).unwrap();
            assert!(h.iter().all(|&x| usize::from(x) <= truth));
            let o = pool(8).install(|| solve(&b, c, fr, &Limits::default(), false));
            assert_eq!(o.length, Some(truth), "scramble {}", format_moves(&scr));
            check_solution(&scr, &o.solutions[0]);
        }
    }
    // admissible: never above the brute-force distance (positions up to 10 moves out)
    let mut tight = 0;
    for _ in 0..300 {
        let n = rng.gen_range(1..=10);
        let c = Cube::solved().apply_seq(
            &(0..n)
                .map(|_| rng.gen_range(0..N_MOVES))
                .collect::<Vec<_>>(),
        );
        let Some(d) = mitm().distance(&c) else {
            continue;
        };
        let h = exact(&b, &Pos::of(&c)).unwrap();
        assert!(
            h.iter().all(|&x| usize::from(x) <= d),
            "bound above distance {d}"
        );
        tight += usize::from(h[3..LargeBounds::K].iter().any(|&x| usize::from(x) == d));
    }
    assert!(tight > 0);
    for _ in 0..100 {
        let p = Pos::of(&Cube::random(&mut rng));
        let h = exact(&b, &p).unwrap();
        let hm = exact(&medium, &p).unwrap();
        // corners x edges and the A-centre/edge lookups are the medium tier's; each A-triple
        // lookup is at least corners x B-centres, each B-triple lookup corners x A-centres
        assert_eq!([h[0], h[1], h[2]], [hm[0], hm[3], hm[4]]);
        assert!(h[3..7].iter().all(|&x| x >= hm[1]));
        assert!(h[7..11].iter().all(|&x| x >= hm[2]));
        for m in 0..N_MOVES {
            let q = mt.apply(&p, m);
            let hq = exact(&b, &q).unwrap();
            for k in 0..LargeBounds::K {
                assert!(h[k].abs_diff(hq[k]) <= 1);
                assert_eq!(decode(h[k], b.table(k).get(b.index(k, &q))), hq[k]);
            }
        }
    }
}

/// Symmetric start: on symmetric positions the filtered search (default) finds the same
/// optimal lengths as the brute force and as the unfiltered search (`-n` keeps every task). The
/// positions: every position within 5 moves (breadth-first, in a fixed order) that some of the
/// 47 non-trivial symmetries fix, at distance 3 or more.
#[test]
fn symmetric_start_keeps_optimal_lengths() {
    use omelettehedron::optimal::sym::{sym48, N_SYM};
    let s48 = sym48();
    let mut seen: HashSet<Cube> = HashSet::from([Cube::solved()]);
    let mut front = vec![(Cube::solved(), vec![])];
    let mut symmetric: Vec<(Cube, Vec<usize>)> = vec![];
    for d in 1..=5 {
        let mut next = vec![];
        for (c, path) in &front {
            for m in 0..N_MOVES {
                let y = c.apply(m);
                if seen.insert(y) {
                    let mut p: Vec<usize> = path.clone();
                    p.push(m);
                    if d >= 3 && (1..N_SYM).any(|k| s48.conj_by_maps(&y, k) == y) {
                        symmetric.push((y, p.clone()));
                    }
                    next.push((y, p));
                }
            }
        }
        front = next;
    }
    assert!(symmetric.len() >= 10, "{}", symmetric.len());
    let (mut filtered_nodes, mut full_nodes) = (0u64, 0u64);
    for (c, path) in symmetric.iter().take(60) {
        let scr = solution_to_physical(Frame::identity(), path, false);
        // (the physical scramble ends in the frame its wide internal moves turned to)
        let (c2, fr) = scramble_to_cube(&scr);
        assert_eq!(c2, *c);
        let o = solve_on(1, *c, fr, &Limits::default());
        assert!(o.symmetries > 1);
        let truth = mitm().distance(c).unwrap();
        assert_eq!(o.length, Some(truth), "{}", format_moves(&scr));
        check_solution(&scr, &o.solutions[0]);
        let full = solve_on(
            1,
            *c,
            fr,
            &Limits {
                max_solutions: Some(1_000_000),
                ..Limits::default()
            },
        );
        assert_eq!(full.length, Some(truth));
        let m = solve_on(8, *c, fr, &Limits::default());
        assert_eq!(m.length, Some(truth));
        filtered_nodes += o.nodes;
        full_nodes += full.nodes;
    }
    eprintln!(
        "{} symmetric positions: {filtered_nodes} nodes, unfiltered (-n) {full_nodes}",
        symmetric.len()
    );
    assert!(filtered_nodes < full_nodes);
}

const SLICE: MoveSet = MoveSet {
    wide: false,
    slice: true,
};

/// The small tier in the metric with slice turns (built in memory once).
fn slice_tables() -> &'static OptTables {
    static T: OnceLock<OptTables> = OnceLock::new();
    T.get_or_init(|| OptTables::generate(OptTier::Small, SLICE, false).unwrap())
}

fn slice_mitm() -> &'static Mitm {
    static M: OnceLock<Mitm> = OnceLock::new();
    M.get_or_init(|| Mitm::with(SLICE))
}

/// Random moves of every layer kind (face, wide and slice turns).
fn random_layer_scramble(rng: &mut impl Rng, n: usize) -> Vec<PhysMove> {
    (0..n)
        .map(|_| PhysMove {
            face: rng.gen_range(0..8),
            layer: [Layer::Face, Layer::Wide, Layer::Slice][rng.gen_range(0..3)],
            cw: rng.gen(),
        })
        .collect()
}

/// With slice turns (`--moves slice`): the brute force's lengths in the slice metric, on 1 and 8
/// threads, for scrambles of face, wide and slice turns; the solutions use face and slice turns
/// only, solve on the sticker model, use slices where they help, and are never longer than in
/// the face metric.
#[test]
fn optimal_matches_brute_force_with_slices() {
    let b = slice_tables().small().unwrap();
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(111);
    let (mut slices_used, mut shorter, mut checked) = (0, 0, 0);
    for len in [0usize, 1, 2, 3, 5, 7, 8, 9, 10] {
        for k in 0..5 {
            let scr = if k % 2 == 0 {
                random_scramble(&mut rng, len)
            } else {
                random_layer_scramble(&mut rng, len)
            };
            let (c, fr) = scramble_to_cube(&scr);
            let Some(truth) = slice_mitm().distance(&c) else {
                continue;
            };
            for threads in [1, 8] {
                let o = pool(threads).install(|| solve(&b, c, fr, &Limits::default(), false));
                assert_eq!(o.length, Some(truth), "scramble {}", format_moves(&scr));
                let sol = &o.solutions[0];
                check_solution(&scr, sol);
                assert!(sol.iter().all(|m| m.layer != Layer::Wide));
                slices_used += sol.iter().filter(|m| m.layer == Layer::Slice).count();
            }
            if k % 2 == 0 {
                let face = solve_on(8, c, fr, &Limits::default());
                assert!(truth <= face.length.unwrap());
                shorter += usize::from(truth < face.length.unwrap());
            }
            checked += 1;
        }
    }
    assert!(checked >= 40 && slices_used > 0 && shorter > 0);
}

/// The slice-metric tables are consistent: decoded child values (by every move, slices too)
/// equal walked exact values, and values change by at most one per move.
#[test]
fn slice_tables_are_consistent() {
    let b = slice_tables().small().unwrap();
    let mt = opt_moves();
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(112);
    for i in 0..60 {
        let c = if i % 2 == 0 {
            Cube::random(&mut rng)
        } else {
            Cube::solved().apply_seq(
                &(0..rng.gen_range(1..8))
                    .map(|_| rng.gen_range(0..MAX_MOVES))
                    .collect::<Vec<_>>(),
            )
        };
        let p = Pos::of(&c);
        let h = exact(&b, &p).unwrap();
        if i % 2 == 1 {
            if let Some(d) = slice_mitm().distance(&c) {
                assert!(h.iter().all(|&x| usize::from(x) <= d));
            }
        }
        for m in 0..MAX_MOVES {
            let q = mt.apply(&p, m);
            let hq = exact(&b, &q).unwrap();
            for k in 0..SmallBounds::K {
                assert!(h[k].abs_diff(hq[k]) <= 1, "lookup {k}, move {m}");
                assert_eq!(decode(h[k], b.table(k).get(b.index(k, &q))), hq[k]);
            }
        }
    }
}

/// Symmetric start with slices: symmetric positions (fixed by some symmetry, within 4 moves of
/// the slice metric) get the brute force's lengths, filtered and unfiltered (`-n`).
#[test]
fn symmetric_start_with_slices() {
    use omelettehedron::optimal::sym::{sym48, N_SYM};
    let s48 = sym48();
    let b = slice_tables().small().unwrap();
    let mut seen: HashSet<Cube> = HashSet::from([Cube::solved()]);
    let mut front = vec![(Cube::solved(), vec![])];
    let mut symmetric: Vec<(Cube, Vec<usize>)> = vec![];
    for d in 1..=4 {
        let mut next = vec![];
        for (c, path) in &front {
            for m in 0..MAX_MOVES {
                let y = c.apply(m);
                if seen.insert(y) {
                    let mut p: Vec<usize> = path.clone();
                    p.push(m);
                    if d >= 3 && (1..N_SYM).any(|k| s48.conj_by_maps(&y, k) == y) {
                        symmetric.push((y, p.clone()));
                    }
                    next.push((y, p));
                }
            }
        }
        front = next;
    }
    assert!(symmetric.len() >= 10, "{}", symmetric.len());
    for (c, path) in symmetric.iter().take(40) {
        let scr = solution_to_physical(Frame::identity(), path, false);
        let (c2, fr) = scramble_to_cube(&scr);
        assert_eq!(c2, *c);
        let truth = slice_mitm().distance(c).unwrap();
        let o = pool(1).install(|| solve(&b, *c, fr, &Limits::default(), false));
        assert!(o.symmetries > 1);
        assert_eq!(o.length, Some(truth), "{}", format_moves(&scr));
        check_solution(&scr, &o.solutions[0]);
        let full = pool(1).install(|| {
            solve(
                &b,
                *c,
                fr,
                &Limits {
                    max_solutions: Some(1_000_000),
                    ..Limits::default()
                },
                false,
            )
        });
        assert_eq!(full.length, Some(truth));
    }
}

/// With slices, the two-phase solver's optimal mode (exact in the slice metric too) gives the
/// same lengths on medium scrambles of face, wide and slice turns.
#[test]
fn slice_optimal_matches_two_phase_optimal_mode() {
    use omelettehedron::search as tp;
    use omelettehedron::tables::{AeSym, Tables};
    let two_phase = Tables::generate(AeSym::new(false), SLICE, false);
    let b = slice_tables().small().unwrap();
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(113);
    for len in [10usize, 11, 12, 12] {
        let scr = random_layer_scramble(&mut rng, len);
        let (c, fr) = scramble_to_cube(&scr);
        let want = pool(8)
            .install(|| tp::solve(&two_phase, c, fr, &tp::Limits::default(), tp::Report::Quiet));
        assert!(want.optimal);
        let o = pool(8).install(|| solve(&b, c, fr, &Limits::default(), false));
        assert_eq!(o.length, want.length, "scramble {}", format_moves(&scr));
        check_solution(&scr, &o.solutions[0]);
    }
}
