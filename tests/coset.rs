//! Tests of the coset solver ([`omelettehedron::coset`]): the B-centres' labels, the moves on
//! the index, which side the sequences and the prepass act on, the elements behind the indices,
//! and exact distances against breadth-first search (small tier, built in memory).

use omelettehedron::coord::corners;
use omelettehedron::coset::enumerate::Phase1;
use omelettehedron::coset::space::{
    labelled_attachment, lift, unlabel, Marks, Space, N_A, N_COSET, UNSET,
};
use omelettehedron::coset::{self, invert, Options, Searched};
use omelettehedron::cube::*;
use omelettehedron::notation::{parse_moves, scramble_to_cube};
use omelettehedron::tables::*;
use rand::{Rng, SeedableRng};
use std::collections::HashMap;
use std::sync::OnceLock;

fn tables() -> &'static Tables<AeSym> {
    static T: OnceLock<Tables<AeSym>> = OnceLock::new();
    T.get_or_init(|| Tables::generate(AeSym::new(false), MoveSet::FACE, false))
}

fn space() -> &'static Space {
    static S: OnceLock<Space> = OnceLock::new();
    S.get_or_init(|| Space::new(&tables().p2))
}

/// The solved puzzle with its B-centres told apart (label i at B-centre position i).
fn labelled_solved() -> Cube {
    Cube {
        xb: std::array::from_fn(|i| i as u8),
        ..Cube::solved()
    }
}

fn scramble(s: &str) -> Cube {
    scramble_to_cube(&parse_moves(s).unwrap()).0
}

fn a_moves() -> Vec<usize> {
    (0..N_MOVES).filter(|&m| is_a_move(m)).collect()
}

/// In H_A the labelled B-centres follow the corners: random class-A words from the labelled
/// solved puzzle end with the labels the table gives for their corner state.
#[test]
fn labelled_b_centres_follow_corners() {
    let att = labelled_attachment().unwrap();
    assert_eq!(att.len(), 960);
    let a = a_moves();
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(1);
    for _ in 0..2000 {
        let mut c = labelled_solved();
        for _ in 0..rng.gen_range(0..30) {
            c = c.apply(a[rng.gen_range(0..a.len())]);
        }
        assert_eq!(att[usize::from(corners().id(&c))], c.xb);
    }
}

/// `lift` labels each colour's centres in position order, and `unlabel` undoes it.
#[test]
fn lift_round_trip() {
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(2);
    for _ in 0..100 {
        let c = Cube::random(&mut rng);
        let l = lift(&c);
        assert_eq!(unlabel(&l), c);
        let mut seen = [false; 12];
        for &x in &l.xb {
            assert!(!seen[usize::from(x)]);
            seen[usize::from(x)] = true;
        }
    }
}

/// The moves on the index are the class-A moves applied to the positions of H_A (every index,
/// every move), and `Phase2::cube` inverts `Phase2::index`.
#[test]
fn index_moves_are_the_right_action() {
    let p2 = &tables().p2;
    let sp = space();
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(3);
    for _ in 0..20_000 {
        let i = rng.gen_range(0..N_COSET);
        let c = p2.cube(i).unwrap();
        assert_eq!(p2.index(&c), Some(i));
        for k in 0..N_A {
            let j = p2.index(&c.apply(sp.moves[k])).unwrap();
            assert_eq!(sp.step(i as u32, k) as usize, j);
        }
    }
    assert!(p2.cube(N_COSET).is_none());
}

/// The trivial coset (the position is in H_A): the marks are the phase-2 distances.
#[test]
fn trivial_coset_is_phase_2() {
    let t = tables();
    let opt = Options {
        target: 12,
        max_search: 0,
        ..Options::default()
    };
    let c = coset::solve(t, &scramble("U F' BR BL'"), &opt, |_| {}).unwrap();
    assert_eq!(c.phase1_depth, Some(0));
    let h = coset::space::histogram(&c.bounds);
    assert_eq!(
        h,
        [1, 8, 48, 288, 1728, 9896, 51808, 220_111, 480_467, 166_276, 2457, 32]
    );
}

/// Breadth-first search over the labelled puzzle (B-centres told apart) from solved, to
/// `depth` moves: the distance of every position within it.
fn labelled_ball(depth: usize) -> HashMap<Cube, u8> {
    let mut dist = HashMap::new();
    dist.insert(labelled_solved(), 0u8);
    let mut fr = vec![labelled_solved()];
    for d in 1..=depth {
        let mut nf = vec![];
        for c in &fr {
            for m in 0..N_MOVES {
                let n = c.apply(m);
                if let std::collections::hash_map::Entry::Vacant(e) = dist.entry(n) {
                    e.insert(d as u8);
                    nf.push(n);
                }
            }
        }
        fr = nf;
    }
    dist
}

/// Every phase-1 length enumerated: the marks are exactly the distances of the labelled puzzle,
/// against breadth-first search, both ways (every element within the ball has its distance as
/// its mark, and every element marked within the ball's radius is in it). Also: the element of
/// each index is in the coset with that index, and the representative word makes the lift.
#[test]
fn marks_are_exact_labelled_distances() {
    const R: usize = 6;
    let t = tables();
    let ball = labelled_ball(R);
    for (k, s) in ["R", "R U'", "B D' L", "BR R' F L"].iter().enumerate() {
        let cube = scramble(s);
        let opt = Options {
            target: R,
            max_search: R,
            stop_when_covered: false,
            keep_paths: true,
            chunks: 1 + k,
            labelled: true,
            relabellings: false,
        };
        let c = coset::solve(t, &cube, &opt, |_| {}).unwrap();
        assert_eq!(c.exact_through, Some(R));
        let lb = Phase1::new(t, cube, None).unwrap().lower_bound();
        assert!(c
            .levels
            .iter()
            .all(|l| l.depth < lb || l.searched == Searched::Full));
        let rep = c.rep.clone().unwrap();
        let lab = labelled_solved().apply_seq(&rep);
        assert_eq!(unlabel(&lab), cube);
        let mut within = 0;
        for p in 0..N_COSET as u32 {
            let (word, z) = c.element(&t.p2, p).unwrap();
            let zl = labelled_solved().apply_seq(&word);
            assert_eq!(unlabel(&zl), z);
            let mark = c.marks.d[p as usize];
            match ball.get(&zl) {
                Some(&d) => {
                    within += 1;
                    assert_eq!(mark, d, "{s}: element {p}");
                }
                None => assert!(mark == UNSET || usize::from(mark) > R, "{s}: element {p}"),
            }
        }
        assert!(within > 0, "{s}: no element within {R}");
    }
}

/// Which side: an element's bound is the shortest solution through the sequences found, each
/// sequence applied to the element itself and completed in phase 2, and the phase-2 length
/// after sequence s equals the distance in H_A between the element's index and the state s
/// ends in (the right action). Every such solution solves the element.
#[test]
fn bounds_match_the_sequences_applied_to_the_elements() {
    let t = tables();
    let sp = space();
    let cube = scramble("R U L' B D BL'");
    let opt = Options {
        target: 13,
        max_search: 11,
        stop_when_covered: false,
        keep_paths: true,
        chunks: 3,
        labelled: true,
        relabellings: false,
    };
    let c = coset::solve(t, &cube, &opt, |_| {}).unwrap();
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(4);
    let seqs: Vec<(Vec<usize>, u32)> = c
        .sequences
        .iter()
        .flat_map(|e| {
            (0..e.index.len()).map(move |k| {
                let s = e.path(k).unwrap().iter().map(|&m| usize::from(m)).collect();
                (s, e.index[k])
            })
        })
        .collect();
    assert!(seqs.len() > 10);
    for _ in 0..40 {
        let p: u32 = rng.gen_range(0..N_COSET as u32);
        let (_, z) = c.element(&t.p2, p).unwrap();
        // distances in H_A from p (class-A moves are closed under inverses)
        let mut m = Marks::new();
        m.mark(p, 0);
        let from_p = m.project(sp);
        let mut best = UNSET as usize;
        for (s, e) in &seqs {
            let end = z.apply_seq(s);
            let d2 =
                t.p2.length(&end)
                    .expect("the sequence takes every element into H_A");
            assert_eq!(d2, from_p[*e as usize]);
            best = best.min(s.len() + usize::from(d2));
        }
        assert_eq!(usize::from(c.bounds[p as usize]), best);
        let sol = c.solution(&t.p2, &z).unwrap();
        assert_eq!(sol.len(), best);
        assert_eq!(z.apply_seq(&sol), Cube::solved());
    }
}

/// With the other relabellings' sequences: each improved bound equals the shortest solution
/// through every sequence kept (labelled or not), applied to the element itself wherever it
/// takes the element into H_A, and is never above the labelled bound.
#[test]
fn relabelled_sequences_improve_soundly() {
    let t = tables();
    let cube = scramble("R U L' B D BL' F R'");
    let opt = Options {
        target: 15,
        max_search: 12,
        stop_when_covered: false,
        keep_paths: true,
        chunks: 2,
        labelled: true,
        relabellings: true,
    };
    let c = coset::solve(t, &cube, &opt, |_| {}).unwrap();
    assert!(c.others.iter().map(|e| e.index.len()).sum::<usize>() > 0);
    let improved = c
        .bounds
        .iter()
        .zip(&c.labelled_bounds)
        .filter(|(a, b)| a < b)
        .count();
    assert!(improved > 0, "no bound improved");
    assert!(c.bounds.iter().zip(&c.labelled_bounds).all(|(a, b)| a <= b));
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(7);
    let mut seen_better = 0;
    for k in 0..60 {
        // half of them among the improved elements
        let p: u32 = if k % 2 == 0 {
            rng.gen_range(0..N_COSET as u32)
        } else {
            loop {
                let p = rng.gen_range(0..N_COSET as u32);
                if c.bounds[p as usize] < c.labelled_bounds[p as usize] {
                    break p;
                }
            }
        };
        let (_, z) = c.element(&t.p2, p).unwrap();
        let sol = c.solution_any(&t.p2, &z).unwrap();
        assert_eq!(sol.len(), usize::from(c.bounds[p as usize]), "element {p}");
        assert_eq!(z.apply_seq(&sol), Cube::solved());
        seen_better += usize::from(c.bounds[p as usize] < c.labelled_bounds[p as usize]);
    }
    assert!(seen_better >= 30);
}

/// Why the B-centres are told apart: marking every sequence that takes the position into H_A
/// as a pattern (`labelled: false`) gives marks that no sequence achieves on the element itself
/// (the optimal solver shows one of them, element 301457 of this coset with mark 10, needs 11
/// moves).
#[test]
fn pattern_marks_are_not_backed_by_the_sequences() {
    let t = tables();
    let opt = Options {
        target: 11,
        max_search: 11,
        stop_when_covered: false,
        keep_paths: true,
        chunks: 1,
        labelled: false,
        relabellings: false,
    };
    let c = coset::solve(t, &scramble("R U L' B"), &opt, |_| {}).unwrap();
    let unbacked = (0..N_COSET as u32)
        .step_by(97)
        .filter(|&p| c.marks.d[p as usize] != UNSET)
        .filter(|&p| {
            let (_, z) = c.element(&t.p2, p).unwrap();
            c.solution_any(&t.p2, &z)
                .is_none_or(|s| s.len() > usize::from(c.marks.d[p as usize]))
        })
        .count();
    assert!(unbacked > 0);
}

/// The inverse of an internal word undoes it.
#[test]
fn invert_undoes() {
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(5);
    let w: Vec<usize> = (0..25).map(|_| rng.gen_range(0..N_MOVES)).collect();
    let c = Cube::random(&mut rng);
    assert_eq!(c.apply_seq(&w).apply_seq(&invert(&w)), c);
}

/// The phase-1 depth of a position is the length of the shortest sequence the enumeration finds,
/// and the batch search agrees with the full enumeration.
#[test]
fn phase1_depth_matches_the_enumeration() {
    let t = tables();
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(6);
    for _ in 0..5 {
        let n = rng.gen_range(3..8);
        let mut cube = Cube::solved();
        for _ in 0..n {
            cube = cube.apply(rng.gen_range(0..N_MOVES));
        }
        let p1 = Phase1::new(t, cube, None).unwrap();
        let (d, _) = p1.depth(12);
        let (len, s) = d.unwrap();
        assert!(t.p2.index(&cube.apply_seq(&s)).is_some());
        assert!(len <= n);
        let stop = std::sync::atomic::AtomicBool::new(false);
        for l in p1.lower_bound()..=len {
            let mut found = 0;
            p1.level(l, false, 2, &stop, |e, _, _| {
                found += e.index.len();
                true
            });
            assert_eq!(found > 0, l == len, "length {l}");
        }
    }
}

/// Turns move pieces in 3-cycles, so the B-centres are always evenly permuted and a position
/// has 6^4 / 2 = 648 lifts (even relabellings of same-coloured B-centres). An element of H_A is
/// fixed only by the even relabellings that permute equal colours within each B-face: the mean
/// number over the 960 corner states is small, so an element can share only a few of the other
/// lifts' phase-1 sequences (see the module documentation of `coset`).
#[test]
fn few_relabellings_fix_an_element_of_h_a() {
    let att = labelled_attachment().unwrap();
    // every internal move permutes the B-centre positions evenly
    for m in 0..N_MOVES {
        let mv = defs().moves[m].xb;
        assert_eq!(perm_parity(&mv), 0, "move {m}");
    }
    let fact = [1u64, 1, 2, 6];
    let mut total = 0u64;
    for xb in att {
        let mut k = 1u64;
        let mut odd = false;
        for f in 0..4 {
            let mut n = [0usize; 4];
            for &l in &xb[3 * f..3 * f + 3] {
                n[usize::from(l / 3)] += 1;
            }
            k *= n.iter().map(|&c| fact[c]).product::<u64>();
            odd |= n.iter().any(|&c| c >= 2);
        }
        total += if odd { k / 2 } else { k };
    }
    let mean = total as f64 / att.len() as f64;
    eprintln!("mean even relabellings fixing an element of H_A: {mean:.3} of 648");
    assert!((1.0..20.0).contains(&mean));
}
