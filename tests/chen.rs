//! End-to-end tests of the block-split solver (`omelettehedron::chen`), with the small tier's
//! phase-1 table and the phase-2 tables built in memory once.

use omelettehedron::chen::relaxed::{relaxed_targets, target_check, RelaxedAe};
use omelettehedron::chen::search::solve;
use omelettehedron::chen::space::space;
use omelettehedron::chen::tables::{ChenTables, TableSet, P2};
use omelettehedron::cube::*;
use omelettehedron::notation::*;
use omelettehedron::search::*;
use omelettehedron::tables::*;
use rand::{Rng, SeedableRng};
use std::sync::OnceLock;
use std::time::Duration;

fn tables() -> &'static (Tables<AeSym>, ChenTables) {
    static T: OnceLock<(Tables<AeSym>, ChenTables)> = OnceLock::new();
    T.get_or_init(|| {
        (
            Tables::generate(AeSym::new(false), MoveSet::FACE, false),
            ChenTables::generate(TableSet { sb: true, px: true }, false),
        )
    })
}

fn check_solution(scr: &[PhysMove], sol: &[PhysMove]) {
    let st = apply_physical(&apply_physical(&solved_stickers(defs()), scr), sol);
    assert!(
        is_solved_physical(&st),
        "solution does not solve: {} / {}",
        format_moves(scr),
        format_moves(sol)
    );
}

fn one_thread<R: Send>(f: impl FnOnce() -> R + Send) -> R {
    rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap()
        .install(f)
}

/// Random-move scrambles of every length up to 12 and random positions: the first solution
/// (with `-n 1`) solves the position, and so does the best after a short time.
#[test]
fn solutions_solve() {
    let (t, ct) = tables();
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(11);
    for len in 0..=12 {
        let scr: Vec<PhysMove> = (0..len)
            .map(|_| PhysMove::face(rng.gen_range(0..8), rng.gen()))
            .collect();
        let (cube, frame) = scramble_to_cube(&scr);
        let limits = Limits {
            max_solutions: Some(1),
            ..Limits::default()
        };
        let o = one_thread(|| {
            solve(
                t,
                ct,
                None,
                cube,
                frame,
                &limits,
                Report::Quiet,
                &Default::default(),
            )
        });
        let sol = o.solution.expect("a solution");
        check_solution(&scr, &sol);
        assert!(!o.optimal);
        if len == 0 {
            assert!(sol.is_empty());
        }
    }
    for i in 0..3 {
        let cube = Cube::random(&mut rng);
        // (two solutions: the second improves on the first)
        let limits = Limits {
            max_solutions: Some(2),
            time: Some(Duration::from_secs(600)),
            ..Limits::default()
        };
        let o = solve(
            t,
            ct,
            None,
            cube,
            Frame::identity(),
            &limits,
            Report::Quiet,
            &Default::default(),
        );
        let sol = o.solution.expect("a solution");
        // (a random position in the fixed frame: the physical scramble is the inverse)
        let st = apply_physical(&cube.to_stickers(), &sol);
        assert!(is_solved_physical(&st), "random position {i}");
        // the events are improvements, ending in the best
        assert!(o.events.windows(2).all(|w| w[0].length > w[1].length));
        assert_eq!(o.events.last().map(|e| e.length), o.length);
    }
}

/// `-l`: the search stops at the first solution of at most that length, and a position in
/// `H_A` (class-A moves only) is solved through phase 1 of length 0.
#[test]
fn max_length_and_short_positions() {
    let (t, ct) = tables();
    // U F BR BL as physical moves are class A
    let scr = parse_moves("U F' BR BL' U'").unwrap();
    let (cube, frame) = scramble_to_cube(&scr);
    let limits = Limits {
        max_length: Some(5),
        ..Limits::default()
    };
    let o = one_thread(|| {
        solve(
            t,
            ct,
            None,
            cube,
            frame,
            &limits,
            Report::Quiet,
            &Default::default(),
        )
    });
    let sol = o.solution.expect("a solution");
    check_solution(&scr, &sol);
    assert!(sol.len() <= 5);
}

/// The phase-2 tables bound the restricted distance (breadth-first search over whole phase-2
/// positions up to 5 moves from solved), and walking down each reaches solved.
#[test]
fn tables_bound_the_restricted_distance() {
    let (_, t) = tables();
    let s = space();
    let start = P2::of(&Cube::solved()).unwrap();
    let mut dist = std::collections::HashMap::from([(start, 0u8)]);
    let mut front = vec![start];
    for d in 1..=5u8 {
        let mut next = vec![];
        for p in &front {
            for m in 0..N_MOVES {
                if let Some(q) = p.apply(s, m) {
                    if let std::collections::hash_map::Entry::Vacant(e) = dist.entry(q) {
                        e.insert(d);
                        next.push(q);
                    }
                }
            }
        }
        front = next;
    }
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(5);
    let all: Vec<(P2, u8)> = dist.into_iter().collect();
    for _ in 0..3000 {
        let (p, d) = all[rng.gen_range(0..all.len())];
        let h = t.exact(&p).unwrap();
        assert!(h.iter().all(|&x| x <= d), "{p:?} {d} {h:?}");
    }
    assert_eq!(t.exact(&start), Some([0; 4]));
    assert!(t.ha_index(s, &start).is_some());
}

/// The relaxed phase-1 table of the small tier's index (built in memory once).
fn relaxed() -> &'static RelaxedAe {
    static R: OnceLock<RelaxedAe> = OnceLock::new();
    R.get_or_init(|| RelaxedAe::generate(&tables().0.ae_ix, false))
}

/// With the relaxed target the solutions solve, a configuration of the phase-2 space is at
/// phase-1 distance 0 (and `H'` among them), and a random position is no farther than from
/// `H'`.
#[test]
fn relaxed_solutions_solve() {
    let (t, ct) = tables();
    let r = relaxed();
    let s = space();
    let f = omelettehedron::fast::fast();
    // every configuration's index is a target (in the small tier's index, which has no edge
    // chirality, some targets also hold positions that are no configuration)
    assert_eq!(r.targets, relaxed_targets(&t.ae_ix));
    let (hits, n) = target_check(&t.ae_ix, &r.targets);
    assert!(hits > 0 && hits <= n);
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(21);
    for &(xa, ep) in s.cfgs.iter().step_by(997) {
        let c = Cube {
            ep,
            xa,
            ..Cube::random(&mut rng)
        };
        assert_eq!(r.exact(&t.ae_ix, &f.of(&c)), Some(0));
    }
    for _ in 0..200 {
        let c = Cube::random(&mut rng);
        let x = f.of(&c);
        assert!(r.exact(&t.ae_ix, &x).unwrap() <= t.exact_fast(&x).0);
    }
    for len in [0, 3, 7, 11] {
        let scr: Vec<PhysMove> = (0..len)
            .map(|_| PhysMove::face(rng.gen_range(0..8), rng.gen()))
            .collect();
        let (cube, frame) = scramble_to_cube(&scr);
        let limits = Limits {
            max_solutions: Some(2),
            time: Some(Duration::from_secs(600)),
            ..Limits::default()
        };
        let o = solve(
            t,
            ct,
            Some(r),
            cube,
            frame,
            &limits,
            Report::Quiet,
            &Default::default(),
        );
        check_solution(&scr, &o.solution.expect("a solution"));
    }
}
