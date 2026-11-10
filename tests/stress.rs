//! Randomised end-to-end checks: every move set, random stopping rules and thread counts.
//!
//! Each solution must solve its scramble on the sticker model, have the reported length and use
//! only the allowed move families; improvement events must shorten strictly; optimal lengths
//! must not depend on the thread count. (Small tier, built in memory once per metric.)

use omelettehedron::cube::*;
use omelettehedron::notation::*;
use omelettehedron::search::*;
use omelettehedron::tables::*;
use rand::{Rng, SeedableRng};
use std::sync::OnceLock;
use std::time::Duration;

const SETS: [MoveSet; 4] = [
    MoveSet::FACE,
    MoveSet {
        wide: true,
        slice: false,
    },
    MoveSet {
        wide: false,
        slice: true,
    },
    MoveSet {
        wide: true,
        slice: true,
    },
];

/// Tables per move set (`wide` only changes the output, so it shares nothing but the flag).
fn tables(moves: MoveSet) -> &'static Tables<AeSym> {
    static T: [OnceLock<Tables<AeSym>>; 4] = [const { OnceLock::new() }; 4];
    let i = usize::from(moves.wide) + 2 * usize::from(moves.slice);
    T[i].get_or_init(|| Tables::generate(AeSym::new(false), moves, false))
}

fn pool(threads: usize) -> &'static rayon::ThreadPool {
    static P: [OnceLock<rayon::ThreadPool>; 3] = [const { OnceLock::new() }; 3];
    let i = [1, 3, 8].iter().position(|&t| t == threads).unwrap();
    P[i].get_or_init(|| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
    })
}

/// Random moves of every layer kind (scrambles may always contain wide and slice turns).
fn scramble(rng: &mut impl Rng, n: usize) -> Vec<PhysMove> {
    (0..n)
        .map(|_| PhysMove {
            face: rng.gen_range(0..8),
            layer: [Layer::Face, Layer::Wide, Layer::Slice][rng.gen_range(0..3)],
            cw: rng.gen(),
        })
        .collect()
}

/// The checks every outcome must pass.
fn check(moves: MoveSet, scr: &[PhysMove], fr: Frame, o: &Outcome, what: &str) {
    let sol = o.solution.as_ref().expect("a solution");
    let ctx = format!("{what}: {} / {}", format_moves(scr), format_moves(sol));
    let st = apply_physical(&apply_physical(&solved_stickers(defs()), scr), sol);
    assert!(is_solved_physical(&st), "does not solve: {ctx}");
    assert_eq!(o.length, Some(sol.len()), "{ctx}");
    for pm in sol {
        let ok = match pm.layer {
            Layer::Face => true,
            Layer::Wide => moves.wide,
            Layer::Slice => moves.slice,
        };
        assert!(ok, "move family not allowed ({moves:?}): {ctx}");
        if moves.wide {
            // wide output keeps the corner the scramble left in front where it is
            assert_eq!(
                fr.to_internal(*pm).1,
                fr,
                "wide output moved the frame: {ctx}"
            );
        }
    }
    for w in o.events.windows(2) {
        assert!(
            w[1].length < w[0].length,
            "events not strictly shorter: {ctx}"
        );
    }
    if let Some(e) = o.events.last() {
        assert_eq!(e.length, sol.len(), "last event is not the result: {ctx}");
    }
}

#[test]
#[ignore = "slow; cargo test --release -- --ignored"]
fn random_limits_threads_and_move_sets() {
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(2611);
    for moves in SETS {
        let t = tables(moves);
        for k in 0..12 {
            let scr = scramble(&mut rng, 40);
            let (c, fr) = scramble_to_cube(&scr);
            let threads = [1, 3, 8][k % 3];
            let limits = match k % 4 {
                0 => Limits {
                    time: Some(Duration::from_millis(rng.gen_range(10..60))),
                    ..Limits::default()
                },
                1 => Limits {
                    max_solutions: Some(rng.gen_range(1..40)),
                    ..Limits::default()
                },
                2 => Limits {
                    max_solutions: Some(1),
                    time: Some(Duration::from_millis(200)),
                    ..Limits::default()
                },
                _ => Limits {
                    // a length cap that the first solution (at most 30 moves here) may not meet
                    max_length: Some(rng.gen_range(22..30)),
                    time: Some(Duration::from_millis(300)),
                    ..Limits::default()
                },
            };
            let o = pool(threads).install(|| solve(t, c, fr, &limits, Report::Quiet));
            let what = format!("{moves:?} {threads}T {limits:?}");
            if o.solution.is_none() {
                // only a time limit may end a search without a solution (the small tier can need
                // seconds for the first one on one thread)
                let limit = limits.time.map(|d| d.as_secs_f64());
                assert!(
                    limit.is_some_and(|l| o.elapsed >= 0.9 * l),
                    "no solution and no time-out: {what} ({} s)",
                    o.elapsed
                );
                assert!(o.length.is_none() && o.events.is_empty(), "{what}");
                continue;
            }
            check(moves, &scr, fr, &o, &what);
            if let Some(n) = limits.max_solutions {
                if limits.time.is_none() {
                    assert!(o.solutions_found >= n, "{what}");
                }
            }
            if let (Some(l), Some(len)) = (limits.max_length, o.length) {
                // the cap stops the search once met; if it was not met the time limit ran out
                assert!(len <= l || o.elapsed >= 0.25, "{what}: {len} > {l}");
            }
        }
    }
}

#[test]
#[ignore = "slow; cargo test --release -- --ignored"]
fn optimal_lengths_do_not_depend_on_threads() {
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(2612);
    for moves in SETS {
        let t = tables(moves);
        for _ in 0..6 {
            let scr = scramble(&mut rng, 7);
            let (c, fr) = scramble_to_cube(&scr);
            let lens: Vec<usize> = [1, 3, 8]
                .iter()
                .map(|&th| {
                    let o = pool(th).install(|| solve(t, c, fr, &Limits::default(), Report::Quiet));
                    check(moves, &scr, fr, &o, &format!("{moves:?} {th}T optimal"));
                    assert!(o.optimal, "{moves:?} {th}T: not proven optimal");
                    o.length.unwrap()
                })
                .collect();
            assert!(
                lens.windows(2).all(|w| w[0] == w[1]),
                "{moves:?}: optimal lengths differ by thread count: {lens:?} ({})",
                format_moves(&scr)
            );
        }
    }
}
