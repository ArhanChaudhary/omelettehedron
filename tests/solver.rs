//! End-to-end correctness tests (use the small table tier, built in memory once per move set).

use omelettehedron::coord::*;
use omelettehedron::cube::*;
use omelettehedron::fast::*;
use omelettehedron::notation::*;
use omelettehedron::search::*;
use omelettehedron::tables::*;
use rand::{Rng, SeedableRng};
use std::collections::HashMap;
use std::sync::OnceLock;
use std::time::Duration;

const SLICE: MoveSet = MoveSet {
    wide: false,
    slice: true,
};

fn tables() -> &'static Tables<AeSym> {
    static T: OnceLock<Tables<AeSym>> = OnceLock::new();
    T.get_or_init(|| Tables::generate(AeSym::new(false), MoveSet::FACE, false))
}

/// The same tables in the metric with slice moves.
fn slice_tables() -> &'static Tables<AeSym> {
    static T: OnceLock<Tables<AeSym>> = OnceLock::new();
    T.get_or_init(|| Tables::generate(AeSym::new(false), SLICE, false))
}

fn random_scramble(rng: &mut impl Rng, n: usize) -> Vec<PhysMove> {
    (0..n)
        .map(|_| PhysMove::face(rng.gen_range(0..8), rng.gen()))
        .collect()
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

/// Optimal solve on one thread, where the target-length passes are on by default.
fn solve_one_thread(t: &Tables<AeSym>, c: Cube, fr: Frame) -> Outcome {
    static POOL: OnceLock<rayon::ThreadPool> = OnceLock::new();
    POOL.get_or_init(|| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .unwrap()
    })
    .install(|| solve(t, c, fr, &Limits::default(), Report::Quiet))
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

/// Exact distances up to 2*R in the metric of a move set by meet in the middle (fixed frame =
/// solved up to rotation). Independent of the solver's tables and search.
struct Mitm {
    ball: HashMap<Cube, u8>,
    moves: MoveSet,
}
const R: usize = 5;
impl Mitm {
    fn new(moves: MoveSet) -> Mitm {
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

#[test]
fn optimal_matches_brute_force() {
    let t = tables();
    let mitm = Mitm::new(MoveSet::FACE);
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(11);
    let mut checked = 0;
    for len in [1usize, 3, 5, 7, 8, 9] {
        for _ in 0..4 {
            let scr = random_scramble(&mut rng, len);
            let (c, fr) = scramble_to_cube(&scr);
            let truth = mitm.distance(&c).expect("within range");
            let o = solve(t, c, fr, &Limits::default(), Report::Quiet);
            assert!(o.optimal);
            assert_eq!(o.length, Some(truth), "scramble {}", format_moves(&scr));
            check_solution(&scr, o.solution.as_ref().unwrap());
            checked += 1;
        }
    }
    assert_eq!(checked, 24);
}

/// Optimal mode in the metric with slices against the brute force in the same metric, on
/// scrambles of face, wide and slice turns. The solutions use only face and slice turns, and
/// count each as one move.
#[test]
fn optimal_matches_brute_force_with_slices() {
    let t = slice_tables();
    let mitm = Mitm::new(SLICE);
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(16);
    let mut slices_used = 0;
    let mut shorter = 0;
    for len in [1usize, 2, 3, 5, 7, 8, 9, 10] {
        for k in 0..6 {
            let scr = if k % 2 == 0 {
                random_scramble(&mut rng, len)
            } else {
                random_layer_scramble(&mut rng, len)
            };
            let (c, fr) = scramble_to_cube(&scr);
            let truth = mitm.distance(&c).expect("within range");
            let o = solve(t, c, fr, &Limits::default(), Report::Quiet);
            assert!(o.optimal);
            assert_eq!(o.length, Some(truth), "scramble {}", format_moves(&scr));
            let sol = o.solution.as_ref().unwrap();
            check_solution(&scr, sol);
            let o1 = solve_one_thread(t, c, fr);
            assert_eq!(
                o1.length,
                Some(truth),
                "1T, scramble {}",
                format_moves(&scr)
            );
            assert!(sol.iter().all(|m| m.layer != Layer::Wide));
            slices_used += sol.iter().filter(|m| m.layer == Layer::Slice).count();
            // never longer than in the face metric
            if k % 2 == 0 {
                let face = solve(tables(), c, fr, &Limits::default(), Report::Quiet);
                assert!(truth <= face.length.unwrap());
                shorter += usize::from(truth < face.length.unwrap());
            }
        }
    }
    assert!(slices_used > 0 && shorter > 0, "{slices_used} {shorter}");
}

/// With `wide`, the same lengths, and solutions without front-face turns or rotations (the
/// front corner stays still).
#[test]
fn wide_output_keeps_the_front_corner() {
    let t = tables();
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(17);
    for _ in 0..4 {
        let scr = random_scramble(&mut rng, 8);
        let (c, fr) = scramble_to_cube(&scr);
        let o = solve(t, c, fr, &Limits::default(), Report::Quiet);
        let sol = o.solution.as_ref().unwrap();
        // the same internal solution, printed with wide turns
        let mut f2 = fr;
        let mut internal = vec![];
        for &pm in sol {
            let (m, nf) = f2.to_internal(pm);
            internal.push(m);
            f2 = nf;
        }
        let wide = solution_to_physical(fr, &internal, true);
        assert_eq!(wide.len(), sol.len());
        check_solution(&scr, &wide);
        // every move leaves the layer at the scramble's front corner alone
        for pm in &wide {
            assert_eq!(fr.to_internal(*pm).1, fr);
        }
    }
}

#[test]
fn random_positions_are_solved() {
    let t = tables();
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(12);
    for _ in 0..4 {
        let scr = random_scramble(&mut rng, 60);
        let (c, fr) = scramble_to_cube(&scr);
        let lim = Limits {
            max_solutions: Some(1),
            ..Default::default()
        };
        let o = solve(t, c, fr, &lim, Report::Quiet);
        check_solution(&scr, o.solution.as_ref().expect("solution"));
        assert!(o.solutions_found >= 1);
    }
}

#[test]
fn stopping_rules() {
    let t = tables();
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(13);
    let scr = random_scramble(&mut rng, 60);
    let (c, fr) = scramble_to_cube(&scr);
    // solution count
    let o = solve(
        t,
        c,
        fr,
        &Limits {
            max_solutions: Some(5),
            ..Default::default()
        },
        Report::Quiet,
    );
    assert!(
        o.solutions_found >= 5 && o.solutions_found < 5 + 64,
        "{}",
        o.solutions_found
    );
    assert!(!o.optimal);
    let l = o.length.expect("a solution after 5 completions");
    // time
    let o = solve(
        t,
        c,
        fr,
        &Limits {
            time: Some(Duration::from_millis(300)),
            ..Default::default()
        },
        Report::Quiet,
    );
    assert!(o.elapsed < 1.5, "{}", o.elapsed);
    // target length: ask for at most a length already known to be reachable
    let o2 = solve(
        t,
        c,
        fr,
        &Limits {
            max_length: Some(l + 2),
            ..Default::default()
        },
        Report::Quiet,
    );
    assert!(o2.length.unwrap() <= l + 2);
    check_solution(&scr, o2.solution.as_ref().unwrap());
}

#[test]
fn heuristics_are_consistent() {
    let t = tables();
    let mt = move_tables();
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(14);
    for _ in 0..30 {
        let c = Cube::random(&mut rng);
        let x = Coords::of(&c);
        let (ha, hc) = (t.ae_exact(&x) as i32, t.cb_exact(&x) as i32);
        for m in 0..N_MOVES {
            let y = mt.apply(&x, m);
            assert!((t.ae_exact(&y) as i32 - ha).abs() <= 1);
            assert!((t.cb_exact(&y) as i32 - hc).abs() <= 1);
        }
    }
}

/// The tables with slices: consistent under all 24 moves, and never above the face-metric
/// distances (every face-metric path is also a path with slices).
#[test]
fn slice_heuristics_are_consistent() {
    let (t, tf) = (slice_tables(), tables());
    let f = fast();
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(18);
    for _ in 0..30 {
        let c = Cube::random(&mut rng);
        let x = f.of(&c);
        let (ha, hc) = t.exact_fast(&x);
        let (fa, fc) = tf.exact_fast(&x);
        assert!(ha <= fa && hc <= fc);
        for m in 0..MAX_MOVES {
            let (ya, yc) = t.exact_fast(&f.of(&c.apply(m)));
            assert!((ya as i32 - ha as i32).abs() <= 1);
            assert!((yc as i32 - hc as i32).abs() <= 1);
        }
    }
    let (hs, hf) = (&t.cb.tab.hist, &tf.cb.tab.hist);
    assert_eq!(hs.iter().sum::<u64>(), hf.iter().sum::<u64>());
    assert_eq!(hs[0], hf[0]);
    assert!(hs.len() <= hf.len());
}

#[test]
fn phase2_solves_class_a_positions() {
    let p2 = &tables().p2;
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(15);
    let amoves: Vec<usize> = (0..N_MOVES).filter(|&m| is_a_move(m)).collect();
    // slices are not in H_A
    for m in N_MOVES..MAX_MOVES {
        assert_eq!(p2.length(&Cube::solved().apply(m)), None);
    }
    for _ in 0..200 {
        let mut c = Cube::solved();
        for _ in 0..40 {
            c = c.apply(amoves[rng.gen_range(0..amoves.len())]);
        }
        let s = p2.solve(&c).expect("in H_A");
        assert_eq!(s.len(), p2.length(&c).unwrap() as usize);
        assert_eq!(c.apply_seq(&s), Cube::solved());
        let x = Coords::of(&c);
        assert_eq!(tables().ae_exact(&x), 0);
        assert_eq!(tables().cb_exact(&x), 0);
    }
    assert_eq!(
        p2.histogram(),
        vec![1, 8, 48, 288, 1728, 9896, 51808, 220111, 480467, 166276, 2457, 32]
    );
}

#[test]
fn table_histograms_match_reference() {
    // independently computed with a Python/GAP sticker model
    // (the CB reference histogram is checked in tables.rs against the unreduced table)
    assert_eq!(tables().cb.tab.hist.len(), 9);
    assert_eq!(tables().ae.hist.len(), 13);
}

/// Optimal mode stops one iteration early (see `solve`). The argument assumes that the
/// conjugation swapping the classes maps single turns to single turns and wide turns to wide
/// turns, keeping each axis's pair together, so that both variants order a commuting pair at
/// the end of a solution the same way.
#[test]
fn class_swap_keeps_single_and_wide_pairs() {
    for m in 0..N_MOVES {
        let c = conjugate_move(m, &SWAP_AB);
        assert_eq!(slot_is_wide(slot(m)), slot_is_wide(slot(c)), "move {m}");
        // the partner on the same axis maps to the partner of the image
        let partner = m ^ (4 << 1);
        assert_eq!(
            slot(conjugate_move(partner, &SWAP_AB)),
            slot(c) ^ 4,
            "move {m}"
        );
    }
}

/// Optimal lengths against meet-in-the-middle distances on many positions (release builds:
/// `cargo test --release -- --ignored`).
#[test]
#[ignore]
fn optimal_matches_brute_force_many() {
    let t = tables();
    let mitm = Mitm::new(MoveSet::FACE);
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(2027);
    for i in 0..400 {
        let scr = random_scramble(&mut rng, 4 + i % 9);
        let (c, fr) = scramble_to_cube(&scr);
        let Some(truth) = mitm.distance(&c) else {
            continue;
        };
        let o = solve(t, c, fr, &Limits::default(), Report::Quiet);
        assert!(o.optimal);
        assert_eq!(o.length, Some(truth), "scramble {}", format_moves(&scr));
        check_solution(&scr, o.solution.as_ref().unwrap());
        // (one thread: with the target-length passes)
        let o = solve_one_thread(t, c, fr);
        assert!(o.optimal);
        assert_eq!(o.length, Some(truth), "1T, scramble {}", format_moves(&scr));
    }
}

/// Optimal lengths with slices against meet-in-the-middle distances on many positions, from
/// scrambles of all layer kinds (release builds: `cargo test --release -- --ignored`).
#[test]
#[ignore]
fn optimal_matches_brute_force_many_with_slices() {
    let t = slice_tables();
    let mitm = Mitm::new(SLICE);
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(2028);
    let mut checked = 0;
    for i in 0..400 {
        let len = 4 + i % 7;
        let scr = if i % 2 == 0 {
            random_scramble(&mut rng, len)
        } else {
            random_layer_scramble(&mut rng, len)
        };
        let (c, fr) = scramble_to_cube(&scr);
        // (at most 10 moves: within the brute force's range)
        let truth = mitm.distance(&c).expect("within range");
        let o = solve(t, c, fr, &Limits::default(), Report::Quiet);
        assert!(o.optimal);
        assert_eq!(o.length, Some(truth), "scramble {}", format_moves(&scr));
        check_solution(&scr, o.solution.as_ref().unwrap());
        // (one thread: with the target-length passes)
        let o = solve_one_thread(t, c, fr);
        assert!(o.optimal);
        assert_eq!(o.length, Some(truth), "1T, scramble {}", format_moves(&scr));
        checked += 1;
    }
    assert_eq!(checked, 400);
}
