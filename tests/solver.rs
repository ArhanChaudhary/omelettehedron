//! End-to-end correctness tests (use the small table tier, built in memory once).

use omelettehedron::coord::*;
use omelettehedron::cube::*;
use omelettehedron::notation::*;
use omelettehedron::search::*;
use omelettehedron::tables::*;
use rand::{Rng, SeedableRng};
use std::collections::HashMap;
use std::sync::OnceLock;
use std::time::Duration;

fn tables() -> &'static Tables<AeSym> {
    static T: OnceLock<Tables<AeSym>> = OnceLock::new();
    T.get_or_init(|| Tables::generate(AeSym::new(false), false))
}

fn random_scramble(rng: &mut impl Rng, n: usize) -> Vec<PhysMove> {
    (0..n)
        .map(|_| PhysMove {
            face: rng.gen_range(0..8),
            cw: rng.gen(),
        })
        .collect()
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

/// Exact distances up to 2*R by meet in the middle (fixed frame = solved up to rotation).
struct Mitm {
    ball: HashMap<Cube, u8>,
}
const R: usize = 5;
impl Mitm {
    fn new() -> Mitm {
        let mut ball = HashMap::new();
        ball.insert(Cube::solved(), 0u8);
        let mut fr = vec![Cube::solved()];
        for d in 1..=R as u8 {
            let mut nf = vec![];
            for c in &fr {
                for m in 0..N_MOVES {
                    let n = c.apply(m);
                    ball.entry(n).or_insert_with(|| {
                        nf.push(n);
                        d
                    });
                }
            }
            fr = nf;
        }
        Mitm { ball }
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
                for m in 0..N_MOVES {
                    if allowed_after(last, m) {
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
    let mitm = Mitm::new();
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

#[test]
fn phase2_solves_class_a_positions() {
    let p2 = &tables().p2;
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(15);
    let amoves: Vec<usize> = (0..N_MOVES).filter(|&m| is_a_move(m)).collect();
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
    let mitm = Mitm::new();
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
    }
}
