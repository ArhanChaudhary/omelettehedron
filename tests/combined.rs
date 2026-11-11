//! The combined solve (both solvers at once, `omelettehedron combined`) and the progress
//! option: a combined answer marked optimal has the length both solvers prove separately, and
//! progress lines do not change the search. Small tiers of both solvers, built in memory once.

use omelettehedron::combined::{self, Config, TwoPhaseMode};
use omelettehedron::cube::*;
use omelettehedron::notation::*;
use omelettehedron::optimal::search as opt;
use omelettehedron::optimal::tables::{OptTables, OptTier};
use omelettehedron::progress::{Progress, Watch};
use omelettehedron::search as tp;
use omelettehedron::tables::{AeSym, Tables};
use rand::{Rng, SeedableRng};
use std::sync::OnceLock;
use std::time::Duration;

fn two_phase() -> &'static Tables<AeSym> {
    static T: OnceLock<Tables<AeSym>> = OnceLock::new();
    T.get_or_init(|| Tables::generate(AeSym::new(false), MoveSet::FACE, false))
}

fn optimal() -> &'static OptTables {
    static T: OnceLock<OptTables> = OnceLock::new();
    T.get_or_init(|| OptTables::generate(OptTier::Small, MoveSet::FACE, false).unwrap())
}

fn pool(threads: usize) -> rayon::ThreadPool {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap()
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

fn config(mode: TwoPhaseMode, time: Option<Duration>) -> Config {
    Config {
        time,
        mode,
        progress: None,
        print: false,
    }
}

/// On scrambles up to 15 moves, in both two-phase modes and with uneven thread splits, the
/// combined answer is proven optimal with the length that each solver proves alone.
#[test]
fn combined_length_matches_both_solvers() {
    let t = two_phase();
    let b = optimal().small().unwrap();
    let p8 = pool(8);
    let (pa, pb, pc) = (pool(4), pool(4), pool(1));
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(2027);
    for (i, len) in [0usize, 1, 3, 8, 10, 11, 12, 12, 13, 13, 14, 14, 15]
        .into_iter()
        .enumerate()
    {
        let scr = random_scramble(&mut rng, len);
        let (c, fr) = scramble_to_cube(&scr);
        let want_tp = p8.install(|| tp::solve(t, c, fr, &tp::Limits::default(), tp::Report::Quiet));
        assert!(want_tp.optimal);
        let want_opt = p8.install(|| opt::solve(&b, c, fr, &opt::Limits::default(), false));
        assert_eq!(want_opt.length, want_tp.length, "{}", format_moves(&scr));
        for mode in [TwoPhaseMode::Improve, TwoPhaseMode::Passes] {
            // (4 + 4 threads, and 1 + 4 / 4 + 1 in turn)
            let (p_tp, p_opt) = match i % 3 {
                0 => (&pa, &pb),
                1 => (&pc, &pb),
                _ => (&pa, &pc),
            };
            let o = combined::solve(t, &b, p_tp, p_opt, c, fr, &config(mode, None)).unwrap();
            assert!(o.optimal, "{mode:?} {}", format_moves(&scr));
            assert_eq!(
                o.upper.map(|u| u.length),
                want_tp.length,
                "{mode:?} {}",
                format_moves(&scr)
            );
            assert_eq!(o.lower.map_or(Some(0), |l| Some(l.length)), want_tp.length);
            let sol = o.solution.unwrap();
            assert_eq!(Some(sol.len()), want_tp.length);
            check_solution(&scr, &sol);
        }
    }
}

/// Under a time limit the bounds stay consistent: the lower bound never exceeds the best
/// length, and an answer is marked optimal only when they meet.
#[test]
fn combined_time_limit_keeps_the_bounds_consistent() {
    let t = two_phase();
    let b = optimal().small().unwrap();
    let (pa, pb) = (pool(2), pool(2));
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(2028);
    for ms in [0u64, 1, 5, 20, 50] {
        let scr = random_scramble(&mut rng, 30);
        let (c, fr) = scramble_to_cube(&scr);
        for mode in [TwoPhaseMode::Improve, TwoPhaseMode::Passes] {
            let cfg = config(mode, Some(Duration::from_millis(ms)));
            let o = combined::solve(t, &b, &pa, &pb, c, fr, &cfg).unwrap();
            let lower = o.lower.map_or(0, |l| l.length);
            if let Some(u) = o.upper {
                assert!(lower <= u.length);
                assert_eq!(o.optimal, lower == u.length);
                let sol = o.solution.unwrap();
                assert_eq!(sol.len(), u.length);
                check_solution(&scr, &sol);
            } else {
                assert!(!o.optimal);
            }
            // the solvers' own answers agree with the combined bounds (every bound they prove
            // is reported)
            assert!(o.optimal_solver.lower_bound <= lower);
            if let Some(l) = o.two_phase.length {
                assert!(lower <= l);
            }
        }
    }
}

/// Progress lines (with a heartbeat thread) change nothing in the search: on one thread the
/// node counts are identical with and without them.
#[test]
fn progress_does_not_change_the_search() {
    let t = two_phase();
    let b = optimal().small().unwrap();
    let p1 = pool(1);
    let prog = Progress::new(Duration::from_millis(1), "test: ");
    let watch = Watch {
        progress: Some(&prog),
        ..Watch::default()
    };
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(2029);
    for len in [9usize, 12, 14] {
        let scr = random_scramble(&mut rng, len);
        let (c, fr) = scramble_to_cube(&scr);
        let limits = tp::Limits::default();
        let a = p1.install(|| tp::solve(t, c, fr, &limits, tp::Report::Quiet));
        let w = p1.install(|| tp::solve_watched(t, c, fr, &limits, tp::Report::Quiet, &watch));
        assert_eq!(
            (a.length, a.nodes, a.optimal),
            (w.length, w.nodes, w.optimal)
        );
        let limits = opt::Limits::default();
        let a = p1.install(|| opt::solve(&b, c, fr, &limits, false));
        let w = p1.install(|| opt::solve_watched(&b, c, fr, &limits, false, &watch));
        assert_eq!((a.length, a.nodes), (w.length, w.nodes));
    }
}

/// A lower bound proven elsewhere (`Watch::lower`), given at the start or raised while the
/// search runs (so that passes start at it, or are left for it), never changes the optimal
/// length of the two-phase solver's optimal mode.
#[test]
fn external_lower_bound_keeps_two_phase_optimal() {
    use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
    let t = two_phase();
    let p4 = pool(4);
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(2030);
    for len in [11usize, 13, 14, 15, 16] {
        let scr = random_scramble(&mut rng, len);
        let (c, fr) = scramble_to_cube(&scr);
        let limits = tp::Limits::default();
        let want = p4.install(|| tp::solve(t, c, fr, &limits, tp::Report::Quiet));
        assert!(want.optimal);
        let opt_len = want.length.unwrap();
        // given at the start: the true optimum, or below it
        for lb in [opt_len.saturating_sub(3), opt_len - 1, opt_len] {
            let lower = AtomicUsize::new(lb);
            let watch = Watch {
                lower: Some(&lower),
                ..Watch::default()
            };
            let o = p4.install(|| tp::solve_watched(t, c, fr, &limits, tp::Report::Quiet, &watch));
            assert!(o.optimal);
            assert_eq!(o.length, Some(opt_len), "{lb} {}", format_moves(&scr));
            check_solution(&scr, o.solution.as_ref().unwrap());
        }
        // raised step by step while the search runs, up to the true optimum
        let lower = AtomicUsize::new(0);
        let done = AtomicBool::new(false);
        let watch = Watch {
            lower: Some(&lower),
            ..Watch::default()
        };
        let o = std::thread::scope(|s| {
            s.spawn(|| {
                for l in opt_len.saturating_sub(4)..=opt_len {
                    if done.load(Ordering::Relaxed) {
                        break;
                    }
                    lower.store(l, Ordering::Relaxed);
                    std::thread::sleep(Duration::from_micros(300));
                }
            });
            let o = p4.install(|| tp::solve_watched(t, c, fr, &limits, tp::Report::Quiet, &watch));
            done.store(true, Ordering::Relaxed);
            o
        });
        assert!(o.optimal);
        assert_eq!(o.length, Some(opt_len), "raised: {}", format_moves(&scr));
    }
}
