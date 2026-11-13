//! Both solvers on one position at once (`omelettehedron solve combined`).
//!
//! The two-phase solver and the optimal solver run at the same time on thread pools of their
//! own and share one state: the upper bound is the shortest solution either has found, the
//! lower bound the largest either has proven (the optimal solver: each completed depth `D`
//! proves at least `D + 1`; the two-phase solver: each completed phase-1 iteration and target
//! pass, see `search::solve`). Every new bound is printed as it comes, and both solvers stop
//! as soon as the bounds meet: the best solution is then optimal, exactly as when either
//! solver proves it alone. They also help each other: the shared stop flag ends the optimal
//! solver's search as soon as the two-phase solver's upper bound `U` is proven (so it never
//! searches deeper than `U - 1`), and the two-phase solver's target-length passes start from,
//! and leave for, the largest lower bound proven by either.

use crate::cube::Cube;
use crate::notation::{format_moves, Frame, PhysMove};
use crate::optimal::search as opt;
use crate::optimal::tables::Bounds;
use crate::progress::{Bound, Progress, Watch};
use crate::search as tp;
use crate::tables::{AeIndex, Tables};
use rayon::ThreadPool;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Mutex, PoisonError};
use std::time::{Duration, Instant};

/// How the two-phase solver searches in a combined run.
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum TwoPhaseMode {
    /// as under a time limit: every better solution as soon as possible, so the upper bound
    /// falls fast; lower bounds come mostly from the optimal solver
    Improve,
    /// its optimal mode: target-length passes, from the largest lower bound proven by either
    /// solver, prove lower bounds too
    Passes,
}

/// Which solver gave a bound.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Source {
    TwoPhase,
    Optimal,
}

impl Source {
    pub fn name(self) -> &'static str {
        match self {
            Source::TwoPhase => "two-phase",
            Source::Optimal => "optimal",
        }
    }
}

/// Settings of a combined solve.
#[derive(Clone, Debug)]
pub struct Config {
    /// stop after this much time (both solvers)
    pub time: Option<Duration>,
    /// stop at the first solution of at most this many moves (the optimal solver also gives up
    /// after this depth)
    pub max_length: Option<usize>,
    pub mode: TwoPhaseMode,
    /// progress lines of both solvers on stderr, with a heartbeat this often
    pub progress: Option<Duration>,
    /// print every new bound on stdout as it comes
    pub print: bool,
}

/// A bound and where it came from.
#[derive(Clone, Copy, Debug)]
pub struct Found {
    pub length: usize,
    pub source: Source,
    /// seconds since the start
    pub time: f64,
}

/// The result of a combined solve.
#[derive(Debug)]
pub struct Outcome {
    /// the shortest solution found (by either solver)
    pub solution: Option<Vec<PhysMove>>,
    /// its length and source
    pub upper: Option<Found>,
    /// the largest proven lower bound (no solution has fewer moves) and its source
    pub lower: Option<Found>,
    /// the bounds met: the solution is optimal
    pub optimal: bool,
    pub elapsed: f64,
    /// the two solvers' own results
    pub two_phase: tp::Outcome,
    pub optimal_solver: opt::Outcome,
}

/// The shared state of a combined solve.
struct State {
    solution: Option<Vec<PhysMove>>,
    upper: Option<Found>,
    lower: Option<Found>,
    /// a lower bound above the upper bound: the tables or a solver are wrong
    conflict: Option<String>,
}

/// Solve one position with both solvers at once: the two-phase solver with tables `t` on
/// `pool_tp`, the optimal solver with bounds `b` on `pool_opt` (both in the same move set).
///
/// # Errors
///
/// If no thread can be started, or a solver fails (tables that are not distance tables, a
/// lower bound above a solution's length).
pub fn solve<I: AeIndex, B: Bounds>(
    t: &Tables<I>,
    b: &B,
    pool_tp: &ThreadPool,
    pool_opt: &ThreadPool,
    cube: Cube,
    frame: Frame,
    cfg: &Config,
) -> Result<Outcome, String> {
    let start = Instant::now();
    let stop = AtomicBool::new(false);
    // the largest lower bound proven by either solver (the two-phase passes start there)
    let lower = AtomicUsize::new(0);
    let state = Mutex::new(State {
        solution: None,
        upper: None,
        lower: None,
        conflict: None,
    });
    let on_bound = |source: Source, bound: Bound| {
        let time = start.elapsed().as_secs_f64();
        let mut s = state.lock().unwrap_or_else(PoisonError::into_inner);
        match bound {
            Bound::Upper(m) => {
                if s.upper.is_none_or(|u| m.len() < u.length) {
                    s.upper = Some(Found {
                        length: m.len(),
                        source,
                        time,
                    });
                    s.solution = Some(m.to_vec());
                    if cfg.print {
                        println!("{}", bound_line(time, source, bound));
                    }
                    if cfg.max_length.is_some_and(|l| m.len() <= l) {
                        stop.store(true, Ordering::Relaxed);
                    }
                }
            }
            Bound::Lower(l) => {
                if s.lower.is_none_or(|x| l > x.length) {
                    s.lower = Some(Found {
                        length: l,
                        source,
                        time,
                    });
                    lower.fetch_max(l, Ordering::Relaxed);
                    if cfg.print {
                        println!("{}", bound_line(time, source, bound));
                    }
                }
            }
        }
        if let Some(u) = s.upper {
            // (no lower bound yet: 0)
            let l = s.lower.map_or(0, |l| l.length);
            if l >= u.length {
                if let (Some(lf), None) = (s.lower.filter(|_| l > u.length), &s.conflict) {
                    s.conflict = Some(format!(
                        "the {} solver found a solution of {} moves, but the {} solver proved \
                         that none has fewer than {l} (are the tables right?)",
                        u.source.name(),
                        u.length,
                        lf.source.name(),
                    ));
                }
                // the bounds met: the best solution is optimal
                stop.store(true, Ordering::Relaxed);
            }
        }
    };
    let on_tp = |b: Bound| on_bound(Source::TwoPhase, b);
    let on_opt = |b: Bound| on_bound(Source::Optimal, b);
    let progress = |label| {
        cfg.progress.map(|every| Progress {
            every,
            label,
            start,
        })
    };
    let (prog_tp, prog_opt) = (progress("two-phase: "), progress("optimal:   "));
    let watch_tp = Watch {
        progress: prog_tp.as_ref(),
        stop: Some(&stop),
        on_bound: Some(&on_tp),
        lower: Some(&lower),
        improve: cfg.mode == TwoPhaseMode::Improve,
    };
    let watch_opt = Watch {
        progress: prog_opt.as_ref(),
        stop: Some(&stop),
        on_bound: Some(&on_opt),
        lower: None,
        improve: false,
    };
    // The time limit of the two-phase solver is a cancellation deadline, which keeps its
    // optimal mode (a stopping rule would turn the target-length passes off).
    let limits_tp = tp::Limits {
        cancel_after: cfg.time,
        ..tp::Limits::default()
    };
    let limits_opt = opt::Limits {
        time: cfg.time,
        max_length: cfg.max_length,
        ..opt::Limits::default()
    };
    let (o_tp, o_opt) = std::thread::scope(|s| {
        let h = std::thread::Builder::new()
            .name("two-phase".into())
            .spawn_scoped(s, || {
                pool_tp.install(|| {
                    tp::solve_watched(t, cube, frame, &limits_tp, tp::Report::Quiet, &watch_tp)
                })
            })
            .map_err(|e| format!("cannot start a thread: {e}"))?;
        let o_opt =
            pool_opt.install(|| opt::solve_watched(b, cube, frame, &limits_opt, false, &watch_opt));
        if o_opt.stop == opt::Stop::BadTables {
            stop.store(true, Ordering::Relaxed);
        }
        let o_tp = h
            .join()
            .map_err(|_| "the two-phase solver failed".to_string())?;
        Ok::<_, String>((o_tp, o_opt))
    })?;
    if o_opt.stop == opt::Stop::BadTables {
        return Err(
            "the optimal solver's tables are not distance tables for this position (delete \
             them to regenerate)"
                .into(),
        );
    }
    let s = state.into_inner().unwrap_or_else(PoisonError::into_inner);
    if let Some(c) = s.conflict {
        return Err(c);
    }
    let optimal = s
        .upper
        .is_some_and(|u| s.lower.map_or(0, |l| l.length) >= u.length);
    Ok(Outcome {
        solution: s.solution,
        upper: s.upper,
        lower: s.lower,
        optimal,
        elapsed: start.elapsed().as_secs_f64(),
        two_phase: o_tp,
        optimal_solver: o_opt,
    })
}

/// One line of the log of bounds that a solve prints as it runs: the time (seconds since the
/// start), the bound and the solver that found it, and for an upper bound its solution, e.g.
/// `    1.076 s  upper 20  two-phase  D R' BR ...` or `    8.207 s  lower 20  optimal`.
pub fn bound_line(time: f64, source: Source, bound: Bound) -> String {
    match bound {
        Bound::Upper(m) => format!(
            "{time:>9.3} s  upper {:>2}  {:<9}  {}",
            m.len(),
            source.name(),
            format_moves(m)
        ),
        Bound::Lower(l) => format!("{time:>9.3} s  lower {l:>2}  {}", source.name()),
    }
}

/// The line after the solution: its length (`upper`) and how far it is proven (`lower`, the
/// largest lower bound; `optimal`: they met), and the time the solve took.
pub fn result_line(
    upper: Option<Found>,
    lower: Option<Found>,
    optimal: bool,
    elapsed: f64,
) -> String {
    let lower_len = lower.map_or(0, |l| l.length);
    let from = |f: Found| format!("{} at {:.3} s", f.source.name(), f.time);
    match (upper, optimal) {
        (Some(u), true) if u.length == 0 => "0 moves: already solved".into(),
        (Some(u), true) => format!(
            "{} moves, optimal (solution: {}; lower bound: {}); {elapsed:.3} s",
            u.length,
            from(u),
            lower.map_or_else(String::new, from),
        ),
        (Some(u), false) => format!(
            "{} moves, not proven optimal: no solution has fewer than {lower_len} moves; \
             {elapsed:.3} s",
            u.length
        ),
        (None, _) => {
            format!("no solution found: none has fewer than {lower_len} moves; {elapsed:.3} s")
        }
    }
}

/// Print the result of a combined solve: the solution, then its length and how far it is
/// proven ([`result_line`]).
pub fn print_outcome(o: &Outcome) {
    if let Some(s) = &o.solution {
        println!("{}", format_moves(s));
    }
    println!("{}", result_line(o.upper, o.lower, o.optimal, o.elapsed));
}
