//! What a solve prints, the same for every method (`solve two-phase`, `optimal`, `combined`):
//!
//! * on stdout, each new bound as it is found ([`combined::bound_line`]): every better solution
//!   (an upper bound on the optimal length) and every larger lower bound proven, with the time
//!   and the solver;
//! * when the search ends, the solution (with `-n`, the optimal solver's solutions), then its
//!   length and how far it is proven ([`combined::result_line`]), and for a pattern the
//!   scramble that makes it.
//!
//! With `--quiet`, only the solution (or `no solution found`). A solve on its own thread
//! (`--batch`) keeps its text and the text is printed in input order.

use super::input::{Job, PatternInput};
use omelettehedron::combined::{self, Found, Source};
use omelettehedron::notation::{format_moves, solves_pure, PhysMove};
use omelettehedron::pattern;
use omelettehedron::progress::Bound;
use std::fmt::Write as _;
use std::sync::{Mutex, PoisonError};
use std::time::Instant;

/// Where the text of one position's solve goes: printed at once, or kept (`--batch`).
pub struct Text {
    keep: bool,
    /// (stdout, stderr), when kept
    kept: Mutex<(String, String)>,
}

impl Text {
    /// Text printed as it comes.
    pub fn printed() -> Text {
        Text {
            keep: false,
            kept: Mutex::new((String::new(), String::new())),
        }
    }

    /// Text kept until [`Text::take`].
    pub fn kept() -> Text {
        Text {
            keep: true,
            kept: Mutex::new((String::new(), String::new())),
        }
    }

    /// A line of stdout.
    pub fn out(&self, line: &str) {
        if self.keep {
            let mut k = self.kept.lock().unwrap_or_else(PoisonError::into_inner);
            k.0 += line;
            k.0.push('\n');
        } else {
            println!("{line}");
        }
    }

    /// A line of stderr.
    pub fn err(&self, line: &str) {
        if self.keep {
            let mut k = self.kept.lock().unwrap_or_else(PoisonError::into_inner);
            k.1 += line;
            k.1.push('\n');
        } else {
            eprintln!("{line}");
        }
    }

    /// The kept text: (stdout, stderr).
    pub fn take(self) -> (String, String) {
        self.kept
            .into_inner()
            .unwrap_or_else(PoisonError::into_inner)
    }
}

/// How much a solve prints.
#[allow(clippy::struct_excessive_bools)] // independent switches of the output
#[derive(Clone, Copy, Debug)]
pub struct Detail {
    /// only the solution
    pub quiet: bool,
    /// whether the solution also solves the pure puzzle (`-v`)
    pub verbose: bool,
    /// the solver's statistics (`-v`, `--progress`)
    pub stats: bool,
    /// a `# label` line before each position
    pub header: bool,
}

/// The bounds of one solver's search as they come: the best solution and the largest lower
/// bound, each with its time, printed as found unless quiet (the bookkeeping that
/// [`combined::solve`] does for two solvers).
pub struct Tracker<'a> {
    start: Instant,
    source: Source,
    text: &'a Text,
    quiet: bool,
    state: Mutex<Bounds>,
}

/// What a [`Tracker`] holds.
#[derive(Default)]
struct Bounds {
    solution: Option<Vec<PhysMove>>,
    upper: Option<Found>,
    lower: Option<Found>,
}

/// The result of a solve, as printed: the solutions (the best first), the bounds and whether
/// they met.
pub struct Summary {
    pub solutions: Vec<Vec<PhysMove>>,
    pub upper: Option<Found>,
    pub lower: Option<Found>,
    pub optimal: bool,
    pub elapsed: f64,
}

impl<'a> Tracker<'a> {
    /// Bounds of `source`'s search starting now, printed to `text` unless `quiet`.
    pub fn new(source: Source, text: &'a Text, quiet: bool) -> Tracker<'a> {
        Tracker {
            start: Instant::now(),
            source,
            text,
            quiet,
            state: Mutex::new(Bounds::default()),
        }
    }

    /// A bound found by the search (the callback of `Watch::on_bound`).
    pub fn bound(&self, bound: Bound) {
        let time = self.start.elapsed().as_secs_f64();
        let mut s = self.state.lock().unwrap_or_else(PoisonError::into_inner);
        let found = |length| {
            Some(Found {
                length,
                source: self.source,
                time,
            })
        };
        let better = match bound {
            Bound::Upper(m) if s.upper.is_none_or(|u| m.len() < u.length) => {
                s.upper = found(m.len());
                s.solution = Some(m.to_vec());
                true
            }
            Bound::Lower(l) if s.lower.is_none_or(|x| l > x.length) => {
                s.lower = found(l);
                true
            }
            _ => false,
        };
        if better && !self.quiet {
            self.text
                .out(&combined::bound_line(time, self.source, bound));
        }
    }

    /// The result: the best solution, or `solutions` if given (the optimal solver's, all of
    /// the optimal length), and the bounds (optimal if they met).
    pub fn finish(self, solutions: Vec<Vec<PhysMove>>) -> Summary {
        let elapsed = self.start.elapsed().as_secs_f64();
        let s = self
            .state
            .into_inner()
            .unwrap_or_else(PoisonError::into_inner);
        let solutions = if solutions.is_empty() {
            s.solution.into_iter().collect()
        } else {
            solutions
        };
        let optimal = s
            .upper
            .is_some_and(|u| s.lower.map_or(0, |l| l.length) >= u.length);
        Summary {
            solutions,
            upper: s.upper,
            lower: s.lower,
            optimal,
            elapsed,
        }
    }
}

/// Print the result of a position's solve: the solutions, the result line and, for a pattern,
/// the scramble that makes it; with `quiet`, the best solution alone (or `no solution found`).
/// The solver's statistics `stats` and, for a scramble, whether the solution also solves the
/// pure puzzle go to stderr if `d` says so.
pub fn print_result(text: &Text, job: &Job, r: &Summary, d: Detail, stats: &[String]) {
    if d.quiet {
        match r.solutions.first() {
            Some(s) => text.out(&format_moves(s)),
            None => text.out("no solution found"),
        }
        return;
    }
    for s in &r.solutions {
        text.out(&format_moves(s));
    }
    text.out(&combined::result_line(
        r.upper, r.lower, r.optimal, r.elapsed,
    ));
    if let (Some(pi), Some(sol)) = (&job.pattern, r.solutions.first()) {
        pattern_scramble(text, pi, sol);
    }
    if d.stats {
        for line in stats {
            text.err(line);
        }
    }
    if d.verbose {
        // the solver works up to centre colours; say whether the solution also restores
        // every centre triangle to its own place
        if let (Some(scr), Some(sol)) = (&job.moves, r.solutions.first()) {
            text.err(&format!(
                "pure FTO (every centre unique): {}",
                if solves_pure(scr, sol) {
                    "also solved"
                } else {
                    "not solved (centres of a colour exchanged)"
                }
            ));
        }
    }
}

/// Print the scramble that makes a solved pattern (and check the solution).
fn pattern_scramble(text: &Text, pi: &PatternInput, solution: &[PhysMove]) {
    let Some(sc) = pattern::scramble(&pi.pattern, solution, pi.usual.as_ref()) else {
        text.err("error: the solution does not solve the pattern (please report this pattern)");
        return;
    };
    let mut out = format!(
        "scramble from solved with {} on U and {} on F: {}",
        sc.start[0],
        sc.start[2],
        format_moves(&sc.moves)
    );
    match sc.from_usual_text() {
        Some((moves, true)) => {
            let _ = write!(out, "\nscramble from your usual hold: {moves}");
        }
        Some((moves, false)) => {
            let _ = write!(
                out,
                "\nscramble from your usual hold (no turns from it make the pattern exactly as \
                 entered: this makes the same puzzle, which matches it once the whole puzzle is \
                 turned): {moves}"
            );
        }
        None if pi.usual.is_some() => {
            out += "\nno scramble from your usual hold: its colours are not a solved state of \
                    this puzzle";
        }
        None => {}
    }
    text.out(&out);
}
