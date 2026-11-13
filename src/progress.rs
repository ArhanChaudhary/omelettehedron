//! What a running search tells the outside world, and what it hears from it.
//!
//! * Progress lines on stderr (`--progress`): each iteration as it ends and, during long
//!   iterations, a heartbeat line with the nodes so far and the tasks started. The search does
//!   no extra work per node for this: the iteration lines are printed between iterations, and
//!   the heartbeat runs on a thread of its own that reads the counters the search already keeps
//!   (the node count is updated once per task).
//! * Links between two searches of the same position (`omelettehedron solve combined`): a stop
//!   flag they share, a callback for every new bound on the optimal length, and a lower bound
//!   proven elsewhere that the two-phase solver's target-length passes start from.
//!
//! [`Watch::default`] changes nothing: the searches run as without it.

use crate::notation::PhysMove;
use std::fmt::Display;
use std::io::Write;
use std::sync::atomic::{AtomicBool, AtomicUsize};
use std::sync::mpsc::{channel, RecvTimeoutError};
use std::time::{Duration, Instant};

/// A bound on the optimal length, found or proven while a search runs.
#[derive(Clone, Copy, Debug)]
pub enum Bound<'a> {
    /// a solution (physical moves): the optimum is at most its length
    Upper(&'a [PhysMove]),
    /// proven: no solution has fewer moves than this
    Lower(usize),
}

/// Progress lines on stderr.
#[derive(Clone, Debug)]
pub struct Progress {
    /// time between two heartbeat lines during an iteration
    pub every: Duration,
    /// start of every line after the time (the solver's name in a combined run)
    pub label: &'static str,
    /// the times printed are measured from here
    pub start: Instant,
}

impl Progress {
    /// Progress lines starting now, with a heartbeat every `every`.
    pub fn new(every: Duration, label: &'static str) -> Progress {
        Progress {
            every,
            label,
            start: Instant::now(),
        }
    }

    /// Print one line: the time since the start, the label and `text`. (One write: stderr is
    /// unbuffered, and `eprintln!` writes each piece of its format separately.)
    pub fn line(&self, text: impl Display) {
        let line = format!(
            "{:>9.3} s  {}{}\n",
            self.start.elapsed().as_secs_f64(),
            self.label,
            text
        );
        let _ = std::io::stderr().lock().write_all(line.as_bytes());
    }
}

/// What a search reports while it runs, and what it hears from outside.
#[derive(Default, Clone, Copy)]
pub struct Watch<'a> {
    /// print progress lines on stderr
    pub progress: Option<&'a Progress>,
    /// stop flag shared with other searches: the search stops soon after it is set (as at a
    /// time limit), and sets it itself when it stops early (a time limit, a proof, a solution
    /// of the optimal solver)
    pub stop: Option<&'a AtomicBool>,
    /// called for every better solution and every larger proven lower bound
    pub on_bound: Option<&'a (dyn Fn(Bound) + Sync)>,
    /// lower bound proven elsewhere (no solution has fewer moves), which may grow while the
    /// search runs. The two-phase solver's optimal mode starts its target-length passes there
    /// and leaves a pass whose target falls below it. Read between tasks only.
    pub lower: Option<&'a AtomicUsize>,
    /// two-phase solver: search as under a time limit (each better solution as soon as
    /// possible, no target-length passes) even without one
    pub improve: bool,
}

impl Watch<'_> {
    /// Report a bound to the callback, if any.
    pub fn bound(&self, b: Bound) {
        if let Some(f) = self.on_bound {
            f(b);
        }
    }

    /// Print a progress line, if progress is on.
    pub fn line(&self, text: impl Display) {
        if let Some(p) = self.progress {
            p.line(text);
        }
    }
}

/// An iteration whose predecessor took less than `1 / SHORT` of the heartbeat interval gets no
/// heartbeat thread: it almost surely ends before the first beat, and starting a thread for
/// each of the many short iterations of a solve cost about 4% of the time of short one-thread
/// solves (12 scrambles of 14-17 moves, 3 rounds: 0.290 -> 0.300-0.305 s).
const SHORT: u32 = 64;

/// Run `f`, an iteration that follows one of `prev` seconds (`None`: the first). With progress
/// on, a thread prints `beat()` every `p.every` while `f` runs (nothing if `f` takes less time;
/// never if `p.every` is zero, or if `prev` is below `p.every / 64`: see [`SHORT`]). If no
/// thread can be started (WebAssembly), `f` runs without it.
pub fn with_heartbeat<R>(
    p: Option<&Progress>,
    prev: Option<f64>,
    beat: impl Fn() -> String + Sync,
    f: impl FnOnce() -> R,
) -> R {
    let short = |p: &Progress| prev.is_some_and(|d| d < (p.every / SHORT).as_secs_f64());
    let Some(p) = p.filter(|p| !p.every.is_zero() && !short(p)) else {
        return f();
    };
    std::thread::scope(|s| {
        let (tx, rx) = channel::<()>();
        let beat = &beat;
        // (the thread ends when `tx` is dropped below; the scope joins it)
        let _ = std::thread::Builder::new()
            .name("progress".into())
            .spawn_scoped(s, move || {
                while let Err(RecvTimeoutError::Timeout) = rx.recv_timeout(p.every) {
                    p.line(beat());
                }
            });
        let r = f();
        drop(tx);
        r
    })
}

/// `n` tasks, e.g. `1 task`, `2688 tasks`.
pub fn tasks(n: usize) -> String {
    format!("{n} task{}", if n == 1 { "" } else { "s" })
}

/// `nodes` in `secs` seconds as a rate, e.g. `2.13e9 nodes/s` (`-` for no time).
pub fn rate(nodes: u64, secs: f64) -> String {
    if secs > 0.0 {
        format!("{:.3e} nodes/s", nodes as f64 / secs)
    } else {
        "- nodes/s".into()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicU64, Ordering};

    #[test]
    fn heartbeat_runs_while_the_work_runs() {
        let p = Progress::new(Duration::from_millis(20), "test: ");
        let beats = AtomicU64::new(0);
        let r = with_heartbeat(
            Some(&p),
            None,
            || {
                beats.fetch_add(1, Ordering::Relaxed);
                "beat".into()
            },
            || {
                std::thread::sleep(Duration::from_millis(150));
                7
            },
        );
        assert_eq!(r, 7);
        assert!(beats.load(Ordering::Relaxed) >= 2);
        // without progress, or after a short iteration, no thread and no beat
        let r = with_heartbeat(None, None, || unreachable!(), || 8);
        assert_eq!(r, 8);
        let r = with_heartbeat(Some(&p), Some(1e-5), || unreachable!(), || 9);
        assert_eq!(r, 9);
        assert_eq!(rate(10, 0.0), "- nodes/s");
        assert_eq!(rate(2_000, 1.0), "2.000e3 nodes/s");
    }
}
