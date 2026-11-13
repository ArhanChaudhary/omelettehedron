#![deny(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::unreachable,
    clippy::todo,
    clippy::unimplemented
)]
#![cfg_attr(
    test,
    allow(
        clippy::unwrap_used,
        clippy::expect_used,
        clippy::panic,
        clippy::unreachable,
        clippy::todo,
        clippy::unimplemented
    )
)]
#![warn(clippy::pedantic)]
// Pedantic lints that would only add noise here:
#![allow(
    // statistics are reported as f64
    clippy::cast_precision_loss
)]

mod cli;

use clap::{Parser, Subcommand};
use cli::Error;
use std::ffi::OsString;
use std::process::ExitCode;

#[derive(Parser)]
#[command(
    name = "omelettehedron",
    version,
    about = "Two-phase and optimal solvers for the face-turning octahedron (FTO)"
)]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Solve positions: with the two-phase solver, the optimal solver, or both at once
    Solve {
        #[command(subcommand)]
        method: cli::solve::Method,
    },
    /// Generate (or check) the pruning tables of a solver
    Gen {
        #[command(subcommand)]
        solver: cli::gen::Solver,
    },
    /// Benchmark a solver on seeded positions
    Bench {
        #[command(subcommand)]
        solver: cli::bench::Solver,
    },
    #[command(long_about = cli::show::LONG_ABOUT)]
    /// Describe positions (scrambles or patterns) without loading any table
    Show(cli::show::Args),
    /// Positions fixed by symmetries of the octahedron
    ///
    /// Without options: the 33 kinds of subgroup of the 48 symmetries (up to conjugation),
    /// largest first, with how many positions each fixes. With --subgroup or --order: the
    /// positions whose symmetries are exactly that subgroup, one of each set of conjugate
    /// positions (conjugates are the same distance from solved), as patterns for `solve
    /// two-phase --patterns` and `solve optimal --patterns`.
    Symmetric(cli::symmetric::Args),
}

fn main() -> ExitCode {
    match run(Cli::parse_from(progress_values(std::env::args_os()))) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

/// `--progress 5` as `--progress=5`. The interval of `--progress` is optional, so clap only
/// takes it after `=` (otherwise `--progress "U R"` would read the scramble as the interval);
/// a number is never a scramble, so a number after `--progress` is its interval.
fn progress_values(args: impl Iterator<Item = OsString>) -> Vec<OsString> {
    let mut out: Vec<OsString> = vec![];
    for a in args {
        let number = a.to_str().is_some_and(|t| t.parse::<f64>().is_ok());
        match out.last_mut() {
            Some(last) if number && last.as_os_str() == "--progress" => {
                last.push("=");
                last.push(&a);
            }
            _ => out.push(a),
        }
    }
    out
}

fn run(cli: Cli) -> Result<(), Error> {
    match cli.cmd {
        Cmd::Solve { method } => cli::solve::run(method),
        Cmd::Gen { solver } => cli::gen::run(solver),
        Cmd::Bench { solver } => cli::bench::run(solver),
        Cmd::Show(args) => cli::show::run(&args),
        Cmd::Symmetric(args) => cli::symmetric::run(&args),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn progress_interval_after_a_space() {
        let args = |v: &[&str]| -> Vec<OsString> { v.iter().map(OsString::from).collect() };
        let fixed = progress_values(args(&["x", "--progress", "1", "U R"]).into_iter());
        assert_eq!(fixed, args(&["x", "--progress=1", "U R"]));
        // a scramble after --progress stays a scramble, and --progress=5 stays as it is
        for v in [&["x", "--progress", "U"][..], &["x", "--progress=5", "1"]] {
            assert_eq!(progress_values(args(v).into_iter()), args(v));
        }
        assert!(
            Cli::try_parse_from(args(&["x", "solve", "two-phase", "--progress=2", "U"])).is_ok()
        );
    }
}
