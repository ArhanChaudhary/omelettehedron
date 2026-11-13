//! `omelettehedron bench two-phase | optimal`: solve seeded positions with a solver and the
//! given limits, and summarise (optionally as JSON).

use super::args::{LimitArgs, TableArgs};
use super::load::{
    optimal_tables, optimal_tier, pool, two_phase_tier, with_bounds, with_two_phase,
};
use super::solve::optimal_limits;
use super::Error;
use omelettehedron::bench;
use omelettehedron::optimal::tables::OptTier;
use omelettehedron::tables::{table_dir, Tier};
use std::path::{Path, PathBuf};

/// The solver to benchmark.
#[derive(clap::Subcommand)]
pub enum Solver {
    /// Benchmark the two-phase solver: for each position the best length over time, the first
    /// solution and whether it was proven optimal; then means, and nodes per second
    ///
    /// Without --time, --max-length or --max-solutions each position is solved to a proof of
    /// optimality with target-length passes (as `solve two-phase --optimal`).
    TwoPhase(TwoPhaseArgs),
    /// Benchmark the optimal solver: each position's optimal length, nodes and time, then a
    /// summary
    Optimal(OptimalArgs),
}

/// The positions of a benchmark.
#[derive(clap::Args)]
#[command(next_help_heading = "Positions")]
pub struct PositionArgs {
    /// Number of positions (default: 20 for two-phase, 10 for optimal)
    #[arg(short = 'N', long, value_name = "N")]
    scrambles: Option<usize>,
    /// Seed of the random positions or scrambles
    #[arg(long, default_value_t = 2026)]
    seed: u64,
    /// Seeded random-move scrambles of this many face turns instead of random positions (e.g.
    /// to benchmark optimal solving: random positions take long to prove)
    #[arg(long, value_name = "MOVES")]
    scramble_len: Option<usize>,
}

/// What a benchmark prints besides its summary.
#[derive(clap::Args)]
#[command(next_help_heading = "Output")]
pub struct BenchOutputArgs {
    /// Write the full results as JSON to this file
    #[arg(long, value_name = "FILE")]
    json: Option<PathBuf>,
    /// Print nothing on stderr but errors (without it: the tables as they are loaded)
    #[arg(short, long, conflicts_with = "verbose")]
    quiet: bool,
    /// Also print each position's result as it ends
    #[arg(short, long)]
    verbose: bool,
}

/// `bench two-phase`.
#[derive(clap::Args)]
pub struct TwoPhaseArgs {
    #[command(flatten)]
    positions: PositionArgs,
    #[command(flatten)]
    limits: LimitArgs,
    /// Stop after this many solutions (every phase-2 completion counts, including ones longer
    /// than the best so far)
    #[arg(short = 'n', long, value_name = "N", help_heading = "Limits")]
    max_solutions: Option<u64>,
    /// Number of worker threads (default: all cores)
    #[arg(short = 'j', long, help_heading = "Search")]
    threads: Option<usize>,
    /// Pruning-table size (as for `solve two-phase --tier`). Without it: the largest tier whose
    /// tables are all in the table directory
    #[arg(long, value_enum, help_heading = "Tables")]
    tier: Option<Tier>,
    #[command(flatten)]
    tables: TableArgs,
    #[command(flatten)]
    output: BenchOutputArgs,
}

/// `bench optimal`.
#[derive(clap::Args)]
pub struct OptimalArgs {
    #[command(flatten)]
    positions: PositionArgs,
    #[command(flatten)]
    limits: LimitArgs,
    /// Find up to this many solutions of the optimal length (all of them if fewer)
    #[arg(short = 'n', long, value_name = "N", help_heading = "Limits")]
    max_solutions: Option<u64>,
    /// Skip the iterations below this depth (as for `solve optimal --start-depth`)
    #[arg(long, value_name = "MOVES", help_heading = "Search")]
    start_depth: Option<usize>,
    /// Number of worker threads (default: all cores)
    #[arg(short = 'j', long, help_heading = "Search")]
    threads: Option<usize>,
    /// Table size (as for `solve optimal --tier`). Without it: the largest tier whose tables
    /// are all in the table directory
    #[arg(long, value_enum, help_heading = "Tables")]
    tier: Option<OptTier>,
    #[command(flatten)]
    tables: TableArgs,
    #[command(flatten)]
    output: BenchOutputArgs,
}

/// Run `bench`.
pub fn run(solver: Solver) -> Result<(), Error> {
    match solver {
        Solver::TwoPhase(a) => two_phase(&a),
        Solver::Optimal(a) => optimal(&a),
    }
}

/// Write `value` as JSON to `path`, if given.
fn write_json(path: Option<&Path>, value: &impl serde::Serialize) -> Result<(), Error> {
    if let Some(path) = path {
        std::fs::write(path, serde_json::to_string_pretty(value)?)
            .map_err(|e| format!("writing {}: {e}", path.display()))?;
    }
    Ok(())
}

/// `bench two-phase`.
fn two_phase(a: &TwoPhaseArgs) -> Result<(), Error> {
    let limits = bench::limits_with_time(a.limits.time, a.limits.max_length, a.max_solutions)
        .map_err(|e| format!("--time: {e}"))?;
    let out = &a.output;
    let threads = pool(a.threads)?;
    let dir = table_dir(a.tables.tables.as_deref());
    let tier = two_phase_tier(a.tier, &dir, a.tables.move_set(), "--tier", out.quiet)?;
    threads.install(|| {
        with_two_phase!(&a.tables, tier, out.quiet, |t| {
            let s = bench::run(
                &t,
                tier.name(),
                a.positions.scrambles.unwrap_or(20),
                a.positions.seed,
                a.positions.scramble_len,
                &limits,
                out.verbose,
            );
            bench::print_summary(&s);
            write_json(out.json.as_deref(), &s)
        })
    })
}

/// `bench optimal`.
fn optimal(a: &OptimalArgs) -> Result<(), Error> {
    let limits = optimal_limits(&a.limits, a.max_solutions, a.start_depth)?;
    let out = &a.output;
    let threads = pool(a.threads)?;
    let dir = table_dir(a.tables.tables.as_deref());
    let tier = optimal_tier(a.tier, &dir, a.tables.move_set(), "--tier", out.quiet)?;
    threads.install(|| {
        let tables = optimal_tables(&a.tables, tier, out.quiet)?;
        with_bounds!(&tables, |b| {
            let s = omelettehedron::optimal::bench::run(
                &b,
                tier.name(),
                a.positions.scrambles.unwrap_or(10),
                a.positions.seed,
                a.positions.scramble_len,
                &limits,
                out.verbose,
            );
            omelettehedron::optimal::bench::print_summary(&s);
            write_json(out.json.as_deref(), &s)
        })
    })
}
