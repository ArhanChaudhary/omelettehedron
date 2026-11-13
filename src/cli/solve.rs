//! `omelettehedron solve two-phase | optimal | combined`: solve positions with the two-phase
//! solver, the optimal solver, or both at once. The three take the same options where they
//! mean the same thing, and print the same way (see [`super::report`]).

use super::args::{InputArgs, LimitArgs, OutputArgs, TableArgs, UsualArgs};
use super::input::{none_read, read_jobs, unread, Entry, Job};
use super::load::{
    optimal_tables, optimal_tier, pool, table_error, two_phase_tier, with_bounds, with_two_phase,
};
use super::report::{print_result, Detail, Summary, Text, Tracker};
use super::Error;
use omelettehedron::combined::{self, Source, TwoPhaseMode};
use omelettehedron::optimal::search::{self as opt, Stop};
use omelettehedron::optimal::tables::{Bounds, OptTier};
use omelettehedron::progress::{Bound, Progress, Watch};
use omelettehedron::search::{self as tp, Report};
use omelettehedron::tables::{table_dir, AeIndex, Tables, Tier};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Mutex, PoisonError};
use std::time::Duration;

/// The solver to solve with.
#[derive(clap::Subcommand)]
pub enum Method {
    /// The two-phase solver: a solution within milliseconds, then shorter ones until the best
    /// is proven optimal
    ///
    /// Each better solution is printed as soon as it is found. Without --time, --max-length or
    /// --max-solutions the search runs until its best solution is proven optimal. With --optimal
    /// it goes straight for an optimal solution, which is faster, but shows no solution until
    /// it reaches the optimal length.
    TwoPhase(TwoPhaseArgs),
    /// The optimal solver: its first solution is optimal
    ///
    /// Single-phase IDA* straight to the solved state, with its own tables (`opt-*.prun`). The
    /// first solution it finds is optimal; until then it proves larger and larger lower
    /// bounds. Positions that some of the octahedron's 48 symmetries fix are searched with a
    /// symmetric start. With --max-length it gives up once no solution of at most that many
    /// moves can exist.
    Optimal(OptimalArgs),
    /// Both solvers at once, until their bounds meet (the solution is optimal)
    ///
    /// Each solver runs on half of the threads: the two-phase solver finds short solutions
    /// (upper bounds on the optimal length), the optimal solver proves lower bounds, and the
    /// search stops when they meet or at a limit.
    Combined(CombinedArgs),
}

/// `solve two-phase`.
#[derive(clap::Args)]
pub struct TwoPhaseArgs {
    #[command(flatten)]
    input: InputArgs,
    #[command(flatten)]
    usual: UsualArgs,
    #[command(flatten)]
    limits: LimitArgs,
    /// Stop after this many solutions (every phase-2 completion counts, including ones longer
    /// than the best so far)
    #[arg(short = 'n', long, value_name = "N", help_heading = "Limits")]
    max_solutions: Option<u64>,
    /// Use the two-phase solver to find an optimal solution: faster than improving step by
    /// step (one thread: about 2.6x; target-length passes), but no solution is shown until the
    /// optimal length is reached
    #[arg(long, conflicts_with_all = ["max_length", "max_solutions"], help_heading = "Search")]
    optimal: bool,
    /// With --optimal, to resume a long proof: no solution has fewer than this many moves, as
    /// an earlier run of the same position printed ("pass N done: no solution of at most N
    /// moves" with --progress gives N + 1). The search finds a first solution, then starts its
    /// passes here. A wrong value can give a wrong answer
    #[arg(
        long,
        value_name = "MOVES",
        requires = "optimal",
        help_heading = "Search"
    )]
    lower_bound: Option<usize>,
    /// Solve the positions concurrently, one thread each (-j of them), instead of one at a time
    /// on all threads: faster for many short solves (with -t or -l). Each position's output is
    /// printed when it ends, in input order
    #[arg(long, conflicts_with_all = ["progress", "lower_bound"], help_heading = "Search")]
    batch: bool,
    /// Number of worker threads (default: all cores)
    #[arg(short = 'j', long, help_heading = "Search")]
    threads: Option<usize>,
    /// Pruning-table size: larger tiers search fewer nodes. Without it: the largest tier whose
    /// tables are all in the table directory
    #[arg(long, value_enum, help_heading = "Tables")]
    tier: Option<Tier>,
    #[command(flatten)]
    tables: TableArgs,
    #[command(flatten)]
    output: OutputArgs,
}

/// `solve optimal`.
#[derive(clap::Args)]
pub struct OptimalArgs {
    #[command(flatten)]
    input: InputArgs,
    #[command(flatten)]
    usual: UsualArgs,
    #[command(flatten)]
    limits: LimitArgs,
    /// Print up to this many solutions of the optimal length (all of them if fewer)
    #[arg(short = 'n', long, value_name = "N", help_heading = "Limits")]
    max_solutions: Option<u64>,
    /// Skip the iterations below this depth (one two below it orders the tasks): the search at
    /// this depth finds any solution of at most this many moves, shorter ones too, but one it
    /// finds is not proven optimal. With -l of the same value: whether the position needs more
    /// than this many moves (e.g. `--start-depth 23 -l 23`), without the iterations below
    #[arg(long, value_name = "MOVES", help_heading = "Search")]
    start_depth: Option<usize>,
    /// Solve the positions concurrently, one thread each (-j of them), instead of one at a time
    /// on all threads: faster for many short solves (with -t or -l). Each position's output is
    /// printed when it ends, in input order
    #[arg(long, conflicts_with = "progress", help_heading = "Search")]
    batch: bool,
    /// Number of worker threads (default: all cores)
    #[arg(short = 'j', long, help_heading = "Search")]
    threads: Option<usize>,
    /// Table size: larger tiers search fewer nodes. Without it: the largest tier whose tables
    /// are all in the table directory
    #[arg(long, value_enum, help_heading = "Tables")]
    tier: Option<OptTier>,
    #[command(flatten)]
    tables: TableArgs,
    #[command(flatten)]
    output: OutputArgs,
}

/// `solve combined`.
#[derive(clap::Args)]
pub struct CombinedArgs {
    #[command(flatten)]
    input: InputArgs,
    #[command(flatten)]
    usual: UsualArgs,
    #[command(flatten)]
    limits: LimitArgs,
    /// The two-phase solver goes straight for an optimal solution (as `solve two-phase
    /// --optimal`: target-length passes, from the largest lower bound either solver has
    /// proven), instead of improving its best solution step by step
    #[arg(long, help_heading = "Search")]
    optimal: bool,
    /// Threads in all (default: all cores), split evenly between the two solvers
    #[arg(short = 'j', long, help_heading = "Search")]
    threads: Option<usize>,
    /// Threads of the two-phase solver (default: half of all)
    #[arg(long, value_name = "N", help_heading = "Search")]
    two_phase_threads: Option<usize>,
    /// Threads of the optimal solver (default: the rest)
    #[arg(long, value_name = "N", help_heading = "Search")]
    optimal_threads: Option<usize>,
    /// The two-phase solver's tables (as for `solve two-phase --tier`). Without it: the largest
    /// tier whose tables are all in the table directory
    #[arg(long, value_enum, help_heading = "Tables")]
    tier: Option<Tier>,
    /// The optimal solver's tables (as for `solve optimal --tier`). Without it: the largest
    /// tier whose tables are all in the table directory
    #[arg(long, value_enum, value_name = "TIER", help_heading = "Tables")]
    optimal_tier: Option<OptTier>,
    #[command(flatten)]
    tables: TableArgs,
    #[command(flatten)]
    output: OutputArgs,
}

/// Run `solve`.
pub fn run(method: Method) -> Result<(), Error> {
    match method {
        Method::TwoPhase(a) => two_phase(&a),
        Method::Optimal(a) => optimal(&a),
        Method::Combined(a) => both(&a),
    }
}

/// How much to print for the positions `jobs`.
fn detail(out: &OutputArgs, input: &InputArgs, jobs: &[Entry]) -> Detail {
    Detail {
        quiet: out.quiet,
        verbose: out.verbose,
        stats: out.verbose || out.progress.is_some(),
        header: !out.quiet && (jobs.len() > 1 || out.verbose || input.patterns),
    }
}

/// The positions to solve. If none can be read, each is reported and the run fails before any
/// table is loaded.
fn positions(input: &InputArgs, usual: &UsualArgs, quiet: bool) -> Result<Vec<Entry>, Error> {
    let jobs = read_jobs(input, usual.usual.as_deref())?;
    if !jobs.iter().any(Result::is_ok) {
        none_read(&jobs, quiet)?;
    }
    Ok(jobs)
}

/// The texts (stdout, stderr) of the positions of `--batch` solved and not yet printed, and the
/// next to print.
type Kept = (Vec<Option<(String, String)>>, usize);

/// Solve every position with `solve` (the text of one position, with a pool of threads
/// installed): one at a time on `pool`, or with `batch` one per thread of `pool`, printed in
/// input order. Positions that could not be read are reported in their turn.
fn each(
    jobs: &[Entry],
    pool: &rayon::ThreadPool,
    batch: bool,
    quiet: bool,
    solve: &(dyn Fn(&Job, &Text) -> Result<(), Error> + Sync),
) -> Result<(), Error> {
    if !batch {
        for entry in jobs {
            match entry {
                Ok(job) => pool.install(|| solve(job, &Text::printed()))?,
                Err(bad) => bad.report(quiet),
            }
        }
        return Ok(());
    }
    let next = AtomicUsize::new(0);
    // the texts solved and not yet printed, and the next to print
    let done: Mutex<Kept> = Mutex::new((jobs.iter().map(|_| None).collect(), 0));
    let failed: Mutex<Option<Error>> = Mutex::new(None);
    std::thread::scope(|scope| {
        for _ in 0..pool.current_num_threads().max(1) {
            scope.spawn(|| {
                let Ok(one) = rayon::ThreadPoolBuilder::new().num_threads(1).build() else {
                    return;
                };
                loop {
                    let i = next.fetch_add(1, Ordering::Relaxed);
                    let Some(entry) = jobs.get(i) else { return };
                    let text = Text::kept();
                    match entry {
                        Ok(job) => {
                            // (the first error is the run's, printed at the end)
                            if let Err(e) = one.install(|| solve(job, &text)) {
                                let mut f = failed.lock().unwrap_or_else(PoisonError::into_inner);
                                f.get_or_insert(e);
                            }
                        }
                        Err(bad) => {
                            if !quiet {
                                text.out(&format!("# {}", bad.label));
                            }
                            text.err(&format!("error: {}: {}", bad.label, bad.message));
                        }
                    }
                    let mut guard = done.lock().unwrap_or_else(PoisonError::into_inner);
                    let (texts, first) = &mut *guard;
                    if let Some(slot) = texts.get_mut(i) {
                        *slot = Some(text.take());
                    }
                    // print every finished position whose turn it is
                    while let Some((out, err)) = texts.get_mut(*first).and_then(Option::take) {
                        print!("{out}");
                        eprint!("{err}");
                        *first += 1;
                    }
                }
            });
        }
    });
    match failed.into_inner().unwrap_or_else(PoisonError::into_inner) {
        Some(e) => Err(e),
        None => Ok(()),
    }
}

/// `solve two-phase`.
fn two_phase(a: &TwoPhaseArgs) -> Result<(), Error> {
    let every = a.output.every()?;
    let time = a.limits.duration()?;
    // Improving: the time limit is a stopping rule (the search then keeps reporting
    // improvements). Passes: it is a cancellation deadline, which keeps the passes.
    let limits = tp::Limits {
        max_solutions: a.max_solutions,
        time: if a.optimal { None } else { time },
        cancel_after: if a.optimal { time } else { None },
        max_length: a.limits.max_length,
        lower_bound: a.lower_bound,
        improve: !a.optimal,
        ..tp::Limits::default()
    };
    let quiet = a.output.quiet;
    let jobs = positions(&a.input, &a.usual, quiet)?;
    let d = detail(&a.output, &a.input, &jobs);
    let p = pool(a.threads)?;
    let dir = table_dir(a.tables.tables.as_deref());
    let tier = two_phase_tier(a.tier, &dir, a.tables.move_set(), "--tier", quiet)?;
    p.install(|| {
        with_two_phase!(&a.tables, tier, quiet, |t| {
            let solve = |job: &Job, text: &Text| {
                two_phase_one(&t, job, &limits, every, d, text);
                Ok(())
            };
            each(&jobs, &p, a.batch, quiet, &solve)
        })
    })?;
    unread(&jobs)
}

/// One position with the two-phase solver.
fn two_phase_one<I: AeIndex>(
    t: &Tables<I>,
    job: &Job,
    limits: &tp::Limits,
    every: Option<Duration>,
    d: Detail,
    text: &Text,
) {
    if d.header {
        text.out(&format!("# {}", job.label));
    }
    let tracker = Tracker::new(Source::TwoPhase, text, d.quiet);
    let o = {
        let on = |b: Bound| tracker.bound(b);
        let progress = every.map(|e| Progress::new(e, ""));
        let watch = Watch {
            progress: progress.as_ref(),
            on_bound: Some(&on),
            ..Watch::default()
        };
        tp::solve_watched(t, job.cube, job.frame, limits, Report::Quiet, &watch)
    };
    let r = tracker.finish(vec![]);
    print_result(text, job, &r, d, &[two_phase_stats(&o)]);
}

/// The two-phase solver's statistics (`-v`).
fn two_phase_stats(o: &tp::Outcome) -> String {
    format!(
        "length {} in {:.3}s, {} solutions, {:.3e} nodes{}",
        o.length.map_or("-".into(), |l| l.to_string()),
        o.elapsed,
        o.solutions_found,
        o.nodes as f64,
        if o.optimal { ", optimal" } else { "" }
    )
}

/// The optimal solver's stopping rules.
pub fn optimal_limits(
    limits: &LimitArgs,
    max_solutions: Option<u64>,
    start_depth: Option<usize>,
) -> Result<opt::Limits, Error> {
    if let (Some(d), Some(l)) = (start_depth, limits.max_length) {
        if d > l {
            return Err(format!("--start-depth {d} is above -l {l}").into());
        }
    }
    Ok(opt::Limits {
        max_solutions,
        time: limits.duration()?,
        max_length: limits.max_length,
        start_depth,
    })
}

/// `solve optimal`.
fn optimal(a: &OptimalArgs) -> Result<(), Error> {
    let every = a.output.every()?;
    let limits = optimal_limits(&a.limits, a.max_solutions, a.start_depth)?;
    let quiet = a.output.quiet;
    let jobs = positions(&a.input, &a.usual, quiet)?;
    let d = detail(&a.output, &a.input, &jobs);
    let p = pool(a.threads)?;
    let dir = table_dir(a.tables.tables.as_deref());
    let tier = optimal_tier(a.tier, &dir, a.tables.move_set(), "--tier", quiet)?;
    p.install(|| {
        let t = optimal_tables(&a.tables, tier, quiet)?;
        with_bounds!(&t, |b| {
            let solve = |job: &Job, text: &Text| optimal_one(&b, job, &limits, every, d, text);
            each(&jobs, &p, a.batch, quiet, &solve)
        })
    })?;
    unread(&jobs)
}

/// One position with the optimal solver.
fn optimal_one<B: Bounds>(
    b: &B,
    job: &Job,
    limits: &opt::Limits,
    every: Option<Duration>,
    d: Detail,
    text: &Text,
) -> Result<(), Error> {
    if d.header {
        text.out(&format!("# {}", job.label));
    }
    let tracker = Tracker::new(Source::Optimal, text, d.quiet);
    let o = {
        let on = |b: Bound| tracker.bound(b);
        let progress = every.map(|e| Progress::new(e, ""));
        let watch = Watch {
            progress: progress.as_ref(),
            on_bound: Some(&on),
            ..Watch::default()
        };
        opt::solve_watched(b, job.cube, job.frame, limits, d.verbose, &watch)
    };
    if o.stop == Stop::BadTables {
        return Err(
            "the tables are not distance tables for this position (delete them to \
                    regenerate)"
                .into(),
        );
    }
    let r = tracker.finish(o.solutions.clone());
    print_result(text, job, &r, d, &[optimal_stats(&o)]);
    Ok(())
}

/// The optimal solver's statistics (`-v`).
fn optimal_stats(o: &opt::Outcome) -> String {
    format!(
        "optimal solver: {} in {:.3}s, {} solution(s), {} nodes",
        o.length.map_or_else(
            || format!("no solution, lower bound {}", o.lower_bound),
            |l| {
                if o.stop == Stop::Found {
                    format!("length {l} (not proven optimal)")
                } else {
                    format!("length {l} (optimal)")
                }
            }
        ),
        o.elapsed,
        o.solutions.len(),
        o.nodes
    )
}

/// `solve combined`.
fn both(a: &CombinedArgs) -> Result<(), Error> {
    use omelettehedron::preflight::{check_memory, Planned};
    let every = a.output.every()?;
    let time = a.limits.duration()?;
    let quiet = a.output.quiet;
    let jobs = positions(&a.input, &a.usual, quiet)?;
    let d = detail(&a.output, &a.input, &jobs);
    // threads: split evenly unless given
    let total = a.threads.unwrap_or_else(rayon::current_num_threads).max(2);
    let n_tp = a
        .two_phase_threads
        .unwrap_or_else(|| {
            a.optimal_threads
                .map_or(total.div_ceil(2), |o| total.saturating_sub(o))
        })
        .max(1);
    let n_opt = a
        .optimal_threads
        .unwrap_or_else(|| total.saturating_sub(n_tp))
        .max(1);
    // the tiers (without --tier / --optimal-tier: the largest whose tables are all there)
    let dir = table_dir(a.tables.tables.as_deref());
    let moves = a.tables.move_set();
    let tier = two_phase_tier(a.tier, &dir, moves, "--tier", quiet)?;
    let opt_tier = optimal_tier(a.optimal_tier, &dir, moves, "--optimal-tier", quiet)?;
    // Both solvers' tables are held at once: check that they fit before loading either
    // (2 bits per entry; a missing table is generated, which `load` checks on its own).
    a.tables.apply_checks();
    let files = tier
        .files(moves)
        .into_iter()
        .chain(opt_tier.files(moves).unwrap_or_default());
    let plan: Vec<Planned> = files
        .map(|f| {
            let path = dir.join(&f.name);
            let bytes = std::fs::metadata(&path).map_or(f.entries.div_ceil(4), |m| m.len());
            Planned {
                path,
                bytes,
                working_set: 0,
            }
        })
        .collect();
    check_memory(&plan).map_err(|e| table_error(&dir, &e))?;
    if a.output.verbose {
        eprintln!(
            "two-phase solver: {n_tp} threads, tier {}; optimal solver: {n_opt} threads, tier \
             {}; tables {:.1} GB in all",
            tier.name(),
            opt_tier.name(),
            plan.iter().map(|p| p.bytes).sum::<u64>() as f64 / 1e9
        );
    }
    let (p_tp, p_opt) = (pool(Some(n_tp))?, pool(Some(n_opt))?);
    let cfg = combined::Config {
        time,
        max_length: a.limits.max_length,
        mode: if a.optimal {
            TwoPhaseMode::Passes
        } else {
            TwoPhaseMode::Improve
        },
        progress: every,
        print: !quiet,
    };
    with_two_phase!(&a.tables, tier, quiet, |t| {
        let ot = optimal_tables(&a.tables, opt_tier, quiet)?;
        with_bounds!(&ot, |b| {
            for entry in &jobs {
                let job = match entry {
                    Ok(j) => j,
                    Err(bad) => {
                        bad.report(quiet);
                        continue;
                    }
                };
                let text = Text::printed();
                if d.header {
                    text.out(&format!("# {}", job.label));
                }
                let o = combined::solve(&t, &b, &p_tp, &p_opt, job.cube, job.frame, &cfg)?;
                let stats = [
                    format!("two-phase solver: {}", two_phase_stats(&o.two_phase)),
                    optimal_stats(&o.optimal_solver),
                ];
                let r = Summary {
                    solutions: o.solution.into_iter().collect(),
                    upper: o.upper,
                    lower: o.lower,
                    optimal: o.optimal,
                    elapsed: o.elapsed,
                };
                print_result(&text, job, &r, d, &stats);
            }
            Ok::<(), Error>(())
        })
    })?;
    unread(&jobs)
}
