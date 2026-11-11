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
    clippy::cast_precision_loss,
    // the command dispatch is one long match
    clippy::too_many_lines
)]

use clap::{Parser, Subcommand};
use omelettehedron::bench;
use omelettehedron::combined::TwoPhaseMode;
use omelettehedron::cube::{Cube, MoveSet};
use omelettehedron::notation::{
    format_moves, invert_moves, parse_moves, pure_stickers, scramble_to_cube, solves_pure, Frame,
    PhysMove,
};
use omelettehedron::optimal::search::{solve_watched as opt_solve, Limits as OptLimits, Stop};
use omelettehedron::optimal::tables::{OptTables, OptTier};
use omelettehedron::progress::{Progress, Watch};
use omelettehedron::search::{solve_watched, Limits, Report};
use omelettehedron::tables::{table_dir, AeIndex, AeRelabel, AeSym, Tables, Tier};
use omelettehedron::{autotier, info, pattern};
use std::io::BufRead;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

type Error = Box<dyn std::error::Error + Send + Sync>;

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

#[derive(clap::Args, Clone)]
struct TableArgs {
    /// Pruning-table size: small (~135 MB), medium (~1.0 GB), large (~8 GB) or huge (~30 GB).
    /// Larger tiers search fewer nodes; large is the fastest. Without it: the largest tier whose
    /// tables are all in the table directory (`gen` needs it)
    #[arg(long, value_enum)]
    tier: Option<Tier>,
    /// Directory for pruning tables (default: $FTO_TABLES or ./tables)
    #[allow(clippy::doc_markdown)] // shown verbatim by --help
    #[arg(long)]
    tables: Option<PathBuf>,
    /// Number of worker threads (default: all cores)
    #[arg(short = 'j', long)]
    threads: Option<usize>,
    /// Print progress information to stderr
    #[arg(short, long)]
    verbose: bool,
    /// Generate missing tables even if the disk or memory check says they do not fit
    #[arg(long)]
    force: bool,
    /// Move families allowed in solutions besides face turns, comma-separated; every move
    /// counts as one. slice: middle-layer turns (own pruning tables, `*-slice.prun`). wide:
    /// two-layer turns; they never shorten a solution (a wide turn is the opposite face turn
    /// plus a rotation), so only the output changes: wide turns of the back faces instead of
    /// front-face turns, with the front corner held still
    #[arg(long, value_enum, value_delimiter = ',')]
    moves: Vec<Family>,
}

/// A move family for `--moves`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
enum Family {
    /// two-layer turns, e.g. `Uw`
    Wide,
    /// middle-layer turns, e.g. `Us` (= `Ds'`)
    Slice,
}

impl TableArgs {
    fn move_set(&self) -> MoveSet {
        MoveSet {
            wide: self.moves.contains(&Family::Wide),
            slice: self.moves.contains(&Family::Slice),
        }
    }
}

#[derive(clap::Args, Clone)]
struct StopArgs {
    /// Stop after this many solutions (every phase-2 completion counts, including ones
    /// longer than the best so far)
    #[arg(short = 'n', long)]
    max_solutions: Option<u64>,
    /// Stop after this many seconds
    #[arg(short = 't', long)]
    time: Option<f64>,
    /// Stop as soon as a solution of at most this many moves is found
    #[arg(short = 'l', long)]
    max_length: Option<usize>,
}

#[derive(Subcommand)]
enum Cmd {
    /// Solve scrambles given as arguments, or one per line on stdin. Without stopping
    /// options the search runs until the solution is proven optimal.
    Solve {
        /// Scramble in FTO notation (U R F L B BR BL D, ' for anticlockwise; a w after the
        /// face for a wide turn, an s for a slice turn: always accepted)
        scramble: Vec<String>,
        #[command(flatten)]
        stop: StopArgs,
        /// best: print only the final best solution; improvements: print each better
        /// solution as soon as it is found
        #[arg(short, long, value_enum, default_value = "best")]
        print: Report,
        /// Solve this many seeded random positions instead of reading scrambles
        #[arg(long)]
        random: Option<usize>,
        /// Seed of the random positions (`--random`)
        #[arg(long, default_value_t = 1)]
        seed: u64,
        /// Resume a long optimal proof: no solution has fewer than this many moves, as an
        /// earlier run of the same position printed ("pass N done: no solution of at most N
        /// moves" with `--progress` gives N + 1). The search finds a first solution, then starts
        /// its target-length passes here. A wrong value can give a wrong answer
        #[arg(long)]
        lower_bound: Option<usize>,
        /// Each given sequence is a solution: solve the position it solves (the sequence
        /// inverted, applied to the solved puzzle) instead of the position it makes
        #[arg(long)]
        inverse: bool,
        /// Solve a pattern instead of a scramble: the 72 sticker colours as you see them, in
        /// the order of `omelettehedron pattern --help`. Also prints a scramble that makes it
        #[arg(long, group = "patterns_in", conflicts_with_all = ["scramble", "random", "inverse"])]
        pattern: Option<String>,
        /// Solve patterns read from stdin, each as with --pattern. A pattern takes one line or
        /// several (joined until they hold 72 stickers; a blank line ends one); `#` starts a
        /// comment. A pattern with a mistake is reported and the others are still solved
        #[arg(long, group = "patterns_in", conflicts_with_all = ["scramble", "random", "inverse"])]
        patterns: bool,
        /// With --pattern(s): your solved puzzle as you usually hold it to scramble, as the
        /// colours of U R F L B BL D BR (8 characters, the letters of the pattern). The
        /// scramble then starts from that hold
        #[arg(long, requires = "patterns_in")]
        usual: Option<String>,
        #[command(flatten)]
        tab: TableArgs,
        #[command(flatten)]
        prog: ProgressArgs,
    },
    /// Describe positions without loading any table
    ///
    /// For each position: the pattern, whether it is solved, its symmetries, the solved pieces,
    /// whether it is in the two-phase solver's subgroups H_A and H_B, and a lower bound.
    /// Scrambles as arguments or one per line on stdin, as for `solve`.
    #[allow(clippy::doc_markdown)] // shown verbatim by --help
    Info {
        /// Scramble in FTO notation
        scramble: Vec<String>,
        #[command(flatten)]
        input: InputArgs,
        /// Print each position as one line of JSON
        #[arg(long)]
        json: bool,
    },
    /// Print the pattern a scramble makes (the input of `solve --pattern`)
    ///
    /// A pattern is the 72 sticker colours, one character each (any 8 different characters;
    /// spaces, `/`, `,` and `|` are ignored). Hold the puzzle with a corner pointing at you and
    /// read the four front faces U, R, F, L (top, right, bottom, left: clockwise from the top),
    /// then turn the whole puzzle around, keeping the top at the top, and read the four back
    /// faces B, BL, D, BR (again clockwise from the top: B is now at the top and BL on the
    /// right). Read each face with its outer edge at the top: the 5 stickers along that edge
    /// from left to right, the 3 of the middle row, then the one at the corner in the middle of
    /// the picture. Printed patterns use the standard colour scheme held the standard way
    /// (white on U, green on F): W R G P / B O Y A (white red green purple / blue orange yellow
    /// grey) for the solved faces U R F L / B BL D BR. Patterns given to `solve --pattern` may
    /// use any letters, and the puzzle may be held any way and have any colour scheme.
    Pattern {
        /// Scramble in FTO notation
        scramble: Vec<String>,
    },
    /// Generate (or verify) the pruning tables for a tier
    Gen {
        #[command(flatten)]
        tab: TableArgs,
    },
    /// Benchmark: solve seeded random positions with the given stopping rule
    Bench {
        /// number of random positions
        #[arg(short = 'N', long, default_value_t = 20)]
        scrambles: usize,
        /// Seed of the random positions or scrambles
        #[arg(long, default_value_t = 2026)]
        seed: u64,
        /// use seeded random-move scrambles of this many face turns instead of random
        /// positions (e.g. to benchmark optimal solving)
        #[arg(long)]
        scramble_len: Option<usize>,
        #[command(flatten)]
        stop: StopArgs,
        /// write the full results as JSON to this file
        #[arg(long)]
        json: Option<PathBuf>,
        #[command(flatten)]
        tab: TableArgs,
    },
    /// The optimal solver: single-phase IDA* straight to the solved state, with its own tables
    /// (`opt-*.prun`); separate from the two-phase solver of `solve`, `gen` and `bench`
    Optimal {
        #[command(subcommand)]
        cmd: OptCmd,
    },
    /// Solve with both solvers at once: the two-phase solver finds short solutions (upper
    /// bounds on the optimal length), the optimal solver proves lower bounds. Prints each new
    /// bound as it comes and stops when they meet (the solution is optimal) or at the time limit
    Combined(CombinedArgs),
}

/// `--progress`: what a search is doing, on stderr while it runs.
#[derive(clap::Args, Clone)]
struct ProgressArgs {
    /// Print what the search is doing to stderr while it runs: each iteration as it ends
    /// (depth, nodes, nodes/s, best length, lower bound proven) and, during a long iteration, a
    /// line every SECONDS (default 1; 0: none) with the nodes so far and the tasks started
    #[arg(
        long,
        value_name = "SECONDS",
        num_args = 0..=1,
        require_equals = true,
        default_missing_value = "1"
    )]
    progress: Option<f64>,
}

impl ProgressArgs {
    /// The heartbeat interval, if progress is on.
    fn every(&self) -> Result<Option<std::time::Duration>, Error> {
        self.progress
            .map(std::time::Duration::try_from_secs_f64)
            .transpose()
            .map_err(|e| format!("--progress: {e}").into())
    }
}

/// Table options of the optimal solver.
#[derive(clap::Args, Clone)]
struct OptTableArgs {
    /// Table size: small (~195 MB), medium (~3.4 GB, fewer nodes) or large (~22.0 GB, fewest
    /// nodes). Without it: the largest tier whose tables are all in the table directory
    /// (`optimal gen` needs it)
    #[arg(long, value_enum)]
    tier: Option<OptTier>,
    /// Directory for the tables (default: $FTO_TABLES or ./tables)
    #[allow(clippy::doc_markdown)] // shown verbatim by --help
    #[arg(long)]
    tables: Option<PathBuf>,
    /// Number of worker threads (default: all cores)
    #[arg(short = 'j', long)]
    threads: Option<usize>,
    /// Print progress to stderr (the bounds of the start position and every iteration)
    #[arg(short, long)]
    verbose: bool,
    /// Generate missing tables even if the disk or memory check says they do not fit
    #[arg(long)]
    force: bool,
    /// Memory (GB) to leave to the rest of the machine: required free on top of the tables
    /// before generating, and the large tier's generation stops (writing nothing) if the
    /// free memory falls below it
    #[arg(long, value_name = "GB", default_value_t = 0.0)]
    min_free_memory: f64,
    /// Move families allowed in solutions besides face turns, comma-separated, as for the
    /// two-phase solver: slice (middle-layer turns; own tables, `*-slice.prun`), wide (output
    /// only: the same search, wide turns of the back faces printed instead of front-face turns)
    #[arg(long, value_enum, value_delimiter = ',')]
    moves: Vec<Family>,
}

impl OptTableArgs {
    fn move_set(&self) -> MoveSet {
        MoveSet {
            wide: self.moves.contains(&Family::Wide),
            slice: self.moves.contains(&Family::Slice),
        }
    }
}

/// Stopping rules of the optimal solver (without any: the first optimal solution).
#[derive(clap::Args, Clone)]
struct OptStopArgs {
    /// Print up to this many solutions of the optimal length (all of them if fewer)
    #[arg(short = 'n', long)]
    max_solutions: Option<u64>,
    /// Stop after this many seconds and report the lower bound proven so far
    #[arg(short = 't', long)]
    time: Option<f64>,
    /// Give up if the optimum is longer than this many moves
    #[arg(short = 'l', long)]
    max_length: Option<usize>,
}

impl OptStopArgs {
    fn limits(&self) -> Result<OptLimits, Error> {
        Ok(OptLimits {
            max_solutions: self.max_solutions,
            time: self
                .time
                .map(std::time::Duration::try_from_secs_f64)
                .transpose()
                .map_err(|e| format!("--time: {e}"))?,
            max_length: self.max_length,
        })
    }
}

#[derive(Subcommand)]
enum OptCmd {
    /// Solve scrambles given as arguments, or one per line on stdin, optimally
    Solve {
        /// Scramble in FTO notation (as for the two-phase `solve`)
        scramble: Vec<String>,
        #[command(flatten)]
        stop: OptStopArgs,
        /// Solve this many seeded random positions instead of reading scrambles
        #[arg(long)]
        random: Option<usize>,
        /// Seed of the random positions (`--random`)
        #[arg(long, default_value_t = 1)]
        seed: u64,
        /// Each given sequence is a solution: solve the position it solves (the sequence
        /// inverted, applied to the solved puzzle) instead of the position it makes
        #[arg(long)]
        inverse: bool,
        /// Solve a pattern instead of a scramble (as for the two-phase `solve --pattern`); also
        /// prints a scramble that makes it
        #[arg(long, group = "patterns_in", conflicts_with_all = ["scramble", "random", "inverse"])]
        pattern: Option<String>,
        /// Solve patterns read from stdin (as for the two-phase `solve --patterns`)
        #[arg(long, group = "patterns_in", conflicts_with_all = ["scramble", "random", "inverse"])]
        patterns: bool,
        /// With --pattern(s): your solved puzzle as you usually hold it, as the colours of U R F
        /// L B BL D BR (as for the two-phase `solve --usual`)
        #[arg(long, requires = "patterns_in")]
        usual: Option<String>,
        #[command(flatten)]
        tab: OptTableArgs,
        #[command(flatten)]
        prog: ProgressArgs,
    },
    /// Generate (or verify) the tables of a tier
    Gen {
        #[command(flatten)]
        tab: OptTableArgs,
    },
    /// Benchmark: solve seeded positions optimally; prints each position's optimal length,
    /// nodes and time, then a summary
    Bench {
        /// number of positions
        #[arg(short = 'N', long, default_value_t = 10)]
        scrambles: usize,
        /// Seed of the scrambles or random positions
        #[arg(long, default_value_t = 2026)]
        seed: u64,
        /// seeded random-move scrambles of this many face turns (the two-phase `bench
        /// --scramble-len` sets); without it, random positions (slow)
        #[arg(long)]
        scramble_len: Option<usize>,
        #[command(flatten)]
        stop: OptStopArgs,
        /// write the full results as JSON to this file
        #[arg(long)]
        json: Option<PathBuf>,
        #[command(flatten)]
        tab: OptTableArgs,
    },
}

#[cfg(not(target_family = "wasm"))]
fn pool(threads: Option<usize>) -> Result<rayon::ThreadPool, Error> {
    let mut b = rayon::ThreadPoolBuilder::new();
    if let Some(t) = threads {
        b = b.num_threads(t);
    }
    Ok(b.build()?)
}

/// WebAssembly (WASI) has no threads: the only worker is the current thread, for this pool and
/// for the global one (which table generation may use).
#[cfg(target_family = "wasm")]
fn pool(_threads: Option<usize>) -> Result<rayon::ThreadPool, Error> {
    // (fails harmlessly if the global pool already exists)
    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .use_current_thread()
        .build_global();
    Ok(rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .use_current_thread()
        .build()?)
}

fn table_error(dir: &Path, e: &std::io::Error) -> Error {
    format!("pruning tables in {}: {e}", dir.display()).into()
}

fn read_scrambles(args: &[String]) -> Vec<String> {
    if !args.is_empty() {
        return vec![args.join(" ")];
    }
    std::io::stdin()
        .lock()
        .lines()
        .map_while(Result::ok)
        .filter(|l| !l.trim().is_empty())
        .collect()
}

/// A pattern to solve (`solve --pattern`).
struct PatternInput {
    pattern: [char; pattern::STICKERS],
    usual: Option<[char; 8]>,
    cube: Cube,
}

fn read_pattern(text: &str, usual: Option<&str>) -> Result<PatternInput, Error> {
    let p = pattern::parse(text).map_err(|e| format!("pattern: {e}"))?;
    let cube = pattern::position(&p).map_err(|e| format!("pattern: {e}"))?;
    Ok(PatternInput {
        pattern: p,
        usual: parse_usual(usual)?,
        cube,
    })
}

/// The colours of `--usual`.
fn parse_usual(usual: Option<&str>) -> Result<Option<[char; 8]>, Error> {
    let Some(u) = usual else { return Ok(None) };
    let v: Vec<char> = u
        .chars()
        .filter(|c| !c.is_whitespace() && *c != '/')
        .collect();
    let u: [char; 8] = v
        .as_slice()
        .try_into()
        .map_err(|_| "--usual: give the 8 colours of U R F L B BL D BR")?;
    Ok(Some(u))
}

/// Print the scramble that makes a solved pattern (and check the solution).
fn print_pattern_scramble(pi: &PatternInput, solution: &[PhysMove]) {
    let Some(sc) = pattern::scramble(&pi.pattern, solution, pi.usual.as_ref()) else {
        eprintln!("error: the solution does not solve the pattern (please report this pattern)");
        return;
    };
    println!(
        "scramble from solved with {} on U and {} on F: {}",
        sc.start[0],
        sc.start[2],
        format_moves(&sc.moves)
    );
    match sc.from_usual_text() {
        Some((text, true)) => println!("scramble from your usual hold: {text}"),
        Some((text, false)) => println!(
            "scramble from your usual hold (no turns from it make the pattern exactly as \
             entered: this makes the same puzzle, which matches it once the whole puzzle is \
             turned): {text}"
        ),
        None if pi.usual.is_some() => println!(
            "no scramble from your usual hold: its colours are not a solved state of this \
             puzzle"
        ),
        None => {}
    }
}

fn run_solve<I: AeIndex>(
    t: &Tables<I>,
    scrambles: &[(Cube, Frame, String)],
    limits: &Limits,
    print: Report,
    verbose: bool,
    input: Option<&PatternInput>,
    progress: Option<std::time::Duration>,
) {
    for (cube, frame, label) in scrambles {
        if scrambles.len() > 1 || verbose {
            println!("# {label}");
        }
        let prog = progress.map(|every| Progress::new(every, ""));
        let watch = Watch {
            progress: prog.as_ref(),
            ..Watch::default()
        };
        let o = solve_watched(t, *cube, *frame, limits, print, &watch);
        if let (Some(pi), Some(sol)) = (input, &o.solution) {
            print_pattern_scramble(pi, sol);
        }
        if verbose || print == Report::Improvements || progress.is_some() {
            eprintln!(
                "length {} in {:.3}s, {} solutions, {:.3e} nodes{}",
                o.length.map_or("-".into(), |l| l.to_string()),
                o.elapsed,
                o.solutions_found,
                o.nodes as f64,
                if o.optimal { ", optimal" } else { "" }
            );
        }
        if o.solution.is_none() {
            println!("no solution found");
        }
        if verbose {
            // the solver works up to centre colours; say whether the solution also restores
            // every centre triangle to its own place
            if let (Ok(scr), Some(sol)) = (parse_moves(label), &o.solution) {
                eprintln!(
                    "pure FTO (every centre unique): {}",
                    if solves_pure(&scr, sol) {
                        "also solved"
                    } else {
                        "not solved (centres of a colour exchanged)"
                    }
                );
            }
        }
    }
}

/// Load the tables of the tier and evaluate `$body` with them; table errors are returned
/// from the enclosing function.
macro_rules! dispatch {
    ($tab:expr, |$t:ident, $name:ident| $body:expr) => {{
        let dir = table_dir($tab.tables.as_deref());
        let moves = $tab.move_set();
        let tier = two_phase_tier(&$tab, &dir)?;
        omelettehedron::preflight::set_force($tab.force);
        // The large tier needs 8 GB and the huge tier's index has 9.1e10 entries: 32-bit
        // targets such as WebAssembly (4 GB) support the small and medium tiers (the medium
        // index has 4.13e9 entries, below 2^32).
        if cfg!(target_pointer_width = "32") && matches!(tier, Tier::Large | Tier::Huge) {
            return Err(
                "only --tier small and medium are supported on 32-bit targets such as \
                        WebAssembly"
                    .into(),
            );
        }
        match tier {
            Tier::Small => {
                let $t = Tables::load(&dir, AeSym::new(false), moves, $tab.verbose)
                    .map_err(|e| table_error(&dir, &e))?;
                let $name = "small";
                $body
            }
            Tier::Medium => {
                let $t = Tables::load(&dir, AeSym::new(true), moves, $tab.verbose)
                    .map_err(|e| table_error(&dir, &e))?;
                let $name = "medium";
                $body
            }
            Tier::Large => {
                let $t = Tables::load_with(&dir, AeSym::new(true), true, moves, $tab.verbose)
                    .map_err(|e| table_error(&dir, &e))?;
                let $name = "large";
                $body
            }
            Tier::Huge => {
                let $t = Tables::load_with(&dir, AeRelabel::new(), true, moves, $tab.verbose)
                    .map_err(|e| table_error(&dir, &e))?;
                let $name = "huge";
                $body
            }
        }
    }};
}

/// Load (generating if needed) the optimal solver's tables.
fn load_opt(tab: &OptTableArgs) -> Result<OptTables, Error> {
    let dir = table_dir(tab.tables.as_deref());
    omelettehedron::preflight::set_force(tab.force);
    // (an `as` cast saturates: negative or NaN -> 0, too large -> u64::MAX)
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let floor = (tab.min_free_memory * 1e9) as u64;
    omelettehedron::preflight::set_memory_floor(floor);
    let tier = if let Some(t) = tab.tier {
        t
    } else {
        let c = autotier::optimal(&dir, tab.move_set())?;
        eprintln!("{}", c.message);
        c.tier
    };
    OptTables::load(&dir, tier, tab.move_set(), tab.verbose).map_err(|e| table_error(&dir, &e))
}

/// Evaluate `$body` with `$b` bound to the lower bounds of the tables' tier.
macro_rules! opt_bounds {
    ($t:expr, |$b:ident| $body:expr) => {{
        let missing = || -> Error { "the loaded tables are not the tier's".into() };
        match $t.tier {
            OptTier::Small => {
                let $b = $t.small().ok_or_else(missing)?;
                $body
            }
            OptTier::Medium => {
                let $b = $t.medium().ok_or_else(missing)?;
                $body
            }
            OptTier::Large => {
                let $b = $t.large().ok_or_else(missing)?;
                $body
            }
        }
    }};
}

fn run_optimal(cmd: OptCmd) -> Result<(), Error> {
    match cmd {
        OptCmd::Solve {
            scramble,
            stop,
            random,
            seed,
            inverse,
            pattern,
            patterns,
            usual,
            tab,
            prog,
        } => {
            let limits = stop.limits()?;
            let every = prog.every()?;
            let jobs = read_jobs(&Inputs {
                scramble: &scramble,
                pattern: pattern.as_deref(),
                patterns,
                usual: usual.as_deref(),
                random,
                seed,
                inverse,
            })?;
            if !jobs.iter().any(Result::is_ok) {
                return none_read(&jobs);
            }
            let p = pool(tab.threads)?;
            p.install(|| {
                let t = load_opt(&tab)?;
                opt_bounds!(t, |b| {
                    for job in &jobs {
                        let Job {
                            cube, frame, label, ..
                        } = match job {
                            Ok(j) => j,
                            Err(bad) => {
                                bad.report();
                                continue;
                            }
                        };
                        if jobs.len() > 1 || tab.verbose {
                            println!("# {label}");
                        }
                        let progress = every.map(|e| Progress::new(e, ""));
                        let watch = Watch {
                            progress: progress.as_ref(),
                            ..Watch::default()
                        };
                        let o = opt_solve(&b, *cube, *frame, &limits, tab.verbose, &watch);
                        for s in &o.solutions {
                            println!("{}", format_moves(s));
                        }
                        if let (
                            Ok(Job {
                                pattern: Some(pi), ..
                            }),
                            Some(sol),
                        ) = (job, o.solutions.first())
                        {
                            print_pattern_scramble(pi, sol);
                        }
                        match o.stop {
                            Stop::Solved => {}
                            Stop::Time => println!(
                                "no solution found (time limit): none has fewer than {} moves",
                                o.lower_bound
                            ),
                            Stop::MaxLength => println!(
                                "no solution of at most {} moves",
                                o.lower_bound.saturating_sub(1)
                            ),
                            Stop::BadTables => {
                                return Err("the tables are not distance tables for this \
                                            position (delete them to regenerate)"
                                    .into());
                            }
                        }
                        if tab.verbose || every.is_some() {
                            eprintln!(
                                "optimal solver: {} in {:.3}s, {} solution(s), {} nodes",
                                o.length.map_or_else(
                                    || format!("no solution, lower bound {}", o.lower_bound),
                                    |l| format!("length {l} (optimal)")
                                ),
                                o.elapsed,
                                o.solutions.len(),
                                o.nodes
                            );
                        }
                        if tab.verbose {
                            if let (Ok(scr), Some(sol)) = (parse_moves(label), o.solutions.first())
                            {
                                eprintln!(
                                    "pure FTO (every centre unique): {}",
                                    if solves_pure(&scr, sol) {
                                        "also solved"
                                    } else {
                                        "not solved (centres of a colour exchanged)"
                                    }
                                );
                            }
                        }
                    }
                    Ok::<(), Error>(())
                })
            })?;
            unread(&jobs)
        }
        OptCmd::Gen { tab } => {
            if tab.tier.is_none() {
                return Err(
                    "optimal gen needs --tier: small, medium or large (the tables to \
                            generate)"
                        .into(),
                );
            }
            let p = pool(tab.threads)?;
            let mut tab = tab;
            tab.verbose = true;
            p.install(|| {
                let t = load_opt(&tab)?;
                eprintln!(
                    "optimal solver, tier {}, moves {}: {:.1} MB of tables",
                    t.tier.name(),
                    t.moves.name(),
                    t.bytes() as f64 / 1e6
                );
                eprintln!("corners x centres depth histogram {:?}", t.cc.hist);
                if let Some(e) = &t.edges {
                    eprintln!("edge permutation depth histogram {:?}", e.hist);
                }
                if let Some((_, ce)) = &t.ce {
                    eprintln!("corners x edges depth histogram {:?}", ce.hist);
                }
                if let Some((_, cbe, _)) = &t.cbe {
                    eprintln!(
                        "corners x B-centres x edge triple depth histogram {:?}",
                        cbe.hist
                    );
                }
                eprintln!(
                    "A-centre/edge (to H_A, two-phase table) depth histogram {:?}",
                    t.ae.tab.hist
                );
                Ok(())
            })
        }
        OptCmd::Bench {
            scrambles,
            seed,
            scramble_len,
            stop,
            json,
            tab,
        } => {
            let limits = stop.limits()?;
            let p = pool(tab.threads)?;
            p.install(|| {
                let t = load_opt(&tab)?;
                opt_bounds!(t, |b| {
                    let s = omelettehedron::optimal::bench::run(
                        &b,
                        t.tier.name(),
                        scrambles,
                        seed,
                        scramble_len,
                        &limits,
                        tab.verbose,
                    );
                    omelettehedron::optimal::bench::print_summary(&s);
                    if let Some(path) = &json {
                        std::fs::write(path, serde_json::to_string_pretty(&s)?)
                            .map_err(|e| format!("writing {}: {e}", path.display()))?;
                    }
                    Ok(())
                })
            })
        }
    }
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

fn run(cli: Cli) -> Result<(), Error> {
    match cli.cmd {
        Cmd::Solve {
            scramble,
            stop,
            print,
            random,
            seed,
            lower_bound,
            inverse,
            pattern,
            patterns,
            usual,
            tab,
            prog,
        } => {
            let progress = prog.every()?;
            let mut limits =
                bench::limits_with_time(stop.time, stop.max_length, stop.max_solutions)
                    .map_err(|e| format!("--time: {e}"))?;
            limits.lower_bound = lower_bound;
            if patterns {
                return solve_pattern_list(&tab, usual.as_deref(), &limits, print, progress);
            }
            let input = match &pattern {
                Some(text) => Some(read_pattern(text, usual.as_deref())?),
                None => None,
            };
            let positions: Vec<(Cube, Frame, String)> = if let Some(pi) = &input {
                vec![(pi.cube, Frame::identity(), pattern::format(&pi.pattern))]
            } else if let Some(n) = random {
                bench::random_positions(n, seed)
                    .into_iter()
                    .enumerate()
                    .map(|(i, c)| {
                        (
                            c,
                            Frame::identity(),
                            format!("random position {i} (seed {seed})"),
                        )
                    })
                    .collect()
            } else {
                let mut v = vec![];
                for s in read_scrambles(&scramble) {
                    match parse_moves(&s) {
                        Ok(m) if inverse => {
                            let m = invert_moves(&m);
                            let (c, f) = scramble_to_cube(&m);
                            v.push((c, f, format_moves(&m)));
                        }
                        Ok(m) => {
                            let (c, f) = scramble_to_cube(&m);
                            v.push((c, f, s));
                        }
                        Err(e) => {
                            eprintln!("error: {e}");
                            std::process::exit(2);
                        }
                    }
                }
                v
            };
            if positions.is_empty() {
                eprintln!("error: no scramble given");
                std::process::exit(2);
            }
            let p = pool(tab.threads)?;
            p.install(|| {
                dispatch!(tab, |t, _name| run_solve(
                    &t,
                    &positions,
                    &limits,
                    print,
                    tab.verbose,
                    input.as_ref(),
                    progress
                ));
                Ok(())
            })
        }
        Cmd::Pattern { scramble } => {
            for s in read_scrambles(&scramble) {
                let m = parse_moves(&s).map_err(|e| e.clone())?;
                println!("{}", pattern::format(&pattern::of_moves(&m)));
            }
            Ok(())
        }
        Cmd::Gen { tab } => {
            if tab.tier.is_none() {
                return Err(
                    "gen needs --tier: small, medium, large or huge (the tables to \
                            generate)"
                        .into(),
                );
            }
            let p = pool(tab.threads)?;
            let mut tab = tab;
            tab.verbose = true;
            p.install(|| {
                dispatch!(tab, |t, name| {
                    eprintln!(
                        "tier {name}, moves {}: AE {:.1} MB, CB {:.1} MB; phase-2 depths {:?}",
                        t.moves.name(),
                        t.ae.bytes() as f64 / 1e6,
                        t.cb.tab.bytes() as f64 / 1e6,
                        t.p2.histogram()
                    );
                    eprintln!("AE depth histogram {:?}", t.ae.hist);
                    eprintln!("CB depth histogram {:?}", t.cb.tab.hist);
                    if let Some(cbe) = &t.cbe {
                        eprintln!("CBE depth histogram {:?}", cbe.hist);
                    }
                    // (read it now: `gen` checks every table of the tier)
                    if let (Some(lazy), Some(cbe)) = (&t.cbs, &t.cbe) {
                        let Some(cbs) = lazy.get(&cbe.ix) else {
                            return Err(format!(
                                "{}: unusable distance-to-solved table",
                                lazy.path().display()
                            )
                            .into());
                        };
                        eprintln!("distance-to-solved depth histogram {:?}", cbs.hist);
                    }
                });
                Ok(())
            })
        }
        Cmd::Bench {
            scrambles,
            seed,
            scramble_len,
            stop,
            json,
            tab,
        } => {
            let limits = bench::limits_with_time(stop.time, stop.max_length, stop.max_solutions)
                .map_err(|e| format!("--time: {e}"))?;
            let p = pool(tab.threads)?;
            p.install(|| {
                dispatch!(tab, |t, name| {
                    let s = bench::run(
                        &t,
                        name,
                        scrambles,
                        seed,
                        scramble_len,
                        &limits,
                        tab.verbose,
                    );
                    bench::print_summary(&s);
                    if let Some(path) = &json {
                        std::fs::write(path, serde_json::to_string_pretty(&s)?)
                            .map_err(|e| format!("writing {}: {e}", path.display()))?;
                    }
                });
                Ok(())
            })
        }
        Cmd::Optimal { cmd } => run_optimal(cmd),
        Cmd::Combined(args) => run_combined(&args),
        Cmd::Info {
            scramble,
            input,
            json,
        } => {
            let jobs = read_jobs(&Inputs {
                scramble: &scramble,
                pattern: input.pattern.as_deref(),
                patterns: input.patterns,
                usual: None,
                random: input.random,
                seed: input.seed,
                inverse: input.inverse,
            })?;
            for job in &jobs {
                match job {
                    Ok(j) => print_info(j, jobs.len() > 1, json)?,
                    // (one line of JSON for it too, in its place)
                    Err(bad) if json => {
                        let out = serde_json::json!({ "input": bad.label, "error": bad.message });
                        println!("{out}");
                        eprintln!("error: {}: {}", bad.label, bad.message);
                    }
                    Err(bad) => bad.report(),
                }
            }
            unread(&jobs)
        }
    }
}

/// The two-phase tier: `--tier`, or the largest whose tables are all in `dir` (said on stderr).
fn two_phase_tier(tab: &TableArgs, dir: &Path) -> Result<Tier, Error> {
    if let Some(t) = tab.tier {
        return Ok(t);
    }
    let c = autotier::two_phase(dir, tab.move_set())?;
    eprintln!("{}", c.message);
    Ok(c.tier)
}

/// The inputs of `info` besides the scrambles.
#[derive(clap::Args, Clone)]
struct InputArgs {
    /// Each given sequence is a solution: the position it solves (the sequence inverted)
    #[arg(long)]
    inverse: bool,
    /// A pattern instead of a scramble (as for `solve --pattern`)
    #[arg(long, conflicts_with_all = ["scramble", "inverse", "random"])]
    pattern: Option<String>,
    /// Patterns read from stdin (as for `solve --patterns`)
    #[arg(long, conflicts_with_all = ["scramble", "inverse", "random", "pattern"])]
    patterns: bool,
    /// This many seeded random positions instead of scrambles
    #[arg(long, conflicts_with_all = ["scramble", "inverse"])]
    random: Option<usize>,
    /// Seed of the random positions (`--random`)
    #[arg(long, default_value_t = 1)]
    seed: u64,
}

/// A position to work on (`optimal solve`, `info`, `solve --patterns`).
struct Job {
    /// printed as `# label` when there are several
    label: String,
    /// a shorter name, if the label holds the pattern (`info` prints the pattern anyway)
    name: Option<String>,
    cube: Cube,
    frame: Frame,
    /// the moves that make it, if given as a scramble
    moves: Option<Vec<PhysMove>>,
    /// the pattern, if given as one: its solution is followed by a scramble that makes it
    pattern: Option<PatternInput>,
}

/// A pattern of `--patterns` with a mistake: reported in its turn, and the run fails at the end.
struct BadInput {
    label: String,
    message: String,
}

impl BadInput {
    fn report(&self) {
        println!("# {}", self.label);
        eprintln!("error: {}: {}", self.label, self.message);
    }
}

/// Where the positions of a command come from.
struct Inputs<'a> {
    scramble: &'a [String],
    pattern: Option<&'a str>,
    patterns: bool,
    usual: Option<&'a str>,
    random: Option<usize>,
    seed: u64,
    inverse: bool,
}

/// The positions to work on, in order: `--patterns` (stdin; a pattern with a mistake is kept
/// as its error), `--pattern`, `--random`, or the scrambles (arguments, or one per line on
/// stdin), inverted with `--inverse`.
fn read_jobs(i: &Inputs) -> Result<Vec<Result<Job, BadInput>>, Error> {
    if i.patterns {
        let usual = parse_usual(i.usual)?;
        let mut text = String::new();
        std::io::Read::read_to_string(&mut std::io::stdin().lock(), &mut text)?;
        let records = pattern::records(&text);
        if records.is_empty() {
            return Err("no pattern on stdin".into());
        }
        return Ok(records
            .iter()
            .enumerate()
            .map(|(k, r)| {
                let name = format!("pattern {} (line {})", k + 1, r.line);
                let read = pattern::parse(&r.text).and_then(|p| Ok((p, pattern::position(&p)?)));
                match read {
                    Ok((p, cube)) => Ok(Job {
                        label: format!("{name}: {}", pattern::format(&p)),
                        name: Some(name),
                        cube,
                        frame: Frame::identity(),
                        moves: None,
                        pattern: Some(PatternInput {
                            pattern: p,
                            usual,
                            cube,
                        }),
                    }),
                    Err(e) => Err(BadInput {
                        label: name,
                        message: e.to_string(),
                    }),
                }
            })
            .collect());
    }
    if let Some(text) = i.pattern {
        let pi = read_pattern(text, i.usual)?;
        return Ok(vec![Ok(Job {
            label: pattern::format(&pi.pattern),
            name: None,
            cube: pi.cube,
            frame: Frame::identity(),
            moves: None,
            pattern: Some(pi),
        })]);
    }
    if let Some(n) = i.random {
        return Ok(bench::random_positions(n, i.seed)
            .into_iter()
            .enumerate()
            .map(|(k, cube)| {
                Ok(Job {
                    label: format!("random position {k} (seed {})", i.seed),
                    name: None,
                    cube,
                    frame: Frame::identity(),
                    moves: None,
                    pattern: None,
                })
            })
            .collect());
    }
    let mut v = vec![];
    for s in read_scrambles(i.scramble) {
        let m = parse_moves(&s)?;
        let (m, label) = if i.inverse {
            let m = invert_moves(&m);
            let label = format_moves(&m);
            (m, label)
        } else {
            (m, s)
        };
        let (cube, frame) = scramble_to_cube(&m);
        v.push(Ok(Job {
            label,
            name: None,
            cube,
            frame,
            moves: Some(m),
            pattern: None,
        }));
    }
    if v.is_empty() {
        return Err("no scramble given".into());
    }
    Ok(v)
}

/// When no pattern could be read: report each, and fail (no table is loaded).
fn none_read(jobs: &[Result<Job, BadInput>]) -> Result<(), Error> {
    for bad in jobs.iter().filter_map(|j| j.as_ref().err()) {
        bad.report();
    }
    unread(jobs)
}

/// After the jobs: an error if some pattern could not be read (each was reported in its turn).
fn unread(jobs: &[Result<Job, BadInput>]) -> Result<(), Error> {
    match jobs.iter().filter(|j| j.is_err()).count() {
        0 => Ok(()),
        bad => Err(format!(
            "{bad} of {} patterns could not be read (see above)",
            jobs.len()
        )
        .into()),
    }
}

/// `solve --patterns`: each pattern of stdin as with `--pattern`; one with a mistake is
/// reported in its turn, the others are solved, and the run fails at the end.
fn solve_pattern_list(
    tab: &TableArgs,
    usual: Option<&str>,
    limits: &Limits,
    print: Report,
    progress: Option<std::time::Duration>,
) -> Result<(), Error> {
    let jobs = read_jobs(&Inputs {
        scramble: &[],
        pattern: None,
        patterns: true,
        usual,
        random: None,
        seed: 0,
        inverse: false,
    })?;
    if !jobs.iter().any(Result::is_ok) {
        return none_read(&jobs);
    }
    let p = pool(tab.threads)?;
    p.install(|| {
        dispatch!(tab, |t, _name| {
            for job in &jobs {
                match job {
                    Ok(j) => {
                        // (`run_solve` prints the label itself with -v)
                        if !tab.verbose {
                            println!("# {}", j.label);
                        }
                        run_solve(
                            &t,
                            &[(j.cube, j.frame, j.label.clone())],
                            limits,
                            print,
                            tab.verbose,
                            j.pattern.as_ref(),
                            progress,
                        );
                    }
                    Err(bad) => bad.report(),
                }
            }
        });
        Ok::<(), Error>(())
    })?;
    unread(&jobs)
}

/// `info`: the facts about one position, as text or one line of JSON.
fn print_info(j: &Job, header: bool, json: bool) -> Result<(), Error> {
    let labels = j.moves.as_deref().map(pure_stickers);
    let facts = info::info(&j.cube, labels.as_ref());
    let pat = match (&j.pattern, &j.moves) {
        (Some(pi), _) => pattern::format(&pi.pattern),
        (None, Some(m)) => pattern::format(&pattern::of_moves(m)),
        (None, None) => pattern::format(&pattern::of_cube(&j.cube)),
    };
    if json {
        #[derive(serde::Serialize)]
        struct Out<'a> {
            input: &'a str,
            pattern: &'a str,
            #[serde(flatten)]
            info: &'a info::Info,
        }
        let out = Out {
            input: j.name.as_deref().unwrap_or(&j.label),
            pattern: &pat,
            info: &facts,
        };
        println!("{}", serde_json::to_string(&out)?);
        return Ok(());
    }
    if header {
        println!("# {}", j.name.as_deref().unwrap_or(&j.label));
    }
    println!("pattern: {pat}");
    for line in facts.lines() {
        println!("{line}");
    }
    Ok(())
}

/// `omelettehedron combined`: both solvers on each position at once.
#[derive(clap::Args)]
struct CombinedArgs {
    /// Scramble in FTO notation (as for `solve`); without one, one scramble per line on stdin
    scramble: Vec<String>,
    /// Stop after this many seconds and print the best solution and the bounds proven
    #[arg(short = 't', long)]
    time: Option<f64>,
    /// Tables of the two-phase solver (as for `solve --tier`; default: the largest tier whose
    /// tables are all there)
    #[arg(long, value_enum)]
    tier: Option<Tier>,
    /// Tables of the optimal solver (as for `optimal solve --tier`; default: the largest tier
    /// whose tables are all there)
    #[arg(long, value_enum)]
    opt_tier: Option<OptTier>,
    /// How the two-phase solver searches. improve: as under a time limit, each better solution
    /// as soon as possible (the optimal solver proves the lower bounds). passes: its optimal
    /// mode, target-length passes that start from the largest lower bound proven by either
    #[arg(long, value_enum, default_value = "improve")]
    two_phase: TwoPhaseMode,
    /// Threads in all (default: all cores), split evenly between the two solvers
    #[arg(short = 'j', long)]
    threads: Option<usize>,
    /// Threads of the two-phase solver (default: half of all)
    #[arg(long)]
    two_phase_threads: Option<usize>,
    /// Threads of the optimal solver (default: the rest)
    #[arg(long)]
    optimal_threads: Option<usize>,
    /// Directory for the tables of both solvers (default: $FTO_TABLES or ./tables)
    #[allow(clippy::doc_markdown)] // shown verbatim by --help
    #[arg(long)]
    tables: Option<PathBuf>,
    /// Move families allowed besides face turns, as for `solve --moves`
    #[arg(long, value_enum, value_delimiter = ',')]
    moves: Vec<Family>,
    /// Solve this many seeded random positions instead of reading scrambles
    #[arg(long)]
    random: Option<usize>,
    /// Seed of the random positions (`--random`)
    #[arg(long, default_value_t = 1)]
    seed: u64,
    /// Each given sequence is a solution: solve the position it solves
    #[arg(long)]
    inverse: bool,
    /// Print table loading and the thread split to stderr
    #[arg(short, long)]
    verbose: bool,
    /// Load or generate the tables even if the memory check says they do not fit
    #[arg(long)]
    force: bool,
    /// Memory (GB) to leave to the rest of the machine (as for `optimal gen`)
    #[arg(long, value_name = "GB", default_value_t = 0.0)]
    min_free_memory: f64,
    #[command(flatten)]
    prog: ProgressArgs,
}

fn run_combined(a: &CombinedArgs) -> Result<(), Error> {
    use omelettehedron::combined;
    use omelettehedron::preflight::{check_memory, set_force, set_memory_floor, Planned};
    let every = a.prog.every()?;
    let time = a
        .time
        .map(std::time::Duration::try_from_secs_f64)
        .transpose()
        .map_err(|e| format!("--time: {e}"))?;
    let jobs = read_jobs(&Inputs {
        scramble: &a.scramble,
        pattern: None,
        patterns: false,
        usual: None,
        random: a.random,
        seed: a.seed,
        inverse: a.inverse,
    })?;
    if !jobs.iter().any(Result::is_ok) {
        return none_read(&jobs);
    }
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
    let mut tab = TableArgs {
        tier: a.tier,
        tables: a.tables.clone(),
        threads: Some(n_tp),
        verbose: a.verbose,
        force: a.force,
        moves: a.moves.clone(),
    };
    let mut otab = OptTableArgs {
        tier: a.opt_tier,
        tables: a.tables.clone(),
        threads: Some(n_opt),
        verbose: a.verbose,
        force: a.force,
        min_free_memory: a.min_free_memory,
        moves: a.moves.clone(),
    };
    // the tiers (without --tier / --opt-tier: the largest whose tables are all there)
    let dir = table_dir(a.tables.as_deref());
    let tier = two_phase_tier(&tab, &dir)?;
    tab.tier = Some(tier);
    let opt_tier = if let Some(t) = a.opt_tier {
        t
    } else {
        let c = autotier::optimal(&dir, otab.move_set())?;
        eprintln!("{}", c.message);
        c.tier
    };
    otab.tier = Some(opt_tier);
    // Both solvers' tables are held at once: check that they fit before loading either
    // (2 bits per entry; a missing table is generated, which `load` checks on its own).
    set_force(a.force);
    // (an `as` cast saturates: negative or NaN -> 0, too large -> u64::MAX)
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let floor = (a.min_free_memory * 1e9) as u64;
    set_memory_floor(floor);
    let files = tier
        .files(tab.move_set())
        .into_iter()
        .chain(opt_tier.files(otab.move_set()).unwrap_or_default());
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
    if a.verbose {
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
        mode: a.two_phase,
        progress: every,
        print: true,
    };
    dispatch!(tab, |t, _name| {
        let ot = load_opt(&otab)?;
        opt_bounds!(ot, |b| {
            for job in &jobs {
                let j = match job {
                    Ok(j) => j,
                    Err(bad) => {
                        bad.report();
                        continue;
                    }
                };
                if jobs.len() > 1 || a.verbose {
                    println!("# {}", j.label);
                }
                let o = combined::solve(&t, &b, &p_tp, &p_opt, j.cube, j.frame, &cfg)?;
                combined::print_outcome(&o);
            }
            Ok::<(), Error>(())
        })
    })?;
    unread(&jobs)
}
