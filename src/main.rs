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
use omelettehedron::cube::{Cube, MoveSet};
use omelettehedron::notation::{
    format_moves, invert_moves, parse_moves, scramble_to_cube, solves_pure, Frame, PhysMove,
};
use omelettehedron::optimal::search::{solve as opt_solve, Limits as OptLimits, Stop};
use omelettehedron::optimal::tables::{OptTables, OptTier};
use omelettehedron::pattern;
use omelettehedron::search::{solve, Limits, Report};
use omelettehedron::tables::{table_dir, AeIndex, AeRelabel, AeSym, Tables, Tier};
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
    /// Pruning-table size (required; there is no default): small (~135 MB), medium (~1.0 GB),
    /// large (~8 GB) or huge (~30 GB). Larger tiers search fewer nodes; large is the fastest.
    #[arg(long, value_enum)]
    tier: Tier,
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
        /// earlier run of the same position logged ("pass N complete: no solution of at most N
        /// moves" with `FTO_DEBUG=1` gives N + 1). The search finds a first solution, then starts
        /// its target-length passes here. A wrong value can give a wrong answer
        #[arg(long)]
        lower_bound: Option<usize>,
        /// Each given sequence is a solution: solve the position it solves (the sequence
        /// inverted, applied to the solved puzzle) instead of the position it makes
        #[arg(long)]
        inverse: bool,
        /// Solve a pattern instead of a scramble: the 72 sticker colours as you see them, in
        /// the order of `omelettehedron pattern --help`. Also prints a scramble that makes it
        #[arg(long, conflicts_with_all = ["scramble", "random", "inverse"])]
        pattern: Option<String>,
        /// With --pattern: your solved puzzle as you usually hold it to scramble, as the colours
        /// of U R F L B BL D BR (8 characters, the letters of the pattern). The scramble then
        /// starts from that hold
        #[arg(long, requires = "pattern")]
        usual: Option<String>,
        #[command(flatten)]
        tab: TableArgs,
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
}

/// Table options of the optimal solver.
#[derive(clap::Args, Clone)]
struct OptTableArgs {
    /// Table size (required; there is no default): small (~195 MB), medium (~3.4 GB, fewer
    /// nodes) or large (~24.4 GB, fewest nodes)
    #[arg(long, value_enum)]
    tier: OptTier,
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
        #[command(flatten)]
        tab: OptTableArgs,
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
    let usual = match usual {
        Some(u) => {
            let v: Vec<char> = u
                .chars()
                .filter(|c| !c.is_whitespace() && *c != '/')
                .collect();
            let u: [char; 8] = v
                .as_slice()
                .try_into()
                .map_err(|_| "--usual: give the 8 colours of U R F L B BL D BR")?;
            Some(u)
        }
        None => None,
    };
    Ok(PatternInput {
        pattern: p,
        usual,
        cube,
    })
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
) {
    for (cube, frame, label) in scrambles {
        if scrambles.len() > 1 || verbose {
            println!("# {label}");
        }
        let o = solve(t, *cube, *frame, limits, print);
        if let (Some(pi), Some(sol)) = (input, &o.solution) {
            print_pattern_scramble(pi, sol);
        }
        if verbose || print == Report::Improvements {
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
        omelettehedron::preflight::set_force($tab.force);
        // The large tier needs 8 GB and the huge tier's index has 9.1e10 entries: 32-bit
        // targets such as WebAssembly (4 GB) support the small and medium tiers (the medium
        // index has 4.13e9 entries, below 2^32).
        if cfg!(target_pointer_width = "32") && matches!($tab.tier, Tier::Large | Tier::Huge) {
            return Err(
                "only --tier small and medium are supported on 32-bit targets such as \
                        WebAssembly"
                    .into(),
            );
        }
        match $tab.tier {
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

/// The positions to solve: `--random` positions, or the given scrambles (or stdin), inverted
/// with `--inverse`.
fn read_positions(
    scramble: &[String],
    random: Option<usize>,
    seed: u64,
    inverse: bool,
) -> Result<Vec<(Cube, Frame, String)>, Error> {
    if let Some(n) = random {
        return Ok(bench::random_positions(n, seed)
            .into_iter()
            .enumerate()
            .map(|(i, c)| {
                (
                    c,
                    Frame::identity(),
                    format!("random position {i} (seed {seed})"),
                )
            })
            .collect());
    }
    let mut v = vec![];
    for s in read_scrambles(scramble) {
        let m = parse_moves(&s)?;
        if inverse {
            let m = invert_moves(&m);
            let (c, f) = scramble_to_cube(&m);
            v.push((c, f, format_moves(&m)));
        } else {
            let (c, f) = scramble_to_cube(&m);
            v.push((c, f, s));
        }
    }
    if v.is_empty() {
        return Err("no scramble given".into());
    }
    Ok(v)
}

/// Load (generating if needed) the optimal solver's tables.
fn load_opt(tab: &OptTableArgs) -> Result<OptTables, Error> {
    let dir = table_dir(tab.tables.as_deref());
    omelettehedron::preflight::set_force(tab.force);
    // (an `as` cast saturates: negative or NaN -> 0, too large -> u64::MAX)
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let floor = (tab.min_free_memory * 1e9) as u64;
    omelettehedron::preflight::set_memory_floor(floor);
    OptTables::load(&dir, tab.tier, tab.verbose).map_err(|e| table_error(&dir, &e))
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
            tab,
        } => {
            let limits = stop.limits()?;
            let positions = read_positions(&scramble, random, seed, inverse)?;
            let p = pool(tab.threads)?;
            p.install(|| {
                let t = load_opt(&tab)?;
                opt_bounds!(t, |b| {
                    for (cube, frame, label) in &positions {
                        if positions.len() > 1 || tab.verbose {
                            println!("# {label}");
                        }
                        let o = opt_solve(&b, *cube, *frame, &limits, tab.verbose);
                        for s in &o.solutions {
                            println!("{}", format_moves(s));
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
                        if tab.verbose {
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
                    Ok(())
                })
            })
        }
        OptCmd::Gen { tab } => {
            let p = pool(tab.threads)?;
            let mut tab = tab;
            tab.verbose = true;
            p.install(|| {
                let t = load_opt(&tab)?;
                eprintln!(
                    "optimal solver, tier {}: {:.1} MB of tables",
                    t.tier.name(),
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
            usual,
            tab,
        } => {
            let mut limits =
                bench::limits_with_time(stop.time, stop.max_length, stop.max_solutions)
                    .map_err(|e| format!("--time: {e}"))?;
            limits.lower_bound = lower_bound;
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
                    input.as_ref()
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
                    if let Some(cbs) = &t.cbs {
                        eprintln!("CBE-to-solved depth histogram {:?}", cbs.hist);
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
    }
}
