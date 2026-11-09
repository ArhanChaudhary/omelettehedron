use clap::{Parser, Subcommand};
use fto_solver::bench;
use fto_solver::cube::Cube;
use fto_solver::notation::*;
use fto_solver::search::*;
use fto_solver::tables::*;
use std::io::BufRead;
use std::path::PathBuf;

#[derive(Parser)]
#[command(
    name = "fto",
    version,
    about = "Two-phase solver for the face-turning octahedron (FTO)"
)]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(clap::Args, Clone)]
struct TableArgs {
    /// Pruning-table size: small (~260 MB), medium (~1.2 GB) or large (~23 GB, fastest)
    #[arg(long, value_enum, default_value = "medium")]
    tier: Tier,
    /// Directory for pruning tables (default: $FTO_TABLES or ./tables)
    #[arg(long)]
    tables: Option<PathBuf>,
    /// Number of worker threads (default: all cores)
    #[arg(short = 'j', long)]
    threads: Option<usize>,
    /// Print progress information to stderr
    #[arg(short, long)]
    verbose: bool,
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
        /// Scramble in FTO notation (U R F L B BR BL D, ' for anticlockwise)
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
        #[arg(long, default_value_t = 1)]
        seed: u64,
        #[command(flatten)]
        tab: TableArgs,
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
        #[arg(long, default_value_t = 2026)]
        seed: u64,
        #[command(flatten)]
        stop: StopArgs,
        /// write the full results as JSON to this file
        #[arg(long)]
        json: Option<PathBuf>,
        #[command(flatten)]
        tab: TableArgs,
    },
}

fn pool(threads: Option<usize>) -> rayon::ThreadPool {
    let mut b = rayon::ThreadPoolBuilder::new();
    if let Some(t) = threads {
        b = b.num_threads(t);
    }
    b.build().expect("thread pool")
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

fn run_solve<I: AeIndex>(
    t: &Tables<I>,
    scrambles: &[(Cube, Frame, String)],
    limits: &Limits,
    print: Report,
    verbose: bool,
) {
    for (cube, frame, label) in scrambles {
        if scrambles.len() > 1 || verbose {
            println!("# {label}");
        }
        let o = solve(t, *cube, *frame, limits, print);
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
    }
}

macro_rules! dispatch {
    ($tab:expr, |$t:ident, $name:ident| $body:expr) => {{
        let dir = table_dir($tab.tables.as_deref());
        match $tab.tier {
            Tier::Small => {
                let $t = Tables::load(&dir, AeSym::new(false), $tab.verbose).expect("tables");
                let $name = "small";
                $body
            }
            Tier::Medium => {
                let $t = Tables::load(&dir, AeSym::new(true), $tab.verbose).expect("tables");
                let $name = "medium";
                $body
            }
            Tier::Large => {
                let $t = Tables::load(&dir, AeRelabel::new(), $tab.verbose).expect("tables");
                let $name = "large";
                $body
            }
        }
    }};
}

fn main() {
    let cli = Cli::parse();
    match cli.cmd {
        Cmd::Solve {
            scramble,
            stop,
            print,
            random,
            seed,
            tab,
        } => {
            let limits = bench::limits_with_time(stop.time, stop.max_length, stop.max_solutions);
            let positions: Vec<(Cube, Frame, String)> = if let Some(n) = random {
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
            let p = pool(tab.threads);
            p.install(|| {
                dispatch!(tab, |t, _name| run_solve(
                    &t,
                    &positions,
                    &limits,
                    print,
                    tab.verbose
                ))
            });
        }
        Cmd::Gen { tab } => {
            let p = pool(tab.threads);
            let mut tab = tab;
            tab.verbose = true;
            p.install(|| {
                dispatch!(tab, |t, name| {
                    eprintln!(
                        "tier {name}: AE {:.1} MB, CB {:.1} MB; phase-2 depths {:?}",
                        t.ae.bytes() as f64 / 1e6,
                        t.cb.tab.bytes() as f64 / 1e6,
                        t.p2.histogram()
                    );
                    eprintln!("AE depth histogram {:?}", t.ae.hist);
                    eprintln!("CB depth histogram {:?}", t.cb.tab.hist);
                })
            });
        }
        Cmd::Bench {
            scrambles,
            seed,
            stop,
            json,
            tab,
        } => {
            let limits = bench::limits_with_time(stop.time, stop.max_length, stop.max_solutions);
            let p = pool(tab.threads);
            p.install(|| {
                dispatch!(tab, |t, name| {
                    let s = bench::run(&t, name, scrambles, seed, &limits, tab.verbose);
                    bench::print_summary(&s);
                    if let Some(path) = &json {
                        std::fs::write(path, serde_json::to_string_pretty(&s).unwrap())
                            .expect("write json");
                    }
                })
            });
        }
    }
}
