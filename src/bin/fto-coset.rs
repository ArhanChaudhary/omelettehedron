#![deny(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::unreachable,
    clippy::todo,
    clippy::unimplemented
)]
#![warn(clippy::pedantic)]
// Pedantic lints that would only add noise here (as in the library):
#![allow(
    // statistics are reported as f64
    clippy::cast_precision_loss,
    // coset indices are below 933,120 and lengths below 32
    clippy::cast_possible_truncation,
    // names follow the maths (t = tables, c = coset, p = index, z = element, ...)
    clippy::many_single_char_names,
    // command-line flags
    clippy::struct_excessive_bools,
    clippy::too_many_lines
)]
//! `fto-coset`: the coset solver for `H_A` = <U, F, BR, BL> ([`omelettehedron::coset`]), as a
//! separate program (the `omelettehedron` command set is fixed).
//!
//! * default: fill the coset of each position to a depth (23 by default), printing per depth
//!   how many of its 933,120 elements are proven within it, then the distribution of the bounds
//!   and the elements not covered;
//! * `--phase1-depth`: only the exact phase-1 depth of each position (the coset's distance from
//!   `H_A`), one position per thread, for sampling many cosets.

use clap::Parser;
use omelettehedron::coset::enumerate::Phase1;
use omelettehedron::coset::space::{histogram, N_COSET, UNSET};
use omelettehedron::coset::{self, Options, Searched};
use omelettehedron::cube::Cube;
use omelettehedron::light::AeLight;
use omelettehedron::notation::{
    apply_physical, format_moves, is_solved_physical, parse_moves, scramble_to_cube,
    solution_to_physical, Frame,
};
use omelettehedron::pattern;
use omelettehedron::tables::{table_dir, AeIndex, AeRelabel, AeSym, Tables, Tier};
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use std::io::{Read, Write};
use std::path::PathBuf;
use std::process::ExitCode;
use std::time::Instant;

type Error = Box<dyn std::error::Error + Send + Sync>;

#[derive(Parser)]
#[command(
    name = "fto-coset",
    version,
    about = "Coset solver for the FTO's tetrahedral subgroup <U, F, BR, BL>",
    long_about = "Coset solver for the FTO's tetrahedral subgroup H_A = <U, F, BR, BL>.\n\nThe coset \
                  of a position is the 933,120 positions that share its phase-1 state. Every \
                  phase-1 sequence of each length (ending in a move outside H_A) marks the \
                  H_A state it ends in; a prepass extends the marks by one H_A move per depth. \
                  Printed per depth: how many elements are within it (exact distances up to the \
                  longest phase-1 length enumerated, upper bounds above)."
)]
struct Args {
    /// A scramble (FTO notation) or a pattern (72 stickers). Without one: one position per
    /// line on stdin (`#` starts a comment, used as the position's name)
    position: Vec<String>,
    /// Read patterns from stdin, each on one line or several (lines are joined until they hold
    /// 72 stickers)
    #[arg(long)]
    patterns: bool,
    /// Use this many seeded random positions (uniformly random cosets)
    #[arg(long)]
    random: Option<usize>,
    /// Use the coset of `H_B` = <R, L, B, D> instead: each position is first conjugated by the
    /// 90-degree rotation about the fixed corner's axis, which swaps the two face classes
    #[arg(long)]
    class_b: bool,
    /// Seed of --random, --check and --samples
    #[arg(long, default_value_t = 1)]
    seed: u64,

    /// Fill the coset to this depth
    #[arg(short, long, default_value_t = 23)]
    depth: usize,
    /// Longest phase-1 sequences to enumerate (default: the depth). Deeper, only the prepass
    #[arg(short = 'S', long)]
    max_search: Option<usize>,
    /// Enumerate every length up to --max-search even after every element is within the depth
    /// (exact distances further up)
    #[arg(long)]
    full: bool,
    /// Groups the tasks of a length are searched in: the coverage is checked after each, so
    /// the last length searched stops sooner with more
    #[arg(long, default_value_t = 16)]
    chunks: usize,

    /// Only the exact phase-1 depth of each position (one thread per position)
    #[arg(long)]
    phase1_depth: bool,
    /// With --phase1-depth: give up above this length
    #[arg(long, default_value_t = 20)]
    max_phase1: usize,

    /// Write the elements not within the depth to this file: pattern and scramble
    #[arg(long)]
    uncovered: Option<PathBuf>,
    /// Print this many random elements of each coset with a solution (keeps the phase-1
    /// sequences; each solution is checked on the sticker model)
    #[arg(long, default_value_t = 0)]
    samples: usize,
    /// Check this many random elements: the bound of each must equal the shortest solution
    /// through the phase-1 sequences found (each sequence plus its phase-2 completion,
    /// computed on the element itself), and those solutions must solve it
    #[arg(long, default_value_t = 0)]
    check: usize,

    /// Mark every phase-1 sequence that ends in `H_A` as a pattern, not only those that end in
    /// it with the B-centres told apart (unsound for the other elements: for showing why)
    #[arg(long)]
    patterns_only: bool,
    /// Write this many random elements of each coset to --elements-out (one pattern per line,
    /// with the index and bound as a comment): input for other solvers
    #[arg(long, default_value_t = 0)]
    elements: usize,
    /// File for --elements (appended to)
    #[arg(long)]
    elements_out: Option<PathBuf>,
    /// Print these elements (by index): bound, pattern, scramble
    #[arg(long, value_delimiter = ',')]
    element: Vec<u32>,
    /// Number of worker threads (default: all cores)
    #[arg(short = 'j', long)]
    threads: Option<usize>,
    /// Two-phase table tier (small, medium, large or huge; default: the largest present)
    #[arg(long, value_enum)]
    tier: Option<Tier>,
    /// Directory of the pruning tables (default: `$FTO_TABLES` or ./tables)
    #[arg(long)]
    tables: Option<PathBuf>,
    /// Less on stderr (no per-depth progress)
    #[arg(short, long)]
    quiet: bool,
}

/// A position to work on.
struct Job {
    label: String,
    cube: Cube,
}

fn main() -> ExitCode {
    match run(&Args::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

fn run(args: &Args) -> Result<(), Error> {
    let mut jobs = read_jobs(args)?;
    if args.class_b {
        for j in &mut jobs {
            j.cube = j
                .cube
                .conjugate(&omelettehedron::cube::SWAP_AB)
                .ok_or("a position could not be conjugated")?;
        }
    }
    let mut b = rayon::ThreadPoolBuilder::new();
    if let Some(t) = args.threads {
        b = b.num_threads(t);
    }
    let pool = b.build()?;
    let dir = table_dir(args.tables.as_deref());
    let moves = omelettehedron::cube::MoveSet::FACE;
    let tier = if let Some(t) = args.tier {
        t
    } else {
        let c = omelettehedron::autotier::two_phase(&dir, moves, "--tier")?;
        if !args.quiet {
            eprintln!("{}", c.message);
        }
        c.tier
    };
    let verbose = !args.quiet;
    let err = |e: std::io::Error| -> Error { format!("tables in {}: {e}", dir.display()).into() };
    match tier {
        Tier::Tiny => {
            let t = Tables::load(&dir, AeLight::new(), moves, verbose).map_err(err)?;
            pool.install(|| work(args, &t, &jobs))
        }
        Tier::Small => {
            let t = Tables::load(&dir, AeSym::new(false), moves, verbose).map_err(err)?;
            pool.install(|| work(args, &t, &jobs))
        }
        Tier::Medium => {
            let t = Tables::load(&dir, AeSym::new(true), moves, verbose).map_err(err)?;
            pool.install(|| work(args, &t, &jobs))
        }
        Tier::Large => {
            let t = Tables::load_with(&dir, AeSym::new(true), true, moves, verbose).map_err(err)?;
            pool.install(|| work(args, &t, &jobs))
        }
        Tier::Huge => {
            let t = Tables::load_with(&dir, AeRelabel::new(), true, moves, verbose).map_err(err)?;
            pool.install(|| work(args, &t, &jobs))
        }
    }
}

/// The positions: from the command line, --random, or stdin.
fn read_jobs(args: &Args) -> Result<Vec<Job>, Error> {
    if let Some(n) = args.random {
        return Ok(omelettehedron::bench::random_positions(n, args.seed)
            .into_iter()
            .enumerate()
            .map(|(k, cube)| Job {
                label: format!("random{k}"),
                cube,
            })
            .collect());
    }
    let mut out = vec![];
    if !args.position.is_empty() {
        let text = args.position.join(" ");
        out.push(Job {
            cube: read_position(&text)?,
            label: text,
        });
        return Ok(out);
    }
    let mut text = String::new();
    std::io::stdin().lock().read_to_string(&mut text)?;
    if args.patterns {
        for (k, r) in pattern::records(&text).iter().enumerate() {
            let p = pattern::parse(&r.text)?;
            let label = text
                .lines()
                .nth(r.line - 1)
                .and_then(|l| l.split_once('#'))
                .map_or_else(
                    || format!("pattern{}", k + 1),
                    |(_, c)| c.trim().to_string(),
                );
            out.push(Job {
                cube: pattern::position(&p)?,
                label,
            });
        }
    } else {
        for (k, line) in text.lines().enumerate() {
            let (body, name) = match line.split_once('#') {
                Some((b, c)) => (b.trim(), Some(c.trim())),
                None => (line.trim(), None),
            };
            if body.is_empty() {
                continue;
            }
            let label = name.map_or_else(|| format!("line{}", k + 1), str::to_string);
            let cube = read_position(body).map_err(|e| format!("line {}: {e}", k + 1))?;
            out.push(Job { label, cube });
        }
    }
    if out.is_empty() {
        return Err("no position given".into());
    }
    Ok(out)
}

/// A scramble, or else a pattern.
fn read_position(text: &str) -> Result<Cube, Error> {
    match parse_moves(text) {
        Ok(m) => Ok(scramble_to_cube(&m).0),
        Err(move_error) => match pattern::parse(text) {
            Ok(p) => Ok(pattern::position(&p)?),
            Err(e) => Err(format!("neither a scramble ({move_error}) nor a pattern ({e})").into()),
        },
    }
}

fn work<I: AeIndex>(args: &Args, t: &Tables<I>, jobs: &[Job]) -> Result<(), Error> {
    if args.phase1_depth {
        return phase1_depths(args, t, jobs);
    }
    for job in jobs {
        fill(args, t, job)?;
    }
    Ok(())
}

/// `--phase1-depth`: one line per position (in input order, printed in groups as they end).
fn phase1_depths<I: AeIndex>(args: &Args, t: &Tables<I>, jobs: &[Job]) -> Result<(), Error> {
    let start = Instant::now();
    let mut out = std::io::stdout().lock();
    writeln!(
        out,
        "# name\tphase1_depth\tlower_bound\tnodes\tseconds\tsequence"
    )?;
    let mut hist = vec![0u64; 32];
    let mut failed = 0u64;
    let group = (rayon::current_num_threads() * 64).max(1);
    for (gi, part) in jobs.chunks(group).enumerate() {
        let rows: Vec<Result<String, String>> = part
            .par_iter()
            .map(|job| {
                let s0 = Instant::now();
                let p1 = Phase1::new(t, job.cube, None)?;
                let lb = p1.lower_bound();
                let (d, stats) = p1.depth(args.max_phase1);
                let secs = s0.elapsed().as_secs_f64();
                Ok(match d {
                    Some((len, seq)) => {
                        // (the sequence must end in H_A: checked in full by the enumeration)
                        let phys = solution_to_physical(Frame::identity(), &seq, false);
                        format!(
                            "{}\t{len}\t{lb}\t{}\t{secs:.3}\t{}",
                            job.label,
                            stats.nodes,
                            format_moves(&phys)
                        )
                    }
                    None => format!(
                        "{}\t>{}\t{lb}\t{}\t{secs:.3}\t",
                        job.label, args.max_phase1, stats.nodes
                    ),
                })
            })
            .collect();
        for r in rows {
            let line = r?;
            match line
                .split('\t')
                .nth(1)
                .and_then(|d| d.parse::<usize>().ok())
            {
                Some(d) => {
                    if let Some(h) = hist.get_mut(d) {
                        *h += 1;
                    }
                }
                None => failed += 1,
            }
            writeln!(out, "{line}")?;
        }
        out.flush()?;
        if !args.quiet {
            eprintln!(
                "{} of {} positions, {:.1} s",
                (gi * group + part.len()).min(jobs.len()),
                jobs.len(),
                start.elapsed().as_secs_f64()
            );
        }
    }
    while hist.last() == Some(&0) {
        hist.pop();
    }
    writeln!(
        out,
        "# phase-1 depth distribution ({} positions):",
        jobs.len()
    )?;
    for (d, &n) in hist.iter().enumerate() {
        if n > 0 {
            writeln!(out, "# {d:>2}  {n}")?;
        }
    }
    if failed > 0 {
        writeln!(out, "# above {}: {failed}", args.max_phase1)?;
    }
    Ok(())
}

/// Fill the coset of one position and report it.
fn fill<I: AeIndex>(args: &Args, t: &Tables<I>, job: &Job) -> Result<(), Error> {
    let opt = Options {
        target: args.depth,
        max_search: args.max_search.unwrap_or(args.depth),
        stop_when_covered: !args.full,
        keep_paths: args.samples > 0 || args.check > 0,
        chunks: args.chunks,
        labelled: !args.patterns_only,
    };
    let start = Instant::now();
    let mut out = std::io::stdout().lock();
    writeln!(out, "coset of {}", job.label)?;
    writeln!(
        out,
        "  pattern {}",
        pattern::format(&pattern::of_cube(&job.cube))
    )?;
    writeln!(
        out,
        "  depth  searched  sequences  as_pattern    from_seq  from_prepass     covered  proj_max         nodes  seconds"
    )?;
    out.flush()?;
    let quiet = args.quiet;
    let c = coset::solve(t, &job.cube, &opt, |l| {
        if !quiet {
            eprintln!(
                "depth {:>2}: {:>7} covered, {} sequences, {:.1} s",
                l.depth, l.covered, l.sequences, l.seconds
            );
        }
    })?;
    for l in &c.levels {
        let s = match l.searched {
            Searched::Full => "full",
            Searched::Partial => "partial",
            Searched::No => "-",
        };
        let pm = l.projected_max.map_or("-".to_string(), |m| m.to_string());
        writeln!(
            out,
            "  {:>5}  {:>8}  {:>9}  {:>10}  {:>10}  {:>12}  {:>10}  {:>8}  {:>12}  {:>7.2}",
            l.depth,
            s,
            l.sequences,
            l.stats.in_h,
            l.from_sequences,
            l.from_prepass,
            l.covered,
            pm,
            l.stats.nodes,
            l.seconds
        )?;
    }
    let total = start.elapsed().as_secs_f64();
    let p1 = c
        .phase1_depth
        .map_or_else(|| "not reached".to_string(), |d| d.to_string());
    let exact = c
        .exact_through
        .map_or_else(|| "-".to_string(), |d| d.to_string());
    writeln!(
        out,
        "  phase-1 depth {p1}; distances exact through {exact}; {total:.1} s"
    )?;
    let h = histogram(&c.bounds);
    let mut cum = 0u64;
    writeln!(
        out,
        "  distribution (exact through {exact}, upper bounds above):"
    )?;
    for (d, &n) in h.iter().enumerate() {
        cum += n;
        if n > 0 {
            writeln!(out, "    {d:>2}  {n:>7}  (<= {d}: {cum})")?;
        }
    }
    let uncovered = c.uncovered(args.depth);
    writeln!(
        out,
        "  not within {}: {} of {}",
        args.depth,
        uncovered.len(),
        N_COSET
    )?;
    if let Some(path) = &args.uncovered {
        write_uncovered(path, t, &c, &uncovered, job)?;
        writeln!(out, "  (listed in {})", path.display())?;
    }
    for &p in &args.element {
        let Some((word, _)) = c.element(&t.p2, p) else {
            return Err("an element could not be made".into());
        };
        let phys = solution_to_physical(Frame::identity(), &word, false);
        writeln!(
            out,
            "  element {p}: bound {}  {}  {}",
            c.bounds.get(p as usize).copied().unwrap_or(UNSET),
            pattern::format(&pattern::of_moves(&phys)),
            format_moves(&phys)
        )?;
    }
    if args.elements > 0 {
        let path = args
            .elements_out
            .as_ref()
            .ok_or("--elements needs --elements-out")?;
        let mut f = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(path)?;
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(args.seed ^ 0xe1e);
        for _ in 0..args.elements {
            let p: u32 = rng.gen_range(0..N_COSET as u32);
            let Some((word, _)) = c.element(&t.p2, p) else {
                return Err("an element could not be made".into());
            };
            let phys = solution_to_physical(Frame::identity(), &word, false);
            writeln!(
                f,
                "{}  # {} element {p} bound {}",
                pattern::format(&pattern::of_moves(&phys)),
                job.label,
                c.bounds[p as usize]
            )?;
        }
    }
    if args.check > 0 {
        check(args, t, &c, &mut out)?;
    }
    if args.samples > 0 {
        samples(args, t, &c, &mut out)?;
    }
    out.flush()?;
    Ok(())
}

/// Append the elements not covered to `path`: bound (or `-`), pattern, scramble.
fn write_uncovered<I: AeIndex>(
    path: &PathBuf,
    t: &Tables<I>,
    c: &coset::Coset,
    uncovered: &[u32],
    job: &Job,
) -> Result<(), Error> {
    let mut f = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)?;
    writeln!(
        f,
        "# coset of {}: {} not covered",
        job.label,
        uncovered.len()
    )?;
    for &p in uncovered {
        let Some((word, _)) = c.element(&t.p2, p) else {
            return Err("an element could not be made".into());
        };
        let b = c.bounds[p as usize];
        let b = if b == UNSET {
            "-".to_string()
        } else {
            b.to_string()
        };
        let phys = solution_to_physical(Frame::identity(), &word, false);
        writeln!(
            f,
            "{b}\t{}\t{}",
            pattern::format(&pattern::of_moves(&phys)),
            format_moves(&phys)
        )?;
    }
    Ok(())
}

/// Whether physical `solution` solves the position that `scramble` makes from solved, on the
/// sticker model (independently of the internal cube model).
fn replays(
    scramble: &[omelettehedron::notation::PhysMove],
    solution: &[omelettehedron::notation::PhysMove],
) -> bool {
    let solved: omelettehedron::cube::Stickers =
        omelettehedron::cube::solved_stickers(omelettehedron::cube::defs());
    is_solved_physical(&apply_physical(
        &apply_physical(&solved, scramble),
        solution,
    ))
}

/// `--check`: random elements' bounds against the per-element computation.
fn check<I: AeIndex>(
    args: &Args,
    t: &Tables<I>,
    c: &coset::Coset,
    out: &mut impl Write,
) -> Result<(), Error> {
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(args.seed);
    let (mut ok, mut bad) = (0u64, 0u64);
    for _ in 0..args.check {
        let p: u32 = rng.gen_range(0..N_COSET as u32);
        let Some((word, z)) = c.element(&t.p2, p) else {
            return Err("an element could not be made".into());
        };
        let b = c.bounds[p as usize];
        let sol = c.solution(&t.p2, &z);
        let len = sol.as_ref().map_or(UNSET, |s| s.len() as u8);
        let solves = sol.as_ref().is_none_or(|s| {
            z.apply_seq(s) == Cube::solved() && {
                let scr = solution_to_physical(Frame::identity(), &word, false);
                // (the solution starts in the frame the scramble ends in)
                let fr = omelettehedron::notation::scramble_to_cube(&scr).1;
                replays(&scr, &solution_to_physical(fr, s, false))
            }
        });
        if len == b && solves {
            ok += 1;
        } else {
            bad += 1;
            writeln!(out, "  CHECK FAILED: element {p}: bound {b}, through the sequences {len}, solves {solves}")?;
        }
    }
    writeln!(out, "  check: {ok} of {} elements agree", ok + bad)?;
    if bad > 0 {
        return Err(format!("{bad} elements disagree").into());
    }
    Ok(())
}

/// `--samples`: random covered elements with their solutions.
fn samples<I: AeIndex>(
    args: &Args,
    t: &Tables<I>,
    c: &coset::Coset,
    out: &mut impl Write,
) -> Result<(), Error> {
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(args.seed ^ 0x5a5a);
    writeln!(
        out,
        "  samples (bound, solution length, pattern, scramble, solution):"
    )?;
    for _ in 0..args.samples {
        let p: u32 = rng.gen_range(0..N_COSET as u32);
        let Some((word, z)) = c.element(&t.p2, p) else {
            return Err("an element could not be made".into());
        };
        let scr = solution_to_physical(Frame::identity(), &word, false);
        let Some(sol) = c.solution(&t.p2, &z) else {
            writeln!(
                out,
                "    -  {}  (no solution through the sequences kept)",
                pattern::format(&pattern::of_moves(&scr))
            )?;
            continue;
        };
        let fr = scramble_to_cube(&scr).1;
        let phys = solution_to_physical(fr, &sol, false);
        let ok = z.apply_seq(&sol) == Cube::solved() && replays(&scr, &phys);
        writeln!(
            out,
            "    {}  {}  {}  {}  {}{}",
            c.bounds[p as usize],
            sol.len(),
            // (the pattern as the scramble leaves the puzzle, which the solution solves)
            pattern::format(&pattern::of_moves(&scr)),
            format_moves(&scr),
            format_moves(&phys),
            if ok { "" } else { "  (DOES NOT SOLVE)" }
        )?;
    }
    Ok(())
}
