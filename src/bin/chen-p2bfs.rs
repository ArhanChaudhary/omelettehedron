#![deny(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::unreachable,
    clippy::todo,
    clippy::unimplemented
)]
#![warn(clippy::pedantic)]
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    // (the GB and minute options are positive)
    clippy::cast_sign_loss,
    clippy::cast_lossless,
    // names follow the library's (c = class, x = B-centres, p = corners, ...)
    clippy::many_single_char_names,
    clippy::similar_names,
    clippy::too_many_lines
)]
//! `chen-p2bfs`: the exact breadth-first search of the block split's phase-2 space
//! ([`omelettehedron::chen::p2bfs`]), and the checks that go with it.
//!
//! * `run`: start or resume a search in a directory;
//! * `status`: the layers made so far, the totals and the `H'` classes;
//! * `reference`, `positions`, `tables`, `orbits`: independent counts to compare a search
//!   with;
//! * `sample`: exact distances of slots of a layer by IDA* with the phase-2 tables;
//! * `antipodes`: the positions of a layer (the last one by default), with solutions.

use anyhow::{bail, Context, Result};
use clap::{Parser, Subcommand, ValueEnum};
use omelettehedron::chen::p2bfs::engine::{
    bench_push, canonical, dense_plan, for_slots, Bfs, Config, Layout, Manifest, Outcome,
};
use omelettehedron::chen::p2bfs::instance::{BView, CornerView, Instance, Kind, NO_EDGE};
use omelettehedron::chen::p2bfs::reference::{position_bfs, slot_bfs};
use omelettehedron::chen::p2bfs::solve::{write, PState, Solver};
use omelettehedron::chen::tables::{ChenTables, TableSet};
use omelettehedron::notation::format_moves;
use rand::{Rng, SeedableRng};
use std::path::{Path, PathBuf};

#[derive(Parser)]
#[command(about = "Exact breadth-first search of the block split's phase-2 space")]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
    /// Worker threads (default: all)
    #[arg(short = 'j', long, global = true)]
    threads: Option<usize>,
}

#[derive(Clone, Copy, ValueEnum)]
enum Corners {
    Full,
    Perm,
    None,
}

#[derive(Clone, Copy, ValueEnum)]
enum Bs {
    Full,
    Unlabelled,
    None,
}

#[derive(clap::Args, Clone, Copy)]
struct InstanceArgs {
    /// Corners kept (full: 11,520 states; perm: 360; none)
    #[arg(long, value_enum, default_value = "full")]
    corners: Corners,
    /// B-centres kept (full: 369,600 colourings; unlabelled: 15,400; none)
    #[arg(long, value_enum, default_value = "full")]
    b: Bs,
}

impl InstanceArgs {
    fn kind(self) -> Kind {
        Kind {
            corners: match self.corners {
                Corners::Full => CornerView::Full,
                Corners::Perm => CornerView::Perm,
                Corners::None => CornerView::None,
            },
            b: match self.b {
                Bs::Full => BView::Full,
                Bs::Unlabelled => BView::Unlabelled,
                Bs::None => BView::None,
            },
        }
    }
}

#[derive(Subcommand)]
enum Cmd {
    /// Start or resume a search
    Run {
        /// Search directory
        #[arg(long)]
        dir: PathBuf,
        #[command(flatten)]
        inst: InstanceArgs,
        /// Memory for a part of a layer, GB (arena and directory; the tables and buffers come
        /// on top, about 1 GB)
        #[arg(long, default_value_t = 55.0)]
        ram_gb: f64,
        /// Stop once this depth is made
        #[arg(long)]
        max_depth: Option<usize>,
        /// Stop after this many parts in this run
        #[arg(long)]
        max_parts: Option<usize>,
        /// Stop after the first part that ends past this many minutes
        #[arg(long)]
        max_minutes: Option<f64>,
        /// Keep the layers of depths 0..=KEEP
        #[arg(long, default_value_t = 0)]
        keep: usize,
        /// Directories for the layer files (repeat for several disks: class c on the
        /// (c mod n)-th); default the search directory. Fixed when the search starts.
        #[arg(long)]
        store: Vec<PathBuf>,
        /// Threads reading layer files in a push
        #[arg(long, default_value_t = 4)]
        readers: usize,
        /// No progress lines
        #[arg(short, long)]
        quiet: bool,
    },
    /// The layers made so far
    Status {
        /// Search directory
        #[arg(long)]
        dir: PathBuf,
    },
    /// Slots per depth of an instance by the in-memory table search (must fit in memory: 2
    /// bits per slot)
    Reference {
        #[command(flatten)]
        inst: InstanceArgs,
        /// Compare with this search directory
        #[arg(long)]
        compare: Option<PathBuf>,
    },
    /// Positions per depth by a plain search without symmetry (the first depths)
    Positions {
        /// Last depth
        #[arg(long, default_value_t = 6)]
        depth: usize,
        /// Compare with this search directory (of the full instance)
        #[arg(long)]
        compare: Option<PathBuf>,
    },
    /// The phase-2 tables' histograms (AC: corners only; PX: B-centres only)
    Tables {
        /// Table directory
        #[arg(long)]
        tables: PathBuf,
        /// Compare AC with this search directory (corners full, B-centres none)
        #[arg(long)]
        ac: Option<PathBuf>,
        /// Compare PX with this search directory (corners none, B-centres full)
        #[arg(long)]
        px: Option<PathBuf>,
    },
    /// Positions up to symmetry of a layer, counted one slot per orbit (against the
    /// manifest's count by fixed points)
    Orbits {
        /// Search directory
        #[arg(long)]
        dir: PathBuf,
        /// Depths (default: all kept)
        #[arg(long)]
        depth: Vec<usize>,
    },
    /// Exact distances of random slots of a layer by IDA* (full instance)
    Sample {
        /// Search directory
        #[arg(long)]
        dir: PathBuf,
        /// Phase-2 table directory
        #[arg(long)]
        tables: PathBuf,
        /// Layer
        #[arg(long)]
        depth: usize,
        /// Slots
        #[arg(short, long, default_value_t = 100)]
        n: usize,
        /// Random seed
        #[arg(long, default_value_t = 1)]
        seed: u64,
    },
    /// The parts of the busiest step for a memory budget, and the layer files they read
    Plan {
        #[command(flatten)]
        inst: InstanceArgs,
        /// Memory for a part, GB
        #[arg(long, default_value_t = 55.0)]
        ram_gb: f64,
    },
    /// Time the push of random blocks (in memory) into a part of all classes
    BenchPush {
        #[command(flatten)]
        inst: InstanceArgs,
        /// Source classes
        #[arg(long, default_value_t = 1611)]
        classes: usize,
        /// Blocks per source class
        #[arg(long, default_value_t = 200)]
        blocks: usize,
        /// Density of the blocks
        #[arg(long, default_value_t = 0.4)]
        density: f64,
        /// The part: classes lo..hi
        #[arg(long, default_value_t = 0)]
        lo: usize,
        #[arg(long, default_value_t = 1611)]
        hi: usize,
    },
    /// Exact distances of uniformly random phase-2 positions by IDA* (an estimate of the
    /// distribution)
    Random {
        /// Phase-2 table directory
        #[arg(long)]
        tables: PathBuf,
        /// Positions
        #[arg(short, long, default_value_t = 20)]
        n: usize,
        /// Random seed
        #[arg(long, default_value_t = 1)]
        seed: u64,
        /// Give up on a position after this many moves
        #[arg(long, default_value_t = 30)]
        max: usize,
    },
    /// The joint check of the two phases (see `omelettehedron::chen::p2bfs::joint`): exact
    /// relaxed phase-1 distances, the geodesic set, membership checks, then every candidate
    /// above the last layer resolved or listed
    Joint {
        /// Search directory (complete, its last layer kept)
        #[arg(long)]
        dir: PathBuf,
        /// Directory of the relaxed phase-1 table and the phase-2 tables
        #[arg(long)]
        tables: PathBuf,
        /// Random positions for the checks of d1 (table, exact array, conjugates)
        #[arg(long, default_value_t = 2000)]
        checks: u64,
        /// Random phase-2 positions solved by IDA* (within 24 moves) that must test outside the
        /// layer
        #[arg(long, default_value_t = 16)]
        outside: usize,
        /// Process only the first N classes of the layer (to time the run)
        #[arg(long)]
        limit: Option<usize>,
        /// Write the unresolved candidates here (tab-separated: pattern of the fixed-frame
        /// position, then the position's coordinates)
        #[arg(long)]
        out: Option<PathBuf>,
        /// Classes of the layer whose first candidates get an explicit solution (phase-1
        /// alternative plus an IDA* phase-2 solution)
        #[arg(long, default_value_t = 100)]
        sample: usize,
        /// Write those solutions here (tab-separated: pattern, phase-1 and phase-2 lengths,
        /// length, solution in face turns)
        #[arg(long)]
        solutions: Option<PathBuf>,
    },
    /// Counts of ascending phase-1 paths and of the candidates they give from the last two
    /// layers (an estimate for a bound one move lower)
    Paths {
        /// Search directory (complete)
        #[arg(long)]
        dir: PathBuf,
        /// Directory of the relaxed phase-1 table
        #[arg(long)]
        tables: PathBuf,
    },
    /// The bad sets above the last layer (positions all of whose optimal phase-1 solutions
    /// end in it), level by level, to a top distance
    Bad {
        /// Search directory (complete)
        #[arg(long)]
        dir: PathBuf,
        /// Directory of the relaxed phase-1 table
        #[arg(long)]
        tables: PathBuf,
        /// Top distance (default: the largest)
        #[arg(long)]
        top: Option<u8>,
        /// Keep only positions on the geodesic set to the largest distance
        #[arg(long)]
        geodesic: bool,
    },
    /// The relaxed phase-1 distance along given solutions (how far each leaves the phase-2
    /// space): a tab-separated file with the pattern in column 5 and the solution in column 9
    Radius {
        /// Directory of the relaxed phase-1 table
        #[arg(long)]
        tables: PathBuf,
        /// The solutions
        #[arg(long)]
        input: PathBuf,
    },
    /// Radius-N phase-2 distances (every position along the solution within relaxed
    /// phase-1 distance N) of the positions of a file (pattern in column 5; a known lower bound,
    /// e.g. the distance with every move allowed, in column 7), by IDA* with the optimal
    /// solver's large tier as lower bounds
    RadiusSolve {
        /// The positions
        #[arg(long)]
        input: PathBuf,
        /// Directory of the relaxed phase-1 table
        #[arg(long)]
        tables: PathBuf,
        /// Directory of the optimal solver's tables (large tier)
        #[arg(long)]
        opt_tables: PathBuf,
        /// Radii (each position is solved for each, in order)
        #[arg(long, value_delimiter = ',', default_value = "1")]
        n: Vec<u8>,
        /// Give up on a length after this many seconds per position and radius (report the
        /// bound proven)
        #[arg(long, default_value_t = 600.0)]
        seconds: f64,
        /// Longest solution searched (the block rule's 25 is always an upper bound)
        #[arg(long, default_value_t = 25)]
        max: usize,
        /// Wait before loading the tables while less memory than this is available, GB (and
        /// before each position while less than this minus what the tables took)
        #[arg(long, default_value_t = 30.0)]
        min_available_gb: f64,
        /// Append the results here (tab-separated: class, x, p, N, lower, upper, seconds,
        /// nodes, solution, relaxed phase-1 distance after each move)
        #[arg(long)]
        out: PathBuf,
    },
    /// Check an independent list of geodesic states and their neighbours' relaxed phase-1
    /// distances (the phase-1 verifier's `geodesic.txt`) against the relaxed table
    VerifyGeodesic {
        /// The list
        #[arg(long)]
        file: PathBuf,
        /// Directory of the relaxed phase-1 table
        #[arg(long)]
        tables: PathBuf,
    },
    /// The positions of a layer (the last by default) up to symmetry, with solutions
    Antipodes {
        /// Search directory
        #[arg(long)]
        dir: PathBuf,
        /// Phase-2 table directory (to solve)
        #[arg(long)]
        tables: Option<PathBuf>,
        /// Layer (default: the last)
        #[arg(long)]
        depth: Option<usize>,
        /// At most this many positions, the first ones (default: all)
        #[arg(long)]
        limit: Option<usize>,
        /// This many classes chosen uniformly at random (with --seed), in place of the first
        #[arg(long)]
        sample: Option<usize>,
        /// Random seed of --sample
        #[arg(long, default_value_t = 1)]
        seed: u64,
        /// List the positions without solving them: class, x, p, depth and the pattern of the
        /// position as held in the fixed frame (each checked to read back as its slot)
        #[arg(long)]
        no_solve: bool,
        /// Write the positions here (tab-separated: class, x, p, depth, pattern, solution,
        /// scramble)
        #[arg(long)]
        out: Option<PathBuf>,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let threads = cli.threads.unwrap_or_else(|| {
        std::thread::available_parallelism().map_or(1, std::num::NonZeroUsize::get)
    });
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .context("thread pool")?;
    match cli.cmd {
        Cmd::Run {
            dir,
            inst,
            ram_gb,
            max_depth,
            max_parts,
            max_minutes,
            keep,
            store,
            readers,
            quiet,
        } => {
            let i = Instance::new(inst.kind());
            let cfg = Config {
                dir: dir.clone(),
                ram: (ram_gb * 1e9) as usize,
                max_depth,
                max_parts,
                max_seconds: max_minutes.map(|m| m * 60.0),
                threads,
                keep,
                stores: store,
                readers,
                verbose: !quiet,
            };
            let mut bfs = Bfs::open(&i, cfg).context("opening the search")?;
            let out = bfs.run().context("search")?;
            status(&dir, &i, bfs.manifest());
            println!(
                "{}",
                match out {
                    Outcome::Complete => "complete",
                    Outcome::Stopped => "stopped (run again to resume)",
                }
            );
        }
        Cmd::Status { dir } => {
            let m = load(&dir)?;
            let i = Instance::new(m.kind);
            status(&dir, &i, &m);
        }
        Cmd::Reference { inst, compare } => {
            let i = Instance::new(inst.kind());
            let t = slot_bfs(&i, true);
            println!("slots per depth: {:?}", t.hist);
            if let Some(dir) = compare {
                let m = load(&dir)?;
                if m.kind != i.kind {
                    bail!("{} is a search of another instance", dir.display());
                }
                compare_counts(
                    "slots",
                    &t.hist,
                    &m.levels.iter().map(|l| l.slots).collect::<Vec<_>>(),
                    m.complete,
                )?;
            }
        }
        Cmd::Positions { depth, compare } => {
            let counts = position_bfs(depth);
            println!("positions per depth: {counts:?}");
            if let Some(dir) = compare {
                let m = load(&dir)?;
                let got: Vec<u64> = m
                    .levels
                    .iter()
                    .map(|l| l.positions)
                    .take(depth + 1)
                    .collect();
                compare_counts("positions", &counts, &got, false)?;
            }
        }
        Cmd::Tables { tables, ac, px } => {
            let t = ChenTables::load(
                &tables,
                TableSet {
                    sb: false,
                    px: true,
                },
                false,
            )
            .context("loading the phase-2 tables")?;
            println!("AC entries per depth: {:?}", t.ac.hist);
            let px_hist = t.px.as_ref().map(|p| p.hist.clone()).unwrap_or_default();
            println!("PX entries per depth: {px_hist:?}");
            for (name, hist, dir) in [("AC", &t.ac.hist, ac), ("PX", &px_hist, px)] {
                if let Some(dir) = dir {
                    let m = load(&dir)?;
                    let got: Vec<u64> = m.levels.iter().map(|l| l.slots).collect();
                    compare_counts(name, hist, &got, m.complete)?;
                }
            }
        }
        Cmd::Orbits { dir, depth } => {
            let m = load(&dir)?;
            let i = Instance::new(m.kind);
            let depths: Vec<usize> = if depth.is_empty() {
                (0..m.levels.len())
                    .filter(|&d| Layout::of(&dir, &m).has(d))
                    .collect()
            } else {
                depth
            };
            let mut bad = false;
            for d in depths {
                let mut n = 0u64;
                for_slots(&Layout::of(&dir, &m), &i, d, 0..i.n_class, |c, x, p| {
                    if canonical(&i, c, x, p) {
                        n += 1;
                    }
                })?;
                let want = m.levels.get(d).map_or(0, |l| l.orbits);
                println!("depth {d}: {n} canonical slots, manifest {want}");
                bad |= n != want;
            }
            if bad {
                bail!("MISMATCH");
            }
        }
        Cmd::Sample {
            dir,
            tables,
            depth,
            n,
            seed,
        } => {
            let m = load(&dir)?;
            let i = Instance::new(m.kind);
            let t = ChenTables::load(
                &tables,
                TableSet {
                    sb: false,
                    px: true,
                },
                false,
            )
            .context("loading the phase-2 tables")?;
            let solver = Solver::new(&t);
            // reservoir sample of the layer's slots
            let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(seed);
            let mut pick: Vec<(usize, usize, usize)> = vec![];
            let mut seen = 0usize;
            for_slots(&Layout::of(&dir, &m), &i, depth, 0..i.n_class, |c, x, p| {
                seen += 1;
                if pick.len() < n {
                    pick.push((c, x, p));
                } else {
                    let j = rng.gen_range(0..seen);
                    if j < n {
                        pick[j] = (c, x, p);
                    }
                }
            })?;
            let mut bad = 0;
            for &(c, x, p) in &pick {
                let Some(st) = PState::of_slot(&i, c, x, p) else {
                    bail!("not the full instance");
                };
                if !solver.indices_agree(&st) {
                    bail!("table indices disagree for slot ({c}, {x}, {p})");
                }
                let sol = solver.solve(&st, depth + 2);
                let len = sol.as_ref().map(Vec::len);
                if len != Some(depth) {
                    bad += 1;
                    println!("slot ({c}, {x}, {p}): IDA* distance {len:?}, layer {depth}");
                }
            }
            println!(
                "{} slots of {seen} at depth {depth}: {} with IDA* distance {depth}",
                pick.len(),
                pick.len() - bad
            );
            if bad > 0 {
                bail!("MISMATCH");
            }
        }
        Cmd::BenchPush {
            inst,
            classes,
            blocks,
            density,
            lo,
            hi,
        } => {
            let i = Instance::new(inst.kind());
            let (moved, touched, secs) =
                bench_push(&i, classes, blocks, density, (lo, hi.min(i.n_class)));
            println!(
                "{moved} slot moves into {touched} blocks in {secs:.2} s: {:.2} ns per slot move \
                 per thread ({threads} threads)",
                secs * 1e9 * threads as f64 / moved.max(1) as f64
            );
        }
        Cmd::Random {
            tables,
            n,
            seed,
            max,
        } => {
            let t = ChenTables::load(
                &tables,
                TableSet {
                    sb: false,
                    px: true,
                },
                false,
            )
            .context("loading the phase-2 tables")?;
            let solver = Solver::new(&t);
            let s = omelettehedron::chen::space::space();
            let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(seed);
            let mut hist = vec![0usize; max + 2];
            let positions: Vec<PState> = (0..n)
                .map(|_| {
                    let mut l = omelettehedron::cube::Cube::solved().xb;
                    rand::seq::SliceRandom::shuffle(&mut l[..], &mut rng);
                    let mut xb = [0u8; 16];
                    xb[..12].copy_from_slice(&l);
                    PState {
                        q: rng.gen_range(0..s.pat.cfg.len()) as u16,
                        pc: rng.gen_range(0..s.corners.next.len()) as u16,
                        xb,
                    }
                })
                .collect();
            // one position per thread, each searched on its thread
            let found: Vec<usize> = {
                use rayon::prelude::*;
                positions
                    .par_iter()
                    .enumerate()
                    .map(|(k, st)| {
                        let t0 = std::time::Instant::now();
                        let h = solver.exact(st).unwrap_or([0; 3]);
                        let d = solver.solve(st, max).map_or(max + 1, |p| p.len());
                        println!(
                            "{k}: distance {d} (tables AC {} PX {} CC {}), {:.1} s",
                            h[0],
                            h[1],
                            h[2],
                            t0.elapsed().as_secs_f64()
                        );
                        d
                    })
                    .collect()
            };
            for d in found {
                hist[d] += 1;
            }
            println!(
                "distances: {:?}",
                hist.iter()
                    .enumerate()
                    .filter(|e| *e.1 > 0)
                    .collect::<Vec<_>>()
            );
        }
        Cmd::Plan { inst, ram_gb } => {
            let i = Instance::new(inst.kind());
            let (parts, reads) = dense_plan(&i, (ram_gb * 1e9) as usize);
            // allowed moves per configuration, on average
            let moves: f64 = (0..i.n_class)
                .map(|c| {
                    let n = i.edges[c].iter().filter(|e| e.0 != NO_EDGE).count();
                    (n as u64 * i.class_size(c)) as f64
                })
                .sum::<f64>()
                / i.n_config as f64;
            println!("allowed moves per configuration: {moves:.3} on average");
            println!(
                "{} slots, {} positions, {} up to symmetry",
                i.n_slots(),
                i.n_positions(),
                i.n_orbits()
            );
            println!(
                "{} parts of up to {} classes; {reads} class files read per step = {:.2} x the \
                 {} classes",
                parts.len(),
                parts.iter().map(|p| p.1 - p.0).max().unwrap_or(0),
                reads as f64 / i.n_class as f64,
                i.n_class
            );
        }
        Cmd::Joint {
            dir,
            tables,
            checks,
            outside,
            limit,
            out,
            sample,
            solutions,
        } => joint(
            &dir,
            &tables,
            checks,
            outside,
            limit,
            out.as_deref(),
            sample,
            solutions.as_deref(),
        )?,
        Cmd::Paths { dir, tables } => paths(&dir, &tables)?,
        Cmd::Radius { tables, input } => radius(&tables, &input)?,
        Cmd::VerifyGeodesic { file, tables } => {
            use omelettehedron::chen::p2bfs::joint::{ae_depths, verify_geodesic};
            use omelettehedron::chen::relaxed::RelaxedAe;
            use omelettehedron::tables::AeSym;
            // the list's positions (its legend): A-centres, then the edges' A stickers
            const A: [[i32; 3]; 12] = [
                [2, 2, 5],
                [2, 5, 2],
                [5, 2, 2],
                [2, -5, -2],
                [2, -2, -5],
                [5, -2, -2],
                [-5, 2, -2],
                [-2, 2, -5],
                [-2, 5, -2],
                [-5, -2, 2],
                [-2, -5, 2],
                [-2, -2, 5],
            ];
            const E: [[i32; 3]; 12] = [
                [1, 4, 4],
                [4, 1, 4],
                [4, 4, 1],
                [1, -4, -4],
                [4, -4, -1],
                [4, -1, -4],
                [-4, 1, -4],
                [-4, 4, -1],
                [-1, 4, -4],
                [-4, -4, 1],
                [-4, -1, 4],
                [-1, -4, 4],
            ];
            let ix = AeSym::new(true);
            let relaxed = RelaxedAe::load(&tables, &ix, false).context("relaxed table")?;
            let dep = ae_depths(&ix, &relaxed.targets, false);
            let text =
                std::fs::read_to_string(&file).with_context(|| file.display().to_string())?;
            let r = verify_geodesic(&text, &ix, &dep, &relaxed.tab, &A, &E)
                .map_err(|e| anyhow::anyhow!(e))?;
            println!(
                "{} states, {} agree; {} neighbours, {} agree; {} distinct neighbours outside \
                 the list; {} entries whose table value is not the distance mod 3",
                r.states, r.states_ok, r.neighbours, r.neighbours_ok, r.outside, r.table_mismatch
            );
            for m in &r.mismatches {
                println!("MISMATCH {m}");
            }
            if !r.mismatches.is_empty() || r.table_mismatch > 0 {
                bail!("MISMATCH");
            }
        }
        Cmd::RadiusSolve {
            input,
            tables,
            opt_tables,
            n,
            seconds,
            max,
            min_available_gb,
            out,
        } => radius_solve(
            &input,
            &tables,
            &opt_tables,
            &n,
            seconds,
            max,
            min_available_gb,
            &out,
        )?,
        Cmd::Bad {
            dir,
            tables,
            top,
            geodesic,
        } => {
            use omelettehedron::chen::p2bfs::joint::{
                ae_depths, bad_sets, geodesic as geo_set, Phase2Slots,
            };
            use omelettehedron::chen::relaxed::RelaxedAe;
            use omelettehedron::tables::AeSym;
            let m = load(&dir)?;
            let last = m.depth().context("no layer")?;
            let ix = AeSym::new(true);
            let relaxed = RelaxedAe::load(&tables, &ix, false).context("relaxed table")?;
            let dep = ae_depths(&ix, &relaxed.targets, false);
            let g = geodesic.then(|| geo_set(&ix, &dep));
            let top = top.unwrap_or((dep.hist.len() - 1) as u8);
            let layer = Phase2Slots::load(&dir, last)?;
            let sizes = bad_sets(&ix, &dep, &layer, top, g.as_ref(), true);
            println!(
                "bad sets above layer {last} (to distance {top}, geodesic {geodesic}): {sizes:?}"
            );
        }
        Cmd::Antipodes {
            dir,
            tables,
            depth,
            limit,
            sample,
            seed,
            no_solve,
            out,
        } => antipodes(
            &dir,
            tables.as_deref(),
            depth,
            limit,
            sample,
            seed,
            no_solve,
            out.as_deref(),
        )?,
    }
    Ok(())
}

/// `antipodes`: a layer's positions up to symmetry (all, the first, or a random sample), with
/// shortest phase-2 solutions or just their patterns.
#[allow(clippy::too_many_arguments)]
fn antipodes(
    dir: &Path,
    tables: Option<&Path>,
    depth: Option<usize>,
    limit: Option<usize>,
    sample: Option<usize>,
    seed: u64,
    no_solve: bool,
    out: Option<&Path>,
) -> Result<()> {
    use omelettehedron::chen::p2bfs::engine::canonical_slot;
    use omelettehedron::chen::p2bfs::joint::Phase2Slots;
    let m = load(dir)?;
    let i = Instance::new(m.kind);
    let d = depth.or(m.depth()).context("no layer")?;
    let mut slots = vec![];
    for_slots(&Layout::of(dir, &m), &i, d, 0..i.n_class, |c, x, p| {
        if canonical(&i, c, x, p) && (sample.is_some() || limit.is_none_or(|l| slots.len() < l)) {
            slots.push((c, x, p));
        }
    })?;
    let total = slots.len();
    if let Some(k) = sample {
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(seed);
        rand::seq::SliceRandom::shuffle(&mut slots[..], &mut rng);
        slots.truncate(k);
        slots.sort_unstable();
    }
    let mut lines = vec![];
    if no_solve {
        let layer = Phase2Slots::load(dir, d)?;
        for &(c, x, p) in &slots {
            let slot = (c as u32, x as u32, p as u32);
            let cube = layer.cube(slot).context("slot without a position")?;
            let pattern = omelettehedron::pattern::format(&omelettehedron::pattern::of_cube(&cube));
            // the pattern read back: the same slot (up to the class's stabiliser)
            let parsed =
                omelettehedron::pattern::parse(&pattern).map_err(|e| anyhow::anyhow!("{e}"))?;
            let back =
                omelettehedron::pattern::position(&parsed).map_err(|e| anyhow::anyhow!("{e}"))?;
            let p2 = omelettehedron::chen::tables::P2::of(&back)
                .context("the pattern is no phase-2 position")?;
            let (c2, x2, q2) = layer.slot_of(&p2);
            let (x2, q2) = canonical_slot(&i, c2 as usize, x2 as usize, q2 as usize);
            if (c2 as usize, x2, q2) != (c, x, p) || layer.contains(&back) != Some(true) {
                bail!("slot ({c}, {x}, {p}): the pattern reads back as ({c2}, {x2}, {q2})");
            }
            let line = format!("{c}\t{x}\t{p}\t{d}\t{pattern}");
            println!("{line}");
            lines.push(line);
        }
        if let Some(out) = out {
            let text = format!("# class\tx\tp\tdepth\tpattern\n{}\n", lines.join("\n"));
            std::fs::write(out, text).with_context(|| out.display().to_string())?;
        }
        println!(
            "{} of the {total} classes at depth {d}; every pattern reads back as its slot",
            slots.len()
        );
        return Ok(());
    }
    let tables = tables.context("--tables is needed to solve (or --no-solve)")?;
    let t = ChenTables::load(
        tables,
        TableSet {
            sb: false,
            px: true,
        },
        false,
    )
    .context("loading the phase-2 tables")?;
    let solver = Solver::new(&t);
    for &(c, x, p) in &slots {
        let Some(st) = PState::of_slot(&i, c, x, p) else {
            bail!("not the full instance");
        };
        let sol = solver
            .solve_par(&st, d + 2)
            .with_context(|| format!("no solution for slot ({c}, {x}, {p})"))?;
        if sol.len() != d {
            bail!(
                "slot ({c}, {x}, {p}): solution of {} moves, layer {d}",
                sol.len()
            );
        }
        let w = write(&st, &sol);
        if !w.consistent {
            bail!("slot ({c}, {x}, {p}): the scramble's pattern is not the position's");
        }
        let line = format!(
            "{c}\t{x}\t{p}\t{d}\t{}\t{}\t{}",
            w.pattern,
            format_moves(&w.solution),
            format_moves(&w.scramble)
        );
        println!("{line}");
        lines.push(line);
    }
    if let Some(out) = out {
        let text = format!(
            "# class\tx\tp\tdepth\tpattern\tsolution\tscramble\n{}\n",
            lines.join("\n")
        );
        std::fs::write(out, text).with_context(|| out.display().to_string())?;
    }
    println!("{} positions at depth {d} (up to symmetry)", slots.len());
    Ok(())
}

fn load(dir: &Path) -> Result<Manifest> {
    Manifest::load(dir)?.with_context(|| format!("{}: no search here", dir.display()))
}

/// Compare two count lists; a search not complete is compared up to its last depth.
fn compare_counts(what: &str, want: &[u64], got: &[u64], complete: bool) -> Result<()> {
    let n = if complete {
        want.len().max(got.len())
    } else {
        got.len()
    };
    let ok = (0..n).all(|d| want.get(d) == got.get(d));
    if ok {
        println!(
            "{what}: the search agrees at depths 0..{}",
            n.saturating_sub(1)
        );
        Ok(())
    } else {
        println!("{what}: reference {want:?}");
        println!("{what}: search    {got:?}");
        bail!("MISMATCH")
    }
}

fn status(dir: &Path, i: &Instance, m: &Manifest) {
    println!("search {} ({:?})", dir.display(), m.kind);
    println!("depth\tslots\tpositions\tup to symmetry\tGB\tseconds\tGB read");
    let (mut slots, mut pos, mut orb) = (0u64, 0u128, 0u64);
    for l in &m.levels {
        println!(
            "{}\t{}\t{}\t{}\t{:.3}\t{:.1}\t{:.3}",
            l.depth,
            l.slots,
            l.positions,
            l.orbits,
            l.bytes as f64 / 1e9,
            l.seconds,
            l.read as f64 / 1e9
        );
        slots += l.slots;
        pos += u128::from(l.positions);
        orb += l.orbits;
    }
    println!(
        "total\t{slots}\t{pos}\t{orb}\t(of {} slots, {} positions){}",
        i.n_slots(),
        i.n_positions(),
        if m.complete { ", complete" } else { "" }
    );
    if m.complete {
        println!(
            "every position reached: {} (positions {pos} of {}; up to symmetry {orb} of {})",
            if pos == i.n_positions() { "yes" } else { "NO" },
            i.n_positions(),
            i.n_orbits()
        );
        let hp_max = (0..i.n_class)
            .filter(|&c| i.hprime[c])
            .filter_map(|c| {
                m.levels
                    .iter()
                    .rev()
                    .find(|l| l.per_class.get(c).copied().unwrap_or(0) > 0)
                    .map(|l| l.depth)
            })
            .max();
        println!(
            "deepest layer: {}; deepest in H': {}",
            m.levels.last().map_or(0, |l| l.depth),
            hp_max.map_or("-".into(), |d| d.to_string())
        );
    }
    // the deepest layer of each class of H'
    let deepest = |c: usize| {
        m.levels
            .iter()
            .rev()
            .find(|l| l.per_class.get(c).copied().unwrap_or(0) > 0)
            .map(|l| l.depth)
    };
    let hp: Vec<(usize, Option<usize>)> = (0..i.n_class)
        .filter(|&c| i.hprime[c])
        .map(|c| (c, deepest(c)))
        .collect();
    println!(
        "H' classes (class: deepest layer so far): {}",
        hp.iter()
            .map(|(c, d)| format!("{c}: {}", d.map_or("-".into(), |d| d.to_string())))
            .collect::<Vec<_>>()
            .join(", ")
    );
    if let Some(p) = &m.progress {
        println!(
            "making depth {}: {} of {} parts done",
            p.depth,
            p.done.iter().filter(|d| d.is_some()).count(),
            p.parts.len()
        );
    }
}

/// `joint`: the joint check, with its checks first.
#[allow(clippy::too_many_arguments)]
fn joint(
    dir: &Path,
    tables: &Path,
    checks: u64,
    outside: usize,
    limit: Option<usize>,
    out: Option<&Path>,
    sample: usize,
    solutions: Option<&Path>,
) -> Result<()> {
    use omelettehedron::chen::p2bfs::joint::{
        ae_depths, check_d1, check_geodesic, geodesic, state_counts, Joint, Phase2Slots,
    };
    use omelettehedron::chen::relaxed::RelaxedAe;
    use omelettehedron::tables::{AeIndex, AeSym};
    use std::time::Instant;
    let t0 = Instant::now();
    let m = load(dir)?;
    if !m.complete {
        bail!("{}: the search is not complete", dir.display());
    }
    let last = m.depth().context("no layer")?;
    let ix = AeSym::new(true);
    let relaxed = RelaxedAe::load(tables, &ix, false).context("loading the relaxed table")?;
    println!(
        "relaxed table: {} entries, {} target entries, entries per distance {:?}",
        ix.size(),
        relaxed.targets.len(),
        relaxed.tab.hist
    );
    // exact distances, checked against the table entry by entry (mod 3) and per distance
    let dep = ae_depths(&ix, &relaxed.targets, true);
    println!(
        "exact distances: entries per distance {:?} ({:.0} s)",
        dep.hist,
        t0.elapsed().as_secs_f64()
    );
    if dep.hist != relaxed.tab.hist {
        bail!("MISMATCH: the exact distances' histogram is not the table's");
    }
    let bad = {
        use rayon::prelude::*;
        (0..ix.size())
            .into_par_iter()
            .filter(|&i| dep.get(i) % 3 != relaxed.tab.get(i))
            .count()
    };
    if bad > 0 {
        bail!("MISMATCH: {bad} entries differ from the table mod 3");
    }
    println!("every entry agrees with the table (mod 3)");
    let n = check_d1(&ix, &relaxed, &dep, checks).map_err(|e| anyhow::anyhow!(e))?;
    println!(
        "d1 of {n} random positions: the table's walk, the exact array and all 24 conjugates agree"
    );
    let geo = geodesic(&ix, &dep);
    println!(
        "geodesic set (entries on a shortest path to distance {}), per distance: {:?} ({:.0} s)",
        geo.top,
        geo.per_level,
        t0.elapsed().as_secs_f64()
    );
    // the layer, and membership checks
    let layer = Phase2Slots::load(dir, last)?;
    println!(
        "layer {last}: {} slots, {} classes",
        layer.len(),
        layer.reps.len()
    );
    let mut missing = 0;
    let s48 = omelettehedron::optimal::sym::sym48();
    for (k, &slot) in layer.reps.iter().enumerate() {
        let c = layer.cube(slot).context("slot without a position")?;
        if layer.contains(&c) != Some(true) {
            missing += 1;
        }
        if k % 97 == 0 {
            for sy in 0..24 {
                let (cc, _) = s48.conj_cube(&c, sy).context("conjugate")?;
                if layer.contains(&cc) != Some(true) {
                    missing += 1;
                }
            }
        }
    }
    if missing > 0 {
        bail!("MISMATCH: {missing} positions of the layer test outside it");
    }
    println!(
        "all {} classes' positions test in the layer, and the conjugates of every 97th",
        layer.reps.len()
    );
    let ct = ChenTables::load(
        tables,
        TableSet {
            sb: false,
            px: true,
        },
        false,
    )
    .context("loading the phase-2 tables")?;
    let solver = Solver::new(&ct);
    let s = omelettehedron::chen::space::space();
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(0x5025);
    let mut lens = vec![];
    for _ in 0..outside {
        let mut l = omelettehedron::cube::Cube::solved().xb;
        rand::seq::SliceRandom::shuffle(&mut l[..], &mut rng);
        let mut xb = [0u8; 16];
        xb[..12].copy_from_slice(&l);
        let st = PState {
            q: rng.gen_range(0..s.pat.cfg.len()) as u16,
            pc: rng.gen_range(0..s.corners.next.len()) as u16,
            xb,
        };
        let sol = solver
            .solve_par(&st, last - 1)
            .context("a random position beyond 24")?;
        let p2 = solver.views().to_p2(&st);
        let (xa, ep) = s.cfgs[p2.cfg as usize];
        let (cp, co) = omelettehedron::coord::corners().id_to_state[p2.c as usize];
        let mut xbf = [0u8; 12];
        xbf.copy_from_slice(&p2.xb[..12]);
        let c = omelettehedron::cube::Cube {
            cp,
            co,
            ep,
            xa,
            xb: xbf,
        };
        if layer.contains(&c) != Some(false) {
            bail!(
                "MISMATCH: a position at distance {} tests in the layer",
                sol.len()
            );
        }
        lens.push(sol.len());
    }
    println!(
        "{outside} random positions solved within {} moves ({lens:?}) test outside",
        last - 1
    );
    // the check
    let mut layer = layer;
    if let Some(l) = limit {
        layer.reps.truncate(l);
    }
    let j = Joint {
        ix: &ix,
        dep: &dep,
        geo: &geo,
        layer: &layer,
    };
    let t1 = Instant::now();
    let r = j.run(true);
    println!(
        "joint: {} starts ({} on the geodesic set), {} candidates, {} resolved by another \
         optimal phase-1 solution, {} unresolved, {} anomalies; {} downward nodes; {:.0} s",
        r.starts,
        r.starts_on_geodesics,
        r.candidates,
        r.resolved,
        r.unresolved.len(),
        r.anomalies,
        r.down_nodes,
        t1.elapsed().as_secs_f64()
    );
    let lines: Vec<String> = r
        .unresolved
        .iter()
        .map(|c| {
            format!(
                "{}\t{:?}",
                omelettehedron::pattern::format(&omelettehedron::pattern::of_cube(c)),
                c
            )
        })
        .collect();
    for l in lines.iter().take(20) {
        println!("unresolved: {l}");
    }
    if let Some(out) = out {
        std::fs::write(out, lines.join("\n") + "\n").with_context(|| out.display().to_string())?;
    }
    if r.anomalies > 0 {
        bail!("{} anomalies", r.anomalies);
    }
    // the geodesic set is exactly the backward closure
    let bad = check_geodesic(&ix, &dep, &geo);
    if bad > 0 {
        bail!("MISMATCH: {bad} entries break the geodesic set's rule");
    }
    println!("the geodesic set is exactly the backward closure (all entries checked)");
    let (all, on) = state_counts(&ix, &dep, &geo);
    println!("A-centre/edge states per distance: {all:?}");
    println!("of them on the geodesic set: {on:?}");
    // explicit solutions of sampled candidates: an alternative optimal phase-1 solution, then
    // a phase-2 solution by IDA* within 24 moves (independent of the layer files)
    let f = omelettehedron::fast::fast();
    let step = (layer.reps.len() / sample.max(1)).max(1);
    let mut lines = vec![];
    let mut longest = 0;
    let mut ends = (0u64, 0u64);
    for &slot in layer.reps.iter().step_by(step).take(sample) {
        let sc = layer.cube(slot).context("slot without a position")?;
        for q in j.candidates_from(&sc, 2) {
            let d1 = relaxed.exact(&ix, &f.of(&q)).context("walk")?;
            let (inside, outside) = j.count_ends(&q);
            if inside == 0 || outside == 0 {
                bail!("MISMATCH: {inside} optimal solutions end inside, {outside} outside");
            }
            ends.0 += inside;
            ends.1 += outside;
            let (w, e) = j
                .alternative(&q)
                .context("a candidate without an alternative")?;
            let mut c = q;
            for &m in &w {
                c = c.apply(m);
            }
            let p2 = omelettehedron::chen::tables::P2::of(&e).context("no phase-2 position")?;
            if c != e || d1 != 11 || w.len() != 11 {
                bail!("MISMATCH: d1 {d1}, path of {} moves", w.len());
            }
            let st = solver.views().of_p2(&p2);
            let sol2 = solver
                .solve_par(&st, last - 1)
                .context("the end needs more than 24 moves")?;
            let sol2f = solver
                .views()
                .to_fixed(&st, &sol2)
                .context("a phase-2 move without a fixed-frame move")?;
            let all: Vec<usize> = w.iter().chain(&sol2f).copied().collect();
            let mut z = q;
            for &m in &all {
                z = z.apply(m);
            }
            if z != omelettehedron::cube::Cube::solved() {
                bail!("MISMATCH: the joined solution does not solve the candidate");
            }
            let phys = omelettehedron::notation::solution_to_physical(
                omelettehedron::notation::Frame::identity(),
                &all,
                false,
            );
            longest = longest.max(phys.len());
            lines.push(format!(
                "{}\t{}\t{}\t{}\t{}",
                omelettehedron::pattern::format(&omelettehedron::pattern::of_cube(&q)),
                w.len(),
                sol2.len(),
                phys.len(),
                format_moves(&phys)
            ));
        }
    }
    println!(
        "{} sampled candidates solved explicitly: 11 + at most {} = at most {longest} moves; \
         their optimal phase-1 solutions end {} times in the layer and {} times outside",
        lines.len(),
        longest - 11,
        ends.0,
        ends.1
    );
    if let Some(p) = solutions {
        std::fs::write(p, lines.join("\n") + "\n").with_context(|| p.display().to_string())?;
    }
    Ok(())
}

/// `paths`: ascending phase-1 path counts per configuration class, weighted by the layers'
/// positions per class.
fn paths(dir: &Path, tables: &Path) -> Result<()> {
    use omelettehedron::chen::p2bfs::joint::{
        ae_depths, geodesic, path_counts, Joint, Phase2Slots,
    };
    use omelettehedron::chen::relaxed::RelaxedAe;
    use omelettehedron::coord::Coords;
    use omelettehedron::tables::{AeIndex, AeSym};
    use std::time::Instant;
    let t0 = Instant::now();
    let m = load(dir)?;
    let last = m.depth().context("no layer")?;
    let i = Instance::new(Kind::FULL);
    let ix = AeSym::new(true);
    let relaxed = RelaxedAe::load(tables, &ix, false).context("loading the relaxed table")?;
    let dep = ae_depths(&ix, &relaxed.targets, false);
    let top = (dep.hist.len() - 1) as u8;
    let s = omelettehedron::chen::space::space();
    // the A-centre/edge entry of each configuration class (its representative, frame 0)
    let entry: Vec<usize> = (0..i.n_class)
        .map(|c| {
            let (xa, ep) = s.cfgs[s.pat.cfg[i.rep[c] as usize] as usize];
            ix.index(&Coords::of(&omelettehedron::cube::Cube {
                xa,
                ep,
                ..omelettehedron::cube::Cube::solved()
            }))
        })
        .collect();
    let pos = |d: usize| -> Vec<f64> {
        m.levels.get(d).map_or(vec![0.0; i.n_class], |l| {
            (0..i.n_class)
                .map(|c| (i.class_size(c) * l.per_class.get(c).copied().unwrap_or(0)) as f64)
                .collect()
        })
    };
    let (p_last, p_prev) = (pos(last), pos(last - 1));
    let all_cfg: f64 = (0..i.n_class).map(|c| i.class_size(c) as f64).sum();
    for l in [top, top - 1] {
        let val = path_counts(&ix, &dep, l);
        let per: Vec<f64> = entry.iter().map(|&e| f64::from(val[e])).collect();
        let uniform: f64 = (0..i.n_class)
            .map(|c| i.class_size(c) as f64 * per[c])
            .sum::<f64>()
            / all_cfg;
        let maxp = per.iter().copied().fold(0.0, f64::max);
        let from = |p: &[f64]| -> (f64, f64) {
            let n: f64 = p.iter().sum();
            let c: f64 = p.iter().zip(&per).map(|(a, b)| a * b).sum();
            (c, c / n.max(1.0))
        };
        let (c_last, a_last) = from(&p_last);
        let (c_prev, a_prev) = from(&p_prev);
        println!(
            "paths of {l} moves up to distance {l}: per configuration at most {maxp:.4e}, \
             {uniform:.4e} on average over all phase-2 positions; from layer {last}: \
             {c_last:.6e} candidates ({a_last:.4e} per position); from layer {}: {c_prev:.6e} \
             ({a_prev:.4e} per position) ({:.0} s)",
            last - 1,
            t0.elapsed().as_secs_f64()
        );
        if l == top {
            // the same from the layer's slots, and per class representative (the joint
            // check's enumeration)
            let layer = Phase2Slots::load(dir, last)?;
            let mut by_slot = 0f64;
            let mut by_rep = 0f64;
            let reps: std::collections::HashSet<_> = layer.reps.iter().copied().collect();
            omelettehedron::chen::p2bfs::engine::for_slots(
                &Layout::of(dir, &m),
                &i,
                last,
                0..i.n_class,
                |c, x, p| {
                    if let Some(cube) = layer.cube((c as u32, x as u32, p as u32)) {
                        let v = f64::from(val[ix.index(&Coords::of(&cube))]);
                        by_slot += i.class_size(c) as f64 * v;
                        if reps.contains(&(c as u32, x as u32, p as u32)) {
                            by_rep += v;
                        }
                    }
                },
            )?;
            println!(
                "check: from the layer's slots {by_slot:.6e} candidates; per class \
                 representative {by_rep:.0} (the joint check enumerated 1,534,416)"
            );
            // the measured rate of the joint check (enumeration and resolution)
            let geo = geodesic(&ix, &dep);
            let j = Joint {
                ix: &ix,
                dep: &dep,
                geo: &geo,
                layer: &layer,
            };
            let t1 = Instant::now();
            let r = j.run(false);
            let secs = t1.elapsed().as_secs_f64();
            println!(
                "joint check: {} candidates in {secs:.3} s = {:.3e} per second ({} threads)",
                r.candidates,
                r.candidates as f64 / secs,
                rayon::current_num_threads()
            );
        }
    }
    Ok(())
}

/// `radius`: the relaxed phase-1 distance after each move of each solution.
fn radius(tables: &Path, input: &Path) -> Result<()> {
    use omelettehedron::chen::relaxed::RelaxedAe;
    use omelettehedron::notation::{parse_moves, Frame};
    use omelettehedron::tables::AeSym;
    let ix = AeSym::new(true);
    let relaxed = RelaxedAe::load(tables, &ix, false).context("loading the relaxed table")?;
    let f = omelettehedron::fast::fast();
    let text = std::fs::read_to_string(input).with_context(|| input.display().to_string())?;
    let mut hist = std::collections::BTreeMap::new();
    for line in text
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
    {
        let c: Vec<&str> = line.split('\t').collect();
        let (Some(pat), Some(sol)) = (c.get(4), c.get(8)) else {
            continue;
        };
        if *sol == "-" {
            continue;
        }
        let parsed = omelettehedron::pattern::parse(pat).map_err(|e| anyhow::anyhow!("{e}"))?;
        let mut cube =
            omelettehedron::pattern::position(&parsed).map_err(|e| anyhow::anyhow!("{e}"))?;
        let moves = parse_moves(sol).map_err(|e| anyhow::anyhow!("{e}"))?;
        let d1 = |x: &omelettehedron::cube::Cube| relaxed.exact(&ix, &f.of(x)).unwrap_or(99);
        let mut fr = Frame::identity();
        let mut profile = vec![d1(&cube)];
        for pm in moves.iter().copied() {
            let (m, nf) = fr.to_internal(pm);
            cube = cube.apply(m);
            fr = nf;
            profile.push(d1(&cube));
        }
        if cube != omelettehedron::cube::Cube::solved() {
            bail!(
                "{}: the solution does not solve the position",
                c.first().unwrap_or(&"?")
            );
        }
        let top = profile.iter().copied().max().unwrap_or(0);
        *hist.entry(top).or_insert(0) += 1;
        println!(
            "{}\t{}\t{}\t{}\tmax d1 {top}\t{}",
            c.first().unwrap_or(&""),
            c.get(1).unwrap_or(&""),
            c.get(2).unwrap_or(&""),
            moves.len(),
            profile
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(" ")
        );
    }
    println!("solutions by the largest relaxed phase-1 distance along them: {hist:?}");
    Ok(())
}

/// Memory available (GB), from `/proc/meminfo` (`None` where there is none).
fn available_gb() -> Option<f64> {
    let text = std::fs::read_to_string("/proc/meminfo").ok()?;
    let line = text.lines().find(|l| l.starts_with("MemAvailable:"))?;
    let kb: f64 = line.split_whitespace().nth(1)?.parse().ok()?;
    Some(kb / 1e6)
}

/// `radius-solve`: radius-N distances of the positions of a file.
#[allow(clippy::too_many_arguments)]
fn radius_solve(
    input: &Path,
    tables: &Path,
    opt_tables: &Path,
    radii: &[u8],
    seconds: f64,
    max: usize,
    min_available_gb: f64,
    out: &Path,
) -> Result<()> {
    use omelettehedron::chen::p2bfs::joint::ae_depths;
    use omelettehedron::chen::p2bfs::radius::{NearSet, Radius};
    use omelettehedron::chen::relaxed::RelaxedAe;
    use omelettehedron::cube::MoveSet;
    use omelettehedron::notation::{solution_to_physical, Frame};
    use omelettehedron::optimal::coords::{opt_moves, Pos};
    use omelettehedron::optimal::tables::{OptTables, OptTier};
    use omelettehedron::tables::AeSym;
    use std::io::Write;
    let ix = AeSym::new(true);
    let near = {
        let relaxed = RelaxedAe::load(tables, &ix, false).context("loading the relaxed table")?;
        let dep = ae_depths(&ix, &relaxed.targets, false);
        NearSet::new(&dep, radii.iter().copied().max().unwrap_or(0))
    };
    eprintln!(
        "{} entries within relaxed phase-1 distance {}",
        near.len(),
        near.n
    );
    // memory: the threshold before the tables are loaded, then the same less what they took
    let wait_for = |gb: f64| {
        while available_gb().is_some_and(|g| g < gb) {
            eprintln!("waiting: less than {gb:.1} GB available");
            std::thread::sleep(std::time::Duration::from_secs(60));
        }
    };
    wait_for(min_available_gb);
    let before = available_gb().unwrap_or(0.0);
    let ot = OptTables::load(opt_tables, OptTier::Large, MoveSet::FACE, false)
        .context("loading the optimal solver's large tier")?;
    let taken = (before - available_gb().unwrap_or(before)).max(0.0);
    let per_position = (min_available_gb - taken).max(4.0);
    let b = ot.large().context("no large-tier bounds")?;
    let text = std::fs::read_to_string(input).with_context(|| input.display().to_string())?;
    // (position, radius) pairs already in the output: skipped, so that a stopped run resumes
    let done: std::collections::HashSet<String> = std::fs::read_to_string(out)
        .unwrap_or_default()
        .lines()
        .map(|l| l.split('\t').take(4).collect::<Vec<_>>().join("\t"))
        .collect();
    let mut f = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(out)
        .with_context(|| out.display().to_string())?;
    for &n in radii {
        let r = Radius {
            b: &b,
            ix: &ix,
            near: &near,
            n,
        };
        for line in text
            .lines()
            .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        {
            let c: Vec<&str> = line.split('\t').collect();
            let (Some(pat), Some(known)) = (c.get(4), c.get(6)) else {
                continue;
            };
            let from: usize = known.parse().unwrap_or(0);
            let key = format!(
                "{}\t{}\t{}\t{n}",
                c.first().unwrap_or(&""),
                c.get(1).unwrap_or(&""),
                c.get(2).unwrap_or(&"")
            );
            if done.contains(&key) {
                continue;
            }
            wait_for(per_position);
            let parsed = omelettehedron::pattern::parse(pat).map_err(|e| anyhow::anyhow!("{e}"))?;
            let cube =
                omelettehedron::pattern::position(&parsed).map_err(|e| anyhow::anyhow!("{e}"))?;
            let t0 = std::time::Instant::now();
            let res = r
                .solve(&cube, from, max, seconds)
                .context("the position is outside the radius or a table does not walk")?;
            let secs = t0.elapsed().as_secs_f64();
            let (upper, sol, profile) = match &res.solution {
                Some(path) => {
                    // the relaxed phase-1 distance along it, and that it solves
                    let mt = opt_moves();
                    let mut p = Pos::of(&cube);
                    let mut prof = vec![near.d1(&ix, &p).unwrap_or(99)];
                    for &m in path {
                        p = mt.apply(&p, m);
                        prof.push(near.d1(&ix, &p).unwrap_or(99));
                    }
                    if p != Pos::of(&omelettehedron::cube::Cube::solved())
                        || prof.iter().any(|&d| d > n)
                    {
                        bail!("a solution that does not solve or leaves the radius");
                    }
                    let phys = solution_to_physical(Frame::identity(), path, false);
                    (
                        path.len().to_string(),
                        format_moves(&phys),
                        prof.iter()
                            .map(ToString::to_string)
                            .collect::<Vec<_>>()
                            .join(" "),
                    )
                }
                None => ("25".to_string(), "-".to_string(), "-".to_string()),
            };
            let out_line = format!(
                "{}\t{}\t{}\t{n}\t{}\t{upper}\t{secs:.1}\t{}\t{sol}\t{profile}",
                c.first().unwrap_or(&""),
                c.get(1).unwrap_or(&""),
                c.get(2).unwrap_or(&""),
                res.lower,
                res.nodes
            );
            println!("{out_line}");
            writeln!(f, "{out_line}").context("writing the results")?;
        }
    }
    Ok(())
}
