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
    /// The positions of a layer (the last by default) up to symmetry, with solutions
    Antipodes {
        /// Search directory
        #[arg(long)]
        dir: PathBuf,
        /// Phase-2 table directory
        #[arg(long)]
        tables: PathBuf,
        /// Layer (default: the last)
        #[arg(long)]
        depth: Option<usize>,
        /// At most this many positions (default: all)
        #[arg(long)]
        limit: Option<usize>,
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
        Cmd::Antipodes {
            dir,
            tables,
            depth,
            limit,
            out,
        } => {
            let m = load(&dir)?;
            let i = Instance::new(m.kind);
            let d = depth.or(m.depth()).context("no layer")?;
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
            let mut slots = vec![];
            for_slots(&Layout::of(&dir, &m), &i, d, 0..i.n_class, |c, x, p| {
                if canonical(&i, c, x, p) && limit.is_none_or(|l| slots.len() < l) {
                    slots.push((c, x, p));
                }
            })?;
            let mut lines = vec![];
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
                std::fs::write(&out, text).with_context(|| out.display().to_string())?;
            }
            println!("{} positions at depth {d} (up to symmetry)", slots.len());
        }
    }
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
