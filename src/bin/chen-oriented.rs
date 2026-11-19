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
    clippy::cast_sign_loss,
    clippy::cast_lossless,
    clippy::many_single_char_names,
    clippy::similar_names,
    clippy::too_many_lines
)]
//! `chen-oriented`: exact phase-1 distances to the oriented targets of the block split
//! ([`omelettehedron::chen::oriented`]: every class-A hexagon solid and unturned, alone or with
//! the configurations one class-A turn away), by a breadth-first search on disk
//! ([`omelettehedron::chen::p1disk`]), and the checks that go with it.

use anyhow::{bail, Context, Result};
use clap::{Parser, Subcommand, ValueEnum};
use omelettehedron::chen::oriented::{P1Space, Target, N_EPERM};
use omelettehedron::chen::p1check::{self, state_key, State};
use omelettehedron::chen::p1disk::{self, Config, Nibbles};
use omelettehedron::cube::Cube;
use omelettehedron::notation::{format_moves, solution_to_physical, Frame};
use rayon::prelude::*;
use std::collections::BTreeMap;
use std::io::Write;
use std::path::PathBuf;

/// Positions per A-centre/edge state: corners up to a whole turn x B-centre colourings.
const PER_STATE: u128 = 960 * 369_600;

#[derive(Parser)]
#[command(about = "Exact phase-1 distances to the oriented targets of the block split")]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
    /// Worker threads (default: all)
    #[arg(short = 'j', long, global = true)]
    threads: Option<usize>,
}

#[derive(Clone, Copy, ValueEnum)]
enum Tgt {
    Oriented,
    T1,
}

impl Tgt {
    fn target(self) -> Target {
        match self {
            Tgt::Oriented => Target::Oriented,
            Tgt::T1 => Target::T1,
        }
    }
}

#[derive(Subcommand)]
enum Cmd {
    /// Sizes, targets and the source rows read per level for a few cache sizes
    Info,
    /// Start or resume a search
    Run {
        #[arg(long)]
        dir: PathBuf,
        #[arg(long, value_enum)]
        target: Tgt,
        /// Memory for the cache of source rows, GB
        #[arg(long, default_value_t = 8.0)]
        cache_gb: f64,
        /// Keep this much memory available (the cache shrinks below it), GB
        #[arg(long, default_value_t = 16.0)]
        floor_gb: f64,
        /// Stop after this level
        #[arg(long)]
        max_level: Option<usize>,
    },
    /// The levels so far: slots, states, positions
    Status {
        #[arg(long)]
        dir: PathBuf,
    },
    /// Build the 4-bit table `D.nib` from the level planes
    Table {
        #[arg(long)]
        dir: PathBuf,
    },
    /// Local consistency of `D.nib` over every slot (targets at 0, every other slot's smallest
    /// neighbour one less)
    Check {
        #[arg(long)]
        dir: PathBuf,
        #[arg(long, value_enum)]
        target: Tgt,
        /// Rows of the 4-bit table to cache (120 MB each)
        #[arg(long, default_value_t = 48)]
        cache_rows: usize,
    },
    /// The states of a level (default: the last): count, classes under the 288 pairs, a list
    /// (`<dir>/level<L>.tsv`), and with `--solve` a shortest solution of each read off the table
    /// (`<dir>/solutions<L>.tsv`, for an independent replay)
    Level {
        #[arg(long)]
        dir: PathBuf,
        #[arg(long, value_enum)]
        target: Tgt,
        #[arg(long)]
        level: Option<usize>,
        #[arg(long)]
        solve: bool,
    },
    /// Meet-in-the-middle lower bound for the states of a level, without the table: the
    /// classes within `--ball` moves of the targets (breadth-first search with cube moves) and
    /// every move sequence of up to `--depth` moves from one state of each class
    Mitm {
        #[arg(long)]
        dir: PathBuf,
        #[arg(long, value_enum)]
        target: Tgt,
        #[arg(long)]
        level: Option<usize>,
        #[arg(long, default_value_t = 8)]
        ball: usize,
        #[arg(long, default_value_t = 6)]
        depth: usize,
    },
    /// The exact strict table (two-phase index and moves) and slot-by-slot comparisons of the
    /// oriented and T1 tables with it and with each other, and of the oriented table mod 3 with
    /// the giant tier's table
    Cross {
        #[arg(long)]
        o: PathBuf,
        #[arg(long)]
        t1: PathBuf,
        #[arg(long)]
        giant: Option<PathBuf>,
    },
    /// The ball of radius R around T1 (the phase-2 configuration set of `chen-p2bfs --t1`)
    /// against the T1 table: every configuration's label is its table value, and every move
    /// out of the ball reaches a value above R
    BallCheck {
        /// The T1 search directory (with `D.nib`)
        #[arg(long)]
        dir: PathBuf,
        #[arg(long, default_value_t = 4)]
        radius: u8,
    },
    /// For every state of the oriented table's last level: a solution of that many moves whose
    /// last move is a class-A turn (a shortest path to T1 read off the T1 table, then the A-turn
    /// to the oriented target), written to `<dir>/solutions<L>-a.tsv` for an independent replay
    AEnd {
        /// The oriented search directory
        #[arg(long)]
        dir: PathBuf,
        /// The T1 search directory (with `D.nib`)
        #[arg(long)]
        t1: PathBuf,
    },
    /// The exact strict table alone (two-phase index and moves): entries and states per
    /// distance
    Strict,
    /// Table-free breadth-first search from the targets over A-centre/edge states (cube moves,
    /// a hash set of states): states per depth
    Plain {
        #[arg(long, value_enum)]
        target: Tgt,
        #[arg(long, default_value_t = 6)]
        depth: usize,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    if let Some(j) = cli.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(j)
            .build_global()
            .context("thread pool")?;
    }
    let sp = P1Space::new().context("tables")?;
    match cli.cmd {
        Cmd::Info => {
            println!(
                "rows {} (A-colouring classes), slots {} ({:.3e}), bit plane {:.1} GB",
                sp.rows(),
                sp.slots(),
                sp.slots() as f64,
                sp.rows() as f64 * omelettehedron::chen::oriented::ROW_STRIDE as f64 / 1e9
            );
            for t in [Target::Oriented, Target::T1] {
                println!(
                    "{}: {} target states, {} target slots",
                    t.name(),
                    sp.target_states(t).len(),
                    sp.target_slots(t).len()
                );
            }
            let order = p1disk::row_order(&sp);
            for cap in [40usize, 80, 160, 250, 400, 700] {
                let r = p1disk::simulate_reads(&sp, &order, cap);
                let g = p1disk::simulate_reads(&sp, &p1disk::greedy_order(&sp, cap), cap);
                println!("cache {cap} rows: source rows read per dense level: breadth-first order {r}, greedy {g}");
            }
        }
        Cmd::Run {
            dir,
            target,
            cache_gb,
            floor_gb,
            max_level,
        } => {
            let cfg = Config {
                dir,
                target: target.target(),
                cache_bytes: (cache_gb * 1e9) as u64,
                floor_bytes: (floor_gb * 1e9) as u64,
                max_level,
            };
            let lv = p1disk::run(&sp, &cfg).context("search")?;
            print_levels(&lv);
        }
        Cmd::Plain { target, depth } => plain(&sp, target.target(), depth),
        Cmd::Table { dir } => {
            let hist = p1disk::build_nibbles(&sp, &dir).context("4-bit table")?;
            println!("D.nib written; slots per value {hist:?}");
        }
        Cmd::Check {
            dir,
            target,
            cache_rows,
        } => {
            let rep =
                p1disk::check_table(&sp, &dir, target.target(), cache_rows).context("check")?;
            println!(
                "check {}: slots per value {:?}; {} bad target slots, {} slots whose smallest \
                 neighbour is not one less; examples {:?}",
                dir.display(),
                rep.hist,
                rep.bad_targets,
                rep.bad_slots,
                rep.examples
            );
            if rep.bad_targets + rep.bad_slots > 0 {
                bail!("the table is not locally consistent");
            }
            println!("LOCALLY CONSISTENT: every slot holds its exact distance");
        }
        Cmd::Level {
            dir,
            target,
            level,
            solve,
        } => level_cmd(&sp, &dir, target.target(), level, solve)?,
        Cmd::Mitm {
            dir,
            target,
            level,
            ball,
            depth,
        } => {
            let l = resolve_level(&sp, &dir, level)?;
            check_listable(&sp, &dir, l)?;
            let states = p1check::expand(&sp, &p1check::level_slots(&sp, &dir, l)?);
            let mut reps: BTreeMap<u64, State> = BTreeMap::new();
            for s in &states {
                reps.entry(p1check::class_key(&s.0, &s.1)).or_insert(*s);
            }
            println!(
                "level {l}: {} states, {} classes; ball of radius {ball} ...",
                states.len(),
                reps.len()
            );
            let b = p1check::class_ball(&sp, target.target(), ball);
            println!(
                "ball: classes per distance {:?}",
                b.iter()
                    .map(std::collections::HashSet::len)
                    .collect::<Vec<_>>()
            );
            let mut lowest = usize::MAX;
            for (i, (k, s)) in reps.iter().enumerate() {
                let r = p1check::mitm(s, &b, depth);
                println!("class {i} (key {k}): shortest through the ball {r:?}");
                lowest = lowest.min(r.unwrap_or(usize::MAX));
            }
            if lowest > ball + depth {
                println!(
                    "LOWER BOUND: no state of level {l} is within {} moves of the target",
                    ball + depth
                );
            } else {
                println!("lowest found {lowest} (<= {})", ball + depth);
            }
        }
        Cmd::BallCheck { dir, radius } => {
            let ball = omelettehedron::chen::radius1::t1_ball(radius).context("ball")?;
            let tab = Nibbles::open(&dir, sp.rows(), false)?;
            let res: Vec<(u64, u64, u64, u64)> = ball
                .cfgs
                .par_iter()
                .enumerate()
                .map(|(c, &(xa, ep))| -> Result<(u64, u64, u64, u64)> {
                    let (row, col) = sp.slot_of(&xa, &ep);
                    let v = tab.get(row, col)?;
                    let bad_in = u64::from(v != ball.cfg_d1(c as u32));
                    let (mut out, mut bad_out) = (0u64, 0u64);
                    let cube = Cube {
                        ep,
                        xa,
                        ..Cube::solved()
                    };
                    for m in 0..omelettehedron::cube::N_MOVES {
                        if ball.next[c][m] == u32::MAX {
                            let n = cube.apply(m);
                            let (r2, c2) = sp.slot_of(&n.xa, &n.ep);
                            out += 1;
                            if tab.get(r2, c2)? != radius + 1 {
                                bad_out += 1;
                            }
                        }
                    }
                    Ok((1, bad_in, out, bad_out))
                })
                .collect::<Result<_>>()?;
            let sum = res.iter().fold((0, 0, 0, 0), |a, b| {
                (a.0 + b.0, a.1 + b.1, a.2 + b.2, a.3 + b.3)
            });
            println!(
                "ball of radius {radius} around T1: {} configurations (fixed frame), {} with a \
                 label other than the table's; {} moves out of the ball, {} not to {}",
                sum.0,
                sum.1,
                sum.2,
                sum.3,
                radius + 1
            );
            if sum.1 + sum.3 > 0 {
                bail!("the ball and the table disagree");
            }
        }
        Cmd::AEnd { dir, t1 } => {
            use omelettehedron::cube::{is_a_move, N_MOVES};
            let l = resolve_level(&sp, &dir, None)?;
            check_listable(&sp, &dir, l)?;
            let states = p1check::expand(&sp, &p1check::level_slots(&sp, &dir, l)?);
            let tab = Nibbles::open(&t1, sp.rows(), false)?;
            let o: std::collections::HashSet<State> =
                sp.target_states(Target::Oriented).into_iter().collect();
            let sols: Vec<Option<Vec<usize>>> = states
                .par_iter()
                .map(|s| -> Result<Option<Vec<usize>>> {
                    let Some(mut sol) = p1check::descend(&sp, &tab, s, (l - 1) as u8)? else {
                        return Ok(None);
                    };
                    let c = Cube {
                        ep: s.1,
                        xa: s.0,
                        ..Cube::solved()
                    }
                    .apply_seq(&sol);
                    let last = (0..N_MOVES).filter(|&m| is_a_move(m)).find(|&m| {
                        let n = c.apply(m);
                        o.contains(&(n.xa, n.ep))
                    });
                    Ok(last.map(|m| {
                        sol.push(m);
                        sol
                    }))
                })
                .collect::<Result<_>>()?;
            let path = dir.join(format!("solutions{l}-a.tsv"));
            let mut f = std::io::BufWriter::new(std::fs::File::create(&path)?);
            writeln!(
                f,
                "# stickers\tsolution (face turns, the last a class-A turn)\tinternal moves"
            )?;
            let mut bad = 0;
            for (s, sol) in states.iter().zip(&sols) {
                match sol {
                    Some(sol) if sol.len() == l => {
                        let phys =
                            format_moves(&solution_to_physical(Frame::identity(), sol, false));
                        writeln!(f, "{}\t{phys}\t{sol:?}", stickers(s))?;
                    }
                    _ => bad += 1,
                }
            }
            f.flush()?;
            println!(
                "level {l}: {} of {} states have a {l}-move solution ending in a class-A move \
                 (a {}-move path to T1 from the T1 table); written to {}",
                states.len() - bad,
                states.len(),
                l - 1,
                path.display()
            );
            if bad > 0 {
                bail!("{bad} states without one");
            }
        }
        Cmd::Strict => {
            let t0 = std::time::Instant::now();
            let (_, strict, hist) = p1check::strict_exact();
            let mut states = vec![0u128; hist.len()];
            let n = omelettehedron::coord::N_ECOSET;
            for (row, &w) in sp.weight.iter().enumerate() {
                let mut h = [0u64; 16];
                for k in 0..n {
                    h[p1check::strict_get(&strict, row * n + k) as usize] += 1;
                }
                for (d, s) in states.iter_mut().enumerate() {
                    *s += u128::from(h[d]) * u128::from(w);
                }
            }
            println!(
                "strict table: entries per distance {hist:?}; states per distance {states:?} \
                 ({:.0} s)",
                t0.elapsed().as_secs_f64()
            );
        }
        Cmd::Cross { o, t1, giant } => {
            let t0 = std::time::Instant::now();
            let (_, strict, hist) = p1check::strict_exact();
            println!(
                "strict table (two-phase index): entries per distance {hist:?} ({:.0} s)",
                t0.elapsed().as_secs_f64()
            );
            if let Some(g) = &giant {
                let gt = p1check::Giant::open(g, sp.slots()).context("giant table")?;
                println!("giant table header: entries per distance {:?}", gt.hist);
            }
            let c = p1check::cross(&sp, &o, &t1, &strict, giant.as_deref()).context("cross")?;
            println!("{c:#?}");
        }
        Cmd::Status { dir } => {
            let lv = p1disk::levels(&sp, &dir).context("levels")?;
            if lv.is_empty() {
                bail!("no levels in {}", dir.display());
            }
            print_levels(&lv);
        }
    }
    Ok(())
}

fn plain(sp: &P1Space, target: Target, depth: usize) {
    use omelettehedron::cube::{Cube, N_MOVES};
    use rayon::prelude::*;
    use std::collections::HashSet;
    let starts = sp.target_states(target);
    let mut seen: HashSet<u64> = starts.iter().map(|(a, e)| state_key(a, e)).collect();
    let mut front: Vec<State> = starts;
    println!("plain {} depth 0: {} states", target.name(), front.len());
    for d in 1..=depth {
        let mut next: Vec<State> = vec![];
        // a million states at a time, so that the candidates stay small
        for part in front.chunks(1 << 20) {
            let mut cand: Vec<(u64, State)> = part
                .par_iter()
                .flat_map_iter(|&(xa, ep)| {
                    let c = Cube {
                        xa,
                        ep,
                        ..Cube::solved()
                    };
                    (0..N_MOVES).map(move |m| {
                        let n = c.apply(m);
                        (state_key(&n.xa, &n.ep), (n.xa, n.ep))
                    })
                })
                .collect();
            cand.par_sort_unstable_by_key(|x| x.0);
            cand.dedup_by_key(|x| x.0);
            next.extend(
                cand.into_iter()
                    .filter(|(k, _)| seen.insert(*k))
                    .map(|(_, s)| s),
            );
        }
        front = next;
        println!("plain {} depth {d}: {} states", target.name(), front.len());
    }
}

fn resolve_level(sp: &P1Space, dir: &std::path::Path, level: Option<usize>) -> Result<usize> {
    let lv = p1disk::levels(sp, dir)?;
    match level {
        Some(l) => Ok(l),
        None => lv.iter().rposition(|l| l.slots > 0).context("no levels"),
    }
}

fn stickers(s: &State) -> String {
    let c = Cube {
        xa: s.0,
        ep: s.1,
        ..Cube::solved()
    };
    c.to_stickers()
        .iter()
        .map(|&x| char::from(b'0' + x))
        .collect()
}

/// Levels are listed in memory (16 bytes a slot, 24 a state): refuse one too large to list.
fn check_listable(sp: &P1Space, dir: &std::path::Path, l: usize) -> Result<()> {
    const MAX_STATES: u128 = 50_000_000;
    let lv = p1disk::levels(sp, dir)?;
    let states = lv.get(l).map_or(0, |x| x.states);
    if states > MAX_STATES {
        bail!("level {l} has {states} states: too many to list (at most {MAX_STATES})");
    }
    Ok(())
}

fn level_cmd(
    sp: &P1Space,
    dir: &std::path::Path,
    target: Target,
    level: Option<usize>,
    solve: bool,
) -> Result<()> {
    let l = resolve_level(sp, dir, level)?;
    check_listable(sp, dir, l)?;
    let slots = p1check::level_slots(sp, dir, l)?;
    let states = p1check::expand(sp, &slots);
    let keys: Vec<u64> = states
        .par_iter()
        .map(|s| p1check::class_key(&s.0, &s.1))
        .collect();
    let mut classes: BTreeMap<u64, usize> = BTreeMap::new();
    for &k in &keys {
        *classes.entry(k).or_default() += 1;
    }
    let held = states.len() as f64 / 12.0;
    println!(
        "level {l}: {} slots, {} states (fixed frame), {held} held-edge states (/12), {} classes \
         under the 288 pairs, {} positions",
        slots.len(),
        states.len(),
        classes.len(),
        states.len() as u128 * PER_STATE
    );
    let mut sizes: BTreeMap<usize, usize> = BTreeMap::new();
    for &n in classes.values() {
        *sizes.entry(n).or_default() += 1;
    }
    println!("classes by size (states: classes): {sizes:?}");
    let path = dir.join(format!("level{l}.tsv"));
    let mut f = std::io::BufWriter::new(std::fs::File::create(&path)?);
    writeln!(
        f,
        "# class_key\tstickers (face index per sticker; corners and B-centres solved)"
    )?;
    for (s, k) in states.iter().zip(&keys) {
        writeln!(f, "{k}\t{}", stickers(s))?;
    }
    f.flush()?;
    println!("listed in {}", path.display());
    if solve {
        let tab = Nibbles::open(dir, sp.rows(), false)?;
        let sols: Vec<Option<Vec<usize>>> = states
            .par_iter()
            .map(|s| p1check::descend(sp, &tab, s, l as u8))
            .collect::<Result<_, _>>()?;
        let ts: std::collections::HashSet<State> = sp.target_states(target).into_iter().collect();
        let path = dir.join(format!("solutions{l}.tsv"));
        let mut f = std::io::BufWriter::new(std::fs::File::create(&path)?);
        writeln!(f, "# stickers\tsolution (face turns)\tinternal moves")?;
        let mut bad = 0;
        for (s, sol) in states.iter().zip(&sols) {
            let Some(sol) = sol else {
                bad += 1;
                continue;
            };
            let c = Cube {
                xa: s.0,
                ep: s.1,
                ..Cube::solved()
            }
            .apply_seq(sol);
            if sol.len() != l || !ts.contains(&(c.xa, c.ep)) {
                bad += 1;
            }
            let phys = format_moves(&solution_to_physical(Frame::identity(), sol, false));
            writeln!(f, "{}\t{phys}\t{sol:?}", stickers(s))?;
        }
        f.flush()?;
        println!(
            "solutions of {} moves for {} of {} states (in the cube model they reach the target); \
             {bad} failures; written to {}",
            l,
            states.len() - bad,
            states.len(),
            path.display()
        );
        if bad > 0 {
            bail!("{bad} states without a solution of {l} moves");
        }
    }
    Ok(())
}

fn print_levels(lv: &[p1disk::Level]) {
    println!("depth  slots  states(fixed frame)  states(held edge, /12)  positions");
    let (mut ts, mut tst) = (0u64, 0u128);
    for l in lv {
        ts += l.slots;
        tst += l.states;
        println!(
            "{:>2} {:>15} {:>18} {:>18} {:>26}",
            l.d,
            l.slots,
            l.states,
            l.states / 12,
            l.states * PER_STATE
        );
    }
    let all = 369_600u128 * N_EPERM as u128;
    println!(
        "total {ts} slots, {tst} states of {all} ({}), {} positions",
        if tst == all { "all" } else { "NOT all" },
        tst * PER_STATE
    );
}
