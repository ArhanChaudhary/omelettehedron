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
    /// Positions at the table's last level on both face classes: side A is the position's
    /// A-centres and edges, side B the A-centres and edges of its `SWAP_AB` conjugate (the
    /// conjugate's edges are a function `g` of the edges, its A-centres a bijection of the
    /// B-centres). Positions are held with the front corner in place (960 corner states), as
    /// the relaxed target's both-sides check: both sides at the maximum needs an edge
    /// arrangement `ep` of a maximal state with `g(ep)` one too.
    BothSides {
        #[arg(long)]
        dir: PathBuf,
        #[arg(long)]
        level: Option<usize>,
        /// Positions built from the maximal states (random corners, B-centres) whose side B is
        /// read off the table
        #[arg(long, default_value_t = 20_000)]
        samples: usize,
    },
    /// Every position with both sides at the last level (each listed triple with each of the
    /// 960 corner states), up to the 24 class-preserving symmetries as `chen-p2bfs
    /// conjugate-solutions` takes orbits: one pattern per class to `--out`, and the orbit sizes,
    /// whose sum is checked against the count
    BothClasses {
        #[arg(long)]
        dir: PathBuf,
        #[arg(long)]
        level: Option<usize>,
        #[arg(long)]
        out: PathBuf,
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
        Cmd::BothSides {
            dir,
            level,
            samples,
        } => both_sides(&sp, &dir, level, samples)?,
        Cmd::BothClasses { dir, level, out } => both_classes(&sp, &dir, level, &out)?,
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

/// The states of level `l`: from its bit plane, or (once the planes are removed) from the list
/// that `level` wrote (`<dir>/level<l>.tsv`), checked against the level's count.
fn level_states(sp: &P1Space, dir: &std::path::Path, l: usize) -> Result<Vec<State>> {
    if dir.join(format!("L{l:02}.bits")).exists() {
        return Ok(p1check::expand(sp, &p1check::level_slots(sp, dir, l)?));
    }
    let path = dir.join(format!("level{l}.tsv"));
    let text =
        std::fs::read_to_string(&path).with_context(|| format!("read {}", path.display()))?;
    let mut out = vec![];
    for line in text.lines().filter(|x| !x.starts_with('#')) {
        let st = line.split('\t').nth(1).context("a line without stickers")?;
        let b = st.as_bytes();
        if b.len() != 72 {
            bail!("{}: a line of {} stickers", path.display(), b.len());
        }
        let s: [u8; 72] = std::array::from_fn(|i| b[i].wrapping_sub(b'0'));
        let c = Cube::from_stickers(&s).context("not a position")?;
        out.push((c.xa, c.ep));
    }
    let want = p1disk::levels(sp, dir)?.get(l).map_or(0, |x| x.states);
    if out.len() as u128 != want {
        bail!(
            "{}: {} states, level {l} has {want}",
            path.display(),
            out.len()
        );
    }
    Ok(out)
}

/// A position's parts (A-centres, edges, B-centres, corners), as a key for ordering.
type PosKey = ([u8; 12], [u8; 12], [u8; 12], [u8; 6], [u8; 6]);

/// A listed position's side-A and side-B distances and its class (key and representative).
type Checked = (u8, u8, Option<(PosKey, Cube)>);

/// A position conjugated by a class-preserving symmetry `m` and held again with the front
/// corner in place (one of the 12 rotations does it).
fn conj_held(c: &Cube, m: &omelettehedron::geometry::Mat, rots: &[Vec<usize>]) -> Option<Cube> {
    use omelettehedron::geometry::{face_of_normal, mat_vec, FACES};
    let d = omelettehedron::cube::defs();
    let st = c.to_stickers();
    let perm = d.sm.whole(m);
    let mut n = [0u8; 72];
    for i in 0..72 {
        n[perm[i]] = face_of_normal(mat_vec(m, FACES[st[i] as usize])) as u8;
    }
    rots.iter().find_map(|r| {
        let mut h = [0u8; 72];
        for i in 0..72 {
            h[r[i]] = n[i];
        }
        Cube::from_stickers(&h).filter(|x| x.cp[0] == 0 && x.co[0] == 0)
    })
}

/// Every (A-centres, edges, B-centres) of a position with both sides at the maximum (corners
/// solved; any of the 960 corner states goes with each), both distances read off the table for
/// each, and their classes under the 24 class-preserving symmetries (conjugation, held again):
/// `<dir>/both<l>.tsv`, one line per class (size, pattern with solved corners).
fn both_list(
    sp: &P1Space,
    dir: &std::path::Path,
    l: usize,
    tab: &Nibbles,
    both: &[[u8; 12]],
    by_ep: &std::collections::HashMap<[u8; 12], Vec<[u8; 12]>>,
) -> Result<()> {
    use omelettehedron::cube::SWAP_AB;
    use omelettehedron::geometry::{class_preserving_group, rotation_group, transpose};
    let d = omelettehedron::cube::defs();
    let rots: Vec<Vec<usize>> = rotation_group().iter().map(|r| d.sm.whole(r)).collect();
    let syms = class_preserving_group();
    let mut triples: Vec<Cube> = vec![];
    for ep in both {
        let gb = Cube {
            ep: *ep,
            ..Cube::solved()
        }
        .conjugate(&SWAP_AB)
        .context("SWAP_AB")?
        .ep;
        for xa in by_ep.get(ep).map_or(&[][..], |v| &v[..]) {
            for xa2 in by_ep.get(&gb).map_or(&[][..], |v| &v[..]) {
                let xb = Cube {
                    xa: *xa2,
                    ..Cube::solved()
                }
                .conjugate(&transpose(&SWAP_AB))
                .context("SWAP_AB^-1")?
                .xb;
                triples.push(Cube {
                    ep: *ep,
                    xa: *xa,
                    xb,
                    ..Cube::solved()
                });
            }
        }
    }
    let key = |c: &Cube| -> PosKey { (c.xa, c.ep, c.xb, c.cp, c.co) };
    let checked: Vec<Checked> = triples
        .par_iter()
        .map(|c| -> Result<_> {
            let (ra, ca) = sp.slot_of(&c.xa, &c.ep);
            let y = c.conjugate(&SWAP_AB).context("SWAP_AB")?;
            let (rb, cb) = sp.slot_of(&y.xa, &y.ep);
            let canon = syms
                .iter()
                .filter_map(|m| conj_held(c, m, &rots))
                .map(|x| (key(&x), x))
                .min_by_key(|x| x.0);
            Ok((tab.get(ra, ca)?, tab.get(rb, cb)?, canon))
        })
        .collect::<Result<_>>()?;
    let bad = checked
        .iter()
        .filter(|x| usize::from(x.0) != l || usize::from(x.1) != l || x.2.is_none())
        .count();
    let mut classes: std::collections::BTreeMap<_, (usize, Cube)> =
        std::collections::BTreeMap::new();
    for (_, _, c) in &checked {
        if let Some((k, x)) = c {
            classes.entry(*k).or_insert((0, *x)).0 += 1;
        }
    }
    let path = dir.join(format!("both{l}.tsv"));
    let mut f = std::io::BufWriter::new(std::fs::File::create(&path)?);
    writeln!(
        f,
        "# positions with both sides at {l}, corners solved (each goes with all 960 corner states), \
         up to the 24 class-preserving symmetries: class size, pattern"
    )?;
    for (n, x) in classes.values() {
        writeln!(
            f,
            "{n}\t{}",
            omelettehedron::pattern::format(&omelettehedron::pattern::of_cube(x))
        )?;
    }
    f.flush()?;
    println!(
        "{} (A-centres, edges, B-centres) with both sides at {l} (x 960 corner states = {} \
         positions); both distances read off the table for each: {} not both {l}; {} classes \
         under the 24 symmetries, listed in {}",
        triples.len(),
        triples.len() * 960,
        bad,
        classes.len(),
        path.display()
    );
    if bad > 0 {
        bail!("{bad} triples are not at {l} on both sides");
    }
    Ok(())
}

/// A position in 120 bits (as `chen-p2bfs`'s), for ordering and distinctness.
fn packed(c: &Cube) -> u128 {
    let fields: [(&[u8], u32); 5] = [(&c.cp, 3), (&c.co, 1), (&c.ep, 4), (&c.xa, 2), (&c.xb, 2)];
    let mut v = 0u128;
    for (values, bits) in fields {
        for &x in values {
            v = v << bits | u128::from(x);
        }
    }
    v
}

fn both_classes(
    sp: &P1Space,
    dir: &std::path::Path,
    level: Option<usize>,
    out: &std::path::Path,
) -> Result<()> {
    use omelettehedron::cube::SWAP_AB;
    use omelettehedron::geometry::transpose;
    use std::collections::HashMap;
    let l = resolve_level(sp, dir, level)?;
    check_listable(sp, dir, l)?;
    let states = level_states(sp, dir, l)?;
    let mut by_ep: HashMap<[u8; 12], Vec<[u8; 12]>> = HashMap::new();
    for (xa, ep) in &states {
        by_ep.entry(*ep).or_default().push(*xa);
    }
    let g = |ep: &[u8; 12]| -> Result<[u8; 12]> {
        Ok(Cube {
            ep: *ep,
            ..Cube::solved()
        }
        .conjugate(&SWAP_AB)
        .context("SWAP_AB")?
        .ep)
    };
    let mut triples: Vec<Cube> = vec![];
    let mut eps: Vec<&[u8; 12]> = by_ep.keys().collect();
    eps.sort_unstable();
    for ep in eps {
        let gb = g(ep)?;
        let Some(bs) = by_ep.get(&gb) else { continue };
        for xa in by_ep.get(ep).map_or(&[][..], |v| &v[..]) {
            for xa2 in bs {
                let xb = Cube {
                    xa: *xa2,
                    ..Cube::solved()
                }
                .conjugate(&transpose(&SWAP_AB))
                .context("SWAP_AB^-1")?
                .xb;
                triples.push(Cube {
                    ep: *ep,
                    xa: *xa,
                    xb,
                    ..Cube::solved()
                });
            }
        }
    }
    let corner_states = &omelettehedron::coord::corners().id_to_state;
    let s48 = omelettehedron::optimal::sym::sym48();
    let err = |e: std::io::Error| anyhow::anyhow!("{}: {e}", out.display());
    let mut w = std::io::BufWriter::new(std::fs::File::create(out).map_err(err)?);
    let (mut positions, mut classes) = (0u64, 0u64);
    let mut sizes: BTreeMap<usize, u64> = BTreeMap::new();
    for chunk in triples.chunks(256) {
        let found: Vec<Vec<(String, usize)>> = chunk
            .par_iter()
            .map(|t| -> Result<Vec<(String, usize)>> {
                let mut reps = vec![];
                for &(cp, co) in corner_states {
                    let c = Cube { cp, co, ..*t };
                    let k = packed(&c);
                    let mut images: Vec<u128> = Vec::with_capacity(24);
                    let mut rep = true;
                    for sy in 0..24 {
                        let (cc, _) = s48.conj_cube(&c, sy).context("no conjugate")?;
                        let kk = packed(&cc);
                        if kk < k {
                            rep = false;
                            break;
                        }
                        if !images.contains(&kk) {
                            images.push(kk);
                        }
                    }
                    if rep {
                        reps.push((
                            omelettehedron::pattern::format(&omelettehedron::pattern::of_cube(&c)),
                            images.len(),
                        ));
                    }
                }
                Ok(reps)
            })
            .collect::<Result<_>>()?;
        for (pat, n) in found.into_iter().flatten() {
            writeln!(w, "{pat}").map_err(err)?;
            positions += n as u64;
            classes += 1;
            *sizes.entry(n).or_default() += 1;
        }
    }
    w.flush().map_err(err)?;
    let want = triples.len() as u64 * corner_states.len() as u64;
    println!(
        "level {l}: {} triples x {} corner states = {want} positions with both sides at {l}; {classes} \
         classes under the 24 symmetries (orbit sizes: classes {sizes:?}) holding {positions} \
         positions; listed in {}",
        triples.len(),
        corner_states.len(),
        out.display()
    );
    if positions != want {
        bail!("the orbits hold {positions} positions, not {want}");
    }
    Ok(())
}

fn both_sides(
    sp: &P1Space,
    dir: &std::path::Path,
    level: Option<usize>,
    samples: usize,
) -> Result<()> {
    use omelettehedron::cube::SWAP_AB;
    use omelettehedron::geometry::transpose;
    use rand::{Rng, SeedableRng};
    use std::collections::{HashMap, HashSet};
    let conj = |c: &Cube| {
        c.conjugate(&SWAP_AB)
            .context("SWAP_AB keeps the front corner")
    };
    let unconj = |c: &Cube| {
        c.conjugate(&transpose(&SWAP_AB))
            .context("SWAP_AB^-1 keeps the front corner")
    };
    let g = |ep: &[u8; 12]| -> Result<[u8; 12]> {
        Ok(conj(&Cube {
            ep: *ep,
            ..Cube::solved()
        })?
        .ep)
    };
    // the conjugate's edges depend only on the edges, its A-centres only on the B-centres
    let mut rng = rand::rngs::StdRng::seed_from_u64(1);
    for _ in 0..100_000 {
        let x = Cube::random(&mut rng);
        let y = conj(&x)?;
        let ya = conj(&Cube {
            xb: x.xb,
            ..Cube::solved()
        })?
        .xa;
        if y.ep != g(&x.ep)? || y.xa != ya || unconj(&y)? != x {
            bail!("the conjugate does not split into edges and centres");
        }
    }
    let l = resolve_level(sp, dir, level)?;
    check_listable(sp, dir, l)?;
    let states = level_states(sp, dir, l)?;
    // A-colourings per edge arrangement among the maximal states
    let mut by_ep: HashMap<[u8; 12], Vec<[u8; 12]>> = HashMap::new();
    for (xa, ep) in &states {
        by_ep.entry(*ep).or_default().push(*xa);
    }
    let eb: HashSet<[u8; 12]> = by_ep.keys().map(&g).collect::<Result<_>>()?;
    let both: Vec<[u8; 12]> = by_ep
        .keys()
        .filter_map(|ep| g(ep).ok().filter(|e| by_ep.contains_key(e)).map(|_| *ep))
        .collect();
    let positions: u128 = both
        .iter()
        .map(|ep| -> Result<u128> {
            let na = by_ep.get(ep).map_or(0, Vec::len) as u128;
            let nb = by_ep.get(&g(ep)?).map_or(0, Vec::len) as u128;
            Ok(na * nb * 960)
        })
        .sum::<Result<u128>>()?;
    println!(
        "level {l}: {} states (fixed frame) with {} edge arrangements EA; g(EA) {}; EA and \
         g^-1(EA) in common: {}; positions (front corner held) with both sides at {l}: {positions}",
        states.len(),
        by_ep.len(),
        eb.len(),
        both.len()
    );
    // side B of positions built from the maximal states, by the table
    let tab = Nibbles::open(dir, sp.rows(), false)?;
    let picks: Vec<(usize, Cube)> = (0..samples)
        .map(|_| {
            let (xa, ep) = states[rng.gen_range(0..states.len())];
            let r = Cube::random(&mut rng);
            (0, Cube { ep, xa, ..r })
        })
        .collect();
    let mut hist: std::collections::BTreeMap<(u8, u8), usize> = std::collections::BTreeMap::new();
    let res: Vec<(u8, u8)> = picks
        .par_iter()
        .map(|(_, c)| -> Result<(u8, u8)> {
            let (ra, ca) = sp.slot_of(&c.xa, &c.ep);
            let y = conj(c)?;
            let (rb, cb) = sp.slot_of(&y.xa, &y.ep);
            Ok((tab.get(ra, ca)?, tab.get(rb, cb)?))
        })
        .collect::<Result<_>>()?;
    for r in res {
        *hist.entry(r).or_default() += 1;
    }
    println!("{samples} positions built from the states at {l} (random corners and B-centres): (side A, side B) {hist:?}");
    if !both.is_empty() {
        both_list(sp, dir, l, &tab, &both, &by_ep)?;
        // a few positions with both sides at the maximum, both distances read off the table
        for ep in both.iter().take(5) {
            let gb = g(ep)?;
            let xa = by_ep
                .get(ep)
                .and_then(|v| v.first())
                .copied()
                .unwrap_or_default();
            let xa2 = by_ep
                .get(&gb)
                .and_then(|v| v.first())
                .copied()
                .unwrap_or_default();
            let xb = unconj(&Cube {
                xa: xa2,
                ..Cube::solved()
            })?
            .xb;
            let c = Cube {
                ep: *ep,
                xa,
                xb,
                ..Cube::solved()
            };
            let (ra, ca) = sp.slot_of(&c.xa, &c.ep);
            let y = conj(&c)?;
            let (rb, cb) = sp.slot_of(&y.xa, &y.ep);
            println!(
                "  both at {l}: {} (side A {}, side B {})",
                omelettehedron::pattern::format(&omelettehedron::pattern::of_cube(&c)),
                tab.get(ra, ca)?,
                tab.get(rb, cb)?
            );
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
