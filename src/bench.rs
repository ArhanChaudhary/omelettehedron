//! Benchmark harness: solve seeded random positions with a time limit and summarise.

use crate::cube::Cube;
use crate::notation::Frame;
use crate::search::*;
use crate::tables::*;
use rand::SeedableRng;
use serde::Serialize;
use std::time::{Duration, TryFromFloatSecsError};

#[derive(Serialize, Clone, Debug)]
pub struct ScrambleResult {
    pub index: usize,
    pub length: Option<usize>,
    pub first_time: Option<f64>,
    pub first_length: Option<usize>,
    pub optimal: bool,
    pub elapsed: f64,
    pub nodes: u64,
    pub solutions: u64,
    pub events: Vec<(f64, usize)>,
}

#[derive(Serialize, Clone, Debug)]
pub struct Summary {
    pub tier: String,
    /// the move set (`face`, `face+slice`, ...)
    pub moves: String,
    pub threads: usize,
    pub scrambles: usize,
    pub time_limit: Option<f64>,
    pub seed: u64,
    /// (seconds, mean best length over positions solved by then, number solved by then)
    pub checkpoints: Vec<(f64, f64, usize)>,
    pub mean_length: f64,
    pub mean_first_time: f64,
    pub mean_first_length: f64,
    pub mean_elapsed: f64,
    pub nodes_per_sec: f64,
    pub solved_optimally: usize,
    pub results: Vec<ScrambleResult>,
}

pub fn random_positions(n: usize, seed: u64) -> Vec<Cube> {
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(seed);
    (0..n).map(|_| Cube::random(&mut rng)).collect()
}

/// Seeded random-move scrambles of `len` face turns (no face twice in a row).
pub fn scrambled_positions(n: usize, len: usize, seed: u64) -> Vec<Cube> {
    use rand::Rng;
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(seed);
    (0..n)
        .map(|_| {
            let mut moves = vec![];
            let mut last = 8;
            while moves.len() < len {
                let face = rng.gen_range(0..8);
                if face != last {
                    moves.push(crate::notation::PhysMove::face(face, rng.gen()));
                    last = face;
                }
            }
            crate::notation::scramble_to_cube(&moves).0
        })
        .collect()
}

fn length_at(events: &[(f64, usize)], t: f64) -> Option<usize> {
    events.iter().filter(|e| e.0 <= t).map(|e| e.1).min()
}

pub fn run<I: AeIndex>(
    tables: &Tables<I>,
    tier: &str,
    n: usize,
    seed: u64,
    scramble_len: Option<usize>,
    limits: &Limits,
    verbose: bool,
) -> Summary {
    let positions = match scramble_len {
        Some(len) => scrambled_positions(n, len, seed),
        None => random_positions(n, seed),
    };
    let mut results = vec![];
    for (i, c) in positions.iter().enumerate() {
        let o = solve(tables, *c, Frame::identity(), limits, Report::Quiet);
        let events: Vec<(f64, usize)> = o.events.iter().map(|e| (e.time, e.length)).collect();
        let r = ScrambleResult {
            index: i,
            length: o.length,
            first_time: events.first().map(|e| e.0),
            first_length: events.first().map(|e| e.1),
            optimal: o.optimal,
            elapsed: o.elapsed,
            nodes: o.nodes,
            solutions: o.solutions_found,
            events,
        };
        if verbose {
            eprintln!(
                "#{i:<4} best {:>3} first {:>3} @ {:>8.4}s  elapsed {:>7.3}s  nodes {:>12}  sols {:>8}{}",
                r.length.map_or("-".into(), |l| l.to_string()),
                r.first_length.map_or("-".into(), |l| l.to_string()),
                r.first_time.unwrap_or(f64::NAN),
                r.elapsed,
                r.nodes,
                r.solutions,
                if r.optimal { "  optimal" } else { "" }
            );
        }
        results.push(r);
    }
    let tl = limits.time.map(|d| d.as_secs_f64());
    let mut cps: Vec<f64> = [0.01, 0.1, 1.0, 10.0, 60.0, 600.0]
        .into_iter()
        .filter(|&c| tl.is_none_or(|t| c < t))
        .collect();
    if let Some(t) = tl {
        cps.push(t);
    }
    let mean = |f: &dyn Fn(&ScrambleResult) -> Option<f64>| -> f64 {
        let v: Vec<f64> = results.iter().filter_map(f).collect();
        if v.is_empty() {
            f64::NAN
        } else {
            v.iter().sum::<f64>() / v.len() as f64
        }
    };
    let checkpoints = cps
        .iter()
        .map(|&c| {
            let solved = results
                .iter()
                .filter(|r| length_at(&r.events, c).is_some())
                .count();
            (
                c,
                mean(&|r: &ScrambleResult| length_at(&r.events, c).map(|l| l as f64)),
                solved,
            )
        })
        .collect();
    let total_nodes: u64 = results.iter().map(|r| r.nodes).sum();
    let total_time: f64 = results.iter().map(|r| r.elapsed).sum();
    Summary {
        tier: tier.to_string(),
        moves: tables.moves.name(),
        threads: rayon::current_num_threads(),
        scrambles: n,
        time_limit: tl,
        seed,
        checkpoints,
        mean_length: mean(&|r| r.length.map(|l| l as f64)),
        mean_first_time: mean(&|r| r.first_time),
        mean_first_length: mean(&|r| r.first_length.map(|l| l as f64)),
        mean_elapsed: mean(&|r| Some(r.elapsed)),
        nodes_per_sec: total_nodes as f64 / total_time,
        solved_optimally: results.iter().filter(|r| r.optimal).count(),
        results,
    }
}

pub fn print_summary(s: &Summary) {
    println!(
        "tier {}  moves {}  threads {}  scrambles {}  seed {}  time limit {}",
        s.tier,
        s.moves,
        s.threads,
        s.scrambles,
        s.seed,
        s.time_limit.map_or("none".into(), |t| format!("{t}s"))
    );
    for (c, l, k) in &s.checkpoints {
        println!(
            "  after {:>6}s: mean best length {:.3} ({}/{} solved)",
            c, l, k, s.scrambles
        );
    }
    println!("  mean final length:        {:.3}", s.mean_length);
    println!(
        "  mean first solution:      {:.3} moves after {:.4}s",
        s.mean_first_length, s.mean_first_time
    );
    println!("  mean elapsed:             {:.3}s", s.mean_elapsed);
    println!(
        "  proven optimal:           {}/{}",
        s.solved_optimally, s.scrambles
    );
    println!("  nodes/s:                  {:.3e}", s.nodes_per_sec);
}

/// Stopping rules from command-line values.
///
/// # Errors
///
/// If `secs` is negative, not finite or too large for a [`Duration`].
pub fn limits_with_time(
    secs: Option<f64>,
    max_len: Option<usize>,
    max_sol: Option<u64>,
) -> Result<Limits, TryFromFloatSecsError> {
    Ok(Limits {
        time: secs.map(Duration::try_from_secs_f64).transpose()?,
        max_length: max_len,
        max_solutions: max_sol,
        ..Limits::default()
    })
}
