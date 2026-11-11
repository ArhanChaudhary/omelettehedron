//! Benchmark harness of the optimal solver: solve seeded positions and summarise.

use super::search::{solve, Limits, Stop};
use super::tables::Bounds;
use crate::cube::Cube;
use crate::notation::Frame;
use serde::Serialize;

#[derive(Serialize, Clone, Debug)]
pub struct PositionResult {
    pub index: usize,
    /// the optimal length, if found
    pub length: Option<usize>,
    /// no solution is shorter (proven)
    pub lower_bound: usize,
    pub stop: Stop,
    pub nodes: u64,
    /// table lookups
    pub lookups: u64,
    /// children by the set of lookups that prune them, and by moves left (see
    /// `search::Outcome`)
    pub prune_masks: Vec<u64>,
    pub children_by_left: Vec<u64>,
    pub elapsed: f64,
    /// nodes of each iteration (depth, nodes)
    pub iterations: Vec<(usize, u64)>,
}

#[derive(Serialize, Clone, Debug)]
pub struct Summary {
    /// always "optimal" (the solver), to tell these results from the two-phase solver's
    pub solver: String,
    pub tier: String,
    pub threads: usize,
    pub positions: usize,
    pub seed: u64,
    /// random-move scrambles of this many face turns (none: random positions)
    pub scramble_len: Option<usize>,
    pub time_limit: Option<f64>,
    pub solved: usize,
    pub mean_length: f64,
    pub total_nodes: u64,
    pub total_lookups: u64,
    pub total_time: f64,
    pub nodes_per_sec: f64,
    pub results: Vec<PositionResult>,
}

/// Solve `n` seeded positions (random-move scrambles of `scramble_len` face turns, the same as
/// the two-phase solver's `bench --scramble-len`, or random positions) and summarise; each
/// position's result is printed as it completes.
pub fn run<B: Bounds>(
    b: &B,
    tier: &str,
    n: usize,
    seed: u64,
    scramble_len: Option<usize>,
    limits: &Limits,
    verbose: bool,
) -> Summary {
    let positions: Vec<Cube> = match scramble_len {
        Some(len) => crate::bench::scrambled_positions(n, len, seed),
        None => crate::bench::random_positions(n, seed),
    };
    let mut results = vec![];
    for (i, c) in positions.iter().enumerate() {
        let o = solve(b, *c, Frame::identity(), limits, verbose);
        let r = PositionResult {
            index: i,
            length: o.length,
            lower_bound: o.lower_bound,
            stop: o.stop,
            nodes: o.nodes,
            lookups: o.lookups,
            prune_masks: o.prune_masks.clone(),
            children_by_left: o.children_by_left.clone(),
            elapsed: o.elapsed,
            iterations: o.iterations.iter().map(|it| (it.depth, it.nodes)).collect(),
        };
        println!(
            "#{i:<3} optimal {:>3}  nodes {:>14}  lookups {:>14}  {:>9.3} s{}",
            r.length
                .map_or_else(|| format!(">={}", r.lower_bound), |l| l.to_string()),
            r.nodes,
            r.lookups,
            r.elapsed,
            match r.stop {
                Stop::Solved => "",
                Stop::Found => "  (start depth: not proven optimal)",
                Stop::Time => "  (time limit)",
                Stop::MaxLength => "  (longer than the length limit)",
                Stop::BadTables => "  (bad tables)",
            }
        );
        results.push(r);
    }
    let solved: Vec<&PositionResult> = results.iter().filter(|r| r.length.is_some()).collect();
    let total_nodes: u64 = results.iter().map(|r| r.nodes).sum();
    let total_lookups: u64 = results.iter().map(|r| r.lookups).sum();
    let total_time: f64 = results.iter().map(|r| r.elapsed).sum();
    Summary {
        solver: "optimal".into(),
        tier: tier.into(),
        threads: rayon::current_num_threads(),
        positions: n,
        seed,
        scramble_len,
        time_limit: limits.time.map(|d| d.as_secs_f64()),
        solved: solved.len(),
        mean_length: solved
            .iter()
            .filter_map(|r| r.length)
            .map(|l| l as f64)
            .sum::<f64>()
            / solved.len().max(1) as f64,
        total_nodes,
        total_lookups,
        total_time,
        nodes_per_sec: total_nodes as f64 / total_time.max(1e-9),
        results,
    }
}

pub fn print_summary(s: &Summary) {
    println!(
        "optimal solver: tier {}  threads {}  positions {}  seed {}  {}  time limit {}",
        s.tier,
        s.threads,
        s.positions,
        s.seed,
        s.scramble_len.map_or_else(
            || "random positions".into(),
            |l| format!("{l}-move scrambles")
        ),
        s.time_limit.map_or("none".into(), |t| format!("{t}s"))
    );
    println!("  solved optimally:  {}/{}", s.solved, s.positions);
    println!("  mean length:       {:.3}", s.mean_length);
    println!("  total nodes:       {}", s.total_nodes);
    println!(
        "  total lookups:     {} ({:.2} per node)",
        s.total_lookups,
        s.total_lookups as f64 / s.total_nodes.max(1) as f64
    );
    println!(
        "  mean nodes:        {:.4e}",
        s.total_nodes as f64 / s.positions.max(1) as f64
    );
    println!(
        "  total time:        {:.3} s ({:.3} s per position)",
        s.total_time,
        s.total_time / s.positions.max(1) as f64
    );
    println!("  nodes/s:           {:.3e}", s.nodes_per_sec);
}
