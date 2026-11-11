//! Search coverage check.
//!
//! With every pruning table zeroed (heuristic = 0 everywhere), the search must reach each
//! phase-1 leaf exactly once on each side. The search runs for two targets at once (side 0:
//! H_A; side 1: H_B, on the position conjugated by SWAP_AB), walking each move sequence once
//! with the sides it is alive on. A leaf of side 0 in the depth-d iteration is a move sequence
//! of length d that follows the move-order rules and ends with a move outside H_A; a leaf of
//! side 1 one that ends with a move outside class B (the moves SWAP_AB maps into class A).
//! Their numbers are counted exactly by a transfer matrix over the last move. A mismatch
//! means that the parallel split, the task scheduling, the batched child evaluation or the
//! per-side bookkeeping drops or repeats part of a side's tree. The tables themselves are not
//! involved.

use omelettehedron::cube::*;
use omelettehedron::mem::HugeBuf;
use omelettehedron::notation::Frame;
use omelettehedron::prune::Table2;
use omelettehedron::search::*;
use omelettehedron::tables::*;
use rand::SeedableRng;

fn zero_table(n: usize) -> Table2 {
    Table2 {
        data: HugeBuf::zeroed(n.div_ceil(32)),
        n,
        hist: vec![1],
    }
}

const SLICE: MoveSet = MoveSet {
    wide: false,
    slice: true,
};

/// Small-tier tables with all pruning values 0, for a move set.
fn null_tables(moves: MoveSet) -> Tables<AeSym> {
    let ae_ix = AeSym::new(false);
    let n = ae_ix.size();
    Tables {
        moves,
        ae_ix,
        ae: zero_table(n),
        // no targets: the exact root walk stops at once and reports distance 0
        ae_targets: vec![],
        cb: CbTable {
            tab: zero_table(N_CB),
            attach: CbTable::attach_targets(),
        },
        p2: Phase2::new(),
        cbe: None,
        cbs: None,
    }
}

/// Allowed move sequences of length `depth` (by the move-order rules) ending with each move.
fn canonical_ways(moves: MoveSet, depth: usize) -> Vec<u64> {
    let n = moves.n_moves();
    // ways[m]: allowed sequences of the current length ending with move m
    let mut ways = vec![1u64; n];
    for _ in 1..depth {
        let mut next = vec![0u64; n];
        for (last, &w) in ways.iter().enumerate() {
            for (m, x) in next.iter_mut().enumerate() {
                if moves.allowed_after(last, m) {
                    *x += w;
                }
            }
        }
        ways = next;
    }
    ways
}

/// Whether `m` ends a phase-1 path on side `side` (0: a move outside class A; 1: a move
/// outside class B, i.e. whose image under `SWAP_AB` is outside class A).
fn ends_path(side: usize, m: usize) -> bool {
    let m = if side == 0 {
        m
    } else {
        conjugate_move(m, &SWAP_AB)
    };
    !is_a_move(m)
}

/// Move sequences of length `depth` allowed by the move-order rules whose last move ends a
/// phase-1 path on `side` (the empty sequence for depth 0).
fn canonical_leaves(moves: MoveSet, side: usize, depth: usize) -> u64 {
    if depth == 0 {
        return 1;
    }
    let ways = canonical_ways(moves, depth);
    (0..moves.n_moves())
        .filter(|&m| ends_path(side, m))
        .map(|m| ways[m])
        .sum()
}

/// Runs the search on `cube` up to depth `d` with each thread count and checks the leaf count
/// of each side.
fn check_leaves(t: &Tables<AeSym>, cube: Cube, expect_solution: bool, d: usize, threads: &[usize]) {
    // both sides (targets H_A and H_B) at every depth 0..=D
    let per_side: [u64; 2] = [0, 1].map(|s| {
        (0..=d)
            .map(|k| canonical_leaves(t.moves, s, k))
            .sum::<u64>()
    });
    let expected = per_side[0] + per_side[1];
    for &threads in threads {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("thread pool");
        let o = pool.install(|| {
            solve(
                t,
                cube,
                Frame::identity(),
                &Limits {
                    max_depth: Some(d),
                    ..Default::default()
                },
                Report::Quiet,
            )
        });
        assert_eq!(o.solution.is_some(), expect_solution);
        assert_eq!(o.leaves_by_side, per_side, "{threads} threads");
        assert_eq!(o.leaves, expected, "{threads} threads");
    }
}

#[test]
fn search_reaches_every_canonical_leaf_once() {
    let t = null_tables(MoveSet::FACE);
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(99);
    // a random position is about 15 moves from H_A, so no leaf is a solution and the search
    // never prunes against a best solution
    check_leaves(&t, Cube::random(&mut rng), false, 6, &[1, 8]);
}

/// The same with slice moves (24 internal moves; depth 5 has about as many leaves as depth 6
/// without slices).
#[test]
fn search_reaches_every_canonical_leaf_once_with_slices() {
    let t = null_tables(SLICE);
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(98);
    check_leaves(&t, Cube::random(&mut rng), false, 5, &[1, 8]);
}

/// Canonical sequences against breadth-first search over all moves, for every length up to
/// `max_k`: every canonical sequence of length k reaches a position at distance exactly k, and
/// every such position is reached (so the move-order rules drop no shortest path: that is what
/// optimality needs). Up to `max_distinct` they also reach pairwise distinct positions (the
/// rules leave no duplicate); returns the number of duplicates at each length.
fn check_canonical(moves: MoveSet, max_k: usize, max_distinct: usize) -> Vec<u64> {
    use std::collections::HashMap;
    let n = moves.n_moves();
    let s0 = Cube::solved();
    let mut dist: HashMap<Cube, usize> = HashMap::from([(s0, 0)]);
    let mut front = vec![s0];
    let mut dups = vec![];
    for k in 1..=max_k {
        let mut next = vec![];
        for c in &front {
            for m in 0..n {
                let y = c.apply(m);
                if let std::collections::hash_map::Entry::Vacant(e) = dist.entry(y) {
                    e.insert(k);
                    next.push(y);
                }
            }
        }
        // all canonical sequences of length k, depth first; the first sequence to each end
        let mut ends: HashMap<Cube, Vec<usize>> = HashMap::with_capacity(next.len());
        let mut count = 0u64;
        let mut examples = vec![];
        let mut stack = vec![(s0, vec![])];
        while let Some((x, path)) = stack.pop() {
            if path.len() == k {
                assert_eq!(dist.get(&x), Some(&k), "{} length {k}", moves.name());
                count += 1;
                if let Some(first) = ends.get(&x) {
                    if examples.len() < 3 {
                        examples.push((first.clone(), path));
                    }
                } else {
                    ends.insert(x, path);
                }
                continue;
            }
            let last = path.last().copied().unwrap_or(usize::MAX);
            for m in 0..n {
                if moves.allowed_after(last, m) {
                    let mut p = path.clone();
                    p.push(m);
                    stack.push((x.apply(m), p));
                }
            }
        }
        assert_eq!(count, canonical_ways(moves, k).iter().sum::<u64>());
        assert_eq!(
            ends.len(),
            next.len(),
            "{} length {k}: missing",
            moves.name()
        );
        let d = count - ends.len() as u64;
        eprintln!(
            "{} length {k}: {count} canonical sequences, {} positions, {d} duplicates {examples:?}",
            moves.name(),
            ends.len()
        );
        if k <= max_distinct {
            assert_eq!(d, 0, "{} length {k}: duplicates", moves.name());
        }
        dups.push(d);
        front = next;
    }
    dups
}

/// Exact up to length 4 (with slices: 24, 440, 8064, 147792 sequences). Without slices, at
/// length 5, 944 of the 448336 positions have two canonical sequences: relations of length
/// 10 across three axes, e.g. BLw' B' BLw' BRw' BLw' = BRw' BLw' BRw' B' BRw' (internal
/// moves 13 1 13 11 13 and 11 13 11 1 11), which rules on consecutive moves cannot see. No
/// position is missed.
#[test]
fn canonical_sequences_are_distinct_positions() {
    check_canonical(MoveSet::FACE, 5, 4);
    check_canonical(SLICE, 4, 4);
}

/// Length 5 with slices (about 2.7 million sequences): `cargo test --release -- --ignored`.
#[test]
#[ignore]
fn canonical_sequences_with_slices_length5() {
    check_canonical(SLICE, 5, 4);
}

/// The same with the bound against the best solution active: the position is in `H_A` with
/// phase-2 length at least 8, so depth 0 already finds a solution of length >= 8, and with all
/// bounds 0 that bound cannot prune anything up to depth 6 (5 with slices).
#[test]
fn search_with_best_bound_reaches_every_canonical_leaf_once() {
    for (moves, depth) in [(MoveSet::FACE, 6), (SLICE, 5)] {
        let t = null_tables(moves);
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(7);
        let a_moves: Vec<usize> = (0..N_MOVES).filter(|&m| is_a_move(m)).collect();
        let cube = loop {
            let mut c = Cube::solved();
            for _ in 0..30 {
                c = c.apply(a_moves[rand::Rng::gen_range(&mut rng, 0..a_moves.len())]);
            }
            if t.p2.length(&c).is_some_and(|d| d >= 8) {
                break c;
            }
        };
        check_leaves(&t, cube, true, depth, &[1, 8]);
    }
}

/// Phase-2 membership must look at the whole position: a position whose edges and corners
/// form an H_A arrangement but whose centres do not is not in H_A.
#[test]
fn phase2_membership_checks_centres() {
    let p2 = Phase2::new();
    let s0 = Cube::solved();
    assert_eq!(p2.length(&s0), Some(0));
    let mut c = s0;
    c.xa.swap(0, 3); // two A-centres of different colours
    assert_eq!(p2.length(&c), None);
    let mut c = s0;
    c.xb.swap(0, 3);
    assert_eq!(p2.length(&c), None);
}

/// Depth 7 also covers the AE-first stage (used with at least 6 moves left). About 1.6e8
/// leaves, so it is not part of the default run: `cargo test --release -- --ignored`.
#[test]
#[ignore]
fn search_reaches_every_canonical_leaf_once_depth7() {
    let t = null_tables(MoveSet::FACE);
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(99);
    check_leaves(&t, Cube::random(&mut rng), false, 7, &[8]);
}
