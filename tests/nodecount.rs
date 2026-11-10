//! Search coverage check.
//!
//! With every pruning table zeroed (heuristic = 0 everywhere), the search must reach each
//! phase-1 leaf exactly once. A leaf of the depth-d iteration is a move sequence of length d
//! that follows the move-order rules and ends with a move outside H_A, so their number is
//! counted exactly by a transfer matrix over the last move. A mismatch means that the
//! parallel split, the task scheduling or the batched child evaluation drops or repeats part of
//! the tree. The tables themselves are not involved.

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

/// Small-tier tables with all pruning values 0.
fn null_tables() -> Tables<AeSym> {
    let ae_ix = AeSym::new(false);
    let n = ae_ix.size();
    Tables {
        moves: MoveSet::FACE,
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
    }
}

/// Move sequences of length `depth` allowed by the move-order rules whose last move is
/// outside H_A (the empty sequence for depth 0).
fn canonical_leaves(depth: usize) -> u64 {
    if depth == 0 {
        return 1;
    }
    // ways[m]: allowed sequences of the current length ending with move m
    let mut ways = [1u64; N_MOVES];
    for _ in 1..depth {
        let mut next = [0u64; N_MOVES];
        for (last, &w) in ways.iter().enumerate() {
            for (m, n) in next.iter_mut().enumerate() {
                if allowed_after(false, last, m) {
                    *n += w;
                }
            }
        }
        ways = next;
    }
    (0..N_MOVES)
        .filter(|&m| !is_a_move(m))
        .map(|m| ways[m])
        .sum()
}

/// Runs the search on `cube` up to depth `d` with each thread count and checks the leaf count.
fn check_leaves(t: &Tables<AeSym>, cube: Cube, expect_solution: bool, d: usize, threads: &[usize]) {
    // two variants (targets H_A and H_B), each searched at every depth 0..=D
    let expected: u64 = 2 * (0..=d).map(canonical_leaves).sum::<u64>();
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
        assert_eq!(o.leaves, expected, "{threads} threads");
    }
}

#[test]
fn search_reaches_every_canonical_leaf_once() {
    let t = null_tables();
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(99);
    // a random position is about 15 moves from H_A, so no leaf is a solution and the search
    // never prunes against a best solution
    check_leaves(&t, Cube::random(&mut rng), false, 6, &[1, 8]);
}

/// The same with the bound against the best solution active: the position is in `H_A` with
/// phase-2 length at least 8, so depth 0 already finds a solution of length >= 8, and with all
/// bounds 0 that bound cannot prune anything up to depth 6.
#[test]
fn search_with_best_bound_reaches_every_canonical_leaf_once() {
    let t = null_tables();
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
    check_leaves(&t, cube, true, 6, &[1, 8]);
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
    let t = null_tables();
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(99);
    check_leaves(&t, Cube::random(&mut rng), false, 7, &[8]);
}
