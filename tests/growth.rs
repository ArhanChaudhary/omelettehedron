//! Analysis: number of distinct positions at each depth vs canonical move sequences
//! (shows how much the canonical-order rule leaves duplicated at short depths).
//! cargo test --release --test growth -- --ignored --nocapture

use fto_solver::cube::*;
use std::collections::HashSet;

#[test]
#[ignore]
fn sphere_vs_sequences() {
    let maxd = 6;
    // canonical sequences: count per last move
    let mut seq = vec![0u64; N_MOVES + 1];
    seq[N_MOVES] = 1;
    let mut seen: HashSet<Cube> = HashSet::new();
    seen.insert(Cube::solved());
    // frontier keeps (cube, last move) for every canonical sequence that reaches a new state
    // first (states only)
    let mut fr: Vec<Cube> = vec![Cube::solved()];
    // all canonical sequences (with multiplicity) whose position is new at that depth
    let mut fr_seq: Vec<(Cube, usize)> = vec![(Cube::solved(), usize::MAX)];
    for d in 1..=maxd {
        let mut nseq = vec![0u64; N_MOVES + 1];
        for last in 0..=N_MOVES {
            for m in 0..N_MOVES {
                let l = if last == N_MOVES { usize::MAX } else { last };
                if allowed_after(l, m) {
                    nseq[m] += seq[last];
                }
            }
        }
        seq = nseq;
        let total: u64 = seq.iter().sum();
        let mut nf = vec![];
        for c in &fr {
            for m in 0..N_MOVES {
                let n = c.apply(m);
                if seen.insert(n) {
                    nf.push(n);
                }
            }
        }
        // canonical sequences of length d that do not end at a new position (duplicates or
        // shorter positions), counted by extending the canonical sequences of length d-1
        // that ended at a new position at depth d-1
        let mut nfs = vec![];
        let mut dup = 0u64;
        let newset: HashSet<Cube> = nf.iter().copied().collect();
        let mut first: HashSet<Cube> = HashSet::new();
        if d <= 5 {
            for (c, last) in &fr_seq {
                for m in 0..N_MOVES {
                    if !allowed_after(*last, m) {
                        continue;
                    }
                    let n = c.apply(m);
                    if newset.contains(&n) && first.insert(n) {
                        nfs.push((n, m));
                    } else {
                        dup += 1;
                    }
                }
            }
        }
        let ext: u64 = if d <= 5 { (nfs.len() as u64) + dup } else { 0 };
        eprintln!(
            "depth {d}: positions {} canonical sequences {} ratio {:.4}; extending unique-prefix sequences: {} of which non-first {} ({:.4})",
            nf.len(),
            total,
            total as f64 / nf.len() as f64,
            ext,
            dup,
            dup as f64 / ext.max(1) as f64
        );
        fr = nf;
        fr_seq = nfs;
    }
}
