//! Can a position be 11 moves (the largest distance) from the relaxed phase-1 target of the
//! block split on both sides at once? Side A splits the hexagons of U, F, BR and BL; side B
//! those of R, L, B and D, which is side A of the position conjugated by `SWAP_AB` (as the
//! solver's side 1).
//!
//! Usage: `cargo run --release --example chen_both_sides [samples [uniform]]` (positions per
//! distance-11 state, default 1,000; uniformly random positions, default 4,000,000) with
//! `FTO_TABLES` holding `ae-sym288-chiral-relaxed.prun`. With `--patterns FILE` instead: every
//! position with both sides at 11 (one per rotation class: the front corner in place), built
//! here, against the patterns of FILE (one per line): the same positions, each with both sides
//! 11 moves away by the table.
//!
//! The side-A distance depends only on the A-centres and the edges, and the side-B distance
//! only on the B-centres and the edges. The four parts of a position (corners, A-centres,
//! B-centres, edges) are independent ([`Cube::random`]). The conjugate's edges depend only on
//! the edges (map `g`), and its A-centres only on the B-centres (map `f`, a bijection of the
//! colourings); both are checked here. So for an edge arrangement `ep`, the number of
//! A-colourings with side A at least `d` moves away is `n_d[coset(ep)]`. The number of
//! B-colourings with side B at least `d` moves away is `n_d[coset(g(ep))]`, where `n_d` counts
//! the states with distance at least `d` in each edge coset. Positions with both sides at least
//! `d` moves away: 960 corner states times the sum over `ep` of the products. Prints:
//!
//! 1. the edge arrangements of the distance-11 states (`EA`), their images under `g`, and the
//!    intersection of `EA` and `g^-1(EA)` (the edges of a position with both sides at 11);
//! 2. the positions with both sides at least 10 moves away, split by the two distances;
//! 3. cross-checks by the table on sampled positions: positions built from each distance-11
//!    state with random corners, B-centres and edge arrangement in its coset (side-B
//!    distances), and uniformly random positions (how often both sides are at least 10 moves
//!    away).

use omelettehedron::chen::relaxed::relaxed_targets;
use omelettehedron::coord::{
    corners, edge_chirality, mset, triple_pieces, Coords, N_ECOSET, N_MSET,
};
use omelettehedron::cube::{perm_parity, Cube, N_MOVES, SWAP_AB};
use omelettehedron::fast::{fast, masks, Fc};
use omelettehedron::geometry::transpose;
use omelettehedron::pattern;
use omelettehedron::prune::Table2;
use omelettehedron::sym::{edge_syms, syms};
use omelettehedron::tables::{ae_tail, table_dir, AeIndex, AeSym};
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use std::collections::{BTreeMap, HashSet};
use std::sync::atomic::{AtomicU32, Ordering};
use std::time::Instant;

/// An A-centre/edge state: (A-colouring rank, triple-labelling rank, 4 chirality bits).
type State = (u32, u32, u8);

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let patterns_file = args
        .iter()
        .position(|a| a == "--patterns")
        .and_then(|i| args.get(i + 1))
        .cloned();
    let numbers: Vec<usize> = args[1..].iter().filter_map(|a| a.parse().ok()).collect();
    let samples = numbers.first().copied().unwrap_or(1000);
    let n_uniform = numbers.get(1).copied().unwrap_or(4_000_000);
    let t0 = Instant::now();
    let dir = table_dir(None);
    let ix = AeSym::new(true);
    let tab = Table2::load(&dir.join("ae-sym288-chiral-relaxed.prun"), ix.size()).expect("table");
    let targets = relaxed_targets(&ix);
    let ms = mset();
    let conj = |c: &Cube| {
        c.conjugate(&SWAP_AB)
            .expect("SWAP_AB keeps the front corner")
    };
    let unconj = |c: &Cube| {
        c.conjugate(&transpose(&SWAP_AB))
            .expect("SWAP_AB^-1 keeps the front corner")
    };
    let dist = |x: &Cube| distance(&ix, &tab, &targets, &fast().of(x)).expect("a distance");
    let g = |ep: &[u8; 12]| {
        conj(&Cube {
            ep: *ep,
            ..Cube::solved()
        })
        .ep
    };

    // the conjugate's edges depend only on the edges, its A-centres only on the B-centres
    let mut rng = rand::rngs::StdRng::seed_from_u64(1);
    for _ in 0..100_000 {
        let x = Cube::random(&mut rng);
        let y = conj(&x);
        assert_eq!(y.ep, g(&x.ep));
        assert_eq!(
            y.xa,
            conj(&Cube {
                xb: x.xb,
                ..Cube::solved()
            })
            .xa
        );
        assert_eq!(unconj(&y), x);
    }
    let images: HashSet<u32> = (0..N_MSET as u32)
        .into_par_iter()
        .map(|b| {
            let x = Cube {
                xb: ms.unrank(b),
                ..Cube::solved()
            };
            ms.rank(&conj(&x).xa)
        })
        .collect();
    assert_eq!(images.len(), N_MSET, "f is a bijection");
    println!(
        "conjugate by SWAP_AB: edges a function g of the edges, A-centres a bijection f of the \
         B-centres (100,000 random positions; all {N_MSET} colourings)"
    );

    // the entries at least 10 moves away, then 11
    let tail = ae_tail(&ix, &targets, 10);
    assert_eq!(tail.len(), 2, "the largest distance is 11");
    println!(
        "entries: {} at least 10 moves away, {} at 11 ({:.0} s)",
        tail[0].len(),
        tail[1].len(),
        t0.elapsed().as_secs_f64()
    );
    let (n11, d11) = per_coset(&ix, &tail[1]);
    if let Some(path) = patterns_file {
        let tp = triple_pieces();
        let ea: HashSet<[u8; 12]> = d11
            .iter()
            .flat_map(|&(_, e, bits)| coset(&tp, e, bits))
            .collect();
        let mut both: Vec<[u8; 12]> = ea
            .iter()
            .filter(|ep| ea.contains(&g(ep)))
            .copied()
            .collect();
        both.sort_unstable();
        // every position with both sides at 11: the A-colourings of the distance-11 states of
        // the edges' coset, the B-colourings whose conjugates are those of the states of the
        // conjugate's coset, any corners (the front corner in place)
        let id = |s: &State| s.1 as usize * 8 + usize::from(s.2 & 7);
        let mut mine: HashSet<Cube> = HashSet::new();
        for ep in &both {
            let (ca, cb) = (coset_id(ep), coset_id(&g(ep)));
            for sa in d11.iter().filter(|s| id(s) == ca) {
                for sb in d11.iter().filter(|s| id(s) == cb) {
                    let xb = unconj(&Cube {
                        xa: ms.unrank(sb.0),
                        ..Cube::solved()
                    })
                    .xb;
                    for &(cp, co) in &corners().id_to_state {
                        mine.insert(Cube {
                            cp,
                            co,
                            ep: *ep,
                            xa: ms.unrank(sa.0),
                            xb,
                        });
                    }
                }
            }
        }
        let text = std::fs::read_to_string(&path).expect("the patterns file");
        let lines: Vec<&str> = text.lines().filter(|l| !l.trim().is_empty()).collect();
        let theirs: Vec<Cube> = lines
            .par_iter()
            .map(|l| {
                let c = pattern::parse(l)
                    .and_then(|p| pattern::position(&p))
                    .expect("a valid pattern");
                assert_eq!((dist(&c), dist(&conj(&c))), (11, 11), "{l}");
                c
            })
            .collect();
        let theirs_set: HashSet<Cube> = theirs.iter().copied().collect();
        println!(
            "positions with both sides at 11 built here: {}; patterns in {path}: {} ({} distinct), \
             every one both sides 11 moves away by the table; the same set: {}",
            mine.len(),
            lines.len(),
            theirs_set.len(),
            theirs_set == mine
        );
        assert_eq!(mine.len(), 138_240);
        assert_eq!(theirs_set, mine);
        return;
    }
    let (n10, _) = per_coset(&ix, &tail[0]);
    let sum = |n: &[u32]| n.iter().map(|&v| u64::from(v)).sum::<u64>();
    // the states per distance from the exact distribution (chen_phase1_dist)
    assert_eq!(sum(&n10), 6_617_546_556 + 96);
    assert_eq!(sum(&n11), 96);
    assert_eq!(d11.len(), 96);
    println!(
        "states: {} at least 10 moves away in {} edge cosets; 96 at 11 in {} edge cosets ({:.0} s)",
        sum(&n10),
        n10.iter().filter(|&&v| v > 0).count(),
        n11.iter().filter(|&&v| v > 0).count(),
        t0.elapsed().as_secs_f64()
    );

    // 1. the edges of the distance-11 states, and of a position with both sides at 11
    let tp = triple_pieces();
    let ea: HashSet<[u8; 12]> = d11
        .iter()
        .flat_map(|&(_, e, bits)| coset(&tp, e, bits))
        .collect();
    let eb: HashSet<[u8; 12]> = ea
        .iter()
        .map(|ep| {
            unconj(&Cube {
                ep: *ep,
                ..Cube::solved()
            })
            .ep
        })
        .collect();
    let mut both: Vec<[u8; 12]> = ea.intersection(&eb).copied().collect();
    both.sort_unstable();
    let both_g = ea.iter().filter(|ep| ea.contains(&g(ep))).count();
    assert_eq!(both.len(), both_g);
    println!(
        "EA: {} edge arrangements (fixed frame); EB = g^-1(EA): {}; EA and EB in common: {}",
        ea.len(),
        eb.len(),
        both.len()
    );
    let both_positions: u128 = both
        .iter()
        .map(|ep| 960 * u128::from(n11[coset_id(ep)]) * u128::from(n11[coset_id(&g(ep))]))
        .sum();
    println!("positions with both sides 11 moves away: {both_positions}");
    for ep in both.iter().take(3) {
        let a = d11
            .iter()
            .find(|s| coset(&tp, s.1, s.2).contains(ep))
            .expect("a state");
        let ga = g(ep);
        let b = d11
            .iter()
            .find(|s| coset(&tp, s.1, s.2).contains(&ga))
            .expect("a state");
        let x = Cube {
            ep: *ep,
            xa: ms.unrank(a.0),
            xb: unconj(&Cube {
                xa: ms.unrank(b.0),
                ..Cube::solved()
            })
            .xb,
            ..Cube::solved()
        };
        // both distances by the table
        assert_eq!((dist(&x), dist(&conj(&x))), (11, 11));
        println!(
            "  example (both sides 11 by the table): {}",
            pattern::format(&pattern::of_cube(&x))
        );
    }

    // 2. positions with both sides at least 10 moves away
    let mut valid = vec![];
    for e in 0..N_MSET as u32 {
        for bits in 0..16u8 {
            let c = coset(&tp, e, bits);
            if perm_parity(&c[0]) == 0 {
                valid.push((e, bits));
            }
        }
    }
    assert_eq!(valid.len(), N_ECOSET);
    // [d_A >= 10 and d_B >= 10, d_A = 11 and d_B >= 10, d_A >= 10 and d_B = 11, both 11]
    let pairs: [u128; 4] = valid
        .par_iter()
        .map(|&(e, bits)| {
            let eps = coset(&tp, e, bits);
            let c = coset_id(&eps[0]);
            assert_eq!(c, e as usize * 8 + usize::from(bits & 7));
            let mut s = [0u128; 4];
            if n10[c] == 0 {
                return s;
            }
            for ep in &eps {
                let cb = coset_id(&g(ep));
                let (a10, a11) = (u128::from(n10[c]), u128::from(n11[c]));
                let (b10, b11) = (u128::from(n10[cb]), u128::from(n11[cb]));
                s[0] += a10 * b10;
                s[1] += a11 * b10;
                s[2] += a10 * b11;
                s[3] += a11 * b11;
            }
            s
        })
        .reduce(
            || [0; 4],
            |a, b| [a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3]],
        );
    let total: u128 = 960 * (N_MSET as u128).pow(2) * (N_ECOSET as u128 * 81);
    let p = |k: usize| 960 * pairs[k];
    println!(
        "positions (of {total}): both sides at least 10 moves away {} ({:.4e} of all); side A 11 \
         and side B at least 10: {}; side A at least 10 and side B 11: {}; both 11: {} ({:.0} s)",
        p(0),
        p(0) as f64 / total as f64,
        p(1),
        p(2),
        p(3),
        t0.elapsed().as_secs_f64()
    );
    println!(
        "so min(side A, side B) is at most {} for every position, and equals 10 for {} positions",
        if p(3) == 0 { 10 } else { 11 },
        p(0) - p(3)
    );

    // 3. cross-checks by the table
    let hist_b: BTreeMap<usize, usize> = d11
        .par_iter()
        .enumerate()
        .map(|(k, &(a, e, bits))| {
            let mut rng = rand::rngs::StdRng::seed_from_u64(1000 + k as u64);
            let eps = coset(&tp, e, bits);
            let mut h = BTreeMap::new();
            for _ in 0..samples {
                let x = Cube {
                    xa: ms.unrank(a),
                    ep: eps[rng.gen_range(0..81)],
                    ..Cube::random(&mut rng)
                };
                assert_eq!(dist(&x), 11);
                *h.entry(dist(&conj(&x))).or_insert(0) += 1;
            }
            h
        })
        .reduce(BTreeMap::new, |mut a, b| {
            for (k, v) in b {
                *a.entry(k).or_insert(0) += v;
            }
            a
        });
    println!(
        "side B of {} positions with side A 11 ({samples} per state: random corners, B-centres, \
         edge arrangement of the coset): {hist_b:?}",
        96 * samples
    );
    let hist_ab: BTreeMap<(usize, usize), usize> = (0..n_uniform)
        .into_par_iter()
        .fold(BTreeMap::new, |mut h, i| {
            let mut rng = rand::rngs::StdRng::seed_from_u64(1 << 40 | i as u64);
            let x = Cube::random(&mut rng);
            let (da, db) = (dist(&x), dist(&conj(&x)));
            if da >= 10 && db >= 10 {
                *h.entry((da, db)).or_insert(0) += 1;
            }
            h
        })
        .reduce(BTreeMap::new, |mut a, b| {
            for (k, v) in b {
                *a.entry(k).or_insert(0) += v;
            }
            a
        });
    let hits: usize = hist_ab.values().sum();
    println!(
        "uniformly random positions: {hits} of {n_uniform} with both sides at least 10 moves \
         away ({hist_ab:?}); expected {:.1} ({:.0} s)",
        n_uniform as f64 * p(0) as f64 / total as f64,
        t0.elapsed().as_secs_f64()
    );
}

/// The states of `entries` counted per edge coset (`e * 8 + chirality & 7`, as in the index),
/// and the states themselves if there are few. Each entry holds one state per A-colouring of
/// its class: the images of the entry's representative state under the 288 symmetries that
/// index back to the entry, one per colouring.
fn per_coset(ix: &AeSym, entries: &[usize]) -> (Vec<u32>, Vec<State>) {
    let (sy, es, ms) = (syms(), edge_syms(), mset());
    let n: Vec<AtomicU32> = (0..N_ECOSET).map(|_| AtomicU32::new(0)).collect();
    let keep = entries.len() <= 1000;
    let states: Vec<State> = entries
        .par_iter()
        .flat_map_iter(|&i| {
            let (rep, e, bits) = ix.decode(i);
            let labels = ms.unrank(rep);
            let mut seen: Vec<State> = Vec::with_capacity(288);
            for &(g, t) in &sy.pairs {
                let (g, t) = (usize::from(g), usize::from(t));
                let a = ms.rank(&sy.apply_a(g, t, &labels));
                let (e2, b2) = es.apply_e(g, t, e, bits);
                let x = Coords {
                    c: 0,
                    a,
                    b: 0,
                    e: e2,
                    ebits: b2,
                };
                if ix.index(&x) == i {
                    seen.push((a, e2, b2));
                }
            }
            // (symmetries fixing the representative state give the same state again)
            seen.sort_unstable();
            seen.dedup();
            for &(_, e2, b2) in &seen {
                n[e2 as usize * 8 + usize::from(b2 & 7)].fetch_add(1, Ordering::Relaxed);
            }
            if keep {
                seen
            } else {
                vec![]
            }
        })
        .collect();
    (n.into_iter().map(AtomicU32::into_inner).collect(), states)
}

/// The edge coset of an arrangement, as in the index.
fn coset_id(ep: &[u8; 12]) -> usize {
    let c = Cube {
        ep: *ep,
        ..Cube::solved()
    };
    mset().rank(&c.etri()) as usize * 8 + usize::from(edge_chirality(ep) & 7)
}

/// The 81 edge arrangements of an edge coset (triple labelling of rank `e`, chirality `bits`).
fn coset(tp: &[[u8; 3]; 4], e: u32, bits: u8) -> Vec<[u8; 12]> {
    let labels = mset().unrank(e);
    let pos: Vec<Vec<usize>> = (0..4)
        .map(|t| (0..12).filter(|&p| labels[p] == t as u8).collect())
        .collect();
    let perms3 = [
        [0, 1, 2],
        [1, 2, 0],
        [2, 0, 1],
        [0, 2, 1],
        [1, 0, 2],
        [2, 1, 0],
    ];
    (0..81)
        .map(|k| {
            let mut ep = [0u8; 12];
            let mut kk = k;
            for t in 0..4 {
                let p = perms3[3 * usize::from(bits >> t & 1) + kk % 3];
                kk /= 3;
                for j in 0..3 {
                    ep[pos[t][j]] = tp[t][p[j]];
                }
            }
            ep
        })
        .collect()
}

/// The exact distance to a target: a walk down the table, each move to a neighbour one closer.
fn distance(ix: &AeSym, tab: &Table2, targets: &[usize], x: &Fc) -> Option<usize> {
    let f = fast();
    let index = |y: &Fc| ix.index_fast(f, y, f.a_rank(y));
    let mut y = *x;
    for d in 0..=tab.hist.len() {
        let i = index(&y);
        if targets.binary_search(&i).is_ok() {
            return Some(d);
        }
        let want = (tab.get(i) + 2) % 3;
        let em = masks(&y.e);
        let m = (0..N_MOVES).find(|&m| tab.get(index(&f.apply_m(&y, &em, m))) == want)?;
        y = f.apply_m(&y, &em, m);
    }
    None
}
