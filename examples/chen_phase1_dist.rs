//! Exact phase-1 distance distributions of the block split, counted over positions, for both
//! targets: `H'` (the hexagons whole; the chiral A-centre/edge table) and the relaxed target
//! (any of the 447,120 configurations of the phase-2 space).
//!
//! Usage: `cargo run --release --example chen_phase1_dist [samples]` with `FTO_TABLES` holding
//! `ae-sym288-chiral.prun` and `ae-sym288-chiral-relaxed.prun` (they are rebuilt in memory and
//! compared word for word). Prints:
//!
//! 1. the weighted histograms: each table entry stands for the A-centre/edge states of its class
//!    (one per A-centre colouring in its symmetry class), and every such state goes with the same
//!    number of whole positions, so this is the distribution over all positions;
//! 2. an exhaustive check of each target set: the states whose index is a target are exactly
//!    the targets' states, and (relaxed) every position whose A-centres and edges are such a
//!    state is a configuration of the phase-2 space (`space().id`);
//! 3. a check of the distances near the targets by brute force over all move sequences, and of
//!    the histograms on uniformly random positions.

use omelettehedron::chen::relaxed::relaxed_targets;
use omelettehedron::chen::space::space;
use omelettehedron::coord::{mset, triple_pieces, Coords, N_ECOSET, N_MSET};
use omelettehedron::cube::{allowed_after, Cube, N_MOVES};
use omelettehedron::fast::{fast, masks, Fc};
use omelettehedron::prune::Table2;
use omelettehedron::sym::{edge_syms, syms};
use omelettehedron::tables::{ae_targets, build_ae_weighted, table_dir, AeIndex, AeSym};
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use std::collections::{BTreeMap, HashSet};
use std::time::Instant;

/// An A-centre/edge state: (A-colouring rank, triple-labelling rank, 4 chirality bits).
type State = (u32, u32, u8);

fn main() {
    let samples: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(10_000_000);
    let dir = table_dir(None);
    let t0 = Instant::now();
    let ix = AeSym::new(true);
    // the weight of each A-centre class: the number of colourings in it
    let mut weight = vec![0u64; ix.cls.rep.len()];
    for &v in ix.cls.sym.iter() {
        weight[(v & 0xffff) as usize] += 1;
    }
    // the number of A-centre/edge states, from the coordinate's definition: A-centre colourings
    // 12! / (3!)^4, times edge cosets (12! / 2 even permutations) / 3^4 (each triple's three
    // edges turned in place)
    let fact12: u128 = (1..=12).product();
    let states = fact12 / 6u128.pow(4) * (fact12 / 2 / 81);
    assert_eq!(states, (N_MSET as u128) * (N_ECOSET as u128));
    assert_eq!(
        weight.iter().map(|&w| u128::from(w)).sum::<u128>(),
        N_MSET as u128
    );
    // positions per A-centre/edge state: corners 960, B-centre colourings, the 81 edge
    // arrangements of a coset
    let per_state: u128 = 960 * N_MSET as u128 * 81;
    println!(
        "A-centre/edge states {states}; positions per state {per_state}; positions {} ({:.4e}); \
         {} A-centre classes ({:.1} s)",
        states * per_state,
        (states * per_state) as f64,
        weight.len(),
        t0.elapsed().as_secs_f64()
    );

    let mut tables = vec![];
    for (name, file, targets) in [
        ("strict (H')", "ae-sym288-chiral.prun", ae_targets(&ix)),
        (
            "relaxed (any configuration)",
            "ae-sym288-chiral-relaxed.prun",
            relaxed_targets(&ix),
        ),
    ] {
        let t1 = Instant::now();
        let (tab, hist) =
            build_ae_weighted(&ix, &targets, ix.class_size(), &weight, "weighted", false);
        let total: u128 = hist.iter().sum();
        assert_eq!(total, states, "{name}: weighted total");
        // the table on disk is the same table
        let disk = Table2::load(&dir.join(file), ix.size()).expect("table file");
        assert!(
            disk.data.iter().eq(tab.data.iter()),
            "{name}: differs from {file}"
        );
        assert_eq!(disk.hist, tab.hist, "{name}: entry histogram");
        let mean = hist
            .iter()
            .enumerate()
            .map(|(d, &c)| d as f64 * c as f64)
            .sum::<f64>()
            / total as f64;
        let emean = tab
            .hist
            .iter()
            .enumerate()
            .map(|(d, &c)| d as f64 * c as f64)
            .sum::<f64>()
            / tab.hist.iter().sum::<u64>() as f64;
        println!(
            "\n## {name}: {} target indices; table identical to {file}; {:.1} s",
            targets.len(),
            t1.elapsed().as_secs_f64()
        );
        println!("depth  states (exact)            fraction        entries");
        for (d, &c) in hist.iter().enumerate() {
            println!(
                "{d:>5}  {c:>22}  {:>14.10}  {:>12}",
                c as f64 / total as f64,
                tab.hist.get(d).copied().unwrap_or(0)
            );
        }
        println!(
            "total  {total:>22}; max {}; mean over positions {mean:.6} (over entries {emean:.6})",
            hist.len() - 1
        );
        tables.push((name, tab, targets, hist));
    }

    // 2. the target sets, exhaustively
    let s = space();
    let s0 = Cube::solved();
    let config_states: HashSet<State> = s
        .cfgs
        .iter()
        .map(|&(xa, ep)| {
            let c = Coords::of(&Cube { ep, xa, ..s0 });
            (c.a, c.e, c.ebits)
        })
        .collect();
    let hprime_states: HashSet<State> = s
        .cfgs
        .iter()
        .enumerate()
        .filter(|&(i, _)| s.in_hprime(i as u32))
        .map(|(_, &(xa, ep))| {
            let c = Coords::of(&Cube { ep, xa, ..s0 });
            (c.a, c.e, c.ebits)
        })
        .collect();
    println!(
        "\nconfigurations {} in {} A-centre/edge states; H' {} configurations in {} states",
        s.cfgs.len(),
        config_states.len(),
        s.cfgs
            .iter()
            .enumerate()
            .filter(|&(i, _)| s.in_hprime(i as u32))
            .count(),
        hprime_states.len()
    );
    for ((name, _, targets, _), want) in tables.iter().zip([&hprime_states, &config_states]) {
        let u = preimage(&ix, targets);
        assert!(u.is_subset(want), "{name}: a target state is no target");
        assert_eq!(&u, want, "{name}: target states");
        println!(
            "{name}: the {} target indices hold exactly the {} target states",
            targets.len(),
            u.len()
        );
    }
    // every position whose A-centres and edges are a target state is a configuration: all 81
    // edge arrangements of each state's coset
    let tp = triple_pieces();
    let mut seen = 0;
    for &(a, e, bits) in &config_states {
        let xa = mset().unrank(a);
        for ep in coset(&tp, e, bits) {
            let c = Cube { ep, xa, ..s0 };
            assert!(s.id(&c).is_some(), "not a configuration: {c:?}");
            let co = Coords::of(&c);
            assert_eq!((co.a, co.e, co.ebits), (a, e, bits));
            seen += 1;
        }
    }
    assert_eq!(seen, s.cfgs.len());
    println!(
        "every one of the {seen} positions (up to corners and B-centres) whose A-centres and \
         edges are a relaxed target state is a configuration"
    );

    // 3a. brute force near the targets: no target within distance - 1, one at the distance
    // (positions a few random moves from a target)
    let f = fast();
    let hprime: Vec<([u8; 12], [u8; 12])> = s
        .cfgs
        .iter()
        .enumerate()
        .filter(|&(i, _)| s.in_hprime(i as u32))
        .map(|(_, &c)| c)
        .collect();
    for (k, (name, tab, targets, _)) in tables.iter().enumerate() {
        let starts = if k == 0 { &hprime } else { &s.cfgs };
        let is_target = |c: &Cube| s.id(c).is_some_and(|i| k == 1 || s.in_hprime(i));
        let mut by_d: BTreeMap<u8, usize> = BTreeMap::new();
        let found: Vec<(u8, bool)> = (0..4000u64)
            .into_par_iter()
            .filter_map(|seed| {
                let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(seed);
                let (xa, ep) = starts[rng.gen_range(0..starts.len())];
                let mut c = Cube {
                    ep,
                    xa,
                    ..Cube::random(&mut rng)
                };
                for _ in 0..rng.gen_range(0..8) {
                    c = c.apply(rng.gen_range(0..N_MOVES));
                }
                let d = walk(&ix, tab, targets, &f.of(&c))?;
                if d > 6 {
                    return None;
                }
                let none_closer = (0..d).all(|l| !reach(&c, l, usize::MAX, &is_target));
                let at_d = reach(&c, d, usize::MAX, &is_target);
                Some((d, none_closer && at_d))
            })
            .collect();
        for &(d, ok) in &found {
            assert!(ok, "{name}: table distance {d} is not the distance");
            *by_d.entry(d).or_default() += 1;
        }
        println!("{name}: brute force agrees on {by_d:?} positions (table distance: count)");
    }

    // 3b. uniformly random positions against the exact fractions
    for (name, tab, targets, hist) in &tables {
        let t1 = Instant::now();
        let chunks = samples.div_ceil(100_000);
        let counts = (0..chunks as u64)
            .into_par_iter()
            .map(|chunk| {
                let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(0xd157 + chunk);
                let mut h = vec![0u64; hist.len()];
                let n = (samples - chunk as usize * 100_000).min(100_000);
                for _ in 0..n {
                    let x = f.of(&Cube::random(&mut rng));
                    let d = walk(&ix, tab, targets, &x).expect("walk") as usize;
                    h[d] += 1;
                }
                h
            })
            .reduce(
                || vec![0u64; hist.len()],
                |a, b| a.iter().zip(&b).map(|(x, y)| x + y).collect(),
            );
        let total: u128 = hist.iter().sum();
        println!(
            "\n{name}: {samples} uniformly random positions ({:.1} s)",
            t1.elapsed().as_secs_f64()
        );
        println!("depth   observed    expected        z");
        let mut zmax = 0f64;
        let mut mean = 0f64;
        for (d, (&o, &c)) in counts.iter().zip(hist).enumerate() {
            let p = c as f64 / total as f64;
            let e = p * samples as f64;
            let z = if e > 0.0 {
                (o as f64 - e) / (e * (1.0 - p)).sqrt()
            } else {
                0.0
            };
            if e >= 5.0 {
                zmax = zmax.max(z.abs());
            }
            mean += d as f64 * o as f64;
            println!("{d:>5}  {o:>9}  {e:>11.1}  {z:>7.2}");
        }
        println!(
            "sample mean {:.5}; largest |z| where at least 5 expected: {zmax:.2}",
            mean / samples as f64
        );
    }
}

/// All A-centre/edge states whose index is one of `targets` (sorted, closed under the index's
/// equivalents): the orbits of the targets' representatives under the index's symmetries,
/// checked to land on targets and counted against the class sizes.
fn preimage(ix: &AeSym, targets: &[usize]) -> HashSet<State> {
    let (sy, es, ms) = (syms(), edge_syms(), mset());
    let mut u = HashSet::new();
    let mut expected = 0usize;
    for &i in targets {
        let (rep, e, bits) = ix.decode(i);
        let labels = ms.unrank(rep);
        let cls = (ix.cls.sym[rep as usize] & 0xffff) as usize;
        expected += ix
            .cls
            .sym
            .iter()
            .filter(|&&v| (v & 0xffff) as usize == cls)
            .count();
        for &(g, t) in &sy.pairs {
            let a = ms.rank(&sy.apply_a(g as usize, t as usize, &labels));
            let (e2, b2) = es.apply_e(g as usize, t as usize, e, bits);
            let x = Coords {
                c: 0,
                a,
                b: 0,
                e: e2,
                ebits: b2,
            };
            assert!(targets.binary_search(&ix.index(&x)).is_ok());
            u.insert((a, e2, b2));
        }
    }
    // each entry stands for one state per colouring of its class: the orbits hold them all
    assert_eq!(u.len(), expected, "preimage size");
    u
}

/// The 81 edge arrangements of an edge coset (triple labelling of rank `e`, chirality `bits`):
/// each triple's three pieces placed on its three positions in an order of the given parity.
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
    let mut out = vec![];
    for k in 0..81 {
        let mut ep = [0u8; 12];
        let mut kk = k;
        for t in 0..4 {
            // the three orders of the parity wanted (even ones first in `perms3`)
            let odd = usize::from(bits >> t & 1);
            let p = perms3[3 * odd + kk % 3];
            kk /= 3;
            for j in 0..3 {
                ep[pos[t][j]] = tp[t][p[j]];
            }
        }
        out.push(ep);
    }
    out
}

/// The table's distance of a position, by walking down the table until a target index.
fn walk(ix: &AeSym, tab: &Table2, targets: &[usize], x: &Fc) -> Option<u8> {
    let f = fast();
    let index = |y: &Fc| ix.index_fast(f, y, f.a_rank(y));
    let mut y = *x;
    for d in 0..=tab.hist.len() {
        let i = index(&y);
        if targets.binary_search(&i).is_ok() {
            return Some(d as u8);
        }
        let want = (tab.get(i) + 2) % 3;
        let em = masks(&y.e);
        y = (0..N_MOVES)
            .map(|m| f.apply_m(&y, &em, m))
            .find(|z| tab.get(index(z)) == want)?;
    }
    None
}

/// Whether some sequence of exactly `left` moves (canonical order) from `c` reaches a target.
fn reach(c: &Cube, left: u8, last: usize, is_target: &dyn Fn(&Cube) -> bool) -> bool {
    if left == 0 {
        return is_target(c);
    }
    (0..N_MOVES)
        .filter(|&m| allowed_after(false, last, m))
        .any(|m| reach(&c.apply(m), left - 1, m, is_target))
}
