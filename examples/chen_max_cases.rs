//! The phase-1 cases at the largest distance, for both targets of the block split: `H'` (14
//! moves) and any configuration of the phase-2 space (11 moves).
//!
//! Usage: `cargo run --release --example chen_max_cases > cases.jsonl` with `FTO_TABLES`
//! holding `ae-sym288-chiral.prun`, `ae-sym288-chiral-relaxed.prun` and the small tier's
//! tables. For each target: the farthest table entries, their A-centre/edge states (the orbits
//! of the entries' representatives), the 81 edge arrangements of each, and their classes under
//! the 12 whole-puzzle rotations (positions turned, colours kept). One case per class: the
//! smallest arrangement of the class with the corners and B-centres solved, an optimal phase-1
//! path (walked down the exact table), a two-phase solution of the whole position, and the
//! scramble that makes it (that solution inverted; it can leave the whole puzzle turned). One JSON line per case on stdout; counts on
//! stderr.

use omelettehedron::chen::relaxed::relaxed_targets;
use omelettehedron::chen::space::space;
use omelettehedron::coord::{mset, triple_pieces, Coords};
use omelettehedron::cube::{Cube, MoveSet, N_MOVES};
use omelettehedron::fast::{fast, masks, Fc};
use omelettehedron::geometry::{face_of_normal, mat_vec, rotation_group, transpose, FACES};
use omelettehedron::notation::{
    format_moves, invert_moves, scramble_to_cube, solution_to_physical, Frame, PhysMove,
};
use omelettehedron::pattern;
use omelettehedron::prune::Table2;
use omelettehedron::search::{solve, Limits, Report};
use omelettehedron::sym::{edge_syms, syms};
use omelettehedron::tables::{ae_farthest, ae_targets, table_dir, AeIndex, AeSym, Tables};
use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::time::Duration;

/// An A-centre/edge state: (A-colouring rank, triple-labelling rank, 4 chirality bits).
type State = (u32, u32, u8);
/// An arrangement of the A-centres and edges: (A-centre colouring, edge permutation).
type Arr = ([u8; 12], [u8; 12]);

fn main() {
    let dir = table_dir(None);
    let ix = AeSym::new(true);
    let strict = Table2::load(&dir.join("ae-sym288-chiral.prun"), ix.size()).expect("table");
    let relaxed =
        Table2::load(&dir.join("ae-sym288-chiral-relaxed.prun"), ix.size()).expect("table");
    let small = Tables::load(&dir, AeSym::new(false), MoveSet::FACE, false).expect("tables");
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(8)
        .build()
        .expect("pool");
    let rots = rotation_group();
    for (kind, tab, targets, want_d) in [
        ("strict", &strict, ae_targets(&ix), 14u8),
        ("relaxed", &relaxed, relaxed_targets(&ix), 11u8),
    ] {
        let (d, far) = ae_farthest(&ix, &targets);
        assert_eq!(d, want_d, "{kind}: largest distance");
        let states = preimage(&ix, &far);
        // every arrangement of each state
        let tp = triple_pieces();
        let arrs: BTreeSet<Arr> = states
            .iter()
            .flat_map(|&(a, e, bits)| {
                let xa = mset().unrank(a);
                coset(&tp, e, bits).into_iter().map(move |ep| (xa, ep))
            })
            .collect();
        // the classes under the 12 rotations: each arrangement's turned versions are
        // arrangements of the set too, at the same distance
        let mut classes: BTreeMap<Arr, Vec<Arr>> = BTreeMap::new();
        for &x in &arrs {
            let turned: Vec<Arr> = rots.iter().map(|m| turn(&x, m)).collect();
            for t in &turned {
                assert!(
                    arrs.contains(t),
                    "{kind}: a turned arrangement is not in the set"
                );
            }
            let distinct: HashSet<&Arr> = turned.iter().collect();
            assert_eq!(
                distinct.len(),
                12,
                "{kind}: an arrangement fixed by a rotation"
            );
            let key = *turned.iter().min().expect("12 rotations");
            classes.entry(key).or_default().push(x);
        }
        assert!(classes.values().all(|v| v.len() == 12));
        eprintln!(
            "{kind}: distance {d}: {} table entries, {} A-centre/edge states, {} arrangements, \
             {} classes under the 12 rotations",
            far.len(),
            states.len(),
            arrs.len(),
            classes.len()
        );
        // the smallest member of each class
        let reps: Vec<Arr> = classes.keys().copied().collect();
        report(kind, d, &reps, tab, &targets, &small, &pool, &rots);
    }
}

/// One JSON line per case: the case's position (A-centres and edges `reps[n]`, the corners and
/// B-centres solved), an optimal phase-1 path, a two-phase solution and the scramble.
#[allow(clippy::too_many_arguments)]
fn report(
    kind: &str,
    d: u8,
    reps: &[Arr],
    tab: &Table2,
    targets: &[usize],
    small: &Tables<AeSym>,
    pool: &rayon::ThreadPool,
    rots: &[omelettehedron::geometry::Mat],
) {
    let want_d = d;
    let s = space();
    let ix = AeSym::new(true);
    let lines: Vec<String> = pool.install(|| {
        reps.iter()
            .enumerate()
            .map(|(n, &(xa, ep))| {
                let cube = Cube {
                    ep,
                    xa,
                    ..Cube::solved()
                };
                let x = fast().of(&cube);
                // an optimal phase-1 path: walk down the exact table
                let path = descend(&ix, tab, targets, &x).expect("walk");
                assert_eq!(path.len(), usize::from(want_d));
                let end = cube.apply_seq(&path);
                let cfg = s.id(&end).expect("a configuration");
                let target = if s.in_hprime(cfg) {
                    "H' (every hexagon whole)"
                } else {
                    "a configuration (not H')"
                };
                assert!(kind == "relaxed" || s.in_hprime(cfg));
                // a solution of the whole position, and the scramble that makes it
                // (stop at a first solution of at most 24 moves, or at the time limit;
                // longer limits until one is found)
                let sol = [500u64, 2000, 8000, 30000, 120_000]
                    .iter()
                    .find_map(|&ms| {
                        let limits = Limits {
                            time: Some(Duration::from_millis(ms)),
                            max_length: Some(24),
                            ..Limits::default()
                        };
                        solve(small, cube, Frame::identity(), &limits, Report::Quiet).solution
                    })
                    .expect("a solution");
                // The solution ends in solved turned as a whole by some rotation r; its
                // inverse with every face renamed by r (or r^-1) makes the position from
                // solved, ending held the standard way (so the corners and B-centres show
                // solved in the pattern).
                let inv = invert_moves(&sol);
                let candidates: Vec<Vec<PhysMove>> = rots
                    .iter()
                    .flat_map(|m| [*m, transpose(m)])
                    .map(|m| {
                        inv.iter()
                            .map(|pm| PhysMove {
                                face: face_of_normal(mat_vec(&m, FACES[pm.face])),
                                ..*pm
                            })
                            .collect::<Vec<PhysMove>>()
                    })
                    .filter(|t| scramble_to_cube(t).0 == cube)
                    .collect();
                let scramble = candidates
                    .iter()
                    .find(|t| scramble_to_cube(t).1 == Frame::identity())
                    .or(candidates.first())
                    .cloned()
                    .expect("a scramble");
                let (_, frame) = scramble_to_cube(&scramble);
                let held = frame == Frame::identity();
                // the phase-1 path as turns of the puzzle as the scramble leaves it
                let phase1 = solution_to_physical(frame, &path, false);
                // the same solution as turns of the puzzle as the scramble leaves it
                let sol = invert_moves(&scramble);
                let pat = pattern::format(&pattern::of_moves(&scramble));
                serde_json::json!({
                    "kind": kind,
                    "case": n + 1,
                    "distance": d,
                    "pattern": pat,
                    "scramble": format_moves(&scramble),
                    "phase1": format_moves(&phase1),
                    "target": target,
                    "solution": format_moves(&sol),
                    "solution_length": sol.len(),
                    "standard_hold": held,
                    "xa": xa,
                    "ep": ep,
                })
                .to_string()
            })
            .collect()
    });
    for l in lines {
        println!("{l}");
    }
}

/// The arrangement with the whole puzzle turned by `m` (stickers moved, colours kept).
fn turn(x: &Arr, m: &omelettehedron::geometry::Mat) -> Arr {
    let sm = &omelettehedron::cube::defs().sm;
    let st = Cube {
        ep: x.1,
        xa: x.0,
        ..Cube::solved()
    }
    .to_stickers();
    let perm = sm.whole(m);
    let mut n = [0u8; 72];
    for i in 0..72 {
        n[perm[i]] = st[i];
    }
    let c = Cube::from_stickers(&n).expect("a turned position");
    (c.xa, c.ep)
}

/// All A-centre/edge states whose index is one of `entries`: the orbits of the entries'
/// representatives under the index's symmetries (each landing on one of the entries).
fn preimage(ix: &AeSym, entries: &[usize]) -> HashSet<State> {
    let (sy, es, ms) = (syms(), edge_syms(), mset());
    let mut u = HashSet::new();
    for &i in entries {
        let (rep, e, bits) = ix.decode(i);
        let labels = ms.unrank(rep);
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
            assert!(entries.contains(&ix.index(&x)));
            u.insert((a, e2, b2));
        }
    }
    u
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

/// An optimal path to a target: each move to a neighbour one closer by the table.
fn descend(ix: &AeSym, tab: &Table2, targets: &[usize], x: &Fc) -> Option<Vec<usize>> {
    let f = fast();
    let index = |y: &Fc| ix.index_fast(f, y, f.a_rank(y));
    let mut y = *x;
    let mut path = vec![];
    for _ in 0..=tab.hist.len() {
        let i = index(&y);
        if targets.binary_search(&i).is_ok() {
            return Some(path);
        }
        let want = (tab.get(i) + 2) % 3;
        let em = masks(&y.e);
        let m = (0..N_MOVES).find(|&m| tab.get(index(&f.apply_m(&y, &em, m))) == want)?;
        y = f.apply_m(&y, &em, m);
        path.push(m);
    }
    None
}
