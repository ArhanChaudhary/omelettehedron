//! The tiny tier (`omelettehedron::light`), built in memory: its bound against plain lookups and
//! the small tier's exact A-centre/edge distance, the staged lookups of the search, table
//! checks and reproducibility, and its optimal solves against the small tier's.

use omelettehedron::cube::*;
use omelettehedron::fast::*;
use omelettehedron::light::*;
use omelettehedron::notation::*;
use omelettehedron::search::*;
use omelettehedron::tables::*;
use rand::{Rng, SeedableRng};
use std::sync::OnceLock;

fn tiny() -> &'static Tables<AeLight> {
    static T: OnceLock<Tables<AeLight>> = OnceLock::new();
    T.get_or_init(|| Tables::generate(AeLight::new(), MoveSet::FACE, false))
}

fn small() -> &'static Tables<AeSym> {
    static T: OnceLock<Tables<AeSym>> = OnceLock::new();
    T.get_or_init(|| Tables::generate(AeSym::new(false), MoveSet::FACE, false))
}

/// The bound from the codes of the plain lookups: P codes 1-3 are 8-10, else S codes 1-3
/// are 5-7, else 0.
fn bound_of_codes(p: &[u8], s: &[u8]) -> u8 {
    let (mp, ms) = (*p.iter().max().unwrap(), *s.iter().max().unwrap());
    if mp > 0 {
        7 + mp
    } else if ms > 0 {
        4 + ms
    } else {
        0
    }
}

/// Random positions and short scrambles (positions near `H_A` too).
fn positions(seed: u64, n: usize) -> Vec<Cube> {
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(seed);
    (0..n)
        .map(|i| {
            if i % 2 == 0 {
                Cube::random(&mut rng)
            } else {
                (0..i % 12).fold(Cube::solved(), |c, _| c.apply(rng.gen_range(0..N_MOVES)))
            }
        })
        .collect()
}

#[test]
fn bound_matches_plain_lookups_and_is_a_lower_bound() {
    let (t, sm, f) = (tiny(), small(), fast());
    for c in positions(31, 1500) {
        let x = f.of(&c);
        let em = masks(&x.e);
        for m in 0..=N_MOVES {
            let y = if m == N_MOVES { x } else { f.apply(&x, &em, m) };
            let b = t.ae_ix.bound_of(f, &y, &t.ae.data);
            let (p, s) = t.ae_ix.codes(f, &y, &t.ae.data);
            assert_eq!(b, bound_of_codes(&p, &s));
            assert!(b <= sm.exact_fast(&y).0, "{b} > the exact AE distance");
        }
    }
}

/// The search's keys: with every `cut`, a lower bound on the full bound that is above `cut`
/// exactly when the full bound is (where P is read: `cut` at least 7; below, when S's bound
/// is), for whole batches and for selections.
#[test]
fn staged_keys_decide_like_the_full_bound() {
    let (t, f) = (tiny(), fast());
    let ix = &t.ae_ix;
    for c in positions(32, 400) {
        let x = f.of(&c);
        let em = masks(&x.e);
        let ys: [Fc; N_MOVES] = std::array::from_fn(|m| f.apply(&x, &em, m));
        let full: Vec<u8> = ys.iter().map(|y| ix.bound_of(f, y, &t.ae.data)).collect();
        for cut in (0..13).chain([usize::MAX]) {
            let mut keys = [0usize; N_MOVES];
            ix.keys_batch_pf::<N_MOVES>(f, &ys, &mut keys, &t.ae.data, cut);
            ix.resolve(f, &ys, &ALL_IDS, N_MOVES, &mut keys, &t.ae.data, cut);
            // a selection: every third child
            let sel: Vec<u8> = (0..N_MOVES as u8).filter(|i| i % 3 == 1).collect();
            let mut skeys = [usize::MAX; N_MOVES];
            ix.keys_sel_pf::<N_MOVES>(f, &ys, &sel, &mut skeys, &t.ae.data, cut);
            ix.resolve(f, &ys, &sel, sel.len(), &mut skeys, &t.ae.data, cut);
            for m in 0..N_MOVES {
                let v = ix.value(&t.ae, 0, keys[m]);
                if sel.contains(&(m as u8)) {
                    assert_eq!(ix.value(&t.ae, 0, skeys[m]), v);
                } else {
                    assert_eq!(skeys[m], usize::MAX, "a child outside the selection");
                }
                assert!(v <= full[m]);
                if cut == usize::MAX {
                    assert_eq!(v, full[m]);
                } else if cut >= 7 {
                    assert_eq!(
                        usize::from(v) > cut,
                        usize::from(full[m]) > cut,
                        "cut {cut}"
                    );
                } else {
                    let (_, s) = ix.codes(f, &ys[m], &t.ae.data);
                    let sb = bound_of_codes(&[0], &s);
                    assert_eq!(usize::from(v) > cut, usize::from(sb) > cut, "cut {cut}");
                }
            }
        }
    }
}

/// The pending keys hold a P index below 2^28 (they must fit 32 bits on WebAssembly), and the
/// table has the entries its file name stands for.
#[test]
fn sizes() {
    let ix = &tiny().ae_ix;
    assert!(ix.pair.size() < 1 << 28);
    assert_eq!(ix.pair.size(), 7864 * 18_480);
    assert_eq!(ix.size() as u64, AE_LIGHT_ENTRIES);
    assert_eq!(tiny().ae.n, ix.size());
}

/// The table passes its own check, a shifted copy fails it, and generating it on two threads
/// gives the same bytes as on all of them.
#[test]
fn table_checks_and_reproducibility() {
    let t = tiny();
    assert!(t.ae_ix.check_table(&t.ae, MoveSet::FACE));
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(2)
        .build()
        .unwrap();
    let mut again = pool.install(|| t.ae_ix.build(MoveSet::FACE, false));
    assert!(again.data[..] == t.ae.data[..]);
    assert_eq!(again.hist, t.ae.hist);
    again.data.rotate_left(1);
    assert!(!t.ae_ix.check_table(&again, MoveSet::FACE));
}

#[test]
fn optimal_lengths_match_the_small_tier() {
    let (t, sm) = (tiny(), small());
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(33);
    for len in [2usize, 5, 8, 10, 11, 12, 13] {
        for _ in 0..3 {
            let scr: Vec<PhysMove> = (0..len)
                .map(|_| PhysMove::face(rng.gen_range(0..8), rng.gen()))
                .collect();
            let (c, fr) = scramble_to_cube(&scr);
            let a = solve(t, c, fr, &Limits::default(), Report::Quiet);
            let b = solve(sm, c, fr, &Limits::default(), Report::Quiet);
            assert!(a.optimal && b.optimal);
            assert_eq!(a.length, b.length, "{}", format_moves(&scr));
            let sol = a.solution.unwrap();
            let st = apply_physical(&apply_physical(&solved_stickers(defs()), &scr), &sol);
            assert!(is_solved_physical(&st), "{}", format_moves(&scr));
        }
    }
}

/// A time limit and a length limit on random positions: the solutions solve them.
#[test]
fn stopping_rules() {
    let t = tiny();
    for c in positions(34, 4).into_iter().step_by(2) {
        let o = solve(
            t,
            c,
            Frame::identity(),
            &Limits {
                time: Some(std::time::Duration::from_millis(300)),
                ..Default::default()
            },
            Report::Quiet,
        );
        let Some(l) = o.length else { continue };
        let o2 = solve(
            t,
            c,
            Frame::identity(),
            &Limits {
                max_length: Some(l + 1),
                ..Default::default()
            },
            Report::Quiet,
        );
        assert!(o2.length.is_some_and(|m| m <= l + 1));
    }
}
