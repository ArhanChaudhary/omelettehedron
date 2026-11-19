//! The disk search of the phase-2 space (`omelettehedron::chen::p2bfs`) against independent
//! searches: the in-memory slot table on small quotient instances (several parts, stopped and
//! resumed), and a plain search without symmetry (positions per depth) on the full instance.

use omelettehedron::chen::p2bfs::engine::{
    canonical, for_slots, Bfs, Config, Layout, Manifest, Outcome,
};
use omelettehedron::chen::p2bfs::instance::{BView, CornerView, Instance, Kind};
use omelettehedron::chen::p2bfs::reference::{position_bfs, slot_bfs};
use std::path::{Path, PathBuf};

fn tmp_dir(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("p2bfs-test-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    d
}

fn config(dir: &Path, ram: usize) -> Config {
    Config {
        dir: dir.to_path_buf(),
        ram,
        max_depth: None,
        max_parts: None,
        max_seconds: None,
        threads: 4,
        keep: 0,
        stores: vec![],
        readers: 2,
        verbose: false,
    }
}

/// Run a search to the end, stopping after every `parts_per_run` parts and resuming in a new
/// `Bfs` (as after an interruption), its layers in `stores` subdirectories (none: in the
/// search directory); return its manifest.
fn run_resuming(
    inst: &Instance,
    dir: &Path,
    ram: usize,
    parts_per_run: Option<usize>,
    stores: usize,
) -> Manifest {
    loop {
        let mut cfg = config(dir, ram);
        cfg.max_parts = parts_per_run;
        cfg.stores = (0..stores).map(|k| dir.join(format!("store{k}"))).collect();
        let mut b = Bfs::open(inst, cfg).unwrap();
        if b.run().unwrap() == Outcome::Complete {
            return b.manifest().clone();
        }
    }
}

/// Every slot of the instance once, the orbit counts equal one canonical slot per orbit, and
/// the slots per depth equal the in-memory table's, for a run in one part and for a run in
/// many small parts that stops after each part and resumes.
fn check_instance(kind: Kind, name: &str) {
    let inst = Instance::new(kind).unwrap();
    let reference = slot_bfs(&inst, false);
    let dir = tmp_dir(name);
    let m = run_resuming(&inst, &dir, 1 << 34, None, 0);
    let slots: Vec<u64> = m.levels.iter().map(|l| l.slots).collect();
    assert_eq!(slots, reference.hist, "{name}: one part");
    assert_eq!(slots.iter().sum::<u64>(), inst.n_slots());
    let positions: u128 = m.levels.iter().map(|l| u128::from(l.positions)).sum();
    assert_eq!(positions, inst.n_positions());
    // the last layer is kept: its canonical slots are its orbits
    let last = m.levels.len() - 1;
    let mut n = 0u64;
    for_slots(
        &Layout::of(&dir, &m),
        &inst,
        last,
        0..inst.n_class,
        |c, x, p| {
            n += u64::from(canonical(&inst, c, x, p));
        },
    )
    .unwrap();
    assert_eq!(n, m.levels[last].orbits, "{name}: orbits of the last layer");
    std::fs::remove_dir_all(&dir).unwrap();
    // small parts, a stop after each
    let dir = tmp_dir(&format!("{name}-parts"));
    let ram = inst.n_x * inst.words * 8 * 40;
    let m2 = run_resuming(&inst, &dir, ram, Some(1), 3);
    let slots2: Vec<u64> = m2.levels.iter().map(|l| l.slots).collect();
    assert_eq!(slots2, reference.hist, "{name}: small parts, resumed");
    let orbits: Vec<u64> = m.levels.iter().map(|l| l.orbits).collect();
    let orbits2: Vec<u64> = m2.levels.iter().map(|l| l.orbits).collect();
    assert_eq!(orbits, orbits2);
    assert!(m2.levels.iter().any(|l| l.depth > 1), "{name}");
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn corner_permutation_instance() {
    check_instance(
        Kind {
            corners: CornerView::Perm,
            b: BView::None,
            radius: 0,
            t1: false,
        },
        "perm",
    );
}

#[test]
fn unlabelled_b_instance() {
    check_instance(
        Kind {
            corners: CornerView::None,
            b: BView::Unlabelled,
            radius: 0,
            t1: false,
        },
        "unlabelled",
    );
}

/// The full instance's first depths: positions per depth as the plain search without symmetry
/// gives them, and Chen Shuang's twsearch counts (1, 16, 160, 1344, 10464, 78152).
#[test]
fn full_instance_first_depths() {
    let inst = Instance::new(Kind::FULL).unwrap();
    let dir = tmp_dir("full");
    let mut cfg = config(&dir, 1 << 33);
    cfg.max_depth = Some(5);
    let mut b = Bfs::open(&inst, cfg).unwrap();
    assert_eq!(b.run().unwrap(), Outcome::Stopped);
    let positions: Vec<u64> = b.manifest().levels.iter().map(|l| l.positions).collect();
    assert_eq!(positions, [1, 16, 160, 1344, 10464, 78152]);
    assert_eq!(position_bfs(5), positions);
    // orbits by one canonical slot each (layers 4 and 5 are kept)
    for d in [4, 5] {
        let mut n = 0u64;
        for_slots(
            &Layout::of(&dir, b.manifest()),
            &inst,
            d,
            0..inst.n_class,
            |c, x, p| {
                n += u64::from(canonical(&inst, c, x, p));
            },
        )
        .unwrap();
        assert_eq!(n, b.manifest().levels[d].orbits);
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

/// The radius-1 instance's first depths: positions per depth, split by the configuration's
/// relaxed phase-1 distance, as a plain search over cube states gives them and as an
/// independent search on the sticker model gave them (1, 16, 208, 2,112, 18,144; at distance
/// 1: 0, 0, 48, 768, 7,680).
#[test]
fn radius1_instance_first_depths() {
    use omelettehedron::chen::p2bfs::engine::positions_by_d1;
    use omelettehedron::chen::p2bfs::reference::cube_bfs;
    use omelettehedron::chen::radius1::radius1;
    let inst = Instance::new(Kind::RADIUS1).unwrap();
    assert_eq!(inst.n_config, 315_900);
    assert_eq!(
        inst.d1.iter().map(|&d| usize::from(d == 0)).sum::<usize>(),
        1611
    );
    let dir = tmp_dir("radius1");
    let mut cfg = config(&dir, 1 << 33);
    cfg.max_depth = Some(4);
    let mut b = Bfs::open(&inst, cfg).unwrap();
    assert_eq!(b.run().unwrap(), Outcome::Stopped);
    let split: Vec<Vec<u64>> = b
        .manifest()
        .levels
        .iter()
        .map(|l| positions_by_d1(&inst, l))
        .collect();
    assert_eq!(
        split,
        [[1, 0], [16, 0], [160, 48], [1344, 768], [10464, 7680]]
    );
    let r = radius1().unwrap();
    assert_eq!(cube_bfs(4, |c| r.id(c).map(|id| r.cfg_d1(id))), split);
    std::fs::remove_dir_all(&dir).unwrap();
}

/// Membership in a layer (`joint::Phase2Slots`): every position of layer 5 of the full
/// instance tests in it, with its conjugates by the 24 symmetries that keep the face classes;
/// every position of layer 4, and positions a few allowed moves from solved, test out; and
/// the conversions between the pattern view and the fixed frame invert each other.
#[test]
fn layer_membership() {
    use omelettehedron::chen::p2bfs::joint::Phase2Slots;
    use omelettehedron::chen::p2bfs::solve::{PState, Views};
    use omelettehedron::chen::tables::P2;
    use omelettehedron::cube::{Cube, N_MOVES};
    use rand::{Rng, SeedableRng};
    let inst = Instance::new(Kind::FULL).unwrap();
    let dir = tmp_dir("member");
    let mut cfg = config(&dir, 1 << 33);
    cfg.max_depth = Some(5);
    let mut b = Bfs::open(&inst, cfg).unwrap();
    assert_eq!(b.run().unwrap(), Outcome::Stopped);
    let l5 = Phase2Slots::load(&dir, 5).unwrap();
    let l4 = Phase2Slots::load(&dir, 4).unwrap();
    assert_eq!(l5.reps.len() as u64, b.manifest().levels[5].orbits);
    let s48 = omelettehedron::optimal::sym::sym48();
    for (k, &slot) in l5.reps.iter().enumerate() {
        let c = l5.cube(slot).unwrap();
        assert_eq!(l5.contains(&c), Some(true));
        assert_eq!(l4.contains(&c), Some(false));
        if k % 50 == 0 {
            for s in 0..24 {
                let (cc, _) = s48.conj_cube(&c, s).unwrap();
                assert_eq!(l5.contains(&cc), Some(true));
            }
        }
    }
    for &slot in &l4.reps {
        let c = l4.cube(slot).unwrap();
        assert_eq!(l4.contains(&c), Some(true));
        assert_eq!(l5.contains(&c), Some(false));
    }
    // positions at most 4 allowed moves from solved are not at distance 5
    let s = omelettehedron::chen::space::space();
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(25);
    let views = Views::new();
    for _ in 0..2000 {
        let mut c = Cube::solved();
        for _ in 0..rng.gen_range(0..5) {
            let m = rng.gen_range(0..N_MOVES);
            if s.id(&c)
                .is_some_and(|id| s.next[id as usize][m] != u32::MAX)
            {
                c = c.apply(m);
            }
        }
        assert_eq!(l5.contains(&c), Some(false));
        // the pattern view and back
        let p2 = P2::of(&c).unwrap();
        let st: PState = views.of_p2(&p2);
        let back = views.to_p2(&st);
        assert_eq!(
            (back.cfg, back.c, &back.xb[..12]),
            (p2.cfg, p2.c, &p2.xb[..12])
        );
        // a pattern-view path as fixed-frame moves
        let mut path = vec![];
        let mut cur = st;
        for _ in 0..6 {
            let m = rng.gen_range(0..N_MOVES);
            if let Some(n) = views.apply(&cur, m) {
                path.push(m);
                cur = n;
            }
        }
        let fixed = views.to_fixed(&st, &path).unwrap();
        let mut q = p2;
        for &m in &fixed {
            q = q.apply(s, m).unwrap();
        }
        let want = views.to_p2(&cur);
        assert_eq!(
            (q.cfg, q.c, &q.xb[..12]),
            (want.cfg, want.c, &want.xb[..12])
        );
    }
    std::fs::remove_dir_all(&dir).unwrap();
}
