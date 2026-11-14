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
    let inst = Instance::new(kind);
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
        },
        "unlabelled",
    );
}

/// The full instance's first depths: positions per depth as the plain search without symmetry
/// gives them, and Chen Shuang's twsearch counts (1, 16, 160, 1344, 10464, 78152).
#[test]
fn full_instance_first_depths() {
    let inst = Instance::new(Kind::FULL);
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
