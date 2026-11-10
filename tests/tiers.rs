//! Cross-check of the symmetry-reduced table against the relabel-only table.
//! Needs both tables on disk: `fto gen --tier medium && fto gen --tier large`, then
//! `cargo test --release -- --ignored`.

use omelettehedron::coord::*;
use omelettehedron::cube::*;
use omelettehedron::tables::*;
use rand::SeedableRng;

#[test]
#[ignore = "needs the medium and huge tables on disk"]
fn medium_and_large_tables_agree() {
    let dir = table_dir(None);
    let med = Tables::load(&dir, AeSym::new(true), MoveSet::FACE, false).unwrap();
    let large = Tables::load(&dir, AeRelabel::new(), MoveSet::FACE, false).unwrap();
    let mt = move_tables();
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(21);
    for _ in 0..3000 {
        let c = Cube::random(&mut rng);
        let x = Coords::of(&c);
        assert_eq!(med.ae_exact(&x), large.ae_exact(&x));
        for m in 0..N_MOVES {
            let y = mt.apply(&x, m);
            assert_eq!(
                med.ae.get(med.ae_ix.index(&y)),
                large.ae.get(large.ae_ix.index(&y))
            );
        }
    }
}
