//! A table of another layout (here: the right size, values moved) must be rejected when loaded,
//! not walked forever or used.

use omelettehedron::cube::*;
use omelettehedron::fast::fast;
use omelettehedron::tables::*;
use rand::SeedableRng;
use std::path::Path;

#[test]
fn walks_reject_a_table_of_another_layout() {
    let mut t = Tables::generate(AeSym::new(false), MoveSet::FACE, false);
    let (ae, cb) = (Path::new("ae.prun"), Path::new("cb.prun"));
    t.check_walks(ae, cb).expect("the generated tables pass");

    // every exact distance stays within the table's largest distance
    let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(5);
    for _ in 0..50 {
        let (ha, hc) = t.exact_fast(&fast().of(&Cube::random(&mut rng)));
        assert!(usize::from(ha) < t.ae.hist.len() && usize::from(hc) < t.cb.tab.hist.len());
    }

    // the AE values moved by one word (the same size and histogram, another "layout")
    t.ae.data.rotate_left(1);
    let e = t
        .check_walks(ae, cb)
        .expect_err("a moved AE table must fail");
    assert!(e.to_string().contains("ae.prun"), "{e}");
    t.ae.data.rotate_right(1);
    t.check_walks(ae, cb).expect("restored");

    // the same for the CB table
    t.cb.tab.data.rotate_left(3);
    let e = t
        .check_walks(ae, cb)
        .expect_err("a moved CB table must fail");
    assert!(e.to_string().contains("cb.prun"), "{e}");
    // and the solver's own walk ends instead of running forever
    let _ = t.exact_fast(&fast().of(&Cube::random(&mut rng)));
}
