//! Table design study: exact values of every optimal-solver lookup on sampled positions.
//!
//! Writes one record per position: the large tier's ten lookups (A-centre/edge on both sides,
//! the eight corner x B-centre x edge-triple lookups), then corners x edges and the two corner x
//! centre lookups of the medium tier (13 bytes). The positions are uniformly random
//! (`Cube::random`), or random-move scrambles of `len` moves with `len > 0`.
//!
//! With a fifth argument (an edge x A-centre table, `optimal::ea`), two more bytes: its lookup
//! on the position and on its `SWAP_AB` conjugate (15 bytes).
//!
//! usage: FTO_TABLES=... cargo run --release --example tabsample -- N seed len out.bin [ea.prun]

use omelettehedron::cube::{Cube, MoveSet};
use omelettehedron::optimal::coords::Pos;
use omelettehedron::optimal::ea::EaIndex;
use omelettehedron::prune::Table2;
use omelettehedron::optimal::tables::{exact, Bounds, OptTables, OptTier};
use rand::SeedableRng;
use rayon::prelude::*;
use std::io::Write;

fn main() -> std::io::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let n: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(100_000);
    let seed: u64 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(1);
    let len: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(0);
    let out = args.get(4).cloned().unwrap_or_else(|| "samples.bin".into());
    let dir = omelettehedron::tables::table_dir(None);
    let lt = OptTables::load(&dir, OptTier::Large, MoveSet::FACE, true)?;
    let mt = OptTables::load(&dir, OptTier::Medium, MoveSet::FACE, true)?;
    let (Some(lb), Some(mb)) = (lt.large(), mt.medium()) else {
        return Err(std::io::Error::other("tiers"));
    };
    let ea = match args.get(5) {
        Some(path) => {
            let ix = EaIndex::new()?;
            let tab = Table2::load(std::path::Path::new(path), ix.size())?;
            let targets = ix.targets();
            Some(EaBounds { ix, tab, targets })
        }
        None => None,
    };
    let width = if ea.is_some() { 15 } else { 13 };
    const CHUNK: usize = 4096;
    let t0 = std::time::Instant::now();
    let recs: Vec<[u8; 15]> = (0..n.div_ceil(CHUNK))
        .into_par_iter()
        .flat_map_iter(|ch| {
            let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(
                seed.wrapping_mul(0x9e37_79b9_7f4a_7c15) ^ ch as u64,
            );
            let m = CHUNK.min(n - ch * CHUNK);
            let mut v = Vec::with_capacity(m);
            for _ in 0..m {
                let c = if len == 0 {
                    Cube::random(&mut rng)
                } else {
                    use rand::Rng;
                    let mut moves = vec![];
                    let mut last = 8;
                    while moves.len() < len {
                        let face = rng.gen_range(0..8);
                        if face != last {
                            moves.push(omelettehedron::notation::PhysMove::face(face, rng.gen()));
                            last = face;
                        }
                    }
                    omelettehedron::notation::scramble_to_cube(&moves).0
                };
                let p = Pos::of(&c);
                let (Some(hl), Some(hm)) = (exact(&lb, &p), exact(&mb, &p)) else {
                    continue;
                };
                let mut r = [0u8; 15];
                r[..10].copy_from_slice(&hl[..10]);
                r[10] = hm[0];
                r[11] = hm[1];
                r[12] = hm[2];
                if let Some(e) = &ea {
                    let Some(he) = exact(e, &p) else {
                        continue;
                    };
                    r[13] = he[0];
                    r[14] = he[1];
                }
                v.push(r);
            }
            v
        })
        .collect();
    eprintln!(
        "{} records in {:.1}s ({} lookups per record: large {}, medium {})",
        recs.len(),
        t0.elapsed().as_secs_f64(),
        width,
        LARGE_K,
        3
    );
    let mut f = std::io::BufWriter::new(std::fs::File::create(&out)?);
    for r in &recs {
        f.write_all(&r[..width])?;
    }
    Ok(())
}

const LARGE_K: usize = <omelettehedron::optimal::tables::LargeBounds<'static> as Bounds>::K;

/// The edge x A-centre table's two lookups, for [`exact`].
struct EaBounds {
    ix: EaIndex,
    tab: Table2,
    targets: Vec<usize>,
}

impl Bounds for EaBounds {
    const K: usize = 2;
    fn index(&self, k: usize, p: &Pos) -> usize {
        self.ix.index(k, p)
    }
    fn table(&self, _k: usize) -> &Table2 {
        &self.tab
    }
    fn is_target(&self, _k: usize, i: usize) -> bool {
        self.targets.contains(&i)
    }
    fn name(&self, k: usize) -> &'static str {
        ["EA", "EA of SWAP_AB"][k & 1]
    }
    fn moves(&self) -> MoveSet {
        MoveSet::FACE
    }
}
