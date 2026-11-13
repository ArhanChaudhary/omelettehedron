//! Table design study: the distance to solved of the whole edge x A-centre projection (without
//! the quotient by `C`, an 83.6 GB table), on the positions that `tabsample` samples, without
//! building it: IDA* in the projection's own graph, bounded by the quotient table (a
//! projection of it). One byte per side per position (capped at `cap`), in the order of
//! `tabsample`'s records.
//!
//! usage: cargo run --release --example fullea -- ea.prun N seed cap out.bin

use omelettehedron::cube::{allowed_after, Cube, N_MOVES};
use omelettehedron::optimal::coords::{opt_moves, Pos};
use omelettehedron::optimal::ea::EaIndex;
use omelettehedron::prune::{decode, Table2};
use omelettehedron::sym::syms;
use rand::SeedableRng;
use rayon::prelude::*;
use std::io::Write;

struct Ctx {
    ix: EaIndex,
    tab: Table2,
    targets: Vec<usize>,
    /// the solved state's 12 rotations (A-colouring, edges)
    solved: Vec<([u8; 16], [u8; 16])>,
}

impl Ctx {
    /// Exact value of the quotient table at `(xa, ep)` (walked down to solved).
    fn exact(&self, xa: &[u8; 16], ep: &[u8; 16]) -> u8 {
        let mt = opt_moves();
        let (mut a, mut e) = (*xa, *ep);
        let mut d = 0u8;
        loop {
            let i = self.ix.index_raw(&a, &e);
            if self.targets.contains(&i) {
                return d;
            }
            let want = (self.tab.get(i) + 2) % 3;
            let Some((a2, e2)) = (0..N_MOVES)
                .map(|m| mt.apply_ea(&a, &e, m))
                .find(|(a2, e2)| self.tab.get(self.ix.index_raw(a2, e2)) == want)
            else {
                return u8::MAX;
            };
            (a, e) = (a2, e2);
            d += 1;
        }
    }

    fn is_solved(&self, xa: &[u8; 16], ep: &[u8; 16]) -> bool {
        self.solved.iter().any(|(a, e)| a == xa && e == ep)
    }

    /// Depth-first search below `bound`: whether a solved rotation is within `bound - g`.
    fn dfs(&self, xa: &[u8; 16], ep: &[u8; 16], h: u8, g: u8, bound: u8, last: usize) -> bool {
        if h == 0 && self.is_solved(xa, ep) {
            return true;
        }
        if g >= bound {
            return false;
        }
        let mt = opt_moves();
        for m in 0..N_MOVES {
            if !allowed_after(false, last, m) {
                continue;
            }
            let (a2, e2) = mt.apply_ea(xa, ep, m);
            let h2 = decode(h, self.tab.get(self.ix.index_raw(&a2, &e2)));
            if g + 1 + h2 <= bound && self.dfs(&a2, &e2, h2, g + 1, bound, m) {
                return true;
            }
        }
        false
    }

    /// The whole projection's distance of `(xa, ep)`, at most `cap`.
    fn full(&self, xa: &[u8; 16], ep: &[u8; 16], cap: u8) -> u8 {
        let h = self.exact(xa, ep);
        (h..cap)
            .find(|&b| self.dfs(xa, ep, h, 0, b, usize::MAX))
            .unwrap_or(cap)
    }
}

fn main() -> std::io::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let path = args.get(1).cloned().unwrap_or_default();
    let n: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(1000);
    let seed: u64 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(7);
    let cap: u8 = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(15);
    let out = args.get(5).cloned().unwrap_or_else(|| "fullea.bin".into());
    let ix = EaIndex::new()?;
    let tab = Table2::load(std::path::Path::new(&path), ix.size())?;
    let targets = ix.targets();
    let s = syms();
    let sp = Pos::solved();
    let mut sl = [0u8; 12];
    sl.copy_from_slice(&sp.xa[..12]);
    let solved = (0..s.mats.len())
        .filter(|&t| omelettehedron::geometry::det(&s.mats[t]) > 0)
        .map(|t| {
            let la = s.apply_a(0, t, &sl);
            let mut a = [0u8; 16];
            a[..12].copy_from_slice(&la);
            let mut e = [0u8; 16];
            for p in 0..12 {
                e[usize::from(s.e_pos[t][p])] = sp.ep[p];
            }
            (a, e)
        })
        .collect();
    let cx = Ctx {
        ix,
        tab,
        targets,
        solved,
    };
    let w = omelettehedron::optimal::sym::sym48();
    let sw = w.swap;
    const CHUNK: usize = 4096;
    let t0 = std::time::Instant::now();
    let recs: Vec<[u8; 4]> = (0..n.div_ceil(CHUNK))
        .into_par_iter()
        .flat_map_iter(|ch| {
            let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(
                seed.wrapping_mul(0x9e37_79b9_7f4a_7c15) ^ ch as u64,
            );
            let m = CHUNK.min(n - ch * CHUNK);
            let mut v = Vec::with_capacity(m);
            for _ in 0..m {
                let c = Cube::random(&mut rng);
                let p = Pos::of(&c);
                let mut r = [0u8; 4];
                r[0] = cx.full(&p.xa, &p.ep, cap);
                r[2] = cx.exact(&p.xa, &p.ep);
                // the SWAP_AB conjugate (B-centres as A-centres)
                let mut xa = [0u8; 16];
                let mut ep = [0u8; 16];
                for q in 0..12 {
                    xa[usize::from(w.cpos[sw][1][q])] = w.ccol[sw][1][usize::from(p.xb[q])];
                    ep[usize::from(w.epos[sw][q])] = w.epos[sw][usize::from(p.ep[q])];
                }
                r[1] = cx.full(&xa, &ep, cap);
                r[3] = cx.exact(&xa, &ep);
                v.push(r);
            }
            v
        })
        .collect();
    eprintln!("{} records in {:.1}s", recs.len(), t0.elapsed().as_secs_f64());
    let mut f = std::io::BufWriter::new(std::fs::File::create(&out)?);
    for r in &recs {
        f.write_all(r)?;
    }
    Ok(())
}
