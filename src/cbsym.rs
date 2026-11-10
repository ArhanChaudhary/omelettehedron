//! Experiment: the CB table reduced by all 24 class-preserving symmetries (conjugation and
//! re-fixing the front corner), ~3.7 MB instead of 89 MB, so that it stays cache-resident.
//! Enabled with FTO_CBSYM=1 (built in memory at load time).

use crate::cbe::conj_renorm_corner_with;
use crate::coord::*;
use crate::cube::*;
use crate::fast::*;
use crate::prune::*;
use crate::sym::*;
use crate::tables::CbTable;
use rayon::prelude::*;
use std::collections::HashMap;

pub struct CbSymIndex {
    /// per corner id (padded to 1024): (class id, symmetry M, reposition R*M)
    kt: Box<[(u16, u8, u8); 1024]>,
    pub rep: Vec<u16>,
    stab: Vec<Vec<(u8, u8)>>,
    bmap: [[u8; 16]; 32],
    bgather: [[u8; 16]; 32],
}

impl CbSymIndex {
    pub fn new() -> CbSymIndex {
        let s = syms();
        let wps: Vec<Vec<usize>> = s.mats.iter().map(|r| defs().sm.whole(r)).collect();
        let conj: Vec<Vec<(u16, usize)>> = (0..N_GEO)
            .into_par_iter()
            .map(|m| {
                (0..N_CORNER as u16)
                    .map(|c| conj_renorm_corner_with(c, m, &wps).expect("re-fix"))
                    .collect()
            })
            .collect();
        let mut rep: Vec<u16> = (0..N_CORNER)
            .map(|c| (0..N_GEO).map(|m| conj[m][c].0).min().unwrap())
            .collect();
        rep.sort();
        rep.dedup();
        let cls_id: HashMap<u16, u16> = rep
            .iter()
            .enumerate()
            .map(|(i, &r)| (r, i as u16))
            .collect();
        let mut kt = Box::new([(0u16, 0u8, 0u8); 1024]);
        for c in 0..N_CORNER {
            let (m, (c2, r)) = (0..N_GEO)
                .map(|m| (m, conj[m][c]))
                .min_by_key(|&(_, (c2, _))| c2)
                .unwrap();
            kt[c] = (cls_id[&c2], m as u8, s.mul[r][m]);
        }
        let stab = rep
            .iter()
            .map(|&r| {
                (1..N_GEO)
                    .filter(|&m| conj[m][r as usize].0 == r)
                    .map(|m| (m as u8, s.mul[conj[m][r as usize].1][m]))
                    .collect()
            })
            .collect();
        let pad_map = |m: usize| {
            let mut v = [0u8; 16];
            if m < N_GEO {
                v[..4].copy_from_slice(&s.b_col[m]);
            }
            v
        };
        let gather = |t: usize| {
            let mut g = [0x80u8; 16];
            if t < N_GEO {
                for p in 0..12 {
                    g[s.b_pos[t][p] as usize] = p as u8;
                }
            }
            g
        };
        CbSymIndex {
            kt,
            rep,
            stab,
            bmap: std::array::from_fn(pad_map),
            bgather: std::array::from_fn(gather),
        }
    }

    pub fn size(&self) -> usize {
        self.rep.len() * N_MSET
    }

    #[inline(always)]
    fn row(&self, f: &FastTables, b: &[u8; 16], m: usize, t: usize) -> usize {
        let b = shuffle(&shuffle(&self.bmap[m & 31], b), &self.bgather[t & 31]);
        f.rank(&masks(&b)) as usize
    }

    #[inline(always)]
    pub fn index(&self, f: &FastTables, x: &Fc) -> usize {
        let (cls, m, t) = self.kt[x.c as usize & 1023];
        cls as usize * N_MSET + self.row(f, &x.b, m as usize, t as usize)
    }

    fn decode(&self, i: usize) -> Fc {
        let mut x = Fc {
            a: [0; 16],
            b: [0; 16],
            e: [0; 16],
            c: self.rep[i / N_MSET],
            ebits: 0,
        };
        x.b[..12].copy_from_slice(&mset().unrank((i % N_MSET) as u32));
        x
    }

    fn equivalents(&self, f: &FastTables, i: usize, out: &mut Vec<usize>) {
        let cls = i / N_MSET;
        if self.stab[cls].is_empty() {
            return;
        }
        let x = self.decode(i);
        for &(m, t) in &self.stab[cls] {
            out.push(cls * N_MSET + self.row(f, &x.b, m as usize, t as usize));
        }
    }

    /// Exact distances (u8) to the CB projection of a target set given by corner ids
    /// (each with its attached B colouring).
    pub fn exact_to(&self, corners_in: &[usize], attach: &[u32]) -> Vec<u8> {
        use std::sync::atomic::{AtomicU8, Ordering};
        let f = fast();
        let mut targets: Vec<usize> = vec![];
        for &c in corners_in {
            let mut x = self.decode(0);
            x.c = c as u16;
            x.b[..12].copy_from_slice(&mset().unrank(attach[c]));
            let i = self.index(f, &x);
            targets.push(i);
            self.equivalents(f, i, &mut targets);
        }
        targets.sort_unstable();
        targets.dedup();
        let ex: Vec<AtomicU8> = (0..self.size()).map(|_| AtomicU8::new(0)).collect();
        bfs_with(
            self.size(),
            &targets,
            |i, out| {
                let x = self.decode(i);
                let em = masks(&x.e);
                for m in 0..N_MOVES {
                    let y = f.apply(&x, &em, m);
                    let j = self.index(f, &y);
                    out.push(j);
                    self.equivalents(f, j, out);
                }
            },
            |i, d| ex[i].store(d as u8, Ordering::Relaxed),
            "cb-route",
            false,
        );
        ex.into_iter().map(AtomicU8::into_inner).collect()
    }

    pub fn build(&self, verbose: bool) -> Table2 {
        let f = fast();
        let attach = CbTable::attach_targets();
        let mut targets: Vec<usize> = vec![];
        for c in 0..N_CORNER {
            let mut x = self.decode(0);
            x.c = c as u16;
            x.b[..12].copy_from_slice(&mset().unrank(attach[c]));
            let i = self.index(f, &x);
            targets.push(i);
            self.equivalents(f, i, &mut targets);
        }
        targets.sort_unstable();
        targets.dedup();
        bfs(
            self.size(),
            &targets,
            |i, out| {
                let x = self.decode(i);
                let em = masks(&x.e);
                for m in 0..N_MOVES {
                    let y = f.apply(&x, &em, m);
                    let j = self.index(f, &y);
                    out.push(j);
                    self.equivalents(f, j, out);
                }
            },
            "cb-sym24",
            verbose,
        )
    }
}
