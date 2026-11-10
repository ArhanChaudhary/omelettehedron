//! Experiment (FTO_CBA=1): joint corner x B-centre x one-A-colour table ("CBA"), the analogue of
//! CBE with the positions of the three centres of A-colour k instead of edge triple k.
//! Built in memory at load time; saturating codes as in CBE (<=6, 7, 8, >=9).

use crate::cbe::{conj_renorm_corner_with, CODE_LB};
use crate::coord::*;
use crate::cube::*;
use crate::fast::*;
use crate::mem::HugeBuf;
use crate::prune::*;
use crate::sym::*;
use crate::tables::{prefetch, CbTable};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::sync::atomic::{AtomicU64, Ordering};

pub const N_A1: usize = 220;
pub const N_ROW: usize = N_MSET * N_A1;

pub struct CbaIndex {
    kt: Box<[[(u16, u8, u8); 4]; 1024]>,
    pub rep: Vec<u16>,
    stab: Vec<Vec<(u8, u8)>>,
    bmap: [[u8; 16]; 32],
    bgather: [[u8; 16]; 32],
    agather: [[u8; 16]; 32],
    c12_unrank: Vec<u16>,
}

impl CbaIndex {
    pub fn new() -> CbaIndex {
        let s = syms();
        let to0 = |k: usize| -> Vec<usize> { (0..N_GEO).filter(|&m| s.a_col[m][k] == 0).collect() };
        let s0 = to0(0);
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
            .map(|c| s0.iter().map(|&m| conj[m][c].0).min().unwrap())
            .collect();
        rep.sort();
        rep.dedup();
        let cls_id: HashMap<u16, u16> = rep
            .iter()
            .enumerate()
            .map(|(i, &r)| (r, i as u16))
            .collect();
        let mut kt = Box::new([[(0u16, 0u8, 0u8); 4]; 1024]);
        for k in 0..4 {
            for c in 0..N_CORNER {
                let (m, (c2, r)) = to0(k)
                    .into_iter()
                    .map(|m| (m, conj[m][c]))
                    .min_by_key(|&(_, (c2, _))| c2)
                    .unwrap();
                kt[c][k] = (cls_id[&c2], m as u8, s.mul[r][m]);
            }
        }
        let stab = rep
            .iter()
            .map(|&r| {
                s0.iter()
                    .filter(|&&m| m != 0 && conj[m][r as usize].0 == r)
                    .map(|&m| (m as u8, s.mul[conj[m][r as usize].1][m]))
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
        let gather = |pos: &[[u8; 12]; N_GEO], t: usize| {
            let mut g = [0x80u8; 16];
            if t < N_GEO {
                for p in 0..12 {
                    g[pos[t][p] as usize] = p as u8;
                }
            }
            g
        };
        let mut c12_unrank = vec![0u16; 220];
        let f = fast();
        for m in 0..4096usize {
            if m.count_ones() == 3 {
                c12_unrank[f.c12[m] as usize] = m as u16;
            }
        }
        CbaIndex {
            kt,
            rep,
            stab,
            bmap: std::array::from_fn(pad_map),
            bgather: std::array::from_fn(|t| gather(&s.b_pos, t)),
            agather: std::array::from_fn(|t| gather(&s.a_pos, t)),
            c12_unrank,
        }
    }

    pub fn size(&self) -> usize {
        self.rep.len() * N_ROW
    }

    #[inline(always)]
    fn row(&self, f: &FastTables, x: &Fc, k: usize, m: usize, t: usize) -> usize {
        let (m, t) = (m & 31, t & 31);
        let b = shuffle(&shuffle(&self.bmap[m], &x.b), &self.bgather[t]);
        let mask = label_mask(&shuffle(&x.a, &self.agather[t]), k as u8) as usize;
        f.rank(&masks(&b)) as usize * N_A1 + f.c12[mask & 0xfff] as usize
    }

    #[inline(always)]
    pub fn index(&self, f: &FastTables, x: &Fc, k: usize) -> usize {
        let (cls, m, t) = self.kt[x.c as usize & 1023][k];
        cls as usize * N_ROW + self.row(f, x, k, m as usize, t as usize)
    }

    fn decode(&self, i: usize) -> Fc {
        let cls = i / N_ROW;
        let r = i % N_ROW;
        let (b, a1) = (r / N_A1, r % N_A1);
        let mask = self.c12_unrank[a1];
        let mut x = Fc {
            a: [0; 16],
            b: [0; 16],
            e: [0; 16],
            c: self.rep[cls],
            ebits: 0,
        };
        x.b[..12].copy_from_slice(&mset().unrank(b as u32));
        for p in 0..12 {
            x.a[p] = if mask >> p & 1 == 1 { 0 } else { 1 };
        }
        x
    }

    fn equivalents(&self, f: &FastTables, i: usize, out: &mut Vec<usize>) {
        let cls = i / N_ROW;
        if self.stab[cls].is_empty() {
            return;
        }
        let x = self.decode(i);
        for &(m, t) in &self.stab[cls] {
            out.push(cls * N_ROW + self.row(f, &x, 0, m as usize, t as usize));
        }
    }
}

pub struct CbaTable {
    pub ix: CbaIndex,
    pub codes: HugeBuf<u64>,
    pub hist: Vec<u64>,
}

impl CbaTable {
    #[inline(always)]
    pub fn code(&self, i: usize) -> u8 {
        ((self.codes[i >> 5] >> ((i & 31) << 1)) & 3) as u8
    }

    pub fn build(verbose: bool) -> CbaTable {
        let ix = CbaIndex::new();
        let f = fast();
        let n = ix.size();
        eprintln!(
            "cba: {} corner classes, {} entries ({:.1} GB per 2-bit array)",
            ix.rep.len(),
            n,
            n as f64 / 4e9
        );
        // targets: every position of H_A (BFS with the class-A moves)
        let amoves: Vec<usize> = (0..N_MOVES).filter(|&m| is_a_move(m)).collect();
        let mut seen: HashSet<Cube> = HashSet::new();
        seen.insert(Cube::solved());
        let mut fr = vec![Cube::solved()];
        while !fr.is_empty() {
            let mut nf = vec![];
            for c in &fr {
                for &m in &amoves {
                    let y = c.apply(m);
                    if seen.insert(y) {
                        nf.push(y);
                    }
                }
            }
            fr = nf;
        }
        let _ = CbTable::attach_targets;
        let mut targets: Vec<usize> = vec![];
        for c in &seen {
            let x = f.of(c);
            let i = ix.index(f, &x, 0);
            targets.push(i);
            ix.equivalents(f, i, &mut targets);
        }
        targets.sort_unstable();
        targets.dedup();
        eprintln!("cba: {} H_A positions, {} target indices", seen.len(), targets.len());
        drop(seen);
        let codes_raw: Vec<AtomicU64> = (0..n.div_ceil(32))
            .into_par_iter()
            .map(|_| AtomicU64::new(0))
            .collect();
        let tab = bfs_with(
            n,
            &targets,
            |i, out| {
                let x = ix.decode(i);
                let em = masks(&x.e);
                for m in 0..N_MOVES {
                    let y = f.apply(&x, &em, m);
                    let j = ix.index(f, &y, 0);
                    out.push(j);
                    ix.equivalents(f, j, out);
                }
            },
            |i, d| {
                if d >= 7 {
                    let code = (d - 6).min(3) as u64;
                    codes_raw[i >> 5].fetch_or(code << ((i & 31) << 1), Ordering::Relaxed);
                }
            },
            "cba",
            verbose,
        );
        let hist = tab.hist.clone();
        drop(tab);
        let mut codes = HugeBuf::<u64>::zeroed(n.div_ceil(32));
        codes
            .par_iter_mut()
            .zip(codes_raw.par_iter())
            .for_each(|(o, a)| *o = a.load(Ordering::Relaxed));
        let tot: u64 = hist.iter().sum();
        let mean = hist.iter().enumerate().map(|(d, &h)| d as f64 * h as f64).sum::<f64>() / tot as f64;
        eprintln!("cba: hist {:?} mean {:.3}", hist, mean);
        CbaTable { ix, codes, hist }
    }

    #[inline(always)]
    pub fn prefetch_all(&self, f: &FastTables, x: &Fc, out: &mut [usize; 4]) {
        for k in 0..4 {
            out[k] = self.ix.index(f, x, k);
            prefetch(&self.codes[out[k] >> 5]);
        }
    }

    #[inline(always)]
    pub fn bound(&self, idx: &[usize; 4]) -> u8 {
        idx.iter()
            .map(|&i| CODE_LB[self.code(i) as usize])
            .max()
            .unwrap()
    }
}
