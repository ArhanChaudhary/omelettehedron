//! Two-bit pruning tables storing distance mod 3 (3 = not yet visited), built by parallel BFS.
//!
//! A search tracks the exact distance of the parent; since one move changes the distance by at
//! most one, the child's exact distance follows from its value mod 3.

use crate::mem::HugeBuf;
use rayon::prelude::*;
use std::io::{Read, Write};
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};

pub struct Table2 {
    pub data: HugeBuf<u64>,
    pub n: usize,
    /// number of entries at each distance
    pub hist: Vec<u64>,
}

/// Decode the exact child distance from the parent distance and the child's mod-3 value.
#[inline(always)]
pub fn decode(parent: u8, v: u8) -> u8 {
    // v == parent-1 (mod 3) -> parent-1 ; v == parent -> parent ; else parent+1
    const T: [[u8; 4]; 3] = [[1, 2, 0, 1], [0, 1, 2, 1], [2, 0, 1, 1]];
    // T[parent % 3][v] = delta+1
    (parent + T[(parent % 3) as usize][v as usize]).wrapping_sub(1)
}

impl Table2 {
    #[inline(always)]
    pub fn get(&self, i: usize) -> u8 {
        // SAFETY-free bounds-checked access; the optimiser removes most checks
        ((self.data[i >> 5] >> ((i & 31) << 1)) & 3) as u8
    }

    pub fn bytes(&self) -> usize {
        self.data.len() * 8
    }

    pub fn save(&self, path: &Path) -> std::io::Result<()> {
        let tmp = path.with_extension("tmp");
        let mut f = std::io::BufWriter::with_capacity(1 << 24, std::fs::File::create(&tmp)?);
        f.write_all(b"FTOPRUN2")?;
        f.write_all(&(self.n as u64).to_le_bytes())?;
        f.write_all(&(self.hist.len() as u64).to_le_bytes())?;
        for h in &self.hist {
            f.write_all(&h.to_le_bytes())?;
        }
        f.write_all(&checksum(&self.data).to_le_bytes())?;
        // SAFETY: u64 slice reinterpreted as bytes for writing
        let bytes = unsafe {
            std::slice::from_raw_parts(self.data.as_ptr() as *const u8, self.data.len() * 8)
        };
        for chunk in bytes.chunks(1 << 26) {
            f.write_all(chunk)?;
        }
        f.flush()?;
        drop(f);
        std::fs::rename(tmp, path)
    }

    pub fn load(path: &Path, n: usize) -> std::io::Result<Table2> {
        let err = |m: &str| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("{}: {m}", path.display()),
            )
        };
        let mut f = std::io::BufReader::with_capacity(1 << 24, std::fs::File::open(path)?);
        let mut b8 = [0u8; 8];
        f.read_exact(&mut b8)?;
        if &b8 != b"FTOPRUN2" {
            return Err(err("bad magic"));
        }
        f.read_exact(&mut b8)?;
        if u64::from_le_bytes(b8) as usize != n {
            return Err(err("size mismatch (table from a different version?)"));
        }
        f.read_exact(&mut b8)?;
        let nh = u64::from_le_bytes(b8) as usize;
        let mut hist = vec![];
        for _ in 0..nh {
            f.read_exact(&mut b8)?;
            hist.push(u64::from_le_bytes(b8));
        }
        f.read_exact(&mut b8)?;
        let sum = u64::from_le_bytes(b8);
        let words = n.div_ceil(32);
        let mut data = alloc_words(words);
        // SAFETY: reading raw bytes into the u64 buffer
        let bytes =
            unsafe { std::slice::from_raw_parts_mut(data.as_mut_ptr() as *mut u8, words * 8) };
        f.read_exact(bytes)?;
        if checksum(&data) != sum {
            return Err(err("checksum mismatch"));
        }
        Ok(Table2 { data, n, hist })
    }
}

fn checksum(d: &[u64]) -> u64 {
    d.par_chunks(1 << 20)
        .enumerate()
        .map(|(i, c)| {
            let mut h = 0x9e3779b97f4a7c15u64 ^ i as u64;
            for &x in c {
                h = (h ^ x).wrapping_mul(0x100000001b3).rotate_left(23);
            }
            h
        })
        .reduce(|| 0, |a, b| a ^ b) // associative: rayon may combine in any order
}

/// Allocate a zeroed, huge-page backed word buffer.
pub fn alloc_words(words: usize) -> HugeBuf<u64> {
    HugeBuf::zeroed(words)
}

/// Breadth-first search over `n` indices. `nbrs(i, out)` appends the neighbours of `i` to
/// `out` (including every index equivalent to a neighbour, if indices are not unique).
pub fn bfs<F>(n: usize, targets: &[usize], nbrs: F, name: &str, verbose: bool) -> Table2
where
    F: Fn(usize, &mut Vec<usize>) + Sync,
{
    bfs_with(n, targets, nbrs, |_, _| {}, name, verbose)
}

/// [`bfs`] that also calls `on_set(index, depth)` for every entry when its depth is found
/// (depth >= 1; targets are depth 0).
pub fn bfs_with<F, S>(
    n: usize,
    targets: &[usize],
    nbrs: F,
    on_set: S,
    name: &str,
    verbose: bool,
) -> Table2
where
    F: Fn(usize, &mut Vec<usize>) + Sync,
    S: Fn(usize, usize) + Sync,
{
    let t0 = std::time::Instant::now();
    let words = n.div_ceil(32);
    let mut raw = alloc_words(words);
    raw.par_chunks_mut(1 << 16).for_each(|c| c.fill(u64::MAX));
    // SAFETY: AtomicU64 has the same size and alignment as u64; `raw` outlives `tab`
    let tab: &[AtomicU64] =
        unsafe { std::slice::from_raw_parts(raw.as_ptr() as *const AtomicU64, raw.len()) };
    let get = |i: usize| ((tab[i >> 5].load(Ordering::Relaxed) >> ((i & 31) << 1)) & 3) as u8;
    let set_if_empty = |i: usize, v: u8| -> bool {
        let w = &tab[i >> 5];
        let sh = (i & 31) << 1;
        let mut cur = w.load(Ordering::Relaxed);
        loop {
            if (cur >> sh) & 3 != 3 {
                return false;
            }
            let new = cur & !(3u64 << sh) | (v as u64) << sh;
            match w.compare_exchange_weak(cur, new, Ordering::Relaxed, Ordering::Relaxed) {
                Ok(_) => return true,
                Err(x) => cur = x,
            }
        }
    };
    let mut nt = 0u64;
    for &t in targets {
        if set_if_empty(t, 0) {
            nt += 1;
        }
    }
    let mut hist = vec![nt];
    let mut done = nt;
    let mut d = 0usize;
    const CH: usize = 1 << 14; // indices per work item (multiple of 32)
    while done < n as u64 {
        let cur = (d % 3) as u8;
        let next = ((d + 1) % 3) as u8;
        let frontier = hist[d];
        let forward = frontier < (n as u64 - done) / 3;
        let newc: u64 = (0..n.div_ceil(CH))
            .into_par_iter()
            .map(|c| {
                let mut nb: Vec<usize> = Vec::with_capacity(64);
                let mut cnt = 0u64;
                let lo = c * CH;
                let hi = ((c + 1) * CH).min(n);
                let mut i = lo;
                while i < hi {
                    let w = tab[i >> 5].load(Ordering::Relaxed);
                    let want = if forward { cur } else { 3 };
                    // skip whole words without a matching entry
                    let pat = 0x5555555555555555u64 * want as u64;
                    let x = !(w ^ pat);
                    let m = x & (x >> 1) & 0x5555555555555555;
                    if m == 0 {
                        i = (i | 31) + 1;
                        continue;
                    }
                    let mut mm = m;
                    while mm != 0 {
                        let b = mm.trailing_zeros() as usize / 2;
                        mm &= mm - 1;
                        let j = (i & !31) + b;
                        if j >= hi || j < lo {
                            continue;
                        }
                        nb.clear();
                        nbrs(j, &mut nb);
                        if forward {
                            for &k in nb.iter() {
                                if set_if_empty(k, next) {
                                    on_set(k, d + 1);
                                    cnt += 1;
                                }
                            }
                        } else if nb.iter().any(|&k| get(k) == cur) && set_if_empty(j, next) {
                            on_set(j, d + 1);
                            cnt += 1;
                        }
                    }
                    i = (i | 31) + 1;
                }
                cnt
            })
            .sum();
        if newc == 0 {
            break;
        }
        hist.push(newc);
        done += newc;
        d += 1;
        if verbose {
            eprintln!(
                "  {name}: depth {d} {newc} ({}) {:.1}s",
                if forward { "fwd" } else { "bwd" },
                t0.elapsed().as_secs_f64()
            );
        }
    }
    if verbose {
        let avg: f64 = hist
            .iter()
            .enumerate()
            .map(|(i, &c)| i as f64 * c as f64)
            .sum::<f64>()
            / done as f64;
        eprintln!(
            "{name}: {n} entries, max depth {}, mean {avg:.3}, {:.1}s",
            hist.len() - 1,
            t0.elapsed().as_secs_f64()
        );
    }
    Table2 { data: raw, n, hist }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decode_table() {
        for p in 0u8..20 {
            for dlt in [-1i32, 0, 1] {
                let c = p as i32 + dlt;
                if c < 0 {
                    continue;
                }
                assert_eq!(decode(p, (c % 3) as u8), c as u8);
            }
        }
    }

    #[test]
    fn bfs_on_cycle() {
        // ring of 1000 nodes, neighbours +-1 (padded with self loops)
        let n = 1000;
        let t = bfs(
            n,
            &[0],
            |i, out| {
                out.push((i + 1) % n);
                out.push((i + n - 1) % n);
            },
            "ring",
            false,
        );
        for i in 0..n {
            assert_eq!(t.get(i) as usize, i.min(n - i) % 3);
        }
        assert_eq!(t.hist.len(), 501);
    }
}
