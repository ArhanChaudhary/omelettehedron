//! Two-bit pruning tables storing distance mod 3 (3 = not yet visited), built by parallel BFS.
//!
//! A search tracks the exact distance of the parent; since one move changes the distance by at
//! most one, the child's exact distance follows from its value mod 3.

use crate::mem::HugeBuf;
use crate::tables::prefetch;
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

    /// Write the table to `path` (through a temporary file, renamed when complete).
    ///
    /// # Errors
    ///
    /// If the file cannot be written.
    pub fn save(&self, path: &Path) -> std::io::Result<()> {
        save_words(path, self.n, &self.hist, &self.data)
    }

    /// Read a table of `n` entries from `path`. (Parallel fused read and checksum, part of the
    /// opt/startup2 stack: large-tier startup 1.07 -> 0.33 s at 48T, 7.9 -> 3.3 s at 1T.)
    ///
    /// # Errors
    ///
    /// If the file cannot be read, or has the wrong format, size or checksum.
    pub fn load(path: &Path, n: usize) -> std::io::Result<Table2> {
        let err = |m: &str| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("{}: {m}", path.display()),
            )
        };
        let mut f = std::io::BufReader::with_capacity(1 << 16, std::fs::File::open(path)?);
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
        let offset = (32 + 8 * nh) as u64;
        let file = f.into_inner();
        // Read in parallel and checksum each piece while it is in cache; the buffer is not
        // zeroed first since every word is overwritten.
        let mut data = HugeBuf::<u64>::for_overwrite(words);
        let sums: std::io::Result<Vec<u64>> = data
            .par_chunks_mut(CK * GROUP)
            .enumerate()
            .map(|(gi, g)| read_and_hash(&file, offset, gi * GROUP, g))
            .collect();
        if sums?.into_iter().fold(0, |a, b| a ^ b) != sum {
            return Err(err("checksum mismatch"));
        }
        Ok(Table2 { data, n, hist })
    }
}

/// Write a table of `n` entries with histogram `hist` and words `data` to `path` (through a
/// temporary file, renamed when complete).
///
/// # Errors
///
/// If the file cannot be written.
pub fn save_words(path: &Path, n: usize, hist: &[u64], data: &[u64]) -> std::io::Result<()> {
    let tmp = path.with_extension("tmp");
    let mut f = std::io::BufWriter::with_capacity(1 << 24, std::fs::File::create(&tmp)?);
    f.write_all(b"FTOPRUN2")?;
    f.write_all(&(n as u64).to_le_bytes())?;
    f.write_all(&(hist.len() as u64).to_le_bytes())?;
    for h in hist {
        f.write_all(&h.to_le_bytes())?;
    }
    f.write_all(&checksum(data).to_le_bytes())?;
    // SAFETY: u64 slice reinterpreted as bytes for writing
    let bytes = unsafe { std::slice::from_raw_parts(data.as_ptr().cast::<u8>(), data.len() * 8) };
    for chunk in bytes.chunks(1 << 26) {
        f.write_all(chunk)?;
    }
    f.flush()?;
    drop(f);
    std::fs::rename(tmp, path)
}

/// Words per checksum chunk.
const CK: usize = 1 << 20;
/// Checksum chunks hashed together (independent dependency chains) by one worker.
const GROUP: usize = 4;
/// Words read per `pread`.
const PIECE: usize = 1 << 17;

/// One chunk's hash step (see [`checksum`]).
#[inline(always)]
fn hash_step(h: u64, x: u64) -> u64 {
    (h ^ x).wrapping_mul(0x0100_0000_01b3).rotate_left(23)
}

/// Read exactly `buf.len()` bytes at `offset` without moving a shared file cursor.
#[cfg(unix)]
fn read_at(file: &std::fs::File, buf: &mut [u8], offset: u64) -> std::io::Result<()> {
    use std::os::unix::fs::FileExt;
    file.read_exact_at(buf, offset)
}

/// Elsewhere (WebAssembly) there is a single thread, so seeking the shared cursor is safe.
#[cfg(not(unix))]
fn read_at(file: &std::fs::File, buf: &mut [u8], offset: u64) -> std::io::Result<()> {
    use std::io::{Read, Seek, SeekFrom};
    let mut f = file;
    f.seek(SeekFrom::Start(offset))?;
    f.read_exact(buf)
}

/// Fill `g` (chunks `c0..` of the table) from the file and return the XOR of their checksums.
fn read_and_hash(
    file: &std::fs::File,
    offset: u64,
    c0: usize,
    g: &mut [u64],
) -> std::io::Result<u64> {
    let nc = g.len().div_ceil(CK);
    let mut h: [u64; GROUP] = std::array::from_fn(|j| 0x9e37_79b9_7f4a_7c15_u64 ^ (c0 + j) as u64);
    let full = g.len() / CK == GROUP;
    let mut pos = 0;
    while pos < CK.min(g.len()) {
        let len = PIECE.min(CK - pos);
        for j in 0..nc {
            let lo = (j * CK + pos).min(g.len());
            let hi = (j * CK + pos + len).min(g.len()).min((j + 1) * CK);
            if lo >= hi {
                continue;
            }
            // SAFETY: reading raw bytes into the u64 buffer
            let bytes = unsafe {
                std::slice::from_raw_parts_mut(g[lo..hi].as_mut_ptr().cast::<u8>(), (hi - lo) * 8)
            };
            read_at(file, bytes, offset + ((c0 * CK + lo) * 8) as u64)?;
        }
        if full {
            // four independent chains: the multiply latency is shared
            for i in pos..pos + len {
                for j in 0..GROUP {
                    h[j] = hash_step(h[j], g[j * CK + i]);
                }
            }
        } else {
            for j in 0..nc {
                let lo = (j * CK + pos).min(g.len());
                let hi = (j * CK + pos + len).min(g.len()).min((j + 1) * CK);
                for &x in &g[lo..hi.max(lo)] {
                    h[j] = hash_step(h[j], x);
                }
            }
        }
        pos += len;
    }
    Ok(h[..nc].iter().fold(0, |a, b| a ^ b))
}

fn checksum(d: &[u64]) -> u64 {
    d.par_chunks(1 << 20)
        .enumerate()
        .map(|(i, c)| {
            let mut h = 0x9e37_79b9_7f4a_7c15_u64 ^ i as u64;
            for &x in c {
                h = (h ^ x).wrapping_mul(0x0100_0000_01b3).rotate_left(23);
            }
            h
        })
        .reduce(|| 0, |a, b| a ^ b) // associative: rayon may combine in any order
}

/// The graph explored by [`Bfs`]: table entries and the moves between them.
pub trait Graph: Sync {
    /// Moves (neighbours) per entry.
    const MOVES: usize;
    /// An entry decoded into whatever [`Graph::step`] needs.
    type Node: Copy;
    /// Per-thread state kept between [`Graph::node`] calls, e.g. data shared by consecutive
    /// entries (entries are decoded mostly in increasing order).
    type Cache: Default;
    /// Decode entry `i`.
    fn node(&self, i: usize, cache: &mut Self::Cache) -> Self::Node;
    /// Index of the neighbour of `x` by move `m` (`m < MOVES`).
    fn step(&self, x: &Self::Node, m: usize) -> usize;
    /// Append every other index that denotes a position equivalent to index `j`, if indices
    /// are not unique. Equivalent indices always get the same value.
    fn equivalents(&self, _j: usize, _out: &mut Vec<usize>) {}
    /// Call `out` with the neighbour of entry `i` by each move (what a forward step needs).
    /// A graph can override it to share work between consecutive entries.
    #[inline(always)]
    fn neighbours(&self, i: usize, cache: &mut Self::Cache, mut out: impl FnMut(usize)) {
        let x = self.node(i, cache);
        for m in 0..Self::MOVES {
            out(self.step(&x, m));
        }
    }
}

/// A level runs forward if frontier * `MOVES` < unvisited * `FORWARD_RATIO`. A forward step
/// makes all neighbours of each frontier entry; a backward step makes those of each unvisited
/// entry up to the first frontier hit (6.3 per entry on average at CBE depth 8). The rule is
/// egg-solver's, re-checked on the levels nearest the boundary (48 threads, opt/gen): AE
/// medium depth 11 (ratio 3.4) forward 4.7 s vs backward 5.4 s, CBE depth 8 (ratio 6.6)
/// backward 15.3 s vs forward 20.6 s; small-tier AE depth 10 (2.4) and CB depth 6 (3.0)
/// differ by at most 0.1 s either way.
const FORWARD_RATIO: u64 = 4;

/// Indices per work item (a multiple of 32, so that work items own whole words).
const CH: usize = 1 << 14;
/// Entries checked in an interleaved way by one thread in a backward step.
const LANES: usize = 16;
/// Repeats the 2-bit pattern 01 over a word (the low bit of every entry).
const LOW: u64 = 0x5555_5555_5555_5555;

#[inline(always)]
fn get(tab: &[AtomicU64], i: usize) -> u8 {
    ((tab[i >> 5].load(Ordering::Relaxed) >> ((i & 31) << 1)) & 3) as u8
}

/// Set entry `i` to `v` if it is unvisited (3). Returns whether this call set it.
///
/// All writers of a level write the same value, so a single `fetch_and` (3 -> v clears the
/// bits of 3 ^ v) replaces a compare-and-swap loop; the plain load first avoids taking the
/// cache line exclusive for the many entries that are already visited. (Alone it measured
/// no faster, AE medium depth 11 7.8 -> 8.1 s at 48T; see `RING`.)
#[inline(always)]
fn set_if_empty(tab: &[AtomicU64], i: usize, v: u8) -> bool {
    let w = &tab[i >> 5];
    let sh = (i & 31) << 1;
    if (w.load(Ordering::Relaxed) >> sh) & 3 != 3 {
        return false;
    }
    let prev = w.fetch_and(!(u64::from(3 ^ v) << sh), Ordering::Relaxed);
    (prev >> sh) & 3 == 3
}

/// Hint that `p` will be written soon.
#[inline(always)]
fn prefetch_write<T>(p: &T) {
    #[cfg(target_arch = "x86_64")]
    // SAFETY: prefetch is a hint and never faults
    unsafe {
        core::arch::x86_64::_mm_prefetch::<{ core::arch::x86_64::_MM_HINT_ET0 }>(
            std::ptr::from_ref(p).cast::<i8>(),
        );
    }
    #[cfg(target_arch = "aarch64")]
    // SAFETY: PRFM is a hint and never faults, whatever the address
    unsafe {
        core::arch::asm!(
            "prfm pstl1keep, [{0}]",
            in(reg) std::ptr::from_ref(p),
            options(nostack, preserves_flags)
        );
    }
    #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
    let _ = p;
}

/// The indices in `lo..hi` whose entry is `code` (`lo` a multiple of 32), each word read when
/// the scan reaches it.
struct Scan<'a> {
    tab: &'a [AtomicU64],
    /// next word to read, and the end
    w: usize,
    end: usize,
    hi: usize,
    pat: u64,
    /// matching entries of the current word (one bit each, at even positions), and its base
    bits: u64,
    base: usize,
}

impl<'a> Scan<'a> {
    fn new(tab: &'a [AtomicU64], lo: usize, hi: usize, code: u8) -> Scan<'a> {
        Scan {
            tab,
            w: lo >> 5,
            end: hi.div_ceil(32),
            hi,
            pat: LOW * code as u64,
            bits: 0,
            base: 0,
        }
    }
}

impl Iterator for Scan<'_> {
    type Item = usize;
    #[inline(always)]
    fn next(&mut self) -> Option<usize> {
        loop {
            if self.bits != 0 {
                let i = self.base + self.bits.trailing_zeros() as usize / 2;
                self.bits &= self.bits - 1;
                if i < self.hi {
                    return Some(i);
                }
                self.bits = 0;
            }
            if self.w >= self.end {
                return None;
            }
            let x = !(self.tab[self.w].load(Ordering::Relaxed) ^ self.pat);
            self.bits = x & (x >> 1) & LOW;
            self.base = self.w << 5;
            self.w += 1;
        }
    }
}

/// Pending writes of a forward step: each index is prefetched when it is pushed and written
/// `RING` pushes later, so that about `RING` table misses are in flight. (With the single
/// `fetch_and`: AE medium depth 11 7.8 -> 6.9 s, medium tier 11.4 -> 10.2 s at 48T when added,
/// opt/gen.)
const RING: usize = 32;

/// Forward step over `lo..hi`: set every unvisited neighbour (and its equivalents) of the
/// entries holding `cur` to `next`. Returns the number of entries set.
fn forward<G: Graph>(g: &G, tab: &[AtomicU64], lo: usize, hi: usize, cur: u8, next: u8) -> u64 {
    let mut ring = [usize::MAX; RING];
    let mut head = 0;
    let mut cnt = 0;
    let write = |k: usize, cnt: &mut u64| {
        if set_if_empty(tab, k, next) {
            *cnt += 1;
        }
    };
    let mut push = |k: usize, cnt: &mut u64| {
        prefetch_write(&tab[k >> 5]);
        let old = std::mem::replace(&mut ring[head], k);
        head = (head + 1) % RING;
        if old != usize::MAX {
            write(old, cnt);
        }
    };
    let mut nb: Vec<usize> = Vec::with_capacity(16);
    let mut cache = G::Cache::default();
    for i in Scan::new(tab, lo, hi, cur) {
        g.neighbours(i, &mut cache, |j| {
            push(j, &mut cnt);
            nb.clear();
            g.equivalents(j, &mut nb);
            for &k in &nb {
                push(k, &mut cnt);
            }
        });
    }
    for k in ring {
        if k != usize::MAX {
            write(k, &mut cnt);
        }
    }
    cnt
}

/// One entry of a backward step in progress: neighbour `j` (by move `m`) of entry `i` has
/// been prefetched.
#[derive(Clone, Copy)]
struct Lane<N> {
    i: usize,
    j: usize,
    m: usize,
    x: N,
}

/// Backward step over `lo..hi`: every unvisited entry with a neighbour holding `cur` is set
/// to `next` (if `set`) and counted; the others are appended to `misses` (if given).
///
/// The neighbours of an entry are made one move at a time and the search stops at the first
/// hit; `LANES` entries are interleaved, each neighbour's word prefetched one round before it
/// is read, so that several misses are in flight. Equivalent indices are not looked up: they
/// hold the same value. (CBE generation 102.7 -> 51.0 s, large tier 121 -> 65 s at 48T when
/// added, opt/gen.)
fn backward<G: Graph>(
    g: &G,
    tab: &[AtomicU64],
    lo: usize,
    hi: usize,
    (cur, next): (u8, u8),
    set: bool,
    mut misses: Option<&mut Vec<usize>>,
) -> u64 {
    let mut cache = G::Cache::default();
    let mut start = |i: usize| -> Lane<G::Node> {
        let x = g.node(i, &mut cache);
        let j = g.step(&x, 0);
        prefetch(&tab[j >> 5]);
        Lane { i, j, m: 0, x }
    };
    let mut scan = Scan::new(tab, lo, hi, 3);
    let mut lanes: Vec<Lane<G::Node>> = scan.by_ref().take(LANES).map(&mut start).collect();
    let mut cnt = 0;
    let mut k = 0;
    while !lanes.is_empty() {
        let l = &mut lanes[k];
        let hit = get(tab, l.j) == cur;
        if hit || l.m + 1 == G::MOVES {
            if hit {
                cnt += 1;
                if set {
                    // Only this work item writes the words of `lo..hi` in a backward step, so
                    // a plain read-modify-write is enough.
                    let w = &tab[l.i >> 5];
                    let sh = (l.i & 31) << 1;
                    let v = w.load(Ordering::Relaxed);
                    w.store(v & !(u64::from(3 ^ next) << sh), Ordering::Relaxed);
                }
            } else if let Some(v) = misses.as_deref_mut() {
                v.push(l.i);
            }
            if let Some(i) = scan.next() {
                *l = start(i);
                k += 1;
            } else {
                lanes.swap_remove(k);
            }
        } else {
            l.m += 1;
            l.j = g.step(&l.x, l.m);
            prefetch(&tab[l.j >> 5]);
            k += 1;
        }
        if k >= lanes.len() {
            k = 0;
        }
    }
    cnt
}

/// Map every 2-bit entry `v` of a word to `map[v]`.
#[inline(always)]
fn remap_word(w: u64, map: [u8; 4]) -> u64 {
    let (lo, hi) = (w & LOW, (w >> 1) & LOW);
    let is = [!hi & !lo & LOW, !hi & lo, hi & !lo, hi & lo];
    let (mut nlo, mut nhi) = (0, 0);
    for v in 0..4 {
        if map[v] & 1 != 0 {
            nlo |= is[v];
        }
        if map[v] & 2 != 0 {
            nhi |= is[v];
        }
    }
    nlo | nhi << 1
}

/// Level-by-level breadth-first search over the 2-bit entries of a table; 3 marks an
/// unvisited entry. Each level reads the entries holding one code (the frontier) and writes
/// another, so a caller can keep the distance mod 3 ([`bfs`]) or saturating codes.
pub struct Bfs<'g, G: Graph> {
    g: &'g G,
    n: usize,
    raw: HugeBuf<u64>,
    /// number of entries at each distance
    pub hist: Vec<u64>,
    done: u64,
    name: &'g str,
    verbose: bool,
    t0: std::time::Instant,
}

impl<'g, G: Graph> Bfs<'g, G> {
    /// All `n` entries unvisited except `targets` (distance 0), which hold `code`.
    pub fn new(
        g: &'g G,
        n: usize,
        targets: &[usize],
        code: u8,
        name: &'g str,
        verbose: bool,
    ) -> Self {
        let t0 = std::time::Instant::now();
        let mut raw = HugeBuf::<u64>::for_overwrite(n.div_ceil(32));
        raw.par_chunks_mut(1 << 16).for_each(|c| c.fill(u64::MAX));
        let mut s = Bfs {
            g,
            n,
            raw,
            hist: vec![],
            done: 0,
            name,
            verbose,
            t0,
        };
        let tab = s.tab();
        let nt = targets
            .iter()
            .filter(|&&t| set_if_empty(tab, t, code))
            .count() as u64;
        s.hist.push(nt);
        s.done = nt;
        s
    }

    fn tab(&self) -> &[AtomicU64] {
        // SAFETY: AtomicU64 has the same size and alignment as u64, and while the returned
        // slice borrows `self`, the words are only accessed through it
        unsafe { std::slice::from_raw_parts(self.raw.as_ptr().cast::<AtomicU64>(), self.raw.len()) }
    }

    /// Entries not reached yet.
    pub fn unvisited(&self) -> u64 {
        self.n as u64 - self.done
    }

    /// One level: every unvisited entry next to an entry holding `cur` gets `next`. Returns
    /// the number of entries reached (also appended to `hist` if non-zero). Older entries may
    /// hold `cur` too (mod-3 distances): a forward step expands them for nothing, a backward
    /// step cannot meet them.
    pub fn level(&mut self, cur: u8, next: u8) -> u64 {
        let tab = self.tab();
        let (g, n) = (self.g, self.n);
        let frontier = self.hist.last().copied().unwrap_or(0);
        let fwd = frontier.saturating_mul(G::MOVES as u64)
            < self.unvisited().saturating_mul(FORWARD_RATIO);
        let newc: u64 = (0..n.div_ceil(CH))
            .into_par_iter()
            .map(|c| {
                let (lo, hi) = (c * CH, ((c + 1) * CH).min(n));
                if fwd {
                    forward(g, tab, lo, hi, cur, next)
                } else {
                    backward(g, tab, lo, hi, (cur, next), true, None)
                }
            })
            .sum();
        if newc > 0 {
            self.hist.push(newc);
            self.done += newc;
        }
        if self.verbose {
            eprintln!(
                "  {}: depth {} {newc} ({}) {:.1}s",
                self.name,
                self.hist.len() - usize::from(newc > 0),
                if fwd { "fwd" } else { "bwd" },
                self.t0.elapsed().as_secs_f64()
            );
        }
        newc
    }

    /// Replace every entry value `v` by `map[v]`.
    pub fn remap(&mut self, map: [u8; 4]) {
        self.raw
            .par_chunks_mut(1 << 16)
            .for_each(|c| c.iter_mut().for_each(|w| *w = remap_word(*w, map)));
    }

    /// Append the distances of the entries not reached yet to `hist`, without changing the
    /// table: a backward pass counts the unvisited entries next to one holding `cur` (the
    /// last level), and a BFS restricted to the entries it does not reach (meant to be few)
    /// counts the rest.
    pub fn count_rest(&mut self, cur: u8) {
        let tab = self.tab();
        let (g, n) = (self.g, self.n);
        let (cnt, mut rest) = (0..n.div_ceil(CH))
            .into_par_iter()
            .map(|c| {
                let mut v = vec![];
                let lo = c * CH;
                let k = backward(
                    g,
                    tab,
                    lo,
                    (lo + CH).min(n),
                    (cur, cur),
                    false,
                    Some(&mut v),
                );
                (k, v)
            })
            .reduce(
                || (0, vec![]),
                |(a, mut va), (b, vb)| {
                    va.extend(vb);
                    (a + b, va)
                },
            );
        if cnt == 0 {
            return;
        }
        let mut counts = vec![cnt];
        // Distances among the rest: the previous level is the unvisited entries outside the
        // rest at first, then the entries found last. (Equivalent indices are in the same
        // sets, so neighbours are made without them.)
        rest.par_sort_unstable();
        let mut last: Option<Vec<usize>> = None;
        loop {
            let in_last = |j: usize| match &last {
                None => get(tab, j) == 3 && rest.binary_search(&j).is_err(),
                Some(l) => l.binary_search(&j).is_ok(),
            };
            let (mut found, left): (Vec<usize>, Vec<usize>) = rest.par_iter().partition(|&&i| {
                let x = g.node(i, &mut G::Cache::default());
                (0..G::MOVES).any(|m| in_last(g.step(&x, m)))
            });
            if found.is_empty() {
                break;
            }
            counts.push(found.len() as u64);
            found.par_sort_unstable();
            last = Some(found);
            rest = left;
        }
        if self.verbose {
            eprintln!(
                "  {}: depths {}.. {counts:?} (counted) {:.1}s",
                self.name,
                self.hist.len(),
                self.t0.elapsed().as_secs_f64()
            );
        }
        self.done += counts.iter().sum::<u64>();
        self.hist.extend(counts);
    }

    /// The table with the distance histogram.
    pub fn finish(self) -> Table2 {
        if self.verbose {
            let avg: f64 = self
                .hist
                .iter()
                .enumerate()
                .map(|(i, &c)| i as f64 * c as f64)
                .sum::<f64>()
                / self.done as f64;
            eprintln!(
                "{}: {} entries, max depth {}, mean {avg:.3}, {:.1}s",
                self.name,
                self.n,
                self.hist.len() - 1,
                self.t0.elapsed().as_secs_f64()
            );
        }
        Table2 {
            data: self.raw,
            n: self.n,
            hist: self.hist,
        }
    }
}

/// Saturating distance codes from `targets`: 0 = at most `base`, 1 = `base + 1`,
/// 2 = `base + 2`, 3 = more (or unreachable). The histogram covers all distances.
///
/// Distances are kept mod 3 up to `base`; then the frontier becomes code 1 and everything
/// older code 0, the next level writes 2 and is relabelled to 1, and the BFS stops after
/// `base + 2`, where the codes are final. Only the histogram of the rest is counted. (Used
/// for CBE: no second array and no copies when building or saving; peak memory of the
/// large-tier generation 14.9 -> 8.1 GB, CBE generation 47.3 -> 44.8 s at 48T when added,
/// opt/gen.)
pub fn saturating<G: Graph>(
    g: &G,
    n: usize,
    targets: &[usize],
    base: u8,
    name: &str,
    verbose: bool,
) -> Table2 {
    let mut b = Bfs::new(g, n, targets, 0, name, verbose);
    for d in 0..base {
        b.level(d % 3, (d + 1) % 3);
    }
    // Distance `base` -> 1, the rest of <= base -> 0. (Distances base - 3, base - 6, ...
    // share the value of `base` and become frontier entries too; that is harmless: all their
    // neighbours are reached.)
    let mut map = [0, 0, 0, 3];
    map[(base % 3) as usize] = 1;
    b.remap(map);
    b.level(1, 2);
    b.remap([0, 0, 1, 3]);
    b.level(1, 2);
    b.count_rest(2);
    b.finish()
}

/// Distance mod 3 of every entry from `targets` (unreachable entries keep 3).
pub fn bfs<G: Graph>(g: &G, n: usize, targets: &[usize], name: &str, verbose: bool) -> Table2 {
    let mut b = Bfs::new(g, n, targets, 0, name, verbose);
    let mut cur = 0;
    while b.unvisited() > 0 && b.level(cur, (cur + 1) % 3) > 0 {
        cur = (cur + 1) % 3;
    }
    b.finish()
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
                let c = u8::try_from(c).unwrap();
                assert_eq!(decode(p, c % 3), c);
            }
        }
    }

    /// Ring of `n` nodes, neighbours +-1.
    struct Ring(usize);

    impl Graph for Ring {
        const MOVES: usize = 2;
        type Node = usize;
        type Cache = ();
        fn node(&self, i: usize, (): &mut ()) -> usize {
            i
        }
        fn step(&self, &i: &usize, m: usize) -> usize {
            if m == 0 {
                (i + 1) % self.0
            } else {
                (i + self.0 - 1) % self.0
            }
        }
    }

    #[test]
    fn bfs_on_cycle() {
        let n = 1000;
        let t = bfs(&Ring(n), n, &[0], "ring", false);
        for i in 0..n {
            assert_eq!(t.get(i) as usize, i.min(n - i) % 3);
        }
        assert_eq!(t.hist.len(), 501);
    }

    /// Circulant graph (steps +-1, +-97, +-9973) on `n` nodes, with `i` and `n - i`
    /// equivalent: several work items, a partial last word, and equivalent indices.
    struct Chords(usize);

    impl Graph for Chords {
        const MOVES: usize = 6;
        type Node = usize;
        type Cache = ();
        fn node(&self, i: usize, (): &mut ()) -> usize {
            i
        }
        fn step(&self, &i: &usize, m: usize) -> usize {
            let d = [1, 97, 9973][m / 2];
            if m.is_multiple_of(2) {
                (i + d) % self.0
            } else {
                (i + self.0 - d) % self.0
            }
        }
        fn equivalents(&self, j: usize, out: &mut Vec<usize>) {
            if j != 0 {
                out.push(self.0 - j);
            }
        }
    }

    /// Exact distances from node 0 and their histogram.
    fn reference(g: &Chords) -> (Vec<u32>, Vec<u64>) {
        let n = g.0;
        let mut dist = vec![u32::MAX; n];
        dist[0] = 0;
        let mut q = std::collections::VecDeque::from([0]);
        while let Some(i) = q.pop_front() {
            for m in 0..Chords::MOVES {
                let j = g.step(&i, m);
                if dist[j] == u32::MAX {
                    dist[j] = dist[i] + 1;
                    q.push_back(j);
                }
            }
        }
        let mut hist = vec![0u64; *dist.iter().max().unwrap() as usize + 1];
        for &d in &dist {
            hist[d as usize] += 1;
        }
        (dist, hist)
    }

    #[test]
    fn bfs_matches_reference() {
        let g = Chords(300_007);
        let (dist, hist) = reference(&g);
        let t = bfs(&g, g.0, &[0], "chords", false);
        assert_eq!(t.hist, hist);
        for i in 0..g.0 {
            assert_eq!(t.get(i) as u32, dist[i] % 3);
        }
    }

    #[test]
    fn saturating_matches_reference() {
        let g = Chords(300_007);
        let (dist, hist) = reference(&g);
        for base in [2u8, 6, 7] {
            let t = saturating(&g, g.0, &[0], base, "chords", false);
            assert_eq!(t.hist, hist);
            for i in 0..g.0 {
                let code = dist[i].saturating_sub(base as u32).min(3);
                assert_eq!(t.get(i) as u32, code, "base {base} index {i}");
            }
        }
    }
}
