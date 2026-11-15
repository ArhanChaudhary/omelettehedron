//! Breadth-first search straight into a packed table ([`Table5`]: distances mod 3, five to a
//! byte), for tables whose 2-bit form ([`crate::prune::Bfs`], which marks "not yet visited"
//! with a fourth value) would not fit in memory: the giant tier's edge x A-centre table takes
//! 83.6 GB at 2 bits and 66.9 GB packed.
//!
//! A packed entry has three values, so "not visited" has none of its own. The unvisited
//! entries share one with a class of visited ones, chosen per level so that the steps can
//! still tell them apart. Level `L` is made in one of three ways:
//! * **Forward** (small frontiers): the unvisited entries hold `(L + 1) mod 3`, and so do the
//!   levels `L - 2`, `L - 5`, ...; those are moved to `(L - 1) mod 3` first (their entries are
//!   listed), so an entry holding `(L + 1) mod 3` is unvisited, and every one next to the
//!   frontier (level `L - 1`) is set to `L mod 3`. Then the unvisited entries move on to
//!   `(L + 2) mod 3` (a pass over the table) and the moved levels go back. (When the frontier
//!   is found by a scan, the older levels holding its value and the moved ones are expanded
//!   too: their neighbours are visited and none holds `(L + 1) mod 3`, so nothing changes.)
//! * **Backward** (large levels): the unvisited entries hold `L mod 3`, and so do the levels
//!   `L - 3`, `L - 6`, ... Every entry holding it (the targets aside) is checked for a
//!   neighbour holding `(L - 1) mod 3`. An unvisited entry has one exactly if it is at level
//!   `L` (its neighbours are at `L - 1` or more), a visited one at level `k >= 1` always has
//!   one (a neighbour at level `k - 1`): both keep their value. The others are unvisited and
//!   move on to `(L + 1) mod 3`. The level's size is the hits less the sizes of the levels
//!   that share the value.
//! * **Listed** (few unvisited left): the entries a backward step left unvisited were listed,
//!   and only they are checked, with no visited entries among them.
//!
//! Lists of levels (for forward steps) and of unvisited entries (for listed steps) are kept in
//! memory or in files of a scratch directory; a list that would exceed [`Bfs5Config::list_cap`]
//! is dropped, and the steps that would need it are made backward instead.

use crate::mem::HugeBuf;
use crate::prune::{prefetch_write, Graph, Table5, FORWARD_RATIO, PER5, TRITS};
use crate::tables::prefetch;
use rayon::prelude::*;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::sync::Mutex;

/// Settings of [`bfs5`].
pub struct Bfs5Config<'a> {
    /// name in messages
    pub name: &'a str,
    pub verbose: bool,
    /// directory for the lists (`None`: kept in memory)
    pub scratch: Option<&'a Path>,
    /// the most entries a list may hold
    pub list_cap: usize,
    /// a level is made forward while its frontier times the moves is below the unvisited
    /// entries times this (and the lists it needs are there)
    pub forward_ratio: u64,
    /// set to give up (the search returns `None`)
    pub stop: &'a AtomicBool,
}

impl<'a> Bfs5Config<'a> {
    /// Lists in memory, up to 2^32 entries; the 2-bit search's forward ratio.
    pub fn new(name: &'a str, verbose: bool, stop: &'a AtomicBool) -> Bfs5Config<'a> {
        Bfs5Config {
            name,
            verbose,
            scratch: None,
            list_cap: 1 << 32,
            forward_ratio: FORWARD_RATIO,
            stop,
        }
    }
}

/// Powers of 3: the weight of each digit of a byte.
const POW3: [u64; 8] = [1, 3, 9, 27, 81, 0, 0, 0];

/// `MATCH[v][b]`: the digits of byte `b` that are `v` (one bit each).
static MATCH: [[u8; 256]; 3] = {
    let mut t = [[0u8; 256]; 3];
    let mut v = 0;
    while v < 3 {
        let mut b = 0;
        while b < 243 {
            let mut x = b;
            let mut j = 0;
            while j < 5 {
                if x % 3 == v {
                    t[v][b] |= 1 << j;
                }
                x /= 3;
                j += 1;
            }
            b += 1;
        }
        v += 1;
    }
    t
};

/// A word with every digit `v`.
const fn fill_word(v: u8) -> u64 {
    let b = (v as u64) * (1 + 3 + 9 + 27 + 81);
    b * 0x0101_0101_0101_0101
}

#[inline(always)]
fn tget(tab: &[AtomicU64], i: usize) -> u8 {
    let q = i / 5;
    let b = (tab[q >> 3].load(Ordering::Relaxed) >> ((q & 7) << 3)) as u8;
    TRITS[usize::from(b)][(i - q * 5) & 7]
}

/// The word of entry `i` and the amount to add to it to turn the entry's digit from `from` to
/// `to` (as a wrapping addition: a digit changes within 0..3 and a byte within 0..243, so
/// nothing carries into the neighbouring digits or bytes).
#[inline(always)]
fn tdelta(i: usize, from: u8, to: u8) -> (usize, u64) {
    let q = i / 5;
    let d = POW3[(i - q * 5) & 7] << ((q & 7) << 3);
    (
        q >> 3,
        u64::from(to).wrapping_sub(u64::from(from)).wrapping_mul(d),
    )
}

/// Set entry `i` from `from` to `to` if it holds `from`; whether this call set it.
#[inline(always)]
fn tcas(tab: &[AtomicU64], i: usize, from: u8, to: u8) -> bool {
    let q = i / 5;
    let (w, delta) = tdelta(i, from, to);
    let (sh, j) = ((q & 7) << 3, (i - q * 5) & 7);
    let word = &tab[w];
    let mut old = word.load(Ordering::Relaxed);
    loop {
        if TRITS[usize::from((old >> sh) as u8)][j] != from {
            return false;
        }
        match word.compare_exchange_weak(
            old,
            old.wrapping_add(delta),
            Ordering::Relaxed,
            Ordering::Relaxed,
        ) {
            Ok(_) => return true,
            Err(x) => old = x,
        }
    }
}

/// Change entry `i` from `from` to `to`, knowing it holds `from` and that no other thread
/// changes it (other digits of its word may change meanwhile).
#[inline(always)]
fn tmove(tab: &[AtomicU64], i: usize, from: u8, to: u8) {
    let (w, delta) = tdelta(i, from, to);
    tab[w].fetch_add(delta, Ordering::Relaxed);
}

/// Words per work item of a scan (`PER5` entries each).
const WCH: usize = 1 << 12;
/// Entries buffered per thread before they are appended to a list.
const BUF: usize = 1 << 16;
/// Entries read from a list file at a time.
const READ: usize = 1 << 22;
/// Pending writes of a forward step (as `crate::prune`'s forward step).
const RING: usize = 32;
/// Entries checked in an interleaved way by one thread in a backward step.
const LANES: usize = 16;

/// A list of entries being written: in memory or appended to a file, dropped (and its file
/// removed) once it would exceed the cap.
struct ListWriter {
    mem: Mutex<Vec<u64>>,
    file: Option<(PathBuf, Mutex<std::io::BufWriter<std::fs::File>>)>,
    len: AtomicUsize,
    cap: usize,
    over: AtomicBool,
    err: Mutex<Option<std::io::Error>>,
}

/// A complete list.
enum List {
    Mem(Vec<u64>),
    File(PathBuf, usize),
}

impl List {
    /// Call `f` on the entries a block at a time (the whole list in memory; blocks of
    /// `READ` entries read in turn from a file).
    fn blocks(&self, mut f: impl FnMut(&[u64])) -> std::io::Result<()> {
        match self {
            List::Mem(v) => {
                f(v);
                Ok(())
            }
            List::File(path, n) => {
                let mut r = std::io::BufReader::with_capacity(1 << 24, std::fs::File::open(path)?);
                let mut left = *n;
                let mut bytes = vec![0u8; READ * 8];
                let mut buf = vec![0u64; READ];
                while left > 0 {
                    let k = left.min(READ);
                    r.read_exact(&mut bytes[..k * 8])?;
                    let (words, _) = bytes[..k * 8].as_chunks::<8>();
                    for (x, e) in buf.iter_mut().zip(words) {
                        *x = u64::from_le_bytes(*e);
                    }
                    f(&buf[..k]);
                    left -= k;
                }
                Ok(())
            }
        }
    }

    /// Call `f` on every entry, in parallel.
    fn for_each(&self, f: impl Fn(usize) + Sync) -> std::io::Result<()> {
        self.blocks(|b| b.par_iter().for_each(|&i| f(i as usize)))
    }

    /// Call `f` on runs of entries, in parallel.
    fn for_runs(&self, f: impl Fn(&mut dyn Iterator<Item = usize>) + Sync) -> std::io::Result<()> {
        self.blocks(|b| {
            b.par_chunks(256)
                .for_each(|c| f(&mut c.iter().map(|&i| i as usize)));
        })
    }

    fn remove(self) {
        if let List::File(path, _) = self {
            let _ = std::fs::remove_file(path);
        }
    }
}

impl ListWriter {
    fn new(scratch: Option<&Path>, name: &str, cap: usize) -> std::io::Result<ListWriter> {
        let file = match scratch {
            Some(dir) => {
                let path = dir.join(name);
                let f = std::fs::File::create(&path)?;
                Some((
                    path,
                    Mutex::new(std::io::BufWriter::with_capacity(1 << 24, f)),
                ))
            }
            None => None,
        };
        Ok(ListWriter {
            mem: Mutex::new(Vec::new()),
            file,
            len: AtomicUsize::new(0),
            cap,
            over: AtomicBool::new(false),
            err: Mutex::new(None),
        })
    }

    /// Append a thread's buffer (and empty it).
    fn flush(&self, buf: &mut Vec<u64>) {
        if buf.is_empty() {
            return;
        }
        let n = buf.len();
        if self.over.load(Ordering::Relaxed)
            || self.len.fetch_add(n, Ordering::Relaxed) + n > self.cap
        {
            self.over.store(true, Ordering::Relaxed);
            buf.clear();
            return;
        }
        match &self.file {
            None => self
                .mem
                .lock()
                .unwrap_or_else(std::sync::PoisonError::into_inner)
                .extend_from_slice(buf),
            Some((_, f)) => {
                let bytes: Vec<u8> = buf.iter().flat_map(|x| x.to_le_bytes()).collect();
                let r = f
                    .lock()
                    .unwrap_or_else(std::sync::PoisonError::into_inner)
                    .write_all(&bytes);
                if let Err(e) = r {
                    self.over.store(true, Ordering::Relaxed);
                    *self
                        .err
                        .lock()
                        .unwrap_or_else(std::sync::PoisonError::into_inner) = Some(e);
                }
            }
        }
        buf.clear();
    }

    /// The list, or `None` if it was dropped.
    fn finish(self) -> std::io::Result<Option<List>> {
        if let Some(e) = self
            .err
            .into_inner()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
        {
            return Err(e);
        }
        let over = self.over.load(Ordering::Relaxed);
        let len = self.len.load(Ordering::Relaxed);
        match self.file {
            None => Ok((!over).then(|| {
                List::Mem(
                    self.mem
                        .into_inner()
                        .unwrap_or_else(std::sync::PoisonError::into_inner),
                )
            })),
            Some((path, f)) => {
                let mut f = f
                    .into_inner()
                    .unwrap_or_else(std::sync::PoisonError::into_inner);
                f.flush()?;
                drop(f);
                if over {
                    std::fs::remove_file(&path)?;
                    Ok(None)
                } else {
                    Ok(Some(List::File(path, len)))
                }
            }
        }
    }
}

/// A thread's buffer for a list writer (flushed when full and when dropped).
struct Appender<'a> {
    w: Option<&'a ListWriter>,
    buf: Vec<u64>,
}

impl Appender<'_> {
    #[inline(always)]
    fn push(&mut self, i: usize) {
        if let Some(w) = self.w {
            self.buf.push(i as u64);
            if self.buf.len() >= BUF {
                w.flush(&mut self.buf);
            }
        }
    }
}

impl Drop for Appender<'_> {
    fn drop(&mut self) {
        if let Some(w) = self.w {
            w.flush(&mut self.buf);
        }
    }
}

fn appender(w: Option<&ListWriter>) -> Appender<'_> {
    Appender {
        w,
        buf: Vec::with_capacity(if w.is_some() { BUF } else { 0 }),
    }
}

/// The entries below `n` holding `v` in words `w..end`, each word read when the scan reaches
/// it.
struct Scan<'a> {
    tab: &'a [AtomicU64],
    w: usize,
    end: usize,
    n: usize,
    v: usize,
    /// the current word, its next byte, and the matching digits of the current byte
    x: u64,
    k: usize,
    m: u8,
    base: usize,
}

impl Iterator for Scan<'_> {
    type Item = usize;
    #[inline(always)]
    fn next(&mut self) -> Option<usize> {
        loop {
            if self.m != 0 {
                let i = self.base + self.m.trailing_zeros() as usize;
                self.m &= self.m - 1;
                return (i < self.n).then_some(i);
            }
            if self.k == 8 {
                if self.w >= self.end {
                    return None;
                }
                self.x = self.tab[self.w].load(Ordering::Relaxed);
                self.k = 0;
                self.w += 1;
            }
            self.m = MATCH[self.v][usize::from((self.x >> (8 * self.k)) as u8)];
            self.base = ((self.w - 1) * 8 + self.k) * 5;
            self.k += 1;
        }
    }
}

/// Call `f` on the entries below `n` holding `v`, in parallel (work items of whole words).
fn scan(tab: &[AtomicU64], n: usize, v: u8, f: impl Fn(&mut dyn Iterator<Item = usize>) + Sync) {
    let words = n.div_ceil(PER5);
    (0..words.div_ceil(WCH)).into_par_iter().for_each(|c| {
        let lo = c * WCH;
        let mut it = Scan {
            tab,
            w: lo,
            end: (lo + WCH).min(words),
            n,
            v: usize::from(v),
            x: 0,
            k: 8,
            m: 0,
            base: 0,
        };
        f(&mut it);
    });
}

/// Turn every digit `from` into `to`, over the whole table.
fn rotate(tab: &[AtomicU64], from: u8, to: u8) {
    let mut lut = [0u8; 256];
    for (b, x) in lut.iter_mut().enumerate() {
        let m = MATCH[usize::from(from)][b];
        let mut nb = b as u64;
        for j in 0..5 {
            if m >> j & 1 != 0 {
                nb = nb + POW3[j] * u64::from(to) - POW3[j] * u64::from(from);
            }
        }
        *x = nb as u8;
    }
    tab.par_chunks(WCH).for_each(|c| {
        for w in c {
            let x = w.load(Ordering::Relaxed);
            let mut y = 0u64;
            for k in 0..8 {
                y |= u64::from(lut[usize::from((x >> (8 * k)) as u8)]) << (8 * k);
            }
            w.store(y, Ordering::Relaxed);
        }
    });
}

/// Forward expansion of `frontier` (calls `push` with each entry to expand): every neighbour
/// (and its equivalents) holding `u` is set to `v` and appended to `out`. Returns how many.
fn expand<G: Graph>(
    g: &G,
    tab: &[AtomicU64],
    (u, v): (u8, u8),
    frontier: &mut dyn Iterator<Item = usize>,
    out: &mut Appender,
) -> u64 {
    let mut ring = [usize::MAX; RING];
    let mut head = 0;
    let mut cnt = 0;
    let mut nb: Vec<usize> = Vec::with_capacity(16);
    let mut cache = G::Cache::default();
    let write = |k: usize, cnt: &mut u64, out: &mut Appender, nb: &mut Vec<usize>| {
        if tcas(tab, k, u, v) {
            *cnt += 1;
            out.push(k);
            nb.clear();
            g.equivalents(k, nb);
            for &e in nb.iter() {
                if tcas(tab, e, u, v) {
                    *cnt += 1;
                    out.push(e);
                }
            }
        }
    };
    for f in frontier {
        g.neighbours(f, &mut cache, |j| {
            prefetch_write(&tab[j / PER5]);
            let old = std::mem::replace(&mut ring[head], j);
            head = (head + 1) % RING;
            if old != usize::MAX {
                write(old, &mut cnt, out, &mut nb);
            }
        });
    }
    for k in ring {
        if k != usize::MAX {
            write(k, &mut cnt, out, &mut nb);
        }
    }
    cnt
}

/// One entry of a backward check in progress (as in `crate::prune`).
#[derive(Clone, Copy)]
struct Lane<N> {
    i: usize,
    j: usize,
    m: usize,
    x: N,
}

/// Backward check of `entries` (holding `u`): an entry with a neighbour holding `f` is a hit
/// and keeps `u`; the others move to `next` and are appended to `misses`. Returns the hits
/// and the misses.
fn check<G: Graph>(
    g: &G,
    tab: &[AtomicU64],
    (f, u, next): (u8, u8, u8),
    entries: &mut dyn Iterator<Item = usize>,
    misses: &mut Appender,
) -> (u64, u64) {
    let mut cache = G::Cache::default();
    let mut start = |i: usize| -> Lane<G::Node> {
        let x = g.node(i, &mut cache);
        let j = g.step(&x, 0);
        prefetch(&tab[j / PER5]);
        Lane { i, j, m: 0, x }
    };
    let mut lanes: Vec<Lane<G::Node>> = Vec::with_capacity(LANES);
    while lanes.len() < LANES {
        match entries.next() {
            Some(i) => lanes.push(start(i)),
            None => break,
        }
    }
    let (mut hits, mut miss) = (0, 0);
    let mut k = 0;
    while !lanes.is_empty() {
        let l = &mut lanes[k];
        let hit = tget(tab, l.j) == f;
        if hit || l.m + 1 == G::MOVES {
            if hit {
                hits += 1;
            } else {
                miss += 1;
                tmove(tab, l.i, u, next);
                misses.push(l.i);
            }
            if let Some(i) = entries.next() {
                *l = start(i);
                k += 1;
            } else {
                lanes.swap_remove(k);
            }
        } else {
            l.m += 1;
            l.j = g.step(&l.x, l.m);
            prefetch(&tab[l.j / PER5]);
            k += 1;
        }
        if k >= lanes.len() {
            k = 0;
        }
    }
    (hits, miss)
}

/// Sum over the work items' counts.
fn add2(a: (u64, u64), b: (u64, u64)) -> (u64, u64) {
    (a.0 + b.0, a.1 + b.1)
}

/// How the next level is made.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Step {
    Forward,
    Backward,
    Listed,
}

/// Breadth-first search from `targets` over the `n` entries of `g`, straight into a packed
/// table of the distances mod 3 (every entry must be reachable). `None` if `cfg.stop` was set.
///
/// # Errors
///
/// If a list file cannot be written or read, or the levels do not add up (a graph whose
/// moves are not symmetric, or with unreachable entries).
pub fn bfs5<G: Graph>(
    g: &G,
    n: usize,
    targets: &[usize],
    cfg: &Bfs5Config,
) -> std::io::Result<Option<Table5>> {
    let t0 = std::time::Instant::now();
    let bad = |m: String| std::io::Error::other(format!("{}: {m}", cfg.name));
    let words = n.div_ceil(PER5);
    let mut data = HugeBuf::<u64>::for_overwrite(words);
    let unseen = fill_word(2);
    data.par_chunks_mut(WCH).for_each(|c| c.fill(unseen));
    let tab = data.as_atomic();
    let mut tg: Vec<usize> = targets.iter().copied().filter(|&i| i < n).collect();
    tg.sort_unstable();
    tg.dedup();
    for &t in &tg {
        tcas(tab, t, 2, 0);
    }
    let mut hist = vec![tg.len() as u64];
    let mut lists: Vec<Option<List>> =
        vec![Some(List::Mem(tg.iter().map(|&t| t as u64).collect()))];
    let mut unvisited = (n - tg.len()) as u64;
    let mut step = Step::Forward;
    let mut rest: Option<List> = None;
    let mut l = 1usize;
    let mv = G::MOVES as u64;
    let name = |what: &str, l: usize| format!("{}-{what}{l}.list", cfg.name);
    while unvisited > 0 {
        if cfg.stop.load(Ordering::Relaxed) {
            for x in lists.into_iter().flatten() {
                x.remove();
            }
            if let Some(r) = rest {
                r.remove();
            }
            return Ok(None);
        }
        let (u, f) = (l % 3, (l + 2) % 3);
        let (u, f) = (u as u8, f as u8);
        let ts = std::time::Instant::now();
        let kind = step;
        let new = match step {
            Step::Forward => {
                // (the forward convention: unvisited entries hold (l + 1) mod 3, and so do the
                // levels in `moved`, which go to the frontier's value first)
                let uf = ((l + 1) % 3) as u8;
                let moved: Vec<usize> = (0..l.saturating_sub(1))
                    .filter(|&k| k % 3 == (l + 1) % 3)
                    .collect();
                for &k in &moved {
                    if let Some(x) = &lists[k] {
                        x.for_each(|i| {
                            tcas(tab, i, uf, f);
                        })?;
                    }
                }
                // record this level if it is likely small enough to be listed
                let record = hist[l - 1].saturating_mul(mv) <= cfg.list_cap as u64;
                let w = if record {
                    Some(ListWriter::new(
                        cfg.scratch,
                        &name("level", l),
                        cfg.list_cap,
                    )?)
                } else {
                    None
                };
                let cnt = AtomicU64::new(0);
                let run = |it: &mut dyn Iterator<Item = usize>| {
                    let mut out = appender(w.as_ref());
                    cnt.fetch_add(expand(g, tab, (uf, u), it, &mut out), Ordering::Relaxed);
                };
                // (the frontier from its list if there is one, else from a scan)
                match &lists[l - 1] {
                    Some(x) => x.for_runs(|it| run(it))?,
                    None => scan(tab, n, f, |it| run(it)),
                }
                let new = cnt.load(Ordering::Relaxed);
                let left = unvisited - new;
                let level = match w {
                    Some(w) => w.finish()?,
                    None => None,
                };
                lists.push(level);
                // the next level forward if its frontier is small and the levels it moves are
                // listed; else backward, whose convention (unvisited hold (l + 1) mod 3) the
                // entries already follow
                let next_moved = (0..l).filter(|&k| k % 3 == (l + 2) % 3);
                let can = next_moved.clone().all(|k| lists[k].is_some());
                let fwd = can && new.saturating_mul(mv) < left.saturating_mul(cfg.forward_ratio);
                if fwd && left > 0 {
                    rotate(tab, uf, ((l + 2) % 3) as u8);
                } else {
                    step = Step::Backward;
                }
                for &k in &moved {
                    if let Some(x) = &lists[k] {
                        x.for_each(|i| {
                            tcas(tab, i, f, uf);
                        })?;
                    }
                }
                new
            }
            Step::Backward => {
                let w = ListWriter::new(cfg.scratch, &name("left", l), cfg.list_cap)?;
                let next = ((l + 1) % 3) as u8;
                let tgs = &tg;
                let (hits, miss) = {
                    let acc = Mutex::new((0u64, 0u64));
                    scan(tab, n, u, |it| {
                        let mut it = it.filter(|i| tgs.binary_search(i).is_err());
                        let mut out = appender(Some(&w));
                        let r = check(g, tab, (f, u, next), &mut it, &mut out);
                        let mut a = acc
                            .lock()
                            .unwrap_or_else(std::sync::PoisonError::into_inner);
                        *a = add2(*a, r);
                    });
                    acc.into_inner()
                        .unwrap_or_else(std::sync::PoisonError::into_inner)
                };
                // the hits include the earlier levels holding the same value
                let shared: u64 = (1..l).filter(|&k| k % 3 == l % 3).map(|k| hist[k]).sum();
                let new = hits
                    .checked_sub(shared)
                    .ok_or_else(|| bad(format!("level {l}: {hits} hits, {shared} shared")))?;
                if new > unvisited || miss != unvisited - new {
                    return Err(bad(format!(
                        "level {l}: {new} new and {miss} left of {unvisited} unvisited"
                    )));
                }
                lists.push(None);
                rest = w.finish()?;
                if rest.is_some() {
                    step = Step::Listed;
                }
                new
            }
            Step::Listed => {
                let Some(r) = rest.take() else {
                    return Err(bad(format!("level {l}: no list of the unvisited entries")));
                };
                let w = ListWriter::new(cfg.scratch, &name("left", l), cfg.list_cap)?;
                let next = ((l + 1) % 3) as u8;
                let acc = Mutex::new((0u64, 0u64));
                r.for_runs(|it| {
                    let mut out = appender(Some(&w));
                    let res = check(g, tab, (f, u, next), it, &mut out);
                    let mut a = acc
                        .lock()
                        .unwrap_or_else(std::sync::PoisonError::into_inner);
                    *a = add2(*a, res);
                })?;
                r.remove();
                let (hits, miss) = acc
                    .into_inner()
                    .unwrap_or_else(std::sync::PoisonError::into_inner);
                if hits + miss != unvisited {
                    return Err(bad(format!(
                        "level {l}: {hits} + {miss} checked of {unvisited} unvisited"
                    )));
                }
                lists.push(None);
                rest = w.finish()?;
                if rest.is_none() {
                    step = Step::Backward;
                }
                hits
            }
        };
        if new == 0 {
            return Err(bad(format!(
                "level {l} is empty with {unvisited} entries unvisited"
            )));
        }
        unvisited -= new;
        hist.push(new);
        if cfg.verbose {
            eprintln!(
                "  {}: depth {l} {new} ({}) {:.1}s",
                cfg.name,
                match kind {
                    Step::Forward => "fwd",
                    Step::Backward => "bwd",
                    Step::Listed => "listed",
                },
                ts.elapsed().as_secs_f64()
            );
        }
        l += 1;
    }
    for x in lists.into_iter().flatten() {
        x.remove();
    }
    if let Some(r) = rest {
        r.remove();
    }
    // (the digits past the last entry: 0, as `Table5::pack` writes them)
    if let Some(last) = data.last_mut() {
        let used = n - (words - 1) * PER5;
        for i in used..PER5 {
            let q = i / 5;
            let d = TRITS[usize::from((*last >> ((q & 7) << 3)) as u8)][i - q * 5];
            *last = last.wrapping_sub(u64::from(d) * (POW3[i - q * 5] << ((q & 7) << 3)));
        }
    }
    if cfg.verbose {
        let mean = hist
            .iter()
            .enumerate()
            .map(|(d, &c)| d as f64 * c as f64)
            .sum::<f64>()
            / n as f64;
        eprintln!(
            "{}: {n} entries, max depth {}, mean {mean:.3}, {:.1}s",
            cfg.name,
            hist.len() - 1,
            t0.elapsed().as_secs_f64()
        );
    }
    Ok(Some(Table5 { data, n, hist }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prune::bfs;

    /// The integers mod an odd `m`, moved by +1, -1, x2 and /2 (diameter about twice log2(m)),
    /// with every state indexed twice (`i` and `i + m` are equivalent; the halving move leads to
    /// the second index): equivalents, levels of several sizes and symmetric moves.
    struct Twin(usize);

    impl Graph for Twin {
        const MOVES: usize = 4;
        type Node = usize;
        type Cache = ();
        fn node(&self, i: usize, (): &mut ()) -> usize {
            i % self.0
        }
        fn step(&self, &x: &usize, k: usize) -> usize {
            let m = self.0;
            match k {
                0 => (x + 1) % m,
                1 => (x + m - 1) % m,
                2 => 2 * x % m,
                // (halving: the inverse of 2 mod an odd m is (m + 1) / 2)
                _ => x * m.div_ceil(2) % m + m,
            }
        }
        fn equivalents(&self, j: usize, out: &mut Vec<usize>) {
            out.push(if j < self.0 { j + self.0 } else { j - self.0 });
        }
    }

    /// The packed search equals the 2-bit one entry for entry, with every kind of step:
    /// forward as long as possible, backward from the second level, no lists (backward scans
    /// only), and the default; lists in memory and in files.
    #[test]
    fn packed_search_equals_two_bit_search() {
        let dir = std::env::temp_dir().join(format!("omelettehedron-bfs5-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let stop = AtomicBool::new(false);
        for m in [7, 41, 1001, 30_011, 200_003, 2_000_001] {
            let g = Twin(m);
            let n = 2 * m;
            let want = bfs(&g, n, &[0, m], "twin", false);
            for (ratio, cap, files) in [
                (u64::MAX / 64, usize::MAX, false),
                (0, usize::MAX, true),
                (FORWARD_RATIO, 0, false),
                (FORWARD_RATIO, 1 << 32, true),
                (FORWARD_RATIO, 40, false),
            ] {
                let cfg = Bfs5Config {
                    name: "twin",
                    verbose: false,
                    scratch: files.then_some(dir.as_path()),
                    list_cap: cap,
                    forward_ratio: ratio,
                    stop: &stop,
                };
                let t = bfs5(&g, n, &[0, m], &cfg).unwrap().unwrap();
                assert_eq!(t.hist, want.hist, "m {m}, ratio {ratio}, cap {cap}");
                assert!(
                    (0..n).all(|i| t.get(i) == want.get(i)),
                    "m {m}, ratio {ratio}, cap {cap}"
                );
                // (the digits past the last entry are 0, as `Table5::pack` writes them)
                let p = Table5::pack(bfs(&g, n, &[0, m], "twin", false));
                assert_eq!(&t.data[..], &p.data[..]);
            }
        }
        // no list file is left behind
        assert_eq!(std::fs::read_dir(&dir).unwrap().count(), 0);
        std::fs::remove_dir_all(&dir).unwrap();
    }

    /// A stop flag set before the search gives up at once.
    #[test]
    fn packed_search_stops() {
        let stop = AtomicBool::new(true);
        let cfg = Bfs5Config::new("twin", false, &stop);
        assert!(bfs5(&Twin(100), 200, &[0, 100], &cfg).unwrap().is_none());
    }

    /// The byte tables: `MATCH` marks the digits of each value, and a move changes one digit.
    #[test]
    fn digit_tables() {
        let tab: Vec<AtomicU64> = (0..3).map(|_| AtomicU64::new(fill_word(2))).collect();
        for i in 0..120 {
            assert_eq!(tget(&tab, i), 2);
        }
        for i in (0..120).step_by(7) {
            assert!(tcas(&tab, i, 2, (i % 3) as u8));
            if i % 3 != 2 {
                // (no longer holding 2)
                assert!(!tcas(&tab, i, 2, 1));
            }
        }
        for i in 0..120 {
            let want = if i % 7 == 0 { (i % 3) as u8 } else { 2 };
            assert_eq!(tget(&tab, i), want, "entry {i}");
        }
        rotate(&tab, 2, 0);
        for i in 0..120 {
            let want = if i % 7 == 0 && i % 3 != 2 {
                (i % 3) as u8
            } else {
                0
            };
            assert_eq!(tget(&tab, i), want, "entry {i}");
        }
        for b in 0..243usize {
            for v in 0..3 {
                let mut x = b;
                for j in 0..5 {
                    assert_eq!(MATCH[v][b] >> j & 1 == 1, x % 3 == v);
                    x /= 3;
                }
            }
        }
    }
}
