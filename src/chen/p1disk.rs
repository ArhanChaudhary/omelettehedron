//! Breadth-first search of the oriented phase-1 space ([`super::oriented`]) on disk, a level at a
//! time, in bit planes.
//!
//! The space has 3.34e11 slots (1,396 rows of 239,500,800 bits, 29.9 MB a row), so one bit
//! plane is 41.8 GB: more than the memory this search may use. Each level is pulled: a slot is
//! new at level `d` if it is not yet visited and one of its 16 neighbours is at level `d - 1`.
//! The neighbours of a row's slots lie in at most 16 rows (one per move), so a destination row
//! needs its 16 source rows of the previous level's plane; those are read from disk into a
//! cache of whole rows (least recently used out), and the destination rows are taken in an
//! order that keeps the sources of consecutive rows close ([`row_order`]).
//!
//! Files in the search directory: `V.bits` (visited), `L<d>.bits` (level `d`), `L<d>.rows` (its
//! slots per row, written when the level is complete) and `levels.tsv`. A level is redone from
//! the start if it was interrupted: a row's update is idempotent (the slots visited before the
//! level are `V & !L<d>`).

use super::oriented::{
    P1Space, Target, GROUP_BYTES, MASK81, N_EPERM, N_TURN, ROW_BYTES, ROW_STRIDE,
};
use crate::coord::N_MSET;
use crate::cube::N_MOVES;
use crate::mem::HugeBuf;
use rayon::prelude::*;
use std::collections::{HashMap, VecDeque};
use std::fs::{File, OpenOptions};
use std::io::Write;
use std::os::unix::fs::FileExt;
#[cfg(target_os = "linux")]
use std::os::unix::fs::OpenOptionsExt;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::time::Instant;

/// A failure of the search: a file, or a check of its data.
#[derive(Debug, thiserror::Error)]
pub enum P1Error {
    #[error("{0}: {1}")]
    Io(String, std::io::Error),
    #[error("{0}")]
    Bad(String),
}

type Res<T> = Result<T, P1Error>;

fn io<T>(what: impl Fn() -> String, r: std::io::Result<T>) -> Res<T> {
    r.map_err(|e| P1Error::Io(what(), e))
}

/// Direct I/O (no page cache) where the system has it (Linux).
fn direct(o: &mut OpenOptions) -> &mut OpenOptions {
    #[cfg(target_os = "linux")]
    o.custom_flags(libc::O_DIRECT);
    o
}

/// A row buffer.
pub type Row = HugeBuf<u8>;

/// A fresh row buffer (contents unspecified).
pub fn row_buf() -> Row {
    HugeBuf::for_overwrite(ROW_STRIDE)
}

/// A file of rows of [`ROW_STRIDE`] bytes, read and written with direct I/O.
pub struct Plane {
    file: File,
    path: PathBuf,
}

impl Plane {
    /// Opens (or creates, sized for `rows`) a plane.
    ///
    /// # Errors
    ///
    /// If the file cannot be opened or sized.
    pub fn open(path: &Path, rows: usize, create: bool) -> Res<Plane> {
        let file = io(
            || format!("open {}", path.display()),
            direct(OpenOptions::new().read(true).write(true).create(create)).open(path),
        )?;
        let len = (rows * ROW_STRIDE) as u64;
        let cur = io(|| format!("stat {}", path.display()), file.metadata())?.len();
        if cur != len {
            if !create {
                return Err(P1Error::Bad(format!(
                    "{}: {cur} bytes, expected {len}",
                    path.display()
                )));
            }
            io(|| format!("size {}", path.display()), file.set_len(len))?;
        }
        Ok(Plane {
            file,
            path: path.to_path_buf(),
        })
    }

    /// Reads a row (the whole stride).
    ///
    /// # Errors
    ///
    /// If the read fails.
    pub fn read(&self, row: usize, buf: &mut [u8]) -> Res<()> {
        io(
            || format!("read row {row} of {}", self.path.display()),
            self.file
                .read_exact_at(&mut buf[..ROW_STRIDE], (row * ROW_STRIDE) as u64),
        )
    }

    /// Writes a row (the whole stride).
    ///
    /// # Errors
    ///
    /// If the write fails.
    pub fn write(&self, row: usize, buf: &[u8]) -> Res<()> {
        io(
            || format!("write row {row} of {}", self.path.display()),
            self.file
                .write_all_at(&buf[..ROW_STRIDE], (row * ROW_STRIDE) as u64),
        )
    }

    /// Flushes the file to the disk.
    ///
    /// # Errors
    ///
    /// If the sync fails.
    pub fn sync(&self) -> Res<()> {
        io(
            || format!("sync {}", self.path.display()),
            self.file.sync_all(),
        )
    }
}

/// The 81-bit block `k` of a row (the row has 16 bytes to spare at its end).
#[inline(always)]
pub fn read81(row: &[u8], k: usize) -> u128 {
    let bit = k * N_TURN;
    let at = bit >> 3;
    let mut b = [0u8; 16];
    b.copy_from_slice(&row[at..at + 16]);
    (u128::from_le_bytes(b) >> (bit & 7)) & MASK81
}

/// The 8 blocks of a group of 81 bytes.
#[inline(always)]
fn unpack(g: &[u8]) -> [u128; 8] {
    let mut buf = [0u8; GROUP_BYTES + 16];
    buf[..GROUP_BYTES].copy_from_slice(&g[..GROUP_BYTES]);
    std::array::from_fn(|b| read81(&buf, b))
}

/// 8 blocks into a group of 81 bytes.
#[inline(always)]
fn pack(blocks: &[u128; 8], g: &mut [u8]) {
    let mut buf = [0u8; GROUP_BYTES + 16];
    for (b, &x) in blocks.iter().enumerate() {
        let bit = b * N_TURN;
        let at = bit >> 3;
        let mut w = [0u8; 16];
        w.copy_from_slice(&buf[at..at + 16]);
        let v = u128::from_le_bytes(w) | (x & MASK81) << (bit & 7);
        buf[at..at + 16].copy_from_slice(&v.to_le_bytes());
    }
    g[..GROUP_BYTES].copy_from_slice(&buf[..GROUP_BYTES]);
}

/// Groups per parallel chunk of a row.
const CHUNK: usize = 1024;

/// Pulls one destination row: `v` holds the visited bits and `l` the level's bits so far
/// (zero unless the level is being redone); both are updated. Returns the new slots.
pub fn pull_row(
    sp: &P1Space,
    row: usize,
    srcs: &[Option<&[u8]>; N_MOVES],
    v: &mut [u8],
    l: &mut [u8],
) -> u64 {
    let nb = &sp.nbr[row];
    v[..ROW_BYTES]
        .par_chunks_mut(GROUP_BYTES * CHUNK)
        .zip(l[..ROW_BYTES].par_chunks_mut(GROUP_BYTES * CHUNK))
        .enumerate()
        .map(|(ci, (vc, lc))| {
            let mut cnt = 0u64;
            for (j, (vg, lg)) in vc
                .chunks_mut(GROUP_BYTES)
                .zip(lc.chunks_mut(GROUP_BYTES))
                .enumerate()
            {
                let e = (ci * CHUNK + j) as u32;
                let vb = unpack(vg);
                let lb = unpack(lg);
                let vold: [u128; 8] = std::array::from_fn(|b| vb[b] & !lb[b]);
                if vold.iter().all(|&x| x == MASK81) {
                    if lb.iter().any(|&x| x != 0) {
                        lg[..GROUP_BYTES].fill(0);
                    }
                    continue;
                }
                let mut acc = [0u128; 8];
                for m in 0..N_MOVES {
                    let Some(src) = srcs[m] else { continue };
                    let n = nb[m];
                    let (g, t) = (n.g as usize, n.t as usize);
                    let se = sp.step_e(m, g, t, e);
                    for b3 in 0..8u8 {
                        if (acc[b3 as usize] | vold[b3 as usize]) == MASK81 {
                            continue;
                        }
                        let (bits3, tr) = sp.step_b(&se, g, sp.bits4(e, b3));
                        let k3 = (se.e3 * 8 + u32::from(bits3 & 7)) as usize;
                        let blk = read81(src, k3);
                        if blk != 0 {
                            acc[b3 as usize] |= sp.permute(g, tr, blk);
                        }
                    }
                    if (0..8).all(|b| (acc[b] | vold[b]) == MASK81) {
                        break;
                    }
                }
                let new: [u128; 8] = std::array::from_fn(|b| acc[b] & !vold[b]);
                let vis: [u128; 8] = std::array::from_fn(|b| vold[b] | new[b]);
                cnt += new.iter().map(|x| u64::from(x.count_ones())).sum::<u64>();
                pack(&new, lg);
                pack(&vis, vg);
            }
            cnt
        })
        .sum()
}

/// Destination rows in breadth-first order of the row graph from the solved row (consecutive
/// rows share many sources).
pub fn row_order(sp: &P1Space) -> Vec<u32> {
    let n = sp.rows();
    let mut seen = vec![false; n];
    let mut out = Vec::with_capacity(n);
    for start in 0..n {
        if seen[start] {
            continue;
        }
        seen[start] = true;
        let mut q = VecDeque::from([start]);
        while let Some(r) = q.pop_front() {
            out.push(r as u32);
            for m in 0..N_MOVES {
                let s = sp.nbr[r][m].row as usize;
                if !seen[s] {
                    seen[s] = true;
                    q.push_back(s);
                }
            }
        }
    }
    out
}

/// Destination rows chosen greedily: next the row with the most sources in a cache of `cap`
/// rows that holds the most recently used sources (ties: the lowest row).
pub fn greedy_order(sp: &P1Space, cap: usize) -> Vec<u32> {
    let n = sp.rows();
    let srcs: Vec<Vec<u32>> = (0..n)
        .map(|r| {
            let mut v: Vec<u32> = (0..N_MOVES).map(|m| sp.nbr[r][m].row).collect();
            v.sort_unstable();
            v.dedup();
            v
        })
        .collect();
    let mut last = vec![0usize; n];
    let mut cached = vec![false; n];
    let mut ncached = 0usize;
    let mut done = vec![false; n];
    let mut out = Vec::with_capacity(n);
    for step in 1..=n {
        let best = (0..n)
            .filter(|&r| !done[r])
            .max_by_key(|&r| {
                let hits = srcs[r].iter().filter(|&&s| cached[s as usize]).count();
                (hits * 64 / srcs[r].len().max(1), std::cmp::Reverse(r))
            })
            .unwrap_or(0);
        done[best] = true;
        out.push(best as u32);
        for &s in &srcs[best] {
            let s = s as usize;
            last[s] = step;
            if !cached[s] {
                cached[s] = true;
                ncached += 1;
            }
        }
        while ncached > cap {
            let victim = (0..n)
                .filter(|&r| cached[r] && last[r] < step)
                .min_by_key(|&r| last[r]);
            match victim {
                Some(v) => {
                    cached[v] = false;
                    ncached -= 1;
                }
                None => break,
            }
        }
    }
    out
}

/// Source rows read for one pass over `order` with a least-recently-used cache of `cap` rows.
pub fn simulate_reads(sp: &P1Space, order: &[u32], cap: usize) -> usize {
    let mut cache: HashMap<u32, usize> = HashMap::new();
    let mut reads = 0;
    for (i, &r) in order.iter().enumerate() {
        for m in 0..N_MOVES {
            let s = sp.nbr[r as usize][m].row;
            if cache.insert(s, i).is_none() {
                reads += 1;
                if cache.len() > cap {
                    if let Some((&old, _)) =
                        cache.iter().filter(|(_, &u)| u < i).min_by_key(|(_, &u)| u)
                    {
                        cache.remove(&old);
                    }
                }
            }
        }
    }
    reads
}

/// `MemAvailable`, bytes.
pub fn mem_available() -> Option<u64> {
    let s = std::fs::read_to_string("/proc/meminfo").ok()?;
    let line = s.lines().find(|l| l.starts_with("MemAvailable:"))?;
    let kb: u64 = line.split_whitespace().nth(1)?.parse().ok()?;
    Some(kb * 1024)
}

/// Search settings.
#[derive(Clone, Debug)]
pub struct Config {
    pub dir: PathBuf,
    pub target: Target,
    /// memory for the cache of source rows, bytes
    pub cache_bytes: u64,
    /// keep at least this much memory available (the cache shrinks below it), bytes
    pub floor_bytes: u64,
    /// stop after this level
    pub max_level: Option<usize>,
}

/// A level's totals.
#[derive(Clone, Copy, Debug, Default)]
pub struct Level {
    pub d: usize,
    pub slots: u64,
    /// slots weighted by class size: A-centre/edge states
    pub states: u128,
}

fn plane_path(dir: &Path, d: usize) -> PathBuf {
    dir.join(format!("L{d:02}.bits"))
}

fn rows_path(dir: &Path, d: usize) -> PathBuf {
    dir.join(format!("L{d:02}.rows"))
}

/// Per-row slot counts of a complete level, if it is complete.
///
/// # Errors
///
/// If the file exists but cannot be read or has the wrong size.
pub fn level_rows(dir: &Path, d: usize, rows: usize) -> Res<Option<Vec<u64>>> {
    let p = rows_path(dir, d);
    if !p.exists() {
        return Ok(None);
    }
    let b = io(|| format!("read {}", p.display()), std::fs::read(&p))?;
    if b.len() != rows * 8 {
        return Err(P1Error::Bad(format!("{}: {} bytes", p.display(), b.len())));
    }
    Ok(Some(
        b.as_chunks::<8>()
            .0
            .iter()
            .map(|c| u64::from_le_bytes(*c))
            .collect(),
    ))
}

fn write_rows(dir: &Path, d: usize, counts: &[u64]) -> Res<()> {
    let p = rows_path(dir, d);
    let tmp = p.with_extension("rows.tmp");
    let b: Vec<u8> = counts.iter().flat_map(|c| c.to_le_bytes()).collect();
    io(
        || format!("write {}", tmp.display()),
        std::fs::write(&tmp, b),
    )?;
    io(
        || format!("rename {}", tmp.display()),
        std::fs::rename(&tmp, &p),
    )
}

/// The levels made so far (complete ones).
///
/// # Errors
///
/// If a level's row counts cannot be read.
pub fn levels(sp: &P1Space, dir: &Path) -> Res<Vec<Level>> {
    let mut out = vec![];
    for d in 0.. {
        let Some(c) = level_rows(dir, d, sp.rows())? else {
            break;
        };
        out.push(Level {
            d,
            slots: c.iter().sum(),
            states: c
                .iter()
                .zip(&sp.weight)
                .map(|(&x, &w)| u128::from(x) * u128::from(w))
                .sum(),
        });
    }
    Ok(out)
}

fn log(dir: &Path, line: &str) {
    let stamp = chrono_like();
    let text = format!("{stamp} {line}");
    println!("{text}");
    if let Ok(mut f) = OpenOptions::new()
        .create(true)
        .append(true)
        .open(dir.join("run.log"))
    {
        let _ = writeln!(f, "{text}");
    }
}

/// UTC time of day as `HH:MM:SS` (no date crate).
pub fn chrono_like() -> String {
    let s = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |d| d.as_secs());
    format!("{:02}:{:02}:{:02}", s / 3600 % 24, s / 60 % 60, s % 60)
}

/// The cache of source rows (one level's plane).
struct Cache {
    cap: usize,
    map: HashMap<u32, (Arc<Row>, u64)>,
    tick: u64,
    reads: u64,
}

impl Cache {
    fn get(&mut self, r: u32) -> Option<Arc<Row>> {
        self.tick += 1;
        let t = self.tick;
        self.map.get_mut(&r).map(|e| {
            e.1 = t;
            e.0.clone()
        })
    }
    /// Evict least recently used rows not in `keep` down to `cap` rows.
    fn evict(&mut self, keep: &[u32]) {
        while self.map.len() > self.cap {
            let victim = self
                .map
                .iter()
                .filter(|(r, _)| !keep.contains(r))
                .min_by_key(|(_, e)| e.1)
                .map(|(&r, _)| r);
            match victim {
                Some(r) => {
                    self.map.remove(&r);
                }
                None => break,
            }
        }
    }
}

/// The sources (nonempty source rows) of a destination row.
fn sources(sp: &P1Space, row: usize, prev: &[u64]) -> Vec<u32> {
    let mut v: Vec<u32> = (0..N_MOVES)
        .map(|m| sp.nbr[row][m].row)
        .filter(|&s| prev[s as usize] > 0)
        .collect();
    v.sort_unstable();
    v.dedup();
    v
}

/// Loads the given rows of `plane` into the cache (those missing), in parallel.
fn load(cache: &Mutex<Cache>, plane: &Plane, rows: &[u32], keep: &[u32]) -> Res<()> {
    let missing: Vec<u32> = {
        let Ok(mut c) = cache.lock() else {
            return Err(P1Error::Bad("cache lock".into()));
        };
        rows.iter()
            .copied()
            .filter(|&r| c.get(r).is_none())
            .collect()
    };
    let loaded: Vec<Res<(u32, Row)>> = missing
        .par_iter()
        .map(|&r| {
            let mut b = row_buf();
            plane.read(r as usize, &mut b)?;
            Ok((r, b))
        })
        .collect();
    let Ok(mut c) = cache.lock() else {
        return Err(P1Error::Bad("cache lock".into()));
    };
    for x in loaded {
        let (r, b) = x?;
        c.tick += 1;
        let t = c.tick;
        c.map.insert(r, (Arc::new(b), t));
        c.reads += 1;
    }
    // shrink the cache when memory runs low
    if let Some(avail) = mem_available() {
        let _ = avail;
    }
    c.evict(keep);
    Ok(())
}

/// A destination row's buffers: visited, level.
struct Dest {
    row: usize,
    v: Row,
    l: Row,
    count: u64,
}

/// Runs (or resumes) the search until the frontier is empty or `max_level`.
///
/// # Errors
///
/// On a file error or a failed check.
pub fn run(sp: &P1Space, cfg: &Config) -> Res<Vec<Level>> {
    let dir = &cfg.dir;
    io(
        || format!("create {}", dir.display()),
        std::fs::create_dir_all(dir),
    )?;
    let nrows = sp.rows();
    let manifest = dir.join("manifest.txt");
    let want = format!(
        "target {}\nrows {nrows}\nrow_stride {ROW_STRIDE}\nslots {}\n",
        cfg.target.name(),
        sp.slots()
    );
    if manifest.exists() {
        let have = io(
            || "read manifest".into(),
            std::fs::read_to_string(&manifest),
        )?;
        if have != want {
            return Err(P1Error::Bad(format!(
                "{} holds another search:\n{have}",
                manifest.display()
            )));
        }
    } else {
        io(|| "write manifest".into(), std::fs::write(&manifest, &want))?;
    }
    let vplane = Plane::open(&dir.join("V.bits"), nrows, true)?;
    // level 0
    if level_rows(dir, 0, nrows)?.is_none() {
        let l0 = Plane::open(&plane_path(dir, 0), nrows, true)?;
        let ts = sp.target_slots(cfg.target);
        let mut counts = vec![0u64; nrows];
        let mut by_row: HashMap<usize, Vec<u32>> = HashMap::new();
        for &(r, c) in &ts {
            by_row.entry(r).or_default().push(c);
        }
        for (&r, cols) in &by_row {
            let mut b = HugeBuf::<u8>::zeroed(ROW_STRIDE);
            for &c in cols {
                b[(c >> 3) as usize] |= 1 << (c & 7);
            }
            counts[r] = cols.len() as u64;
            l0.write(r, &b)?;
            vplane.write(r, &b)?;
        }
        l0.sync()?;
        vplane.sync()?;
        write_rows(dir, 0, &counts)?;
        log(dir, &format!("level 0: {} target slots", ts.len()));
    }
    let row_bytes = ROW_STRIDE as u64;
    let order = greedy_order(
        sp,
        ((cfg.cache_bytes / row_bytes) as usize).max(2 * N_MOVES + 2),
    );
    for d in 1.. {
        if cfg.max_level.is_some_and(|m| d > m) {
            break;
        }
        if level_rows(dir, d, nrows)?.is_some() {
            continue;
        }
        let Some(prev) = level_rows(dir, d - 1, nrows)? else {
            return Err(P1Error::Bad(format!("level {} missing", d - 1)));
        };
        if prev.iter().all(|&c| c == 0) {
            break;
        }
        let t0 = Instant::now();
        let ppath = plane_path(dir, d);
        let resume = ppath.exists();
        let lprev = Plane::open(&plane_path(dir, d - 1), nrows, false)?;
        let lcur = Plane::open(&ppath, nrows, true)?;
        let cap = ((cfg.cache_bytes / row_bytes) as usize).max(2 * N_MOVES + 2);
        let cache = Mutex::new(Cache {
            cap,
            map: HashMap::new(),
            tick: 0,
            reads: 0,
        });
        let mut counts = vec![0u64; nrows];
        // the rows that get anything: those with a nonempty source
        let work: Vec<u32> = order
            .iter()
            .copied()
            .filter(|&r| !sources(sp, r as usize, &prev).is_empty())
            .collect();
        log(
            dir,
            &format!(
                "level {d}: {} destination rows, cache {cap} rows{}",
                work.len(),
                if resume { " (redoing the level)" } else { "" }
            ),
        );
        // pipeline: compute row i while writing row i-1 and reading row i+1
        let read_dest = |row: usize| -> Res<Dest> {
            let mut v = row_buf();
            vplane.read(row, &mut v)?;
            let mut l = row_buf();
            if resume {
                lcur.read(row, &mut l)?;
            } else {
                l[..ROW_STRIDE].fill(0);
            }
            Ok(Dest {
                row,
                v,
                l,
                count: 0,
            })
        };
        // (a row with nothing new changes neither file: had the level been begun before, its
        // slots would be found again)
        let write_dest = |x: &Dest| -> Res<()> {
            if x.count == 0 {
                return Ok(());
            }
            lcur.write(x.row, &x.l)?;
            vplane.write(x.row, &x.v)
        };
        let mut done: Option<Dest> = None;
        let mut next: Option<Dest> = None;
        let mut last_log = Instant::now();
        let mut low_mem = 0u64;
        for i in 0..work.len() {
            let row = work[i] as usize;
            let srcs_i = sources(sp, row, &prev);
            if next.is_none() {
                load(&cache, &lprev, &srcs_i, &srcs_i)?;
                next = Some(read_dest(row)?);
            }
            let Some(mut cur) = next.take() else {
                return Err(P1Error::Bad("pipeline".into()));
            };
            // the source rows of this row, pinned for the computation
            let pinned: HashMap<u32, Arc<Row>> = {
                let Ok(mut c) = cache.lock() else {
                    return Err(P1Error::Bad("cache lock".into()));
                };
                srcs_i
                    .iter()
                    .filter_map(|&s| c.get(s).map(|b| (s, b)))
                    .collect()
            };
            if pinned.len() != srcs_i.len() {
                return Err(P1Error::Bad(format!("row {row}: sources not loaded")));
            }
            let srcs: [Option<&[u8]>; N_MOVES] = std::array::from_fn(|m| {
                let s = sp.nbr[row][m].row;
                pinned.get(&s).map(|b| &b[..])
            });
            let nxt = work.get(i + 1).map(|&r| r as usize);
            let prev_done = done.take();
            // memory guard: shrink the cache when MemAvailable is below the floor
            if let Some(avail) = mem_available() {
                if avail < cfg.floor_bytes {
                    if let Ok(mut c) = cache.lock() {
                        let shrink = ((cfg.floor_bytes - avail) / row_bytes + 1) as usize;
                        c.cap = c.cap.saturating_sub(shrink).max(2 * N_MOVES + 2);
                        low_mem += 1;
                    }
                }
            }
            let (count, io_res) = rayon::join(
                || pull_row(sp, row, &srcs, &mut cur.v, &mut cur.l),
                || -> Res<Option<Dest>> {
                    if let Some(x) = &prev_done {
                        write_dest(x)?;
                    }
                    match nxt {
                        Some(r) => {
                            let s = sources(sp, r, &prev);
                            let mut keep = s.clone();
                            keep.extend_from_slice(&srcs_i);
                            load(&cache, &lprev, &s, &keep)?;
                            Ok(Some(read_dest(r)?))
                        }
                        None => Ok(None),
                    }
                },
            );
            next = io_res?;
            cur.count = count;
            counts[row] = count;
            done = Some(cur);
            if last_log.elapsed().as_secs() >= 300 {
                last_log = Instant::now();
                let reads = cache.lock().map_or(0, |c| c.reads);
                log(
                    dir,
                    &format!(
                        "level {d}: {}/{} rows, {} new so far, {reads} source rows read, {:.0} s",
                        i + 1,
                        work.len(),
                        counts.iter().sum::<u64>(),
                        t0.elapsed().as_secs_f64()
                    ),
                );
            }
        }
        if let Some(x) = done.take() {
            write_dest(&x)?;
        }
        lcur.sync()?;
        vplane.sync()?;
        write_rows(dir, d, &counts)?;
        let slots: u64 = counts.iter().sum();
        let states: u128 = counts
            .iter()
            .zip(&sp.weight)
            .map(|(&x, &w)| u128::from(x) * u128::from(w))
            .sum();
        let reads = cache.lock().map_or(0, |c| c.reads);
        log(
            dir,
            &format!(
                "level {d} done: {slots} slots, {states} states; {reads} source rows read; \
                 {low_mem} low-memory shrinks; {:.0} s",
                t0.elapsed().as_secs_f64()
            ),
        );
        if slots == 0 {
            break;
        }
    }
    levels(sp, dir)
}

/// Bytes of a row of the 4-bit table: two slots a byte, whole pages.
pub const NIB_STRIDE: usize = (N_EPERM / 2).div_ceil(4096) * 4096;

/// The value of slot `col` in a row of the 4-bit table.
#[inline(always)]
pub fn nib(row: &[u8], col: usize) -> u8 {
    row[col >> 1] >> ((col & 1) * 4) & 15
}

/// Builds the 4-bit table `D.nib` (distance of every slot) from the level planes, and checks
/// that every slot has exactly one level.
///
/// # Errors
///
/// On a file error, a slot in two levels or in none, or a level above 15.
pub fn build_nibbles(sp: &P1Space, dir: &Path) -> Res<Vec<u64>> {
    let nrows = sp.rows();
    let lv = levels(sp, dir)?;
    let top = lv.iter().rposition(|l| l.slots > 0).unwrap_or(0);
    if top > 15 {
        return Err(P1Error::Bad(format!("level {top} does not fit in 4 bits")));
    }
    let planes: Vec<Plane> = (0..=top)
        .map(|d| Plane::open(&plane_path(dir, d), nrows, false))
        .collect::<Res<_>>()?;
    let row_counts: Vec<Vec<u64>> = (0..=top)
        .map(|d| level_rows(dir, d, nrows).map(Option::unwrap_or_default))
        .collect::<Res<_>>()?;
    let path = dir.join("D.nib");
    let out = io(
        || format!("open {}", path.display()),
        direct(OpenOptions::new().read(true).write(true).create(true)).open(&path),
    )?;
    io(
        || "size D.nib".into(),
        out.set_len((nrows * NIB_STRIDE) as u64),
    )?;
    let mut hist = vec![0u64; top + 1];
    let mut nb = HugeBuf::<u8>::for_overwrite(NIB_STRIDE);
    let mut b = row_buf();
    let mut seen = row_buf();
    for row in 0..nrows {
        nb[..NIB_STRIDE].fill(0);
        seen[..ROW_STRIDE].fill(0);
        for d in 0..=top {
            if row_counts[d][row] == 0 {
                continue;
            }
            planes[d].read(row, &mut b)?;
            let dd = d as u8;
            let clash = nb[..N_EPERM / 2]
                .par_chunks_mut(4096)
                .zip(b[..ROW_BYTES].par_chunks(2048))
                .zip(seen[..ROW_BYTES].par_chunks_mut(2048))
                .map(|((nc, bc), sc)| {
                    let mut clash = 0u64;
                    for (i, (&x, s)) in bc.iter().zip(sc.iter_mut()).enumerate() {
                        if x == 0 {
                            continue;
                        }
                        clash += u64::from((x & *s).count_ones());
                        *s |= x;
                        for k in 0..8 {
                            if x >> k & 1 == 1 {
                                let c = 8 * i + k;
                                nc[c >> 1] |= dd << ((c & 1) * 4);
                            }
                        }
                    }
                    clash
                })
                .sum::<u64>();
            if clash > 0 {
                return Err(P1Error::Bad(format!(
                    "row {row}: {clash} slots in two levels"
                )));
            }
            hist[d] += row_counts[d][row];
        }
        let full: u64 = seen[..ROW_BYTES]
            .par_iter()
            .map(|x| u64::from(x.count_ones()))
            .sum();
        if full != N_EPERM as u64 {
            return Err(P1Error::Bad(format!(
                "row {row}: {} slots in no level",
                N_EPERM as u64 - full
            )));
        }
        io(
            || format!("write row {row} of D.nib"),
            out.write_all_at(&nb[..NIB_STRIDE], (row * NIB_STRIDE) as u64),
        )?;
    }
    io(|| "sync D.nib".into(), out.sync_all())?;
    Ok(hist)
}

/// The 4-bit table of a search directory, read a row at a time.
pub struct Nibbles {
    file: File,
}

impl Nibbles {
    /// Opens `D.nib`.
    ///
    /// # Errors
    ///
    /// If it cannot be opened or has the wrong size.
    pub fn open(dir: &Path, rows: usize, use_direct: bool) -> Res<Nibbles> {
        let path = dir.join("D.nib");
        let mut o = OpenOptions::new();
        o.read(true);
        if use_direct {
            direct(&mut o);
        }
        let file = io(|| format!("open {}", path.display()), o.open(&path))?;
        let len = io(|| "stat D.nib".into(), file.metadata())?.len();
        if len != (rows * NIB_STRIDE) as u64 {
            return Err(P1Error::Bad(format!("D.nib: {len} bytes")));
        }
        Ok(Nibbles { file })
    }

    /// Reads a row.
    ///
    /// # Errors
    ///
    /// If the read fails.
    pub fn read(&self, row: usize, buf: &mut [u8]) -> Res<()> {
        io(
            || format!("read row {row} of D.nib"),
            self.file
                .read_exact_at(&mut buf[..NIB_STRIDE], (row * NIB_STRIDE) as u64),
        )
    }

    /// One slot's value (a small read; for walks, not for passes).
    ///
    /// # Errors
    ///
    /// If the read fails.
    pub fn get(&self, row: usize, col: u32) -> Res<u8> {
        let mut b = [0u8; 1];
        io(
            || format!("read slot {row}/{col} of D.nib"),
            self.file
                .read_exact_at(&mut b, (row * NIB_STRIDE + (col >> 1) as usize) as u64),
        )?;
        Ok(b[0] >> ((col & 1) * 4) & 15)
    }
}

/// A part of [`check_table`]'s sums: slots per value, bad targets, bad slots, examples.
type CheckPart = (Vec<u64>, u64, u64, Vec<(usize, u32, u8, u8)>);

/// The result of [`check_table`].
#[derive(Clone, Debug, Default)]
pub struct CheckReport {
    /// slots per value
    pub hist: Vec<u64>,
    /// slots at 0 that are no target, or targets not at 0
    pub bad_targets: u64,
    /// slots above 0 whose smallest neighbour is not one less
    pub bad_slots: u64,
    /// examples of failures: (row, col, value, smallest neighbour)
    pub examples: Vec<(usize, u32, u8, u8)>,
}

/// Local consistency of the 4-bit table over every slot: the slots at 0 are exactly the target
/// slots, and every other slot's smallest neighbour value is its value less one. The neighbour
/// of each slot is found slot by slot ([`P1Space::step`] and the digit map), not through the
/// search's block permutation.
///
/// # Errors
///
/// On a file error.
pub fn check_table(
    sp: &P1Space,
    dir: &Path,
    target: Target,
    cache_rows: usize,
) -> Res<CheckReport> {
    let nrows = sp.rows();
    let tab = Nibbles::open(dir, nrows, true)?;
    let mut targets: HashMap<usize, Vec<u32>> = HashMap::new();
    for (r, c) in sp.target_slots(target) {
        targets.entry(r).or_default().push(c);
    }
    let order = greedy_order(sp, cache_rows);
    let mut cache: HashMap<u32, (Arc<HugeBuf<u8>>, u64)> = HashMap::new();
    let mut tick = 0u64;
    let mut rep = CheckReport {
        hist: vec![0; 16],
        ..CheckReport::default()
    };
    let t0 = Instant::now();
    for (i, &row) in order.iter().enumerate() {
        let row = row as usize;
        let mut need: Vec<u32> = (0..N_MOVES).map(|m| sp.nbr[row][m].row).collect();
        need.push(row as u32);
        need.sort_unstable();
        need.dedup();
        let missing: Vec<u32> = need
            .iter()
            .copied()
            .filter(|r| !cache.contains_key(r))
            .collect();
        let loaded: Vec<Res<(u32, HugeBuf<u8>)>> = missing
            .par_iter()
            .map(|&r| {
                let mut b = HugeBuf::<u8>::for_overwrite(NIB_STRIDE);
                tab.read(r as usize, &mut b)?;
                Ok((r, b))
            })
            .collect();
        for x in loaded {
            let (r, b) = x?;
            cache.insert(r, (Arc::new(b), 0));
        }
        for r in &need {
            tick += 1;
            if let Some(e) = cache.get_mut(r) {
                e.1 = tick;
            }
        }
        while cache.len() > cache_rows.max(need.len()) {
            let victim = cache
                .iter()
                .filter(|(r, _)| !need.contains(r))
                .min_by_key(|(_, e)| e.1)
                .map(|(&r, _)| r);
            match victim {
                Some(r) => {
                    cache.remove(&r);
                }
                None => break,
            }
        }
        let get = |r: u32| cache.get(&r).map(|e| e.0.clone());
        let Some(me) = get(row as u32) else {
            return Err(P1Error::Bad("row not loaded".into()));
        };
        let srcs: Vec<Arc<HugeBuf<u8>>> = (0..N_MOVES)
            .map(|m| get(sp.nbr[row][m].row))
            .collect::<Option<_>>()
            .ok_or_else(|| P1Error::Bad("source not loaded".into()))?;
        let tcols = targets.get(&row).cloned().unwrap_or_default();
        let part: CheckPart = (0..N_MSET)
            .into_par_iter()
            .fold(
                || (vec![0u64; 16], 0u64, 0u64, Vec::new()),
                |mut a, e| {
                    let e = e as u32;
                    let mut mins = [15u8 + 1; 8 * N_TURN];
                    for m in 0..N_MOVES {
                        let n = sp.nbr[row][m];
                        let (g, t) = (n.g as usize, n.t as usize);
                        let src = &srcs[m][..];
                        for b3 in 0..8u8 {
                            let (e3, bits3, tr) = sp.step(m, g, t, e, sp.bits4(e, b3));
                            let base = (e3 * 8 + u32::from(bits3 & 7)) as usize * N_TURN;
                            for r in 0..N_TURN {
                                let r3 = sp.add[sp.lin[g][r] as usize][tr as usize] as usize;
                                let v = nib(src, base + r3);
                                let mm = &mut mins[b3 as usize * N_TURN + r];
                                if v < *mm {
                                    *mm = v;
                                }
                            }
                        }
                    }
                    for b3 in 0..8u32 {
                        for r in 0..N_TURN as u32 {
                            let col = (e * 8 + b3) * N_TURN as u32 + r;
                            let v = nib(&me, col as usize);
                            a.0[v as usize] += 1;
                            let mn = mins[(b3 * N_TURN as u32 + r) as usize];
                            let is_t = tcols.contains(&col);
                            if (v == 0) != is_t {
                                a.1 += 1;
                                if a.3.len() < 10 {
                                    a.3.push((row, col, v, mn));
                                }
                            } else if v > 0 && mn != v - 1 {
                                a.2 += 1;
                                if a.3.len() < 10 {
                                    a.3.push((row, col, v, mn));
                                }
                            }
                        }
                    }
                    a
                },
            )
            .reduce(
                || (vec![0u64; 16], 0u64, 0u64, Vec::new()),
                |mut a, b| {
                    for (x, y) in a.0.iter_mut().zip(&b.0) {
                        *x += y;
                    }
                    a.1 += b.1;
                    a.2 += b.2;
                    a.3.extend(b.3);
                    a.3.truncate(10);
                    a
                },
            );
        for (x, y) in rep.hist.iter_mut().zip(&part.0) {
            *x += y;
        }
        rep.bad_targets += part.1;
        rep.bad_slots += part.2;
        rep.examples.extend(part.3);
        rep.examples.truncate(20);
        if (i + 1) % 100 == 0 {
            log(
                dir,
                &format!(
                    "check: {}/{nrows} rows, {} bad targets, {} bad slots, {:.0} s",
                    i + 1,
                    rep.bad_targets,
                    rep.bad_slots,
                    t0.elapsed().as_secs_f64()
                ),
            );
        }
    }
    Ok(rep)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pack_unpack() {
        let mut g = [0u8; GROUP_BYTES];
        let blocks: [u128; 8] =
            std::array::from_fn(|b| (0x1234_5678_9abc_def0_1357u128 * (b as u128 + 3)) & MASK81);
        pack(&blocks, &mut g);
        assert_eq!(unpack(&g), blocks);
    }
}
