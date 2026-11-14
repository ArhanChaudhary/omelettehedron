//! The search: one layer (all slots at one distance) per step, on disk, a part of the next
//! layer at a time in memory.
//!
//! The move graph is undirected (every allowed move is undone by an allowed move), so the
//! layer at distance `d + 1` is the set of children of layer `d` that are in neither layer `d`
//! nor layer `d - 1`; no set of all visited slots is kept. A step:
//!
//! 1. **Plan**: the classes are split into *parts*, consecutive in [`Instance::order`] (the
//!    classes joined by the always-allowed moves together), whose share of the next layer fits
//!    the memory budget (an upper bound per class from the blocks of layer `d` that move into
//!    it).
//! 2. **Push**, per part: one thread reads the layer-`d` files of the classes with a move into
//!    the part, in order, chunk by chunk; the workers decode each block, and for each allowed
//!    move into the part OR the block's bits, permuted, into the child's block in memory
//!    (blocks allocated on first touch from an arena, 64-bit atomic ORs).
//! 3. **Close** under the stabilisers: a class whose representative is fixed by symmetries
//!    gets the images of its new slots, so every layer is closed under them.
//! 4. **Mask and write**: per class of the part, the new blocks minus layers `d` and `d - 1`
//!    are written as layer `d + 1`, counted (slots, positions, orbits under the 24
//!    symmetries) and the part is recorded in the manifest.
//!
//! A step's parts are recorded one by one in `manifest.json` (written atomically), so a
//! stopped search resumes at the first part not recorded; a layer's files are written under a
//! temporary name and renamed once its part is recorded.

use super::codec::{decode_chunk, Blk, ClassReader, ClassWriter, Summary, CHUNK};
use super::instance::{Instance, Kind, NO_EDGE, N_SYM};
use crate::cube::N_MOVES;
use crate::tables::prefetch_at;
use serde::{Deserialize, Serialize};
use std::fs;
use std::io;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicU32, AtomicU64, AtomicUsize, Ordering};
use std::sync::mpsc::sync_channel;
use std::sync::{Mutex, PoisonError};
use std::time::Instant;

/// The manifest's format version.
const VERSION: u32 = 3;

/// Why a search cannot go on.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum BfsError {
    /// Reading or writing a file failed (a damaged block file included).
    #[error("{path}")]
    Io {
        /// the file
        path: String,
        /// what went wrong
        source: io::Error,
    },
    /// The manifest cannot be read, or belongs to another instance.
    #[error("{path}: {what}")]
    Manifest {
        /// the file
        path: String,
        /// what is wrong
        what: String,
    },
    /// A part's arena was too small (the planning bound is wrong).
    #[error("part {part}: more blocks touched than the plan allowed ({cap})")]
    Overflow {
        /// the part
        part: usize,
        /// its arena's blocks
        cap: usize,
    },
    /// A count that must be whole is not (orbits under a stabiliser).
    #[error("class {class}: orbit count {sum} not divisible by {group}")]
    Orbits {
        /// the class
        class: usize,
        /// the sum of fixed slots
        sum: u64,
        /// the stabiliser's order
        group: u64,
    },
}

fn io_err(path: &Path) -> impl FnOnce(io::Error) -> BfsError + '_ {
    move |source| BfsError::Io {
        path: path.display().to_string(),
        source,
    }
}

/// The counts of one layer.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Level {
    /// distance
    pub depth: usize,
    /// slots (the search's index entries)
    pub slots: u64,
    /// positions (each configuration of a class counted)
    pub positions: u64,
    /// positions up to the 24 symmetries
    pub orbits: u64,
    /// slots per class
    pub per_class: Vec<u64>,
    /// nonempty blocks per class
    pub nonempty: Vec<u64>,
    /// bytes of the layer's files
    pub bytes: u64,
    /// seconds to make the layer
    pub seconds: f64,
    /// bytes read to make it
    pub read: u64,
}

/// One class of a part's result.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ClassCount {
    /// the class
    pub class: u32,
    /// slots, nonempty blocks, file bytes
    pub sum: (u64, u64, u64),
    /// orbits under the class's stabiliser
    pub orbits: u64,
}

/// A recorded part.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PartDone {
    /// per class
    pub classes: Vec<ClassCount>,
    /// seconds
    pub seconds: f64,
    /// bytes read
    pub read: u64,
}

/// The step in progress.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Progress {
    /// the depth being made
    pub depth: usize,
    /// the parts: ranges of positions in [`Instance::order`]
    pub parts: Vec<(u32, u32)>,
    /// the parts recorded
    pub done: Vec<Option<PartDone>>,
}

/// The state of a search directory.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Manifest {
    /// format version
    pub version: u32,
    /// the instance
    pub kind: Kind,
    /// layers kept at the start (files of depths `0..=keep`), besides the last two
    pub keep: usize,
    /// the directories of the layer files (empty: the search directory)
    pub stores: Vec<String>,
    /// the layers made, from depth 0
    pub levels: Vec<Level>,
    /// whether the last layer is the last (the next one is empty)
    pub complete: bool,
    /// the step in progress
    pub progress: Option<Progress>,
}

/// The settings of a run.
#[derive(Clone, Debug)]
pub struct Config {
    /// the search directory
    pub dir: PathBuf,
    /// memory for a part (arena and directory), bytes
    pub ram: usize,
    /// stop once this depth is made
    pub max_depth: Option<usize>,
    /// stop after this many parts in this run (to test resuming)
    pub max_parts: Option<usize>,
    /// stop after the first part that ends after this many seconds
    pub max_seconds: Option<f64>,
    /// worker threads
    pub threads: usize,
    /// layers to keep at the start
    pub keep: usize,
    /// directories for the layer files (empty: the search directory); fixed when a search
    /// starts
    pub stores: Vec<PathBuf>,
    /// threads reading layer files during a push (several requests in flight)
    pub readers: usize,
    /// print progress
    pub verbose: bool,
}

/// How a run ended.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Outcome {
    /// the whole space is done
    Complete,
    /// stopped at a limit; run again to resume
    Stopped,
}

/// Where the layer files are: the directories of the stores, class `c`'s files in store
/// `c % stores` (`<store>/f<depth>/c<class>.blk`), so that the data and its traffic are spread
/// over the disks.
#[derive(Clone, Debug)]
pub struct Layout {
    /// the store directories
    pub stores: Vec<PathBuf>,
}

impl Layout {
    /// The layout of a search: its manifest's stores, or the search directory itself.
    pub fn of(dir: &Path, m: &Manifest) -> Layout {
        let stores: Vec<PathBuf> = m.stores.iter().map(PathBuf::from).collect();
        Layout {
            stores: if stores.is_empty() {
                vec![dir.to_path_buf()]
            } else {
                stores
            },
        }
    }

    fn store(&self, c: usize) -> &Path {
        self.stores
            .get(c % self.stores.len().max(1))
            .map_or(Path::new("."), PathBuf::as_path)
    }

    /// The directories of layer `d` (one per store).
    pub fn dirs(&self, d: usize) -> Vec<PathBuf> {
        self.stores
            .iter()
            .map(|s| s.join(format!("f{d:03}")))
            .collect()
    }

    /// The file of class `c` in layer `d`.
    pub fn file(&self, d: usize, c: usize) -> PathBuf {
        self.store(c)
            .join(format!("f{d:03}"))
            .join(format!("c{c:04}.blk"))
    }

    fn tmp(&self, d: usize, c: usize) -> PathBuf {
        self.store(c)
            .join(format!("f{d:03}"))
            .join(format!("c{c:04}.blk.tmp"))
    }

    /// Whether layer `d` is on disk.
    pub fn has(&self, d: usize) -> bool {
        self.dirs(d).iter().all(|p| p.exists())
    }

    fn create(&self, d: usize) -> Result<(), BfsError> {
        for p in self.dirs(d) {
            fs::create_dir_all(&p).map_err(io_err(&p))?;
        }
        Ok(())
    }

    fn remove(&self, d: usize) -> Result<(), BfsError> {
        for p in self.dirs(d) {
            if p.exists() {
                fs::remove_dir_all(&p).map_err(io_err(&p))?;
            }
        }
        Ok(())
    }

    fn sync(&self, d: usize) -> Result<(), BfsError> {
        for p in self.dirs(d) {
            if let Ok(f) = fs::File::open(&p) {
                f.sync_all().map_err(io_err(&p))?;
            }
        }
        Ok(())
    }
}

impl Manifest {
    fn path(dir: &Path) -> PathBuf {
        dir.join("manifest.json")
    }

    /// Read the manifest of a search directory (`None` if there is none).
    ///
    /// # Errors
    ///
    /// If it cannot be read or parsed.
    pub fn load(dir: &Path) -> Result<Option<Manifest>, BfsError> {
        let path = Self::path(dir);
        let text = match fs::read_to_string(&path) {
            Ok(t) => t,
            Err(e) if e.kind() == io::ErrorKind::NotFound => return Ok(None),
            Err(e) => return Err(io_err(&path)(e)),
        };
        let m: Manifest = serde_json::from_str(&text).map_err(|e| BfsError::Manifest {
            path: path.display().to_string(),
            what: e.to_string(),
        })?;
        if m.version != VERSION {
            return Err(BfsError::Manifest {
                path: path.display().to_string(),
                what: format!("version {} (this build reads {VERSION})", m.version),
            });
        }
        Ok(Some(m))
    }

    /// Write the manifest atomically (temporary file, sync, rename, sync the directory).
    ///
    /// # Errors
    ///
    /// If writing fails.
    pub fn save(&self, dir: &Path) -> Result<(), BfsError> {
        let path = Self::path(dir);
        let tmp = dir.join("manifest.json.tmp");
        let text = serde_json::to_string(self).map_err(|e| BfsError::Manifest {
            path: path.display().to_string(),
            what: e.to_string(),
        })?;
        (|| -> io::Result<()> {
            let mut f = fs::File::create(&tmp)?;
            io::Write::write_all(&mut f, text.as_bytes())?;
            f.sync_all()?;
            fs::rename(&tmp, &path)?;
            fs::File::open(dir)?.sync_all()
        })()
        .map_err(io_err(&path))
    }

    /// The depth of the last layer made.
    pub fn depth(&self) -> Option<usize> {
        self.levels.last().map(|l| l.depth)
    }
}

/// Blocks allocated on first touch, written by many threads with atomic ors.
struct Arena {
    words: usize,
    data: Vec<AtomicU64>,
    cap: usize,
    next: AtomicUsize,
    /// per block of the part (`(class - lo) * n_x + x`): its arena index, [`EMPTY`] or [`BUSY`]
    dir: Vec<AtomicU32>,
    overflow: AtomicBool,
}

const EMPTY: u32 = u32::MAX;
const BUSY: u32 = u32::MAX - 1;

// `AtomicU64` and `u64` have the same size and alignment on the targets this builds for, so a
// zeroed `Vec<u64>` can become a `Vec<AtomicU64>` (see `Arena::zeroed`).
const _: () = assert!(std::mem::align_of::<AtomicU64>() == std::mem::align_of::<u64>());
const _: () = assert!(std::mem::size_of::<AtomicU64>() == std::mem::size_of::<u64>());
const _: () = assert!(std::mem::align_of::<AtomicU32>() == std::mem::align_of::<u32>());

/// A vector of zeroed atomics, from zeroed memory the system gives lazily.
fn zeroed_u64(n: usize) -> Vec<AtomicU64> {
    let mut v = std::mem::ManuallyDrop::new(vec![0u64; n]);
    // SAFETY: AtomicU64 has the size and alignment of u64 (checked above at compile time) and
    // the same bit validity, so the allocation is a valid Vec<AtomicU64> of the same length and
    // capacity, freed with the same layout.
    unsafe { Vec::from_raw_parts(v.as_mut_ptr().cast::<AtomicU64>(), v.len(), v.capacity()) }
}

fn filled_u32(n: usize, value: u32) -> Vec<AtomicU32> {
    let mut v = std::mem::ManuallyDrop::new(vec![value; n]);
    // SAFETY: AtomicU32 has the size and alignment of u32 (checked above at compile time) and
    // the same bit validity, so this is a valid Vec<AtomicU32> freed with the same layout.
    unsafe { Vec::from_raw_parts(v.as_mut_ptr().cast::<AtomicU32>(), v.len(), v.capacity()) }
}

/// Ask for transparent huge pages on the 2 MiB-aligned inside of a buffer (its blocks are
/// touched at random: huge pages save most of the TLB misses).
#[cfg(all(target_os = "linux", not(miri)))]
fn advise_huge<T>(v: &[T]) {
    const HUGE: usize = 2 << 20;
    let start = v.as_ptr() as usize;
    let end = start + std::mem::size_of_val(v);
    let a = start.next_multiple_of(HUGE);
    if end > a + HUGE {
        let len = (end - a) / HUGE * HUGE;
        // SAFETY: an advisory call on whole pages inside memory this buffer owns; it does
        // not change the contents
        let p = v.as_ptr().cast::<u8>().wrapping_add(a - start).cast_mut();
        let _ = unsafe { libc::madvise(p.cast::<libc::c_void>(), len, libc::MADV_HUGEPAGE) };
    }
}

#[cfg(not(all(target_os = "linux", not(miri))))]
fn advise_huge<T>(_v: &[T]) {}

impl Arena {
    fn new(blocks_in_part: usize, cap: usize, words: usize) -> Arena {
        let data = zeroed_u64(cap * words);
        advise_huge(&data);
        let dir = filled_u32(blocks_in_part, EMPTY);
        Arena {
            words,
            data,
            cap,
            next: AtomicUsize::new(0),
            dir,
            overflow: AtomicBool::new(false),
        }
    }

    /// The arena index of block `i` of the part, allocated if new (`None` on overflow).
    #[inline]
    fn block(&self, i: usize) -> Option<usize> {
        let d = self.dir.get(i)?;
        let v = d.load(Ordering::Acquire);
        if v < BUSY {
            return Some(v as usize);
        }
        if v == EMPTY
            && d.compare_exchange(EMPTY, BUSY, Ordering::AcqRel, Ordering::Acquire)
                .is_ok()
        {
            let b = self.next.fetch_add(1, Ordering::Relaxed);
            if b >= self.cap {
                self.overflow.store(true, Ordering::Relaxed);
                // (left BUSY: the part is abandoned)
                return None;
            }
            d.store(b as u32, Ordering::Release);
            return Some(b);
        }
        loop {
            let v = d.load(Ordering::Acquire);
            if v < BUSY {
                return Some(v as usize);
            }
            if self.overflow.load(Ordering::Relaxed) {
                return None;
            }
            std::hint::spin_loop();
        }
    }

    /// The arena index of block `i` if it was touched.
    #[inline]
    fn get(&self, i: usize) -> Option<usize> {
        let v = self.dir.get(i)?.load(Ordering::Acquire);
        (v < BUSY).then_some(v as usize)
    }

    #[inline]
    fn or(&self, b: usize, w: usize, v: u64) {
        if let Some(a) = self.data.get(b * self.words + w) {
            a.fetch_or(v, Ordering::Relaxed);
        }
    }

    fn load(&self, b: usize, out: &mut [u64]) {
        let base = b * self.words;
        for (w, o) in out.iter_mut().enumerate() {
            *o = self
                .data
                .get(base + w)
                .map_or(0, |a| a.load(Ordering::Relaxed));
        }
    }
}

/// The search over one directory.
pub struct Bfs<'a> {
    inst: &'a Instance,
    cfg: Config,
    m: Manifest,
    lay: Layout,
    /// per class: the classes with a move into it
    pred: Vec<Vec<u32>>,
    /// per symmetry: the B-centre states and the corner-state mask it fixes
    fixed: Vec<(Vec<u32>, Vec<u64>)>,
}

impl<'a> Bfs<'a> {
    /// Open a search directory: resume its manifest, or start a new search (layer 0: the
    /// solved slot).
    ///
    /// # Errors
    ///
    /// If the manifest belongs to another instance or a file cannot be written.
    pub fn open(inst: &'a Instance, cfg: Config) -> Result<Bfs<'a>, BfsError> {
        fs::create_dir_all(&cfg.dir).map_err(io_err(&cfg.dir))?;
        let m = match Manifest::load(&cfg.dir)? {
            Some(m) => {
                if m.kind != inst.kind {
                    return Err(BfsError::Manifest {
                        path: Manifest::path(&cfg.dir).display().to_string(),
                        what: format!("a search of {:?}, not {:?}", m.kind, inst.kind),
                    });
                }
                m
            }
            None => start(inst, &cfg)?,
        };
        if !cfg.stores.is_empty() {
            let want: Vec<String> = cfg.stores.iter().map(|p| p.display().to_string()).collect();
            if want != m.stores {
                return Err(BfsError::Manifest {
                    path: Manifest::path(&cfg.dir).display().to_string(),
                    what: format!("a search with the stores {:?}, not {want:?}", m.stores),
                });
            }
        }
        let lay = Layout::of(&cfg.dir, &m);
        let mut pred = vec![vec![]; inst.n_class];
        for c in 0..inst.n_class {
            for &(c2, _) in &inst.edges[c] {
                if c2 != NO_EDGE && !pred[c2 as usize].contains(&(c as u32)) {
                    pred[c2 as usize].push(c as u32);
                }
            }
        }
        let fixed = (0..N_SYM)
            .map(|k| {
                let xs: Vec<u32> = (0..inst.n_x as u32)
                    .filter(|&x| inst.xs[k * inst.n_x + x as usize] == x)
                    .collect();
                let mut mask = vec![0u64; inst.words];
                for p in 0..inst.n_p {
                    if inst.ps[p * N_SYM + k] as usize == p {
                        mask[p / 64] |= 1 << (p % 64);
                    }
                }
                (xs, mask)
            })
            .collect();
        Ok(Bfs {
            inst,
            cfg,
            m,
            lay,
            pred,
            fixed,
        })
    }

    /// The manifest.
    pub fn manifest(&self) -> &Manifest {
        &self.m
    }

    fn log(&self, msg: &str) {
        if self.cfg.verbose {
            eprintln!("[{}] {msg}", clock());
        }
    }

    /// Run until the search is complete or a limit is reached.
    ///
    /// # Errors
    ///
    /// If a file cannot be read or written, or a check fails.
    pub fn run(&mut self) -> Result<Outcome, BfsError> {
        let t_run = Instant::now();
        let mut parts_run = 0usize;
        loop {
            if self.m.complete {
                return Ok(Outcome::Complete);
            }
            let d = self.m.depth().unwrap_or(0);
            if self.m.progress.is_none() {
                self.cleanup(d)?;
                if self.cfg.max_depth.is_some_and(|max| d >= max) {
                    return Ok(Outcome::Stopped);
                }
                let parts = self.plan(d);
                self.log(&format!(
                    "depth {}: {} part(s) {:?}",
                    d + 1,
                    parts.len(),
                    if parts.len() <= 8 {
                        parts.clone()
                    } else {
                        vec![]
                    }
                ));
                self.m.progress = Some(Progress {
                    depth: d + 1,
                    done: vec![None; parts.len()],
                    parts,
                });
                self.m.save(&self.cfg.dir)?;
            }
            let Some(prog) = self.m.progress.clone() else {
                continue;
            };
            self.lay.create(prog.depth)?;
            // finish the renames of parts recorded before a stop
            for (j, done) in prog.done.iter().enumerate() {
                if done.is_some() {
                    self.rename_part(prog.depth, prog.parts[j])?;
                }
            }
            for j in 0..prog.parts.len() {
                if prog.done[j].is_some() {
                    continue;
                }
                if self.cfg.max_parts.is_some_and(|max| parts_run >= max)
                    || self
                        .cfg
                        .max_seconds
                        .is_some_and(|s| t_run.elapsed().as_secs_f64() >= s)
                {
                    return Ok(Outcome::Stopped);
                }
                let r = self.part(prog.depth, j, prog.parts[j])?;
                if let Some(p) = self.m.progress.as_mut() {
                    p.done[j] = Some(r);
                }
                self.m.save(&self.cfg.dir)?;
                self.rename_part(prog.depth, prog.parts[j])?;
                parts_run += 1;
            }
            self.finish_level()?;
        }
    }

    /// Remove the layer directories no step needs once layer `d` is the last made (left by a
    /// stop between the manifest and the removal): all but `d - 1`, `d` and the kept ones.
    fn cleanup(&self, d: usize) -> Result<(), BfsError> {
        for store in &self.lay.stores {
            let rd = fs::read_dir(store).map_err(io_err(store))?;
            for e in rd.flatten() {
                let name = e.file_name();
                let Some(k) = name
                    .to_str()
                    .and_then(|n| n.strip_prefix('f'))
                    .and_then(|n| n.parse::<usize>().ok())
                else {
                    continue;
                };
                if (k > self.m.keep && k + 1 < d) || k > d {
                    let path = e.path();
                    fs::remove_dir_all(&path).map_err(io_err(&path))?;
                }
            }
        }
        Ok(())
    }

    /// Upper bound of the blocks of layer `d + 1` per class: the blocks of layer `d` that a
    /// move takes into the class (times the class's stabiliser, for the closure).
    fn estimate(&self, d: usize) -> Vec<usize> {
        let inst = self.inst;
        let mut from = vec![0u64; inst.n_class];
        if let Some(l) = self.m.levels.get(d) {
            for c in 0..inst.n_class {
                let ne = l.nonempty.get(c).copied().unwrap_or(0);
                for &(c2, _) in &inst.edges[c] {
                    if c2 != NO_EDGE {
                        from[c2 as usize] += ne;
                    }
                }
            }
        }
        from.iter()
            .enumerate()
            .map(|(c, &f)| {
                let bound = f.saturating_mul(1 + inst.stab[c].len() as u64);
                bound.min(inst.n_x as u64) as usize
            })
            .collect()
    }

    /// Split the classes into parts that fit the memory budget: ranges of positions in
    /// [`Instance::order`].
    fn plan(&self, d: usize) -> Vec<(u32, u32)> {
        let inst = self.inst;
        let est = self.estimate(d);
        let per_block = inst.words * 8;
        let mut parts = vec![];
        let mut lo = 0usize;
        let mut used = 0usize;
        for (i, &c) in inst.order.iter().enumerate() {
            let cost = est[c as usize] * per_block + inst.n_x * 4;
            if i > lo && used + cost > self.cfg.ram {
                parts.push((lo as u32, i as u32));
                lo = i;
                used = 0;
            }
            used += cost;
        }
        parts.push((lo as u32, inst.n_class as u32));
        parts
    }

    /// Make part `j` of depth `depth` (classes `lo..hi`): push, close, mask, write.
    fn part(&self, depth: usize, j: usize, (lo, hi): (u32, u32)) -> Result<PartDone, BfsError> {
        let t0 = Instant::now();
        let inst = self.inst;
        let (lo, hi) = (lo as usize, hi as usize);
        let est = self.estimate(depth - 1);
        let classes = &inst.order[lo..hi];
        let cap: usize = classes.iter().map(|&c| est[c as usize]).sum();
        let arena = Arena::new((hi - lo) * inst.n_x, cap.max(1), inst.words);
        let read = AtomicU64::new(0);
        // sources: the classes with a move into the part and a nonempty layer
        let ne = &self.m.levels[depth - 1].nonempty;
        let mut sources: Vec<usize> = classes
            .iter()
            .flat_map(|&c| self.pred[c as usize].iter().map(|&s| s as usize))
            .filter(|&s| ne.get(s).copied().unwrap_or(0) > 0)
            .collect();
        sources.sort_unstable();
        sources.dedup();
        self.push(depth - 1, &sources, lo, hi, &arena, &read)?;
        if arena.overflow.load(Ordering::Relaxed) {
            return Err(BfsError::Overflow { part: j, cap });
        }
        let t_push = t0.elapsed().as_secs_f64();
        self.close(lo, hi, &arena);
        if arena.overflow.load(Ordering::Relaxed) {
            return Err(BfsError::Overflow { part: j, cap });
        }
        let counts = self.mask_write(depth, lo, hi, &arena, &read)?;
        let slots: u64 = counts.iter().map(|c| c.sum.0).sum();
        self.log(&format!(
            "depth {depth} part {j} ({} classes): {} sources, {} blocks touched, {slots} \
             slots, push {t_push:.1} s, total {:.1} s, read {:.2} GB",
            hi - lo,
            sources.len(),
            arena.next.load(Ordering::Relaxed),
            t0.elapsed().as_secs_f64(),
            read.load(Ordering::Relaxed) as f64 / 1e9
        ));
        Ok(PartDone {
            classes: counts,
            seconds: t0.elapsed().as_secs_f64(),
            read: read.load(Ordering::Relaxed),
        })
    }

    /// Push the blocks of layer `d` of `sources` into the part's arena.
    fn push(
        &self,
        d: usize,
        sources: &[usize],
        lo: usize,
        hi: usize,
        arena: &Arena,
        read: &AtomicU64,
    ) -> Result<(), BfsError> {
        let inst = self.inst;
        let err: Mutex<Option<BfsError>> = Mutex::new(None);
        let set_err = |e: BfsError| {
            let mut g = err.lock().unwrap_or_else(PoisonError::into_inner);
            g.get_or_insert(e);
        };
        let failed = AtomicBool::new(false);
        let (tx, rx) = sync_channel::<(usize, usize, Vec<u8>)>(2 * self.cfg.threads.max(1));
        let rx = Mutex::new(rx);
        std::thread::scope(|s| {
            // the readers: each its share of the files in order, chunk by chunk
            let readers = self.cfg.readers.clamp(1, sources.len().max(1));
            for r in 0..readers {
                let tx = tx.clone();
                let failed = &failed;
                let set_err = &set_err;
                s.spawn(move || {
                    for &c in sources.iter().skip(r).step_by(readers) {
                        let path = self.lay.file(d, c);
                        let rd = match ClassReader::open(&path, inst.n_x, inst.n_p) {
                            Ok(rd) => rd,
                            Err(e) => {
                                set_err(io_err(&path)(e));
                                failed.store(true, Ordering::Relaxed);
                                return;
                            }
                        };
                        for i in 0..rd.chunks() {
                            if failed.load(Ordering::Relaxed) {
                                return;
                            }
                            let mut buf = Vec::new();
                            if let Err(e) = rd.read_chunk(i, &mut buf) {
                                set_err(io_err(&path)(e));
                                failed.store(true, Ordering::Relaxed);
                                return;
                            }
                            read.fetch_add(buf.len() as u64, Ordering::Relaxed);
                            if tx.send((c, i, buf)).is_err() {
                                return;
                            }
                        }
                    }
                });
            }
            drop(tx);
            for _ in 0..self.cfg.threads.max(1) {
                s.spawn(|| {
                    let mut bits: Vec<u16> = Vec::with_capacity(inst.n_p);
                    let mut words = vec![0u64; inst.words];
                    let mut buf = vec![0u64; inst.words];
                    loop {
                        let job = {
                            let g = rx.lock().unwrap_or_else(PoisonError::into_inner);
                            g.recv()
                        };
                        let Ok((c, i, bytes)) = job else { return };
                        if failed.load(Ordering::Relaxed) {
                            continue;
                        }
                        let x0 = i * CHUNK;
                        let n = CHUNK.min(inst.n_x - x0);
                        let ok = decode_chunk(&bytes, x0, n, inst.words, false, |x, blk| {
                            let src = match blk {
                                Blk::Full => Src::Full,
                                Blk::Raw(_) => {
                                    blk.expand(&mut words, inst.n_p);
                                    Src::Words(&words)
                                }
                                Blk::Empty | Blk::Sparse(_) => {
                                    bits.clear();
                                    blk.bits(&mut bits, inst.n_p);
                                    Src::Bits(&bits)
                                }
                            };
                            push_block(inst, c, x, src, lo, hi, arena, &mut buf);
                        });
                        if ok.is_none() {
                            set_err(BfsError::Io {
                                path: self.lay.file(d, c).display().to_string(),
                                source: io::Error::new(
                                    io::ErrorKind::InvalidData,
                                    format!("chunk {i} malformed"),
                                ),
                            });
                            failed.store(true, Ordering::Relaxed);
                        }
                    }
                });
            }
        });
        match err.into_inner().unwrap_or_else(PoisonError::into_inner) {
            Some(e) => Err(e),
            None => Ok(()),
        }
    }

    /// Close the part's new blocks under each class's stabiliser.
    fn close(&self, lo: usize, hi: usize, arena: &Arena) {
        let inst = self.inst;
        let classes: Vec<(usize, usize)> = (lo..hi)
            .map(|i| (i, inst.order[i] as usize))
            .filter(|&(_, c)| !inst.stab[c].is_empty())
            .collect();
        rayon::scope(|s| {
            for &(i, c) in &classes {
                s.spawn(move |_| {
                    let mut blk = vec![0u64; inst.words];
                    let mut out = vec![0u64; inst.words];
                    for x in 0..inst.n_x {
                        let Some(b) = arena.get((i - lo) * inst.n_x + x) else {
                            continue;
                        };
                        arena.load(b, &mut blk);
                        for &k in &inst.stab[c] {
                            let k = k as usize;
                            let x2 = inst.xs[k * inst.n_x + x] as usize;
                            out.fill(0);
                            for (w, &v) in blk.iter().enumerate() {
                                let mut v = v;
                                while v != 0 {
                                    let p = w * 64 + v.trailing_zeros() as usize;
                                    v &= v - 1;
                                    let t = inst.ps[p * N_SYM + k] as usize;
                                    out[t / 64] |= 1 << (t % 64);
                                }
                            }
                            let Some(b2) = arena.block((i - lo) * inst.n_x + x2) else {
                                return;
                            };
                            for (w, &v) in out.iter().enumerate() {
                                if v != 0 {
                                    arena.or(b2, w, v);
                                }
                            }
                        }
                    }
                });
            }
        });
    }

    /// Write layer `depth` of the part's classes: the new blocks minus layers `depth - 1` and
    /// `depth - 2`; count them.
    fn mask_write(
        &self,
        depth: usize,
        lo: usize,
        hi: usize,
        arena: &Arena,
        read: &AtomicU64,
    ) -> Result<Vec<ClassCount>, BfsError> {
        use rayon::prelude::*;
        let inst = self.inst;
        (lo..hi)
            .into_par_iter()
            .map(|at| -> Result<ClassCount, BfsError> {
                // the class at position `at` of the order; its blocks in the arena from
                // `(at - lo) * n_x`
                let c = inst.order[at] as usize;
                let base = (at - lo) * inst.n_x;
                let open = |d: usize| -> Result<Option<ClassReader>, BfsError> {
                    let path = self.lay.file(d, c);
                    ClassReader::open(&path, inst.n_x, inst.n_p)
                        .map(Some)
                        .map_err(io_err(&path))
                };
                let prev = open(depth - 1)?;
                let prev2 = if depth >= 2 { open(depth - 2)? } else { None };
                let out = self.lay.tmp(depth, c);
                let mut w = ClassWriter::create(&out, inst.n_x, inst.n_p).map_err(io_err(&out))?;
                let mut mask = vec![0u64; CHUNK * inst.words];
                let mut blk = vec![0u64; inst.words];
                let mut scratch = vec![0u64; inst.words];
                let mut buf = vec![];
                let mut fixed_sum = 0u64;
                let n_chunks = inst.n_x.div_ceil(CHUNK);
                for i in 0..n_chunks {
                    let x0 = i * CHUNK;
                    let n = CHUNK.min(inst.n_x - x0);
                    // any new block in this chunk?
                    let any = (x0..x0 + n).any(|x| arena.get(base + x).is_some());
                    if !any {
                        w.push_empty(n).map_err(io_err(&out))?;
                        continue;
                    }
                    mask[..n * inst.words].fill(0);
                    for (dd, r) in [(depth - 1, &prev), (depth.wrapping_sub(2), &prev2)] {
                        let Some(r) = r else { continue };
                        let path = self.lay.file(dd, c);
                        r.for_chunk(i, &mut buf, false, |x, b| {
                            let s = (x - x0) * inst.words;
                            b.expand(&mut scratch, inst.n_p);
                            for (m, v) in mask[s..s + inst.words].iter_mut().zip(&scratch) {
                                *m |= v;
                            }
                        })
                        .map_err(io_err(&path))?;
                        read.fetch_add(buf.len() as u64, Ordering::Relaxed);
                    }
                    for x in x0..x0 + n {
                        match arena.get(base + x) {
                            None => w.push_empty(1).map_err(io_err(&out))?,
                            Some(b) => {
                                arena.load(b, &mut blk);
                                let s = (x - x0) * inst.words;
                                for (v, m) in blk.iter_mut().zip(&mask[s..s + inst.words]) {
                                    *v &= !m;
                                }
                                for &k in &inst.stab[c] {
                                    let (xs, pm) = &self.fixed[k as usize];
                                    if xs.binary_search(&(x as u32)).is_ok() {
                                        fixed_sum += blk
                                            .iter()
                                            .zip(pm)
                                            .map(|(a, b)| u64::from((a & b).count_ones()))
                                            .sum::<u64>();
                                    }
                                }
                                w.push(&blk).map_err(io_err(&out))?;
                            }
                        }
                    }
                }
                let sum: Summary = w.finish().map_err(io_err(&out))?;
                let group = 1 + inst.stab[c].len() as u64;
                let total = sum.bits + fixed_sum;
                if !total.is_multiple_of(group) {
                    return Err(BfsError::Orbits {
                        class: c,
                        sum: total,
                        group,
                    });
                }
                Ok(ClassCount {
                    class: c as u32,
                    sum: (sum.bits, sum.nonempty, sum.bytes),
                    orbits: total / group,
                })
            })
            .collect()
    }

    /// Rename a recorded part's files to their final names (when not done yet), and delete
    /// the part's files of layer `depth - 2` (no step needs them any more, unless kept).
    fn rename_part(&self, depth: usize, (lo, hi): (u32, u32)) -> Result<(), BfsError> {
        for &c in &self.inst.order[lo as usize..hi as usize] {
            let c = c as usize;
            let tmp = self.lay.tmp(depth, c);
            if tmp.exists() {
                let fin = self.lay.file(depth, c);
                fs::rename(&tmp, &fin).map_err(io_err(&fin))?;
            }
            if depth >= 2 && depth - 2 > self.m.keep {
                let old = self.lay.file(depth - 2, c);
                if old.exists() {
                    fs::remove_file(&old).map_err(io_err(&old))?;
                }
            }
        }
        self.lay.sync(depth)
    }

    /// All parts recorded: make the level, drop the layer no longer needed.
    fn finish_level(&mut self) -> Result<(), BfsError> {
        let inst = self.inst;
        let Some(prog) = self.m.progress.take() else {
            return Ok(());
        };
        let mut lv = Level {
            depth: prog.depth,
            per_class: vec![0; inst.n_class],
            nonempty: vec![0; inst.n_class],
            ..Level::default()
        };
        for p in prog.done.iter().flatten() {
            lv.seconds += p.seconds;
            lv.read += p.read;
            for cc in &p.classes {
                let c = cc.class as usize;
                lv.per_class[c] = cc.sum.0;
                lv.nonempty[c] = cc.sum.1;
                lv.bytes += cc.sum.2;
                lv.slots += cc.sum.0;
                lv.positions += inst.class_size(c) * cc.sum.0;
                lv.orbits += cc.orbits;
            }
        }
        self.log(&format!(
            "depth {}: {} slots, {} positions, {} up to symmetry, {:.2} GB, {:.1} s",
            lv.depth,
            lv.slots,
            lv.positions,
            lv.orbits,
            lv.bytes as f64 / 1e9,
            lv.seconds
        ));
        let last = lv.slots == 0;
        if last {
            // the previous layer was the last: keep it (the antipodes)
            self.m.complete = true;
        } else {
            self.m.levels.push(lv);
        }
        self.m.save(&self.cfg.dir)?;
        if last {
            self.lay.remove(prog.depth)?;
        }
        // layer depth - 2 is no longer needed (unless kept)
        if prog.depth >= 2 && prog.depth - 2 > self.m.keep {
            self.lay.remove(prog.depth - 2)?;
        }
        Ok(())
    }
}

/// A source block for [`push_block`].
#[derive(Clone, Copy)]
enum Src<'a> {
    /// every bit set
    Full,
    /// the indices of the bits set (a sparse block)
    Bits(&'a [u16]),
    /// the words of a dense block
    Words(&'a [u64]),
}

/// Push one block (class `c`, B-state `x`) into the part (positions `lo..hi` of the order).
/// A dense block of the full instance is permuted a 32-bit half at a time
/// ([`super::instance::WordPerm`]), and the children's directory entries and blocks are
/// prefetched before the ORs: with huge pages for the arena, 1.96 ns per slot move per thread
/// instead of 4.29 ns one bit at a time (blocks of density 0.4, `bench-push`, 48 threads; about
/// 9 us per block and move, mostly the 180 atomic ORs into a block not in cache).
#[inline]
#[allow(clippy::too_many_arguments)]
fn push_block(
    inst: &Instance,
    c: usize,
    x: usize,
    src: Src,
    lo: usize,
    hi: usize,
    arena: &Arena,
    buf: &mut [u64],
) {
    let xrow = inst.xm[x];
    // the children in the part: (move, symmetry, directory index), their directory entries
    // prefetched; then their blocks, prefetched; then the permutations (the cache misses of
    // the children overlap instead of following each other)
    let mut kids = [(0u8, 0u8, 0usize); N_MOVES];
    let mut n = 0;
    for m in 0..N_MOVES {
        let (c2, k) = inst.edges[c][m];
        if c2 == NO_EDGE {
            continue;
        }
        // the child's class's position in the part order
        let p2 = inst.pos[c2 as usize] as usize;
        if p2 < lo || p2 >= hi {
            continue;
        }
        let x2 = inst.xs[k as usize * inst.n_x + xrow[m] as usize] as usize;
        let i = (p2 - lo) * inst.n_x + x2;
        prefetch_at(&arena.dir, i);
        kids[n] = (m as u8, k, i);
        n += 1;
    }
    let mut blocks = [0usize; N_MOVES];
    for j in 0..n {
        let Some(b) = arena.block(kids[j].2) else {
            return;
        };
        blocks[j] = b;
        for w in (0..inst.words).step_by(8) {
            prefetch_at(&arena.data, b * inst.words + w);
        }
    }
    for j in 0..n {
        let (m, k, _) = kids[j];
        let (m, k, b) = (m as usize, k as usize, blocks[j]);
        let perm = inst.fiber_perm(m, k);
        match src {
            Src::Full => {
                let last = inst.words - 1;
                for w in 0..inst.words {
                    let v = if w == last {
                        inst.last_mask()
                    } else {
                        u64::MAX
                    };
                    arena.or(b, w, v);
                }
                continue;
            }
            Src::Bits(bits) => {
                for &i in bits {
                    let t = perm[i as usize] as usize;
                    buf[t / 64] |= 1 << (t % 64);
                }
            }
            Src::Words(words) => match &inst.word {
                Some(wp) if wp.ok[m * N_SYM + k] => wp.apply(m * N_SYM + k, words, buf),
                _ => {
                    for (w, &v) in words.iter().enumerate() {
                        let mut v = v;
                        while v != 0 {
                            let t = perm[w * 64 + v.trailing_zeros() as usize] as usize;
                            buf[t / 64] |= 1 << (t % 64);
                            v &= v - 1;
                        }
                    }
                }
            },
        }
        for (w, v) in buf.iter_mut().enumerate() {
            if *v != 0 {
                arena.or(b, w, *v);
                *v = 0;
            }
        }
    }
}

/// A new search: layer 0 is the solved slot (closed under its stabiliser: it is fixed).
fn start(inst: &Instance, cfg: &Config) -> Result<Manifest, BfsError> {
    let (c0, x0, p0) = inst.solved;
    let stores: Vec<String> = cfg.stores.iter().map(|p| p.display().to_string()).collect();
    let lay = Layout {
        stores: if cfg.stores.is_empty() {
            vec![cfg.dir.clone()]
        } else {
            cfg.stores.clone()
        },
    };
    lay.create(0)?;
    let mut lv = Level {
        depth: 0,
        per_class: vec![0; inst.n_class],
        nonempty: vec![0; inst.n_class],
        ..Level::default()
    };
    // the solved slot and its images under the stabiliser
    let mut slots = vec![(x0 as usize, p0 as usize)];
    for &k in &inst.stab[c0 as usize] {
        let k = k as usize;
        let s = (
            inst.xs[k * inst.n_x + x0 as usize] as usize,
            inst.ps[p0 as usize * N_SYM + k] as usize,
        );
        if !slots.contains(&s) {
            slots.push(s);
        }
    }
    slots.sort_unstable();
    for c in 0..inst.n_class {
        let path = lay.file(0, c);
        let mut w = ClassWriter::create(&path, inst.n_x, inst.n_p).map_err(io_err(&path))?;
        if c == c0 as usize {
            let mut x = 0;
            let mut blk = vec![0u64; inst.words];
            for &(sx, _) in &slots {
                if sx < x {
                    continue;
                }
                w.push_empty(sx - x).map_err(io_err(&path))?;
                blk.fill(0);
                for &(_, p) in slots.iter().filter(|s| s.0 == sx) {
                    blk[p / 64] |= 1 << (p % 64);
                }
                w.push(&blk).map_err(io_err(&path))?;
                x = sx + 1;
            }
            w.push_empty(inst.n_x - x).map_err(io_err(&path))?;
        } else {
            w.push_empty(inst.n_x).map_err(io_err(&path))?;
        }
        let sum = w.finish().map_err(io_err(&path))?;
        lv.per_class[c] = sum.bits;
        lv.nonempty[c] = sum.nonempty;
        lv.bytes += sum.bytes;
        lv.slots += sum.bits;
        lv.positions += inst.class_size(c) * sum.bits;
    }
    lv.orbits = 1;
    let m = Manifest {
        version: VERSION,
        kind: inst.kind,
        keep: cfg.keep,
        stores,
        levels: vec![lv],
        complete: false,
        progress: None,
    };
    m.save(&cfg.dir)?;
    Ok(m)
}

/// The time of day (UTC) for progress lines.
fn clock() -> String {
    let t = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |d| d.as_secs());
    format!("{:02}:{:02}:{:02}", t / 3600 % 24, t / 60 % 60, t % 60)
}

/// Call `f(c, x, p)` for every slot set in layer `d` of the classes `classes`.
///
/// # Errors
///
/// If a file cannot be read.
pub fn for_slots(
    lay: &Layout,
    inst: &Instance,
    d: usize,
    classes: impl IntoIterator<Item = usize>,
    mut f: impl FnMut(usize, usize, usize),
) -> Result<(), BfsError> {
    let mut buf = vec![];
    let mut bits: Vec<u16> = vec![];
    for c in classes {
        let path = lay.file(d, c);
        let r = ClassReader::open(&path, inst.n_x, inst.n_p).map_err(io_err(&path))?;
        if r.summary().bits == 0 {
            continue;
        }
        for i in 0..r.chunks() {
            r.for_chunk(i, &mut buf, false, |x, b| {
                bits.clear();
                b.bits(&mut bits, inst.n_p);
                for &p in &bits {
                    f(c, x, p as usize);
                }
            })
            .map_err(io_err(&path))?;
        }
    }
    Ok(())
}

/// Whether slot `(c, x, p)` is the smallest of its images under its class's stabiliser (one
/// slot per position up to symmetry).
pub fn canonical(inst: &Instance, c: usize, x: usize, p: usize) -> bool {
    inst.stab[c].iter().all(|&k| {
        let k = k as usize;
        let img = (
            inst.xs[k * inst.n_x + x] as usize,
            inst.ps[p * N_SYM + k] as usize,
        );
        img >= (x, p)
    })
}

/// The parts of a step whose next layer touches every block, for a memory budget, and the
/// class files read for them (the classes with a move into each part, summed over the parts):
/// the step at the peak of the search.
pub fn dense_plan(inst: &Instance, ram: usize) -> (Vec<(u32, u32)>, usize) {
    let per_class = inst.n_x * (inst.words * 8 + 4);
    let per_part = (ram / per_class.max(1)).max(1);
    let mut parts = vec![];
    let mut lo = 0;
    while lo < inst.n_class {
        let hi = (lo + per_part).min(inst.n_class);
        parts.push((lo as u32, hi as u32));
        lo = hi;
    }
    let reads = parts
        .iter()
        .map(|&(lo, hi)| {
            (0..inst.n_class)
                .filter(|&c| {
                    inst.edges[c]
                        .iter()
                        .any(|&(c2, _)| c2 != NO_EDGE && (lo..hi).contains(&inst.pos[c2 as usize]))
                })
                .count()
        })
        .sum();
    (parts, reads)
}

/// Time the push of random blocks of density `density` from the first `n_src` classes into
/// a part of the classes `lo..hi` (in memory, all threads): returns (slot moves pushed, blocks
/// touched, seconds). For estimates of the busiest steps.
pub fn bench_push(
    inst: &Instance,
    n_src: usize,
    blocks_per_class: usize,
    density: f64,
    (lo, hi): (usize, usize),
) -> (u64, usize, f64) {
    use rayon::prelude::*;
    let cap = (hi - lo) * inst.n_x;
    let arena = Arena::new((hi - lo) * inst.n_x, cap, inst.words);
    #[allow(clippy::cast_sign_loss, reason = "the density is clamped to 0..1")]
    let threshold = (density.clamp(0.0, 1.0) * f64::from(u32::MAX)) as u64;
    let jobs: Vec<(usize, usize)> = (0..n_src.min(inst.n_class))
        .flat_map(|c| {
            (0..blocks_per_class.min(inst.n_x))
                .map(move |i| (c, i * (inst.n_x / blocks_per_class.max(1)).max(1)))
        })
        .collect();
    // 64 random blocks (made before the clock starts)
    let mut h: u64 = 0x9e37_79b9_7f4a_7c15;
    let patterns: Vec<Vec<u16>> = (0..64)
        .map(|_| {
            (0..inst.n_p as u16)
                .filter(|_| {
                    h ^= h << 13;
                    h ^= h >> 7;
                    h ^= h << 17;
                    h & 0xffff_ffff < threshold
                })
                .collect()
        })
        .collect();
    let pattern_words: Vec<Vec<u64>> = patterns
        .iter()
        .map(|b| {
            let mut w = vec![0u64; inst.words];
            for &i in b {
                w[i as usize / 64] |= 1 << (i % 64);
            }
            w
        })
        .collect();
    let pass = || -> u64 {
        jobs.par_iter()
            .map_init(
                || vec![0u64; inst.words],
                |buf, &(c, x)| {
                    let j = (c * 31 + x) % 64;
                    let bits = &patterns[j];
                    // dense blocks as words, as the search decodes them
                    let src = if bits.len() >= 4 * inst.words - 2 {
                        Src::Words(&pattern_words[j])
                    } else {
                        Src::Bits(bits)
                    };
                    push_block(inst, c, x, src, lo, hi, &arena, buf);
                    let moves = inst.edges[c]
                        .iter()
                        .filter(|&&(c2, _)| {
                            c2 != NO_EDGE && (lo..hi).contains(&(inst.pos[c2 as usize] as usize))
                        })
                        .count() as u64;
                    bits.len() as u64 * moves
                },
            )
            .sum()
    };
    // a first pass touches the blocks (page faults), the second is timed
    pass();
    let t0 = Instant::now();
    let moved = pass();
    (
        moved,
        arena.next.load(Ordering::Relaxed),
        t0.elapsed().as_secs_f64(),
    )
}
