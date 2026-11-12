//! omelettehedron as a WebAssembly module: a small C ABI that `omelettehedron.mjs` wraps for
//! JavaScript (browsers and Bun).
//!
//! Strings go in as (pointer, length) of UTF-8 in buffers from `fto_alloc`. Results come out
//! as text in an output buffer (`fto_output`, `fto_output_len`): JSON after a call that
//! returns 0, an error message after one that returns -1. Table files (the format of the
//! native `gen`, `Table2::save`) are copied in by the caller into buffers from
//! `fto_table_buffer`, then checked and installed by `fto_load`; `fto_export` hands out the
//! generated tables in the same format.
//!
//! Three sets of tables (the native tiers of the same names), chosen by a number: 0 for the
//! small tier (the 129 MB A-centre/edge table without edge chirality), 1 for the medium tier
//! (the full 1.03 GB A-centre/edge table, 4.13e9 entries: its indices fit 32 bits) and 2 for
//! the tiny tier (the 37 MB pair and split tables of `omelettehedron::light`), each with the
//! 4.7 MB corner/B-centre table. The large tier needs 8 GB, beyond wasm32's 4 GB.
//!
//! Two builds (scripts/build-wasm.sh): `wasm32-wasip1` runs on one thread anywhere, and
//! `wasm32-wasip1-threads` runs rayon's pool on wasi-threads (in a browser: a cross-origin
//! isolated page). In both, the pool is rayon's global one, built once by `fto_init` with the
//! calling thread as one of its workers. It is never dropped, so its threads never exit: a
//! thread that exits can race with wasi-libc's thread exit.

#![cfg_attr(
    not(test),
    deny(
        clippy::unwrap_used,
        clippy::expect_used,
        clippy::panic,
        clippy::unreachable,
        clippy::todo,
        clippy::unimplemented
    )
)]
#![warn(clippy::pedantic)]
// sizes and counts are reported as f64 / u64 (wasm32: usize is 32 bits, so the casts are exact)
#![allow(clippy::cast_precision_loss, clippy::cast_possible_truncation)]

use omelettehedron::bench;
use omelettehedron::cube::{defs, solved_stickers, MoveSet};
use omelettehedron::light::AeLight;
use omelettehedron::mem::HugeBuf;
use omelettehedron::notation::{
    apply_physical, format_moves, is_solved_physical, parse_moves, scramble_to_cube, Frame,
    PhysMove,
};
use omelettehedron::pattern::{self, PatternError};
use omelettehedron::progress::{Bound, Watch};
use omelettehedron::prune::{checksum, Table2};
use omelettehedron::search::{solve, solve_watched, Limits, Report};
use omelettehedron::tables::{ae_targets, AeIndex, AeSym, CbTable, Phase2, Tables, Tier, N_CB};
use serde_json::json;
use std::sync::{Mutex, MutexGuard, PoisonError};
use std::time::{Duration, Instant};

/// Whether this build has threads (`wasm32-wasip1-threads`, see build.rs).
const THREADS: bool = cfg!(wasm_threads);

/// The cancellation deadline of every solve (`Limits::cancel_after`). The JavaScript wrapper
/// cancels a solve by moving the module's monotonic clock forward (see `cancel` in
/// omelettehedron.mjs) past it; it is long enough never to matter otherwise, and unlike a time
/// limit it does not change the search (an optimal solve keeps the one-thread target-length
/// passes).
const NO_LIMIT: Duration = Duration::from_secs(1_000_000);

/// The WASI reactor initializer, which a Rust cdylib does not get from wasi-libc (its
/// `crt1-reactor`) before Rust 1.100 (`cfg(wasm_own_initialize)`, see build.rs; from 1.100 on
/// the toolchain links wasi-libc's, which runs the constructors and, in its newer wasi-libc,
/// needs no thread-pointer call): the host calls it once, on the main thread, before any other
/// export (`omelettehedron.mjs` does). In the threads build it first sets up the main thread's
/// pthread record (without it, wasi-libc's thread list grew from a zeroed record and threads
/// could hang in their exit).
#[cfg(all(target_family = "wasm", wasm_own_initialize))]
#[unsafe(no_mangle)]
pub extern "C" fn _initialize() {
    static DONE: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);
    unsafe extern "C" {
        #[cfg(wasm_threads)]
        fn __wasi_init_tp();
        fn __wasm_call_ctors();
    }
    if DONE.swap(true, std::sync::atomic::Ordering::Relaxed) {
        return;
    }
    // SAFETY: wasi-libc's and the linker's initializers, run once, first.
    unsafe {
        #[cfg(wasm_threads)]
        __wasi_init_tp();
        __wasm_call_ctors();
    }
}

/// The tier of number `code` (0: small, 1: medium, 2: tiny).
fn tier_of(code: u32) -> Result<Tier, String> {
    match code {
        0 => Ok(Tier::Small),
        1 => Ok(Tier::Medium),
        2 => Ok(Tier::Tiny),
        _ => Err(format!("tables {code}: 0 (small), 1 (medium) or 2 (tiny)")),
    }
}

/// The A-centre/edge index of the small or medium tier.
fn ae_index(tier: Tier) -> AeSym {
    AeSym::new(tier != Tier::Small)
}

/// The tables of a tier: the small or medium tier's A-centre/edge index, or the tiny tier's.
enum AnyTables {
    Sym(Box<Tables<AeSym>>),
    Light(Box<Tables<AeLight>>),
}

/// Evaluate `$body` with `$t` bound to the tables of `$any` (an `&AnyTables`), whatever
/// their index (the search is compiled for each).
macro_rules! with_tables {
    ($any:expr, |$t:ident| $body:expr) => {
        match $any {
            AnyTables::Sym(b) => {
                let $t = &**b;
                $body
            }
            AnyTables::Light(b) => {
                let $t = &**b;
                $body
            }
        }
    };
}

impl AnyTables {
    /// The name of the tier (as the native `--tier`).
    fn tier(&self) -> &'static str {
        match self {
            AnyTables::Sym(t) if t.ae_ix.chiral => "medium",
            AnyTables::Sym(_) => "small",
            AnyTables::Light(_) => "tiny",
        }
    }
}

/// A table file being copied in by the caller: its bytes, in a buffer of whole u64 words (the
/// table's words follow a header of whole words, so the table is used where it lies).
struct Incoming {
    words: HugeBuf<u64>,
    len: usize,
}

/// The module's state between calls.
struct State {
    /// Threads of the global pool (0 before `fto_init`).
    threads: usize,
    tables: Option<AnyTables>,
    /// Table files being copied in by the caller (`fto_table_buffer`): AE, then CB.
    incoming: [Option<Incoming>; 2],
    /// The header of the table file last exported (`fto_export`).
    header: Vec<u8>,
    /// The last call's result (JSON) or error message.
    output: String,
}

static STATE: Mutex<State> = Mutex::new(State {
    threads: 0,
    tables: None,
    incoming: [None, None],
    header: Vec::new(),
    output: String::new(),
});

fn state() -> MutexGuard<'static, State> {
    STATE.lock().unwrap_or_else(PoisonError::into_inner)
}

/// Run `f` on the state; its `Ok` text or `Err` message becomes the output. Returns 0 or -1.
fn with_state(f: impl FnOnce(&mut State) -> Result<String, String>) -> i32 {
    let mut state = state();
    let (output, status) = match f(&mut state) {
        Ok(text) => (text, 0),
        Err(message) => (message, -1),
    };
    state.output = output;
    status
}

/// A string argument.
///
/// # Safety
/// If `len` is nonzero, `ptr` addresses one live allocation containing at least
/// `len <= isize::MAX` bytes, immutable for the returned borrow. Zero length permits a null
/// pointer.
unsafe fn text_argument<'a>(ptr: *const u8, len: usize) -> Result<&'a str, String> {
    if len == 0 {
        return Ok("");
    }
    // SAFETY: as the caller promises.
    let bytes = unsafe { std::slice::from_raw_parts(ptr, len) };
    std::str::from_utf8(bytes).map_err(|_| "the text is not UTF-8".to_string())
}

/// The global pool, built now with `threads` threads (1 without the threads build) if it is
/// not yet. Returns its size.
fn ensure_pool(state: &mut State, threads: usize) -> Result<usize, String> {
    if state.threads == 0 {
        let n = if THREADS { threads.max(1) } else { 1 };
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .use_current_thread()
            .build_global()
            .map_err(|e| format!("cannot start {n} threads: {e}"))?;
        state.threads = n;
    }
    Ok(state.threads)
}

/// The tables, or an error if there are none yet.
fn tables(state: &State) -> Result<&AnyTables, String> {
    state
        .tables
        .as_ref()
        .ok_or_else(|| "no tables: generate or load them first".to_string())
}

/// Allocate `len` bytes for the caller (strings passed in); 0 if out of memory. Free with
/// `fto_free`.
#[unsafe(no_mangle)]
pub extern "C" fn fto_alloc(len: usize) -> *mut u8 {
    let Ok(layout) = std::alloc::Layout::from_size_align(len.max(1), 1) else {
        return std::ptr::null_mut();
    };
    // SAFETY: a non-zero size.
    unsafe { std::alloc::alloc(layout) }
}

/// Free a buffer from `fto_alloc`.
///
/// # Safety
/// `ptr` and `len` are those of a `fto_alloc` buffer not yet freed.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn fto_free(ptr: *mut u8, len: usize) {
    if let Ok(layout) = std::alloc::Layout::from_size_align(len.max(1), 1) {
        // SAFETY: as the caller promises.
        unsafe { std::alloc::dealloc(ptr, layout) };
    }
}

/// The output of the last call: JSON or an error message (UTF-8).
#[unsafe(no_mangle)]
pub extern "C" fn fto_output() -> *const u8 {
    state().output.as_ptr()
}

/// The length in bytes of `fto_output`.
#[unsafe(no_mangle)]
pub extern "C" fn fto_output_len() -> usize {
    state().output.len()
}

/// 1 in the threads build, else 0.
#[unsafe(no_mangle)]
pub extern "C" fn fto_has_threads() -> u32 {
    u32::from(THREADS)
}

/// Start the thread pool with `threads` threads, this one included (the threads build; the
/// other has one). Call it once, before the other functions (which otherwise start one
/// thread); the host must be ready to run `threads - 1` new threads. Output:
/// `{"threads":..}`.
#[unsafe(no_mangle)]
pub extern "C" fn fto_init(threads: u32) -> i32 {
    with_state(|state| {
        let n = ensure_pool(state, threads as usize)?;
        Ok(json!({ "threads": n }).to_string())
    })
}

/// `{"tier":..,"ae_bytes":..,"cb_bytes":..}` of the tables, plus `extra`'s fields.
fn tables_json(t: &AnyTables, extra: &serde_json::Value) -> String {
    let (ae, cb) = with_tables!(t, |t| (t.ae.bytes(), t.cb.tab.bytes()));
    let mut v = json!({
        "tier": t.tier(),
        "ae_bytes": ae,
        "cb_bytes": cb,
    });
    if let (Some(v), Some(extra)) = (v.as_object_mut(), extra.as_object()) {
        v.extend(extra.clone());
    }
    v.to_string()
}

/// Memory that generating needs besides the A-centre/edge table (the move, symmetry and
/// coordinate tables, the corner/B-centre table and the generator's work; measured: the
/// single-thread module grows to 352 MiB for the small tier's 123 MiB table).
const GENERATION_EXTRA: usize = 256 << 20;

/// [`GENERATION_EXTRA`] for the tiny tier, whose generation needs neither the edge symmetry
/// tables nor the 288-symmetry class table (measured: with no reserve, the single-thread module
/// grows to 216 MiB generating its 36 MiB of tables).
const GENERATION_EXTRA_TINY: usize = 128 << 20;

/// An error unless a table of `n` entries and `extra` bytes fit in the module's memory,
/// so that a module without room says so instead of trapping when generating. In the threads
/// build, whose memory is fixed, the sum must not exceed it; and the memory is tried by
/// allocating it (the single-thread build grows its memory for it, and keeps it for the
/// tables). The allocation's address goes through `black_box`: an allocation that is freed
/// unused may be removed by the compiler, success assumed.
fn ensure_room(n: usize, extra: usize, tier: &str) -> Result<(), String> {
    let bytes = n.div_ceil(32) * 8 + extra;
    let err = || {
        format!(
            "not enough memory for the {tier} tables (about {} MiB; this module has {} MiB)",
            bytes >> 20,
            module_memory() >> 20
        )
    };
    if THREADS && module_memory() < bytes as u64 {
        return Err(err());
    }
    let probe = HugeBuf::<u8>::try_for_overwrite(bytes);
    if std::hint::black_box(probe.as_ref().map(|b| b.as_ptr())).is_none() {
        return Err(err());
    }
    Ok(())
}

/// The size of the module's memory in bytes (u64: 4 GiB is 2^32).
fn module_memory() -> u64 {
    #[cfg(target_arch = "wasm32")]
    let pages = core::arch::wasm32::memory_size(0) as u64;
    #[cfg(not(target_arch = "wasm32"))]
    let pages = 0;
    pages << 16
}

/// Generate the tables of `tier` (0: small, 1: medium, 2: tiny) in memory, replacing the current ones
/// (progress lines go to stderr). Output: `{"tier":..,"ae_bytes":..,"cb_bytes":..,
/// "seconds":..}`.
#[unsafe(no_mangle)]
pub extern "C" fn fto_generate(tier: u32) -> i32 {
    with_state(|state| {
        let tier = tier_of(tier)?;
        ensure_pool(state, 1)?;
        state.tables = None;
        state.incoming = [None, None];
        let start = Instant::now();
        let t = if tier == Tier::Tiny {
            let ae_ix = AeLight::new();
            ensure_room(ae_ix.size(), GENERATION_EXTRA_TINY, tier.name())?;
            AnyTables::Light(Box::new(Tables::generate(ae_ix, MoveSet::FACE, true)))
        } else {
            let ae_ix = ae_index(tier);
            ensure_room(ae_ix.size(), GENERATION_EXTRA, tier.name())?;
            AnyTables::Sym(Box::new(Tables::generate(ae_ix, MoveSet::FACE, true)))
        };
        let info = tables_json(&t, &json!({ "seconds": start.elapsed().as_secs_f64() }));
        state.tables = Some(t);
        Ok(info)
    })
}

/// A buffer of `len` bytes for the file of table `which` (0: the A-centre/edge table,
/// `ae-sym288-triples.prun`, `ae-sym288-chiral.prun` or `ae-pair48-split96.prun`; 1:
/// `cb-sym24.prun`), which the caller
/// fills before `fto_load`. Drops the current tables first. Returns 0 if out of memory.
#[unsafe(no_mangle)]
pub extern "C" fn fto_table_buffer(which: u32, len: usize) -> *mut u8 {
    let mut state = state();
    state.tables = None;
    let Some(slot) = state.incoming.get_mut(which as usize) else {
        return std::ptr::null_mut();
    };
    *slot = None;
    let Some(mut words) = HugeBuf::<u64>::try_for_overwrite(len.div_ceil(8).max(1)) else {
        return std::ptr::null_mut();
    };
    // (the caller writes `len` bytes: the rest of the last word is defined too)
    if let Some(last) = words.last_mut() {
        *last = 0;
    }
    let ptr = words.as_mut_ptr().cast::<u8>();
    *slot = Some(Incoming { words, len });
    ptr
}

/// A table from its file, as copied in (`fto_table_buffer`). The format of `Table2::save`
/// (magic, entries, histogram length, histogram, checksum, then the 2-bit entries as
/// little-endian u64 words). `n` is the expected number of entries; the checksum is always
/// checked (a damaged table would give wrong bounds: solutions that are not shortest, or
/// "optimal" claims that are wrong). The header is removed in place: no second copy of the
/// table.
fn parse_table(file: Option<Incoming>, n: usize, name: &str) -> Result<Table2, String> {
    let err = |m: &str| format!("{name}: {m}");
    let Incoming { mut words, len } = file.ok_or_else(|| err("no file given"))?;
    words.iter_mut().for_each(|w| *w = u64::from_le(*w));
    let word = |i: usize| words.get(i).copied().ok_or_else(|| err("truncated header"));
    if word(0)? != u64::from_le_bytes(*b"FTOPRUN2") {
        return Err(err("not a table file (bad magic)"));
    }
    if usize::try_from(word(1)?).ok() != Some(n) {
        return Err(err("size mismatch (a table of another tier or version?)"));
    }
    let nh = usize::try_from(word(2)?)
        .ok()
        .filter(|&h| (1..=64).contains(&h))
        .ok_or_else(|| err("bad histogram length"))?;
    let hist = (0..nh)
        .map(|i| word(3 + i))
        .collect::<Result<Vec<u64>, _>>()?;
    let sum = word(3 + nh)?;
    let data_bytes = 8 * n.div_ceil(32);
    let header = 8 * (4 + nh);
    if len != header + data_bytes {
        return Err(err(&format!(
            "{} data bytes instead of {data_bytes}",
            len.saturating_sub(header),
        )));
    }
    words.drop_front(4 + nh);
    if checksum(&words) != sum {
        return Err(err("checksum mismatch (a damaged file?)"));
    }
    Ok(Table2 {
        data: words,
        n,
        hist,
    })
}

/// The tables of an index, from its A-centre/edge file and the corner/B-centre table.
fn tables_of<I: AeIndex>(
    ae_ix: I,
    ae_file: Option<Incoming>,
    cb: CbTable,
) -> Result<Tables<I>, String> {
    let ae_name = format!("{}.prun", ae_ix.name());
    let ae = parse_table(ae_file, ae_ix.size(), &ae_name)?;
    Ok(Tables {
        ae_targets: ae_targets(&ae_ix),
        ae_ix,
        ae,
        cb,
        p2: Phase2::new(),
        cbe: None,
        cbs: None,
        moves: MoveSet::FACE,
    })
}

/// The entry count in the header of a table file being copied in, if it has one.
fn entries(file: Option<&Incoming>) -> Option<usize> {
    file.and_then(|f| f.words.get(1))
        .and_then(|&w| usize::try_from(u64::from_le(w)).ok())
}

/// Install the table files copied into the `fto_table_buffer`s (both must be filled), after
/// checking their format, size and checksums; the tier (small, medium or tiny) is that of the
/// A-centre/edge table's size. Output: `{"tier":..,"ae_bytes":..,"cb_bytes":..,"seconds":..}`.
#[unsafe(no_mangle)]
pub extern "C" fn fto_load() -> i32 {
    with_state(|state| {
        ensure_pool(state, 1)?;
        state.tables = None;
        let start = Instant::now();
        let [ae_file, cb_file] = std::mem::take(&mut state.incoming);
        let n = entries(ae_file.as_ref());
        let cb = parse_table(cb_file, N_CB, "cb-sym24.prun")?;
        let cb = CbTable {
            tab: cb,
            attach: CbTable::attach_targets(),
        };
        let t = if n == usize::try_from(omelettehedron::tables::AE_LIGHT_ENTRIES).ok() {
            AnyTables::Light(Box::new(tables_of(AeLight::new(), ae_file, cb)?))
        } else {
            let ae_ix = [Tier::Small, Tier::Medium]
                .into_iter()
                .map(ae_index)
                .find(|ix| Some(ix.size()) == n)
                .ok_or(
                    "the A-centre/edge table file is none of ae-sym288-triples.prun, \
                     ae-sym288-chiral.prun and ae-pair48-split96.prun of this version",
                )?;
            AnyTables::Sym(Box::new(tables_of(ae_ix, ae_file, cb)?))
        };
        let info = tables_json(&t, &json!({ "seconds": start.elapsed().as_secs_f64() }));
        state.tables = Some(t);
        Ok(info)
    })
}

/// The file of table `which` (0: AE, 1: CB) in the format of `fto_load` and the native `gen`,
/// as two pieces for the caller to copy out: the header (in a module buffer, valid until the
/// next `fto_export`) and the table's data (valid while the tables are). Output:
/// `{"name":..,"header_ptr":..,"header_len":..,"data_ptr":..,"data_len":..}`.
#[unsafe(no_mangle)]
pub extern "C" fn fto_export(which: u32) -> i32 {
    with_state(|state| {
        let t = tables(state)?;
        let (name, tab) = with_tables!(t, |t| match which {
            0 => (format!("{}.prun", t.ae_ix.name()), &t.ae),
            1 => ("cb-sym24.prun".to_string(), &t.cb.tab),
            _ => return Err(format!("table {which}: 0 (AE) or 1 (CB)")),
        });
        // (as `prune::save_words` writes it)
        let mut h = Vec::with_capacity(32 + 8 * tab.hist.len());
        h.extend_from_slice(b"FTOPRUN2");
        h.extend_from_slice(&(tab.n as u64).to_le_bytes());
        h.extend_from_slice(&(tab.hist.len() as u64).to_le_bytes());
        for x in &tab.hist {
            h.extend_from_slice(&x.to_le_bytes());
        }
        h.extend_from_slice(&checksum(&tab.data).to_le_bytes());
        let (data_ptr, data_len) = (tab.data.as_ptr() as usize, tab.data.len() * 8);
        state.header = h;
        Ok(json!({
            "name": name,
            "header_ptr": state.header.as_ptr() as usize,
            "header_len": state.header.len(),
            "data_ptr": data_ptr,
            "data_len": data_len,
        })
        .to_string())
    })
}

/// True if `solution` solves the position that `scramble` makes from solved, checked on the
/// sticker model (up to a whole-puzzle rotation).
fn solves(scramble: &[PhysMove], solution: &[PhysMove]) -> bool {
    let scrambled = apply_physical(&solved_stickers(defs()), scramble);
    is_solved_physical(&apply_physical(&scrambled, solution))
}

/// The limits of a page solve. `seconds` > 0: a time limit (a stopping rule); 0: optimal mode;
/// < 0: the improving search with no time limit (`Limits::improve`: each better solution as it
/// is found, until the last is proven optimal). Every solve can be cancelled (by moving the
/// clock past `NO_LIMIT`) without changing the search.
fn page_limits(seconds: f64, max_length: u32) -> Result<Limits, String> {
    let time = if seconds > 0.0 {
        Some(Duration::try_from_secs_f64(seconds).map_err(|e| format!("time limit: {e}"))?)
    } else {
        None
    };
    Ok(Limits {
        time,
        max_length: (max_length > 0).then_some(max_length as usize),
        cancel_after: Some(NO_LIMIT),
        improve: seconds < 0.0,
        ..Limits::default()
    })
}

/// A page solve. With `live`, stdout gets each better solution (`<seconds>s  <length>  <moves>`,
/// printed by the search) and each larger proven lower bound (`<seconds>s  lower <n>`: no
/// solution has fewer moves), so that the page can show the search's progress.
fn page_solve<I: AeIndex>(
    t: &Tables<I>,
    cube: omelettehedron::cube::Cube,
    frame: Frame,
    limits: &Limits,
    live: bool,
) -> omelettehedron::search::Outcome {
    let start = Instant::now();
    let lower = |b: Bound| {
        if let Bound::Lower(n) = b {
            println!("{:>9.3}s  lower {n}", start.elapsed().as_secs_f64());
        }
    };
    let watch = Watch {
        on_bound: live.then_some(&lower as &(dyn Fn(Bound) + Sync)),
        ..Watch::default()
    };
    let report = if live {
        Report::Improvements
    } else {
        Report::Quiet
    };
    solve_watched(t, cube, frame, limits, report, &watch)
}

/// Solve a scramble (FTO notation, e.g. `U R' BL F`). `seconds` > 0: stop after that time;
/// 0: run until the solution is proven optimal. `max_length` > 0: stop at the first solution
/// of at most that many moves. `improvements` != 0: print each better solution to stdout as
/// it is found (`<seconds>s  <length>  <moves>`). Output: `{"solution":..,"length":..,
/// "optimal":..,"verified":..,"seconds":..,"nodes":..,"solutions":..,"events":[[s,len],..]}`
/// (`solution` and `length` null if none was found in time; `verified`: the solution solves
/// the scramble on the sticker model).
///
/// # Safety
/// If `len` is nonzero, `text` points into one live allocation containing `len <= isize::MAX`
/// bytes, immutable until this call returns.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn fto_solve(
    text: *const u8,
    len: usize,
    seconds: f64,
    max_length: u32,
    improvements: u32,
) -> i32 {
    // SAFETY: as the caller promises.
    let text = unsafe { text_argument(text, len) };
    with_state(|state| {
        let moves = parse_moves(text?)?;
        ensure_pool(state, 1)?;
        let t = tables(state)?;
        let limits = page_limits(seconds, max_length)?;
        let (cube, frame) = scramble_to_cube(&moves);
        let o = with_tables!(t, |t| page_solve(
            t,
            cube,
            frame,
            &limits,
            improvements != 0
        ));
        let events: Vec<(f64, usize)> = o.events.iter().map(|e| (e.time, e.length)).collect();
        Ok(json!({
            "solution": o.solution.as_deref().map(format_moves),
            "length": o.length,
            "optimal": o.optimal,
            "verified": o.solution.as_deref().is_some_and(|s| solves(&moves, s)),
            "seconds": o.elapsed,
            "nodes": o.nodes,
            "solutions": o.solutions_found,
            "events": events,
        })
        .to_string())
    })
}

/// A random-state scramble: a uniformly random position (`seed`), solved up to the first
/// solution, and that solution inverted. (Any solution would do; the first one takes from
/// milliseconds to seconds, and has 22-25 moves.) Output: `{"scramble":..,"length":..}`.
#[unsafe(no_mangle)]
pub extern "C" fn fto_random_scramble(seed: u32) -> i32 {
    with_state(|state| {
        ensure_pool(state, 1)?;
        let t = tables(state)?;
        let cube = bench::random_positions(1, u64::from(seed))
            .into_iter()
            .next()
            .ok_or("no random position")?;
        let limits = Limits {
            max_solutions: Some(1),
            cancel_after: Some(NO_LIMIT),
            ..Limits::default()
        };
        let o = with_tables!(t, |t| solve(
            t,
            cube,
            Frame::identity(),
            &limits,
            Report::Quiet
        ));
        let s = o.solution.ok_or("cancelled")?;
        let scramble: Vec<PhysMove> = s
            .iter()
            .rev()
            .map(|m| PhysMove {
                face: m.face,
                layer: m.layer,
                cw: !m.cw,
            })
            .collect();
        Ok(json!({ "scramble": format_moves(&scramble), "length": scramble.len() }).to_string())
    })
}

/// Whether a move sequence, applied to the solved puzzle, leaves it solved (up to a
/// whole-puzzle rotation), on the sticker model; e.g. a scramble followed by its solution.
/// Output: `{"solved":..,"moves":..}`.
///
/// # Safety
/// As `fto_solve`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn fto_is_solved(text: *const u8, len: usize) -> i32 {
    // SAFETY: as the caller promises.
    let text = unsafe { text_argument(text, len) };
    with_state(|_| {
        let moves = parse_moves(text?)?;
        Ok(json!({ "solved": solves(&[], &moves), "moves": moves.len() }).to_string())
    })
}

/// The native `bench`: `count` seeded random positions (or, with `scramble_len` > 0,
/// random-move scrambles of that length) solved one at a time with a time limit of `seconds`
/// (0: optimal). Output: the native `--json` summary (tier, threads, checkpoints,
/// `mean_length`, `nodes_per_sec`, `results`, ...).
#[unsafe(no_mangle)]
pub extern "C" fn fto_bench(count: u32, seed: u32, seconds: f64, scramble_len: u32) -> i32 {
    with_state(|state| {
        ensure_pool(state, 1)?;
        let t = tables(state)?;
        let time = if seconds > 0.0 {
            Some(Duration::try_from_secs_f64(seconds).map_err(|e| format!("time limit: {e}"))?)
        } else {
            None
        };
        let limits = Limits {
            time,
            ..Limits::default()
        };
        let s = with_tables!(t, |tt| bench::run(
            tt,
            t.tier(),
            count as usize,
            u64::from(seed),
            (scramble_len > 0).then_some(scramble_len as usize),
            &limits,
            false,
        ));
        serde_json::to_string(&s).map_err(|e| e.to_string())
    })
}

/// The pattern a scramble makes (as the native `pattern`): the 72 sticker colours in picture
/// order (src/pattern.rs), with the letters of the solved faces, W R G P / B O Y A for U R F L
/// / B BL D BR. Output: `{"pattern":..}`.
///
/// # Safety
/// As `fto_solve`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn fto_pattern(text: *const u8, len: usize) -> i32 {
    // SAFETY: as the caller promises.
    let text = unsafe { text_argument(text, len) };
    with_state(|_| {
        let moves = parse_moves(text?)?;
        Ok(json!({ "pattern": pattern::format(&pattern::of_moves(&moves)) }).to_string())
    })
}

/// What is wrong with a pattern, for the page: `{"message":.., "kind":.., "stickers": the
/// picture positions to point out, "counts": [[colour, stickers], ..] for the colours not on
/// 9 stickers}`.
fn pattern_error(e: &PatternError) -> serde_json::Value {
    let kind = match e {
        PatternError::Length(_) => "length",
        PatternError::Colours(_) => "colours",
        PatternError::Counts(_) => "counts",
        PatternError::Pieces(_) => "pieces",
        PatternError::Duplicates(_) => "duplicates",
        PatternError::Unreachable(_) => "unreachable",
    };
    let counts: Vec<(String, usize)> = match e {
        PatternError::Counts(c) => c.iter().map(|&(s, n)| (s.to_string(), n)).collect(),
        _ => vec![],
    };
    json!({
        "message": e.to_string(),
        "kind": kind,
        "stickers": e.stickers(),
        "counts": counts,
    })
}

/// Solve a pattern (src/pattern.rs) as `fto_solve` solves a scramble. `usual`: the colours of
/// U R F L B BL D BR on the solved puzzle as usually held (8 letters of the pattern), or empty.
/// Output: as `fto_solve` (`verified`: the solution leaves every face one colour), plus
/// `"scramble": {"start": the colours of U R F L B BL D BR on the solved puzzle it starts
/// from, "moves":.., "usual": the scramble from the usual hold or null, "exact": whether that
/// makes the pattern as entered (otherwise the same puzzle, held another way)}`. A wrong
/// pattern gives `{"error": ..}` (see `pattern_error`) instead of a solve.
///
/// # Safety
/// As `fto_solve`, for both strings.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn fto_solve_pattern(
    text: *const u8,
    len: usize,
    usual: *const u8,
    usual_len: usize,
    seconds: f64,
    max_length: u32,
    improvements: u32,
) -> i32 {
    // SAFETY: as the caller promises.
    let (text, usual) = unsafe { (text_argument(text, len), text_argument(usual, usual_len)) };
    with_state(|state| {
        let p = match pattern::parse(text?).and_then(|p| pattern::position(&p).map(|c| (p, c))) {
            Ok(x) => x,
            Err(e) => return Ok(json!({ "error": pattern_error(&e) }).to_string()),
        };
        let (p, cube) = p;
        let usual: Option<[char; 8]> = match usual? {
            "" => None,
            u => Some(
                u.chars()
                    .collect::<Vec<_>>()
                    .try_into()
                    .map_err(|_| "usual: 8 colours".to_string())?,
            ),
        };
        ensure_pool(state, 1)?;
        let t = tables(state)?;
        let limits = page_limits(seconds, max_length)?;
        let o = with_tables!(t, |t| page_solve(
            t,
            cube,
            Frame::identity(),
            &limits,
            improvements != 0
        ));
        let events: Vec<(f64, usize)> = o.events.iter().map(|e| (e.time, e.length)).collect();
        let scramble = o
            .solution
            .as_deref()
            .and_then(|s| pattern::scramble(&p, s, usual.as_ref()))
            .map(|sc| {
                let from_usual = sc.from_usual_text();
                json!({
                    "start": sc.start.iter().collect::<String>(),
                    "moves": format_moves(&sc.moves),
                    "usual": from_usual.as_ref().map(|(t, _)| t.clone()),
                    "exact": from_usual.is_some_and(|(_, exact)| exact),
                })
            });
        Ok(json!({
            "solution": o.solution.as_deref().map(format_moves),
            "length": o.length,
            "optimal": o.optimal,
            "verified": o.solution.as_deref().is_some_and(|s| pattern::solves(&p, s)),
            "seconds": o.elapsed,
            "nodes": o.nodes,
            "solutions": o.solutions_found,
            "events": events,
            "scramble": scramble,
        })
        .to_string())
    })
}
