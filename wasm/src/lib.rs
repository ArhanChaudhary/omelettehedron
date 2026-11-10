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
//! Two builds (scripts/build-wasm.sh): `wasm32-wasip1` runs on one thread anywhere, and
//! `wasm32-wasip1-threads` runs rayon's pool on wasi-threads (in a browser: a cross-origin
//! isolated page). In both, the pool is rayon's global one, built once by `fto_init` with the
//! calling thread as one of its workers. It is never dropped, so its threads never exit: a
//! thread that exits can race with wasi-libc's thread exit (egg-solver, wasm/src/threads.rs).
//!
//! Only the small tier: the larger tiers' indices and tables exceed wasm32's address space.
//!
//! Adapted from egg-solver's wasm/src/lib.rs (same author): the ABI conventions (`fto_alloc`,
//! `fto_free`, the output buffer, `with_state`, `text_argument`), the reactor initializer
//! `_initialize` and the threads build's set-up.

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
use omelettehedron::mem::HugeBuf;
use omelettehedron::notation::{
    apply_physical, format_moves, is_solved_physical, parse_moves, scramble_to_cube, Frame,
    PhysMove,
};
use omelettehedron::prune::{checksum, Table2};
use omelettehedron::search::{solve, Limits, Report};
use omelettehedron::tables::{ae_targets, AeIndex, AeSym, CbTable, Phase2, Tables, N_CB};
use serde_json::json;
use std::sync::{Mutex, MutexGuard, PoisonError};
use std::time::{Duration, Instant};

/// Whether this build has threads (`wasm32-wasip1-threads`, see build.rs).
const THREADS: bool = cfg!(wasm_threads);

/// The time limit of an optimal solve. The JavaScript wrapper cancels a solve by moving the
/// module's monotonic clock forward (see `cancel` in omelettehedron.mjs), which only stops a
/// search that has a time limit; this one is long enough never to matter otherwise. (The time
/// limit has no other effect on the search: the optimality proof is the same without one.)
const NO_LIMIT: Duration = Duration::from_secs(1_000_000);

/// The WASI reactor initializer, which a Rust cdylib does not get from wasi-libc (its
/// `crt1-reactor`): the host calls it once, on the main thread, before any other export
/// (`omelettehedron.mjs` does). In the threads build it first sets up the main thread's
/// pthread record (without it, wasi-libc's thread list grew from a zeroed record and threads
/// could hang in their exit; egg-solver). Adapted from egg-solver's wasm/src/lib.rs.
#[cfg(target_family = "wasm")]
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

/// The module's state between calls.
struct State {
    /// Threads of the global pool (0 before `fto_init`).
    threads: usize,
    tables: Option<Tables<AeSym>>,
    /// Table files being copied in by the caller (`fto_table_buffer`): AE, then CB.
    incoming: [Vec<u8>; 2],
    /// The header of the table file last exported (`fto_export`).
    header: Vec<u8>,
    /// The last call's result (JSON) or error message.
    output: String,
}

static STATE: Mutex<State> = Mutex::new(State {
    threads: 0,
    tables: None,
    incoming: [Vec::new(), Vec::new()],
    header: Vec::new(),
    output: String::new(),
});

fn state() -> MutexGuard<'static, State> {
    STATE.lock().unwrap_or_else(PoisonError::into_inner)
}

/// Run `f` on the state; its `Ok` text or `Err` message becomes the output. Returns 0 or -1.
/// (From egg-solver's wasm/src/lib.rs.)
fn with_state(f: impl FnOnce(&mut State) -> Result<String, String>) -> i32 {
    let mut state = state();
    let (output, status) = match f(&mut state) {
        Ok(text) => (text, 0),
        Err(message) => (message, -1),
    };
    state.output = output;
    status
}

/// A string argument. (From egg-solver's wasm/src/lib.rs.)
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
fn tables(state: &State) -> Result<&Tables<AeSym>, String> {
    state
        .tables
        .as_ref()
        .ok_or_else(|| "no tables: generate or load them first".to_string())
}

/// Allocate `len` bytes for the caller (strings passed in); 0 if out of memory. Free with
/// `fto_free`. (From egg-solver's wasm/src/lib.rs.)
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

/// `{"ae_bytes":..,"cb_bytes":..}` of the tables, plus `extra`'s fields.
fn tables_json(t: &Tables<AeSym>, extra: &serde_json::Value) -> String {
    let mut v = json!({ "ae_bytes": t.ae.bytes(), "cb_bytes": t.cb.tab.bytes() });
    if let (Some(v), Some(extra)) = (v.as_object_mut(), extra.as_object()) {
        v.extend(extra.clone());
    }
    v.to_string()
}

/// Generate the small tier's tables in memory, replacing the current ones (progress lines go
/// to stderr). Output: `{"ae_bytes":..,"cb_bytes":..,"seconds":..}`.
#[unsafe(no_mangle)]
pub extern "C" fn fto_generate() -> i32 {
    with_state(|state| {
        ensure_pool(state, 1)?;
        state.tables = None;
        state.incoming = [Vec::new(), Vec::new()];
        let start = Instant::now();
        let t = Tables::generate(AeSym::new(false), MoveSet::FACE, true);
        let info = tables_json(&t, &json!({ "seconds": start.elapsed().as_secs_f64() }));
        state.tables = Some(t);
        Ok(info)
    })
}

/// Name and entry count of table `which` (0: AE, 1: CB) of the small tier.
fn table_kind(which: u32) -> Result<(&'static str, usize), String> {
    match which {
        0 => Ok(("ae-sym288-triples.prun", AeSym::new(false).size())),
        1 => Ok(("cb-sym24.prun", N_CB)),
        _ => Err(format!("table {which}: 0 (AE) or 1 (CB)")),
    }
}

/// A buffer of `len` bytes for the file of table `which` (0: `ae-sym288-triples.prun`, 1:
/// `cb-sym24.prun`), which the caller fills before `fto_load`. Drops the current tables
/// first. Returns 0 if out of memory.
#[unsafe(no_mangle)]
pub extern "C" fn fto_table_buffer(which: u32, len: usize) -> *mut u8 {
    let mut state = state();
    state.tables = None;
    let Some(slot) = state.incoming.get_mut(which as usize) else {
        return std::ptr::null_mut();
    };
    *slot = Vec::new();
    if slot.try_reserve_exact(len.max(1)).is_err() {
        return std::ptr::null_mut();
    }
    slot.resize(len, 0);
    slot.as_mut_ptr()
}

/// Little-endian u64 `i` of `bytes`, if there is one.
fn word(bytes: &[u8], i: usize) -> Option<u64> {
    let b = bytes.get(8 * i..8 * i + 8)?;
    Some(u64::from_le_bytes(b.try_into().ok()?))
}

/// A table of `n` entries from the bytes of its file. The format of `Table2::save` (magic,
/// entries, histogram length, histogram, checksum, then the 2-bit entries as little-endian
/// u64 words); the checksum is always checked (a damaged table would give wrong bounds:
/// solutions that are not shortest, or "optimal" claims that are wrong).
fn parse_table(bytes: &[u8], n: usize, name: &str) -> Result<Table2, String> {
    let err = |m: &str| format!("{name}: {m}");
    if bytes.get(..8) != Some(b"FTOPRUN2".as_slice()) {
        return Err(err("not a table file (bad magic)"));
    }
    if word(bytes, 1).and_then(|x| usize::try_from(x).ok()) != Some(n) {
        return Err(err("size mismatch (a table of another tier or version?)"));
    }
    let nh = word(bytes, 2)
        .and_then(|x| usize::try_from(x).ok())
        .filter(|&h| (1..=64).contains(&h))
        .ok_or_else(|| err("bad histogram length"))?;
    let hist = (0..nh)
        .map(|i| word(bytes, 3 + i))
        .collect::<Option<Vec<u64>>>()
        .ok_or_else(|| err("truncated header"))?;
    let sum = word(bytes, 3 + nh).ok_or_else(|| err("truncated header"))?;
    let words = n.div_ceil(32);
    let body = bytes.get(8 * (4 + nh)..).unwrap_or_default();
    if body.len() != words * 8 {
        return Err(err(&format!(
            "{} data bytes instead of {}",
            body.len(),
            words * 8
        )));
    }
    let mut data = HugeBuf::<u64>::for_overwrite(words);
    for (w, b) in data.iter_mut().zip(body.as_chunks::<8>().0) {
        *w = u64::from_le_bytes(*b);
    }
    if checksum(&data) != sum {
        return Err(err("checksum mismatch (a damaged file?)"));
    }
    Ok(Table2 { data, n, hist })
}

/// Install the table files copied into the `fto_table_buffer`s (both must be filled), after
/// checking their format, size and checksums. Output: `{"ae_bytes":..,"cb_bytes":..,
/// "seconds":..}`.
#[unsafe(no_mangle)]
pub extern "C" fn fto_load() -> i32 {
    with_state(|state| {
        ensure_pool(state, 1)?;
        state.tables = None;
        let start = Instant::now();
        let [ae_file, cb_file] = std::mem::take(&mut state.incoming);
        let ae_ix = AeSym::new(false);
        let (ae_name, ae_n) = table_kind(0)?;
        let (cb_name, cb_n) = table_kind(1)?;
        // (each file is freed once parsed: less memory at the peak)
        let cb = parse_table(&cb_file, cb_n, cb_name)?;
        drop(cb_file);
        let ae = parse_table(&ae_file, ae_n, ae_name)?;
        drop(ae_file);
        let t = Tables {
            ae_targets: ae_targets(&ae_ix),
            ae_ix,
            ae,
            cb: CbTable {
                tab: cb,
                attach: CbTable::attach_targets(),
            },
            p2: Phase2::new(),
            cbe: None,
            moves: MoveSet::FACE,
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
        let (name, _) = table_kind(which)?;
        let t = tables(state)?;
        let tab = if which == 0 { &t.ae } else { &t.cb.tab };
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
        let time = if seconds > 0.0 {
            Duration::try_from_secs_f64(seconds).map_err(|e| format!("time limit: {e}"))?
        } else {
            NO_LIMIT
        };
        let limits = Limits {
            time: Some(time),
            max_length: (max_length > 0).then_some(max_length as usize),
            ..Limits::default()
        };
        let (cube, frame) = scramble_to_cube(&moves);
        let report = if improvements == 0 {
            Report::Quiet
        } else {
            Report::Improvements
        };
        let o = solve(t, cube, frame, &limits, report);
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
            time: Some(NO_LIMIT),
            ..Limits::default()
        };
        let o = solve(t, cube, Frame::identity(), &limits, Report::Quiet);
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
        let s = bench::run(
            t,
            "small",
            count as usize,
            u64::from(seed),
            (scramble_len > 0).then_some(scramble_len as usize),
            &limits,
            false,
        );
        serde_json::to_string(&s).map_err(|e| e.to_string())
    })
}
