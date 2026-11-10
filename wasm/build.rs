//! For the threads build (wasm32-wasip1-threads): `cfg(wasm_threads)` (its `atomics` target
//! feature is not visible to `cfg(target_feature)`), and all of the module's memory from the
//! start, so that it never grows.
//!
//! Adapted from egg-solver's wasm/build.rs (same author), which found why: memory that grows
//! while threads run is not safe in V8. A bulk memory operation (`memory.fill`,
//! `memory.copy`) checks its bounds against the size the thread's instance last saw, which a
//! grow updates only later, so a `memory.fill` of a block fresh from `malloc` could trap
//! ("memory access out of bounds"). With the initial memory equal to the maximum, wasi-libc's
//! `malloc` has it all without growing; the engine maps it, and untouched pages cost no RAM.
//!
//! The small tier needs ~300 MB (the 129 MB A-centre/edge table, the 4.7 MB corner/B-centre
//! table, the move and symmetry tables, and while loading a table file its bytes too), and
//! each thread a 2 MiB stack: 1 GiB leaves room for 256 threads.

/// Memory of the threads build, in bytes (a multiple of the 64 KiB page).
const MEMORY: u64 = 1 << 30;

fn main() {
    println!("cargo::rustc-check-cfg=cfg(wasm_threads)");
    println!("cargo::rerun-if-changed=build.rs");
    if std::env::var("TARGET").is_ok_and(|target| target.ends_with("-threads")) {
        println!("cargo::rustc-cfg=wasm_threads");
        println!("cargo::rustc-cdylib-link-arg=--initial-memory={MEMORY}");
        println!("cargo::rustc-cdylib-link-arg=--max-memory={MEMORY}");
    }
}
