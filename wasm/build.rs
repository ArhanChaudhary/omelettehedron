//! For the threads build (wasm32-wasip1-threads): `cfg(wasm_threads)` (its `atomics` target
//! feature is not visible to `cfg(target_feature)`), and all of the module's memory from the
//! start, so that it never grows.
//!
//! The reason: memory that grows while threads run is not safe in V8. A bulk memory operation (`memory.fill`,
//! `memory.copy`) checks its bounds against the size the thread's instance last saw, which a
//! grow updates only later, so a `memory.fill` of a block fresh from `malloc` could trap
//! ("memory access out of bounds"). With the initial memory equal to the maximum, wasi-libc's
//! `malloc` has it all without growing; the engine maps it, and untouched pages cost no RAM.
//!
//! The small tier needs ~300 MB (the 129 MB A-centre/edge table, the 4.7 MB corner/B-centre
//! table, the move and symmetry tables, and while loading a table file its bytes too), and
//! each thread a 2 MiB stack: 1 GiB leaves room for 256 threads.
//!
//! The medium tier needs ~1.3 GB (its A-centre/edge table is 1.03 GB; a table file is loaded
//! where it is copied in, without a second copy). It gets its own module, built with the
//! feature `medium` (`omelettehedron-threads-medium.wasm`, see scripts/build-wasm.sh), with
//! 2 GiB: that leaves room for over 300 threads, and the small tier's module keeps its 1 GiB,
//! which more devices can reserve. (The memory cannot simply be larger than the module says:
//! `malloc` uses only the memory up to the module's initial size, `__heap_end`, and then
//! grows it.)

/// Memory of the threads build, in bytes (a multiple of the 64 KiB page).
const MEMORY: u64 = 1 << 30;

/// Memory of the threads build with the feature `medium`.
const MEDIUM_MEMORY: u64 = 2 << 30;

fn main() {
    println!("cargo::rustc-check-cfg=cfg(wasm_threads)");
    println!("cargo::rustc-check-cfg=cfg(wasm_own_initialize)");
    println!("cargo::rerun-if-changed=build.rs");
    if rustc_minor().is_some_and(|minor| minor < 100) {
        println!("cargo::rustc-cfg=wasm_own_initialize");
    }
    if std::env::var("TARGET").is_ok_and(|target| target.ends_with("-threads")) {
        let memory = if std::env::var_os("CARGO_FEATURE_MEDIUM").is_some() {
            MEDIUM_MEMORY
        } else {
            MEMORY
        };
        println!("cargo::rustc-cfg=wasm_threads");
        println!("cargo::rustc-cdylib-link-arg=--initial-memory={memory}");
        println!("cargo::rustc-cdylib-link-arg=--max-memory={memory}");
    }
}

/// The minor version of the compiler (`rustc 1.100.0 (...)` -> 100), if it can be read.
///
/// Up to 1.99 a Rust cdylib for WASI does not link wasi-libc's `crt1-reactor.o`, so the module
/// defines the reactor initializer `_initialize` itself (`cfg(wasm_own_initialize)`); from 1.100
/// on it does, and a second `_initialize` is a duplicate symbol at link time.
fn rustc_minor() -> Option<u32> {
    let rustc = std::env::var_os("RUSTC")?;
    let out = std::process::Command::new(rustc)
        .arg("--version")
        .output()
        .ok()?;
    let text = String::from_utf8(out.stdout).ok()?;
    text.split_whitespace()
        .nth(1)?
        .split('.')
        .nth(1)?
        .parse()
        .ok()
}
