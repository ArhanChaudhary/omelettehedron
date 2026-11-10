# omelettehedron 🥚

A fast two-phase solver for the face-turning octahedron (FTO), written in Rust.

* **Phase 1** brings the puzzle into H_A, the positions that can be solved by turning
  only the four class-A faces U, F, BR and BL.
* **Phase 2** finishes with those faces only. Its 933,120 positions are fully tabulated, so it
  is an O(length) lookup, and its solutions are optimal even among all-move solutions.

Phase 1 is an IDA* search guided by pruning tables. See [docs/PHASES.md](docs/PHASES.md)
for why this split was chosen and how it compares with alternatives,
[docs/RESULTS.md](docs/RESULTS.md) for measured speed and solution lengths, and
[docs/OPTIMIZATION_LOG.md](docs/OPTIMIZATION_LOG.md) for how it got there.

On 48 threads a random position gets its first solution in about 6 ms (23 moves on average),
21.5 moves after 1 s and 21.1 after 5 s; a proof of optimality takes from a few seconds to about
a minute (optimal lengths are around 20-21).

## Build

```sh
cargo build --release
```

On x86-64 the build targets `x86-64-v3` (AVX2/BMI2, any x86-64 CPU from about 2015 on); see
`.cargo/config.toml`. Other targets build too: aarch64 (e.g. Apple Silicon) uses NEON for
the byte shuffles, wasm32 uses simd128 when enabled (`-C target-feature=+simd128`), and
anything else falls back to scalar code.

## Notation

Hold the puzzle with a corner pointing at you:

| Side | Faces |
|---|---|
| Front | U (top), R (right), F (bottom), L (left) |
| Back | B (top), BR, D (bottom), BL |

`X` turns face X 120° clockwise as seen from outside, and `X'` turns it anticlockwise. `X2` is
accepted as two turns. Opposite pairs: U–D, F–B, R–BL, L–BR.

Face turns can rotate the whole puzzle, so a solution ends in the solved state in *some*
orientation. Internally the front corner is held fixed.

## Usage

```sh
# optimal solution (default when no stopping rule is given)
target/release/omelettehedron solve --tier large "U R' BL F D B' L BR"

# stop after 2 seconds, printing every improvement as it is found
target/release/omelettehedron solve --tier large -t 2 -p improvements "R U' L F' BR B' ..."

# stop at the first solution of at most 22 moves
target/release/omelettehedron solve --tier large -l 22 "..."

# stop after 1000 solutions (every phase-2 completion counts, even if longer than the best)
target/release/omelettehedron solve --tier large -n 1000 "..."

# one scramble per line on stdin
target/release/omelettehedron solve --tier large -t 1 < scrambles.txt
```

Stopping rules can be combined; the search stops when any of them triggers. `-j N` sets the
thread count (default: all cores).

### Pruning tables

`--tier` is required (there is no default size): choose the largest tier that fits
comfortably in RAM. Tables are generated on first use. They are written to `./tables`, or to `--tables DIR` or
`$FTO_TABLES` if given. `omelettehedron gen --tier T` generates them ahead of time. Before generating a
table, the solver checks that the disk and the free memory can hold it, and stops at once
with a message if not.

| `--tier` | Memory | Generation (48 cores) | Contents |
|---|---|---|---|
| `small` | ~260 MB | ~1 s | A-centre/edge table without edge chirality |
| `medium` | ~1.2 GB | ~7 s | full A-centre/edge table (288-fold symmetry) |
| `large` | ~8 GB | ~30 s | medium + joint corner/B-centre/edge-triple table (about 2x fewer nodes) |
| `huge` | ~31 GB | ~3 min | relabel-only A-centre/edge table (cheaper lookups) + the joint table |

With `FTO_HUGE1G=1`, tables of 256 MB or more are placed on 1 GiB pages when the system has
reserved them (`echo N > /sys/kernel/mm/hugepages/hugepages-1048576kB/nr_hugepages` as root; the
large tier needs 9): optimal proofs about 4% faster with the A-centre/edge table on a 1 GiB page.
Without reserved pages the tables use transparent huge pages.

All tiers also use the 89 MB corner/B-centre table and the in-memory phase-2 table.
Generation needs about the memory of the tier (peak 8.1 GB for `large`, 29 GB for `huge`); on
one core it takes ~30 s for `small` and ~3 min for `medium`. Writing the 30 GB of `huge` to
disk is a large part of its time.

## Benchmarks

```sh
target/release/omelettehedron bench --tier large -N 40 -t 5       # 40 seeded random positions, 5 s each
target/release/omelettehedron bench --tier large -N 20 -t 5 -j 1  # single-threaded
scripts/bench.sh <label> [tier]                 # the standard set used during optimisation
scripts/bench-opt.sh <label> [tier]             # optimal mode: time to prove optimality
scripts/bench-randopt.sh <label> <tier>         # optimal mode on random positions (slow)
scripts/compare.py <before> <after>             # markdown comparison of two runs
```

`bench` reports these figures:
* mean best length after 0.01 s, 0.1 s, 1 s, … and at the time limit;
* mean time to the first solution and its length;
* how many positions were proven optimal;
* search speed in nodes per second.

## WebAssembly demo

`wasm/` holds a browser demo: the solver as a WebAssembly module, a JavaScript wrapper and a
static page. It supports the small tier only (wasm32 has a 4 GB address space).

```sh
rustup target add wasm32-wasip1 wasm32-wasip1-threads
scripts/build-wasm.sh        # wasm/dist/omelettehedron.wasm and omelettehedron-threads.wasm
bun wasm/serve.mjs           # then open http://localhost:8000/
bun wasm/test.mjs            # headless test in Bun (FTO_TABLES=DIR also compares the tables)
bun wasm/bench.mjs [--threads N] [--tables DIR] [--count 10] [--seconds 5]
```

* `omelettehedron.wasm` (`wasm32-wasip1`) runs on one thread anywhere.
  `omelettehedron-threads.wasm` (`wasm32-wasip1-threads`) runs rayon's thread pool on
  wasi-threads. It needs `SharedArrayBuffer`, so a browser page must be cross-origin isolated
  (the COOP/COEP headers that `wasm/serve.mjs` sends; GitHub Pages cannot send them). Its
  memory is fixed at 1 GiB from the start: growing the shared memory while threads run is
  unsafe in V8 (see `wasm/build.rs`).
* `wasm/omelettehedron.mjs` provides the WASI functions the module imports and a worker pool
  for its threads. `FtoSolver` generates or loads the tables, solves scrambles optimally or
  with a time limit, makes random-state scrambles and checks solutions on the sticker model.
  A solve can be cancelled from another thread (`cancel`), which moves the module's clock past
  its time limit; this needs shared memory.
* The page (`wasm/index.html`, run in a worker by `wasm/demo-worker.mjs`) starts with every core
  when it is cross-origin isolated, else with one thread. It generates the tables, warms up the
  engine and keeps the tables in IndexedDB for later visits. The page shows each better solution
  while solving. **Stop** cancels a solve. Without shared memory, it restarts the worker and
  reloads the stored tables.
* Tables can also be loaded from the native `gen --tier small` files: the formats are the same,
  and the module generates byte-identical files.

Measured in Bun 1.4.4 (JavaScriptCore, Safari's engine) on the 48-core EPYC VM:

| | native | wasm, 1 thread | wasm, 8 threads | wasm, 48 threads |
|---|---|---|---|---|
| table generation | 23.6 s (1T), ~1 s (48T) | 32.7 s | 31.0 s | 16.4 s |
| table load (both files) | | 0.54 s | 0.5 s | |
| search speed (10 random positions, 5 s each) | 58 / 399 / 2162 M nodes/s (1/8/48T) | 23 M nodes/s | 36 M nodes/s | 74 M nodes/s |
| random positions with any solution after 5 s | 9 / 10 / 10 of 10 | 4 | 7 | 10 |

The modules are 0.55 MB each. The single-thread module grows to 352 MiB while generating and
428 MiB while loading table files. The optimal solves of short scrambles (6-14 moves) take
milliseconds, with the same node counts as natively on one thread.

Threads help little in JavaScriptCore. Its WebAssembly atomic loads behave like locked
read-modify-writes. In a micro-benchmark, eight threads reading shared words got 7 M loads/s
each, against 226 M/s for one thread. The search and the table generator read atomics often:
the `OnceLock` statics such as `fast()`, and the stop and best-length flags. V8 (Chrome) and
SpiderMonkey (Firefox) compile atomic loads on x86 to plain loads, so they should scale like
the native build. No browser was available to measure this. Before a solve, the page runs a
short warm-up so that V8 compiles the search with its optimizing tier.

The JavaScript layer is adapted from [egg-solver](https://github.com/ArhanChaudhary/egg-solver)'s
WebAssembly build (same author). That includes the WASI shim, the thread pool, the reactor
initializer and its lessons: initial memory equals maximum memory, threads never exit, and a
warm-up solve runs for V8.

## Tests

```sh
cargo test --release
scripts/check.sh      # fmt, clippy (pedantic, -D warnings), tests, no runtime asserts
```

Non-test code does not panic: `unwrap`, `expect`, `panic!` and friends are denied by clippy,
and `scripts/check.sh` rejects runtime `assert!` outside `#[cfg(test)]` code.

The tests cover:
* the geometry, against values from an independent sticker model checked with GAP;
* notation and frame handling, via physical sticker simulation;
* coordinates and symmetry reduction;
* pruning-table consistency;
* optimal mode, against a meet-in-the-middle brute force;
* search coverage: with every pruning table zeroed, each allowed phase-1 move sequence is
  reached exactly once (checked against exact counts, at 1 and 8 threads);
* phase-2 membership, including the centres;
* every stopping rule.
