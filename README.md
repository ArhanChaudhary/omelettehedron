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

`Xw` turns face X together with the middle layer (a wide turn). `Xs` turns only the middle
layer between X and its opposite face, in the sense of X, so `Us` = `Ds'` (a slice turn).
Scrambles may always contain wide and slice turns; solutions name slices by their front face
(`Us`, `Rs`, `Fs`, `Ls`).

Face turns can rotate the whole puzzle, so a solution ends in the solved state in *some*
orientation. Internally the front corner is held fixed.

### Move sets

Solutions use face turns, plus the families given with `--moves` (comma-separated). Every
allowed move counts as one.

| `--moves` | Solutions use | Effect |
|---|---|---|
| (none) | face turns | the face-turn metric |
| `slice` | face and slice turns | shorter solutions; own pruning tables (`*-slice.prun`) |
| `wide` | face and wide turns | same lengths and search; see below |
| `wide,slice` | all three | as `slice`, printed as with `wide` |

A wide turn is a turn of the opposite face plus a whole-puzzle rotation, and the solver treats
every orientation as solved, so wide turns never shorten a solution. With `wide` the search is
unchanged and only the output differs: the corner piece between U, R, F and L stays where the
scramble left it, and turns of the faces around it are printed as wide turns of the opposite
faces (`Dw` for `U` if the scramble left that piece in front).

With `slice` the phase-1 search has 24 internal moves instead of 16, and slices are not
phase-2 moves; see [docs/PHASES.md](docs/PHASES.md#8-move-sets---moves) for the move-order
rules, the tables and why optimal mode stays exact.

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

# optimal solution with slice turns allowed (generates the slice tables on first use)
target/release/omelettehedron solve --tier large --moves slice "U R' BL F D B' L BR"
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

All tiers also use the 89 MB corner/B-centre table and the in-memory phase-2 table. With
`--moves slice` every tier has its own set of tables of the same sizes (names ending in
`-slice`); `gen --tier T --moves slice` generates them ahead of time.
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
* optimal mode, against a meet-in-the-middle brute force, with and without slice turns;
* the move-order rules: per axis against the axis group, and the canonical sequences of each
  length up to 5 (4 with slices) against breadth-first search (every position reached; exactly
  once up to length 4);
* search coverage: with every pruning table zeroed, each allowed phase-1 move sequence is
  reached exactly once (checked against exact counts, at 1 and 8 threads, with and without
  slices);
* phase-2 membership, including the centres;
* every stopping rule.
