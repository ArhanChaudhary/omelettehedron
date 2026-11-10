# fto-solver

A fast two-phase solver for the face-turning octahedron (FTO), written in Rust.

* **Phase 1** brings the puzzle into H_A, the positions that can be solved by turning
  only the four class-A faces U, F, BR and BL.
* **Phase 2** finishes with those faces only. Its 933,120 positions are fully tabulated, so it
  is an O(length) lookup, and its solutions are optimal even among all-move solutions.

Phase 1 is an IDA* search guided by two pruning tables. See [docs/PHASES.md](docs/PHASES.md)
for why this split was chosen and how it compares with alternatives.

## Build

```sh
cargo build --release
```

The build targets `x86-64-v3` (AVX2/BMI2, any x86-64 CPU from about 2015 on); see
`.cargo/config.toml`.

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
target/release/fto solve --tier large "U R' BL F D B' L BR"

# stop after 2 seconds, printing every improvement as it is found
target/release/fto solve --tier large -t 2 -p improvements "R U' L F' BR B' ..."

# stop at the first solution of at most 22 moves
target/release/fto solve --tier large -l 22 "..."

# stop after 1000 solutions (every phase-2 completion counts, even if longer than the best)
target/release/fto solve --tier large -n 1000 "..."

# one scramble per line on stdin
target/release/fto solve --tier large -t 1 < scrambles.txt
```

Stopping rules can be combined; the search stops when any of them triggers. `-j N` sets the
thread count (default: all cores).

### Pruning tables

`--tier` is required (there is no default size): choose the largest tier that fits
comfortably in RAM. Tables are generated on first use. They are written to `./tables`, or to `--tables DIR` or
`$FTO_TABLES` if given. `fto gen --tier T` generates them ahead of time.

| `--tier` | Memory | Generation (48 cores) | Contents |
|---|---|---|---|
| `small` | ~260 MB | ~5 s | A-centre/edge table without edge chirality |
| `medium` | ~1.2 GB | ~20 s | full A-centre/edge table (288-fold symmetry) |
| `large` | ~8 GB | ~2.5 min | medium + joint corner/B-centre/edge-triple table (about 2x fewer nodes) |
| `huge` | ~31 GB | ~5 min | relabel-only A-centre/edge table (cheaper lookups) + the joint table |

All tiers also use the 89 MB corner/B-centre table and the in-memory phase-2 table.

## Benchmarks

```sh
target/release/fto bench --tier large -N 40 -t 5       # 40 seeded random positions, 5 s each
target/release/fto bench --tier large -N 20 -t 5 -j 1  # single-threaded
scripts/bench.sh <label> [tier]                 # the standard set used during optimisation
scripts/bench-opt.sh <label> [tier]             # optimal mode: time to prove optimality
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
* optimal mode, against a meet-in-the-middle brute force;
* every stopping rule.
