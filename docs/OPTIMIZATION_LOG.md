# Optimisation log

Each change was made on its own branch, checked (fmt, clippy, tests; identical search tree or
identical optimal lengths where applicable), benchmarked before/after and committed with the
benchmark table in its commit message. Branches that improved the benchmarks were merged; the
others are kept unmerged for reference as `neg/` (rejected) and `archive/` branches, each listed
with its result in [history/UNMERGED_EXPERIMENTS.md](history/UNMERGED_EXPERIMENTS.md).

Time-limited benchmark: seeded random positions; 48 threads 40 x 5 s and 1 thread 20 x 5 s.
Optimal benchmark: proven-optimal solves of seeded random-move scrambles (48 threads: 10 x 19
moves; 1 thread: 10 x 17 moves). The large tier is used from the joint-table change onwards.

| Change | Branch | Result | Merged |
|---|---|---|---|
| AE table checked before CB with 7+ moves left | `opt/ae-first` | +7% (48T) / +10% (1T) nodes/s | yes |
| Huge pages actually used (2 MiB-aligned buffers; madvise was failing) | `opt/hugepages` | +18% / +16% nodes/s; first solution 0.77x | yes |
| Batched child evaluation with staged software prefetch | `opt/prefetch` | 2.0x nodes/s; 5 s mean length -0.125 | yes |
| Cache-resident coordinates (label vectors, pshufb, mask ranks) | `opt/fastcoord` | +35% / +42% nodes/s; 1T first solution 0.71x | yes |
| Shallow parallel task splitting | `opt/split` | 48T first solution 0.50x (51 -> 25 ms) | yes |
| Joint corner x B-centre x edge-triple table (large tier) | `opt/cbe` | 2.05x fewer nodes; first solution 0.75x / 0.65x | yes |
| Class-B bound against the best solution (from egg's two-sided bounds) | `opt/classb` | optimal mode 0.29x (48T) / 0.67x (1T) time | yes |
| Task ordering by previous-iteration statistics (egg's job order) | `opt/order` | 5 s mean length -0.225 (48T); optimal 1T 0.69x | yes |
| Staged child evaluation (earlier AE prefetch, AE from 6 left, CBE before CB) | `opt/stages` | +10% / +14% nodes/s; optimal 0.89x / 0.88x | yes |
| Memoized split tree, per-node statistics, adaptive task count | `opt/split2` | 48T first solution 0.51x (17.8 -> 9.1 ms) | yes |
| Class-B lookups issued with the CB/CBE stages | `opt/bpath` | +6% / +4% nodes/s; optimal 0.95x | yes |
| Faster startup (parallel table read, orbit-based symmetry classes) | `opt/startup2` | large-tier startup 1.07 -> 0.33 s (48T), 7.9 -> 3.3 s (1T); search neutral | yes |
| Cheaper per-node compute (pext where fast, padded tables, deferred CB reads) | `opt/micro-compute` | 1.29x / 1.28x nodes/s; optimal 0.78x | yes |
| Direct distance-to-solved bound (edges 12!/2 + corners), from a competing solver | `neg/direct` | no node pruned that the existing bounds keep; optimal 1.10x / 1.12x time | no |
| Joint table also on the class-B side | `neg/bcbe` | neutral (optimal 1.02x / 1.00x) | no |
| Children visited in order of smallest bound | `neg/childorder2` | neutral (optimal 1.01x / 1.00x; 48T first solution 1.04x) | no |
| Lane-interleaved search (egg's lanes) | `neg/lanes` | slower (search is compute-bound) | no |
| Sibling pruning in the AE stage | `neg/sibling` | slower (exposes a second DRAM latency) | no |
| Optimality proven one iteration early, and at once on a solution of length depth + 1 | `opt/finalvar2` | random positions 419 -> 92 s, 481 -> 95 s, 39 -> 7.6 s (48T); optimal 0.61x / 0.78x | yes |
| Cross-check of the early proof (`FTO_FULL_PROOF`, `scripts/check-optimal.sh`) | `test/proofcheck` | all proven lengths agree on ten 19-move scrambles and a random position | yes |
| Corner/B-centre table reduced by the 24 class-preserving symmetries | `opt/cbsym` | 89 MB -> 4.7 MB; speed unchanged (1.01x / 1.00x nodes/s, two interleaved rounds) | yes |
| Greedy probe for a quick first solution | `neg/probe` | never finds a solution within 20000 nodes; 48T first solution 1.26x | no |
| Sibling pipeline in the search (expand / stage2 / finish) | `opt/micro3` | 1.34x nodes/s (48T and 1T); optimal 0.74x / 0.77x; first solution 0.88x / 0.79x | yes |
| Target-length passes in optimal mode (default with one thread) | `opt/target2` | optimal 1T 0.39x, 48T 0.98x | yes |
| Optional 1 GiB pages for large tables (`FTO_HUGE1G=1`) | `opt/huge1g` | A-centre/edge table on a 1 GiB page: optimal 0.96x / 0.97x, nodes/s 1.01x / 1.03x (opt-in) | yes |
| Faster table generation (egg-solver's techniques) | `opt/gen` | large tier 121 -> 30 s and peak 14.9 -> 8.1 GB; medium 16.0 -> 6.6 s (48T); tables byte-identical | yes |
| Start-up move and symmetry tables from label masks (6-bit compress/expand tables, `MaskPerm`) | `opt/startup` | start-up (2-move solve) 1T 2.66 -> 1.38 s large, 1.69 -> 0.41 s small; 48T 0.29 -> 0.26 s; WebAssembly table load (1 thread, Bun) 2.08 -> 0.54 s; search unchanged (nodes/s 1.01x / 1.01x, optimal 1.00x / 0.99x) | yes |
| CB symmetry index held in `FastTables` (no `OnceLock` read per child) | `opt/cbref` | WebAssembly 8 threads (Bun) 33.1 -> 36.7 M nodes/s (1.11x), 5 of 6 positions solved within 5 s instead of 4; native nodes/s 1.00x (48T) / 1.02x (1T) | yes |
| Both variants in one traversal (a node alive in both is evaluated once) | `opt/merged` | optimal 0.69x / 0.77x time (slices 0.85x / 0.67x; random positions 0.83x / 0.94x / 0.94x); time-limited nodes/s 0.93x / 0.94x, best after 5 s +0.09 at 48T | yes |
| `FastTables` read from `Shared` in the pipeline stages (no `fast()` `OnceLock` per stage and node) | `opt/atomics` | WebAssembly 8 threads (Bun) 34.7 -> 58.8 M nodes/s (1.70x); native nodes/s 0.98x / 1.00x (within noise) | yes |
| With slices, last proof iteration restricted to paths ending in a slice | `neg/slicelast` | neutral: few nodes reach the last level (48T node counts within 0.1%); slice optimal 0.97x (48T) / 1.05x (1T) | no |
| Profile-guided build (LLVM PGO, training on other seeds) | `neg/pgo` | nodes/s 0.987x (48T) / 0.947x (1T); optimal 1.002x / 0.978x time | no |

Cumulative, first working version -> main of 2026-11-10 ~13:30 (r10 and rc runs; 48 threads, 40
random positions, 5 s): mean best length 21.60 -> 21.10; first solution 86 ms -> 6.1 ms; 1.70e9
nodes/s. Optimal mode, 10 x 19-move scrambles at 48 threads: 15.1 s -> 1.21 s mean per position;
random positions (seeds 1-3) 419 / 481 / 39 s (morning of 2026-11-10) -> 65 / 65 / 5.3 s. Round 2
alone (main-r2 -> micro-compute): 1.52x nodes/s at 48T and 1T, first solution 0.36x (48T) /
0.64x (1T), optimal proofs 0.65x / 0.67x time.

## Table generation (`opt/gen`)

Generation time and peak memory of `fto gen --tier T` into an empty directory (wall time and
peak RSS from `/usr/bin/time -v`); every file byte-identical to the previous tables. Each
step was measured on its own build, at 48 threads unless noted.

| Step | small | medium | large (peak RSS) |
|---|---|---|---|
| main | 3.92 s | 16.04 s | 121.1 s (14.9 GB) |
| 1. backward step stops at the first frontier hit, 16 interleaved lanes with prefetch, no equivalents | 2.27 s | 11.41 s | 64.9 s |
| 2. forward writes: one `fetch_and` and a 32-entry prefetching write ring | 2.20 s | 10.22 s | 60.1 s |
| 3. CBE saturating codes made directly (no second array, no copies) | | | 56.3 s (8.1 GB) |
| 4. AE generation on the cache-resident coordinates of the search | 1.46 s | 8.66 s | 55.2 s |
| 5. label table instead of `Mset::unrank` when decoding entries | 1.24 s | 8.13 s | 47.0 s |
| 6. A part of the AE neighbours tabulated per class and move | 1.16 s | 6.69 s | 46.2 s |
| 7. CBE neighbours made per group (corner class and B colouring shared) | 1.16 s | 6.91 s | 29.5 s |
| 8. AE forward neighbours shared by the 8 entries that differ only in chirality | 1.24 s | 6.58 s | 29.7 s (8.1 GB) |

One thread, main -> step 8: small 107.2 -> 27.4 s, medium 384.6 -> 167.0 s (step 8 alone:
185.7 -> 167.0 s). Huge tier: 3 min 14 s, peak 29.4 GB, files identical (measured once, with
the build of the rejected step below, which does not change this tier; not measured on main).

Rejected: AE equivalents with the cache-resident symmetry transform (and `Mset::unrank` for the
labels) instead of the 35 MB tables of `Syms::apply_e`: medium 6.58 -> 7.29 s at 48T (the
unchanged CBE table was 6% slower in the same run too), 167.0 -> 216.7 s at 1T (AE depth 11
118.5 -> 164.0 s).
