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
| Class-B bound against the best solution (two-sided bounds) | `opt/classb` | optimal mode 0.29x (48T) / 0.67x (1T) time | yes |
| Task ordering by previous-iteration statistics | `opt/order` | 5 s mean length -0.225 (48T); optimal 1T 0.69x | yes |
| Staged child evaluation (earlier AE prefetch, AE from 6 left, CBE before CB) | `opt/stages` | +10% / +14% nodes/s; optimal 0.89x / 0.88x | yes |
| Memoized split tree, per-node statistics, adaptive task count | `opt/split2` | 48T first solution 0.51x (17.8 -> 9.1 ms) | yes |
| Class-B lookups issued with the CB/CBE stages | `opt/bpath` | +6% / +4% nodes/s; optimal 0.95x | yes |
| Faster startup (parallel table read, orbit-based symmetry classes) | `opt/startup2` | large-tier startup 1.07 -> 0.33 s (48T), 7.9 -> 3.3 s (1T); search neutral | yes |
| Cheaper per-node compute (pext where fast, padded tables, deferred CB reads) | `opt/micro-compute` | 1.29x / 1.28x nodes/s; optimal 0.78x | yes |
| Direct distance-to-solved bound (edges 12!/2 + corners), from a competing solver | `neg/direct` | no node pruned that the existing bounds keep; optimal 1.10x / 1.12x time | no |
| Joint table also on the class-B side | `neg/bcbe` | neutral (optimal 1.02x / 1.00x) | no |
| Children visited in order of smallest bound | `neg/childorder2` | neutral (optimal 1.01x / 1.00x; 48T first solution 1.04x) | no |
| Lane-interleaved search (several paths advanced in turn) | `neg/lanes` | slower (search is compute-bound) | no |
| Sibling pruning in the AE stage | `neg/sibling` | slower (exposes a second DRAM latency) | no |
| Optimality proven one iteration early, and at once on a solution of length depth + 1 | `opt/finalvar2` | random positions 419 -> 92 s, 481 -> 95 s, 39 -> 7.6 s (48T); optimal 0.61x / 0.78x | yes |
| Cross-check of the early proof (`FTO_FULL_PROOF`, `scripts/check-optimal.sh`) | `test/proofcheck` | all proven lengths agree on ten 19-move scrambles and a random position | yes |
| Corner/B-centre table reduced by the 24 class-preserving symmetries | `opt/cbsym` | 89 MB -> 4.7 MB; speed unchanged (1.01x / 1.00x nodes/s, two interleaved rounds) | yes |
| Greedy probe for a quick first solution | `neg/probe` | never finds a solution within 20000 nodes; 48T first solution 1.26x | no |
| Sibling pipeline in the search (expand / stage2 / finish) | `opt/micro3` | 1.34x nodes/s (48T and 1T); optimal 0.74x / 0.77x; first solution 0.88x / 0.79x | yes |
| Target-length passes in optimal mode (default with one thread) | `opt/target2` | optimal 1T 0.39x, 48T 0.98x | yes |
| Optional 1 GiB pages for large tables (`FTO_HUGE1G=1`) | `opt/huge1g` | A-centre/edge table on a 1 GiB page: optimal 0.96x / 0.97x, nodes/s 1.01x / 1.03x (opt-in) | yes |
| Faster table generation (backward steps stopping at the first frontier hit, prefetching lanes, in-place joint table) | `opt/gen` | large tier 121 -> 30 s and peak 14.9 -> 8.1 GB; medium 16.0 -> 6.6 s (48T); tables byte-identical | yes |
| Start-up move and symmetry tables from label masks (6-bit compress/expand tables, `MaskPerm`) | `opt/startup` | start-up (2-move solve) 1T 2.66 -> 1.38 s large, 1.69 -> 0.41 s small; 48T 0.29 -> 0.26 s; WebAssembly table load (1 thread, Bun) 2.08 -> 0.54 s; search unchanged (nodes/s 1.01x / 1.01x, optimal 1.00x / 0.99x) | yes |
| CB symmetry index held in `FastTables` (no `OnceLock` read per child) | `opt/cbref` | WebAssembly 8 threads (Bun) 33.1 -> 36.7 M nodes/s (1.11x), 5 of 6 positions solved within 5 s instead of 4; native nodes/s 1.00x (48T) / 1.02x (1T) | yes |
| Both variants in one traversal (a node alive in both is evaluated once) | `opt/merged` | optimal 0.69x / 0.77x time (slices 0.85x / 0.67x; random positions 0.83x / 0.94x / 0.94x); time-limited nodes/s 0.93x / 0.94x, best after 5 s +0.09 at 48T | yes |
| `FastTables` read from `Shared` in the pipeline stages (no `fast()` `OnceLock` per stage and node) | `opt/atomics` | WebAssembly 8 threads (Bun) 34.7 -> 58.8 M nodes/s (1.70x); native nodes/s 0.98x / 1.00x (within noise) | yes |
| Stop flag and class-B limit cached per worker (refreshed every 16384 nodes, per task and after each leaf) | `opt/stopcache` | WebAssembly 8 threads (Bun) 58.8 -> 144 M nodes/s (2.4x); native nodes/s 1.01x / 1.02x, optimal 0.99x / 0.96x time | yes |
| Cancellation without a time limit (`Limits::cancel_after`): WebAssembly optimal solves keep the one-thread target-length passes | `fix/cancel` | demo test, 14-move scramble 467,558 -> 169,652 nodes, 11-move 582 -> 486 (one thread; natively the passes take 0.39x time) | yes |
| Target-length passes on by default at every thread count (gap 2; was one thread only) | `opt/gap48` | 48T optimal: 19-move set 0.954x, random positions 1-3 0.99x / 0.97x / 0.95x time | yes |
| With slices, last proof iteration restricted to paths ending in a slice | `neg/slicelast` | neutral: few nodes reach the last level (48T node counts within 0.1%); slice optimal 0.97x (48T) / 1.05x (1T) | no |
| Profile-guided build (LLVM PGO, training on other seeds) | `neg/pgo` | nodes/s 0.987x (48T) / 0.947x (1T); optimal 1.002x / 0.978x time | no |
| Sibling pipeline three deep (first lookups three siblings ahead) | `neg/pipe4` | nodes/s 1.008x / 0.997x; optimal 1.009x / 0.983x: the table-read stall is not prefetch timing | no |
| `-l L` searched as a target pass at L once a first solution exists | `neg/lpass` | time to a solution of at most L 1.01x (48T, -l 21) / 1.10x (1T, -l 22): the search already prunes against its best | no |
| Target-length passes under a time limit (three switching rules) | `neg/tlpasses` | proof within the limit 1T 0.80 -> 0.63 s, but best after 5 s 21.45 -> 22.20 (48T), 22.00 -> 22.55 (1T): a time limit keeps the improving search | no |
| Tables mapped from the files instead of read into huge-page memory (start-up) | (analysis) | start-up at 48T is 0.26 s, mostly kernel time (faults, zeroing, copying 8 GB); but file mappings get only 4 KB pages here (FilePmdMapped 0, THPeligible 0), and huge pages were worth +16-18% nodes/s | no |
| Class-B bound of one-side nodes from the CB table only (skip the other side's AE read) | (analysis) | the other side's AE values are mod 3 and decoded from the parent's exact value, so skipping one level loses the chain below it; AE is also the stronger half of that bound (mean 10.9 vs 6.1) | no |
| `AE_FIRST` 5 or 7 instead of 6 (after the merged traversal) | (no branch: constant) | optimal 1.03x / 1.01x (5) and 1.04x / 1.01x (7) time; nodes/s 0.97-1.02x | no |
| `FTO_TASKN` 128 or 2048 instead of 512 | (no branch: env knob) | optimal 19-move set 1.04x / 1.03x, random position 3 0.98x / 0.97x | no |
| Distance-to-solved bound in optimal mode: 1-bit table (distance >= 11) on the joint table's index, 8 lookups (4 triples x 2 sides), read with at most 10 moves left; large and huge tiers, +3.5 GB, read on first use (estimate first: `archive/tosolved-est`) | `opt/tosolved` | 1T nodes 0.743x / 0.771x (17 / 19-move sets); optimal 0.788x (1T) / 0.853x (48T) time with the table read on first use (0.784x / 0.841x read at start-up), random positions 0.82-0.90x, edge_hard 0.85x, corner 6-flip 626.5 -> 117.2 s; time-limited and start-up unchanged | yes |
| Distance-to-solved lookups in two waves: triple 2 of both sides first (issued in `stage2`), the other six only for the children they keep (27.8% on the 17-move set; table reads 0.459x) | `opt/tosolved-order` | optimal 0.982x (1T) / 0.985x (48T) time, 4 rounds, identical node counts; both waves in `finish`: 1.007x / 1.008x; second wave issued before the CB / class-B stages: 0.985x / 0.986x, no better (`neg/tosolved-wave2`) | yes |
| Joint-table (CBE) entries prefetched inside the AE/CB loop as soon as a child passes, instead of in a loop of their own before the reads (profile: 8.1% / 17.9% of the samples on that read, optimal / time-limited) | `neg/cbe-nextpf` | optimal 0.999x (1T) / 1.006x (48T) time, identical node counts; time limit 5 s: mean length unchanged (21.450 at 48T, 22.000 at 1T), nodes/s 1.000x / 0.980x | no |
| Distance-to-solved bound under a time limit and `-l` too (it prunes against the best, as admissible as a target pass; only counting solutions turns it off) | `opt/tl-cbs` | corner 6-flip `-t 1800` 624 -> 84.0 s; edge_hard 67.0 -> 56.4 s; random positions `-t 1200` 0.85x / 0.80x / 0.80x; best lengths at 0.1-3 s unchanged (48T, 1T, 2 rounds) | yes |
| Edge label masks computed once per child and shared by the A-centre/edge, joint-table and distance-to-solved indices (kept in the frame) | `neg/shared-em` | identical node counts; optimal 1.032x (1T) / 1.032x (48T) time, time limit 5 s nodes/s 0.974x (48T) / 0.973x (1T); computed only where every child needs them: 1.020x / 1.039x, nodes/s 0.972x / 0.975x; mean lengths unchanged | no |

Cumulative, first working version -> main 97452e0 of 2026-11-10 (48 threads, 40 random
positions, 5 s): mean best length 21.60 -> 21.175; first solution 86 ms -> 6.0 ms. Optimal mode,
10 x 19-move scrambles at 48 threads: 15.1 s -> 0.80 s mean per position; random positions
(seeds 1-3) 419 / 481 / 39 s (morning of 2026-11-10) -> 55 / 64 / 5.3 s (single runs vary by a
few seconds). Round 2 alone
(main-r2 -> micro-compute): 1.52x nodes/s at 48T and 1T, first solution 0.36x (48T) / 0.64x
(1T), optimal proofs 0.65x / 0.67x time.

Per-node work on the two-phase solver closed on 2026-11-11. The profile of the large tier (perf,
task-clock sampling, the VM has no hardware counters; 1 thread), optimal mode on the 19-move set
/ time limit 1 s: A-centre/edge and CB index arithmetic 16.0% / 13.1% of the samples, the stall
on their table reads 12.8% / 14.2%, prefetch instructions (line fill buffers full) 12.1% / 9.7%,
the joint table's reads 8.1% / 17.9%, the rest thinly spread (move application, the stack). The
two changes aimed at it, prefetching the joint table's entries earlier and computing the edge
masks once per child, did not pay, like the earlier `neg/pipe4` and `neg/sibling`: the read
stall does not shrink with more slack, and the arithmetic left is spread over many small steps.

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

## Optimal solver

The separate optimal solver (`src/optimal/`, [OPTIMAL.md](OPTIMAL.md)), branch `optimal/v0`.
Measured by node counts on one thread (deterministic) while the machine was busy; timings
follow when it is free.

| Change | Branch | Result (nodes, 10 x 15-move scrambles, seed 1717) | Merged |
|---|---|---|---|
| v0: IDA* to solved; corners x centres (both orbits, one table) + edge permutation (small tier) | `optimal/v0` | 8.73e7 per position | (branch) |
| Corners x edges table, 48 symmetries (medium tier) | `optimal/v0` | 2.87e6 (0.033x) | (branch) |
| Symmetric start (keep the first task of each class of conjugate positions) | `optimal/v0` | random scrambles unchanged (no symmetry); hard positions with 48 / 24 / 16 / 8 symmetries 0.021x / 0.042x / 0.063x / 0.125x nodes (iterations to 15, 1 thread) | (branch) |
| The two-phase solver's AE-to-`H_A` table as two more lookups (the position and its `SWAP_AB` conjugate; small tier: the triples table, medium: with chirality) | `optimal/ae` | small 1.35e6 (0.0154x), medium 1.16e5 (0.0405x); 17-move set, medium: 1.48e7 per position (0.0383x); same lengths; the 6-flip unchanged (AE 0 at its root) | (branch) |
| The two-phase joint corner/B-centre/edge-triple table (CBE) on both sides as a saturating second stage, `large` tier (+7 GB) | `neg/optimal-cbe` | 0.9998x (15-move set), 0.978x (17-move set) on top of medium + AE; 0.0002 CBE lookups per node (hardly any node has fewer than 9 moves left) | no |
| Lookups in the order of their pruning power with early exits (AE, its twin, corners x edges, corners x centres; from per-child pruning statistics) | `optimal/order` | nodes unchanged; lookups 0.305x / 0.292x (medium / small, 15-move set), 0.302x (17-move set): 1.52 / 1.46 / 1.51 per node instead of 5; time 0.818x (1T, 17-move set) / 0.811x (48T, 19-move set), medium, quiet machine | (branch) |
| Exact checks near the leaves (one move left: is the child one of the 16 neighbours of solved; or a ball of positions within a few moves) | (analysed, not pursued) | the pruning statistics show 92% of the children with 9 or 10 moves left and fewer than 0.01% with at most 6 (about 130 per 1.16 million with one left), so no measurable change in nodes | no |
| Slice turns (`--moves slice`): 24 moves, slice-metric tables (own and two-phase AE) | `optimal/slice` | 15-move set: mean optimal 13.6 (14.5 without slices; the two-phase solver's lengths on all ten); small 1.30e9 nodes (96x the face metric), medium 1.40e8 (120x), 1.34 / 1.37 lookups per node; 17-move set, medium: mean 15.5 (16.5), 4.70e10 nodes (318x); face-metric time against main (with round 3) 1.007x (1T, 17-move set) / 1.001x (48T, 19-move set), same nodes | (branch) |
| Corner x B-centre x edge-triple table to solved (ordered positions of an A-face's edges, 6 symmetries, 21.0 GB), 8 lookups (4 A-triples, 4 B-triples through `SWAP_AB`), `large` tier (24.4 GB); chosen by an exact estimate on the recorded search (0.34x, against 0.62x for positions as a set + chirality at 7.0 GB) | `optimal/cbe-solved` | nodes 0.549x / 0.515x (15 / 17-move sets, first solution; full trees 0.289x / 0.336x), 0.249x (19-move set, 48T); time 0.885x (1T, 17-move set), 0.470x (48T, 19-move set); corner 6-flip 214x fewer nodes, 29.0 -> 0.50 s (48T) | (branch) |
| Optimal solver (large tier) against the two-phase solver's optimal mode (large tier), same positions, two rounds, quiet machine | `optimal/vs-two-phase` | the two-phase solver 2.5x faster on the 17-move set (1T), 1.8x on the 19-move set (48T), 1.7x / 3.5x / 21x on random positions of 21 / 21 / 20 moves, 4.8x on edge_hard; the optimal solver 1,310x faster on the corner 6-flip (0.48 s against 626.5 s) | (comparison) |
| Next-stage entries prefetched as soon as a child passes a stage (profile: 36-44% of `dfs` samples stalled on the table read) | `optimal/r3-prefetch` | nodes unchanged; time 0.98x / 0.99x (medium 1T / 48T), 0.77-0.82x / 0.82-0.84x (large) | (branch) |
| Sibling pipeline: the next surviving child's children generated and their first-stage entries prefetched before the search descends into the current one | `optimal/r3-prefetch` | nodes unchanged; with the next-stage prefetch 0.928x / 0.945x (medium 1T / 48T), 0.818x / 0.875x (large); alone 0.953x / 0.940x, 0.890x / 0.973x | (branch) |
| Piece positions from SIMD compares instead of inverse permutations (A-centre/edge chirality bits, triple table index; 15% of the medium tier's samples) | `optimal/r3-masks` | nodes unchanged; time 0.885x / 0.885x (medium 1T / 48T), 0.907x / 0.925x (large) against round 3 | (branch) |
| The first stage's A-centre/edge class entries prefetched while the sibling pipeline generates the children (Bounds::pre / index_with) | `optimal/r3-aepre` | nodes unchanged; time 0.898x / 0.919x (medium 1T / 48T), 0.946x / 0.941x (large) on top of the SIMD positions | (branch) |
| Per-node lookup order: the lookups that cannot prune a node's children last, those at its limit first | [`neg/optimal-r3-dynorder`](history/UNMERGED_EXPERIMENTS.md#ref-neg-optimal-r3-dynorder) | lookups 0.97x (medium) / 0.91x (large), same nodes; time mixed: medium 0.993x / 0.987x, large 1.029x / 0.959x (1T / 48T) | no |
| The large tier without corners x edges (it prunes 0.3% of the children after the others): 22.0 GB instead of 24.4 | `optimal/large-noce` | nodes 1.028x / 1.001x (17 / 19-move sets), time 0.979x (1T) / 0.954x (48T) | (branch) |
| Lanes: several tasks per worker as explicit stacks, interleaved a node at a time (with and without the sibling pipeline inside) | [`neg/optimal-lanes`](history/UNMERGED_EXPERIMENTS.md#ref-neg-optimal-lanes) | completed-iteration nodes unchanged; time 0.97-1.04x (medium), 1.00-1.04x (large) for 2 or 3 lanes: the first-stage stall is not hidden by more work in flight | no |
| Hot loops aligned to 64 bytes (`-C llvm-args=-align-loops=64`), both solvers (main 8834eed, the same flags otherwise) | (measured, not adopted) | optimal solver: medium 0.997x / 0.998x, large 0.995x / 0.996x (1T 17-move / 48T 19-move set); two-phase optimal mode, large tier: 0.999x (1T, 17-move set) / 1.005x (48T, 19-move set); three interleaved rounds each: neutral, so the build flags stay | no |
| 1 GiB pages for the tables (`FTO_HUGE1G=1`, already supported by the shared table loader) | (measured, no change) | medium with the whole working set on 1 GiB pages: 0.996x (1T) / 0.983x (48T); large with only the A-centre/edge table on one (5 pages could be reserved): 0.992x / 0.995x; the first-stage read stall is not address translation | no |
| Prefetch pass for the A-centre/edge lookups' symmetry-class entries (Bounds::pre) | [`neg/optimal-r3-aeclass`](history/UNMERGED_EXPERIMENTS.md#ref-neg-optimal-r3-aeclass) | 1.104x / 1.120x (medium 1T / 48T), 1.007x / 1.009x (large): the extra pass costs more than the stall | no |
