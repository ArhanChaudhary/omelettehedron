# Unmerged experiments

[Optimisation log](../OPTIMIZATION_LOG.md)

Every optimisation is tried on its own branch, benchmarked before and after (48 threads and 1
thread) and committed with the benchmark tables in its commit message, even when it loses.
Branches that win are merged into `main` with `--no-ff`, so their commits stay in the history;
merged branch names are then deleted. Branches that do not win are kept under a name saying
why, and listed here with their date and measured result, so the record survives even if the
branches are later dropped (for example when the history is rebuilt for publication):

* `neg/<name>`: rejected or neutral experiments;
* `archive/<name>`: superseded work in progress and analysis instrumentation;
* `work/<name>`: anything else that was not merged.

Results use different baselines and workloads (see each branch's commit message), so
percentages are not additive. "Optimal" is the time to prove optimality (48 threads, 19-move
scrambles / 1 thread, 17-move scrambles).

## Rejected and neutral experiments

| Experiment | Date | Summary |
|---|---|---|
| <a id="ref-neg-wasm-plainload"></a>`neg/wasm-plainload` | 2026-11-11 | Table generation in WebAssembly with threads: table words read with plain loads instead of atomic ones (a racy read, experiment only), to test whether atomic loads of the table made threads slower than one thread: no change (6 threads 47.3 s, 2 threads 43.3 s, against 48.6 s and 51.3 s); the cause was the per-entry `OnceLock` reads (`opt/wasm-gen`) |
| <a id="ref-neg-fwd-ratio"></a>`neg/fwd-ratio` | 2026-11-11 | Table generation: `FORWARD_RATIO` 2 instead of 4 (small-tier AE depth 10 and CB depth 6 backward instead of forward), for one WebAssembly thread: much slower (WebAssembly 1 thread 29.8-30.8 -> 51.9-52.1 s, 6 threads 6.7-7.2 -> 15.1-15.9 s; native 1 thread 21.3 -> 37.1 s, 48 threads 1.08 -> 1.27 s; AE depth 10 backward 44.8 s against 23.6 s forward) |
| <a id="ref-neg-shared-em"></a>`neg/shared-em` | 2026-11-11 | Two-phase solver: each child's edge label masks (`fast::masks`) computed once and kept in the frame for the A-centre/edge index, the joint table's four triple indices and the distance-to-solved lookups, which each computed them themselves: slower, optimal 1.032x / 1.032x time at 1T / 48T, time limit 5 s nodes/s 0.974x (48T) / 0.973x (1T); computed in `expand` only when the AE stage comes first and otherwise for the CB survivors: 1.020x / 1.039x, nodes/s 0.972x / 0.975x (three rounds each against main 9a9262f; identical node counts, same mean lengths). Storing and reloading the masks costs more than the four SIMD compares that compute them where they are used |
| <a id="ref-neg-cbe-nextpf"></a>`neg/cbe-nextpf` | 2026-11-11 | Two-phase solver: each child's joint corner/B-centre/edge-triple (CBE) entries prefetched inside the AE/CB loop as soon as it passes, where `finish` reads that table itself (budget <= 8 and not prefetched in `stage2`), instead of in a loop of their own right before the reads (the profile put 8.1% of the samples in optimal mode and 17.9% under a time limit on that read): optimal 0.999x / 1.006x time at 1T / 48T (three rounds against main 9a9262f), time limit 5 s: same mean lengths, first solution 0.0057 -> 0.0059 s (48T) and 0.237 -> 0.241 s (1T), nodes/s 1.000x / 0.980x; identical node counts. As in the optimal solver, the read stall does not shrink when the prefetch comes earlier |
| <a id="ref-neg-tosolved-wave2"></a>`neg/tosolved-wave2` | 2026-11-11 | Two-phase solver: the second wave of distance-to-solved lookups issued before the CB and class-B stages (the first wave read at the start of `finish`), so that its latency overlaps them: 0.985x / 0.986x time at 1T / 48T against main, the same as the first wave alone in `stage2` (0.982x / 0.985x; rounds mixed), identical node counts |
| <a id="ref-neg-lanes"></a>`neg/lanes` | 2026-11-10 | Lane-interleaved phase-1 search (several paths advanced in turn): slower; the search is compute-bound (~110 cycles per child). 48T nodes/s 0.94x, first solution 2.70x; optimal 1.07x |
| <a id="ref-neg-sibling"></a>`neg/sibling` | 2026-11-10 | Sibling pruning in the AE stage: slower, it exposes a second DRAM latency. Nodes/s 0.95x (48T and 1T) |
| <a id="ref-neg-direct"></a>`neg/direct` | 2026-11-10 | Direct distance-to-solved bound (edge permutation 12!/2 + corners, from a competing solver's design): prunes nothing beyond max(distance to H_A, distance to H_B). Optimal 1.10x / 1.12x, nodes/s 0.91x |
| <a id="ref-neg-bcbe"></a>`neg/bcbe` | 2026-11-10 | Joint corner/B-centre/edge-triple table also on the class-B side: neutral (optimal 1.02x / 1.00x) |
| <a id="ref-neg-childorder2"></a>`neg/childorder2` | 2026-11-10 | Children visited in order of their smallest exact bound: neutral (optimal 1.01x / 1.00x; 48T first solution 1.04x) |
| <a id="ref-neg-probe"></a>`neg/probe` | 2026-11-10 | Greedy depth-first probe (20000 nodes, children by bound) for a quick first solution before the iterations: it practically never finds a phase-1 path, and only delays the iterations (48T first solution 1.26x, positions solved after 0.01 s 31 -> 26) |
| <a id="ref-neg-optimal-cbe"></a>`neg/optimal-cbe` | 2026-11-11 | Optimal solver: the two-phase joint corner/B-centre/edge-triple table on both sides as a saturating second stage (`large` tier, +7 GB): 0.978x nodes on the 17-move set (0.9998x on the 15-move set) on top of medium + A-centre/edge bound; hardly any node has fewer than 9 moves left, where its codes can prune (0.0002 lookups per node). Not worth 7 GB |
| <a id="ref-neg-optimal-lanes"></a>`neg/optimal-lanes` | 2026-11-11 | Optimal solver: 2-8 searches per worker interleaved a node at a time (explicit stacks; one lane = today's recursion), with and without the sibling pipeline inside each lane: same trees, but 0.97-1.04x time on medium and 1.00-1.04x on large for 2 or 3 lanes, worse with more. The read stall is almost all in the first stage (26% / 35% of the samples, medium / large); a node takes 40-62 ns, less than a DRAM miss, yet 2-8 nodes of slack do not hide it, and 48 threads run 43x one thread's nodes (not bandwidth). Address translation or the lanes' own working set are the likely causes; this VM has no hardware counters to tell |
| <a id="ref-neg-optimal-r3-dynorder"></a>`neg/optimal-r3-dynorder` | 2026-11-11 | Optimal solver: each node's own lookup order, the lookups that cannot prune its children (`h[k] + 1 <= budget`) last and those at its limit (`h[k] == budget + 1`) first: 3% / 9% fewer lookups (medium / large), same nodes, but mixed time (medium 0.993x / 0.987x, large 1.029x / 0.959x at 1T / 48T; the deferral alone 1.002x / 1.008x, 1.008x / 0.961x) |
| <a id="ref-neg-optimal-r3-aeclass"></a>`neg/optimal-r3-aeclass` | 2026-11-11 | Optimal solver: a pass over each stage's children that computes their A-centre rank and prefetches its symmetry-class entry before the A-centre/edge indices (the profile showed 6.3% of the medium tier's samples on that 1.5 MB table): slower, 1.104x / 1.120x time (medium, 1T / 48T), 1.007x / 1.009x (large); the extra pass costs more than the stall it hides |

## Archived experiments

| Experiment | Date | Summary |
|---|---|---|
| <a id="ref-archive-tosolved-est"></a>`archive/tosolved-est` | 2026-11-11 | Two-phase solver: exact estimate of a distance-to-solved bound on the joint table's index before building it (feature `est`: ball of radius 8 + radius-9 bit set around the solved projection, canonical search above, per-candidate subtree accounting; also cross-checks the built table: 0 mismatches in 1,684,303 children of the 17-move searches). 17-move set: exact distances 0.760x nodes, one side only 0.820x, 2-bit codes of base 8 / 9 / 10 / 11: 0.800x / 0.760x / 0.761x / 0.935x. Led to `opt/tosolved` |
| <a id="ref-archive-childorder-wip"></a>`archive/childorder-wip` | 2026-11-10 | First child-ordering attempt (its key used deferred CB values); superseded by `neg/childorder2` |
| <a id="ref-neg-slicelast"></a>`neg/slicelast` | 2026-11-10 | With slices, the last proof iteration (and depth T of a target pass) restricted to phase-1 paths ending in a slice: correct (check-optimal large slice) but neutral, since few nodes reach the last level (48T node counts within 0.1%; slice optimal 17-move set 0.97x at 48T, 15-move set 1.05x at 1T, time-limited unchanged) |
| <a id="ref-neg-pgo"></a>`neg/pgo` | 2026-11-10 | Profile-guided build (`scripts/build-pgo.sh`): neutral to worse (time-limited nodes/s 0.987x at 48T, 0.947x at 1T; optimal 1.002x / 0.978x time; two interleaved rounds) |
| <a id="ref-neg-pipe4"></a>`neg/pipe4` | 2026-11-10 | First lookups issued three siblings ahead (4 frames per level): neutral (nodes/s 1.008x / 0.997x; optimal 1.009x / 0.983x; identical tree) |
| <a id="ref-neg-lpass"></a>`neg/lpass` | 2026-11-10 | `-l L` as a target pass at L after the first solution: neutral to worse (48T -l 21 16.73 -> 16.84 s, 1T -l 22 10.09 -> 11.07 s mean time to a solution of at most L) |
| <a id="ref-neg-tlpasses"></a>`neg/tlpasses` | 2026-11-10 | Target-length passes under a time limit (after the first solution / after a stalled iteration / when a forecast says the proof fits): proofs that fit 1.27x faster at 1T, but answers 0.55-1.0 moves worse when they do not; the forecast never triggers |
| <a id="ref-archive-nodecount-v1"></a>`archive/nodecount-v1` | 2026-11-10 | Search-coverage test and full phase-2 membership check before the no-panic merge; superseded by the rebased `test/nodecount2` (benchmarked neutral) |
| <a id="ref-archive-analysis-round3"></a>`archive/analysis-round3` | 2026-11-10 | Round-3 analysis instrumentation: `stats` feature (evaluations and survivors by remaining moves), FTO_DELAY memory-latency probe, early-stop and target-length-pass experiments, CB reduced by 24 symmetries (+1.3% / +1.0% nodes/s, 4.7 MB), CB x A-colour joint table (rejected: time +2-4.5%), forced-route bounds (<= 0.1% of nodes), duplicate-sequence growth (0.62% excess at depth 6) |
| <a id="ref-archive-finalvar-v1"></a>`archive/finalvar-v1` | 2026-11-10 | Early end of optimality proofs before the rebase onto the no-panic and test merges; its commits were rebased into the merged `opt/finalvar2` (random positions 4.5-5.1x faster) |
| <a id="ref-archive-target-v1"></a>`archive/target-v1` | 2026-11-10 | Target-length passes before the sibling pipeline (tuning runs for gaps 0/2/3); re-applied as the merged `opt/target2` (1T optimal 0.39x) |
| <a id="ref-work-pure"></a>`work/pure` | 2026-11-10 | Start of a pure-FTO mode (labelled centres): needs an exact pure-H_A membership test (75.6 M states) before the search can use it; waiting for a decision |
