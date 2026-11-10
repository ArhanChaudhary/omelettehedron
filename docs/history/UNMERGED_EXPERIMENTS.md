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
| <a id="ref-neg-lanes"></a>`neg/lanes` | 2026-11-10 | Lane-interleaved phase-1 search (egg-solver's lanes): slower; the search is compute-bound (~110 cycles per child). 48T nodes/s 0.94x, first solution 2.70x; optimal 1.07x |
| <a id="ref-neg-sibling"></a>`neg/sibling` | 2026-11-10 | Sibling pruning in the AE stage (egg-solver): slower, it exposes a second DRAM latency. Nodes/s 0.95x (48T and 1T) |
| <a id="ref-neg-direct"></a>`neg/direct` | 2026-11-10 | Direct distance-to-solved bound (edge permutation 12!/2 + corners, from a competing solver's design): prunes nothing beyond max(distance to H_A, distance to H_B). Optimal 1.10x / 1.12x, nodes/s 0.91x |
| <a id="ref-neg-bcbe"></a>`neg/bcbe` | 2026-11-10 | Joint corner/B-centre/edge-triple table also on the class-B side: neutral (optimal 1.02x / 1.00x) |
| <a id="ref-neg-childorder2"></a>`neg/childorder2` | 2026-11-10 | Children visited in order of their smallest exact bound: neutral (optimal 1.01x / 1.00x; 48T first solution 1.04x) |
| <a id="ref-neg-probe"></a>`neg/probe` | 2026-11-10 | Greedy depth-first probe (20000 nodes, children by bound) for a quick first solution before the iterations: it practically never finds a phase-1 path, and only delays the iterations (48T first solution 1.26x, positions solved after 0.01 s 31 -> 26) |

## Archived experiments

| Experiment | Date | Summary |
|---|---|---|
| <a id="ref-archive-childorder-wip"></a>`archive/childorder-wip` | 2026-11-10 | First child-ordering attempt (its key used deferred CB values); superseded by `neg/childorder2` |
| <a id="ref-archive-nodecount-v1"></a>`archive/nodecount-v1` | 2026-11-10 | Search-coverage test and full phase-2 membership check before the no-panic merge; superseded by the rebased `test/nodecount2` (benchmarked neutral) |
| <a id="ref-archive-analysis-round3"></a>`archive/analysis-round3` | 2026-11-10 | Round-3 analysis instrumentation: `stats` feature (evaluations and survivors by remaining moves), FTO_DELAY memory-latency probe, early-stop and target-length-pass experiments, CB reduced by 24 symmetries (+1.3% / +1.0% nodes/s, 4.7 MB), CB x A-colour joint table (rejected: time +2-4.5%), forced-route bounds (<= 0.1% of nodes), duplicate-sequence growth (0.62% excess at depth 6) |
