# Optimisation log

Each change was made on its own branch, checked (fmt, clippy, tests; identical search tree or
identical optimal lengths where applicable), benchmarked before/after and committed with the
benchmark table in its commit message. Branches that improved the benchmarks were merged; the
others are kept unmerged for reference.

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
| Direct distance-to-solved bound (edges 12!/2 + corners), from a competing solver | `opt/direct` | no node pruned that the existing bounds keep; optimal 1.10x / 1.12x time | no |
| Joint table also on the class-B side | `opt/bcbe` | neutral (optimal 1.02x / 1.00x) | no |
| Children visited in order of smallest bound | `opt/childorder2` | neutral (optimal 1.01x / 1.00x; 48T first solution 1.04x) | no |
| Lane-interleaved search (egg's lanes) | `opt/lanes` | slower (search is compute-bound) | no |
| Sibling pruning in the AE stage | `opt/sibling` | slower (exposes a second DRAM latency) | no |

Cumulative (48 threads, 40 random positions, 5 s): mean best length 21.60 -> 21.20; first
solution 86 ms -> 6.8 ms. Optimal mode (19-move scrambles, 48 threads): 15.1 s -> 2.76 s mean.
Last round alone (main-r2 -> micro-compute): 1.52x nodes/s at 48T and 1T, first solution 0.36x
(48T) / 0.64x (1T), optimal proofs 0.65x / 0.67x time.
