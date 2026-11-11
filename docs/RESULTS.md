# Results

[Documentation home](../README.md) · [Optimisation log](OPTIMIZATION_LOG.md) ·
[Phase design](PHASES.md)

Measured on a 48-core AMD EPYC VM (94 GB RAM, no hardware performance counters), large tier
(A-centre/edge 1.03 GB + corner/B-centre + joint 7.0 GB tables), tables loaded before timing
(loading takes 0.33 s with 48 threads). Positions are seeded random states unless stated
otherwise; the numbers below are from the standard suites (`scripts/bench.sh`,
`scripts/bench-opt.sh`, `scripts/bench-randopt.sh`): means of two interleaved runs of main
784a24a (7278262 for the 1-thread optimal row) unless stated otherwise.

## Two-phase solving (time limit)

Mean best length, 40 random positions (48 threads) / 20 (1 thread):

| | 48 threads | 1 thread |
|---|---|---|
| first solution | 23.1 moves @ 5.7 ms | 23.1 moves @ 0.16 s |
| after 0.01 s | 22.31 (36 of 40 solved) | (none solved yet) |
| after 0.1 s | 21.75 | 22.67 (9 of 20 solved) |
| after 1 s | 21.50 | 22.30 |
| after 5 s | 21.10 | 22.00 |
| search speed | 1.71e9 nodes/s | 3.8e7 nodes/s |

Since the merged traversal (`opt/merged`; a node alive in both variants now counts once), final
main 97452e0 of 2026-11-10 (single run of the same suites): first solution 22.98 moves @ 6.0 ms
(48T) / 23.25 @ 0.17 s (1T); best after 0.1 s 21.75 / 22.70, after 1 s 21.525 / 22.35, after
5 s 21.175 / 22.00 (two of the 40 positions reach their 5 s best later than before);
1.52e9 / 3.54e7 nodes/s.

## Proven optimal solutions

| set | threads | time to prove (total) | mean optimal length |
|---|---|---|---|
| 10 x 19-move scrambles | 48 | 8.0 s (0.80 s per position, main 97452e0; 12.5 s before `opt/merged`) | 18.1 |
| 10 x 17-move scrambles | 1 | 2.9 s (0.29 s per position; target-length passes; 3.6 s before) | 16.5 |

Random positions (48 threads). "Proof": main a1f665d with the merged traversal (mean of two
runs); in brackets main 66ac3cb before it, and the morning of 2026-11-10 before the early proof
end and the sibling pipeline. The last column proves the same positions with
`FTO_FULL_PROOF=1` (no early end, no target passes) on main d970d72 (in brackets on 66ac3cb,
before the merged traversal): the proven lengths agree.

| seed | optimal | first solution | proof | full proof |
|---|---|---|---|---|
| 1 | 21 | 25 @ 12 ms | 54.0 s (64.8 s; 419 s) | 181.3 s (275.6 s) |
| 2 | 21 | 24 @ 3 ms | 61.8 s (64.6 s; 481 s) | 216.0 s (316.2 s) |
| 3 | 20 | 23 @ 6 ms | 4.96 s (5.2 s; 39 s) | 17.0 s (26.1 s) |

### With the distance-to-solved bound (`opt/tosolved`)

The large tier's 3.5 GB distance-to-solved table (optimal mode only; see
[PHASES.md](PHASES.md#4-search)), main c1dfbd3 against the branch, quiet machine, each run
alone under the bench lock, two interleaved rounds (random positions: one run each):

| | main | with the bound | ratio |
|---|---|---|---|
| 10 x 17-move scrambles, 1 thread (time per position) | 0.297 s | 0.233 s | 0.784x (nodes 0.743x) |
| the same, table read on first use (its reading inside the first solves) | 0.301 s | 0.237 s | 0.788x (nodes 0.777x) |
| 10 x 19-move scrambles, 48 threads | 0.833 s | 0.701 s | 0.841x (nodes 0.806x) |
| the same, table read on first use | 0.838 s | 0.715 s | 0.853x (nodes 0.812x) |
| random position, seed 1 (21) | 71.5 s | 64.2 s | 0.897x |
| random position, seed 2 (21) | 69.4 s | 58.2 s | 0.838x |
| random position, seed 3 (20) | 5.90 s | 4.83 s | 0.819x |
| time-limited search, nodes/s (48 / 1 threads) | | | 0.992x / 1.001x (unchanged: not used) |

The other rows read the table with the other tables. All proven lengths are identical. (These
main times are above the table above, measured on an earlier main and a quieter day.)

## Move sets (`--moves`)

Large tier; wide turns never shorten a solution, so only slices are measured.

| | face turns | with slices |
|---|---|---|
| optimal length, 10 × 19-move scrambles (mean) | 18.10 | 16.80 |
| time to prove optimality (48 threads, mean) | 1.29 s | 39.4 s |
| random positions, best after 5 s (48 threads) | 21.10 | 18.93 |
| random positions, best after 5 s (1 thread) | 22.00 | 19.80 |
| large-tier table generation (48 threads) | 30 s | 43 s |

(Measured before `opt/merged`; with it, slice proofs take 0.85x the time on the 17-move set at 48
threads and 0.67x on the 15-move set at 1 thread.) With slices the early end of the proof is off (see docs/PHASES.md section 8) and the branching
factor rises from 12.9 to about 18.3, so proofs take about 20 times as long as the full proof
without slices.

## Hard positions (`bench/hard.tsv`)

Positions supplied by the user: each is the position its listed algorithm solves (the
algorithm inverted, applied to the solved puzzle). Solved optimally in the normal FTO (centres of
a colour interchangeable), large tier, 48 threads. (A first run took the algorithms themselves as
scrambles, i.e. the inverse positions, which in the normal FTO need not have the same length:
21 / 21 / 21 / 22 there. Those numbers do not apply to the listed positions.)

| position | listed | proven optimum (normal FTO) | proof | optimal solver, large tier |
|---|---|---|---|---|
| tri_adj+edge_opp | 23 | 23 | 8,866 s (2.5 h, 1.4e13 nodes) | **23 in 266.5 s** (1.2e11 nodes) |
| tri_adj+edge_adj | 23 | 23 | 11,877 s (3.3 h, 1.9e13 nodes) | **23 in 803.0 s** (4.0e11 nodes) |
| tri_adj-[Uf_Fu;Db_Bd] | 23 | 23 | 10,942 s (3.0 h, 1.7e13 nodes) | **23 in 1,155.1 s** (5.7e11 nodes) |
| tri_adj+[UB_FD,NL_MR] | 23 | 23 | 8,397 s (2.3 h, 1.3e13 nodes) | **23 in 447.8 s** (2.1e11 nodes) |
| edge_hard | 21 | 21 | 74 s (64.3 s with the distance-to-solved bound) | 21 in 344.7 s |
| pure_corn_6flip | >= 24 (pure) | 20 | 661 s (`--lower-bound 19`, after an earlier run excluded 18); 117.2 s from scratch with the distance-to-solved bound | **20 in 0.48 s** |

The corner 6-flip undoes itself, so it is the same position either way. Its normal-FTO optimum is
20 (no solution of at most 19 moves after 620 s, a 20-move solution at 661 s; verified on the
sticker model); the solution exchanges centres of a colour, so it does not bear on the pure
count of at least 24. With edges and centres solved, the A-centre/edge table bounds it by 0 and
only the corner tables prune: each target pass costs about 9 times the previous one. The
distance-to-solved bound of optimal mode (`opt/tosolved`, large tier: corners, B-centres and an
edge triple together, to solved) proves it from scratch in 117.2 s at 48 threads (6.45e10
nodes; passes 18 and 19 complete at 8.2 s and 90.4 s), against 626.5 s (7.97e11 nodes) without
it, a 5.3x speed-up; edge_hard 75.6 -> 64.3 s (1.18e11 -> 9.41e10 nodes), same solutions.

All four positions listed at 23 are 23 in the normal FTO too (two-phase solver, optimal mode,
2026-11-10/11; every solution checked: scramble plus solution leaves each face one colour):

- tri_adj+edge_adj: BL' F BL' U F U L' D B L D' B' BL U' BR BL' F' R' L' R' L' D B
- tri_adj-[Uf_Fu;Db_Bd]: F L B' BL D' R L' D L F BL' F' L' B R' U BL U BR F' U BL' BR'
- tri_adj+[UB_FD,NL_MR]: B' BR BL L D' R' L' BR' U' BL BR' BL L' B' R' BR U' L R D' B D' L
- tri_adj+edge_opp: F BR' U' BL' B' U' F L' R B' D L' BR U BR F R' BR R BR' U F BL'

Each proof takes 2.3-3.3 hours on 48 threads (1.3-1.9e13 nodes), most of it in the last target
pass (no solution of at most 22 moves).

The optimal solver's large tier (`omelettehedron optimal solve --tier large`, 48 threads, quiet
machine, main 5106ccf) proves tri_adj+edge_opp 23 in **266.5 s** (1.23e11 nodes: depth 22
complete, no solution, at 152.8 s; a 23-move solution at 266.5 s), against 8,866 s for the
two-phase solver; solution `B D' L R' D L' B U L BL' F' BL D U BR B' R' BR L' B' U' BL D'`
(exchanges centres of a colour; scramble plus solution leaves each face one colour). It proves
the corner 6-flip 20 in 0.48 s (9.6e7 nodes), against 626.5 s for the two-phase solver without
a lower bound, tri_adj+edge_adj 23 in 803.0 s (depth 22 complete at 345.1 s; the two-phase
solver: 11,877 s; solution `B' D F' U' D R BR' D' U' F BL' D' L' F BL' D B' L BL U' L BL' U'`,
checked on the sticker model), and edge_hard 21 in 344.7 s (its corners are solved, so the triple table adds
little; the two-phase solver: 71-74 s). tri_adj+[UB_FD,NL_MR] 23 in 447.8 s (depth 22 complete at 439.3 s; solution `B D' R F D B' D' R'
BR' R L' B BL BR' B' BL F B' BL' BR L B' BL'`), and tri_adj-[Uf_Fu;Db_Bd] 23 in 1,155.1 s (depth 22
complete at 1,144.0 s; solution `BR U' BL L BL R BR B' U BL L' D B BL D' R' U BR' L' U' R' D R`;
light one-thread work ran beside this proof). All four 23-move positions in 44.5 minutes, against
11.1 hours for the two-phase solver; every solution checked on the sticker model. (Optimal
solver times on the quiet machine, 48 threads; the two-phase proofs ran on the same machine.)

## Start-up

Wall time of a solve of a 2-move scramble, tables already on disk (mostly loading the tables and
building the move and symmetry tables; `bench` leaves it out of its per-position times):

| tier | 1 thread | 48 threads |
|---|---|---|
| small | 0.41 s | 0.13 s |
| large | 1.38 s | 0.26 s |
| large with the distance-to-solved table (`opt/tosolved`), read at start-up / read on first use (main in the same sessions: 1.71 / 0.27 s, 1.72 / 0.28 s) | 2.13 / 1.72 s | 0.34 / 0.28 s |

## Tables

| tier | memory while solving | generation (48 threads) | generation (1 thread) | peak memory generating |
|---|---|---|---|---|
| small | 134 MB | 1.2 s | 27 s | 0.38 GB |
| medium | 1.04 GB | 6.6 s | 167 s | 1.26 GB |
| large | 11.6 GB (8.1 GB before `opt/tosolved`) | 30 s + 20 s (64 s from scratch beside a build) | — | 18.3 GB |
| huge | 33.5 GB | 3 min 14 s (incl. writing 30 GB) + about 20 s (not measured) | — | 29.4 GB before `opt/tosolved` (not measured since) |

The huge tier holds the same A-centre/edge values as medium with a cheaper index and is measured
slower than large (0.85x nodes/s), so large is the fastest tier. Tables are generated on first use
and loaded from disk afterwards; all tiers include the 4.7 MB symmetry-reduced corner/B-centre
table and the in-memory phase-2 table.

## WebAssembly (small and medium tiers)

Measured in Bun (JavaScriptCore) on the same machine; see the README's WebAssembly section.

| | native small tier | wasm, 1 thread | wasm, 48 threads |
|---|---|---|---|
| table generation | 27 s (1T), 1.2 s (48T) | 32.7 s | 16.4 s |
| search speed | 58 M nodes/s (1T) | 23 M nodes/s | 74 M nodes/s |
| random positions with any solution after 5 s | 9 of 10 (1T) | 4 of 10 | 10 of 10 |

The modules are 0.63 MB (one thread) and 0.67 MB (each threads module). Threads scale poorly
in JavaScriptCore because its WebAssembly atomic loads are slow under contention.

The medium tier (the demo's Faster choice), measured on 2026-11-11 while proofs kept 43 of the
48 cores busy (lowest priority; for one thread the process CPU time stands in for the time on
an idle machine):

| | small | medium |
|---|---|---|
| generation, wasm, 16 threads (wall time under load) | 257 s | 872 s (files identical to the native ones) |
| table load, wasm, 1 thread (CPU time) | 0.8 s | 1.6 s |
| memory of the one-thread module after loading | 358 MiB | 1220 MiB |
| 10 optimal solves of 16-move scrambles, wasm, 1 thread | 53.2 M nodes, 7.3 s CPU | 5.9 M nodes, 1.6 s CPU |
| 10 optimal solves of 17-move scrambles, native, 1 thread | 37.2 s CPU (1.1e9 nodes) | 3.3 s CPU (1.0e8 nodes) |

With the medium tier the module needs about 9-11 times fewer nodes; a node costs about twice
as much in WebAssembly (natively the same), so proofs take about a fifth of the time there.

## Optimal solver (v0)

The separate optimal solver (`omelettehedron optimal`, single-phase IDA* to the solved state;
[OPTIMAL.md](OPTIMAL.md)). Node counts only (one thread unless noted; the machine was busy, so
no timings yet). Same optimal lengths as the two-phase solver's optimal mode.

| set | tier | mean optimal length | mean nodes per position |
|---|---|---|---|
| 10 x 15-move scrambles, seed 1717 | small | 14.5 | 8.73e7 |
| 10 x 15-move scrambles, seed 1717 | medium | 14.5 | 2.87e6 |
| 10 x 17-move scrambles, seed 1717 | medium | 16.5 | 3.85e8 |
| 10 x 19-move scrambles, seed 1919 (48 threads) | medium | 18.1 | 3.76e10 |
| 10 x 15-move scrambles, seed 1717, with the two-phase AE table (`optimal/ae`) | small | 14.5 | 1.35e6 |
| 10 x 15-move scrambles, seed 1717, with the two-phase AE table | medium | 14.5 | 1.16e5 |
| 10 x 17-move scrambles, seed 1717, with the two-phase AE table | medium | 16.5 | 1.48e7 |

Hard positions (`bench/hard.tsv`; 48 threads beside the two-phase proofs, medium tier,
symmetric start): `pure_corn_6flip` proven optimal at 20 in 232 s (3.05e10 nodes; the
two-phase solver: 661 s with `--lower-bound 19`). The others reach lower bounds 19-20 in 300 s
(details in [OPTIMAL.md](OPTIMAL.md#symmetric-start-on-the-hard-positions)).
