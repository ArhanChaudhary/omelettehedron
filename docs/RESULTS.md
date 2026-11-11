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

| position | listed | proven optimum (normal FTO) | proof |
|---|---|---|---|
| tri_adj+edge_opp | 23 | (running) | |
| tri_adj+edge_adj | 23 | 23 | 11,877 s (3.3 h, 1.9e13 nodes) |
| tri_adj-[Uf_Fu;Db_Bd] | 23 | (running) | |
| tri_adj+[UB_FD,NL_MR] | 23 | (running) | |
| edge_hard | 21 | 21 | 74 s |
| pure_corn_6flip | >= 24 (pure) | 20 | 661 s (`--lower-bound 19`, after an earlier run excluded 18) |

The corner 6-flip undoes itself, so it is the same position either way. Its normal-FTO optimum is
20 (no solution of at most 19 moves after 620 s, a 20-move solution at 661 s; verified on the
sticker model); the solution exchanges centres of a colour, so it does not bear on the pure
count of at least 24. With edges and centres solved, the A-centre/edge table bounds it by 0 and
only the corner tables prune: each target pass costs about 9 times the previous one.

## Start-up

Wall time of a solve of a 2-move scramble, tables already on disk (mostly loading the tables and
building the move and symmetry tables; `bench` leaves it out of its per-position times):

| tier | 1 thread | 48 threads |
|---|---|---|
| small | 0.41 s | 0.13 s |
| large | 1.38 s | 0.26 s |

## Tables

| tier | memory while solving | generation (48 threads) | generation (1 thread) | peak memory generating |
|---|---|---|---|---|
| small | 134 MB | 1.2 s | 27 s | 0.38 GB |
| medium | 1.04 GB | 6.6 s | 167 s | 1.26 GB |
| large | 8.1 GB | 30 s | — | 8.1 GB |
| huge | 30 GB | 3 min 14 s (incl. writing 30 GB) | — | 29.4 GB |

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
