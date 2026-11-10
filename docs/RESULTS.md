# Results

[Documentation home](../README.md) · [Optimisation log](OPTIMIZATION_LOG.md) ·
[Phase design](PHASES.md)

Measured on a 48-core AMD EPYC VM (94 GB RAM, no hardware performance counters), large tier
(A-centre/edge 1.03 GB + corner/B-centre + joint 7.0 GB tables), tables loaded before timing
(loading takes 0.33 s with 48 threads). Positions are seeded random states unless stated
otherwise; the numbers below are from the standard suites (`scripts/bench.sh`,
`scripts/bench-opt.sh`, `scripts/bench-randopt.sh`): means of two interleaved runs of main
5a2deed (14551aa for the 1-thread optimal row) unless stated otherwise.

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
main b0da038 of 2026-11-10 (single run of the same suites): first solution 22.98 moves @ 6.0 ms
(48T) / 23.25 @ 0.17 s (1T); best after 0.1 s 21.75 / 22.70, after 1 s 21.525 / 22.35, after
5 s 21.175 / 22.00 (two of the 40 positions reach their 5 s best later than before);
1.52e9 / 3.54e7 nodes/s.

## Proven optimal solutions

| set | threads | time to prove (total) | mean optimal length |
|---|---|---|---|
| 10 x 19-move scrambles | 48 | 8.0 s (0.80 s per position, main b0da038; 12.5 s before `opt/merged`) | 18.1 |
| 10 x 17-move scrambles | 1 | 2.9 s (0.29 s per position; target-length passes; 3.6 s before) | 16.5 |

Random positions (48 threads). "Proof": main ddcb5ba with the merged traversal (mean of two
runs); in brackets main 496b389 before it, and the morning of 2026-11-10 before the early proof
end and the sibling pipeline. The last column proves the same positions with
`FTO_FULL_PROOF=1` (no early end, no target passes) on main ef10b90 (in brackets on 496b389,
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

Positions supplied by the user, solved optimally in the normal FTO (centres of a colour
interchangeable), large tier, 48 threads, main b677856. Every solution was checked on the sticker
model; the `edge_hard` scramble gives an edges-only position here, as listed, so the notation
agrees with the source's.

| position | listed | proven optimum | proof | solves the pure FTO too |
|---|---|---|---|---|
| tri_adj+edge_opp | 23 | 21 | 76 s | no |
| tri_adj+edge_adj | 23 | 21 | 122 s | no |
| tri_adj-[Uf_Fu;Db_Bd] | 23 | 21 | 76 s | no |
| tri_adj+[UB_FD,NL_MR] | 23 | 22 | 654 s | no |
| edge_hard | 21 | 21 | 74 s | no (a different 21 does, per the source) |
| pure_corn_6flip | >= 24 (pure) | 18-23 (normal FTO) | stopped | no |

The corner 6-flip leaves edges and centres solved, so the A-centre/edge table bounds it by 0 and
only the corner tables prune: each target pass costs about 9 times the previous one (no
solution of at most 16 moves after 14.8 s, at most 17 after 128.7 s), so excluding 22 moves would
take months at this rate. Its normal-FTO optimum lies between 18 and 23 (a 23-move solution is
found in 1 s); the pure optimum is at least that.

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

## WebAssembly (small tier)

Measured in Bun (JavaScriptCore) on the same machine; see the README's WebAssembly section.

| | native small tier | wasm, 1 thread | wasm, 48 threads |
|---|---|---|---|
| table generation | 27 s (1T), 1.2 s (48T) | 32.7 s | 16.4 s |
| search speed | 58 M nodes/s (1T) | 23 M nodes/s | 74 M nodes/s |
| random positions with any solution after 5 s | 9 of 10 (1T) | 4 of 10 | 10 of 10 |

The modules are 0.55 MB each. Threads scale poorly in JavaScriptCore because its WebAssembly
atomic loads are slow under contention.
