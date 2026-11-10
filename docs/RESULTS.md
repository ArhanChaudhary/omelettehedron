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

## Proven optimal solutions

| set | threads | time to prove (total) | mean optimal length |
|---|---|---|---|
| 10 x 19-move scrambles | 48 | 12.5 s (1.25 s per position) | 18.1 |
| 10 x 17-move scrambles | 1 | 3.6 s (0.36 s per position; target-length passes) | 16.5 |

Random positions (48 threads, single runs; in brackets the same proof before the early proof end
and the sibling pipeline). The last column proves the same positions with `FTO_FULL_PROOF=1`
(no early end, no target passes), main 496b389: the proven lengths agree.

| seed | optimal | first solution | proof | full proof |
|---|---|---|---|---|
| 1 | 21 | 25 @ 12 ms | 64.8 s (419 s on 2026-11-10 morning) | 275.6 s |
| 2 | 21 | 24 @ 3 ms | 64.6 s (481 s) | 316.2 s |
| 3 | 20 | 23 @ 6 ms | 5.2 s (39 s) | 26.1 s |

## Move sets (`--moves`)

Large tier; wide turns never shorten a solution, so only slices are measured.

| | face turns | with slices |
|---|---|---|
| optimal length, 10 × 19-move scrambles (mean) | 18.10 | 16.80 |
| time to prove optimality (48 threads, mean) | 1.29 s | 39.4 s |
| random positions, best after 5 s (48 threads) | 21.10 | 18.93 |
| random positions, best after 5 s (1 thread) | 22.00 | 19.80 |
| large-tier table generation (48 threads) | 30 s | 43 s |

With slices the early end of the proof is off (see docs/PHASES.md section 8) and the branching
factor rises from 12.9 to about 18.3, so proofs take about 20 times as long as the full proof
without slices.

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
