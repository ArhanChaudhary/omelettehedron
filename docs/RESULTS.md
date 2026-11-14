# Results

[Documentation home](../README.md) · [Optimisation log](OPTIMIZATION_LOG.md) ·
[Phase design](PHASES.md)

The reference run of 2026-11-11 comes first; the sections after it are earlier measurements, each
labelled with its commit (the previous reference run: main 97452e0 of 2026-11-10).

Measured on a 48-core AMD EPYC VM (94 GB RAM, no hardware performance counters), large tier
(A-centre/edge 1.03 GB + corner/B-centre + joint 7.0 GB tables), tables loaded before timing
(loading takes 0.33 s with 48 threads). Positions are seeded random states unless stated
otherwise; the numbers below are from the standard suites (`scripts/bench.sh`,
`scripts/bench-opt.sh`, `scripts/bench-randopt.sh`): means of two interleaved runs of main
784a24a (7278262 for the 1-thread optimal row) unless stated otherwise.

## Reference run, 2026-11-11 (main 9a9262f, tag v0.1.0)

Quiet machine, each item alone under the bench lock, tables loaded before timing; large tier for
the two-phase solver unless stated otherwise. From now on comparisons use the quick set
(`scripts/bench-quick.sh`, under 5 minutes; see the README), so its numbers are the baseline.

### What changed since the previous reference run

"Now" is main 9a9262f on the same sets as "before" (means of two runs, one run of each single
position; see the last table of this section for the full suites).

| | before | now | change | merged on 2026-11-11 |
|---|---|---|---|---|
| two-phase, optimal mode, 10 x 19 moves, 48 threads (per position) | 0.80 s (main 97452e0) | 0.667 s | 0.83x | `opt/tosolved`, `opt/tosolved-order` |
| two-phase, optimal mode, 10 x 17 moves, 1 thread (per position) | 0.286 s | 0.227 s | 0.79x | the same |
| two-phase, edge_hard (21), 48 threads | 74 s | 62.8 s | 0.85x | `opt/tosolved` |
| two-phase, corner 6-flip (20) from scratch, 48 threads | 626.5 s | 117.2 s (measured on the branch, not rerun) | 0.19x | `opt/tosolved` |
| two-phase, time limit 5 s, 48 threads: first solution, best after 5 s | 22.98 @ 6.0 ms, 21.175 | 22.98 @ 5.8 ms, 21.19 | unchanged (nodes/s 1.03x) | (the new bound is used in optimal mode only) |
| two-phase, time limit 5 s, 1 thread: first solution, best after 5 s | 23.25 @ 0.168 s, 22.00 | 23.25 @ 0.167 s, 22.00 | unchanged (nodes/s 1.02x) | |
| optimal solver, 10 x 17 moves, 1 thread (per position) | 1.03 s (first timing: medium, main 3deb77e) | medium 0.585 s, large 0.482 s | 0.57x / 0.47x | `optimal/order`, `optimal/cbe-solved` (the large tier), `optimal/r3-prefetch`, `optimal/r3-masks`, `optimal/r3-aepre`, `optimal/large-noce` |
| optimal solver, 10 x 19 moves, 48 threads (per position) | 4.12 s (medium, main 3deb77e) | medium 2.37 s, large 0.98 s | 0.58x / 0.24x | the same |
| optimal solver, large tier, start of round 3 (main 4c96995) -> now: 17 / 19-move sets | 0.758 / 1.528 s | 0.482 / 0.982 s | 0.64x / 0.64x | round 3, `optimal/large-noce` |
| optimal solver, random position seed 1 (21), large, 48 threads | 91.7 s (main 4c96995) | 66.4 s | 0.72x | the same |
| optimal solver, random positions seeds 2 (21) / 3 (20), large, 48 threads | | 143.4 s (9.16e10 nodes) / 73.3 s (measured later, main 47ca15d) | | (the two-phase solver under `-t 1200`: 48.6 / 4.0 s) |
| optimal solver, edge_hard (21), large, 48 threads | 344.7 s (main 5106ccf) | 226.5 s (main 47ca15d) | 0.66x | round 3, `optimal/large-noce` (the two-phase solver: 56.4 s) |
| optimal solver, corner 6-flip (20), 48 threads | 232 s (medium, v0, beside other work) | large 0.34 s, medium 21.6 s | | the large tier; round 3 |
| optimal solver, the four 23-move positions of `bench/hard.tsv` | (two-phase solver: 11.1 h) | 44.5 min (large, measured earlier today, not rerun) | 0.067x | |

### Quick set (`scripts/bench-quick.sh --opt-tier medium,large --rounds 2`)

Two rounds (each figure the mean of the two; the rounds differ by up to 4%). Two-phase optimal
mode with the distance-to-solved table read at start-up (`FTO_CBS_EAGER=1`).

| item | threads | result |
|---|---|---|
| two-phase, optimal mode, 10 x 17 moves (seed 1717) | 1 | 2.24 s (0.224 s per position), 7.86e7 nodes, mean 16.50 |
| two-phase, optimal mode, 10 x 19 moves (seed 1919) | 48 | 7.02 s (0.702 s per position), 1.00e10 nodes, mean 18.10 |
| two-phase, time limit 1 s, 10 random positions (seed 2026) | 48 | best 21.80, first solution 23.00 @ 6.9 ms, 1.45e9 nodes/s |
| two-phase, time limit 1 s, 10 random positions (seed 7) | 1 | best 22.20, first solution 23.10 @ 0.211 s, 3.56e7 nodes/s |
| optimal solver, medium, 10 x 17 moves | 1 | 5.88 s, 1.48e8 nodes |
| optimal solver, medium, 10 x 19 moves | 48 | 24.4 s, 2.52e10 nodes |
| optimal solver, large, 10 x 17 moves | 1 | 5.02 s, 7.82e7 nodes |
| optimal solver, large, 10 x 19 moves | 48 | 9.99 s, 6.28e9 nodes |

Single runs (48 threads): the two-phase solver proves edge_hard 21 in 62.8 s (9.43e10 nodes);
the optimal solver proves the corner 6-flip 20 in 0.34 s with the large tier (9.6e7 nodes) and
in 21.6 s with medium (2.06e10 nodes).

### The previous suites, once more

Before the switch to the quick set, the suites of the previous reference run ran once more on
the same main (two runs of `scripts/bench.sh` and `bench-opt.sh`, one of `bench-randopt.sh`),
so they line up with the numbers below:

| | main 97452e0 (2026-11-10) | main 9a9262f |
|---|---|---|
| 48 threads, 40 positions, 5 s: first solution | 22.975 @ 6.0 ms | 22.975 @ 5.8 ms |
| after 0.01 s / 0.1 s / 1 s / 5 s | 22.30 (33 of 40 solved) / 21.75 / 21.525 / 21.175 | 22.36 (36 of 40) / 21.75 / 21.50 / 21.19 |
| search speed | 1.52e9 nodes/s | 1.56e9 nodes/s |
| 1 thread, 20 positions, 5 s: first solution | 23.25 @ 0.168 s | 23.25 @ 0.167 s |
| after 0.1 s / 1 s / 5 s | 22.70 (10 of 20) / 22.35 / 22.00 | 22.70 (10 of 20) / 22.325 / 22.00 |
| search speed | 3.54e7 nodes/s | 3.62e7 nodes/s |
| optimal mode, 10 x 19 moves, 48 threads | 8.0 s | 6.68 s (table read on first use) |
| optimal mode, 10 x 17 moves, 1 thread | 2.86 s | 2.27 s |
| random positions 1 / 2 / 3, `-t 1200` | 54.6 / 64.3 / 5.3 s | 53.2 / 60.5 / 4.9 s |

`bench-randopt.sh` caps each position with `-t 1200`, and with a time limit the search is the
time-limited one, which on this main did not use the distance-to-solved bound: the random
positions, and edge_hard (67.9 s) and the corner 6-flip (624 s) under `-t 1800`, show the speed
without it. Without `-t` (optimal mode) edge_hard takes 62.8 s and the 6-flip 117 s. Since
`opt/tl-cbs` time-limited runs use it too; see [Under a time limit](#under-a-time-limit-opttl-cbs).

## The optimal solver's huge tier (2026-11-13, `optimal/ea-huge`; not in the reference run)

The `huge` tier (27.9 GB: an edge x A-centre table to solved with 288 symmetries on both sides,
and corners x centres; [OPTIMAL.md](OPTIMAL.md#optimisation-round-4-what-the-memory-buys-the-huge-tier-optimalea-huge))
against the large tier (22.0 GB), the same binary, quiet machine, bench lock, two rounds:

| set | threads | large | huge | huge / large |
|---|---|---|---|---|
| 10 x 17 moves (seed 1717) | 1 | 4.83 s, 7.82e7 nodes | 0.74, 0.78 s, 1.79e7 nodes | 0.153x, 0.160x time, 0.229x nodes |
| 10 x 19 moves (seed 1919) | 48 | 9.79, 9.82 s, 6.28e9 nodes | 2.17, 2.20 s, 2.17e9 nodes | 0.221x, 0.223x time, 0.345x nodes |

Random positions (48 threads, one position per seed as in `bench-randopt.sh`, proven optimal):
huge 18.3 s (seed 1, 21 moves, 1.95e10 nodes), 27.5 s (seed 2, 21 moves, 2.91e10 nodes), 3.2 s
(seed 3, 20 moves, 3.33e9 nodes); the large tier took 66.4 / 143.4 / 73.3 s and the two-phase
solver's optimal mode about 53 / 49 / 4 s (earlier mains, see above). With the huge tier the
optimal solver proves these positions and the two sets faster than the two-phase solver's
optimal mode (17-move set 0.75 s against 2.24 s on one thread, 19-move set 2.2 s against 7.0 s on
48).

## Progress and both solvers at once

Measured on 2026-11-11 (`feat/progress-combined`, 48-core EPYC, large tiers, under the bench
lock).

### `--progress`

Without `--progress`, against main bbb1c7d (`scripts/bench-quick.sh --compare`, 2 interleaved
rounds; one core was busy with another session's job, the same for both builds):

| item | B / A (medians) | per round | B better | nodes (1 thread) |
|---|---|---|---|---|
| two-phase, optimal mode, 10 x 17 moves, 1 thread | time 1.002x | 1.009, 0.996 | 1/2 | same |
| two-phase, optimal mode, 10 x 19 moves, 48 threads | time 0.967x | 0.996, 0.940 | 2/2 | |
| two-phase, 1 s, 10 random positions, 48 threads | nodes/s 1.010x | 1.007, 1.014 | 2/2 | |
| two-phase, 1 s, 10 random positions, 1 thread | nodes/s 1.026x | 1.009, 1.043 | 2/2 | |
| optimal solver, 10 x 17 moves, 1 thread | time 1.007x | 0.995, 1.018 | 1/2 | same |
| optimal solver, 10 x 19 moves, 48 threads | time 0.986x | 0.971, 1.001 | 1/2 | |

With `--progress` against without (the solvers' own times, 3 interleaved rounds; on / off per
round):

| set | before the fix below | with it |
|---|---|---|
| two-phase, 15 scrambles of 14-19 moves, 48 threads | 1.034, 1.012, 0.999 | 1.000, 1.004, 0.991 |
| two-phase, random position seed 3 (20), 48 threads | 1.001, 0.984, 0.987 | 0.991, 1.001, 0.981 |
| two-phase, 12 scrambles of 14-17 moves, 1 thread (24 ms each) | 1.034, 1.055, 1.048 | 1.025, 1.017, 1.010 |
| optimal solver, the 15 scrambles, 48 threads | 1.000, 1.049, 0.984 | 1.009, 1.005, 0.986 |

The heartbeat thread was first started for every iteration, and a solve has many short ones;
now an iteration that follows one shorter than 1/64 of the heartbeat interval gets none. What
remains on the short one-thread solves (about 0.4 ms per solve) is printing their 50 or so
lines. One-thread node counts are identical with and without progress lines (test
`progress_does_not_change_the_search`).

### `combined`

Seconds on 48 threads in all (combined: 24 for each solver), each run the solvers' own time
(tables loaded before), one run each. The 17- to 19-move scrambles are seeded random face-turn
sequences; the random position is `--random 1 --seed 3`. The edge_hard runs shared the machine
with one busy core of another session.

| position | optimal | two-phase alone (optimal mode) | optimal solver alone | combined, `improve` | combined, `passes` |
|---|---|---|---|---|---|
| `BL BR U F' R' D R' F' B' U BR' L BR U F BR' D'` | 15 | 0.002 | 0.002 | 0.002 | 0.002 |
| `U' B' U L BL' F' BL L' F' BR L B' L B' U' L F' BL'` | 18 | 0.109 | 0.521 | 0.274 (solution 0.272, lower bound 0.116) | 0.222 |
| `L' F L B' BL B' BR' U BL' B F' BR' D' F L D' F BR' L` | 19 | 1.564 | 2.970 | 1.397 | 1.348 |
| random position seed 3 | 20 | 4.655 | 73.6 | 8.47 (solution 1.04) | 9.11 (solution 9.11) |
| edge_hard | 21 | 62.8 (above) | 344.7 (hard positions) | 112.8 (solution 0.54) | 123.3 (solution 123.3) |
| corner 6-flip | 20 | 117.2 (above) | 0.337 | 0.632 | 0.698 |

Each solver has half of the threads, so the combined solve takes about twice the time of the
faster solver alone (1.8x on the random position and edge_hard, 1.9x on the 6-flip), but is
far faster than the slower one; on the 19-move scramble it beat both. The optimal solver's
help shows where the two-phase solver finds the optimum early: the optimal solver then needs
only the depths below it (on the random position, alone, its depth 19 ended after about 6 s
of its 74 s). `improve` and `passes` prove optimality about equally fast, and `improve` has the
optimal solution almost at once, while `passes` finds a first solution only when it reaches
the optimal length, so `improve` is the default.

Correctness: 12 seeded scrambles of 14-17 moves (optimal 12-17) give the same length with the
two-phase solver, the optimal solver and both combined modes (`improve` 24 + 24 threads,
`passes` 8 + 8); edge_hard and the corner 6-flip give 21 and 20, as each solver alone. The
tests (`tests/combined.rs`, small tiers) compare the combined lengths with both solvers' on
scrambles of up to 15 moves with even and uneven thread splits, check the bounds under time
limits, and solve with a lower bound given from outside, also raised during the search.

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

The large tier's 3.5 GB distance-to-solved table (optimal mode only until `opt/tl-cbs`; see
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

### Under a time limit (`opt/tl-cbs`)

Until `opt/tl-cbs` a time limit (`-t`) or `-l` turned the bound off, although it prunes against
the best length found, which is just as admissible as a target pass. Main 1cce546 against the
branch, quiet machine, each run alone under the bench lock (wall times including the table
reading; the hard positions and random positions one run each):

| | main | with the bound | ratio |
|---|---|---|---|
| corner 6-flip (20), `-t 1800`, 48 threads | 624 s (on 9a9262f) | 84.0 s | 0.135x |
| edge_hard (21), `-t 1800`, 48 threads | 67.0 s | 56.4 s | 0.84x |
| random position seed 1 (21), `-t 1200` | 53.5 s (8.61e10 nodes) | 45.5 s (7.27e10) | 0.85x |
| random position seed 2 (21), `-t 1200` | 60.8 s (9.69e10) | 48.6 s (7.66e10) | 0.80x |
| random position seed 3 (20), `-t 1200` | 5.0 s (7.88e9) | 4.0 s (6.26e9) | 0.80x |
| 20 random positions, 3 s, 48 threads: mean best after 0.1 / 0.3 / 1 / 3 s (2 rounds) | 21.95 / 21.85 / 21.65 / 21.50 | 21.95 / 21.80-21.85 / 21.65 / 21.50 | unchanged; nodes/s 0.994x, 0.967x |
| 8 random positions, 3 s, 1 thread (2 rounds) | 22.43 (0.3 s) / 21.875-22.00 (1 s) / 21.50 (3 s) | 22.43 / 21.875 / 21.50 | unchanged; nodes/s 1.013x, 1.003x |

Under a time limit the 6-flip is now faster than in optimal mode (117 s): the first solution
found is already 20, and pruning against it amounts to pass 19 alone, while optimal mode also
runs passes 18 and 20. Memory: the table adds 3.5 GB to a long time-limited run (peak 11.4
against 8.0 GB); short runs never read it.

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

## Tiny tier

The `tiny` tier (`src/light.rs`, the demo's Light) is for phones: the user reported that the
134 MB of the small tier took too long to prepare on an iPhone.

### Choosing the tables (estimate before building)

Projections of the A-centre/edge coordinate merge labels of the A-centres (A) or of the edge
triples (E): all four told apart (Full, 369,600 labellings), two told apart and two merged
(Pair, 18,480), two pairs (Split, 924) or one against three (Single, 220). Each projection's
exact distance to the projection of `H_A` was computed (a throwaway program, not kept), read
on seeded random positions through relabellings (4M for the 1.49x-1.99x rows and the codes
below, 100k-1M for the others), and the maximum with the corner/B-centre table compared with
the small tier by Korf's formula (relative nodes = sum over k of 12.93^-k P(h <= k)):

| bound | entries | nodes vs small (estimate) |
|---|---|---|
| small tier (Full x Full, 288 symmetries) | 516M | 1.00x |
| A Full x E Pair (48 symmetries), 6 lookups | 145M | 1.49x |
| the same, 4 / 3 / 2 / 1 lookups | | 1.70x / 1.99x / 2.90x / 16x |
| A Pair x E Full, 6 lookups | 145M | 1.48x |
| A Pair x E Pair (6), Full x Split (3), Split x Full (3), Full x Single (4), Single x Full (4), all together | 17M | 9.5x |
| A Full x E Split, 3 lookups | 3.6M | 24x |

Nothing lies between: the next smaller projections have about 7M entries. A maximum of
lookups cannot be decoded from distances mod 3 (that needs each lookup's parent value), so the
tables hold direct 2-bit codes: P (Full x Pair) 0 = at most 7, 1 = 8, 2 = 9, 3 = 10 or more,
which alone is 81x (code 0 says nothing); with S (Full x Split, codes at most 4 / 5 / 6 / 7 or
more) read where every P lookup says 0, 1.50x (1.72x with 4 P lookups). Without the
corner/B-centre table the estimate is far worse (2.2x even with exact S values, and the leaves
of the search would not be found), so the tier keeps it.

### Measured (2026-11-11, the 48-core VM; one thread unless noted)

The same 10 random positions (seed 2026) and 10 random 15-move scrambles as `bench`; the
WebAssembly figures in Bun 1.4.4 (JavaScriptCore), loaded from the native files.

| | small | tiny |
|---|---|---|
| tables (A-centre/edge + corner/B-centre) | 129.0 + 4.7 MB | 37.3 + 4.7 MB |
| generation, native, 1 / 48 threads | 21-27 s / 1.1-1.2 s | 5.0 s / 0.28 s (peak 0.20 GB) |
| generation, WebAssembly, 1 / 8 threads | 29.6-30.8 s / 5.7 s | 6.8 s / 1.4 s |
| table load, WebAssembly | 0.47 s | 0.29 s |
| one-thread module after loading / after generating | 314 / 352 MiB (352: before the edge symmetry tables were made in place) | 190 / 216 MiB |
| search speed, native / WebAssembly | 57.5 / 20.9 M nodes/s | 31.6 / 13.4 M nodes/s |
| random positions with a solution within 1 s / 5 s, native | 3 / 9 of 10 | 0 / 5 of 10 |
| first solution (mean over those found), native | 2.00 s, 23.6 moves | 2.79 s, 22.8 moves |
| random positions with a solution within 5 s, WebAssembly | 5 of 10 (first after 2.63 s) | 2 of 10 (first after 4.17 s) |
| optimal, 15-move scrambles, native (mean per scramble) | 7.5 ms | 42.6 ms (about 2.8x the nodes) |
| optimal, the demo test's 14-move scramble | 169,652 nodes | 315,594 nodes |

Per node the tiny tier costs about twice as much (natively 0.55x the speed, in WebAssembly
0.64x): a child reads up to six P words instead of one A-centre/edge word, and at most 7 moves
left, three S words. Its tables are read in stages to keep that low (for children the first P
lookup does not cut, the other five; S first where it can cut, P not at all below 7, where it cut
0.2% of what S kept). The first-solution times of random positions on one thread are long for
both tiers; the tiny tier finds about half as many within a few seconds. Optimal solves of short
scrambles need more nodes than the estimate for random positions suggests: there the A-centre
and edge distances are small, where the codes are coarse.

K, the P lookups per position: with six instead of four (one thread, the same positions with 10 s
each), the first solution came in 0.94-0.98x the time on each of the 6 positions both solved,
and 8 of 10 positions were solved instead of 6, at 31.4 against 33.4 M nodes/s; so six.

Existing tiers: the search's hooks for the tiny tier sit behind the constant
`AeIndex::DIRECT`, so its hot functions compile for the other indices byte for byte as before
(`expand_any`, `stage2_any`, `finish`; objdump, addresses aside), and one-thread node counts are
unchanged. The quick set still came out 0.1-0.9% slower on the 48-thread 19-move optimal set (3
rounds: 1.005x, 1.009x, 1.001x; the 17-move set on one thread 1.011x, 0.994x, 1.005x), about the
same as before the gate (0.1-1.4%), so it is not the search code; the binary around it changed.

## Optimal solver (v0)

The separate optimal solver (`omelettehedron solve optimal`, single-phase IDA* to the solved state;
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
