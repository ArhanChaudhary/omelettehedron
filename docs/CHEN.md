# Chen Shuang's block split

[Documentation home](../README.md) · [Two-phase design](PHASES.md) · [Results](RESULTS.md)

On 2026-11-12 Chen Shuang proposed a new two-phase split for the FTO: phase 1 solves the
tetraminx centres, and phase 2 solves the rest with moves that keep each face's six
tetraminx-centre pieces in at most two same-coloured blocks of three, much like Square-1's
constraint. We implemented it as a second solver (`--method chen`) and measured it against
omelettehedron's two-phase solver. It is not in the released solver; this page explains why.

The code is on the branch [`neg/chen-twophase`](https://github.com/ArhanChaudhary/omelettehedron/tree/neg/chen-twophase)
(`src/chen/`, with a longer technical write-up in its `docs/CHEN.md`).

## In short

* **The numbers check out.** Our model reproduces the 460 × 3⁴ = 37,260 block configurations,
  the phase-2 space of 1.59 × 10¹⁴ positions (the "158T"), and the phase-2 state counts for
  depths 0 to 5.
* **Faster first solution, but a longer one.** With the smaller table sets it finds a first
  solution up to 7 times sooner. That solution is 3 to 5 moves longer.
* **Longer solutions after the first moment.** From 0.1 s on 48 threads, or 1 s on one
  thread, the two-phase solver's mean solution is shorter at every table size, though at those
  points it has not always solved as many positions (tiny tier on 48 threads at 0.1 s: 5 of 10;
  tiny and small tiers on one thread at 1 s: 3 and 5 of 10). At 10 s it is 0.2 to 1.5 moves
  shorter.
* **More memory.** It needs 158 MB of extra tables.
* **The cause.** Phase 2 is a search over 1.59 × 10¹⁴ positions. The largest tables that fit
  estimate about 11 moves for phase 2 where about 16 are needed, so the search finds little in
  the nodes it can afford. The two-phase solver's phase 2 is read off an exact table, with no
  search.

## The split

The faces of the FTO split into two classes of four, and no two faces of a class share an
edge. In our labelling, class A is U, F, BR, BL and class B is R, L, B, D. Each class-A face
has a *hexagon*: its 3 centres and the class-A stickers of its 3 edges. Together the four
hexagons are the "tetraminx centres".

Chen Shuang's description has the restricted moves on U, F, BL, BR. With the hexagons on the
class-A faces as here, the restricted moves are R, L, B and D instead. The two are the same
split, up to symmetry.

* **Phase 1** solves the four hexagons: each holds its own 3 centres and 3 edges, possibly
  turned by a multiple of 120°. omelettehedron already has a table for this (its A-centre/edge
  coordinate, exact from the medium tier up). It needs at most 14 moves, 10.9 on average,
  matching the 14 moves quoted from the FTO group.
* **Phase 2** solves the rest with constrained moves. U, F, BR and BL turns are always allowed.
  A turn of R, L, B or D cuts three hexagons into 3 stickers inside the turning layer and 3
  outside. It is allowed only if, on each of the three, the 3 stickers inside are one colour
  and the 3 outside are one colour.

**What we reproduced.** From solved, the allowed moves reach 37,260 hexagon configurations: 460
colour shapes with 3⁴ = 81 block turns each, the "460 × 3⁴ cases". Each case has 12!/6⁴ ×
11,520 positions (369,600 B-centre colourings × 11,520 corner states). That makes
37,260 × 11,520 × 369,600 = 1.59 × 10¹⁴ positions in phase 2. A breadth-first search in our
Python model gives 1, 16, 160, 1,344, 10,464 and 78,152 positions at depths 0 to 5, the same as
the twsearch distribution in the proposal.

## How it was implemented

The block-split solver reuses everything except phase 2. Phase 1 uses the tier's own
A-centre/edge table and runs IDA* in both face classes, as the two-phase solver does. Each
phase-1 solution then gets a phase-2 IDA* with the allowed moves.

A phase-2 path ends once it reaches the positions solvable with U, F, BR and BL turns alone
(933,120 of them; the two-phase solver keeps an exact table of these), and that table gives the
rest. Every allowed path to solved ends this way, so this finds the shortest allowed
completion. Phase 2 is pruned by the largest of three tables, the best combinations of two
parts that fit in memory:

| Table | Size | Max | Mean |
|---|---:|---:|---:|
| block configuration × corners (24 symmetries) | 4.6 MB | 13 | 10.05 |
| block configuration × B-centres | 148.9 MB | 15 | 11.97 |
| corners × B-centres (24 symmetries) | 4.7 MB | 11 | 8.35 |

Each phase-1 solution gets at most 2¹⁸ phase-2 nodes before the search moves on. We tried 2¹⁶,
2¹⁸, 2²⁰ and no limit. The three limits give about the same lengths after 1 s; 2¹⁸ needed the
fewest nodes to the first solution and gave the fastest first solution on 48 threads. Without a
limit, the first phase-1 solutions hold the search for a long time. We also tried ordering
phase-1 solutions or phase-2 children by their bound, and a smaller 52 MB table set; none did
better.

## Measurements

These runs were made on 2026-11-12 on a 48-core AMD EPYC 9V45 virtual machine. Each run had the
machine to itself: it held an exclusive lock, and we checked that no other timed work
overlapped. Both solvers come from the same build of the branch (9dcb91c; the block split is
`--method chen`, the two-phase solver `--method two-phase`). Each run uses the same 10 seeded
random positions (`bench -N 10 --seed 1 -t 10`). The tiers are omelettehedron's table sizes:
tiny (42 MB of tables), small (134 MB) and medium (about 1 GB).

### Random positions, 1 and 48 threads

Each cell is the mean best length over the positions solved by then; a number in brackets is
how many of the 10 are solved, when not all. "First" is the mean first solution and when it
came. A dash means no position was solved yet.

| Tier | Threads | Solver | First | 0.01 s | 0.1 s | 1 s | 10 s |
|---|---:|---|---|---|---|---|---|
| tiny | 1 | two-phase | 22.9 after 3.30 s (9) | - | - | 23.0 (3) | 22.7 (9) |
| tiny | 1 | block split | 27.5 after 0.75 s | - | - | 24.9 (7) | 24.2 |
| tiny | 48 | two-phase | 23.3 after 0.095 s | - | 22.4 (5) | 22.1 | 21.5 |
| tiny | 48 | block split | 26.1 after 0.014 s | 25.6 (8) | 24.6 | 23.9 | 22.7 |
| small | 1 | two-phase | 23.0 after 1.92 s | - | - | 22.6 (5) | 22.4 |
| small | 1 | block split | 27.5 after 0.76 s | - | - | 24.8 (8) | 22.9 |
| small | 48 | two-phase | 23.0 after 0.037 s | 23.5 (2) | 22.4 | 21.6 | 21.5 |
| small | 48 | block split | 26.1 after 0.0064 s | 25.2 (9) | 23.7 | 22.6 | 21.7 |
| medium | 1 | two-phase | 23.3 after 0.35 s | - | 23.0 (2) | 22.7 (9) | 21.7 |
| medium | 1 | block split | 27.5 after 0.78 s | - | 25.0 (1) | 24.3 (7) | 22.6 |
| medium | 48 | two-phase | 23.5 after 0.0081 s | 22.9 (7) | 22.1 | 21.5 | 20.8 |
| medium | 48 | block split | 26.1 after 0.0070 s | 24.8 (8) | 22.9 | 22.0 | 21.5 |

On the tiny and small tiers the block split finds its first solution 2.5 to 6.8 times sooner,
and at 0.01 s on 48 threads it has solved more positions. On one thread at 1 s it has also
solved more (tiny 7 and small 8 of 10, against 3 and 5). On
the medium tier the two-phase solver's first solution comes about as soon (48 threads) or sooner
(one thread). Wherever both have solved everything, the two-phase solver's solutions are shorter.

### 16 threads, with memory

These runs were made the same day on the same machine, from the branch's final build (d72fe13).
Memory is the peak resident size of the whole process.

| Tier | Solver | Memory | First | 0.1 s | 1 s | 10 s |
|---|---|---:|---|---|---|---|
| small | two-phase | 0.30 GB | 23.0 after 0.090 s | 22.3 (6) | 22.3 | 21.5 |
| small | block split | 0.56 GB | 27.1 after 0.032 s | 24.8 | 22.7 | 22.0 |
| medium | two-phase | 1.21 GB | 23.5 after 0.017 s | 22.7 | 21.5 | 21.4 |
| medium | block split | 1.47 GB | 27.0 after 0.028 s | 24.0 | 22.6 | 21.6 |

### Hard positions

`solve -t 20` on the 24-symmetric position and on the hard positions in
[`bench/hard.tsv`](../bench/hard.tsv) (optimal lengths in brackets). Each cell is the first
solution and the best within 20 s, as length @ seconds.

| Position | small 1T two-phase | small 1T block split | small 48T two-phase | small 48T block split | tiny 1T two-phase | tiny 1T block split |
|---|---|---|---|---|---|---|
| 24-symmetric (22) | none | none | 23 @ 5.3, 22 @ 18.3 | 25 @ 15.4 | none | none |
| tri_adj+edge_opp (23) | none | 29 @ 1.1, 26 @ 3.3 | 26 @ 3.1, 24 @ 7.4 | 29 @ 0.04, 26 @ 0.13 | none | 29 @ 1.2, 27 @ 1.3 |
| tri_adj+edge_adj (23) | none | none | 23 @ 2.9 | 29 @ 15.5, 26 @ 17.3 | none | none |
| tri_adj-[Uf_Fu;Db_Bd] (23) | 23 @ 0.68 | 29 @ 3.4, 23 @ 17.0 | 23 @ 0.017 | 29 @ 0.15, 23 @ 0.45 | 23 @ 2.4 | 29 @ 3.6, 26 @ 16.7 |
| tri_adj+[UB_FD,NL_MR] (23) | 25 @ 4.9 | 29 @ 8.8, 27 @ 19.9 | 25 @ 0.14, 23 @ 1.0 | 29 @ 0.18, 23 @ 8.9 | none | 29 @ 9.3, 28 @ 19.0 |
| edge_hard (21) | 23 @ 3.6, 22 @ 5.2 | 28 @ 0.23, 23 @ 1.9 | 23 @ 0.06, 21 @ 4.3 | 27 @ 0.007, 21 @ 7.9 | 23 @ 8.7, 22 @ 12.8 | 27 @ 0.33, 24 @ 7.1 |
| corner 6-flip (20) | none | none | 23 @ 8.7 | 24 @ 6.1, 23 @ 15.1 | none | none |

On one thread the block split solves more of these within 20 s (4 of 7 against 3 on the small
tier, 4 against 2 on the tiny tier), always with a long solution at first. On 48 threads the
two-phase solver solves all 7, and its best is as short or shorter on every one.

## Why it loses

A two-phase solver is fast when phase 1 ends somewhere phase 2 is cheap. omelettehedron's
two-phase solver puts all the search into phase 1. Phase 1 brings the puzzle into the
positions solvable with U, F, BR and BL turns alone, and phase 2 is read off the exact table of
those 933,120 positions, with no search. Every phase-1 solution is at once a complete solution, and
the time goes into covering as many phase-1 paths as possible, pruned by large exact tables.

The block split moves the search into phase 2. Its phase 1 is short and cheap: at most 14
moves, with an exact table. But phase 2 is a search in 1.59 × 10¹⁴ positions, and the only
tables that fit cover two of its three parts at a time. At a typical phase-2 start they say 10
or 11 moves, while the phase 2 actually found is about 16. IDA* with a bound 5 moves short does
a great deal of extra work:

* on the first random position, a single phase-1 solution's phase-2 iteration at bound 16 takes
  1 to 2 × 10⁵ nodes;
* most phase-1 solutions use up their 2¹⁸ phase-2 nodes without finding anything;
* the first solution of a random position takes about 10⁶ to 10⁷ nodes (1 thread, small tier:
  0.9, 1.7, 3.7 and 8.2 × 10⁶ on four positions).

An exact phase-2 table would remove the phase-2 search. But at 1.59 × 10¹⁴ entries it would
take about 40 TB at 2 bits per position, or about 1.7 TB with 24 symmetries. omelettehedron's
largest table tier is about 30 GB.

The proposal estimated that phase 2 always finishes in about 18 to 20 moves; the phase 2 of the
solutions found here is about 16, in line with that. The split itself is sound; the cost lies in
pruning a phase 2 this large.

## What might change this

* A much better phase-2 estimate, such as a table over the block configuration, the corners and
  part of the B-centres together. That table would be far larger than the 149 MB one tried here.
* Settings where a quick first solution matters more than its length. The block split's
  first-solution lead is real on the small and tiny tiers, but it is gone by 0.1 s on 48 threads.

The God's-number side of the proposal (an upper bound of 14 plus the phase-2 maximum) needs the
largest phase-2 distance from the 3.4 × 10¹¹ phase-1-solved positions. A search for it covers
the whole phase-2 space; we did not run one.

## Reproducing

```sh
git fetch origin neg/chen-twophase && git checkout neg/chen-twophase
cargo build --release
B=target/release/omelettehedron
$B gen --tier small && $B gen --method chen --tier small
$B bench --method chen --tier small -N 10 --seed 1 -t 10 -j 48
$B bench --method two-phase --tier small -N 10 --seed 1 -t 10 -j 48
```

Thanks to Chen Shuang for the idea and for the state counts that made it easy to check.
