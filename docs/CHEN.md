# Chen Shuang's block split

[Documentation home](../README.md) · [Two-phase design](PHASES.md) · [Results](RESULTS.md)

On 2026-11-12 Chen Shuang proposed a new two-phase split for the FTO: phase 1 solves the
tetraminx centres, and phase 2 solves the rest with moves that keep each face's six
tetraminx-centre pieces in at most two same-coloured blocks of three, much like Square-1's
constraint. We implemented it as a second solver (`--method chen`) and measured it against
omelettehedron's two-phase solver. It is not in the released solver; this page explains why.

The code is on the branch [`neg/chen-twophase`](https://github.com/ArhanChaudhary/omelettehedron/tree/neg/chen-twophase)
(`src/chen/`, with a longer technical write-up in its `docs/CHEN.md`). Chen Shuang then
suggested letting phase 1 stop earlier; that variant is measured in
[A relaxed phase 1](#a-relaxed-phase-1).

## In short

* **The numbers check out.** Our model reproduces the 460 × 3⁴ = 37,260 block configurations,
  the phase-2 space of 1.59 × 10¹⁴ positions (the "158T"), and the phase-2 state counts for
  depths 0 to 5.
* **Faster first solution, but a longer one.** On the tiny and small table tiers it finds a
  first solution 3 to 10 times sooner. That solution is 3 to 4 moves longer: 26 to 27 moves,
  against 23.
* **Longer solutions after the first moment.** At 10 s the two-phase solver's mean solution is
  shorter at every table size and thread count, by 0.32 to 1.49 moves, measured on 100 random
  positions. Every gap is well outside its 95% interval (the smallest is 0.32 ± 0.11), and the
  block split is shorter on at most 4 of the 100 positions. The block split solves more
  positions early on: on one thread at 1 s, 83 of 100 on the tiny tier against 30.
* **More memory.** It needs 158 MB of extra tables, or 0.16 to 0.31 GB more at peak.
* **A relaxed phase 1 helps, but not enough.** Chen Shuang suggested letting phase 1 stop at any
  block configuration. That makes the block split's first solution about 1 move shorter, and
  its best at 10 s 0.6 to 0.7 moves shorter on one thread (about equal on 16 and 48 threads).
  On the small tier it draws level with the two-phase solver (at 10 s on one thread, and up to
  1 s on 16 and 48 threads). Elsewhere the two-phase solver is still 0.21 to 0.45 moves shorter at 10 s.
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

These runs were made on 2026-11-13 on a 48-core AMD EPYC 9V45 virtual machine, whose cores come
in 6 groups of 8 that each share a 32 MB cache. Both solvers come from the same build of the
branch (d72fe13; the block split is `--method chen`, the two-phase solver `--method
two-phase`). Every run uses the same 100 seeded random positions, with 10 s each (`bench -N 100
--seed 1 -t 10`). Each batch of runs held an exclusive lock, so no other timed work ran
alongside it:

* the six one-thread runs (three tiers, two solvers) ran side by side, each pinned to one core
  of its own cache group;
* the 16-thread runs ran in pairs, the two solvers at once on cores 0 to 15 and 16 to 31;
* each 48-thread run ran alone.

The tiers are omelettehedron's table sizes: tiny (42 MB of tables), small (134 MB) and medium
(about 1 GB). The first version of this page used 10 positions. The conclusions have not
changed, and the differences below now come with 95% intervals.

### Mean solution length

Each cell is the mean best length over the positions solved by then; a number in brackets is
how many of the 100 are solved, when not all. "First" is the mean first solution and when it
came. A dash means no position was solved yet. (On one thread on the tiny tier, the two-phase
solver solved 90 of the 100 within 10 s.)

| Tier | Threads | Solver | First | 0.01 s | 0.1 s | 1 s | 10 s |
|---|---:|---|---|---|---|---|---|
| tiny | 1 | two-phase | 22.92 after 3.19 s (90) | - | 21.67 (3) | 22.70 (30) | 22.61 (90) |
| tiny | 1 | block split | 27.02 after 0.55 s | 25.00 (10) | 25.88 (24) | 25.25 (83) | 24.12 |
| tiny | 48 | two-phase | 23.09 after 0.097 s | 21.80 (5) | 22.48 (63) | 22.08 | 21.50 |
| tiny | 48 | block split | 26.23 after 0.0093 s | 25.50 (76) | 24.36 | 23.58 | 22.76 |
| small | 1 | two-phase | 23.01 after 1.56 s | 20.00 (1) | 21.67 (6) | 22.60 (50) | 22.39 |
| small | 1 | block split | 27.03 after 0.51 s | 25.00 (10) | 25.44 (25) | 24.67 (88) | 22.92 |
| small | 16 | two-phase | 23.02 after 0.098 s | 21.80 (5) | 22.46 (63) | 22.17 | 21.56 |
| small | 16 | block split | 26.75 after 0.029 s | 25.27 (48) | 24.21 (94) | 22.74 | 21.95 |
| small | 48 | two-phase | 23.03 after 0.038 s | 22.70 (20) | 22.57 (98) | 21.80 | 21.40 |
| small | 48 | block split | 26.23 after 0.0064 s | 24.86 (86) | 23.39 | 22.50 | 21.72 |
| medium | 1 | two-phase | 23.06 after 0.27 s | 21.60 (5) | 22.58 (38) | 22.61 (97) | 21.70 |
| medium | 1 | block split | 27.05 after 0.51 s | 25.00 (10) | 25.12 (24) | 24.16 (85) | 22.53 |
| medium | 16 | two-phase | 23.09 after 0.016 s | 22.71 (49) | 22.33 (99) | 21.62 | 21.19 |
| medium | 16 | block split | 26.70 after 0.028 s | 25.45 (56) | 23.67 (93) | 22.43 | 21.66 |
| medium | 48 | two-phase | 23.11 after 0.0079 s | 22.59 (69) | 21.93 | 21.45 | 20.96 |
| medium | 48 | block split | 26.08 after 0.0056 s | 24.48 (88) | 22.79 | 21.91 | 21.43 |

### Paired differences

The block split's best length minus the two-phase solver's on the same position, averaged over
the positions both have solved by then, with a 95% interval (1.96 standard errors). Positive
means the two-phase solver is shorter. A number in brackets is how many positions are counted,
when fewer than 100; a dash means fewer than 10. The last column counts the positions at 10 s
where the two-phase solver is shorter, where they are equal, and where the block split is
shorter.

| Tier | Threads | 0.1 s | 1 s | 10 s | At 10 s: two-phase shorter / equal / block split shorter |
|---|---:|---|---|---|---|
| tiny | 1 | - | +2.15 ± 0.77 (26) | +1.49 ± 0.25 (90) | 72 / 15 / 3 |
| tiny | 48 | +1.75 ± 0.35 (63) | +1.50 ± 0.23 | +1.26 ± 0.22 | 72 / 28 / 0 |
| small | 1 | - | +1.66 ± 0.60 (44) | +0.53 ± 0.19 | 39 / 57 / 4 |
| small | 16 | +1.55 ± 0.43 (60) | +0.57 ± 0.19 | +0.39 ± 0.14 | 29 / 71 / 0 |
| small | 48 | +0.81 ± 0.23 (98) | +0.70 ± 0.18 | +0.32 ± 0.11 | 28 / 72 / 0 |
| medium | 1 | - | +1.55 ± 0.40 (82) | +0.83 ± 0.18 | 60 / 39 / 1 |
| medium | 16 | +1.34 ± 0.29 (92) | +0.81 ± 0.18 | +0.47 ± 0.13 | 39 / 61 / 0 |
| medium | 48 | +0.86 ± 0.19 | +0.46 ± 0.13 | +0.47 ± 0.13 | 39 / 61 / 0 |

Every difference is positive, and no interval reaches zero. Where fewer than 100 positions are
counted, they are the ones the two-phase solver had solved by then, which may be the easier
ones.

On the tiny and small tiers the block split finds its first solution 3 to 10 times sooner, and
early on it has solved more positions:

* at 0.01 s on 48 threads: 76 and 86 of 100, against 5 and 20;
* at 1 s on one thread: 83 and 88, against 30 and 50.

On the medium tier the two-phase solver's first solution comes sooner on 1 and 16 threads, and
the block split's comes 1.4 times sooner on 48. Wherever both have solved a position, the
two-phase solver's solution is shorter on average.

### Memory

Peak resident size of the whole process, in GB, as two-phase / block split.

| Tier | 1 thread | 16 threads | 48 threads |
|---|---|---|---|
| tiny | 0.21 / 0.37 | not run | 0.21 / 0.45 |
| small | 0.30 / 0.53 | 0.31 / 0.58 | 0.31 / 0.62 |
| medium | 1.20 / 1.44 | 1.21 / 1.48 | 1.22 / 1.52 |

### Hard positions

These runs are from 2026-11-12 (build 9dcb91c), each alone on the machine: `solve -t 20` on the
24-symmetric position and on the hard positions in
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

## A relaxed phase 1

After reading the first version of this page, Chen Shuang suggested a change. Phase 1 can stop
at any of the 460 × 3⁴ block configurations, not only at those with every hexagon intact,
because phase 2 can start from any of them. We implemented it as `--method chen-relaxed` on the
branch [`work/chen-relaxed`](https://github.com/ArhanChaudhary/omelettehedron/tree/work/chen-relaxed).
Phase 2 is unchanged: its tables already cover every block configuration.

* **Shorter phase 1.** Phase 1 now ends in any of 447,120 configurations (37,260 for each of the
  12 ways the hexagons can sit) instead of 972. On the medium tier its new table is exact: at
  most 11 moves, 8.45 on average, against 14 and 10.94. On the small tier the table is a lower
  bound (at most 10, 7.61 on average).
* **More memory.** The new phase-1 table takes 1.03 GB on the medium tier and 129 MB on the
  small tier. Peak memory is 2.47 to 2.56 GB on medium (strict block split 1.44 to 1.52,
  two-phase 1.20 to 1.22) and 0.66 to 0.74 GB on small. The tiny tier has no table to build it
  from.

These runs were made on 2026-11-13 the same way as above, on the same 100 positions (build
badc27e). The strict block split, rerun from this build beside the relaxed runs on 1 and 16
threads, gave the same lengths as the earlier build at 10 s on all 100 positions. So the
comparisons below use the earlier runs of the other two solvers.

### Mean solution length, relaxed

Same format as the [main table](#mean-solution-length).

| Tier | Threads | First | 0.01 s | 0.1 s | 1 s | 10 s |
|---|---:|---|---|---|---|---|
| small | 1 | 25.93 after 0.4 s | 24.00 (8) | 24.64 (28) | 23.40 (92) | 22.22 |
| small | 16 | 25.49 after 0.022 s | 24.22 (51) | 22.91 (95) | 22.18 | 21.79 |
| small | 48 | 25.14 after 0.0046 s | 23.74 (89) | 22.35 | 21.88 | 21.61 |
| medium | 1 | 25.78 after 0.38 s | 24.00 (8) | 24.48 (31) | 23.05 (93) | 21.95 |
| medium | 16 | 25.49 after 0.022 s | 24.12 (56) | 22.55 (94) | 21.88 | 21.54 |
| medium | 48 | 25.16 after 0.0039 s | 23.46 (91) | 22.20 | 21.79 | 21.41 |

### Against the strict block split

The relaxed split's best length minus the strict split's on the same position, as in
[Paired differences](#paired-differences): negative means the relaxed split is shorter.

| Tier | Threads | 0.1 s | 1 s | 10 s | At 10 s: strict shorter / equal / relaxed shorter |
|---|---:|---|---|---|---|
| small | 1 | - | -1.19 ± 0.42 (81) | -0.70 ± 0.22 | 10 / 36 / 54 |
| small | 16 | -1.27 ± 0.28 (89) | -0.56 ± 0.24 | -0.16 ± 0.19 | 19 / 50 / 31 |
| small | 48 | -1.04 ± 0.22 | -0.62 ± 0.22 | -0.11 ± 0.18 | 17 / 57 / 26 |
| medium | 1 | - | -1.10 ± 0.47 (79) | -0.58 ± 0.23 | 11 / 38 / 51 |
| medium | 16 | -1.01 ± 0.26 (87) | -0.55 ± 0.22 | -0.12 ± 0.16 | 15 / 62 / 23 |
| medium | 48 | -0.59 ± 0.24 | -0.12 ± 0.19 | -0.02 ± 0.14 | 15 / 66 / 19 |

The relaxed target is never worse. Its first solution is 0.9 to 1.3 moves shorter and comes in
0.70 to 0.78 of the time. At 0.1 s and 1 s it is clearly shorter, except at 1 s on the medium
tier on 48 threads. At 10 s it is 0.58 to 0.70 moves shorter on one thread; on 16 and 48
threads the two are about equal by then.

### Against the two-phase solver

The relaxed split's best length minus the two-phase solver's: positive means the two-phase
solver is shorter.

| Tier | Threads | 0.1 s | 1 s | 10 s | At 10 s: two-phase shorter / equal / relaxed shorter |
|---|---:|---|---|---|---|
| small | 1 | - | +0.44 ± 0.46 (48) | -0.17 ± 0.23 | 22 / 44 / 34 |
| small | 16 | +0.27 ± 0.34 (62) | +0.01 ± 0.23 | +0.23 ± 0.18 | 30 / 53 / 17 |
| small | 48 | -0.23 ± 0.25 (98) | +0.08 ± 0.21 | +0.21 ± 0.15 | 26 / 63 / 11 |
| medium | 1 | +1.00 ± 1.31 (10) | +0.46 ± 0.37 (91) | +0.25 ± 0.20 | 34 / 48 / 18 |
| medium | 16 | +0.32 ± 0.24 (93) | +0.26 ± 0.18 | +0.35 ± 0.14 | 37 / 55 / 8 |
| medium | 48 | +0.27 ± 0.24 | +0.34 ± 0.17 | +0.45 ± 0.14 | 41 / 55 / 4 |

The relaxed target closes part of the gap. On the small tier the two are about equal up to 1 s
on 16 and 48 threads, and at 10 s on one thread (−0.17 ± 0.23). There, by 1 s the relaxed split
has solved 92 of the 100 positions against the two-phase solver's 50. This is the one setting
measured where the block split keeps up. Everywhere else the two-phase solver is shorter at
10 s, by 0.21 to 0.45 moves, and no interval reaches zero.

On 48 threads, in a short screening run, the relaxed split searched fewer phase-1 nodes per
second than the strict one (3.6 × 10⁷ against 6.0 × 10⁷). A likely cause is the many more
phase-1 solutions, each of which starts a phase-2 search; we have not looked into it.

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
* A different partition: one with few enough phase-1 positions for an exact table, and a phase 2
  that still has a strong estimate. The relaxed target moves in that direction.
* Settings where a quick first solution matters more than its length. The block split's
  first-solution lead is real on the small and tiny tiers. On 48 threads it is gone by about
  0.1 s on the small tier and 1 s on the tiny tier; on one thread it lasts longer.

The God's-number side of the proposal (an upper bound of 14 plus the phase-2 maximum) needs the
largest phase-2 distance from the 3.4 × 10¹¹ phase-1-solved positions. A search for it covers
the whole phase-2 space; we did not run one.

## Reproducing

```sh
git fetch origin neg/chen-twophase && git checkout neg/chen-twophase
cargo build --release
B=target/release/omelettehedron
$B gen --tier small && $B gen --method chen --tier small
$B bench --method chen --tier small -N 100 --seed 1 -t 10 -j 48
$B bench --method two-phase --tier small -N 100 --seed 1 -t 10 -j 48

# the relaxed phase 1
git fetch origin work/chen-relaxed && git checkout work/chen-relaxed
cargo build --release
$B gen --method chen-relaxed --tier small
$B bench --method chen-relaxed --tier small -N 100 --seed 1 -t 10 -j 48
```

Thanks to Chen Shuang for the idea, for the relaxed target, and for the state counts that made it easy to check.
