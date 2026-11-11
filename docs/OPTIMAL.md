# The optimal solver

[Documentation home](../README.md) · [Two-phase design](PHASES.md) · [Results](RESULTS.md)

`omelettehedron optimal` is a second solver, separate from the two-phase solver: a single-phase
IDA* straight to the solved state, in the style of the optimal 3x3x3 solvers (pattern databases
of distances to solved). Everything on this page is about the **optimal solver**; the two-phase
solver (`omelettehedron solve|gen|bench`, [PHASES.md](PHASES.md)) is unchanged.

Status: **v0** (correct first, fast later). The two-phase solver's optimal mode is much faster
today; this solver is the base for the optimisation loop.

## Usage

```sh
# optimal solution (the first one found; IDA* finds the shortest first)
omelettehedron optimal solve --tier medium "U R' BL F D B' L BR"
# all optimal solutions, up to 20
omelettehedron optimal solve --tier medium -n 20 "..."
# give up after 60 s: prints the lower bound proven so far
omelettehedron optimal solve --tier medium -t 60 "..."
# give up if the optimum is longer than 18 moves
omelettehedron optimal solve --tier medium -l 18 "..."
# a position given as the algorithm that solves it
omelettehedron optimal solve --tier medium --inverse "U F BL B' ..."
# per-iteration progress (depth, nodes, tasks, time) on stderr
omelettehedron optimal solve --tier medium -v "..."
# generate the tables ahead of time
omelettehedron optimal gen --tier medium
# benchmark: seeded random-move scrambles (the two-phase bench sets), one line per position
omelettehedron optimal bench --tier medium -N 10 --scramble-len 17 --seed 1717 -j 1
```

`--tier` is required. `--tables DIR` / `$FTO_TABLES` / `./tables` hold the tables, as for the
two-phase solver; the optimal solver's files start with `opt-`. Before generating, the disk
space and free memory are checked for all tables of the tier (`--force` skips the check).
`--moves slice` adds the slice turns (each counts one move), with tables of their own in that
metric (`*-slice.prun`); `--moves wide` only changes the output, as for the two-phase solver.

## Search

IDA* on the whole position (fixed-corner frame, 16 internal moves, see
[PHASES.md](PHASES.md#7-implementation-notes-differences-from-the-original-proposal)):

* Iteration `D = h(root), h(root) + 1, ...` is a depth-first search over the canonical move
  sequences (`cube::allowed_after`: never the same face twice in a row, single before wide on
  an axis) that skips a child whose lower bound exceeds the moves left. It finds exactly the
  solutions of `D` moves. All shorter iterations found none, so the first solution is optimal,
  and each completed iteration proves "no solution of at most `D` moves" (`-t` reports that
  bound).
* A node generates its children (a move is one byte shuffle per vector and a corner-table
  lookup), computes all their table indices and prefetches them, then reads the values, so the
  cache misses of the siblings overlap. Values are exact distances mod 3; each node carries its
  exact distances and decodes its children's (`prune::decode`). The children of a node with no
  move left after them are leaves: only the solved state passes, and it is recognised without
  lookups.
* Threads: the top levels (at most 4 moves) are split into tasks that the workers take in order
  from a shared counter. Tasks are ordered by a statistic of the previous iteration (sum of 8^-h
  over the surviving nodes of the task's subtree), so that the last iteration tends to meet a
  solution early. Small iterations run on the calling thread.
* **Symmetric start**: when some of the 48 symmetries fix the position, each task level keeps
  only the first task (in move-prefix order) of each class of conjugate positions (key: the
  least conjugate under the position's symmetry group). Conjugate positions are equally far from
  solved; the first optimal solution is not lost because the smallest prefix of a canonical
  optimal solution is always the first of its class (a smaller conjugate task would give a
  canonical optimal solution with a smaller prefix: at the junction a wide turn followed by its
  axis's single turn is reordered single first, which only makes the prefix smaller, and a merge
  would make the solution shorter). Off with `-n` (all optimal solutions are wanted) and with
  `FTO_OPT_SYMSTART=0`. Random positions have no symmetry; the user's hard positions
  (`bench/hard.tsv`) have 48, 24, 8, 16, 1 and 48, and the filter divides their node counts by
  about the order of their symmetry group (below).
* Stopping rules: by default the first solution ends the search; `-n N` collects up to `N`
  solutions of the optimal length (the iteration that finds the first one runs to its end, so
  with fewer than `N` the output is all of them); `-t` stops and reports the proven lower bound;
  `-l L` gives up above `L` moves. Every solution is checked on the cubie model before it is
  reported (and on the sticker model by the tests). `-v` prints the root's bounds, how many of
  the 48 symmetries fix the position, and each iteration (depth, nodes, tasks, time).

The normal FTO has no inverse trick: same-coloured centres are interchangeable, so a position
and its inverse can have different lengths, and no table is looked up on the inverse.

## Tables

Every table is the exact distance of a **projection** of the position to the projection of the
solved state, generated by breadth-first search (`prune::bfs`, 2 bits per entry, distance
mod 3) from that single state. Every solution projects to a path of the same length, so each
table is a lower bound, and the bound of a position is the **maximum** of its lookups. Each
table is consistent (its value changes by at most one per move).

| table | projection | entries | bytes | max | mean (entries) |
|---|---|---|---|---|---|
| `opt-cc-sym24.prun` | corners x B-centre colouring, 24 symmetries | 51 x 369,600 = 18,849,600 | 4.7 MB | 11 | 8.35 |
| `opt-edges.prun` | edge permutation (even) | 12!/2 = 239,500,800 | 60 MB | 11 | 8.47 |
| `opt-ce-sym48.prun` | corners x edge permutation, 48 symmetries | 40 x 239,500,800 = 9,580,032,000 | 2.4 GB | 13 | 10.56 |
| `ae-sym288-triples.prun` (two-phase) | A-centres x edge A-triples, to `H_A` | 515,961,600 | 129 MB | 12 | 10.10 |
| `ae-sym288-chiral.prun` (two-phase) | A-centres x edge A-triples x chirality, to `H_A` | 4,127,692,800 | 1.03 GB | 14 | 10.94 |

The last two are the two-phase solver's A-centre/edge tables (its small and medium tiers'
files, shared, read through its public index `tables::AeSym`). They hold the distance to
`H_A` = <U, F, BR, BL> of the A-centres and edges, not to solved; solved is in `H_A`, so they
are lower bounds on the distance to solved too, and so is the same table on the `SWAP_AB`
conjugate (the position's distance to `H_B`, which holds solved as well). Each is still an
exact BFS distance (mod 3) that changes by at most one per move, so the search decodes it from
the parent's value like the others. The optimal solver computes the two-phase coordinates of a
position and of its conjugate from its own byte vectors (`coords::FcMaps`, checked against
`fast().of`).

Tiers (`--tier`):

| tier | lookups per position | memory |
|---|---|---|
| `small` | edges, corners x B-centres, corners x A-centres, A-centres/edge triples on both sides | 195 MB |
| `medium` | corners x edges, corners x B-centres, corners x A-centres, A-centres/edges (with chirality) on both sides | 3.4 GB |

(The corner x edge table subsumes the edge table: every corner x edge path is an edge path.)

Histograms (entries at depth 0, 1, ...; the reduced tables count indices, not positions):

* edge permutation: 1, 16, 184, 1,992, 20,814, 207,312, 1,939,311, 16,114,044, 92,853,054,
  124,720,121, 3,643,521, 430;
* corners x B-centres (reduced): 1, 2, 17, 132, 1,677, 19,138, 207,478, 1,765,748,
  8,311,698, 8,199,439, 344,207, 63;
* corners x edges (reduced): 1, 2, 17, 112, 1,614, 19,862, 253,507, 3,007,751, 36,398,797,
  408,573,186, 3,386,024,615, 5,654,317,203, 91,435,214, 119.

Generation (`optimal gen`, release build, on the busy machine: the two-phase proofs held all 48
cores, so these ran on a few cores at `nice 10`): small tier 18 s, corners x edges 181 s
(peak memory: the tables themselves).

### Symmetry

The distance to solved is unchanged by conjugating a position with any of the **48**
symmetries of the octahedron (stickers moved and recoloured; a reflection mirrors the turns
onto turns) followed by the whole-puzzle rotation `R` that brings the front corner piece home
again (`optimal/sym.rs`). `R` depends on the corner state, so a table reduced this way must
include the corners. On the coordinates, symmetry `s` sends corner state `c` to
`conj[s][c]`, and relabels centres/edges by `s` and repositions them by `t = R s`.

These are not the two-phase solver's 288 AE symmetries (independent relabelling and
repositioning, which preserve `H_A` but not the solved state).

* **Corners x centres**: the 24 class-preserving symmetries leave 51 corner classes. The index
  maps the corner state to its class representative and applies the same symmetry to the
  B-colouring (two byte shuffles and a multiset rank). Representatives with a stabiliser have
  several indices per position; the BFS fills all of them (`equivalents`).
* **The A-centres use the same table**: `SWAP_AB` (the 90° rotation about the fixed corner's
  axis) fixes the solved state and the front corner and swaps the face classes, so the
  corners x A-centres distance of a position is the corners x B-centres distance of its
  `SWAP_AB` conjugate. The index composes the swap into the symmetry (class-swapping relabel
  and reposition maps), so no second position is carried down the search.
* **Corners x edges**: an edge has one face of each class, so all 48 symmetries apply: 40
  corner classes (960 / 48 = 20 without stabilisers).
* An AE-to-solved table (A-centres x edge triples x chirality) can use only the symmetries
  that fix the front corner (4 class-preserving ones): it includes no corners, so no `R`.
  That is 2.7e11 entries (68 GB), not built.

### Checks

* `optimal/sym.rs`: the coordinate maps of all 48 symmetries equal the sticker-model
  conjugation (random positions); `SWAP_AB` equals `Cube::conjugate`.
* `optimal/coords.rs`: the A lookup equals the B lookup of the `SWAP_AB` conjugate; both
  indices are invariant (up to equivalent indices) under their symmetries; ranks round-trip.
* `optimal/tables.rs`: the reduced corner x centre table equals an unreduced 89 MB BFS on
  3,000 random positions and all their neighbours, through both centre orbits.
* On load, every table is walked down from random positions to the solved entry (a file of
  another layout fails), besides the size and checksum checks of `prune::Table2::load`.

## Tests (`tests/optimal.rs`)

* Optimal lengths equal the meet-in-the-middle brute force on 65 scrambles of 0-12 moves, on 1
  and 8 threads; every solution solves on the sticker model. (`--ignored`: 300 more.)
* Optimal lengths equal the **two-phase solver's optimal mode** (both exact) on 6 scrambles of
  11-15 moves.
* `-n`: the solutions are exactly the canonical optimal solutions found by brute-force
  enumeration (including a position with two, from a length-10 relation), on 1 and 8 threads;
  `-n 1` and `-n 2` cap the count. `-l` gives up below the optimum and finds it at or above;
  `-t` stops promptly with the root bound as lower bound.
* Node counts on one thread are deterministic, and the iterations' counts add up.
* Symmetric start: the 112 positions within 5 moves (distance 3 or more) that a symmetry fixes
  get the brute force's lengths, filtered (default) and unfiltered, on 1 and 8 threads.
* Tables: decoded child values equal walked exact values, values change by at most one per
  move, only the solved state is 0 in every lookup, the bound never exceeds the true distance;
  the edge table's first levels equal an independent BFS.

## Measured node counts

Nodes = children generated (counted as in the two-phase solver), on one thread unless noted
(deterministic: the same tree every run). The sets are the two-phase solver's optimal-mode sets
(`bench --scramble-len L --seed S`, seeded random-move scrambles). Optimal solver v0, commit
0a9f8cd; times are not meaningful yet (busy machine, `nice 10`).

| set | tier | mean optimal length | total nodes | mean nodes per position |
|---|---|---|---|---|
| 10 x 15 moves, seed 1717 | small | 14.5 | 873,050,832 | 8.73e7 |
| 10 x 15 moves, seed 1717 | medium | 14.5 | 28,682,242 | 2.87e6 (0.033x small) |
| 10 x 17 moves, seed 1717 | medium | 16.5 | 3,852,034,640 | 3.85e8 |
| 10 x 19 moves, seed 1919 (48 threads) | medium | 18.1 | 376,022,447,738 | 3.76e10 |

The optimal lengths equal the two-phase solver's (17-move set mean 16.5, 19-move set 18.1,
docs/RESULTS.md). The 19-move set ran on 48 threads at `nice 10` beside the two-phase proofs
(2,642 s in all, 1.4e8 nodes/s); with several threads the node count of the last iteration
depends on when the first solution is met, so it varies a little between runs. Per position
(19-move set): 19 moves 5.18e+10, 19 moves 9.60e+10, 17 moves 4.45e+08, 17 moves 5.45e+08, 19 moves 4.03e+10, 19 moves 1.98e+10, 16 moves 3.60e+07, 17 moves 8.95e+08, 19 moves 8.72e+10, 19 moves 7.90e+10.

Per position (17-move set, medium tier, one thread):

| # | optimal | nodes | last three iterations (depth: nodes) |
|---|---|---|---|
| 0 | 17 | 672,765,370 | 15: 8,102,382, 16: 102,549,812, 17: 561,422,804 |
| 1 | 17 | 133,651,458 | 15: 9,098,336, 16: 112,530,458, 17: 11,218,896 |
| 2 | 15 | 7,242,568 | 13: 44,278, 14: 580,858, 15: 6,612,910 |
| 3 | 17 | 143,809,664 | 15: 6,514,098, 16: 87,421,710, 17: 49,361,098 |
| 4 | 17 | 401,202,730 | 15: 6,472,742, 16: 87,436,488, 17: 306,795,874 |
| 5 | 16 | 13,158,570 | 14: 623,006, 15: 8,005,202, 16: 4,480,432 |
| 6 | 15 | 5,219,336 | 13: 43,954, 14: 571,680, 15: 4,600,120 |
| 7 | 17 | 1,225,451,128 | 15: 8,367,906, 16: 105,847,592, 17: 1,110,512,724 |
| 8 | 17 | 269,803,928 | 15: 7,119,924, 16: 94,725,400, 17: 167,399,566 |
| 9 | 17 | 979,729,888 | 15: 7,190,244, 16: 95,168,084, 17: 876,789,836 |

Each iteration costs about 12.4 times the previous one (the canonical branching factor is
12.93): the bounds are still weak compared with the depth (corners x edges averages 10.56,
the optima of these sets are 15-19). The last iteration stops at the first solution; the task
ordering of the previous iteration's statistics puts it early in some positions (#1: 1.1e7 of
an estimated 1.4e9) and not in others (#0).

### Symmetric start on the hard positions

`bench/hard.tsv` positions, medium tier, one thread, iterations up to depth 15 (`-l 15`; none has
a solution that short), with and without the filter (`FTO_OPT_SYMSTART=0`):

| position | symmetries | nodes, filtered | nodes, unfiltered | ratio |
|---|---|---|---|---|
| tri_adj+edge_opp | 48 | 172,138 | 8,239,136 | 0.021 (47.9x fewer) |
| tri_adj+edge_adj | 24 | 287,114 | 6,875,800 | 0.042 (23.9x) |
| tri_adj-[Uf_Fu;Db_Bd] | 8 | 19,992,502 | 159,640,000 | 0.125 (8.0x) |
| tri_adj+[UB_FD,NL_MR] | 16 | 705,154 | 11,215,152 | 0.063 (15.9x) |
| pure_corn_6flip | 48 | 685,776 | 32,842,576 | 0.021 (47.9x) |
| edge_hard | 1 | (no symmetry: unchanged) | | |

The hard positions on 48 threads (at `nice 10` beside the two-phase proofs, about 1.35e8
nodes/s), medium tier with symmetric start (commit df9e167), 300 s each, the 6-flip without a
limit:

| position | known | root bound | symmetries | result | nodes |
|---|---|---|---|---|---|
| tri_adj+edge_opp | 23 | 12 | 48 | no solution below 20 moves (300 s) | 41,139,142,554 |
| tri_adj+edge_adj | 23 | 11 | 24 | no solution below 20 (300 s) | 40,986,460,294 |
| tri_adj-[Uf_Fu;Db_Bd] | 23 | 10 | 8 | no solution below 19 (300 s) | 41,302,196,336 |
| tri_adj+[UB_FD,NL_MR] | 23 | 11 | 16 | no solution below 20 (300 s) | 42,527,118,126 |
| edge_hard | 21 | 12 | 1 | no solution below 19 (300 s) | 39,977,083,142 |
| pure_corn_6flip | 20 (normal FTO) | 9 | 48 | **20, optimal, in 232 s** | 30,526,471,644 |

The 6-flip's optimum agrees with the two-phase solver's (20, docs/RESULTS.md: 661 s with
`--lower-bound 19` after an earlier run had excluded 18). Its solution
`D B BR L' BL U' L' D B' L' D B L' U D BL' BR' D' L B'` exchanges centres of a colour (not a
pure-FTO solution), like the two-phase solver's. Each further iteration of the others costs
about 12.5 times the last (about 9 minutes for depth 20 of tri_adj+edge_opp, 2 hours for 21).

## Optimisation round 1

Node counts on one thread (deterministic), lookups = table entries read. The v0 rows are the
measurements above (v0 counted no lookups: every child of a node with moves left after it read
all K = 3 tables).

**A-centre/edge bound (`optimal/ae`).** The two-phase AE table on the position and on its
`SWAP_AB` conjugate, two more lookups per child (K = 5). Same optimal lengths on every set:

| set | tier | v0 nodes | with AE: nodes | ratio | lookups per node |
|---|---|---|---|---|---|
| 10 x 15 moves, seed 1717 | small | 873,050,832 | 13,456,478 | 0.0154 (65x fewer) | 5.00 |
| 10 x 15 moves, seed 1717 | medium | 28,682,242 | 1,162,092 | 0.0405 (24.7x) | 5.00 |
| 10 x 17 moves, seed 1717 | medium | 3,852,034,640 | 147,723,722 | 0.0383 (26.1x) | 5.00 |

AE averages 10.9 (with chirality) against 10.56 for corners x edges, and it measures the
centres, which no other table of the tiers combines with the edges. On the corner 6-flip it
adds nothing (its edges and centres are solved: AE is 0 at the root on both sides; iterations
to 16 within 0.4% of v0's nodes), so that proof stays as in v0.

**Joint corner/B-centre/edge-triple table (`neg/optimal-cbe`, not merged).** The two-phase
CBE table (7.0 GB, distance to `H_A` as a saturating code: at most 6, 7, 8, at least 9) for all
four A-faces' triples, on the position and on its `SWAP_AB` conjugate, as a second stage (a
child that passes the exact tables is checked against it when it has fewer than 9 moves left),
in a `large` tier of 10.4 GB. Same lengths; on top of medium + AE: 1,162,092 -> 1,161,858 nodes
on the 15-move set (0.9998x), 147,723,722 -> 144,410,408 on the 17-move set (0.978x), with
0.0002 CBE lookups per node. With bounds around 11 hardly any node has fewer than 9 moves left,
so the table, which halved the two-phase solver's phase-1 tree near its leaves, prunes little
here: 2% fewer nodes for 7 GB more.

**Lookup order and early exits (`optimal/order`).** Which lookups prune a child does not depend
on the order they are read in, so the tree (and the node count) is the same in any order; only
the lookups change. A child is dropped at its first value that exceeds its moves left, and each
stage computes and prefetches the indices of all children still alive before reading any. The
order comes from one run that reads everything (`FTO_OPT_PRUNE_STATS=1`): for every child the
set of lookups that alone would prune it (bench JSON `prune_masks`), which gives the lookups of
every order. Medium tier, 15-move set: 7.9% of the children survive; alone, AE prunes 74%, its
`SWAP_AB` twin 72%, corners x edges 61%, each corner x centre lookup 1.4%. Orders: medium AE,
AE twin, corners x edges, corners x B, corners x A; small AE, AE twin, corners x B, corners x A,
edges (the corner x centre table is cache-resident, the edge table is not). Same lengths and
the same node count per position on every set:

| set | tier | lookups before | after | per node | per node by table (in order) |
|---|---|---|---|---|---|
| 10 x 15 moves, seed 1717 | medium | 5,809,810 | 1,769,241 (0.305x) | 1.52 | AE 1.000, AE twin 0.260, CE 0.105, CB 0.079, CA 0.079 |
| 10 x 15 moves, seed 1717 | small | 67,281,740 | 19,660,878 (0.292x) | 1.46 | AE 1.000, AE twin 0.220, CB 0.082, CA 0.080, edges 0.078 |
| 10 x 17 moves, seed 1717 | medium | 738,617,940 | 223,268,625 (0.302x) | 1.51 | AE 1.000, AE twin 0.252, CE 0.105, CB 0.078, CA 0.077 |

The tables that miss the cache (AE 1.03 GB, corners x edges 2.4 GB) are read 1.36 times per node
instead of 3. The stages add dependent memory rounds (up to five per node instead of one), so
the time gain needs measuring once the machine is free. The statistics also show where the
nodes are: 92% of the children have 9 or 10 moves left and fewer than 0.01% at most 6 (17-move
set: 1,700 of 1.5e8), so exact checks near the leaves (for example one move left: one of the 16
neighbours of solved) cannot change the node count measurably; not done.

**Slice turns (`--moves slice`, `optimal/slice`).** 24 internal moves (the slices are the middle
layers of the four axes, turned like the back face; they move no corner) with the slice metric's
canonical order (`cube::allowed_after`: on an axis only X Xw and X' Xw'). Every table of the
tiers is a distance in the metric of the move set, so the slice metric has its own files
(`opt-*-slice.prun`, and the two-phase A-centre/edge tables `ae-sym288-*-slice.prun`); the
symmetry reductions hold, as every symmetry maps slices to slices. The search is compiled once
per number of moves (16 or 24), so the face metric runs the same code. The symmetric start
still keeps the first optimal solution: the canonical form of a two-move element of an axis
starts with its single turn, the smallest move index of the axis, so the reordering at the end
of a task's prefix only makes the prefix smaller. Tests: the brute force in the slice metric on
1 and 8 threads (scrambles of face, wide and slice turns; the solutions use slices and are never
longer than in the face metric), the two-phase solver's optimal mode with slices, consistency of
the slice tables, the symmetric start on symmetric positions of the slice metric.

Slice-metric tables (entries; max, mean): corners x centres 18,849,600 (10, 7.52), edges
239,500,800 (9, 7.17), A-centre/edge triples to `H_A` 515,961,600 (11, 8.83). Small tier, one
thread, 10 x 15-move face-turn scrambles (seed 1717): mean optimal length 13.6 with slices
(14.5 without), the same lengths as the two-phase solver's optimal mode with slices on all ten,
1,297,054,007 nodes (96x the face metric's 13,456,478: about 18 children per node instead of
13, and shallower tables), 1.34 lookups per node (the face-metric lookup order is also within
0.1% of the best one here, by the pruning statistics of a 13-move set).

## Missing (v0)

* Timings: the machine was busy with long two-phase proofs, so only node counts are measured.
* Printing improvements (`-p`) has no use here (the first solution is optimal); resuming a
  proof from a known lower bound.
* Lanes / micro-threading and huge-page tuning: not done.

## Next optimisations (ranked)

1. **A bound on the A-centres and edges (AE) and its `SWAP_AB` twin.** Done in round 1
   (`optimal/ae`, above): the two-phase solver's AE-to-`H_A` table, 24.7x / 26.1x fewer nodes
   on the 15 / 17-move sets. A real AE-to-solved table would be stronger, but in the fixed
   frame it keeps only the 4 class-preserving symmetries that fix the front corner (it has no
   corners, so no re-fixing rotation): 2.7e11 entries (68 GB) with chirality, 3.4e10 (8.5 GB)
   without.
2. **Corners x B-centres x one A-face's edges, to solved** ("CBE-to-solved"): 960 x 369,600 x
   1,320 ordered positions of 3 edges, reduced by the 24 class-preserving symmetries (about 2e10
   entries, 4.9 GB), looked up four times (each A-face mapped to face 0 by a symmetry) and four
   more times through `SWAP_AB` (A-centres with a B-face's edges). The two-phase solver's
   joint table of this shape halved its nodes.
3. **Timing, once the machine is free** (node counts unchanged): lanes (several depth-first
   stacks per worker, so the misses of different nodes overlap, not only those of siblings) or
   the two-phase solver's sibling pipeline; staged lookups with an early exit (corners x edges
   first, the corner x centre lookups only for children it keeps); no lookups with one move left
   (a child must then be one of the 16 neighbours of solved); 1 GiB pages for the corners x
   edges table (`FTO_HUGE1G`); a SIMD even-permutation rank.
4. **A centres table** (A x B colourings to solved; 8 front-corner-fixing symmetries, 1.7e10 to
   3.4e10 entries): complementary to corners x edges.

Not applicable here: sibling pruning with bounds learned from subtrees (with consistent tables
the bound learned from a subtree never exceeds the node's own: for a descendant at distance j,
`max_k h_k(D) - j <= max_k h_k(node)`), and anything on the inverse position (the normal FTO's
inverse can have another length).
