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
space and free memory are checked for all tables of the tier (`--force` skips the check);
`--min-free-memory GB` leaves that much memory to the rest of the machine (added to the check,
and the large tier's 21 GB table stops generating, writing nothing, if the free memory falls
below it). `--moves slice` adds the slice turns (each counts one move), with tables of their own
in that metric (`*-slice.prun`; small and medium tiers: the large tier is face turns only);
`--moves wide` only changes the output, as for the two-phase solver.

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
| `large` | corners x edges, A-centres/edges (with chirality) on both sides, and corners x B-centres x edge triple to solved for the 4 A-triples and the 4 B-triples (no corner x centre lookups: a triple lookup is never below the corner x centre value of its orbit) | 24.4 GB |

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
* **Corners x B-centres x edge triple** (large tier, `optimal/cbe.rs`): the table follows the
  three edges of A-face 0 in order, so only the 6 class-preserving symmetries that fix A-face 0
  reduce it: 172 corner classes. A lookup for A-face `k` conjugates with a symmetry that maps
  A-face `k` onto A-face 0 (the 6 such symmetries give the same classes), and a lookup for a
  B-face composes `SWAP_AB` in front, as for the A-centres above: eight lookups, one table.
  The index maps the corner state to its representative, the B- (or A-) colouring by two
  shuffles and a multiset rank, and the three edges' positions by the symmetry's position map
  (the pieces that the symmetry sends onto triple 0, in triple 0's order).

### Checks

* `optimal/sym.rs`: the coordinate maps of all 48 symmetries equal the sticker-model
  conjugation (random positions); `SWAP_AB` equals `Cube::conjugate`.
* `optimal/coords.rs`: the A lookup equals the B lookup of the `SWAP_AB` conjugate; both
  indices are invariant (up to equivalent indices) under their symmetries; ranks round-trip.
* `optimal/tables.rs`: the reduced corner x centre table equals an unreduced 89 MB BFS on
  3,000 random positions and all their neighbours, through both centre orbits.
* `optimal/cbe.rs`: every lookup's symmetry maps the corner state onto its representative
  and the lookup's triple onto triple 0; the eight lookups of a position and of each of its
  48 conjugates agree (up to order and equivalent indices); a decoded entry has its own index;
  the B-face lookups equal the A-face lookups of the `SWAP_AB` conjugate. **Brute force:** a
  sparse BFS over the table's own indices (the generator's steps and equivalents) to depth 5
  gives every position within 5 moves (all of them for lookup 0, an eighth for the others) the
  distance of its raw projection, computed independently from that ball of positions.
  (`--ignored`, with the 21 GB table: 3,000 sampled entries walk down to solved, equivalent
  entries agree, every entry was reached.)
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
* Large tier (`--ignored`, needs its tables): lengths equal the brute force on 36 scrambles of
  0-12 moves, every bound is at most the brute-force distance (300 positions up to 10 moves
  out), the first five lookups equal the medium tier's, and the tables are consistent.

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
instead of 3. The stages add dependent memory rounds (up to five per node instead of one), but
the fewer misses win on wall time. Timing on the quiet machine (after the two-phase proofs,
under the bench lock), medium tier, against main 3deb77e, three interleaved A/B rounds (the
order alternating), `scripts/bench-optimal.sh` sets:

| set | threads | main 3deb77e (s, median) | optimal/order (s, median) | ratio (rounds) | nodes/s |
|---|---|---|---|---|---|
| 10 x 17 moves, seed 1717 | 1 | 10.33 | 8.50 | 0.818 (0.818, 0.828, 0.818) | 1.42e7 -> 1.74e7 |
| 10 x 19 moves, seed 1919 | 48 | 41.19 | 33.42 | 0.811 (0.805, 0.818, 0.811) | 6.07e8 -> 7.54e8 |

Same lengths (means 16.50 / 18.10) and, on one thread, the same nodes (147,723,722); on 48
threads the node counts differ by under 0.01% between runs (the parallel split), and the
lookups drop to 0.308x. The statistics also show where the
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
0.1% of the best one here, by the pruning statistics of a 13-move set). Medium tier with slices
(corners x edges in the slice metric: 9,580,032,000 entries, max 11, mean 9.27; generated in 47
minutes on the busy machine): 139,550,571 nodes on the same set (0.108x the small tier, 120x
the face metric's 1,162,092), same lengths, 1.37 lookups per node. On the 17-move set (medium,
one thread): mean optimal length 15.5 (16.5 in face turns), 47,025,092,761 nodes (318x the face
metric's 147,723,722; two positions of 17 slice-metric moves take 1.7e10 and 2.8e10), 1.38
lookups per node.

What the slices cost face turns: the search is compiled once per move count, and with 16
moves it applies face moves only (`OptMoves::apply_n`), reading the face rows of the corner
table (the first 32 KB; the slice rows, a slice keeping the corner id, are a second block), so
the face metric has no branch on the move kind and the same table footprint as before.
Face-metric time against main (4c96995), medium tier, quiet machine, interleaved rounds: 48
threads, 19-move set 0.999x (0.999, 0.997, 0.999); 1 thread, 17-move set 1.015x (1.016,
1.016, 0.999) and in six more rounds 1.031x (1.009 to 1.058); same nodes. Earlier layouts: a
branch in the corner move 1.083x / 1.059x (1T / 48T), one 64-byte row of 32 moves per corner
1.012x / 1.024x, the two blocks without `apply_n` 1.025x / 1.024x. With main's round 3 merged
in (main 7b93c60, three rounds): 1 thread 1.007x (0.966, 1.003, 1.024), 48 threads 1.001x
(1.010, 0.995, 0.996): at parity. (The earlier 1T residue came with identical `dfs` machine
code, so it was code placement, which every merge shifts.) Stacked on `optimal/r3-aepre`:
1.015x (1.017, 1.001, 1.018) / 1.004x (1.000, 1.010, 1.004) against it.

## Optimisation round 2: a bound for the nodes with 9 to 11 moves left

**Estimate before building.** With medium + AE, the expanded nodes of the 13-move set (10
scrambles, seed 1717: 18,582 nodes, 1,425 expansions, one thread) sit at 10 and 11 moves left
(498 and 682 expansions; 38 at 9). Each candidate table was evaluated exactly on that recorded
search, without building it: a breadth-first ball of radius 6 around the solved projection,
and for every expansion a depth-first search (bounded by the tier's own tables) that decides
whether the projection's distance fits in the moves left. The search order does not change
and the bounds are admissible, so a node is kept iff every ancestor passes, and the new node
count is exact on the set (`optimal/estimate`, an experiment branch: a recording hook in the
search and `examples/estimate.rs`).

| candidate (to solved) | entries | size | nodes | cut at 9 / 10 / 11 moves left |
|---|---|---|---|---|
| today (medium + AE) | | | 18,582 | of 38 / 498 / 682 expansions |
| A-centres x edge A-triples x chirality, both sides | 2.75e11 | 68.7 GB | 0.340x | 22 / 442 / 477 |
| the same without chirality | 3.43e10 | 8.6 GB | 0.679x | 22 / 348 / 89 |
| both centre orbits (A x B colourings), 8 symmetries | 1.37e10 | 4.3 GB | 0.880x | |
| **corners x B-centres x ordered positions of an A-face's edges, 4 triples x 2 sides** | 8.39e10 | 21.0 GB | **0.342x** | 28 / 483 / 430 |
| the same with the positions as a set + chirality (the two-phase CBE layout) | 2.80e10 | 7.0 GB | 0.616x | 28 / 474 / 47 |
| the same with two of the three edges, ordered | 8.39e9 | 2.1 GB | 0.706x | 27 / 389 / 4 |
| A-centres of two colours x edge A-triples x chirality, both sides | 1.37e10 | 3.4 GB | 0.808x | 15 / 245 / 13 |

(The corner x B-centre x triple tables keep the 6 symmetries that fix A-face 0: 172 corner
classes x 369,600 colourings x 1,320 ordered positions (440 as a set + chirality); one table
serves all eight lookups. Without the corners, an A-centre/edge table keeps only the 4
frame-fixing symmetries.) Dropping the chirality loses most of the A-centre/edge table's
gain at 11 moves left (477 -> 89 cuts), and the table that keeps it does not fit next to the
proofs; the ordered triple table is as strong at a third of the size. The set + chirality
variant (7.0 GB, 0.62x) is the fallback for smaller 64-bit machines (with the medium tier,
10.4 GB in all; not built). WebAssembly stays on the small tier: wasm32 has 4 GB of address
space, and every table beyond the small tier has more than 2^32 entries. A check on two
positions of the 15-move set (83,960 nodes) gave 0.351x for the ordered triple table, as on the
13-move set.

**The `large` tier (`optimal/cbe-solved`).** The ordered triple table (`opt-cbe-sym6.prun`,
`optimal/cbe.rs`): 83,913,984,000 entries, 20,978,496,152 bytes, generated in 98.7 s on the
free machine (48 threads; peak memory the table plus the medium tier), every entry reached, max
depth 14, mean 10.98; depths 0..14: 1, 4, 42, 461, 6,187, 77,829, 988,289, 12,226,937,
149,316,933, 1,699,817,301, 14,901,628,662, 50,036,333,989, 17,080,092,594, 33,494,770, 1.
The tier reads corners x edges, AE and its twin, and the eight triple lookups (K = 11): the
corner x centre lookups are dominated (each triple lookup is at least the corner x centre value
of its orbit; with every lookup read they prune no child). Lookup order AE, B-triple 3, AE twin,
the other triples, corners x edges (greedy on the pruning masks of the 17-move set).

Nodes (one thread; the 48-thread counts vary slightly with the split), same lengths everywhere:

| set | medium + AE | large | ratio | completed iterations | all optimal solutions (`-n`) |
|---|---|---|---|---|---|
| 10 x 15 moves, seed 1717 | 1,162,092 | 638,270 | 0.549x | 0.285x | 4,461,034 -> 1,290,732 (0.289x) |
| 10 x 17 moves, seed 1717 | 147,723,722 | 76,064,980 | 0.515x | 0.350x | 626,281,636 -> 210,553,326 (0.336x) |
| 10 x 19 moves, seed 1919 (48 threads) | 25,165,160,444 | 6,272,933,698 | 0.249x | | |

The full trees shrink as estimated; to the first solution the gain is smaller on the shorter
sets because the last iteration's task order (from the previous iteration's statistics)
changes with the bounds and puts the solution elsewhere (position 3 of the 15-move set: 9,344
-> 216,266 nodes in its last iteration). Lookups per node: 2.36 on the 17-move set (2.53 on
the 15-move set) against 1.51, 1.08 of them in the 21 GB table.

Time on the quiet machine (bench lock, three interleaved rounds, the same binary for both
tiers, medians):

| set | threads | medium (s) | large (s) | ratio (rounds) |
|---|---|---|---|---|
| 10 x 15 moves, seed 1717 | 1 | 0.069 | 0.066 | 0.956 (0.955, 0.970, 0.975) |
| 10 x 15 moves, seed 1717 | 48 | 0.021 | 0.018 | 0.861 |
| 10 x 17 moves, seed 1717 | 1 | 8.077 | 7.147 | 0.885 (0.890, 0.885, 0.888) |
| 10 x 17 moves, seed 1717 | 48 | 0.236 | 0.209 | 0.888 |
| 10 x 19 moves, seed 1919 | 48 | 31.399 | 14.753 | **0.470** (0.459, 0.471, 0.468) |

A node costs about 1.7 times as much (more lookups, and the triple table misses the cache and
the TLB), so the gain grows with the depth: 0.89x at 17 moves, 0.47x at 19.

**The corner 6-flip** (`bench/hard.tsv`, 48 threads, quiet machine): medium 20,555,511,618 nodes
in 29.0 s; large **96,270,740 nodes in 0.497 s** (214x fewer nodes, 58x faster), the same
20-move solution. Its edges and centres are solved, so AE is 0 and every root bound is 9, but
deeper the triple table sees the corners together with the edges and centres they disturb
(depth 19: 4.17e9 -> 5.82e7 nodes).

## The optimal solver against the two-phase solver's optimal mode

Both solvers prove optimality; the same positions, quiet machine, bench lock, two interleaved
rounds (`opt-vs-2p.sh`): the two-phase solver (`bench`/`solve --tier large`, main 4c96995)
against this solver's large tier. Times to a proven optimum (both rounds):

| positions | threads | two-phase (s) | optimal solver (s) | ratio | nodes, two-phase / optimal |
|---|---|---|---|---|---|
| 10 x 17 moves, seed 1717 (mean 16.5) | 1 | 2.82, 2.87 | 7.06, 7.28 | 2.51x, 2.54x | 1.06e8 / 7.61e7 |
| 10 x 19 moves, seed 1919 (mean 18.1) | 48 | 7.95, 8.01 | 14.23, 14.93 | 1.79x, 1.86x | 1.24e10 / 6.27e9 |
| random position, seed 1 (21) | 48 | 53.8, 53.5 | 91.7, 91.4 | 1.70x, 1.71x | 8.60e10 / 3.85e10 |
| random position, seed 2 (21) | 48 | 61.2, 60.4 | 210.9, 214.6 | 3.45x, 3.55x | 9.69e10 / 9.10e10 |
| random position, seed 3 (20) | 48 | 4.94, 4.97 | 106.0, 105.5 | 21.4x, 21.2x | 7.88e9 / 4.52e10 |
| edge_hard (21, `bench/hard.tsv`) | 48 | 71.3 | 344.7 | 4.84x | 1.08e11 / 1.67e11 |
| corner 6-flip (20, normal FTO) | 48 | 626.5 | **0.478** | **0.00076x (1,310x faster)** | 7.97e11 / 9.60e7 |

The two-phase solver stays faster everywhere except on the corner 6-flip, where the large
tier's triple table sees the flipped corners together with the edges and centres they would
disturb (every root bound is 9, but its depth-19 iteration has 98.6% fewer nodes than the
medium tier's). The optimal solver usually
searches fewer nodes (0.4x to 0.9x on the sets and the 21-move positions), but each costs 3.5
times as much (1.08e7 against 3.75e7 nodes/s on one thread, 4.4e8 against 1.57e9 on 48; the
large tier reads 2.4 to 2.6 table entries per node, nearly half of them in the 21 GB table). It
is weakest where the optimum is reached early: on the
random position of length 20 the two-phase solver has its 20 within 0.5 s and only has to
exclude 19 (7.9e9 nodes), while IDA* must finish the whole depth-19 iteration (4.5e10).

## Optimisation round 3: the cost per node

**Profile** (perf with task-clock sampling, the VM has no hardware counters; 1 thread, 17-move
set, main 4c96995). Over 90% of the time is in `dfs` (everything is inlined into it). The
largest single cost is the stall on the table read: 36.5% of the medium tier's `dfs` samples
and 44% of the large tier's sit on the instruction after the table load, although every
entry is prefetched before it is read (a stage computed all its children's indices and
prefetched them right before reading them). Then 6.3% on the A-centre/edge index's u32
symmetry-class table (1.5 MB, misses the 1 MB L2), 3.6% (medium) / 10% (large) on the
prefetch instructions themselves (line fill buffers full); the rest is index arithmetic,
thinly spread, and move application (the shuffles) is 2-4%.

Changes (`optimal/r3-prefetch`), each with the same tree, nodes and lookups:

* **Next-stage prefetch:** a child that passes stage k gets its stage-k+1 index computed and its
  entry prefetched at once, inside stage k's read loop, so the miss overlaps the rest of stage
  k as well.
* **Sibling pipeline:** before the search descends into a surviving child, the next surviving
  child's children are generated and their first-stage entries prefetched (two frames per
  level), so those misses overlap the current child's subtree. A node's children are counted
  when the search enters it.

Time against main 4c96995 (quiet machine, bench lock, three rounds with the order rotating;
1 thread on the 17-move set, 48 threads on the 19-move set):

| variant | medium 1T | medium 48T | large 1T | large 48T |
|---|---|---|---|---|
| next-stage prefetch | 0.976, 0.982 | 0.996, 0.986 | 0.772, 0.815 | 0.821, 0.835 |
| sibling pipeline | 0.953 | 0.940 | 0.890 | 0.973 |
| **both** | **0.928** | **0.945** | **0.818** | **0.875** |
| AE class-entry prefetch pass (`neg/optimal-r3-aeclass`) | 1.104 | 1.120 | 1.007 | 1.009 |

(Two numbers: two sessions.) The large tier gains most from the next-stage prefetch (1.36
lookups per node after the first stage, against 0.51 for medium), the medium tier from the
pipeline (its first stage is most of its lookups). The AE prefetch pass computes the
A-centre rank of every child of a stage first and prefetches its class entry: the extra pass
costs more than the stall it hides.

**Piece positions from SIMD compares (`optimal/r3-masks`).** The profile after the prefetch
changes put 15% of the medium tier's samples on the A-centre/edge chirality bits
(`FcMaps::ebits`: an inverse permutation scattered into memory, then the parity of four
position triples). It and the triple table's index (which also built an inverse permutation
to find three pieces) now take each piece's position as a one-bit mask from a SIMD compare
(`fast::label_mask`); comparing two masks compares the positions, so a parity is three mask
comparisons. Against round 3: medium 0.885x (1T) / 0.885x (48T), large 0.907x / 0.925x; same
nodes.

**The first stage's A-centre/edge class entries prefetched while the children are generated
(`optimal/r3-aepre`).** The A-centre/edge index first reads a u32 symmetry-class entry (1.5
MB table, misses L2; 9.6% of the medium tier's samples after the prefetch changes) and only
then knows the table entry. When the sibling pipeline generates a node's children, it now
also computes each child's A-centre rank for the first stage (`Bounds::pre`) and prefetches
its class entry; the first stage's indices finish from that rank (`Bounds::index_with`), with
a whole sibling subtree of slack. (Doing this as an extra pass before every stage was slower,
`neg/optimal-r3-aeclass`.) On top of the SIMD positions: medium 0.898x (1T) / 0.919x (48T),
large 0.946x / 0.941x; same nodes.

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
2. **Corners x B-centres x one A-face's edges, to solved** ("CBE-to-solved"): round 2 (the
   `large` tier, above). Only the 6 symmetries that fix A-face 0 reduce it (the others move the
   triple), so it has 172 x 369,600 x 1,320 = 8.39e10 entries (21.0 GB), looked up four times
   (each A-face mapped to face 0 by a symmetry) and four more times through `SWAP_AB`.
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
