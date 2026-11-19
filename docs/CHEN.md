# The block split (`--method chen`)

[Documentation home](../README.md) · [Two-phase design](PHASES.md) · [Results](RESULTS.md)

A second two-phase solver with a different split, proposed by Chen Shuang, in `src/chen/`. It
runs with `solve`, `bench` and `gen` and `--method chen`; everything else is as for the two-phase
solver (`--tier`, `-t`, `-l`, `-n`, `-j`, `--random`, patterns).

```sh
omelettehedron gen --method chen --tier small                 # make the phase-2 tables
omelettehedron solve --method chen --tier small -t 5 -p improvements "..."
omelettehedron bench --method chen --tier small -N 10 -t 10 -j 1
```

## The split

Each class-A face (U, F, BR, BL) has a hexagon: its 3 centres and the A-stickers of its 3 edges.

* **Phase 1** solves the A-centres and the edges: every hexagon holds its own 3 centres and 3
  edges, turned by any multiple of 120° (up to a whole-puzzle turn). That is "the A-centre/edge
  part solved", so the tier's A-centre/edge table is the phase-1 table (exact for medium and
  larger, a lower bound for small and tiny). At most 14 moves, 10.9 on average.
* **Phase 2** solves the rest with fewer moves. U, F, BR and BL turns are always allowed. A turn
  of R, L, B or D cuts three hexagons into 3 stickers inside its layer and 3 outside; it is
  allowed only if, on each, the 3 inside are one colour and the 3 outside are one colour. A
  phase-2 path ends once it reaches `H_A` (with an R, L, B or D turn, or at once): the exact
  `H_A` table gives the rest. Every allowed path to solved ends that way, so this finds the
  shortest allowed completion (within the node budget below).

As in the two-phase solver, phase 1 targets both face classes, phase-1 depths grow by one (on
several threads a thread with no task left starts the next depth), and each phase-1 solution
gets a phase-2 search for anything shorter than the best so far. A
phase-1 solution gets at most 2^18 phase-2 nodes before the search moves on (`FTO_CHEN_BUDGET`,
0: no limit). Phase 2 is restricted, so the search cannot prove a solution optimal: without
`-t`, `-l` or `-n` it runs until no phase-1 depth below the best is left.

## The phase-2 space

Without a fixed corner (the pattern view: the 8 face turns), the allowed moves reach 37,260
hexagon configurations from solved, in 460 colour shapes of 81 block turns each, and never move
a block to another face. In the solver's frame (front corner held) a wide turn also turns the
whole puzzle, so each configuration is seen in 12 frames: 447,120 configurations, no two of
them a turn of each other. The phase-2 space is therefore 447,120 x 960 corner states x 369,600
B-centre colourings = 1.59 x 10^14 positions (37,260 x 11,520 x 369,600 in the pattern view), 12
times the 1.3 x 10^13 that a closed set of configurations would give.

## Breadth-first search of phase 2 (`chen-p2bfs`)

`src/bin/chen-p2bfs.rs` and `src/chen/p2bfs/`: an exact breadth-first search of the whole phase-2
space, for its distance distribution, its antipodes and the deepest position of `H'`.

**Size.** Each pattern-view configuration goes with all 11,520 corner states, not 960: the fixed
frame's 12 frames are the 12 places of corner piece 0 (6 corners, 2 orientations) in the pattern
view, and 8 of the 16 face turns move that corner, so no frame is closed. Taking the corners up
to a whole turn (960 states) would identify solved with 11 other positions (solved with the
corners turned as a whole) and give distances to the nearest of the 12. So the space is 37,260 x
11,520 x 369,600 = 1.586 x 10^14 positions, 6,610,224,733,560 up to the 24 symmetries; the
search's index has 1,611 classes x 369,600 x 11,520 = 6.86 x 10^12 slots, 857 GB at a bit each.

**Search.** A *block* is the 11,520 corner states of one (configuration class, B-centre
colouring): a move takes a block to one block, its bits permuted by a table (the move and the
symmetry back to the class representative). The move graph is undirected, so layer d + 1 is the
children of layer d minus layers d and d - 1: only layers are stored, one file per class,
blocks compressed (runs of empty or full blocks, lists of bits, plain blocks) in checksummed
chunks of 4,096. A step splits the classes into parts that fit the memory budget (classes joined
by U, F, BR and BL turns kept together: those moves never leave a part), and per part one thread
reads the layer-d files of the classes with a move into it while the workers OR each block's
children into the part's blocks in memory; the stabilisers' images are added, the new blocks
minus layers d and d - 1 are written and counted (slots, positions, positions up to symmetry by
fixed points). Each part is recorded in `manifest.json` when written, so a stopped search
resumes at the next part. `antipodes` lists a layer's positions up to symmetry with a shortest
solution (IDA* with the AC, PX and CC tables, in face turns) and the scramble that is its inverse.

A dense block is permuted 32 bits at a time: the corner states are numbered permutation x 32 +
orientations, and each (move, symmetry) maps the 32 orientations of a permutation to those of
one permutation by an affine map of the bit index (byte lookups, then a flip). Layer files can
be spread over several disks (`--store`, class c on the (c mod n)-th), and several threads read
them during a push (`--readers`).

**Cost.** Uniformly random phase-2 positions (480, IDA*) are at distance 14-23, 19.6 on average
(20: 26%, 19 and 21: 21% each; 11 of the 480 not finished after 12 minutes, at 22 or more), so
about six layers are dense (a bit in most blocks: 857 GB each as stored) and two of them are on
disk at once at the peak: the search needs about 1.7 TB of layer space and writes about 5.7 TB
in all (an estimate from that distribution).

```sh
chen-p2bfs run --dir DIR --ram-gb 55          # start or resume (Ctrl-C or --max-minutes: resume later)
chen-p2bfs status --dir DIR                   # layers so far, totals, the H' classes' deepest layers
chen-p2bfs antipodes --dir DIR --tables T --out antipodes.tsv
```

**Result** (complete run, 2026-11-15): every one of the 1.586 x 10^14 positions is reached;
the distance is at most 25. Positions per distance (up to symmetry in brackets): 13: 9.8 x
10^10, 14: 4.0 x 10^11, 15: 1.5 x 10^12, 16: 4.4 x 10^12, 17: 1.14 x 10^13, 18: 2.30 x 10^13, 19:
3.59 x 10^13, 20: 4.00 x 10^13, 21: 2.80 x 10^13, 22: 1.18 x 10^13, 23: 2.0 x 10^12, 24: 5.3 x
10^10, 25: 1,180,056 (49,804). The deepest position of `H'` is at 24.

**Joint check** (`chen-p2bfs joint`, `src/chen/p2bfs/joint.rs`): every position is within 35
moves. The relaxed phase 1 needs at most 11 moves (10 or fewer: 10 + 25 = 35). A position at 11
needs more than 35 with this split only if every optimal phase-1 solution ends in the last
layer S25; every such position is `s v` with `s` in S25 and `v` 11 moves along which the
phase-1 distance rises by one each step. The check computes exact phase-1 distances of all
4,127,692,800 entries of the relaxed index (equal to the table entry by entry, mod 3), the
entries on a shortest path to the 3 entries at distance 11 (25, 95, 158, 171, 185, 182, 176,
165, 157, 64, 8, 3 per distance; checked to be exactly the backward closure; as A-centre/edge
states 5,412, 25,152, 41,184, 45,792, 47,232, 47,808, 46,368, 45,504, 41,472, 16,560, 1,536 and
96, the weights giving every distance's exact state count), and from one
position of each of the 49,804 classes of S25 (the 24 symmetries keep the phase-1 distance,
checked on 20,000 random positions, and S25) enumerates the 1,534,416 candidates `s v`: each
has another optimal phase-1 solution ending outside S25, so none needs an explicit solution.
Sampled candidates (60) were solved explicitly, phase 1 then phase 2 by IDA*, in at most 33
moves, replayed on the independent sticker model; their optimal phase-1 solutions end 135 times
in S25 (at least once each) and 95,001 times outside. With all its checks: 2 min 19 s at 48
threads, 6.0 GB.

**Checks** (`tests/p2bfs.rs`, and on the command line): the positions per depth of the full
space are 1, 16, 160, 1,344, 10,464, 78,152 (Chen Shuang's twsearch counts), 560,100, 3,812,363
(the sticker model of `/data3/scratch/chen`, and a plain search without symmetry), 24,234,475 at
depth 8; the instances with the corners or the B-centres dropped give the AC and PX tables'
histograms exactly; a quotient instance (corner permutation x colourings up to renaming the
colours, 8.9 x 10^9 slots, depth 17) gives the in-memory table search's counts at every depth,
in one part and in many small parts stopped and resumed after every few; the orbit counts equal a
count of one canonical slot per orbit; sampled slots of layers 3 to 8 have that IDA* distance;
solutions of layer-8 positions replay to solved, with allowed moves only, on the independent
sticker model, and `omelettehedron show` gives the stated pattern for each scramble.

## The radius-1 phase-2 space (`chen-p2bfs --radius 1`)

The block split's phase-2 space is small because the move rule is strict: a shortest solution of
a depth-25 position with any moves takes 18-21 moves (100 sampled classes). The *radius-1 space*
relaxes the rule as little as possible for an exact search: every position whose relaxed phase-1
distance is at most 1, with every move that stays among them. The distance depends on the
A-centres and edges alone and is the same in every frame, so the space is a set of pattern-view
configurations with every corner state and B-centre colouring: 315,900 configurations (the
37,260 at distance 0 and 278,640 one face turn away), 13,245 symmetry classes (1,611 and
11,634), 1,345,036,492,800,000 positions, 56,043,189,900,216 up to symmetry, 8.5 times the block
split's space. Only 5.51 moves per configuration stay in the set on average. Phase 1 then has to
reach distance 1 only (10 moves at most), so the bound is 10 + the deepest layer; for positions
whose configuration is at distance 0 (the block split's phase-2 positions) the deepest layer
containing one is the gain over Chen's 25.

`src/chen/radius1.rs` builds the set in the fixed frame from the 447,120 configurations and
their neighbours under the 16 moves (3,790,800), each configuration's frame following its
parent's by the move's frame change (one rotation per frame and move, read off the block split's
space), and checks the construction: 12 frames of 315,900, one frame per configuration however
it is reached, every configuration turned by its frame's rotation in frame 0, frame 0 closed
under the 24 symmetries, every move undone by its inverse, the block split's move table agreed
with. The search is the same engine with these configurations in place of the 37,260
(`Kind::radius`), and a class order of greedily grown clusters (`src/chen/p2bfs/order.rs`: no
move is always allowed here, so the A-component order does not apply): a dense step reads 2.48
times the class files (2.0 in the block split's search).

**Checks.** Positions per depth 0-7 (1, 16, 208, 2,112, 18,144, 142,568, 1,063,524, 7,564,483;
at distance 1: 0, 0, 48, 768, 7,680, 64,320, 501,696, 3,732,344) agree between the search, a
plain search over whole cube states in the fixed frame (`positions --radius1`) and an independent
search on the Python sticker model; depth 5 already holds 78,248 distance-0 positions against
the block split's 78,152 (shortcuts through distance 1). The 100 sampled depth-25 classes' radius-1
distances (`radius-solve`, section above) are looked up in the layers as the search makes them.

**Run.** A dense layer is 7.05 TB (a bit per slot); two layers at the peak, 14.1 TB of a 17 TB
disk; 42 GB for a part (78 classes, about 170 parts per dense level, 44 GB in all); a dense level
moves about 39 TB, 5-6 hours at the disk's 2-2.6 GB/s. Started 2026-11-17 05:34 UTC.

**Result** (complete, 2026-11-19 01:12 UTC): every one of the 1.345 x 10^15 positions is
reached; the deepest layer is 26 (24 positions, one class, at distance 1). Block positions
(distance 0) per depth at the end: 23: 307,366,026,849, 24: 91,479,353, 25: 3, none deeper.

**Certificates for 33** (phase 1 to a block configuration in at most 10 moves on the better
side, then the radius-1 distance: 10 + 23 once every block position at depth 24 or 25 has a
solution of at most 23 moves, any moves, on both face classes):

* `layer-positions --distance 0 --expand`: depth 24's block classes only (1,611 of the 13,245
  classes read): 3,815,395 positions up to symmetry, 4,197,343 slots (each class's count equal
  to the level's), every pattern read back as its slot; their orbits give 91,479,353 positions,
  all distinct, each in the layer.
* Depth 25 (its layer deleted by the search): `neighbours --distance 0` of the 24 depth-26
  positions gives exactly 3 positions, one orbit (class 1578), in none of layers 0, 24 and 26.
  A depth-26 position's block neighbours are at depth 25, and no block position is deeper than
  25, so these are the 3. Optimal distance 22 (`solve optimal --tier huge`: no solution of 21
  moves), e.g. `B U F' BR' U' BR BL F' L BL' L R BR R' B' BL B U' BR F U BR`.
* Each class is solved by the two-phase solver (`-l 23`, large tier; about 3.5 CPU-seconds a
  class). `conjugate-solutions` and `verify-side-b --depth 24` then give the certificates of
  both face classes. Each conjugated solution is checked and replayed on the Python sticker
  model's tables. The positions of the B side are checked to have their `SWAP_AB` conjugate in
  layer 24, at relaxed distance 0. Both tools now stream their files a chunk at a time and
  check that the positions are distinct at the end, so 91 million lines fit in memory.

## Oriented phase 1 (`chen-oriented`)

`src/chen/oriented.rs`, `p1disk.rs`, `p1check.rs` and `src/bin/chen-oriented.rs`: exact phase-1
distances to targets finer than `H'`, every class-A hexagon one colour *and unturned*:

* **oriented**: the solved configuration (A-centres and edges) in each of the 12 frames;
* **T1**: those and the 8 configurations one class-A turn (U, U', F, F', BR, BR', BL, BL') from
  each, 9 per frame (108).

Neither is closed under turning a hexagon in place, so a state is the whole edge permutation:
369,600 A-centre colourings x 239,500,800 even edge permutations = 8.85 x 10^13 states, each
with 960 x 369,600 positions (corners up to a whole turn, B-centres). Both targets are closed
under the two-phase table's 288 pairs (`src/sym.rs`; a pair is a rotation times a conjugation
that fixes solved and maps class-A turns to class-A turns), so the index has one row per
A-colouring class (1,396) of 239,500,800 slots: 3.34 x 10^11 slots, 41.8 GB a bit plane (the
optimal solver's giant-tier table has the same rows and holds the oriented distance mod 3). A
column is (edge coset of the two-phase table, turn of each A-triple in `Z_3^4`): a move or a pair
maps a coset's 81 slots by an affine map of `Z_3^4`, so they move as one 81-bit block.

**Search.** A level at a time on disk: a slot is new at level d if it is unvisited and a
neighbour is at level d - 1 (pull). A destination row reads at most 16 source rows of the
previous level's plane, whole (29.9 MB) into a least-recently-used cache of 6-8 GB; the rows go
in a greedy order that keeps the sources cached. Bit planes `V.bits` and `L<d>.bits` with
direct I/O, then `D.nib` (4 bits a slot, 167 GB). A dense level took 6-8 minutes with both
searches running on a loaded machine; each search about an hour.

**Result** (2026-11-19; states in the fixed frame, 960 x 369,600 positions each; divide by 12
for states with one edge held):

| depth | oriented: states | slots | T1: states | slots |
|---:|---:|---:|---:|---:|
| 0 | 12 | 1 | 108 | 9 |
| 1 | 192 | 9 | 1,152 | 33 |
| 2 | 2,208 | 36 | 13,056 | 92 |
| 3 | 24,576 | 138 | 144,960 | 594 |
| 4 | 272,448 | 1,072 | 1,602,816 | 6,090 |
| 5 | 3,005,376 | 11,472 | 17,636,064 | 67,351 |
| 6 | 32,993,280 | 125,821 | 193,129,704 | 736,130 |
| 7 | 360,490,332 | 1,369,696 | 2,104,985,484 | 7,982,882 |
| 8 | 3,920,792,028 | 14,850,724 | 22,823,672,916 | 86,363,283 |
| 9 | 42,424,555,428 | 160,482,714 | 245,348,265,276 | 927,788,548 |
| 10 | 452,961,167,244 | 1,712,577,206 | 2,543,499,292,020 | 9,614,511,442 |
| 11 | 4,558,966,848,552 | 17,229,723,946 | 21,391,629,875,244 | 80,821,227,795 |
| 12 | 32,844,356,972,184 | 124,069,059,872 | 58,294,300,125,552 | 220,143,305,683 |
| 13 | 49,515,312,979,248 | 186,991,529,581 | 6,019,575,659,448 | 22,741,121,307 |
| 14 | 1,101,155,571,960 | 4,163,384,427 | 1,276,200 | 5,561 |
| 15 | 4,932 | 85 | 0 | 0 |
| all | 88,519,495,680,000 | 334,343,116,800 | 88,519,495,680,000 | 334,343,116,800 |

* **Oriented**: at most **15** moves; 4,932 states at 15 (**411** with one edge held, 85 slots,
  1.75 x 10^12 positions); mean 12.520877 over positions. Every count over 12 equals the
  held-edge phase-1 distribution published with the 36-move upper bound (411 at 15), and the
  slot counts equal the giant tier's table's histogram.
* **T1**: at most **14** moves; 1,276,200 states at 14 (106,350 with one edge held, 5,561
  slots, 4.53 x 10^14 positions); mean 11.759395. So phase 1 to the ball of radius 4 around T1
  takes at most 10 moves.
* The strict target (`H'`) for comparison: at most 14, mean 10.938721.

**Checks** (logs and lists in `/data3/oriented4/stageA/`):

* Local consistency of every one of the 3.34 x 10^11 slots of each `D.nib` (`check`): the slots
  at 0 are exactly the target slots, and every other slot's smallest neighbour value is its
  value less one (the neighbours found slot by slot with the digit map, not the search's block
  permutation). Both tables pass, so each holds the exact distance given the move model.
* The move model: the index's steps equal cube moves (tests over random slots and every 997th
  coset), and the 288 pairs map each state's neighbours to its image's neighbours up to a
  rotation and fix both target sets (tests).
* Table-free searches from the targets over whole A-centre/edge states (`plain`, a hash set):
  oriented depths 0-6 and T1 depths 0-6 equal the tables'; the Python sticker model
  (`p1py.py bfs`): oriented 0-5, T1 0-4.
* The oriented distribution over 12 equals the published held-edge distribution (every depth),
  and its slot counts the giant tier's table's histogram (every depth).
* Every state at the maximum has a solution read off its table: 4,932 of 15 moves (oriented) and
  1,276,200 of 14 (T1), each replayed on the Python sticker model to the target. With the T1
  table, every oriented state at 15 also has a 15-move solution ending in a class-A turn (a
  14-move path to T1, then the turn), replayed the same way.
* Lower bounds without a table (`mitm`): the classes within 8 moves of the targets by
  breadth-first search over classes (cube moves, the 288 pairs as class keys), and every move
  sequence from one state of each class at the maximum: none of the oriented 29 classes is
  within 14 moves (8 + 6) of the target, and none of T1's 4,670 classes at 14 is within 13
  (8 + 5).
* Slot by slot over all 3.34 x 10^11 slots (`cross`), against an exact strict table made with
  the two-phase index and moves (its states per distance equal the earlier strict
  distribution): oriented - strict is 0 to 4 everywhere (slots 1.71e10, 1.31e11, 1.62e11,
  2.35e10, 8.02e8), oriented - T1 is 0 or 1 everywhere (7.97e10, 2.55e11), T1 is never below
  strict, and the smallest oriented value and the smallest T1 value over a coset's 81 turns both
  equal the strict value on every coset. The oriented table mod 3 equals the giant tier's table on
  every slot.

## The radius-4 space around T1 (`chen-p2bfs --t1 --radius 4`)

Phase 1 to within 4 moves of T1 takes at most 14 - 4 = 10 moves; phase 2 is then an exact
search of every position whose configuration is within 4 moves of T1 (any frame), with every
move that stays in that set, so the bound is 10 + its deepest layer. `t1_ball` in
`src/chen/radius1.rs` builds the set like the radius-1 set: from the 108 fixed-frame T1
configurations (in the block split's space, with their frames) by breadth-first search with the
16 moves, each configuration's frame following its parent's. Checks: 146,841 configurations in
each of the 12 frames (the Python sticker model's ball: 9, 105, 1,193, 13,273, 146,841 within 0
to 4), one frame per configuration however it is reached, every configuration turned by its
frame in frame 0 at the same distance, frame 0 closed under the 24 symmetries, every move undone
by its inverse, no move changing the distance by 2. With the T1 table (`chen-oriented
ball-check`): every one of the 1,762,092 configurations' distance is its table value and every
one of the 21,635,712 moves out of the set reaches a value of 5.

The space: 146,841 configurations (6,134 classes) x 11,520 corner states x 369,600 B-centre
colourings = 6.25 x 10^14 positions, 2.61 x 10^13 slots (3.3 TB a dense layer); 3.72 moves per
configuration stay in the set on average. Positions per depth 0-7 (1, 16, 208, 2,688, 34,752,
278,456, 845,016, 6,406,223), split by distance to T1, agree between the search, a plain search
over whole cube states (`positions --t1-radius 4`) and an independent search on the Python
sticker model.

**Run** (started 2026-11-19 07:28 UTC with 36 GB for a part, resumed at depth 15 with 30 GB to
leave the machine 15 GB): `chen-p2bfs -j 48 run --t1 --radius 4 --ram-gb 30`,
layers in `/data4/oriented4/store`; a watcher (`keep-layers.py`) hard-links each completed layer
of at most 1.5 TB into `keep/` and keeps the 4 deepest, since the search deletes layer d - 2
while it makes layer d.

## Tables

| Table | File | Entries | Size | Max | Mean |
|---|---|---:|---:|---:|---:|
| AC: configuration x corners, pattern view, 24 symmetries (1,611 x 11,520) | `chen-acs.prun` | 18,558,720 | 4.6 MB | 13 | 10.05 |
| PX: pattern-view configuration x B-centres, 24 symmetries (1,611 classes) | `chen-px.prun` | 595,425,600 | 148.9 MB | 15 | 11.97 |
| CC: corners x B-centres, every move allowed (the optimal solver's, 24 symmetries) | `opt-cc-sym24.prun` | 18,849,600 | 4.7 MB | 11 | 8.35 |
| SB: pattern-view shape x B-centres (not used by default: PX is never lower) | `chen-sb.prun` | 170,016,000 | 42.5 MB | 13 | 9.72 |

AC, PX and SB hold distances with the allowed moves (mod 3), so a search decodes each child's
value from its parent's; the means are over table entries. The tier's two-phase tables come on
top (phase 1 and the in-memory `H_A` table).

## Search choices (node counts)

Measured on 4 seeded random positions (`--random 4 --seed 1`, small tier, one thread), as the
nodes to the first solution (`-n 1`; all but a few thousand are phase-2 nodes; on one thread
they repeat, up to the few thousand searched before the stop is seen):

| Choice | Nodes to the first solution | Lengths |
|---|---|---|
| PX, 2^18 nodes per phase-1 solution (default) | 1.7e6, 3.7e6, 8.2e6, 0.9e6 | 27, 28, 28, 28 |
| no limit per phase-1 solution | 1.1e7, none in 1.6e7, 1.7e7, 1.2e8 | 29, -, 29, 31 |
| 2^16 / 2^20 per phase-1 solution | 1.6e7 / 1.7e7 in all | 27-28 |
| SB instead of PX | 1.2e7, none in 1.5e7, none in 2.3e7, 6.9e6 | 25, -, -, 27 |
| an iteration's phase-1 solutions searched in order of their phase-2 bound | 1.5e7 in all | 27-28 |
| phase-2 children searched in order of their bound | 1.4e7 in all | 27-29 |

Timed (10 random positions, small tier, 10 s each), the limit of 2^16, 2^18 or 2^20 nodes gives
on one thread a first solution after 0.53, 0.76 and 1.17 s and 22.8, 22.9 and 22.8 moves at 10 s;
on 48 threads the first after 0.011, 0.0064 and 0.0065 s and the same lengths from 1 s on.

Without a limit per phase-1 solution, the first ones hold the search for long: their phase 2
is about 16 moves, while the tables say about 11 (on the first random position the leaves start
at a bound of 10 or 11, and the iteration at bound 16 needs 1-2 x 10^5 nodes).

## Against the two-phase solver

Measured on 2026-11-12, every run holding the exclusive bench lock (no other timed work) on the
48-core machine, with build 9dcb91c (the phase-2 prefetch and the overlapping depths, the
fixed-frame corner table; the later changes keep the node counts: the corner table under
symmetry costs 1.03x the time). The two-phase rows are the same build's `--method two-phase`.

### Random positions

`bench -N 10 --seed 1 -t 10`: the mean best length over the positions solved by then (how
many), and the mean first solution.

| Tier | Threads | Solver | First solution | 0.01 s | 0.1 s | 1 s | 10 s |
|---|---:|---|---|---|---|---|---|
| tiny (Light) | 1 | two-phase | 22.9 after 3.30 s (9 of 10) | - | - | 23.0 (3) | 22.7 (9) |
| tiny (Light) | 1 | block split | 27.5 after 0.75 s | - | - | 24.9 (7) | 24.2 (10) |
| tiny (Light) | 48 | two-phase | 23.3 after 0.095 s | - | 22.4 (5) | 22.1 | 21.5 |
| tiny (Light) | 48 | block split | 26.1 after 0.014 s | 25.6 (8) | 24.6 | 23.9 | 22.7 |
| small (Standard) | 1 | two-phase | 23.0 after 1.92 s | - | - | 22.6 (5) | 22.4 |
| small (Standard) | 1 | block split | 27.5 after 0.76 s | - | - | 24.8 (8) | 22.9 |
| small (Standard) | 48 | two-phase | 23.0 after 0.037 s | 23.5 (2) | 22.4 | 21.6 | 21.5 |
| small (Standard) | 48 | block split | 26.1 after 0.0064 s | 25.2 (9) | 23.7 | 22.6 | 21.7 |
| medium (Faster) | 1 | two-phase | 23.3 after 0.35 s | - | 23.0 (2) | 22.7 (9) | 21.7 |
| medium (Faster) | 1 | block split | 27.5 after 0.78 s | - | 25.0 (1) | 24.3 (7) | 22.6 |
| medium (Faster) | 48 | two-phase | 23.5 after 0.0081 s | 22.9 (7) | 22.1 | 21.5 | 20.8 |
| medium (Faster) | 48 | block split | 26.1 after 0.0070 s | 24.8 (8) | 22.9 | 22.0 | 21.5 |

(A dash: no position solved by then.) The block split finds a first solution sooner on the
small and tiny tiers (2.5x and 4.4x sooner on one thread, 5.8x and 6.8x on 48), but 3 to 5
moves longer. From 0.1 s (48 threads) or 1 s (one thread) on, the two-phase solver's solutions
are shorter at every tier, by 0.2 to 1.5 moves at 10 s (on one thread at 1 s it has solved
fewer of the positions on the small and tiny tiers: 5 and 3 of 10, against 8 and 7). With the
medium tier the two-phase solver's first solution comes about as soon (48 threads) or sooner
(one thread).

### Hard positions

`solve -t 20 -p improvements` on the 24-symmetric position (optimal 22) and on `bench/hard.tsv`
(known lengths in brackets; the corner 6-flip needs 20 moves on the normal FTO): first solution
and best within 20 s, as length @ seconds.

| Position | small 1T two-phase | small 1T block split | small 48T two-phase | small 48T block split | tiny 1T two-phase | tiny 1T block split |
|---|---|---|---|---|---|---|
| 24-symmetric (22) | none | none | 23 @ 5.3, 22 @ 18.3 | 25 @ 15.4 | none | none |
| tri_adj+edge_opp (23) | none | 29 @ 1.1, 26 @ 3.3 | 26 @ 3.1, 24 @ 7.4 | 29 @ 0.04, 26 @ 0.13 | none | 29 @ 1.2, 27 @ 1.3 |
| tri_adj+edge_adj (23) | none | none | 23 @ 2.9 | 29 @ 15.5, 26 @ 17.3 | none | none |
| tri_adj-[Uf_Fu;Db_Bd] (23) | 23 @ 0.68 | 29 @ 3.4, 23 @ 17.0 | 23 @ 0.017 | 29 @ 0.15, 23 @ 0.45 | 23 @ 2.4 | 29 @ 3.6, 26 @ 16.7 |
| tri_adj+[UB_FD,NL_MR] (23) | 25 @ 4.9 | 29 @ 8.8, 27 @ 19.9 | 25 @ 0.14, 23 @ 1.0 | 29 @ 0.18, 23 @ 8.9 | none | 29 @ 9.3, 28 @ 19.0 |
| edge_hard (21) | 23 @ 3.6, 22 @ 5.2 | 28 @ 0.23, 23 @ 1.9 | 23 @ 0.06, 21 @ 4.3 | 27 @ 0.007, 21 @ 7.9 | 23 @ 8.7, 22 @ 12.8 | 27 @ 0.33, 24 @ 7.1 |
| corner 6-flip (20) | none | none | 23 @ 8.7 | 24 @ 6.1, 23 @ 15.1 | none | none |

On one thread the block split finds a solution for more of these positions (small tier: 4 of
7 against 3; tiny tier: 4 against 2), always a long one; where both find one, the two-phase
solver's is as short or shorter. On 48 threads the two-phase solver solves all 7, and at 20 s
its solution is as short or shorter on every one (where they tie, it got there sooner); the
block split's first solution comes sooner on three of them. The 24-symmetric position and the
corner 6-flip are already in `H'` (phase-1 length 0), so the block split's whole work there is
phase 2 (the 24-symmetric position has the hexagons and corners solved and only the B-centres
to place), where its tables give lower bounds of only 11 and 9 moves.

### Memory and generation

| Table | Size | Generation, 48 threads | One thread |
|---|---:|---:|---:|
| AC (under symmetry) | 4.6 MB | 0.1 s | 1.4 s |
| PX | 148.9 MB | 1.4 s | 47.2 s |
| SB (optional) | 42.5 MB | 0.2 s | 8.3 s |
| CC (the optimal solver's) | 4.7 MB | | |

Building the configuration space at start-up takes 0.27 s (one thread). The default phase-2
tables (AC, PX, CC) are 158 MB on top of the tier: Light 42 MB becomes 200 MB, Standard 134 MB
becomes 292 MB. With SB instead of PX (52 MB; build 9f5c88f) the first-solution advantage is
gone (one thread, `FTO_CHEN_TABLES=sb`): small tier, first solution 26.2 moves after 3.22 s
(two-phase: 23.0 after 1.92 s), 23.6 moves at 10 s (22.4); tiny tier, 26.2 after 3.16 s
(two-phase: 22.9 after 3.30 s, 9 of 10 positions), 24.8 at 10 s (22.7).

## Verdict

The block split is not a better two-phase solver for this puzzle. Its phase 2 is a search over
1.59 x 10^14 positions in which the tables (the best pairwise projections that fit: AC, PX, CC)
fall about 5 moves short of the true distance (leaves at a bound of 10 or 11 whose phase 2 is 16
moves), so a phase-1 solution takes up to its 2.6 x 10^5 phase-2 nodes, mostly for nothing, and
the first solution of a random position about 10^6 to 10^7 nodes, where the two-phase solver's phase 2
is one table lookup. It gives a first solution sooner on the small and tiny tiers, and on one
thread it solves more of the hard positions within 20 s, but its solutions are 3 to 5 moves
longer at first; it needs 158 MB more memory for that, and from about 0.1 s (48 threads) or 1 s
(one thread) on its solutions are longer at every tier; on the hard positions that both solve,
its best is never shorter.
