# Symmetric positions

How far from solved are the positions that symmetries of the octahedron fix? The aim is a
position that needs 24 moves (face-turn metric, centres of a colour alike). The known hardest
positions (23 moves, [bench/hard.tsv](../bench/hard.tsv)) are all symmetric. A symmetric
position also has few neighbours up to symmetry, so it is cheap to search: the optimal solver
searches it with a symmetric start, about as many times faster as the symmetries that fix it.

## The positions

A symmetry `M` fixes a position `P` if `M P M^-1 = P` up to a whole-puzzle turn (the notion of
`show` and of the optimal solver). The symmetries fixing a position form a subgroup of the 48.
`omelettehedron symmetric` enumerates, for each of the 33 subgroups up to conjugation, every
position it fixes: in the frame with the front corner held, the subgroup fixes the corner state
and permutes the edge places and centre places; an edge (centre colour) chosen at one place
fixes the edges (colours) on the whole orbit of that place, so a backtracking search over the
orbits lists them all (`src/symsearch.rs`). Conjugate positions are the same distance from
solved, so each set of conjugates is listed once.

Checks (unit tests): every listed position is fixed by its subgroup on the sticker model and
is reachable; the counts agree with a count from cycle types for the 47 cyclic subgroups and
with generate-and-test from a cyclic subgroup for the others; the count of all positions up to
symmetry by Burnside's lemma equals the sum of the classes; the 48- and 24-symmetric positions
of bench/hard.tsv are in the lists.

| # | order | group | symmetries | fixed | up to symmetry |
|---|---:|---|---|---:|---:|
| 1 | 48 | Oh | E 8C3 3C2 6m' 6S4 6C2' 6C4 i 8S6 3m | 8 | 8 |
| 2 | 24 | Th | E 8C3 3C2 i 8S6 3m | 32 | 12 |
| 3 | 24 | O | E 8C3 3C2 6C2' 6C4 | 16 | 4 |
| 4 | 24 | Td | E 8C3 3C2 6m' 6S4 | 16 | 4 |
| 5 | 16 | D4h | E 3C2 2m' 2S4 2C2' 2C4 i 3m | 256 | 248 |
| 6 | 12 | D3d | E 2C3 3m' 3C2' i 2S6 | 48 | 40 |
| 7 | 12 | T | E 8C3 3C2 | 384 | 84 |
| 8 | 8 | C4h | E C2 2S4 2C4 i m | 4,096 | 1,920 |
| 9 | 8 | D2d | E C2 2S4 2C2' 2m | 4,096 | 1,920 |
| 10 | 8 | C4v | E C2 2m' 2C4 2m | 2,048 | 896 |
| 11 | 8 | D2h | E C2 2m' 2C2' i m | 2,048 | 896 |
| 12 | 8 | D2h | E 3C2 i 3m | 32,768 | 5,332 |
| 13 | 8 | D4 | E 3C2 2C2' 2C4 | 4,096 | 1,916 |
| 14 | 8 | D2d | E 3C2 2m' 2S4 | 4,096 | 1,916 |
| 15 | 6 | S6 | E 2C3 i 2S6 | 15,552 | 7,740 |
| 16 | 6 | D3 | E 2C3 3C2' | 7,776 | 3,860 |
| 17 | 6 | C3v | E 2C3 3m' | 864 | 404 |
| 18 | 4 | C2v | E m' C2' m | 163,840 | 80,896 |
| 19 | 4 | C2h | E m' C2' i | 163,840 | 80,856 |
| 20 | 4 | C2v | E C2 2m | 1,310,720 | 318,080 |
| 21 | 4 | C2h | E C2 i m | 1,310,720 | 318,080 |
| 22 | 4 | C4 | E C2 2C4 | 983,040 | 243,328 |
| 23 | 4 | D2 | E C2 2C2' | 655,360 | 161,408 |
| 24 | 4 | S4 | E C2 2S4 | 3,145,728 | 783,488 |
| 25 | 4 | C2v | E C2 2m' | 524,288 | 129,152 |
| 26 | 4 | D2 | E 3C2 | 12,582,912 | 1,043,900 |
| 27 | 3 | C3 | E 2C3 | 136,048,896 | 34,006,116 |
| 28 | 2 | Cs | E m | 27,249,868,800 | 3,405,709,824 |
| 29 | 2 | Ci | E i | 136,249,344,000 | 5,676,851,876 |
| 30 | 2 | C2 | E C2' | 136,249,344,000 | 34,062,087,404 |
| 31 | 2 | Cs | E m' | 6,291,456,000 | 1,572,651,628 |
| 32 | 2 | C2 | E C2 | 603,979,776,000 | 75,494,920,960 |
| 33 | 1 | C1 | E | 31,408,133,379,194,880,000,000 | 654,336,112,006,441,760,108 |

"fixed": the positions the group fixes (it may fix more symmetries too); "up to symmetry": the
positions whose symmetries are exactly a copy of the group, counting turned or mirrored copies
of a position once. Kinds of symmetry: `C3` a third of a turn about a face axis; `C4`, `C2` a
quarter and a half turn about a corner axis; `C2'` a half turn about an edge axis; `i` the point
reflection; `m` a mirror through four corners, `m'` through two; `S4`, `S6` a turn and a mirror.
The 8 positions with all 48 symmetries combine corners solved or all flipped, edges solved or
one 48-symmetric edge pattern, and centres solved or one 48-symmetric centre pattern.

## Looking for 24

No position with at least 6 symmetries needs 24 moves: each of the 27,200 (one of each set of
conjugates, subgroups 1 to 17) has a solution of at most 23 moves. They are all in
[bench/symmetric.tsv](../bench/symmetric.tsv), and each one was replayed on its pattern exactly as
`symmetric --order N` lists it with an independent sticker model (a separate Python
implementation of the geometry), which found every one solved.

A position is ruled out by any solution of at most 23 moves, so each list was filtered with the
cheapest tool first (large tier, 48 threads):

1. the two-phase solver, `solve two-phase --batch -l 23` (one thread per position, stop at the first
   solution of at most 23 moves) with 5 s each, then 60 s and 600 s for the positions left;
2. for the rest, the optimal solver, `--start-depth 23 -l 23`: one iteration at depth 21 to
   order the tasks, then depth 23, which finds any solution of at most 23 moves (capped at
   20 minutes; a position with none would have needed at least 24, then an uncapped run and
   independent checks).

Orders 6 and 8 (26,800 positions): 9,772 were left after 5 s, 2,425 after 60 s and 207 after
600 s. The optimal solver found a solution of at most 23 moves for all 207, in at most 14
minutes each (11 hours on 48 threads in all). The two-phase solver is a weak filter for these
positions: in 600 s on one thread it missed solutions that the optimal solver finds in
seconds. Orders 12 to 48 (400 positions) were filtered first with `-l 22 -t 3` on 48 threads
per position; 37 were left, the optimal solver found a solution of at most 23 moves for 35 of
them, and the other 2 are known 23s ([bench/hard.tsv](../bench/hard.tsv)).

The two-phase solver's lengths can be far from optimal here: 4.1 (24 symmetries: corners, edges
and the A-centres solved, the B-centres in a pattern), for which it found nothing shorter than
22 moves, is exactly 19 (`solve optimal`, 0.06 s on 48 threads:
`B BR B' BL D U BL D' BL BR' D' BR' D BL' U B' BL' B U`).

Shortest solution found, by subgroup ("23": the first solution of at most 23 moves found; most
are not proven optimal):

| # | group | positions | at most 22 | 23 |
|---:|---|---:|---:|---:|
| 1 | Oh | 8 | 6 | 2 |
| 2 | Th | 12 | 8 | 4 |
| 3 | O | 4 | 1 | 3 |
| 4 | Td | 4 | 4 | 0 |
| 5 | D4h | 248 | 139 | 109 |
| 6 | D3d | 40 | 17 | 23 |
| 7 | T | 84 | 55 | 29 |
| 8 | C4h | 1,920 | 489 | 1,431 |
| 9 | D2d | 1,920 | 518 | 1,402 |
| 10 | C4v | 896 | 280 | 616 |
| 11 | D2h | 896 | 256 | 640 |
| 12 | D2h | 5,332 | 1,596 | 3,736 |
| 13 | D4 | 1,916 | 517 | 1,399 |
| 14 | D2d | 1,916 | 677 | 1,239 |
| 15 | S6 | 7,740 | 2,106 | 5,634 |
| 16 | D3 | 3,860 | 1,050 | 2,810 |
| 17 | C3v | 404 | 124 | 280 |

Positions of subgroups 18 to 33 (4 symmetries or fewer) are not covered: subgroups 18 to 26
alone have about 3.2 million positions up to symmetry, far more than a day of this machine at
the rates above.
