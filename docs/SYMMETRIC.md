# Symmetric positions

How far from solved are the positions that symmetries of the octahedron fix? The aim is a
position that needs 24 moves (face-turn metric, centres of a colour alike). The known hardest
positions (23 moves, [bench/hard.tsv](../bench/hard.tsv)) are all symmetric. A symmetric
position also has few neighbours up to symmetry, so it is cheap to search: the optimal solver
searches it with a symmetric start, about as many times faster as the symmetries that fix it.

## The positions

A symmetry `M` fixes a position `P` if `M P M^-1 = P` up to a whole-puzzle turn (the notion of
`info` and of the optimal solver). The symmetries fixing a position form a subgroup of the 48.
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
