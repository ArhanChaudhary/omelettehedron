# Upper bounds through the A-block set: why the planned coset method fails

**Status: no bound is claimed.** The plan was to bound God's number by 14 + T, where T is the
largest distance of a position whose class-A hexagons are intact, with a coset ("block") method
over that set. The method's central step does not hold for this puzzle, because centres of a
colour are identical. This page records the failure, the checks behind it, the structure that
does hold, and what a correct version would cost. The checks are `tests/bound.rs`.

## The sets

* `H_A` = ⟨U, F, BR, BL⟩, the phase-2 group: 933,120 patterns, every corner state (960) times
  972 edge states. Its A-centres follow from the edges, and its B-centres copy the corner
  sticker next to them.
  - Checked: with every centre labelled it still has 933,120 elements. So it acts faithfully
    on corners and edges, and as a group it is C × E_A.
* `H'`, the relaxed set: positions whose A-centre/edge (AE) coordinate is solved.
  - Checked: these are exactly "each A hexagon (three A-centres and three edges) a block on its
    own face, any block rotation, up to whole-puzzle rotation". The 972 A-centre/edge
    configurations of `H_A` are the 12 rotations times the 81 hexagon rotations.
  - So `H'` = 972 configurations × 960 corner states × 369,600 B-colourings = 3.45 × 10¹¹
    positions (front corner held).
  - The medium AE table (`ae-sym288-chiral.prun`) has these as its only target (one index
    entry, closed under its 288 symmetries).
  - Its distances equal a brute force that knows only the 972 configurations, on 60 positions
    (distances 0–12; sequences of up to 7 moves against a ball of radius 5 around `H'`).
  - Entries per distance (symmetry-reduced index entries, not positions): 1, 1, 4, 50, 524,
    5849, 63563, 686316, 7385536, 77576093, 703614668, 2707776227, 630579953, 4014, 1 at
    0–14.

## The planned method

A block is the set of positions with one phase-1 coordinate. A phase-1 walk `s0` of the block
maps it onto `H_A` by `y ↦ y·s0`. Any other walk `s` then acts as `y·s = (y·s0)·g` with
`g = s0⁻¹s`. On corners and edges this is the right translation by `h_g ∈ C × E_A`. The
method marks `h_g⁻¹` at length `|s|` and grows the covered set by left multiplication with
A-moves.

## Where it fails

`y·s = (y·s0)·g` is right on corners and edges, which are distinct pieces. But it also
assumes that `y·s` is in `H_A`, i.e. that its B-centres copy its corner stickers. `g` can
permute identical B-centres:

* It maps `y0·s0` to `y0·s` for the position `y0` the walks were found from.
* For the other elements of the block it can take them out of `H_A`.

**Smallest example.** Take the block of the solved position (`s0` empty; the block is the
`H_A` set itself).

* `W = B F' BR D U' BL' U' BL D' U BR' U F B'` (14 moves) 3-cycles the three centres of face R
  and changes nothing else. So it keeps the solved pattern: it is a walk of this block, with
  `h_W` = identity.
* The method therefore credits every element `z` of `H_A` with the solution `W` + (phase 2 of
  `z`).
* For `z` = U that solution is `W U'`, which does not solve U. `U·W` keeps U's corners and
  edges, but its R centres no longer copy their corners, so it is not in `H_A`.
* 9,372 of 9,620 sampled elements of `H_A` leave `H_A` under `W`.
* There is such a 3-cycle for every B face (16–18 internal moves each).

**Real walks.** For random positions `y0` of `H'` we took all phase-1 walks of lengths d1 and
d1 + 1 (exhaustive search), with `s0` the first.

* For each other walk `s` and 20 random elements `y` of the planned block, we tested whether
  `y·s` lands in `H_A`. It mostly does not: 252 of 260, 1,248 of 1,260, 1,637 of 1,800 and 529
  of 760 land outside, on four positions with d1 = 11, 12, 12 and 12. (A fifth position has
  no walk of at most 13 moves; see Step 1 below.)

* **False certificates.** The method marks one element of the block as solved by a walk `s`
  alone. On the checked blocks that element is not solved by `s`. Exhaustive search finds no
  solution of `|s|` moves for it at all. On three of the four positions the method would credit
  it with 12, 13 and 13 moves; it needs more (its best solutions through walks of that length
  have 16, 17 and 17 moves). On the fourth, the pair of walks happens to agree.
* **The phase-1 distance is not constant on a planned block.** Four random elements of each:
  10, 10, 11, 10 (d1 of `y0` = 11); 7, 12, 10, 11; 10, 11, 11, 10; 12, 9, 11, 11.

## The structure that does hold

**Walks depend only on corners and B-centres.** Take a position of `H'`. Whether a walk takes
it into `H_A` depends only on its corner state and B-colouring (its *CB value*). Which of the
972 A-centre/edge configurations it has does not matter.

* Why: the AE and CB coordinates move independently, and the 972 configurations share their
  walks (`H_A` maps the A-centre and edge places of each A face to those of an A face).
* Checked: every walk of the sampled positions also works for 27 other configurations each.

**Different CB values never share all their walks.**

1. Suppose positions `y` and `y'` of `H'` have the same walks. Take a walk `s` of `y`, and
   the A-moves `v` that solve `y·s`.
2. Then the solved position and `u = y'·s·v` (in `H_A`) have the same walks.
3. The four B-face 3-cycles are walks of the solved position, so they keep `u` in `H_A`. That
   needs each B face of `u` to be one colour.
4. Only the solved corner state does that (checked over all 960). So `u` has the solved CB
   value, and then `y'` has the CB value of `y`.

So a block whose elements share their walks has at most 972 elements. `H'` has 354,816,000
such blocks, at least 1.48 × 10⁷ (354,816,000 / 24) up to the 24 symmetries (the CB table has
18,849,600 entries).

## Consequences

* **Step 1 (Q = max over `H'` of the distance to `H_A`)** needs every CB value. The planned
  30,800 representatives with solved corners are not enough: the distance depends on the
  corner state too (it varies within a planned block, above), and a position with solved
  corners stands only for its own CB value.
  * Distances to `H_A` of 32 random CB values (with their walks): 7 (1 value), 9 (2), 10 (2),
    11 (6), 12 (11), 13 (10).
  * **A position of `H'` needs 14 moves to reach `H_A`.** It is
    `OOOOGOOOW PRRYBBBBR AAAAAAAAG BYPYBRRPP / YPYPPYBRR WWWWWWWWA RBYRYBPPY OGGGOGGGG`
    (pattern format of `--pattern`), found among random positions by the real-walk test.
    - Exhaustive search finds no walk of at most 13 moves, with the CB distances taken from an
      independent breadth-first search over all 354,816,000 (corner state, B-colouring) pairs
      (`deep_phase1_position`; that search agrees with the CB table on 20,000 random pairs).
    - A search that uses neither table finds none either: CB distances from that independent
      search, and A-centre/edge distances only from an exact ball of radius 5 around `H'`
      (with at most 5 moves left; `deep_phase1_position_without_ae_table`, 8.7 × 10⁸ nodes
      at 13 moves).
    - It finds 48 walks of 14 moves, e.g. `B U' R L BR' F' BL' F BR BL L' R' U B'`.
    - So Q ≥ 14, more than the 13 moves the planned bound used for this step, and
      14 + Q + 11 ≥ 39.
    - Its distance to `H_B` (the class-B group, through the conjugation by the class-swapping
      rotation) is also 14, so reaching either group takes 14 moves.
* **Per-element bounds over `H'`** with valid blocks: one exhaustive walk enumeration per CB
  value, 1.48 × 10⁷ of them. Sample of 32 random blocks, walks of lengths d1 and d1 + 1
  (`block_bounds_sample`):

  | d1 | blocks | mean nodes | largest bound |
  |---:|---:|---:|---:|
  | 7 | 1 | 3.1e4 | 15 |
  | 9 | 2 | 5.7e5 | 18 |
  | 10 | 2 | 2.4e6 | 19 |
  | 11 | 6 | 1.5e7 | 20 |
  | 12 | 11 | 1.0e8 | 21 |
  | 13 | 10 | 5.5e8 | 22 |

  * Mean: 2.1 × 10⁸ nodes per block. 1.48 × 10⁷ blocks × 2.1 × 10⁸ ≈ 3 × 10¹⁵ nodes.
  * At the two-phase solver's 1.67 × 10⁹ nodes/s (48 threads, idle machine) that is about 22
    days. The sample's own search is simpler and ran at about 10⁷ nodes/s on the shared
    machine.
  * The bounds over the 31,104 sampled elements: at most 22 (69 elements at 22, 1,172 at 21).
  * A bound of 23 alone needs less, because phase 2 needs at most 11 moves. When d1 ≤ 12, any
    one walk of at most 12 moves bounds the whole block by 23. A walk of 13 moves covers the
    elements whose phase-2 length after it is at most 10, so blocks with d1 = 13 need a few
    walks. Walks of 14 cover phase-2 lengths up to 9 only (blocks with d1 = 14 exist, above).
* **Labelling the B-centres** makes the method valid, but then it bounds a harder puzzle:
  every B-centre must return to its own place. By counting, that costs about
  log₁₂.₉₃ 1,296 ≈ 2.8 moves (1,296 = (3!)⁴ labelled arrangements per solved pattern).
* **A decomposition without this problem:** phase 1 solves the edges and one centre orbit
  exactly. That is a pattern coordinate with a single target value, and the phase-2 bounds are
  then per element. This is the decomposition of the published proof of 36
  (github.com/antipodeaddict/ftoupper).

## Reproducing

* `cargo test --release --test bound`: structure and the smallest example (no tables).
* With the medium tier in `FTO_TABLES`: `cargo test --release --test bound -- --ignored
  --nocapture`. This runs:
  - the AE brute force;
  - the real walks;
  - a check of the walk search against an unpruned search;
  - the block sample (about 10 minutes at 48 threads).
* Logs of the runs above are kept with the scratch files of the computation.
