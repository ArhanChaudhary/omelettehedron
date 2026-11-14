# Coset solver for H_A

`fto-coset` (src/bin/fto-coset.rs, library `src/coset/`) follows the coset method that proved
the Rubik's Cube's diameter is 20 (Rokicki, Kociemba, Davidson and Dethridge, and their
`hcoset` program), with `H_A` = <U, F, BR, BL> as the subgroup. It was built to look for a
position that needs 24 moves in the cosets of `H_A` furthest from it.

**Result in one line:** the method is exact for the FTO with its B-centres told apart, but not
for the real FTO, whose same-coloured centres are interchangeable; there it gives upper bounds
only, and a "deep" coset is not made of deep positions. Details below.

## How it works

A coset is the 933,120 positions that one phase-1 sequence takes into `H_A` together. The
solver takes one position `y` of it and:

1. finds every phase-1 sequence of each length `n` that ends in a move outside `H_A` (IDA* with
   the two-phase solver's tables, cutting only on lower bounds, so none is missed);
2. marks the state of `H_A` that each one ends in (`y s`), at depth `n`: the coset's elements
   are indexed by these states, as in `hcoset`;
3. before each depth, the *prepass* applies each of the 8 class-A moves to the states marked
   at the depth before (a sequence ending in a class-A move is a shorter one plus that move).

The mark of an element is then the length of its shortest solution whose phase-1 part has been
enumerated: its exact distance if every length up to the mark was enumerated, an upper bound
otherwise. One byte per element (912 KB) holds the whole distribution; a prepass takes about
10 ms. Every step uses the ordinary right action of moves on positions (the derivation, and why
a left action would not do, are in `src/coset/mod.rs`).

```sh
fto-coset --tier large "<scramble or pattern>"          # phase 1 until all are within 23
fto-coset --tier large -S 19 --uncovered left.txt ...    # phase 1 up to 19 moves, list the rest
fto-coset --tier large --phase1-depth --random 50000      # distance to H_A of many positions
```

## Identical centres

In `H_A` the B-centres travel with the corners, so inside `H_A` a B-face's three centres
usually have three different colours. A permutation of same-coloured B-centres leaves the
solved puzzle as it is but not the other positions of `H_A`. So positions that share one
phase-1 sequence do not share the others, and marking every sequence that takes `y` into `H_A`
is wrong: in the coset of `R U L' B`, element 301457 got the mark 10, and the optimal solver
proves it needs 11.

With the 12 B-centres told apart (A-centres stay interchangeable: within `H_A` each A-face's
centres are one colour, so that does no harm) the cosets are real cosets and everything above
holds exactly. That is what `fto-coset` does: it labels the B-centres of `y` as the first
two-phase solution it finds solves them, and keeps only the sequences that take this labelled
position into `H_A` with every label in place. Its marks are distances in that labelled puzzle,
which are upper bounds for the real one.

The price: turns move pieces in 3-cycles, so a position has 6^4 / 2 = 648 labellings (lifts),
and a sequence that takes `y` into `H_A` as a pattern does so with the labels in place for only
one of them. The labelled coset gets about 1/648 of the sequences (measured: 1/400 to 1/580),
which takes about 2.5 more moves of phase 1 for the same coverage. A sequence that ends in
`H_A` up to a relabelling still serves the elements whose own B-colouring that relabelling
preserves (`Coset::improve`, on by default): on average 4.9 of the 648 relabellings, over the
960 corner states of `H_A`. And an element's own phase-1 sequences are mostly not `y`'s at
all, so its bound can stay well above its distance.

## Checks

* `tests/coset.rs`: the labelled B-centres follow the corners in all of `H_A`; the moves on the
  index are the class-A moves; the trivial coset gives the phase-2 distribution; on four
  cosets the marks equal a breadth-first search of the labelled puzzle to 6 moves, both ways;
  on a coset filled to 13, each bound equals the shortest solution through the sequences
  applied to the element itself, and the phase-2 length after each sequence equals the
  distance in `H_A` between the element's index and the sequence's end state; with the other
  relabellings' sequences, each improved bound equals the shortest solution through every
  sequence kept that takes the element itself into `H_A`; the batch phase-1 depth agrees
  with the full enumeration.
* `--check K` repeats the per-element comparison on K random elements of any coset (with the
  other relabellings' sequences too), and replays every solution on the sticker model: 300 of
  300 on each of three cosets. `--samples K` prints elements with their solutions, which an
  independent Python sticker model replays (300 of 300 on three cosets).
* Against the optimal solver, 100 random elements of the coset of `R U L' B D BL'` (exact
  labelled distances through 14): 99 bounds equal the optimal distance, one is 2 above (a
  shorter solution exists through another labelling of its B-centres).

## Measurements

Large tier, 48 threads.

**Distance to `H_A`** (`--phase1-depth`, the depth of the coset with the natural labelling):

| positions | <= 13 | 14 | 15 | 16 | 17 | 18 | mean |
|---|---:|---:|---:|---:|---:|---:|---:|
| 100,000 random | 450 | 5,132 | 41,702 | 52,453 | 263 | 0 | 15.47 |
| 100,000 of the 3.16 M order-4 symmetric positions | 1,163 | 3,942 | 22,929 | 61,140 | 10,825 | 1 | 15.76 |
| 27,200 symmetric (orders 6-48), `H_A` | 785 | 1,322 | 5,314 | 15,433 | 4,339 | 7 | 15.75 |
| the same, `H_B` | 822 | 1,324 | 5,360 | 15,403 | 4,286 | 5 | 15.74 |

Of the four
known 23s, two are 17 from both `H_A` and `H_B`, one 16 and one 15. The deepest coset found is
that of the symmetric position D4h 5.165, 18 from both.

**Enumeration cost**: about 1.5e9 nodes/s; each length costs about 13 times the one before,
whatever the coset's depth (random coset of depth 15: lengths 18, 19 in 8 s, 98 s; the depth-17
coset of a known 23: lengths 19, 20 in 67 s, 886 s).

**Coverage at 23** (labelled sequences only): the random coset after length 19, 84.5% of its
elements; the depth-17 coset after length 20, 17%. Covering a depth 17-18 coset at 23 this
way needs lengths 21-23: hours to weeks of the 48 cores.

**Are a deep coset's elements deep?** Their distance to `H_A` is the minimum over their 648
labellings, and only one of these is in the coset. Over 2,000 random elements: random coset
(depth 15) mean 14.90; the depth-17 coset 14-16 for 99.25% of them, mean 15.64, against 15.47
for random positions.

**Are a deep coset's elements hard?** They are the deep position with up to 11 class-A moves
made *before* its scramble, so some are. The two-phase solver (`-l 23 -t 3`) left 31 of 2,000
random elements of the depth-17 coset without a solution of at most 23 moves (none of 2,000 for
the random coset); with 60 s all 31 got one (22 or 23 moves).

**The deepest coset found**, that of D4h 5.165 (18 from `H_A`), after phase 1 to 19 moves
(86 s): exact labelled distances 18: 2, 19: 86; bounds 20: 621, 21: 3,618, 22: 19,965,
23: 96,490, 24: 322,669, 25: 409,730, 26: 79,831, 27: 108. So 120,782 elements are proven
within 23 (47,593 by the labelled sequences alone) and 812,338 are not covered. Of 5,000
random elements, the two-phase solver left 164 (3.3%) without a solution of at most 23 moves
within 3 s, more than for the depth-17 coset (1.55%).

## Cost of a coset at 23

| | random coset (depth 15) | depth 17-18 coset |
|---|---|---|
| labelled phase 1 to 19 | 107 s, 84.5% within 23 | 67-86 s, 12-17% within 23 |
| to 20 | ~15 min | 15 min (depth 17: 17% within 23) |
| to 21 / 22 / 23 (estimated, x13 each) | ~3 h / - / - | ~3 h / ~40 h / ~3 weeks |
| two-phase solver on every element (`-l 23 -t 3`, 0.8 CPU-s each) | ~4.3 h | ~4.3 h, plus ~1.5% at 60 s (~5 h) |

## Conclusion

The coset method of the cube's proof does not carry over to the FTO with interchangeable
centres: the 933,120 positions of a coset share their phase-1 sequences only as positions of
the puzzle with the B-centres told apart, where they are 2.5 moves harder, and a coset's
distance from `H_A` says little about its elements' distances from it. A coset can be covered
at 23 only for shallow cosets, in hours; for the deep cosets the per-element check by the
two-phase solver is cheaper than the enumeration, and that is plain sampling of positions
near a deep one. No position needing 24 moves was found.

