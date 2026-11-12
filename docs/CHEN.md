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
  shortest allowed completion.

As in the two-phase solver, phase 1 targets both face classes, phase-1 depths grow by one, and
each phase-1 solution gets a phase-2 search for anything shorter than the best so far. A
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

## Tables

| Table | File | Entries | Size | Max | Mean |
|---|---|---:|---:|---:|---:|
| AC: configuration x corners | `chen-ac.prun` | 429,235,200 | 107.3 MB | 13 | 10.05 |
| PX: pattern-view configuration x B-centres, 24 symmetries (1,611 classes) | `chen-px.prun` | 595,425,600 | 148.9 MB | 15 | 11.97 |
| CC: corners x B-centres, every move allowed (the optimal solver's) | `opt-cc-sym24.prun` | 18,849,600 | 4.7 MB | 11 | |
| SB: pattern-view shape x B-centres (not used by default: PX is never lower) | `chen-sb.prun` | 170,016,000 | 42.5 MB | 13 | 9.72 |

AC, PX and SB hold distances with the allowed moves (mod 3), so a search decodes each child's
value from its parent's; the means are over table entries. The tier's two-phase tables come on
top (phase 1 and the in-memory `H_A` table).

## Search choices (node counts)

Measured on 4 seeded random positions (`--random 4 --seed 1`, small tier, one thread), as the
phase-2 nodes to the first solution (they repeat exactly on one thread):

| Choice | Phase-2 nodes to the first solution | Lengths |
|---|---|---|
| PX, 2^18 nodes per phase-1 solution (default) | 1.7e6, 3.7e6, 8.2e6, 0.9e6 | 27, 28, 28, 28 |
| no limit per phase-1 solution | 1.1e7, none in 1.6e7, 1.7e7, 1.2e8 | 29, -, 29, 31 |
| 2^16 / 2^20 per phase-1 solution | 1.6e7 / 1.7e7 in all | 27-28 |
| SB instead of PX | 1.2e7, none in 1.5e7, none in 2.3e7, 6.9e6 | 25, -, -, 27 |
| an iteration's phase-1 solutions searched in order of their phase-2 bound | 1.5e7 in all | 27-28 |
| phase-2 children searched in order of their bound | 1.4e7 in all | 27-29 |

Without a limit per phase-1 solution, the first ones hold the search for long: their phase 2
is about 16 moves, while the tables say about 11 (on the first random position the leaves start
at a bound of 10 or 11, and the iteration at bound 16 needs 1-2 x 10^5 nodes).
