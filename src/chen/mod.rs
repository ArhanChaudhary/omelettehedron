//! A second two-phase solver with another split: the *block split*, proposed by Chen Shuang.
//!
//! Phase 1 only solves the A-centres and edges (each class-A face's centres and edges form an
//! intact hexagon, a *block*, turned by any multiple of 120°); phase 2 solves the rest with a
//! restricted move set that never breaks a colour group of those hexagons further, and ends in
//! the two-phase solver's `H_A` table. Phase 1 is short (at most 14 moves, 10.9 on average),
//! phase 2 is long. Separate from the two-phase solver of [`crate::search`]; it reads that
//! solver's A-centre/edge tables (phase 1) and `H_A` table, the optimal solver's corner x
//! B-centre table, and two tables of its own.
//!
//! * [`space`]: the block configurations and the move rule;
//! * [`tables`]: the phase-2 pruning tables;
//! * [`search`]: the two-phase loop.

pub mod search;
pub mod space;
pub mod tables;
