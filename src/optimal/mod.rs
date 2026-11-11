//! The optimal solver: single-phase IDA* straight to the solved state.
//!
//! Separate from the two-phase solver (`crate::search`, `crate::tables`, `crate::cbe`,
//! `crate::cbsym`), it shares only the foundations: the cubie model and moves
//! ([`crate::cube`]), coordinates ([`crate::coord`], [`crate::fast`]), the breadth-first table
//! generator and table files ([`crate::prune`]), memory and preflight checks, and notation.
//!
//! * [`sym`]: the 48 symmetries of the octahedron as they act on positions (they fix the solved
//!   state);
//! * [`coords`]: positions as byte vectors, and the indices of the tables;
//! * [`tables`]: the pattern databases (distance to solved of a projection) and the tiers;
//! * [`cbe`]: the large tier's corner x B-centre x edge-triple table;
//! * [`search`]: IDA* with canonical move order and a parallel task split;
//! * [`bench`]: the benchmark harness.
//!
//! See docs/OPTIMAL.md.

pub mod bench;
pub mod cbe;
pub mod coords;
pub mod search;
pub mod sym;
pub mod tables;
