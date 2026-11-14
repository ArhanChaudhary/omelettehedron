//! An exact breadth-first search of the block split's whole phase-2 space (the `chen-p2bfs`
//! binary): its distance distribution, its antipodes, and the greatest distance within `H'`.
//!
//! * [`instance`]: the slot index (configuration class x B-centre colouring x corner state, in
//!   the pattern view under the 24 symmetries) and its move tables, or a quotient of it;
//! * [`codec`]: the layer files (one per class, compressed blocks in checksummed chunks);
//! * [`engine`]: the search, one layer per step, a part of the next layer at a time in memory,
//!   resumable after any part;
//! * [`reference`]: independent searches to check it against;
//! * [`solve`]: positions of slots, and their shortest phase-2 solutions by IDA*.

pub mod codec;
pub mod engine;
pub mod instance;
pub mod reference;
pub mod solve;
