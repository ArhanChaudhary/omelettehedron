//! An exact breadth-first search of the block split's whole phase-2 space (the `chen-p2bfs`
//! binary): its distance distribution, its antipodes, and the greatest distance within `H'`.
//!
//! * [`instance`]: the slot index (configuration class x B-centre colouring x corner state, in
//!   the pattern view under the 24 symmetries) and its move tables, or a quotient of it;
//! * [`order`]: the orders of the classes that make the parts;
//! * [`lookup`]: positions looked up in the layers on disk;
//! * [`codec`]: the layer files (one per class, compressed blocks in checksummed chunks);
//! * [`engine`]: the search, one layer per step, a part of the next layer at a time in memory,
//!   resumable after any part;
//! * [`reference`]: independent searches to check it against;
//! * [`solve`]: positions of slots, and their shortest phase-2 solutions by IDA*;
//! * [`joint`]: the joint check of the relaxed phase 1 and the last phase-2 layer;
//! * [`radius`]: the radius-N phase-2 distance (blocks broken within N moves of a
//!   configuration), by IDA*.

pub mod codec;
pub mod engine;
pub mod instance;
pub mod joint;
pub mod lookup;
pub mod order;
pub mod radius;
pub mod reference;
pub mod solve;
