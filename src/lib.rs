#![deny(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::unreachable,
    clippy::todo,
    clippy::unimplemented
)]
#![cfg_attr(
    test,
    allow(
        clippy::unwrap_used,
        clippy::expect_used,
        clippy::panic,
        clippy::unreachable,
        clippy::todo,
        clippy::unimplemented
    )
)]
#![warn(clippy::pedantic)]
// Pedantic lints that would only add noise here:
#![allow(
    // narrowing and widening `as` casts are deliberate in the index/coordinate arithmetic
    clippy::cast_possible_truncation,
    clippy::cast_lossless,
    // statistics and timings are reported as f64
    clippy::cast_precision_loss,
    // names follow the maths (m = move, p = position, l = labels, x/y = positions, ...)
    clippy::many_single_char_names,
    clippy::similar_names,
    // almost every function is a pure helper; #[must_use] on all of them adds nothing
    clippy::must_use_candidate,
    // the hot-path helpers are force-inlined on purpose (measured)
    clippy::inline_always,
    // modules share one coordinate/geometry vocabulary through glob imports
    clippy::wildcard_imports,
    // the staged search loop and the table builders are long by design
    clippy::too_many_lines
)]
#![allow(
    clippy::needless_range_loop,
    clippy::too_many_arguments,
    clippy::new_without_default
)]
//! Two-phase solver for the face-turning octahedron.
pub mod bench;
pub mod cba;
pub mod cbe;
pub mod cbsym;
pub mod coord;
pub mod cube;
pub mod fast;
pub mod geometry;
pub mod mem;
pub mod notation;
pub mod prune;
pub mod search;
#[cfg(feature = "stats")]
pub mod stats;
pub mod sym;
pub mod tables;
