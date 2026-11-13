//! The command line of the `omelettehedron` binary (not part of the library): one module per
//! command, and the parts they share.
//!
//! * [`args`]: the options several commands share (input, limits, tables, output), each set
//!   flattened into every command that takes it, so that an option has the same name and help
//!   everywhere.
//! * [`input`]: the positions to work on: scrambles, patterns or random positions.
//! * [`load`]: thread pools, the choice of tier and table loading.
//! * [`report`]: what a solve prints, the same for every method.
//! * [`solve`], [`gen`], [`bench`], [`show`], [`symmetric`]: the commands.

pub mod args;
pub mod bench;
pub mod gen;
pub mod input;
pub mod load;
pub mod report;
pub mod show;
pub mod solve;
pub mod symmetric;

/// The error of a command, printed as `error: ...`.
pub type Error = Box<dyn std::error::Error + Send + Sync>;
