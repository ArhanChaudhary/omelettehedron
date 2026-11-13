//! Thread pools, the choice of tier (`--tier`, or the largest whose tables are all there) and
//! table loading, for every command that loads tables. What is loaded or generated is said on
//! stderr unless the command is quiet.

use super::args::TableArgs;
use super::Error;
use omelettehedron::autotier;
use omelettehedron::cube::MoveSet;
use omelettehedron::optimal::tables::{OptTables, OptTier};
use omelettehedron::tables::Tier;
use std::path::Path;

#[cfg(not(target_family = "wasm"))]
pub fn pool(threads: Option<usize>) -> Result<rayon::ThreadPool, Error> {
    let mut b = rayon::ThreadPoolBuilder::new();
    if let Some(t) = threads {
        b = b.num_threads(t);
    }
    Ok(b.build()?)
}

/// WebAssembly (WASI) has no threads: the only worker is the current thread, for this pool and
/// for the global one (which table generation may use).
#[cfg(target_family = "wasm")]
pub fn pool(_threads: Option<usize>) -> Result<rayon::ThreadPool, Error> {
    // (fails harmlessly if the global pool already exists)
    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .use_current_thread()
        .build_global();
    Ok(rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .use_current_thread()
        .build()?)
}

/// A table error, with the directory.
pub fn table_error(dir: &Path, e: &std::io::Error) -> Error {
    format!("pruning tables in {}: {e}", dir.display()).into()
}

/// The two-phase tier: `tier`, or the largest whose tables are all in `dir` (said on stderr
/// unless `quiet`). `flag` is the option that sets it, for the messages.
pub fn two_phase_tier(
    tier: Option<Tier>,
    dir: &Path,
    moves: MoveSet,
    flag: &str,
    quiet: bool,
) -> Result<Tier, Error> {
    if let Some(t) = tier {
        return Ok(t);
    }
    let c = autotier::two_phase(dir, moves, flag)?;
    if !quiet {
        eprintln!("{}", c.message);
    }
    Ok(c.tier)
}

/// The optimal solver's tier, as [`two_phase_tier`].
pub fn optimal_tier(
    tier: Option<OptTier>,
    dir: &Path,
    moves: MoveSet,
    flag: &str,
    quiet: bool,
) -> Result<OptTier, Error> {
    if let Some(t) = tier {
        return Ok(t);
    }
    let c = autotier::optimal(dir, moves, flag)?;
    if !quiet {
        eprintln!("{}", c.message);
    }
    Ok(c.tier)
}

/// Load (generating if needed) the optimal solver's tables of `tier`, saying what is loaded or
/// generated on stderr unless `quiet`.
pub fn optimal_tables(tab: &TableArgs, tier: OptTier, quiet: bool) -> Result<OptTables, Error> {
    let dir = omelettehedron::tables::table_dir(tab.tables.as_deref());
    tab.apply_checks();
    OptTables::load(&dir, tier, tab.move_set(), !quiet).map_err(|e| table_error(&dir, &e))
}

/// Load (generating if needed) the two-phase tables of `$tier` with the options `$tab`, saying
/// what is loaded or generated on stderr unless `$quiet`, and evaluate `$body` with `$t` bound
/// to them (each tier has its own index type); table errors are returned from the enclosing
/// function.
macro_rules! with_two_phase {
    ($tab:expr, $tier:expr, $quiet:expr, |$t:ident| $body:expr) => {{
        use ::omelettehedron::light::AeLight;
        use ::omelettehedron::tables::{AeRelabel, AeSym, Tables, Tier};
        let tab: &$crate::cli::args::TableArgs = $tab;
        let tier: Tier = $tier;
        let verbose = !$quiet;
        let dir = ::omelettehedron::tables::table_dir(tab.tables.as_deref());
        let moves = tab.move_set();
        tab.apply_checks();
        let err = |e: std::io::Error| $crate::cli::load::table_error(&dir, &e);
        // The large tier needs 8 GB and the huge tier's index has 9.1e10 entries: 32-bit
        // targets such as WebAssembly (4 GB) support the tiny, small and medium tiers (the
        // medium index has 4.13e9 entries, below 2^32).
        if cfg!(target_pointer_width = "32") && matches!(tier, Tier::Large | Tier::Huge) {
            return Err(
                "only --tier tiny, small and medium are supported on 32-bit targets such \
                        as WebAssembly"
                    .into(),
            );
        }
        match tier {
            Tier::Tiny => {
                let $t = Tables::load(&dir, AeLight::new(), moves, verbose).map_err(err)?;
                $body
            }
            Tier::Small => {
                let $t = Tables::load(&dir, AeSym::new(false), moves, verbose).map_err(err)?;
                $body
            }
            Tier::Medium => {
                let $t = Tables::load(&dir, AeSym::new(true), moves, verbose).map_err(err)?;
                $body
            }
            Tier::Large => {
                let $t =
                    Tables::load_with(&dir, AeSym::new(true), true, moves, verbose).map_err(err)?;
                $body
            }
            Tier::Huge => {
                let $t =
                    Tables::load_with(&dir, AeRelabel::new(), true, moves, verbose).map_err(err)?;
                $body
            }
        }
    }};
}
pub(crate) use with_two_phase;

/// Evaluate `$body` with `$b` bound to the lower bounds of the optimal solver's tables `$t`
/// (each tier has its own type).
macro_rules! with_bounds {
    ($t:expr, |$b:ident| $body:expr) => {{
        use ::omelettehedron::optimal::tables::OptTier;
        let t: &::omelettehedron::optimal::tables::OptTables = $t;
        let missing = || -> $crate::cli::Error { "the loaded tables are not the tier's".into() };
        match t.tier {
            OptTier::Small => {
                let $b = t.small().ok_or_else(missing)?;
                $body
            }
            OptTier::Medium => {
                let $b = t.medium().ok_or_else(missing)?;
                $body
            }
            OptTier::Large => {
                let $b = t.large().ok_or_else(missing)?;
                $body
            }
        }
    }};
}
pub(crate) use with_bounds;
