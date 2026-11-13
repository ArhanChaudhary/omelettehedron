//! `omelettehedron gen two-phase | optimal`: generate (or check) the tables of a tier ahead of
//! time. Solves generate missing tables on first use too; `gen` reads every table of the tier,
//! so it also checks them.

use super::args::TableArgs;
use super::load::{optimal_tables, pool, with_two_phase};
use super::Error;
use omelettehedron::optimal::tables::OptTier;
use omelettehedron::tables::Tier;

/// The solver whose tables to generate.
#[derive(clap::Subcommand)]
pub enum Solver {
    /// Generate (or check) the two-phase solver's tables of a tier
    TwoPhase(TwoPhaseArgs),
    /// Generate (or check) the optimal solver's tables of a tier (`opt-*.prun`)
    Optimal(OptimalArgs),
}

/// `gen two-phase`.
#[derive(clap::Args)]
pub struct TwoPhaseArgs {
    /// The tier whose tables to generate
    #[arg(long, value_enum, help_heading = "Tables")]
    tier: Tier,
    #[command(flatten)]
    tables: TableArgs,
    #[command(flatten)]
    run: GenArgs,
}

/// `gen optimal`.
#[derive(clap::Args)]
pub struct OptimalArgs {
    /// The tier whose tables to generate
    #[arg(long, value_enum, help_heading = "Tables")]
    tier: OptTier,
    #[command(flatten)]
    tables: TableArgs,
    #[command(flatten)]
    run: GenArgs,
}

/// The options of both `gen` commands besides the tables.
#[derive(clap::Args)]
pub struct GenArgs {
    /// Number of worker threads (default: all cores)
    #[arg(short = 'j', long, help_heading = "Generation")]
    threads: Option<usize>,
    /// Print nothing but errors (without it: each table as it is loaded or generated, then the
    /// tables' depth histograms)
    #[arg(short, long, help_heading = "Output")]
    quiet: bool,
}

/// Run `gen`.
pub fn run(solver: Solver) -> Result<(), Error> {
    match solver {
        Solver::TwoPhase(a) => two_phase(&a),
        Solver::Optimal(a) => optimal(&a),
    }
}

/// `gen two-phase`.
fn two_phase(a: &TwoPhaseArgs) -> Result<(), Error> {
    let quiet = a.run.quiet;
    let p = pool(a.run.threads)?;
    p.install(|| {
        with_two_phase!(&a.tables, a.tier, quiet, |t| {
            if !quiet {
                eprintln!(
                    "tier {}, moves {}: AE {:.1} MB, CB {:.1} MB; phase-2 depths {:?}",
                    a.tier.name(),
                    t.moves.name(),
                    t.ae.bytes() as f64 / 1e6,
                    t.cb.tab.bytes() as f64 / 1e6,
                    t.p2.histogram()
                );
                eprintln!("AE depth histogram {:?}", t.ae.hist);
                eprintln!("CB depth histogram {:?}", t.cb.tab.hist);
                if let Some(cbe) = &t.cbe {
                    eprintln!("CBE depth histogram {:?}", cbe.hist);
                }
            }
            // (read it now: `gen` checks every table of the tier)
            if let (Some(lazy), Some(cbe)) = (&t.cbs, &t.cbe) {
                let Some(cbs) = lazy.get(&cbe.ix) else {
                    return Err(format!(
                        "{}: unusable distance-to-solved table",
                        lazy.path().display()
                    )
                    .into());
                };
                if !quiet {
                    eprintln!("distance-to-solved depth histogram {:?}", cbs.hist);
                }
            }
            Ok(())
        })
    })
}

/// `gen optimal`.
fn optimal(a: &OptimalArgs) -> Result<(), Error> {
    let quiet = a.run.quiet;
    let p = pool(a.run.threads)?;
    p.install(|| {
        let t = optimal_tables(&a.tables, a.tier, quiet)?;
        if quiet {
            return Ok(());
        }
        eprintln!(
            "optimal solver, tier {}, moves {}: {:.1} MB of tables",
            t.tier.name(),
            t.moves.name(),
            t.bytes() as f64 / 1e6
        );
        eprintln!("corners x centres depth histogram {:?}", t.cc.hist);
        if let Some(e) = &t.edges {
            eprintln!("edge permutation depth histogram {:?}", e.hist);
        }
        if let Some((_, ce)) = &t.ce {
            eprintln!("corners x edges depth histogram {:?}", ce.hist);
        }
        if let Some((_, cbe, _)) = &t.cbe {
            eprintln!(
                "corners x B-centres x edge triple depth histogram {:?}",
                cbe.hist
            );
        }
        eprintln!(
            "A-centre/edge (to H_A, two-phase table) depth histogram {:?}",
            t.ae.tab.hist
        );
        Ok(())
    })
}
