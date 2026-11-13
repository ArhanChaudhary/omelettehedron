//! `omelettehedron symmetric`: the positions fixed by symmetries of the octahedron (the notion
//! of `show`): the table of subgroups, the positions of some of them, or where given positions
//! are listed ([`omelettehedron::symsearch`]).

use super::args::InputArgs;
use super::input::{read_jobs, unread, Job};
use super::Error;
use omelettehedron::{pattern, symsearch};

/// `symmetric`.
#[derive(clap::Args)]
pub struct Args {
    /// List the positions of this subgroup (its number in the table)
    #[arg(long, value_name = "K", conflicts_with = "order")]
    subgroup: Option<usize>,
    /// List the positions with exactly this many symmetries (of every subgroup of that order)
    #[arg(long, value_name = "N")]
    order: Option<usize>,
    /// List at most this many positions of a subgroup (counted first): a subgroup with more is
    /// only counted
    #[arg(long, value_name = "N", default_value_t = 1_000_000)]
    max: u128,
    /// Say where positions are listed (their subgroup and number, e.g. 2.4): scrambles or
    /// patterns as arguments or one per line on stdin, or patterns with --patterns
    #[arg(long, conflicts_with_all = ["subgroup", "order"])]
    identify: bool,
    /// The positions for --identify: a scramble or a pattern
    #[arg(value_name = "SCRAMBLE|PATTERN", requires = "identify")]
    position: Vec<String>,
    /// --identify reads patterns of one or more lines from stdin (as `solve two-phase
    /// --patterns`)
    #[arg(long, requires = "identify")]
    patterns: bool,
}

/// Run `symmetric`.
pub fn run(a: &Args) -> Result<(), Error> {
    if a.identify {
        let input = InputArgs {
            position: a.position.clone(),
            inverse: false,
            patterns: a.patterns,
            random: None,
            seed: 1,
        };
        let jobs = read_jobs(&input, None)?;
        let census = symsearch::census();
        for entry in &jobs {
            match entry {
                Ok(j) => identify(j, &census, a.max),
                Err(bad) => bad.report(false),
            }
        }
        return unread(&jobs);
    }
    list(a.subgroup, a.order, a.max)
}

/// `symmetric --identify`: the subgroup of a position's symmetries and its number in the list
/// (if that subgroup has at most `max` positions up to symmetry, as when listing).
fn identify(j: &Job, census: &[(symsearch::SubgroupClass, symsearch::Counts)], max: u128) {
    let (k, rep) = symsearch::identify(&j.cube);
    println!("# {}", j.title());
    let Some((c, n)) = census.get(k) else { return };
    if c.order == 1 {
        println!("no symmetry: not listed");
        return;
    }
    let name = symsearch::schoenflies(c.rep);
    let at = if n.classes <= max {
        symsearch::exact(c.rep)
            .and_then(|l| l.iter().position(|x| *x == rep))
            .map_or_else(String::new, |i| format!(": position {}.{}", k + 1, i + 1))
    } else {
        String::new()
    };
    println!(
        "subgroup {}: {name}, {} symmetries ({}){at}",
        k + 1,
        c.order,
        symsearch::kinds_text(c.rep)
    );
}

/// The table of subgroups, or the positions of some of them.
fn list(subgroup: Option<usize>, order: Option<usize>, max: u128) -> Result<(), Error> {
    let census = symsearch::census();
    if subgroup.is_none() && order.is_none() {
        println!(
            "{:>2} {:>5} {:>4} {:>4} {:>24} {:>24} {:>24}  symmetries",
            "#", "order", "name", "conj", "fixed", "exactly", "up to symmetry"
        );
        for (i, (c, n)) in census.iter().enumerate() {
            println!(
                "{:>2} {:>5} {:>4} {:>4} {:>24} {:>24} {:>24}  {}",
                i + 1,
                c.order,
                symsearch::schoenflies(c.rep),
                c.conjugates,
                n.fixed,
                n.exact,
                n.classes,
                symsearch::kinds_text(c.rep)
            );
        }
        return Ok(());
    }
    if let Some(k) = subgroup {
        if k == 0 || k > census.len() {
            return Err(format!("--subgroup: give a number from 1 to {}", census.len()).into());
        }
    }
    let mut listed = 0;
    for (i, (c, n)) in census.iter().enumerate() {
        if subgroup.is_some_and(|k| k != i + 1) || order.is_some_and(|o| o != c.order) {
            continue;
        }
        let name = symsearch::schoenflies(c.rep);
        println!(
            "# subgroup {}: {name}, {} symmetries ({}): {} positions up to symmetry",
            i + 1,
            c.order,
            symsearch::kinds_text(c.rep),
            n.classes
        );
        if n.classes > max || c.order == 1 {
            println!("# (not listed: more than --max {max})");
            continue;
        }
        let list = symsearch::exact(c.rep).unwrap_or_default();
        if list.len() as u128 != n.classes {
            return Err(format!(
                "subgroup {}: {} positions listed, {} counted (please report this)",
                i + 1,
                list.len(),
                n.classes
            )
            .into());
        }
        for (j, cube) in list.iter().enumerate() {
            println!(
                "{}  # {name} {}.{}",
                pattern::format(&pattern::of_cube(cube)),
                i + 1,
                j + 1
            );
        }
        listed += 1;
    }
    if listed == 0 && subgroup.is_none() {
        if let Some(o) = order {
            if !census.iter().any(|(c, _)| c.order == o) {
                return Err(format!("--order: no subgroup has {o} symmetries").into());
            }
        }
    }
    Ok(())
}
