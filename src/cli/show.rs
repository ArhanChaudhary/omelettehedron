//! `omelettehedron show`: what can be said about positions without loading any table: the
//! moves (for a scramble), the pattern, the symmetries (how many, and their group as listed by
//! `omelettehedron symmetric`), the pieces in place, membership of the two-phase solver's
//! subgroups and a lower bound ([`omelettehedron::info`]).

use super::args::InputArgs;
use super::input::{read_jobs, unread, Job};
use super::Error;
use omelettehedron::notation::pure_stickers;
use omelettehedron::{info, pattern, symsearch};

/// The long help of `show`: what a pattern is.
pub const LONG_ABOUT: &str = "\
Describe positions (scrambles or patterns) without loading any table

For each position: the moves (for a scramble), the pattern it makes, whether it is solved, its \
symmetries and their group (numbered as in `omelettehedron symmetric`), the solved pieces, \
whether it is in the two-phase solver's subgroups H_A and H_B, and a lower bound on its \
distance from solved. It takes well under a second.

A pattern is the 72 sticker colours, one character each (any 8 different characters; spaces, \
`/`, `,` and `|` are ignored). Hold the puzzle with a corner pointing at you and read the four \
front faces U, R, F, L (top, right, bottom, left: clockwise from the top), then turn the whole \
puzzle around, keeping the top at the top, and read the four back faces B, BL, D, BR (again \
clockwise from the top: B is now at the top and BL on the right). Read each face with its outer \
edge at the top: the 5 stickers along that edge from left to right, the 3 of the middle row, \
then the one at the corner in the middle of the picture. Printed patterns use the standard \
colour scheme held the standard way (white on U, green on F): W R G P / B O Y A (white red green \
purple / blue orange yellow grey) for the solved faces U R F L / B BL D BR. Patterns given to \
`solve` and `show` may use any letters, and the puzzle may be held any way and have any colour \
scheme.";

/// `show`.
#[derive(clap::Args)]
pub struct Args {
    #[command(flatten)]
    input: InputArgs,
    /// Print each position as one line of JSON
    #[arg(long, help_heading = "Output")]
    json: bool,
}

/// Run `show`.
pub fn run(a: &Args) -> Result<(), Error> {
    let jobs = read_jobs(&a.input, None)?;
    for entry in &jobs {
        match entry {
            Ok(j) => show(j, jobs.len() > 1, a.json)?,
            // (one line of JSON for it too, in its place)
            Err(bad) if a.json => {
                let out = serde_json::json!({ "input": bad.label, "error": bad.message });
                println!("{out}");
                eprintln!("error: {}: {}", bad.label, bad.message);
            }
            Err(bad) => bad.report(false),
        }
    }
    unread(&jobs)
}

/// The group of symmetries fixing a position: its number in `omelettehedron symmetric`'s table
/// (from 1), its name and its kinds of symmetry.
#[derive(serde::Serialize)]
struct Group {
    subgroup: usize,
    name: &'static str,
    kinds: String,
}

fn group(j: &Job) -> Group {
    let (k, _) = symsearch::identify(&j.cube);
    let stab = symsearch::stabiliser(&j.cube);
    Group {
        subgroup: k + 1,
        name: symsearch::schoenflies(stab),
        kinds: symsearch::kinds_text(stab),
    }
}

/// The number of turns in a scramble as written (`X2` is one; a rotation `Xv` is none).
fn turns(text: &str) -> usize {
    text.split(|c: char| c.is_whitespace() || c == ',')
        .filter(|t| {
            let name = t.trim_end_matches(|c: char| c == '\'' || c == '’' || c.is_ascii_digit());
            !name.is_empty() && !name.ends_with(['v', 'V'])
        })
        .count()
}

/// The facts about one position, as text or one line of JSON.
fn show(j: &Job, header: bool, json: bool) -> Result<(), Error> {
    let labels = j.moves.as_deref().map(pure_stickers);
    let facts = info::info(&j.cube, labels.as_ref());
    let pat = match (&j.pattern, &j.moves) {
        (Some(pi), _) => pattern::format(&pi.pattern),
        (None, Some(m)) => pattern::format(&pattern::of_moves(m)),
        (None, None) => pattern::format(&pattern::of_cube(&j.cube)),
    };
    let g = group(j);
    if json {
        #[derive(serde::Serialize)]
        struct Out<'a> {
            input: &'a str,
            moves: Option<&'a str>,
            pattern: &'a str,
            group: &'a Group,
            #[serde(flatten)]
            info: &'a info::Info,
        }
        let out = Out {
            input: j.title(),
            moves: j.moves.is_some().then_some(j.label.as_str()),
            pattern: &pat,
            group: &g,
            info: &facts,
        };
        println!("{}", serde_json::to_string(&out)?);
        return Ok(());
    }
    if header {
        println!("# {}", j.title());
    }
    if j.moves.is_some() {
        match turns(&j.label) {
            0 => println!("moves: none"),
            n => println!(
                "moves: {} ({n} move{})",
                j.label,
                if n == 1 { "" } else { "s" }
            ),
        }
    }
    println!("pattern: {pat}");
    for line in facts.lines() {
        let symmetries = line.starts_with("symmetries: ");
        println!("{line}");
        if symmetries {
            println!(
                "symmetry group: {}, subgroup {} in `omelettehedron symmetric` ({})",
                g.name, g.subgroup, g.kinds
            );
        }
    }
    Ok(())
}
