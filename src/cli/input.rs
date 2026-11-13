//! The positions a command works on: a scramble or a pattern as arguments, one per line on
//! stdin, patterns of one or more lines on stdin (`--patterns`), or seeded random positions
//! (`--random`). A position that cannot be read is kept as its error, reported in its turn; the
//! others are still worked on, and the command fails at the end.

use super::args::InputArgs;
use super::Error;
use omelettehedron::bench::random_positions;
use omelettehedron::cube::Cube;
use omelettehedron::notation::{
    format_moves, invert_moves, parse_moves, scramble_to_cube, Frame, PhysMove,
};
use omelettehedron::pattern;
use std::io::{BufRead, Read};

/// A position given as a pattern: its solution is followed by a scramble that makes it.
pub struct PatternInput {
    pub pattern: [char; pattern::STICKERS],
    /// the colours of `--usual`
    pub usual: Option<[char; 8]>,
}

/// A position to work on.
pub struct Job {
    /// printed as `# label` before the position's output when there are several
    pub label: String,
    /// a shorter name, if the label holds the pattern (`show` prints the pattern anyway)
    pub name: Option<String>,
    pub cube: Cube,
    pub frame: Frame,
    /// the moves that make it, if given as a scramble (inverted with `--inverse`)
    pub moves: Option<Vec<PhysMove>>,
    /// the pattern, if given as one
    pub pattern: Option<PatternInput>,
}

impl Job {
    /// What `# ...` lines call the position.
    pub fn title(&self) -> &str {
        self.name.as_deref().unwrap_or(&self.label)
    }
}

/// A position that could not be read: reported in its turn, and the run fails at the end.
pub struct BadInput {
    pub label: String,
    pub message: String,
}

impl BadInput {
    /// Report it: its `# label` line on stdout (unless `quiet`), the error on stderr.
    pub fn report(&self, quiet: bool) {
        if !quiet {
            println!("# {}", self.label);
        }
        eprintln!("error: {}: {}", self.label, self.message);
    }
}

/// A position or the reason it could not be read.
pub type Entry = Result<Job, BadInput>;

/// The colours of `--usual`.
pub fn parse_usual(usual: Option<&str>) -> Result<Option<[char; 8]>, Error> {
    let Some(u) = usual else { return Ok(None) };
    let v: Vec<char> = u
        .chars()
        .filter(|c| !c.is_whitespace() && *c != '/')
        .collect();
    let u: [char; 8] = v
        .as_slice()
        .try_into()
        .map_err(|_| "--usual: give the 8 colours of U R F L B BL D BR")?;
    Ok(Some(u))
}

/// The positions of `input`, in order: `--patterns` (stdin), `--random`, the arguments (one
/// position: the arguments joined), or one position per line of stdin. A scramble is inverted
/// with `--inverse`; `usual` goes with the positions given as patterns.
///
/// # Errors
///
/// If stdin cannot be read, `--usual` is not 8 colours, or there is no position at all.
pub fn read_jobs(input: &InputArgs, usual: Option<&str>) -> Result<Vec<Entry>, Error> {
    let usual = parse_usual(usual)?;
    if input.patterns {
        let mut text = String::new();
        std::io::stdin().lock().read_to_string(&mut text)?;
        let records = pattern::records(&text);
        if records.is_empty() {
            return Err("no pattern on stdin".into());
        }
        return Ok(records
            .iter()
            .enumerate()
            .map(|(k, r)| {
                let name = format!("pattern {} (line {})", k + 1, r.line);
                match read_pattern(&r.text, usual) {
                    Ok(mut job) => {
                        job.label = format!("{name}: {}", job.label);
                        job.name = Some(name);
                        Ok(job)
                    }
                    Err(e) => Err(BadInput {
                        label: name,
                        message: e.to_string(),
                    }),
                }
            })
            .collect());
    }
    if let Some(n) = input.random {
        return Ok(random_positions(n, input.seed)
            .into_iter()
            .enumerate()
            .map(|(k, cube)| {
                Ok(Job {
                    label: format!("random position {k} (seed {})", input.seed),
                    name: None,
                    cube,
                    frame: Frame::identity(),
                    moves: None,
                    pattern: None,
                })
            })
            .collect());
    }
    let texts: Vec<String> = if input.position.is_empty() {
        std::io::stdin()
            .lock()
            .lines()
            .map_while(Result::ok)
            .map(|l| l.split('#').next().unwrap_or("").trim().to_string())
            .filter(|l| !l.is_empty())
            .collect()
    } else {
        vec![input.position.join(" ")]
    };
    if texts.is_empty() {
        return Err("no scramble given".into());
    }
    Ok(texts
        .into_iter()
        .map(|text| {
            read_position(&text, input.inverse, usual).map_err(|message| BadInput {
                label: text,
                message,
            })
        })
        .collect())
}

/// One position given as text: a scramble (inverted if `inverse`) or else a pattern.
fn read_position(text: &str, inverse: bool, usual: Option<[char; 8]>) -> Result<Job, String> {
    let move_error = match parse_moves(text) {
        Ok(m) => {
            let (m, label) = if inverse {
                let m = invert_moves(&m);
                let label = format_moves(&m);
                (m, label)
            } else {
                (m, text.to_string())
            };
            let (cube, frame) = scramble_to_cube(&m);
            return Ok(Job {
                label,
                name: None,
                cube,
                frame,
                moves: Some(m),
                pattern: None,
            });
        }
        Err(e) => e,
    };
    match pattern::parse(text) {
        // 72 stickers: meant as a pattern
        Ok(_) if inverse => Err("--inverse takes scrambles, not patterns".into()),
        Ok(_) => read_pattern(text, usual).map_err(|e| e.to_string()),
        Err(e) => Err(format!(
            "neither a scramble ({move_error}) nor a pattern ({e})"
        )),
    }
}

/// A position given as a pattern.
fn read_pattern(text: &str, usual: Option<[char; 8]>) -> Result<Job, pattern::PatternError> {
    let p = pattern::parse(text)?;
    let cube = pattern::position(&p)?;
    Ok(Job {
        label: pattern::format(&p),
        name: None,
        cube,
        frame: Frame::identity(),
        moves: None,
        pattern: Some(PatternInput { pattern: p, usual }),
    })
}

/// After the positions: an error if some could not be read (each was reported in its turn).
///
/// # Errors
///
/// If any position could not be read.
pub fn unread(jobs: &[Entry]) -> Result<(), Error> {
    match jobs.iter().filter(|j| j.is_err()).count() {
        0 => Ok(()),
        bad => Err(format!(
            "{bad} of {} positions could not be read (see above)",
            jobs.len()
        )
        .into()),
    }
}

/// When no position could be read: report each, and fail (before any table is loaded).
///
/// # Errors
///
/// Always (when some position could not be read, which is the only case it is called in).
pub fn none_read(jobs: &[Entry], quiet: bool) -> Result<(), Error> {
    for bad in jobs.iter().filter_map(|j| j.as_ref().err()) {
        bad.report(quiet);
    }
    unread(jobs)
}
