//! The tier to use when `--tier` is not given: the largest one whose table files are all in
//! the table directory and complete. Files are judged from their sizes and headers only
//! ([`crate::prune::file_state`]): no table is read, and nothing is generated.

use crate::cube::MoveSet;
use crate::optimal::tables::OptTier;
use crate::prune::{file_state, FileState, TableFile};
use crate::tables::Tier;
use std::path::Path;

/// The chosen tier and the line that says why.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Choice<T> {
    pub tier: T,
    pub message: String,
}

/// The first problem among `files` in `dir` (`None` if they are all there and complete).
fn problem(dir: &Path, files: &[TableFile]) -> Option<String> {
    files
        .iter()
        .find_map(|f| match file_state(&dir.join(&f.name), f.entries) {
            FileState::Complete => None,
            FileState::Missing => Some(format!("{} is missing", f.name)),
            FileState::Bad(why) => Some(format!(
                "{} is incomplete or not this version's ({why})",
                f.name
            )),
        })
}

/// The largest of `tiers` (smallest first: name, files, or `None` if the tier cannot be used
/// here) whose files are all complete in `dir`. `flag` is the option that sets the tier (e.g.
/// `--tier`), and `gen` gives the command that generates the tables of a tier (by name), for
/// the message when there is none.
fn pick<T: Copy>(
    dir: &Path,
    tiers: &[(T, &str, Option<Vec<TableFile>>)],
    flag: &str,
    gen: &dyn Fn(&str) -> String,
) -> Result<Choice<T>, String> {
    // why each larger tier was not chosen, largest first
    let mut larger: Vec<String> = vec![];
    for (tier, name, files) in tiers.iter().rev() {
        let Some(files) = files else { continue };
        match problem(dir, files) {
            None => {
                // (the next larger tier's reason is the one that matters)
                let why = larger.last().map_or(String::new(), |l| format!(" ({l})"));
                return Ok(Choice {
                    tier: *tier,
                    message: format!(
                        "no {flag} given: using {flag} {name}, the largest tier whose tables are \
                         all in {}{why}",
                        dir.display()
                    ),
                });
            }
            Some(p) => larger.push(format!("{name}: {p}")),
        }
    }
    let smallest = tiers
        .iter()
        .find(|t| t.2.is_some())
        .map_or("small", |t| t.1);
    Err(format!(
        "no {flag} given, and no tier has all its tables in {} ({}). Pass {flag} (missing tables \
         are then generated), or generate them first with `{}` (or a larger tier)",
        dir.display(),
        larger.last().map_or("", String::as_str),
        gen(smallest)
    ))
}

/// The command (`omelettehedron gen two-phase` or `omelettehedron gen optimal`: `solver`) that
/// generates the tables of a tier in the metric of `moves`.
fn gen_command(solver: &'static str, moves: MoveSet) -> impl Fn(&str) -> String {
    let slice = if moves.slice { " --moves slice" } else { "" };
    move |tier| format!("omelettehedron gen {solver} --tier {tier}{slice}")
}

/// The largest two-phase tier whose tables (in the metric of `moves`) are all complete in
/// `dir`. On 32-bit targets only the small and medium tiers count. `flag` is the option that
/// sets the tier, for the messages.
///
/// # Errors
///
/// If no tier is complete: the message says what to do.
pub fn two_phase(dir: &Path, moves: MoveSet, flag: &str) -> Result<Choice<Tier>, String> {
    let tiers: Vec<(Tier, &str, Option<Vec<TableFile>>)> = Tier::ALL
        .iter()
        .map(|&t| {
            let usable =
                !(cfg!(target_pointer_width = "32") && matches!(t, Tier::Large | Tier::Huge));
            (t, t.name(), usable.then(|| t.files(moves)))
        })
        .collect();
    pick(dir, &tiers, flag, &gen_command("two-phase", moves))
}

/// The largest tier of the optimal solver whose tables (in the metric of `moves`) are all
/// complete in `dir` (the large tier has no slice tables). `flag` is the option that sets the
/// tier, for the messages.
///
/// # Errors
///
/// If no tier is complete: the message says what to do.
pub fn optimal(dir: &Path, moves: MoveSet, flag: &str) -> Result<Choice<OptTier>, String> {
    let tiers: Vec<(OptTier, &str, Option<Vec<TableFile>>)> = OptTier::ALL
        .iter()
        .map(|&t| (t, t.name(), t.files(moves)))
        .collect();
    pick(dir, &tiers, flag, &gen_command("optimal", moves))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    /// A file with the header of a complete table of `entries` entries and the right length
    /// (sparse: nothing but the header is written, so even the huge tier costs no disk).
    fn fake(path: &Path, entries: u64) {
        let mut f = std::fs::File::create(path).unwrap();
        f.write_all(b"FTOPRUN2").unwrap();
        f.write_all(&entries.to_le_bytes()).unwrap();
        f.write_all(&0u64.to_le_bytes()).unwrap();
        f.set_len(32 + entries.div_ceil(32) * 8).unwrap();
    }

    fn fake_tier(dir: &Path, files: &[TableFile]) {
        for f in files {
            fake(&dir.join(&f.name), f.entries);
        }
    }

    const SLICE: MoveSet = MoveSet {
        wide: false,
        slice: true,
    };

    // (files on a real file system: not on WebAssembly)
    #[cfg(unix)]
    #[test]
    fn picks_the_largest_complete_tier() {
        let dir =
            std::env::temp_dir().join(format!("omelettehedron-autotier-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        // nothing yet: an error that says what to do
        let e = two_phase(&dir, MoveSet::FACE, "--tier").unwrap_err();
        // (the smallest tier, quickest to generate)
        assert!(
            e.contains("--tier") && e.contains("omelettehedron gen two-phase --tier tiny"),
            "{e}"
        );
        let e = optimal(&dir, SLICE, "--tier").unwrap_err();
        assert!(
            e.contains("omelettehedron gen optimal --tier small --moves slice"),
            "{e}"
        );

        fake_tier(&dir, &Tier::Tiny.files(MoveSet::FACE));
        let c = two_phase(&dir, MoveSet::FACE, "--tier").unwrap();
        assert_eq!(c.tier, Tier::Tiny);
        assert!(
            c.message
                .contains("small: ae-sym288-triples.prun is missing"),
            "{}",
            c.message
        );
        fake_tier(&dir, &Tier::Small.files(MoveSet::FACE));
        let c = two_phase(&dir, MoveSet::FACE, "--tier").unwrap();
        assert_eq!(c.tier, Tier::Small);
        assert!(c.message.contains("--tier small"), "{}", c.message);
        assert!(
            c.message
                .contains("medium: ae-sym288-chiral.prun is missing"),
            "{}",
            c.message
        );
        // the slice tables are other files
        assert!(two_phase(&dir, SLICE, "--tier").is_err());
        // the wide move set reads the face-turn tables
        let wide = MoveSet {
            wide: true,
            slice: false,
        };
        assert_eq!(two_phase(&dir, wide, "--tier").unwrap().tier, Tier::Small);

        // the large tier's files, but the joint table cut short: medium (whose files are a
        // subset of large's) is the largest complete one
        fake_tier(&dir, &Tier::Large.files(MoveSet::FACE));
        let cbe = dir.join("cbe.prun");
        let len = std::fs::metadata(&cbe).unwrap().len();
        std::fs::OpenOptions::new()
            .write(true)
            .open(&cbe)
            .unwrap()
            .set_len(len - 8)
            .unwrap();
        let c = two_phase(&dir, MoveSet::FACE, "--tier").unwrap();
        assert_eq!(c.tier, Tier::Medium);
        assert!(
            c.message.contains("cbe.prun is incomplete"),
            "{}",
            c.message
        );
        fake(&cbe, crate::tables::CBE_ENTRIES);
        assert_eq!(
            two_phase(&dir, MoveSet::FACE, "--tier").unwrap().tier,
            Tier::Large
        );
        // a table of another size (another version's index) does not count
        fake(&cbe, crate::tables::CBE_ENTRIES + 32);
        assert_eq!(
            two_phase(&dir, MoveSet::FACE, "--tier").unwrap().tier,
            Tier::Medium
        );
        fake(&cbe, crate::tables::CBE_ENTRIES);
        // a temporary file of an unfinished generation is not the table
        std::fs::rename(dir.join("cbe-solved11.prun"), dir.join("cbe-solved11.tmp")).unwrap();
        assert_eq!(
            two_phase(&dir, MoveSet::FACE, "--tier").unwrap().tier,
            Tier::Medium
        );
        std::fs::rename(dir.join("cbe-solved11.tmp"), dir.join("cbe-solved11.prun")).unwrap();

        if cfg!(target_pointer_width = "64") {
            assert_eq!(
                two_phase(&dir, MoveSet::FACE, "--tier").unwrap().tier,
                Tier::Large
            );
            fake_tier(&dir, &Tier::Huge.files(MoveSet::FACE));
            let c = two_phase(&dir, MoveSet::FACE, "--tier").unwrap();
            assert_eq!(c.tier, Tier::Huge);
            assert!(!c.message.contains('('), "{}", c.message);
        }

        // the optimal solver's tiers: the two-phase small tier gives it one file of three
        let e = optimal(&dir, MoveSet::FACE, "--tier").unwrap_err();
        assert!(e.contains("opt-cc-sym24.prun is missing"), "{e}");
        for t in OptTier::ALL {
            fake_tier(&dir, &t.files(SLICE).unwrap_or_default());
        }
        // (the large tier has no slice tables: medium is the largest with slices)
        assert_eq!(
            optimal(&dir, SLICE, "--tier").unwrap().tier,
            if cfg!(target_pointer_width = "64") {
                OptTier::Medium
            } else {
                OptTier::Small
            }
        );
        assert!(optimal(&dir, MoveSet::FACE, "--tier").is_err());
        fake_tier(&dir, &OptTier::Small.files(MoveSet::FACE).unwrap());
        assert_eq!(
            optimal(&dir, MoveSet::FACE, "--tier").unwrap().tier,
            OptTier::Small
        );
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
