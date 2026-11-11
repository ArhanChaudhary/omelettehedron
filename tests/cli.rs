//! The command line: patterns read from stdin (`solve --patterns`, `optimal solve --pattern(s)`),
//! the tier chosen from the table files when `--tier` is not given, and `info`.
//!
//! The solves use the small tiers: the files in `$FTO_TABLES` if both small tiers are complete
//! there (nothing is ever written to it), otherwise tables generated once in the temporary
//! directory (`gen --tier small`, `optimal gen --tier small`).

use omelettehedron::cube::MoveSet;
use omelettehedron::optimal::tables::OptTier;
use omelettehedron::prune::{file_state, FileState, TableFile};
use omelettehedron::tables::Tier;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::OnceLock;

/// Exit status, stdout and stderr of the binary run with `args` and `stdin`.
fn run(args: &[&str], stdin: &str) -> (bool, String, String) {
    let mut child = Command::new(env!("CARGO_BIN_EXE_omelettehedron"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(stdin.as_bytes())
        .unwrap();
    let out = child.wait_with_output().unwrap();
    (
        out.status.success(),
        String::from_utf8(out.stdout).unwrap(),
        String::from_utf8(out.stderr).unwrap(),
    )
}

/// The small tiers' files (two-phase and optimal).
fn small_files() -> Vec<TableFile> {
    let mut v = Tier::Small.files(MoveSet::FACE);
    // (the A-centre/edge table is shared)
    for f in OptTier::Small.files(MoveSet::FACE).unwrap() {
        if !v.contains(&f) {
            v.push(f);
        }
    }
    v
}

fn complete(dir: &Path) -> bool {
    small_files()
        .iter()
        .all(|f| file_state(&dir.join(&f.name), f.entries) == FileState::Complete)
}

/// A directory with both small tiers.
fn small_tables() -> &'static Path {
    static DIR: OnceLock<PathBuf> = OnceLock::new();
    DIR.get_or_init(|| {
        if let Some(d) = std::env::var_os("FTO_TABLES").map(PathBuf::from) {
            if complete(&d) {
                return d;
            }
        }
        let d = std::env::temp_dir().join("omelettehedron-cli-small");
        if !complete(&d) {
            let dir = d.to_str().unwrap();
            for gen in [&["gen"][..], &["optimal", "gen"]] {
                let args: Vec<&str> = gen
                    .iter()
                    .copied()
                    .chain(["--tier", "small", "--tables", dir])
                    .collect();
                let (ok, _, err) = run(&args, "");
                assert!(ok, "{err}");
            }
        }
        d
    })
}

/// The pattern a scramble makes (`pattern`).
fn pattern_of(scramble: &str) -> String {
    let (ok, out, err) = run(&["pattern", scramble], "");
    assert!(ok, "{err}");
    out.trim().to_string()
}

/// Lines of `out` that follow a line starting with `head`, up to the next `# ` line.
fn section<'a>(out: &'a str, head: &str) -> Vec<&'a str> {
    out.lines()
        .skip_while(|l| !l.starts_with(head))
        .skip(1)
        .take_while(|l| !l.starts_with("# "))
        .collect()
}

/// The scramble printed after "scramble from your usual hold: " in `lines`.
fn usual_scramble<'a>(lines: &[&'a str]) -> &'a str {
    lines
        .iter()
        .find_map(|l| l.strip_prefix("scramble from your usual hold: "))
        .unwrap()
}

/// Three patterns on stdin, the second with a mistake (a sticker of the wrong colour), the third
/// over two lines: the good ones are solved (the scrambles printed make them), the bad one is
/// reported in plain words, and the run fails at the end. For both solvers.
#[test]
fn patterns_on_stdin_with_a_bad_one() {
    let dir = small_tables().to_str().unwrap();
    let a = pattern_of("R U' BL F D B'");
    let b = pattern_of("L BR' D");
    // a sticker entered with the wrong colour
    let good = pattern_of("U F");
    let wrong = if good.starts_with('W') { "R" } else { "W" };
    let bad = format!("{wrong}{}", &good[1..]);
    let (front, back) = b.split_once(" / ").unwrap();
    let stdin = format!("# three patterns\n{a}\n\n{bad}\n{front}\n{back}  # over two lines\n");
    for solver in [&["solve"][..], &["optimal", "solve"]] {
        let args: Vec<&str> = solver
            .iter()
            .copied()
            .chain([
                "--tier",
                "small",
                "--tables",
                dir,
                "-j",
                "2",
                "--patterns",
                "--usual",
                "WRGPBOYA",
            ])
            .collect();
        let (ok, out, err) = run(&args, &stdin);
        assert!(!ok, "{out}\n{err}");
        assert!(out.contains(&format!("# pattern 1 (line 2): {a}")), "{out}");
        assert!(out.contains("# pattern 2 (line 4)\n"), "{out}");
        assert!(out.contains(&format!("# pattern 3 (line 5): {b}")), "{out}");
        assert!(
            err.contains("error: pattern 2 (line 4): each colour should be on exactly 9 stickers"),
            "{err}"
        );
        assert!(err.contains("1 of 3 patterns could not be read"), "{err}");
        // the scrambles make the patterns (from the standard hold, which --usual gives)
        for (head, want) in [("# pattern 1 ", &a), ("# pattern 3 ", &b)] {
            let lines = section(&out, head);
            assert_eq!(lines.len(), 3, "{out}");
            assert_eq!(&pattern_of(usual_scramble(&lines)), want, "{out}");
        }
    }
    // one pattern to the optimal solver
    let (ok, out, err) = run(
        &[
            "optimal",
            "solve",
            "--tier",
            "small",
            "--tables",
            dir,
            "-j",
            "2",
            "--pattern",
            &a,
            "--usual",
            "WRGPBOYA",
        ],
        "",
    );
    assert!(ok, "{err}");
    let lines: Vec<&str> = out.lines().collect();
    assert_eq!(lines.len(), 3, "{out}");
    assert!(lines[1].starts_with("scramble from solved with "), "{out}");
    assert_eq!(pattern_of(usual_scramble(&lines)), a);
    // nothing on stdin, and a pattern list whose only pattern is wrong
    let (ok, _, err) = run(&["solve", "--tier", "small", "--patterns"], "# none\n");
    assert!(!ok && err.contains("no pattern on stdin"), "{err}");
    for solver in [&["solve"][..], &["optimal", "solve"]] {
        let args: Vec<&str> = solver
            .iter()
            .copied()
            .chain(["--tier", "small", "--tables", "/nonexistent", "--patterns"])
            .collect();
        let (ok, out, err) = run(&args, &bad);
        assert!(!ok);
        assert_eq!(out, "# pattern 1 (line 1)\n");
        assert!(
            err.contains("error: pattern 1 (line 1): each colour"),
            "{err}"
        );
        assert!(err.contains("1 of 1 patterns could not be read"), "{err}");
    }
}

/// Without `--tier`, the largest tier whose files are all there is used, and said on stderr;
/// with none, the run stops with a message and generates nothing.
#[test]
fn tier_is_chosen_from_the_files() {
    let dir = std::env::temp_dir().join(format!("omelettehedron-cli-tier-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let d = dir.to_str().unwrap();
    for solver in [&["solve"][..], &["optimal", "solve"]] {
        let args: Vec<&str> = solver
            .iter()
            .copied()
            .chain(["--tables", d, "U R"])
            .collect();
        let (ok, out, err) = run(&args, "");
        assert!(!ok, "{out}");
        // (the smallest tier of each solver)
        let smallest = if solver[0] == "optimal" {
            "small"
        } else {
            "tiny"
        };
        assert!(
            err.contains("no tier has all its tables")
                && err.contains(&format!("--tier {smallest}")),
            "{err}"
        );
    }
    assert_eq!(
        std::fs::read_dir(&dir).unwrap().count(),
        0,
        "nothing generated"
    );
    // `gen` still needs --tier
    let (ok, _, err) = run(&["gen", "--tables", d], "");
    assert!(!ok && err.contains("gen needs --tier"), "{err}");
    // the small tiers' files (and a medium-tier table cut short, which does not count)
    let small = small_tables();
    for f in small_files() {
        let from = small.join(&f.name);
        #[cfg(unix)]
        std::os::unix::fs::symlink(&from, dir.join(&f.name)).unwrap();
        #[cfg(not(unix))]
        std::fs::copy(&from, dir.join(&f.name)).unwrap();
    }
    std::fs::write(dir.join("ae-sym288-chiral.prun"), b"FTOPRUN2").unwrap();
    for solver in [&["solve"][..], &["optimal", "solve"]] {
        let args: Vec<&str> = solver
            .iter()
            .copied()
            .chain(["--tables", d, "-j", "2", "U R"])
            .collect();
        let (ok, out, err) = run(&args, "");
        assert!(ok, "{err}");
        assert_eq!(out.trim(), "R' U'");
        let line: Vec<&str> = err.lines().collect();
        assert_eq!(line.len(), 1, "{err}");
        assert!(
            line[0].starts_with("no --tier given: using --tier small")
                && line[0].contains("ae-sym288-chiral.prun is incomplete"),
            "{err}"
        );
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

/// `info` needs no table (the table directory does not even exist): the solved state, a
/// symmetric position of bench/hard.tsv and random positions, as text and as JSON.
#[test]
fn info_without_tables() {
    let info = |args: &[&str], stdin: &str| {
        let mut c = Command::new(env!("CARGO_BIN_EXE_omelettehedron"));
        c.arg("info")
            .args(args)
            .env("FTO_TABLES", "/nonexistent/omelettehedron-tables")
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped());
        let start = std::time::Instant::now();
        let mut child = c.spawn().unwrap();
        child
            .stdin
            .take()
            .unwrap()
            .write_all(stdin.as_bytes())
            .unwrap();
        let out = child.wait_with_output().unwrap();
        // well under a second when the machine is quiet; this only catches table loading
        assert!(start.elapsed().as_secs_f64() < 10.0);
        (
            out.status.success(),
            String::from_utf8(out.stdout).unwrap(),
            String::from_utf8(out.stderr).unwrap(),
        )
    };
    let (ok, out, err) = info(&[""], "");
    assert!(ok, "{err}");
    for want in [
        "pattern: WWWWWWWWW RRRRRRRRR GGGGGGGGG PPPPPPPPP / BBBBBBBBB OOOOOOOOO YYYYYYYYY AAAAAAAAA",
        "solved: yes",
        "symmetries: 48 of 48",
        "corners: 6 of 6 solved",
        "edges: 12 of 12 solved",
        "centres: 24 of 24",
        "each of the 24 is its own piece",
        "in H_A (solvable with U F BR BL alone): yes, in 0 moves",
        "lower bound: 0 moves",
    ] {
        assert!(out.contains(want), "{want}\n{out}");
    }

    // tri_adj+edge_opp of bench/hard.tsv: fixed by all 48 symmetries, corners solved
    let hard = include_str!("../bench/hard.tsv");
    let scramble = hard
        .lines()
        .find(|l| l.starts_with("tri_adj+edge_opp\t"))
        .unwrap()
        .split('\t')
        .nth(3)
        .unwrap();
    let (ok, out, err) = info(&[scramble], "");
    assert!(ok, "{err}");
    assert!(out.contains("solved: no"), "{out}");
    assert!(out.contains("symmetries: 48 of 48"), "{out}");
    assert!(out.contains("corners: 6 of 6 solved"), "{out}");
    assert!(
        out.contains("in H_A (solvable with U F BR BL alone): no"),
        "{out}"
    );
    // its pattern is the one `pattern` prints, and given as a pattern it says the same, but
    // cannot tell which centre is which
    let p = pattern_of(scramble);
    assert!(out.contains(&format!("pattern: {p}")), "{out}");
    let (ok, out2, err) = info(&["--pattern", &p], "");
    assert!(ok, "{err}");
    assert!(out2.contains("symmetries: 48 of 48"), "{out2}");
    assert!(out2.contains("not known"), "{out2}");

    // random positions (one per line of JSON), and scrambles on stdin
    let (ok, out, err) = info(&["--random", "3", "--seed", "7", "--json"], "");
    assert!(ok, "{err}");
    let rows: Vec<serde_json::Value> = out
        .lines()
        .map(|l| serde_json::from_str(l).unwrap())
        .collect();
    assert_eq!(rows.len(), 3);
    for (k, r) in rows.iter().enumerate() {
        assert_eq!(r["input"], format!("random position {k} (seed 7)"));
        assert_eq!(r["solved"], false);
        assert_eq!(r["symmetries"]["count"], 1);
        assert_eq!(r["pure"], serde_json::Value::Null);
        assert!(r["lower_bound"]["moves"].as_u64().unwrap() >= 1, "{r}");
        // 72 stickers, faces set apart, " / " between the halves
        assert_eq!(r["pattern"].as_str().unwrap().len(), 72 + 9);
    }
    let (ok, out, err) = info(&[], "U\n\nR' D\n");
    assert!(ok, "{err}");
    assert!(out.starts_with("# U\n"), "{out}");
    assert!(out.contains("\n# R' D\n"), "{out}");
    assert!(
        out.contains("in H_A (solvable with U F BR BL alone): yes, in 1 move\n"),
        "{out}"
    );
    assert!(out.contains("lower bound: 1 move ("), "{out}");

    // patterns on stdin, one with a mistake: reported in its place, and the run fails
    let stdin = format!("{p}\nWWW\n");
    let (ok, out, err) = info(&["--patterns", "--json"], &stdin);
    assert!(!ok);
    let rows: Vec<serde_json::Value> = out
        .lines()
        .map(|l| serde_json::from_str(l).unwrap())
        .collect();
    assert_eq!(rows[0]["symmetries"]["count"], 48);
    assert_eq!(rows[1]["input"], "pattern 2 (line 2)");
    assert!(rows[1]["error"].as_str().unwrap().contains("72 stickers"));
    assert!(err.contains("1 of 2 patterns could not be read"), "{err}");
}

/// Every way tables get generated checks first that they fit, and writes nothing when they do
/// not: `optimal gen` and the optimal solver's tables on first use (`--min-free-memory` asks for
/// more memory than any machine has), and the two-phase loader behind `gen`, `solve` and `bench`
/// (the same floor, set in this process).
#[cfg(target_os = "linux")]
#[test]
fn nothing_written_when_tables_do_not_fit() {
    let dir = std::env::temp_dir().join(format!("omelettehedron-cli-fit-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    let d = dir.to_str().unwrap();
    for args in [
        &["optimal", "gen", "--tier", "small"][..],
        &["optimal", "solve", "--tier", "small", "U R"],
    ] {
        let args: Vec<&str> = args
            .iter()
            .copied()
            .chain(["--tables", d, "--min-free-memory", "1e9"])
            .collect();
        let (ok, _, err) = run(&args, "");
        assert!(!ok && err.contains("not enough memory"), "{err}");
    }
    omelettehedron::preflight::set_memory_floor(u64::MAX / 4);
    for (tier, joint) in [(false, false), (true, true)] {
        let r = omelettehedron::tables::Tables::load_with(
            &dir,
            omelettehedron::tables::AeSym::new(tier),
            joint,
            MoveSet::FACE,
            false,
        );
        let e = r.err().unwrap();
        assert!(e.to_string().contains("not enough memory"), "{e}");
    }
    // (the check may create the directory, but no file in it)
    let files = std::fs::read_dir(&dir).map_or(0, Iterator::count);
    assert_eq!(files, 0);
    let _ = std::fs::remove_dir_all(&dir);
}
