//! The command line: what `solve` prints (the same for its three methods), patterns read from
//! stdin (`solve ... --patterns`) or given as arguments, the tier chosen from the table files
//! when `--tier` is not given, `show` and `symmetric`.
//!
//! The solves use the small tiers: the files in `$FTO_TABLES` if both small tiers are complete
//! there (nothing is ever written to it), otherwise tables generated once in the temporary
//! directory (`gen two-phase --tier small`, `gen optimal --tier small`).

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
            for gen in [&["gen", "two-phase"][..], &["gen", "optimal"]] {
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

/// The pattern a scramble makes (`show`).
fn pattern_of(scramble: &str) -> String {
    let (ok, out, err) = run(&["show", scramble], "");
    assert!(ok, "{err}");
    out.lines()
        .find_map(|l| l.strip_prefix("pattern: "))
        .unwrap()
        .to_string()
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

/// `solve <method>` on the small tiers with two threads, with `extra` options.
fn solve(method: &str, extra: &[&str], stdin: &str) -> (bool, String, String) {
    let dir = small_tables().to_str().unwrap();
    let mut args = vec![
        "solve", method, "--tier", "small", "--tables", dir, "-j", "2",
    ];
    if method == "combined" {
        args.extend(["--optimal-tier", "small"]);
    }
    args.extend_from_slice(extra);
    run(&args, stdin)
}

/// A line of the log of bounds (`    0.012 s  upper 21  two-phase  ...`).
fn is_bound(line: &str) -> bool {
    line.split_once(" s  ").is_some_and(|(t, rest)| {
        t.trim().parse::<f64>().is_ok()
            && (rest.starts_with("upper ") || rest.starts_with("lower "))
    })
}

/// The output of one position: bound lines, the solution, the result line (and, for a
/// pattern, the scramble lines after it). Returns the solution's length.
fn check_output(lines: &[&str], source: &str) -> usize {
    let n_bounds = lines.iter().take_while(|l| is_bound(l)).count();
    let rest = &lines[n_bounds..];
    assert!(rest.len() >= 2, "{lines:?}");
    let solution = rest[0];
    let n = solution.split_whitespace().count();
    assert!(
        rest[1].starts_with(&format!("{n} moves, optimal (solution: ")),
        "{lines:?}"
    );
    // the best upper bound is the solution, and the lower bound met it
    let uppers: Vec<&str> = lines[..n_bounds]
        .iter()
        .filter(|l| l.contains(" s  upper "))
        .copied()
        .collect();
    let last = uppers.last().unwrap();
    assert!(last.ends_with(solution), "{lines:?}");
    assert!(last.contains(&format!("upper {n:>2}  ")), "{lines:?}");
    assert!(
        lines[..n_bounds]
            .iter()
            .any(|l| l.contains(&format!("lower {n:>2}"))),
        "{lines:?}"
    );
    if source != "combined" {
        assert!(
            lines[..n_bounds].iter().all(|l| l.contains(source)),
            "{lines:?}"
        );
    }
    n
}

/// The three `solve` methods print alike: each new bound as it comes, then the solution and the
/// result line; with `-q` the solution alone; with `--batch` (two-phase, optimal) the same per
/// position, in input order. The tables loaded are said on stderr unless `-q`.
#[test]
fn the_three_methods_print_alike() {
    let scramble = "U R' BL F D B' L BR";
    let mut lengths = vec![];
    for method in ["two-phase", "optimal", "combined"] {
        let (ok, out, err) = solve(method, &[scramble], "");
        assert!(ok, "{err}");
        let lines: Vec<&str> = out.lines().collect();
        lengths.push(check_output(&lines, method));
        assert!(err.lines().any(|l| l.starts_with("loaded ")), "{err}");
        // -q: the solution alone, and nothing on stderr
        let (ok, quiet, err) = solve(method, &["-q", scramble], "");
        assert!(ok, "{err}");
        assert_eq!(err, "");
        assert_eq!(quiet.lines().count(), 1, "{quiet}");
        assert_eq!(quiet.split_whitespace().count(), lengths[0], "{quiet}");
    }
    assert!(lengths.iter().all(|&n| n == lengths[0]), "{lengths:?}");

    // several positions on stdin, one of them a pattern and one not a position at all
    let pattern = pattern_of("R U' BL F D B'");
    let stdin = format!("{scramble}\n# a comment\nU X\n{pattern}\n");
    for method in ["two-phase", "optimal"] {
        for batch in [false, true] {
            let extra: &[&str] = if batch { &["--batch"] } else { &[] };
            let (ok, out, err) = solve(method, extra, &stdin);
            assert!(!ok, "{out}");
            assert!(err.contains("error: U X: neither a scramble"), "{err}");
            assert!(err.contains("1 of 3 positions could not be read"), "{err}");
            let first = section(&out, &format!("# {scramble}"));
            assert_eq!(check_output(&first, method), lengths[0], "{out}");
            assert!(out.contains("\n# U X\n# "), "{out}");
            let third = section(&out, &format!("# {pattern}"));
            check_output(&third, method);
            assert!(
                third
                    .iter()
                    .any(|l| l.starts_with("scramble from solved with ")),
                "{out}"
            );
            // -q: one line per position read
            let mut args = vec!["-q"];
            args.extend_from_slice(extra);
            let (_, out, _) = solve(method, &args, &stdin);
            assert_eq!(out.lines().count(), 2, "{out}");
        }
    }
}

/// The user's case that printed nothing: a random position with no limit. The two-phase solver
/// improves step by step by default, so solutions come at once (with target-length passes the
/// first one would come only at the optimal length, minutes later).
#[test]
fn improvements_come_without_a_limit() {
    use std::io::BufRead;
    let dir = small_tables().to_str().unwrap();
    let mut child = Command::new(env!("CARGO_BIN_EXE_omelettehedron"))
        .args([
            "solve",
            "two-phase",
            "--tier",
            "small",
            "--tables",
            dir,
            "-j",
            "2",
            "--random",
            "1",
            "--seed",
            "123789",
            "-v",
            "--progress",
        ])
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .unwrap();
    let stdout = child.stdout.take().unwrap();
    let (tx, rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        for line in std::io::BufReader::new(stdout)
            .lines()
            .map_while(Result::ok)
        {
            if tx.send(line).is_err() {
                return;
            }
        }
    });
    // (the first solution takes milliseconds; the margin is for a loaded machine)
    let deadline = std::time::Instant::now() + std::time::Duration::from_secs(60);
    let mut seen = vec![];
    let found = loop {
        let left = deadline.saturating_duration_since(std::time::Instant::now());
        match rx.recv_timeout(left) {
            Ok(line) if line.contains(" s  upper ") => break true,
            Ok(line) => seen.push(line),
            Err(_) => break false,
        }
    };
    let _ = child.kill();
    let _ = child.wait();
    assert!(found, "no solution printed: {seen:?}");
    assert_eq!(
        seen.first().map(String::as_str),
        Some("# random position 0 (seed 123789)")
    );
}

/// Three patterns on stdin, the second with a mistake (a sticker of the wrong colour), the third
/// over two lines: the good ones are solved (the scrambles printed make them), the bad one is
/// reported in plain words, and the run fails at the end. For both solvers.
#[test]
fn patterns_on_stdin_with_a_bad_one() {
    let a = pattern_of("R U' BL F D B'");
    let b = pattern_of("L BR' D");
    // a sticker entered with the wrong colour
    let good = pattern_of("U F");
    let wrong = if good.starts_with('W') { "R" } else { "W" };
    let bad = format!("{wrong}{}", &good[1..]);
    let (front, back) = b.split_once(" / ").unwrap();
    let stdin = format!("# three patterns\n{a}\n\n{bad}\n{front}\n{back}  # over two lines\n");
    for method in ["two-phase", "optimal"] {
        let (ok, out, err) = solve(method, &["--patterns", "--usual", "WRGPBOYA"], &stdin);
        assert!(!ok, "{out}\n{err}");
        assert!(out.contains(&format!("# pattern 1 (line 2): {a}")), "{out}");
        assert!(out.contains("# pattern 2 (line 4)\n"), "{out}");
        assert!(out.contains(&format!("# pattern 3 (line 5): {b}")), "{out}");
        assert!(
            err.contains("error: pattern 2 (line 4): each colour should be on exactly 9 stickers"),
            "{err}"
        );
        assert!(err.contains("1 of 3 positions could not be read"), "{err}");
        // the scrambles make the patterns (from the standard hold, which --usual gives)
        for (head, want) in [("# pattern 1 ", &a), ("# pattern 3 ", &b)] {
            let lines = section(&out, head);
            check_output(&lines, method);
            assert_eq!(&pattern_of(usual_scramble(&lines)), want, "{out}");
        }
    }
    // one pattern as the argument, to the optimal solver
    let (ok, out, err) = solve("optimal", &[&a, "--usual", "WRGPBOYA"], "");
    assert!(ok, "{err}");
    let lines: Vec<&str> = out.lines().collect();
    check_output(&lines, "optimal");
    assert!(
        lines
            .iter()
            .any(|l| l.starts_with("scramble from solved with ")),
        "{out}"
    );
    assert_eq!(pattern_of(usual_scramble(&lines)), a);
    // nothing on stdin, and a pattern list whose only pattern is wrong
    let (ok, _, err) = run(
        &["solve", "two-phase", "--tier", "small", "--patterns"],
        "# none\n",
    );
    assert!(!ok && err.contains("no pattern on stdin"), "{err}");
    for method in ["two-phase", "optimal"] {
        let args = [
            "solve",
            method,
            "--tier",
            "small",
            "--tables",
            "/nonexistent",
            "--patterns",
        ];
        let (ok, out, err) = run(&args, &bad);
        assert!(!ok);
        assert_eq!(out, "# pattern 1 (line 1)\n");
        assert!(
            err.contains("error: pattern 1 (line 1): each colour"),
            "{err}"
        );
        assert!(err.contains("1 of 1 positions could not be read"), "{err}");
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
    for method in ["two-phase", "optimal"] {
        let (ok, out, err) = run(&["solve", method, "--tables", d, "U R"], "");
        assert!(!ok, "{out}");
        // (the smallest tier of each solver)
        let smallest = if method == "optimal" { "small" } else { "tiny" };
        assert!(
            err.contains("no tier has all its tables")
                && err.contains(&format!("omelettehedron gen {method} --tier {smallest}")),
            "{err}"
        );
    }
    assert_eq!(
        std::fs::read_dir(&dir).unwrap().count(),
        0,
        "nothing generated"
    );
    // `gen` needs --tier
    let (ok, _, err) = run(&["gen", "two-phase", "--tables", d], "");
    assert!(!ok && err.contains("--tier"), "{err}");
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
    for method in ["two-phase", "optimal"] {
        let (ok, out, err) = run(&["solve", method, "--tables", d, "-j", "2", "U R"], "");
        assert!(ok, "{err}");
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(check_output(&lines, method), 2, "{out}");
        assert!(lines.contains(&"R' U'"), "{out}");
        // the tier chosen, then the tables loaded
        let line: Vec<&str> = err.lines().collect();
        assert!(
            line[0].starts_with("no --tier given: using --tier small")
                && line[0].contains("ae-sym288-chiral.prun is incomplete"),
            "{err}"
        );
        assert!(line[1..].iter().all(|l| l.starts_with("loaded ")), "{err}");
        // -q: nothing but the solution
        let (ok, out, err) = run(
            &["solve", method, "--tables", d, "-j", "2", "-q", "U R"],
            "",
        );
        assert!(ok, "{err}");
        assert_eq!((out.as_str(), err.as_str()), ("R' U'\n", ""));
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

/// `show` needs no table (the table directory does not even exist): the solved state, a
/// symmetric position of bench/hard.tsv and random positions, as text and as JSON.
#[test]
fn show_without_tables() {
    let show = |args: &[&str], stdin: &str| {
        let mut c = Command::new(env!("CARGO_BIN_EXE_omelettehedron"));
        c.arg("show")
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
    let (ok, out, err) = show(&[""], "");
    assert!(ok, "{err}");
    for want in [
        "moves: none",
        "pattern: WWWWWWWWW RRRRRRRRR GGGGGGGGG PPPPPPPPP / BBBBBBBBB OOOOOOOOO YYYYYYYYY AAAAAAAAA",
        "solved: yes",
        "symmetries: 48 of 48",
        "symmetry group: Oh, subgroup 1 in `omelettehedron symmetric` (E 8C3 3C2 6m' 6S4 6C2' 6C4 i 8S6 3m)",
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
    let (ok, out, err) = show(&[scramble], "");
    assert!(ok, "{err}");
    // (the scramble starts with a whole-puzzle rotation, which is no move)
    assert!(
        out.starts_with(&format!("moves: {scramble} (23 moves)\n")),
        "{out}"
    );
    assert!(out.contains("solved: no"), "{out}");
    assert!(out.contains("symmetries: 48 of 48"), "{out}");
    assert!(out.contains("symmetry group: Oh, subgroup 1 "), "{out}");
    assert!(out.contains("corners: 6 of 6 solved"), "{out}");
    assert!(
        out.contains("in H_A (solvable with U F BR BL alone): no"),
        "{out}"
    );
    // given as its pattern it says the same, but cannot tell which centre is which
    let p = pattern_of(scramble);
    assert!(out.contains(&format!("pattern: {p}")), "{out}");
    let (ok, out2, err) = show(&[&p], "");
    assert!(ok, "{err}");
    assert!(!out2.contains("moves:"), "{out2}");
    assert!(out2.contains("symmetries: 48 of 48"), "{out2}");
    assert!(out2.contains("not known"), "{out2}");
    // neither a scramble nor a pattern
    let (ok, _, err) = show(&["U X"], "");
    assert!(
        !ok && err.contains("neither a scramble (unknown move `X`) nor a pattern"),
        "{err}"
    );

    // random positions (one per line of JSON), and scrambles on stdin
    let (ok, out, err) = show(&["--random", "3", "--seed", "7", "--json"], "");
    assert!(ok, "{err}");
    let rows: Vec<serde_json::Value> = out
        .lines()
        .map(|l| serde_json::from_str(l).unwrap())
        .collect();
    assert_eq!(rows.len(), 3);
    for (k, r) in rows.iter().enumerate() {
        assert_eq!(r["input"], format!("random position {k} (seed 7)"));
        assert_eq!(r["moves"], serde_json::Value::Null);
        assert_eq!(r["solved"], false);
        assert_eq!(r["symmetries"]["count"], 1);
        assert_eq!(r["group"]["name"], "C1");
        assert_eq!(r["group"]["subgroup"], 33);
        assert_eq!(r["pure"], serde_json::Value::Null);
        assert!(r["lower_bound"]["moves"].as_u64().unwrap() >= 1, "{r}");
        // 72 stickers, faces set apart, " / " between the halves
        assert_eq!(r["pattern"].as_str().unwrap().len(), 72 + 9);
    }
    let (ok, out, err) = show(&[], "U\n\nR' D  # two moves\n");
    assert!(ok, "{err}");
    assert!(out.starts_with("# U\nmoves: U (1 move)\n"), "{out}");
    assert!(out.contains("\n# R' D\n"), "{out}");
    assert!(
        out.contains("in H_A (solvable with U F BR BL alone): yes, in 1 move\n"),
        "{out}"
    );
    assert!(out.contains("lower bound: 1 move ("), "{out}");

    // patterns on stdin, one with a mistake: reported in its place, and the run fails
    let stdin = format!("{p}\nWWW\n");
    let (ok, out, err) = show(&["--patterns", "--json"], &stdin);
    assert!(!ok);
    let rows: Vec<serde_json::Value> = out
        .lines()
        .map(|l| serde_json::from_str(l).unwrap())
        .collect();
    assert_eq!(rows[0]["symmetries"]["count"], 48);
    assert_eq!(rows[1]["input"], "pattern 2 (line 2)");
    assert!(rows[1]["error"].as_str().unwrap().contains("72 stickers"));
    assert!(err.contains("1 of 2 positions could not be read"), "{err}");
}

/// Every way tables get generated checks first that they fit, and writes nothing when they do
/// not: `gen` and the tables of a solve on first use (`--min-free-memory` asks for more memory
/// than any machine has), and the two-phase loader behind `gen`, `solve` and `bench` (the same
/// floor, set in this process).
#[cfg(target_os = "linux")]
#[test]
fn nothing_written_when_tables_do_not_fit() {
    let dir = std::env::temp_dir().join(format!("omelettehedron-cli-fit-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    let d = dir.to_str().unwrap();
    for args in [
        &["gen", "optimal", "--tier", "small"][..],
        &["gen", "two-phase", "--tier", "small"],
        &["solve", "optimal", "--tier", "small", "U R"],
        &["solve", "two-phase", "--tier", "small", "U R"],
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

/// `symmetric`: the table of subgroups, and the positions of one listed as patterns that `show`
/// reads back with that many symmetries (no table needed).
#[test]
fn symmetric_positions() {
    let (ok, out, err) = run(&["symmetric"], "");
    assert!(ok, "{err}");
    // a header and the 33 classes, the full group first and the trivial one last
    assert_eq!(out.lines().count(), 34, "{out}");
    let first: Vec<&str> = out.lines().nth(1).unwrap().split_whitespace().collect();
    assert_eq!(&first[..7], ["1", "48", "Oh", "1", "8", "8", "8"], "{out}");
    let last: Vec<&str> = out.lines().last().unwrap().split_whitespace().collect();
    assert_eq!(
        &last[..5],
        ["33", "1", "C1", "1", "31408133379194880000000"]
    );

    // the 24-symmetry positions: 20 (three subgroups), each with 24 symmetries; read back as
    // patterns of one line each (the `#` comments are skipped) or with --patterns
    let (ok, out, err) = run(&["symmetric", "--order", "24"], "");
    assert!(ok, "{err}");
    let patterns: Vec<&str> = out.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(patterns.len(), 20, "{out}");
    for args in [&["show"][..], &["show", "--patterns"]] {
        let (ok, show, err) = run(args, &out);
        assert!(ok, "{err}");
        assert_eq!(show.matches("symmetries: 24 of 48").count(), 20, "{show}");
    }
    // tri_adj+edge_adj of bench/hard.tsv is one of them (up to symmetry: some listed
    // position has the same pattern up to turning, mirroring and renaming colours, so check
    // the count of one subgroup and that the scramble's pattern has 24 symmetries)
    assert!(out.contains("# subgroup 2: Th, 24 symmetries"), "{out}");
    let first = patterns[0].split('#').next().unwrap().trim();
    let (ok, show, err) = run(&["show", first], "");
    assert!(ok, "{err}");
    assert!(show.contains("symmetry group: Th, subgroup 2 "), "{show}");

    // a subgroup with too many positions is only counted
    let (ok, out, err) = run(&["symmetric", "--subgroup", "33"], "");
    assert!(ok, "{err}");
    assert!(out.contains("not listed"), "{out}");
    let (ok, _, err) = run(&["symmetric", "--order", "5"], "");
    assert!(!ok);
    assert!(err.contains("no subgroup has 5 symmetries"), "{err}");
}

/// `solve optimal --start-depth`: with `-l` of the same value, whether a position needs more
/// moves (one iteration); above the optimum, a solution that is not proven optimal.
#[test]
fn optimal_start_depth() {
    let scramble = "U R' BL F D B' L BR";
    let opt = |extra: &[&str]| {
        let mut args = extra.to_vec();
        args.push(scramble);
        let (ok, out, err) = solve("optimal", &args, "");
        assert!(ok, "{err}");
        out
    };
    let out = opt(&["-q"]);
    let optimum = out.split_whitespace().count();
    let below = (optimum - 1).to_string();
    let out = opt(&["--start-depth", &below, "-l", &below]);
    assert!(
        out.lines().last().unwrap().starts_with(&format!(
            "no solution found: none has fewer than {optimum} moves"
        )),
        "{out}"
    );
    assert_eq!(
        opt(&["-q", "--start-depth", &below, "-l", &below]),
        "no solution found\n"
    );
    let above = (optimum + 2).to_string();
    let out = opt(&["--start-depth", &above]);
    let last = out.lines().last().unwrap();
    assert!(
        last.contains("moves") && !last.starts_with("no solution"),
        "{out}"
    );
    let sol = opt(&["-q", "--start-depth", &above]);
    let n = sol.split_whitespace().count();
    assert!(n >= optimum && n <= optimum + 2, "{out}");
    let solved = pattern_of(&format!("{scramble} {}", sol.trim()));
    // (solved up to a whole-puzzle turn: every face one colour)
    let faces: Vec<&str> = solved.split_whitespace().filter(|f| *f != "/").collect();
    assert_eq!(faces.len(), 8, "{solved}");
    assert!(
        faces.iter().all(|f| f.chars().all(|c| f.starts_with(c))),
        "{solved}"
    );
    let (ok, _, err) = run(
        &[
            "solve",
            "optimal",
            "--start-depth",
            "9",
            "-l",
            "8",
            scramble,
        ],
        "",
    );
    assert!(!ok);
    assert!(err.contains("--start-depth"), "{err}");
}
