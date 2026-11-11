//! Checks before generating tables: fail at once, not after minutes of work, when the tables
//! cannot be held in memory or saved to disk. All tables of a run are checked together
//! before any is generated, so a later, larger table cannot fail after an earlier one has
//! already been built.

use std::io;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};

/// Head room kept free on disk and in memory beyond the tables themselves.
const MARGIN: u64 = 256 << 20;

static FORCE: AtomicBool = AtomicBool::new(false);

/// Skip the checks (`--force`).
pub fn set_force(force: bool) {
    FORCE.store(force, Ordering::Relaxed);
}

/// A table a run needs.
pub struct Planned {
    pub path: PathBuf,
    /// size in memory and on disk
    pub bytes: u64,
    /// extra memory needed while it is generated
    pub working_set: u64,
}

/// Bytes available to unprivileged users on the file system holding `dir`.
#[cfg(unix)]
fn free_disk(dir: &Path) -> io::Result<u64> {
    use std::os::unix::ffi::OsStrExt;
    let mut c = dir.as_os_str().as_bytes().to_vec();
    c.push(0);
    // SAFETY: statvfs is plain old data, for which all-zero bytes are a valid value
    let mut st: libc::statvfs = unsafe { std::mem::zeroed() };
    // SAFETY: `c` is a NUL-terminated path and `st` is a valid, writable statvfs
    let r = unsafe { libc::statvfs(c.as_ptr().cast(), &raw mut st) };
    if r != 0 {
        return Err(io::Error::last_os_error());
    }
    // (the field widths differ between platforms: u64 on Linux, u32 on macOS)
    #[allow(clippy::useless_conversion)]
    Ok(u64::from(st.f_bavail).saturating_mul(u64::from(st.f_frsize)))
}

/// Unknown elsewhere: do not block generation.
#[cfg(not(unix))]
#[allow(clippy::unnecessary_wraps)] // same signature as the Unix version
fn free_disk(_dir: &Path) -> io::Result<u64> {
    Ok(u64::MAX)
}

/// `MemAvailable` from /proc/meminfo, if it can be read.
fn available_memory() -> Option<u64> {
    let s = std::fs::read_to_string("/proc/meminfo").ok()?;
    let line = s.lines().find(|l| l.starts_with("MemAvailable:"))?;
    let kb: u64 = line.split_whitespace().nth(1)?.parse().ok()?;
    Some(kb * 1024)
}

/// Room left under this process's cgroup v2 memory limit, if it has one.
fn cgroup_room() -> Option<u64> {
    let cg = std::fs::read_to_string("/proc/self/cgroup").ok()?;
    let rel = cg.lines().find_map(|l| l.strip_prefix("0::"))?;
    let base = Path::new("/sys/fs/cgroup").join(rel.trim_start_matches('/'));
    // "max" (no limit) does not parse
    let max: u64 = std::fs::read_to_string(base.join("memory.max"))
        .ok()?
        .trim()
        .parse()
        .ok()?;
    let cur: u64 = std::fs::read_to_string(base.join("memory.current"))
        .ok()?
        .trim()
        .parse()
        .ok()?;
    Some(max.saturating_sub(cur))
}

fn gb(b: u64) -> f64 {
    b as f64 / 1e9
}

/// Make sure all tables of a run fit: the ones still to be generated (those whose files do
/// not exist yet) on disk, and every table plus the largest generation working set in memory.
///
/// # Errors
///
/// If a table directory cannot be created or inspected, or the disk or the free memory is too
/// small (unless the checks are switched off with [`set_force`]).
pub fn check_plan(tables: &[Planned]) -> io::Result<()> {
    let missing: Vec<&Planned> = tables.iter().filter(|t| !t.path.exists()).collect();
    if missing.is_empty() || FORCE.load(Ordering::Relaxed) {
        return Ok(());
    }
    let mut need_disk: Vec<(PathBuf, u64)> = vec![];
    for t in &missing {
        let dir = t.path.parent().unwrap_or(Path::new(".")).to_path_buf();
        std::fs::create_dir_all(&dir)?;
        match need_disk.iter_mut().find(|(d, _)| *d == dir) {
            Some((_, b)) => *b = b.saturating_add(t.bytes),
            None => need_disk.push((dir, t.bytes)),
        }
    }
    for (dir, bytes) in &need_disk {
        let disk = free_disk(dir)?;
        if disk < bytes.saturating_add(MARGIN) {
            return Err(io::Error::other(format!(
                "not enough disk space for the tables to generate in {}: they need {:.2} GB, \
                 {:.2} GB is free. Free some space or choose another directory with --tables \
                 (or skip this check with --force).",
                dir.display(),
                gb(bytes.saturating_add(MARGIN)),
                gb(disk)
            )));
        }
    }
    let resident: u64 = tables.iter().map(|t| t.bytes).sum();
    let working: u64 = missing.iter().map(|t| t.working_set).max().unwrap_or(0);
    let need_mem = resident.saturating_add(working).saturating_add(MARGIN);
    let mem = match (available_memory(), cgroup_room()) {
        (Some(a), Some(c)) => Some(a.min(c)),
        (a, c) => a.or(c),
    };
    if let Some(mem) = mem {
        if mem < need_mem {
            return Err(io::Error::other(format!(
                "not enough memory to generate the tables: this tier needs {:.2} GB while \
                 generating, {:.2} GB is available. Use a smaller --tier (or skip this check \
                 with --force).",
                gb(need_mem),
                gb(mem)
            )));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    // (without statvfs and /proc, e.g. on WebAssembly, the limits are unknown and nothing is
    // refused)
    #[cfg(unix)]
    #[test]
    fn rejects_tables_that_do_not_fit() {
        let dir = std::env::temp_dir().join("omelettehedron-preflight-test");
        let plan = |bytes| {
            vec![Planned {
                path: dir.join("missing.prun"),
                bytes,
                working_set: 0,
            }]
        };
        assert!(check_plan(&plan(1 << 20)).is_ok());
        let err = check_plan(&plan(u64::MAX / 4)).expect_err("an exabyte table must not fit");
        assert!(err.to_string().contains("not enough"));
        // a table whose file exists is loaded, not generated: nothing to check
        std::fs::create_dir_all(&dir).expect("temp dir");
        let present = dir.join("present.prun");
        std::fs::write(&present, b"x").expect("temp file");
        let loaded = [Planned {
            path: present,
            bytes: 1 << 20,
            working_set: u64::MAX / 4,
        }];
        assert!(check_plan(&loaded).is_ok());
    }
}
