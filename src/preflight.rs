//! Checks before generating a table: fail at once, not after minutes of work, when the table
//! cannot be held in memory or saved to disk.

use std::io;
use std::path::Path;

/// Head room kept free on disk and in memory beyond the table itself.
const MARGIN: u64 = 256 << 20;

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
    // (the field widths differ between platforms)
    Ok(u64::from(st.f_bavail).saturating_mul(u64::from(st.f_frsize)))
}

/// Unknown elsewhere: do not block generation.
#[cfg(not(unix))]
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

fn gb(b: u64) -> f64 {
    b as f64 / 1e9
}

/// Make sure a table of `bytes` can be generated in memory and saved at `path`.
pub fn check_table(path: &Path, bytes: u64) -> io::Result<()> {
    let dir = path.parent().unwrap_or(Path::new("."));
    std::fs::create_dir_all(dir)?;
    let disk = free_disk(dir)?;
    if disk < bytes.saturating_add(MARGIN) {
        return Err(io::Error::other(format!(
            "not enough disk space to save {}: it needs {:.2} GB, {} has {:.2} GB free. \
             Free some space or choose another directory with --tables.",
            path.display(),
            gb(bytes.saturating_add(MARGIN)),
            dir.display(),
            gb(disk)
        )));
    }
    if let Some(mem) = available_memory() {
        if mem < bytes.saturating_add(MARGIN) {
            return Err(io::Error::other(format!(
                "not enough memory to generate {}: it needs {:.2} GB, {:.2} GB is available. \
                 Use a smaller --tier.",
                path.display(),
                gb(bytes.saturating_add(MARGIN)),
                gb(mem)
            )));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_tables_that_do_not_fit() {
        let dir = std::env::temp_dir().join("fto-preflight-test");
        let path = dir.join("t.prun");
        assert!(check_table(&path, 1 << 20).is_ok());
        let err = check_table(&path, u64::MAX / 4).expect_err("an exabyte table must not fit");
        assert!(err.to_string().contains("not enough"));
    }
}
