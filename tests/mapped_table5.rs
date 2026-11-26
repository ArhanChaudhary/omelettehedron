//! The packed table (`Table5`, the optimal solver's large tables) mapped read-only from a
//! sealed memfd equals the table read into memory, and unsealed, corrupt or truncated files
//! are refused (see `Table5::load_mapped`).
#![cfg(all(target_os = "linux", not(miri)))]

use omelettehedron::prune::{checksum, Table5};
use std::fs::File;
use std::io::Write;
use std::os::fd::{AsRawFd, FromRawFd};
use std::os::unix::fs::FileExt;
use std::path::PathBuf;

const SEALED: i32 = libc::F_SEAL_WRITE | libc::F_SEAL_GROW | libc::F_SEAL_SHRINK;

fn memfd(bytes: &[u8], seals: i32) -> File {
    // SAFETY: valid static name, documented flags; a successful fd is owned exactly once.
    let fd = unsafe { libc::memfd_create(c"fto-table5-test".as_ptr(), libc::MFD_ALLOW_SEALING) };
    assert!(fd >= 0, "{}", std::io::Error::last_os_error());
    // SAFETY: this is our freshly returned, open descriptor.
    let mut file = unsafe { File::from_raw_fd(fd) };
    file.write_all(bytes).unwrap();
    if seals != 0 {
        // SAFETY: adding the supplied seals to the valid owned descriptor.
        assert_eq!(unsafe { libc::fcntl(fd, libc::F_ADD_SEALS, seals) }, 0);
    }
    file
}

fn path(file: &File) -> PathBuf {
    PathBuf::from(format!("/proc/self/fd/{}", file.as_raw_fd()))
}

/// A table file of `n` entries (40 per word) with histogram [4, 5, 6].
fn fixture() -> (Vec<u8>, Vec<u64>, usize) {
    let n: usize = 100_003;
    let data: Vec<u64> = (0..n.div_ceil(40))
        .map(|i| (i as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15))
        .collect();
    let mut bytes = b"FTOPRUN5".to_vec();
    for value in [n as u64, 3, 4, 5, 6, checksum(&data)] {
        bytes.extend(value.to_le_bytes());
    }
    for value in &data {
        bytes.extend(value.to_le_bytes());
    }
    (bytes, data, n)
}

#[test]
fn mapped_table5_matches_the_read_table() {
    let (bytes, data, n) = fixture();
    let file = memfd(&bytes, SEALED);
    let heap = Table5::load(&path(&file), n).unwrap();
    let mapped = Table5::load_mapped(&path(&file), n).unwrap();
    assert_eq!(&*mapped.data, &*heap.data);
    assert_eq!(&*mapped.data, &data);
    assert_eq!(mapped.hist, vec![4, 5, 6]);
    assert!(file.write_at(&[0], 56).is_err());
    drop(file);
    for i in 0..n {
        assert_eq!(mapped.get(i), heap.get(i));
    }
}

#[test]
fn mapped_table5_refuses_bad_files() {
    let (bytes, _, n) = fixture();
    for seals in [0, libc::F_SEAL_SHRINK | libc::F_SEAL_GROW] {
        let f = memfd(&bytes, seals);
        assert!(Table5::load_mapped(&path(&f), n).is_err());
    }
    let mut corrupt = bytes.clone();
    *corrupt.last_mut().unwrap() ^= 1;
    let f = memfd(&corrupt, SEALED);
    assert!(Table5::load_mapped(&path(&f), n).is_err());
    let f = memfd(&bytes[..bytes.len() - 1], SEALED);
    assert!(Table5::load_mapped(&path(&f), n).is_err());
    let f = memfd(&bytes, SEALED);
    assert!(Table5::load_mapped(&path(&f), n + 1).is_err());
}
