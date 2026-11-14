//! Block files: the blocks of one class (one bitmap of `words` words per B-centre state, in
//! order), compressed, in chunks of [`CHUNK`] blocks that decode independently.
//!
//! A chunk is a byte stream of codes:
//!
//! * `0` + varint `n`: `n` empty blocks;
//! * `1` + varint `n`: `n` full blocks (every corner state);
//! * `2` + `8 * words` bytes: a block as it is (little-endian words);
//! * `3` + varint `k` + `k` little-endian `u16`: a block with `k` bits set, their indices in
//!   increasing order (used while that is shorter than the plain block).
//!
//! After the chunks: an index (per chunk its offset, length and hash) and a footer of nine
//! `u64` (magic, `n_x`, `words`, chunk size, chunks, index offset, bits set, nonempty blocks,
//! the hash of the index and the footer). A reader checks the footer, that hash and each
//! chunk's hash.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

/// Blocks per chunk.
pub const CHUNK: usize = 4096;
/// The footer's magic number ("P2BFBLK1").
const MAGIC: u64 = u64::from_le_bytes(*b"P2BFBLK1");
/// Footer words.
const FOOTER: usize = 9;

const T_EMPTY: u8 = 0;
const T_FULL: u8 = 1;
const T_RAW: u8 = 2;
const T_SPARSE: u8 = 3;

/// A 64-bit hash of a byte string (multiply-rotate over 8-byte words), to detect corruption.
pub fn hash(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0x9e37_79b9_7f4a_7c15 ^ bytes.len() as u64;
    let (words, rest) = bytes.as_chunks::<8>();
    for c in words {
        h = (h ^ u64::from_le_bytes(*c))
            .wrapping_mul(0xbf58_476d_1ce4_e5b9)
            .rotate_left(31);
    }
    for (i, &b) in rest.iter().enumerate() {
        h = (h ^ (u64::from(b) << (8 * i))).wrapping_mul(0x94d0_49bb_1331_11eb);
    }
    h ^ h >> 29
}

/// The hash of the index and the footer's first words (all but the hash itself).
fn tail_hash(idx: &[u8], footer: &[u64]) -> u64 {
    let mut b = idx.to_vec();
    for w in footer.iter().take(FOOTER - 1) {
        b.extend_from_slice(&w.to_le_bytes());
    }
    hash(&b)
}

fn put_varint(out: &mut Vec<u8>, mut v: u64) {
    while v >= 0x80 {
        out.push((v as u8) | 0x80);
        v >>= 7;
    }
    out.push(v as u8);
}

fn get_varint(b: &[u8], i: &mut usize) -> Option<u64> {
    let mut v = 0u64;
    let mut shift = 0;
    loop {
        let byte = *b.get(*i)?;
        *i += 1;
        v |= u64::from(byte & 0x7f) << shift;
        if byte < 0x80 {
            return Some(v);
        }
        shift += 7;
        if shift > 63 {
            return None;
        }
    }
}

/// The bad-data error of a file.
fn bad(path: &Path, what: &str) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        format!("{}: {what}", path.display()),
    )
}

/// Bits set and nonempty blocks of a file.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Summary {
    /// bits set
    pub bits: u64,
    /// blocks with a bit set
    pub nonempty: u64,
    /// bytes of the file
    pub bytes: u64,
}

/// Writes the blocks of one class in order.
pub struct ClassWriter {
    out: BufWriter<File>,
    path: PathBuf,
    n_x: usize,
    words: usize,
    n_p: usize,
    /// the chunk being encoded
    buf: Vec<u8>,
    /// the run being collected: its kind (empty or full) and length
    run: Option<(u8, u64)>,
    index: Vec<(u64, u64, u64)>,
    pos: u64,
    x: usize,
    sum: Summary,
    scratch: Vec<u16>,
}

impl ClassWriter {
    /// Create the file (truncating it).
    ///
    /// # Errors
    ///
    /// If the file cannot be created.
    pub fn create(path: &Path, n_x: usize, n_p: usize) -> io::Result<ClassWriter> {
        let f = File::create(path)?;
        Ok(ClassWriter {
            out: BufWriter::with_capacity(8 << 20, f),
            path: path.to_path_buf(),
            n_x,
            words: n_p.div_ceil(64),
            n_p,
            buf: Vec::new(),
            run: None,
            index: Vec::new(),
            pos: 0,
            x: 0,
            sum: Summary::default(),
            scratch: Vec::new(),
        })
    }

    fn flush_run(&mut self) {
        if let Some((t, n)) = self.run.take() {
            self.buf.push(t);
            put_varint(&mut self.buf, n);
        }
    }

    fn add_run(&mut self, t: u8, n: u64) {
        match &mut self.run {
            Some((t0, n0)) if *t0 == t => *n0 += n,
            _ => {
                self.flush_run();
                self.run = Some((t, n));
            }
        }
    }

    /// After a block: close the chunk if it is complete.
    fn advance(&mut self, n: usize) -> io::Result<()> {
        self.x += n;
        if self.x.is_multiple_of(CHUNK) || self.x == self.n_x {
            self.flush_run();
            let h = hash(&self.buf);
            self.out.write_all(&self.buf)?;
            self.index.push((self.pos, self.buf.len() as u64, h));
            self.pos += self.buf.len() as u64;
            self.buf.clear();
        }
        Ok(())
    }

    /// The next block. `block` holds `words` words (bits past the last state clear).
    ///
    /// # Errors
    ///
    /// If writing fails, or if all blocks are written already.
    pub fn push(&mut self, block: &[u64]) -> io::Result<()> {
        if self.x >= self.n_x {
            return Err(bad(&self.path, "too many blocks"));
        }
        let bits: u32 = block.iter().map(|w| w.count_ones()).sum();
        if bits == 0 {
            self.add_run(T_EMPTY, 1);
        } else if bits as usize == self.n_p {
            self.sum.bits += u64::from(bits);
            self.sum.nonempty += 1;
            self.add_run(T_FULL, 1);
        } else {
            self.sum.bits += u64::from(bits);
            self.sum.nonempty += 1;
            self.flush_run();
            if (bits as usize) < 4 * self.words - 2 {
                self.buf.push(T_SPARSE);
                put_varint(&mut self.buf, u64::from(bits));
                self.scratch.clear();
                for (w, &word) in block.iter().enumerate() {
                    let mut v = word;
                    while v != 0 {
                        self.scratch
                            .push((w * 64 + v.trailing_zeros() as usize) as u16);
                        v &= v - 1;
                    }
                }
                for &i in &self.scratch {
                    self.buf.extend_from_slice(&i.to_le_bytes());
                }
            } else {
                self.buf.push(T_RAW);
                for &w in block {
                    self.buf.extend_from_slice(&w.to_le_bytes());
                }
            }
        }
        self.advance(1)
    }

    /// `n` empty blocks.
    ///
    /// # Errors
    ///
    /// If writing fails, or if that is more blocks than are left.
    pub fn push_empty(&mut self, mut n: usize) -> io::Result<()> {
        if self.x + n > self.n_x {
            return Err(bad(&self.path, "too many blocks"));
        }
        while n > 0 {
            let room = CHUNK - self.x % CHUNK;
            let k = n.min(room).min(self.n_x - self.x);
            self.add_run(T_EMPTY, k as u64);
            self.advance(k)?;
            n -= k;
        }
        Ok(())
    }

    /// Write the index and the footer, and sync the file to disk.
    ///
    /// # Errors
    ///
    /// If writing fails, or if not every block was written.
    pub fn finish(mut self) -> io::Result<Summary> {
        if self.x != self.n_x {
            return Err(bad(&self.path, "blocks missing"));
        }
        let mut idx = Vec::with_capacity(self.index.len() * 24);
        for &(o, l, h) in &self.index {
            idx.extend_from_slice(&o.to_le_bytes());
            idx.extend_from_slice(&l.to_le_bytes());
            idx.extend_from_slice(&h.to_le_bytes());
        }
        let mut footer: [u64; FOOTER] = [
            MAGIC,
            self.n_x as u64,
            self.words as u64,
            CHUNK as u64,
            self.index.len() as u64,
            self.pos,
            self.sum.bits,
            self.sum.nonempty,
            0,
        ];
        footer[FOOTER - 1] = tail_hash(&idx, &footer);
        self.out.write_all(&idx)?;
        for w in footer {
            self.out.write_all(&w.to_le_bytes())?;
        }
        self.out.flush()?;
        let f = self
            .out
            .into_inner()
            .map_err(io::IntoInnerError::into_error)?;
        f.sync_all()?;
        self.sum.bytes = self.pos + idx.len() as u64 + (FOOTER * 8) as u64;
        Ok(self.sum)
    }
}

/// A block as stored.
#[derive(Clone, Copy, Debug)]
pub enum Blk<'a> {
    /// no bit set
    Empty,
    /// every bit set
    Full,
    /// the words, little-endian
    Raw(&'a [u8]),
    /// the indices of the bits set, little-endian `u16`
    Sparse(&'a [u8]),
}

impl Blk<'_> {
    /// Write the block into `out` (`words` words; bits past `n_p` stay clear).
    pub fn expand(&self, out: &mut [u64], n_p: usize) {
        match *self {
            Blk::Empty => out.fill(0),
            Blk::Full => {
                out.fill(u64::MAX);
                if !n_p.is_multiple_of(64) {
                    if let Some(last) = out.last_mut() {
                        *last = (1u64 << (n_p % 64)) - 1;
                    }
                }
            }
            Blk::Raw(b) => {
                for (w, c) in out.iter_mut().zip(b.as_chunks::<8>().0) {
                    *w = u64::from_le_bytes(*c);
                }
            }
            Blk::Sparse(b) => {
                out.fill(0);
                for c in b.as_chunks::<2>().0 {
                    let i = usize::from(u16::from_le_bytes(*c));
                    if let Some(w) = out.get_mut(i / 64) {
                        *w |= 1 << (i % 64);
                    }
                }
            }
        }
    }

    /// The indices of the bits set, appended to `out` (not for [`Blk::Full`]: all of them).
    pub fn bits(&self, out: &mut Vec<u16>, n_p: usize) {
        match *self {
            Blk::Empty => {}
            Blk::Full => out.extend(0..n_p as u16),
            Blk::Raw(b) => {
                for (w, c) in b.as_chunks::<8>().0.iter().enumerate() {
                    let mut v = u64::from_le_bytes(*c);
                    while v != 0 {
                        out.push((w * 64 + v.trailing_zeros() as usize) as u16);
                        v &= v - 1;
                    }
                }
            }
            Blk::Sparse(b) => {
                out.extend(b.as_chunks::<2>().0.iter().map(|c| u16::from_le_bytes(*c)));
            }
        }
    }
}

/// Decodes one chunk: calls `f(x, block)` for each block of the chunk that starts at block
/// `x0`, runs of empty blocks skipped unless `with_empty`. `None` if the chunk is malformed.
pub fn decode_chunk(
    bytes: &[u8],
    x0: usize,
    n_blocks: usize,
    words: usize,
    with_empty: bool,
    mut f: impl FnMut(usize, Blk),
) -> Option<()> {
    let mut i = 0;
    let mut x = x0;
    let end = x0 + n_blocks;
    while i < bytes.len() {
        let t = bytes[i];
        i += 1;
        match t {
            T_EMPTY | T_FULL => {
                let n = get_varint(bytes, &mut i)? as usize;
                if x + n > end {
                    return None;
                }
                if t == T_FULL || with_empty {
                    let b = if t == T_FULL { Blk::Full } else { Blk::Empty };
                    for y in x..x + n {
                        f(y, b);
                    }
                }
                x += n;
            }
            T_RAW => {
                let b = bytes.get(i..i + 8 * words)?;
                i += 8 * words;
                if x >= end {
                    return None;
                }
                f(x, Blk::Raw(b));
                x += 1;
            }
            T_SPARSE => {
                let k = get_varint(bytes, &mut i)? as usize;
                let b = bytes.get(i..i + 2 * k)?;
                i += 2 * k;
                if x >= end {
                    return None;
                }
                f(x, Blk::Sparse(b));
                x += 1;
            }
            _ => return None,
        }
    }
    (x == end).then_some(())
}

/// Reads a class file.
pub struct ClassReader {
    file: File,
    path: PathBuf,
    /// per chunk: offset, length, hash
    index: Vec<(u64, u64, u64)>,
    n_x: usize,
    words: usize,
    sum: Summary,
}

impl ClassReader {
    /// Open a file and check its footer and index against the instance's sizes.
    ///
    /// # Errors
    ///
    /// If the file cannot be read or is not a block file of these sizes.
    pub fn open(path: &Path, n_x: usize, n_p: usize) -> io::Result<ClassReader> {
        let file = File::open(path)?;
        let len = file.metadata()?.len();
        let fl = (FOOTER * 8) as u64;
        if len < fl {
            return Err(bad(path, "too short"));
        }
        let mut fb = vec![0u8; FOOTER * 8];
        read_at(&file, &mut fb, len - fl)?;
        let w: Vec<u64> = fb
            .as_chunks::<8>()
            .0
            .iter()
            .map(|c| u64::from_le_bytes(*c))
            .collect();
        let words = n_p.div_ceil(64);
        let n_chunks = n_x.div_ceil(CHUNK);
        if w[0] != MAGIC
            || w[1] != n_x as u64
            || w[2] != words as u64
            || w[3] != CHUNK as u64
            || w[4] != n_chunks as u64
        {
            return Err(bad(path, "not a block file of this instance"));
        }
        let idx_off = w[5];
        if idx_off + (n_chunks * 24) as u64 + fl != len {
            return Err(bad(path, "wrong length"));
        }
        let mut ib = vec![0u8; n_chunks * 24];
        read_at(&file, &mut ib, idx_off)?;
        if tail_hash(&ib, &w) != w[8] {
            return Err(bad(path, "index checksum"));
        }
        let index = ib
            .as_chunks::<24>()
            .0
            .iter()
            .map(|c| {
                let g = |k: usize| {
                    let mut a = [0u8; 8];
                    a.copy_from_slice(&c[8 * k..8 * k + 8]);
                    u64::from_le_bytes(a)
                };
                (g(0), g(1), g(2))
            })
            .collect();
        Ok(ClassReader {
            file,
            path: path.to_path_buf(),
            index,
            n_x,
            words,
            sum: Summary {
                bits: w[6],
                nonempty: w[7],
                bytes: len,
            },
        })
    }

    /// Bits set, nonempty blocks and bytes.
    pub fn summary(&self) -> Summary {
        self.sum
    }

    /// Chunks of the file.
    pub fn chunks(&self) -> usize {
        self.index.len()
    }

    /// The first block and the number of blocks of chunk `i`.
    pub fn chunk_range(&self, i: usize) -> (usize, usize) {
        let x0 = i * CHUNK;
        (x0, CHUNK.min(self.n_x.saturating_sub(x0)))
    }

    /// Read chunk `i` into `buf` and check its hash.
    ///
    /// # Errors
    ///
    /// If reading fails or the hash does not match.
    pub fn read_chunk(&self, i: usize, buf: &mut Vec<u8>) -> io::Result<()> {
        let &(o, l, h) = self
            .index
            .get(i)
            .ok_or_else(|| bad(&self.path, "no such chunk"))?;
        buf.resize(l as usize, 0);
        read_at(&self.file, buf, o)?;
        if hash(buf) != h {
            return Err(bad(&self.path, &format!("chunk {i} checksum")));
        }
        Ok(())
    }

    /// Decode chunk `i` (read into `buf`), calling `f` as [`decode_chunk`] does.
    ///
    /// # Errors
    ///
    /// If reading fails or the chunk is malformed.
    pub fn for_chunk(
        &self,
        i: usize,
        buf: &mut Vec<u8>,
        with_empty: bool,
        f: impl FnMut(usize, Blk),
    ) -> io::Result<()> {
        self.read_chunk(i, buf)?;
        let (x0, n) = self.chunk_range(i);
        decode_chunk(buf, x0, n, self.words, with_empty, f)
            .ok_or_else(|| bad(&self.path, &format!("chunk {i} malformed")))
    }
}

/// Read exactly `buf.len()` bytes at `offset`.
fn read_at(file: &File, buf: &mut [u8], offset: u64) -> io::Result<()> {
    #[cfg(unix)]
    {
        std::os::unix::fs::FileExt::read_exact_at(file, buf, offset)
    }
    #[cfg(not(unix))]
    {
        use std::io::{Read, Seek, SeekFrom};
        let mut f = file.try_clone()?;
        f.seek(SeekFrom::Start(offset))?;
        f.read_exact(buf)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{Rng, SeedableRng};

    /// Random blocks of every kind (empty, full, sparse, dense) and runs across chunk
    /// boundaries come back as written; the summary counts them.
    #[test]
    fn blocks_round_trip() {
        let dir = std::env::temp_dir().join(format!("p2bfs-codec-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(5);
        for &(n_x, n_p) in &[
            (10_000usize, 11_520usize),
            (9_000, 360),
            (5_000, 1),
            (100, 100),
        ] {
            let words = n_p.div_ceil(64);
            let mut blocks = vec![vec![0u64; words]; n_x];
            let mut x = 0;
            while x < n_x {
                let run = rng.gen_range(1..3000).min(n_x - x);
                let kind = rng.gen_range(0..4);
                for b in &mut blocks[x..x + run] {
                    for i in 0..n_p {
                        let on = match kind {
                            0 => false,
                            1 => true,
                            2 => rng.gen_range(0..200) == 0,
                            _ => rng.gen_bool(0.5),
                        };
                        if on {
                            b[i / 64] |= 1 << (i % 64);
                        }
                    }
                }
                x += run;
            }
            let path = dir.join("c.blk");
            let mut w = ClassWriter::create(&path, n_x, n_p).unwrap();
            let mut i = 0;
            while i < n_x {
                if blocks[i].iter().all(|&v| v == 0) && rng.gen_bool(0.5) {
                    let mut j = i;
                    while j < n_x && blocks[j].iter().all(|&v| v == 0) {
                        j += 1;
                    }
                    w.push_empty(j - i).unwrap();
                    i = j;
                } else {
                    w.push(&blocks[i]).unwrap();
                    i += 1;
                }
            }
            let sum = w.finish().unwrap();
            let r = ClassReader::open(&path, n_x, n_p).unwrap();
            assert_eq!(r.summary(), sum);
            let bits: u64 = blocks
                .iter()
                .flatten()
                .map(|w| u64::from(w.count_ones()))
                .sum();
            assert_eq!(sum.bits, bits);
            let mut got = vec![vec![u64::MAX; words]; n_x];
            let mut buf = vec![];
            for c in 0..r.chunks() {
                r.for_chunk(c, &mut buf, true, |x, b| b.expand(&mut got[x], n_p))
                    .unwrap();
            }
            assert_eq!(got, blocks);
            // the bit lists agree with the blocks
            for c in 0..r.chunks() {
                r.for_chunk(c, &mut buf, false, |x, b| {
                    let mut l = vec![];
                    b.bits(&mut l, n_p);
                    let want: Vec<u16> = (0..n_p as u16)
                        .filter(|&i| blocks[x][i as usize / 64] >> (i % 64) & 1 == 1)
                        .collect();
                    assert_eq!(l, want);
                })
                .unwrap();
            }
            // a damaged byte anywhere is caught (by the footer, the index or a chunk)
            assert!(ClassReader::open(&path, n_x + 1, n_p).is_err());
            let good = std::fs::read(&path).unwrap();
            for _ in 0..20 {
                let mut bytes = good.clone();
                let i = rng.gen_range(0..bytes.len());
                bytes[i] ^= 1 << rng.gen_range(0..8);
                std::fs::write(&path, &bytes).unwrap();
                let caught = match ClassReader::open(&path, n_x, n_p) {
                    Err(_) => true,
                    Ok(r) => (0..r.chunks()).any(|c| r.read_chunk(c, &mut buf).is_err()),
                };
                assert!(caught, "byte {i}");
            }
        }
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
