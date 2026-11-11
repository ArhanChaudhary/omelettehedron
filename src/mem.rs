//! Large, 2 MiB-aligned buffers backed by transparent huge pages.
//!
//! Random lookups into tables of hundreds of megabytes are dominated by cache and TLB misses;
//! huge pages remove most of the TLB misses. `madvise` only works on page-aligned ranges, so
//! the allocation itself must be 2 MiB aligned (a plain `Vec` is not). (+18% / +16% nodes/s,
//! first solution 0.77x at 48T when added, opt/hugepages.)

use rayon::prelude::*;
use std::alloc::{alloc, dealloc, Layout};
use std::ops::{Deref, DerefMut};
use std::ptr::NonNull;

const HUGE: usize = 2 << 20;

/// Fixed-size buffer of plain-old-data elements, zero-initialised, huge-page aligned.
pub struct HugeBuf<T: Copy> {
    ptr: NonNull<T>,
    len: usize,
    layout: Option<Layout>,
    /// bytes of a 1 GiB-page mapping (0: allocated with `layout`)
    #[cfg_attr(not(all(target_os = "linux", not(miri))), allow(dead_code))]
    mapped: usize,
}

// SAFETY: HugeBuf owns its memory exclusively, like Vec<T>.
unsafe impl<T: Copy + Send> Send for HugeBuf<T> {}
// SAFETY: shared access only hands out &[T].
unsafe impl<T: Copy + Sync> Sync for HugeBuf<T> {}

/// Ask the kernel to back the `bytes` bytes at `p` with transparent huge pages.
#[cfg(all(target_os = "linux", not(miri)))]
fn advise_huge(p: *mut u8, bytes: usize) {
    use std::sync::atomic::{AtomicBool, Ordering};
    static WARNED: AtomicBool = AtomicBool::new(false);
    // SAFETY: advisory call on memory the caller owns; the range is 2 MiB aligned
    let r = unsafe { libc::madvise(p.cast::<libc::c_void>(), bytes, libc::MADV_HUGEPAGE) };
    if r != 0 && !WARNED.swap(true, Ordering::Relaxed) {
        eprintln!("warning: madvise(MADV_HUGEPAGE) failed; tables use normal pages");
    }
}

/// Elsewhere (and under Miri) the buffer simply uses the system's normal pages.
#[cfg(not(all(target_os = "linux", not(miri))))]
fn advise_huge(_p: *mut u8, _bytes: usize) {}

/// Map `bytes` (rounded up to whole 1 GiB pages) from the hugetlb pool, if `FTO_HUGE1G=1`, the
/// buffer is large and the pool has the pages (reserved with
/// `echo N > /sys/kernel/mm/hugepages/hugepages-1048576kB/nr_hugepages`). The kernel hands the
/// memory out zeroed. (A-centre/edge table on a 1 GiB page: optimal proofs 0.96x / 0.97x time,
/// nodes/s 1.01x / 1.03x at 48T / 1T, two interleaved rounds, when added.)
#[cfg(all(target_os = "linux", not(miri)))]
fn map_1g(bytes: usize) -> Option<(NonNull<u8>, usize)> {
    const GIB: usize = 1 << 30;
    let on = std::env::var("FTO_HUGE1G").is_ok_and(|v| v == "1");
    if !on || bytes < GIB / 4 {
        return None;
    }
    let len = bytes.checked_next_multiple_of(GIB)?;
    let flags = libc::MAP_PRIVATE | libc::MAP_ANONYMOUS | libc::MAP_HUGETLB | libc::MAP_HUGE_1GB;
    // SAFETY: an anonymous mapping of fresh memory; the result is checked below
    let p = unsafe {
        libc::mmap(
            std::ptr::null_mut(),
            len,
            libc::PROT_READ | libc::PROT_WRITE,
            flags,
            -1,
            0,
        )
    };
    let ok = p != libc::MAP_FAILED;
    if std::env::var_os("FTO_DEBUG").is_some() {
        eprintln!(
            "1 GiB pages for {len} bytes: {}",
            if ok { "yes" } else { "no (THP)" }
        );
    }
    if !ok {
        return None;
    }
    NonNull::new(p.cast::<u8>()).map(|p| (p, len))
}

#[cfg(not(all(target_os = "linux", not(miri))))]
fn map_1g(_bytes: usize) -> Option<(NonNull<u8>, usize)> {
    None
}

/// Layout of a buffer of `len` elements of `T`, rounded up to whole huge pages and aligned to
/// them; `None` if the size is beyond the address space.
fn try_huge_layout<T>(len: usize) -> Option<Layout> {
    len.checked_mul(std::mem::size_of::<T>())
        .and_then(|b| b.max(1).checked_next_multiple_of(HUGE))
        .and_then(|b| Layout::from_size_align(b, HUGE.max(std::mem::align_of::<T>())).ok())
}

/// [`try_huge_layout`]; a size beyond the address space is handled like a failed allocation
/// (the process aborts, as for any allocation the system cannot satisfy).
fn huge_layout<T>(len: usize) -> Layout {
    let Some(layout) = try_huge_layout::<T>(len) else {
        eprintln!(
            "error: cannot allocate {len} elements of {} bytes",
            std::mem::size_of::<T>()
        );
        std::process::abort()
    };
    layout
}

impl<T: Copy> HugeBuf<T> {
    /// Zero-filled buffer of `len` elements. `T` must be valid when all bits are zero.
    pub fn zeroed(len: usize) -> HugeBuf<T> {
        let layout = huge_layout::<T>(len);
        let bytes = layout.size();
        if let Some((p, mapped)) = map_1g(bytes) {
            // (1 GiB pages are zeroed by the kernel)
            return HugeBuf {
                ptr: p.cast::<T>(),
                len,
                layout: None,
                mapped,
            };
        }
        // Allocate without touching the memory (an aligned zeroed allocation would memset it,
        // faulting in small pages before madvise), then advise, then zero in parallel.
        // SAFETY: layout has non-zero size
        let p = unsafe { alloc(layout) };
        let Some(ptr) = NonNull::new(p.cast::<T>()) else {
            std::alloc::handle_alloc_error(layout)
        };
        advise_huge(p, bytes);
        // SAFETY: the allocation is `bytes` long; zero bytes are valid for T (documented)
        let raw = unsafe { std::slice::from_raw_parts_mut(p, bytes) };
        // (sequential under Miri: rayon's crossbeam-epoch uses a container-of cast that Stacked
        // Borrows rejects, which would hide checks of this module's own unsafe code)
        #[cfg(not(miri))]
        raw.par_chunks_mut(HUGE).for_each(|c| c.fill(0));
        #[cfg(miri)]
        raw.fill(0);
        HugeBuf {
            ptr,
            len,
            layout: Some(layout),
            mapped: 0,
        }
    }

    /// Buffer of `len` elements whose contents are unspecified (any bit pattern must be valid
    /// for `T`); the caller overwrites all of it. Skips the zeroing pass of [`Self::zeroed`].
    pub fn for_overwrite(len: usize) -> HugeBuf<T> {
        let layout = huge_layout::<T>(len);
        Self::alloc_for_overwrite(len, layout)
            .unwrap_or_else(|| std::alloc::handle_alloc_error(layout))
    }

    /// [`Self::for_overwrite`], or `None` if the memory cannot be allocated (where the system
    /// reports that, e.g. WebAssembly's fixed or full memory) instead of aborting.
    pub fn try_for_overwrite(len: usize) -> Option<HugeBuf<T>> {
        Self::alloc_for_overwrite(len, try_huge_layout::<T>(len)?)
    }

    /// A buffer of `len` elements with `layout` (their huge-page layout), contents unspecified.
    fn alloc_for_overwrite(len: usize, layout: Layout) -> Option<HugeBuf<T>> {
        let bytes = layout.size();
        if let Some((p, mapped)) = map_1g(bytes) {
            return Some(HugeBuf {
                ptr: p.cast::<T>(),
                len,
                layout: None,
                mapped,
            });
        }
        // SAFETY: layout has non-zero size
        let p = unsafe { alloc(layout) };
        let ptr = NonNull::new(p.cast::<T>())?;
        advise_huge(p, bytes);
        Some(HugeBuf {
            ptr,
            len,
            layout: Some(layout),
            mapped: 0,
        })
    }

    /// Remove the first `n` elements (all of them if there are fewer), moving the rest to the
    /// front; the memory stays allocated. (The WebAssembly module reads a table file into one
    /// buffer and removes its header this way, without a second copy of the table.)
    pub fn drop_front(&mut self, n: usize) {
        let n = n.min(self.len);
        self.copy_within(n.., 0);
        self.len -= n;
    }

    /// Copy a slice into a new huge-page buffer.
    pub fn from_slice(s: &[T]) -> HugeBuf<T> {
        let mut b = Self::zeroed(s.len());
        b.copy_from_slice(s);
        b
    }
}

impl<T: Copy> Deref for HugeBuf<T> {
    type Target = [T];
    #[inline(always)]
    fn deref(&self) -> &[T] {
        // SAFETY: ptr is valid for len elements for the lifetime of self
        unsafe { std::slice::from_raw_parts(self.ptr.as_ptr(), self.len) }
    }
}

impl<T: Copy> DerefMut for HugeBuf<T> {
    #[inline(always)]
    fn deref_mut(&mut self) -> &mut [T] {
        // SAFETY: unique access through &mut self
        unsafe { std::slice::from_raw_parts_mut(self.ptr.as_ptr(), self.len) }
    }
}

impl<T: Copy> Drop for HugeBuf<T> {
    fn drop(&mut self) {
        #[cfg(all(target_os = "linux", not(miri)))]
        if self.mapped > 0 {
            // SAFETY: mapped with exactly this length in `map_1g`
            unsafe { libc::munmap(self.ptr.as_ptr().cast::<libc::c_void>(), self.mapped) };
            return;
        }
        if let Some(l) = self.layout {
            // SAFETY: allocated with exactly this layout in `zeroed`
            unsafe { dealloc(self.ptr.as_ptr().cast::<u8>(), l) };
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn huge_buf_contents() {
        let z: HugeBuf<u64> = HugeBuf::zeroed(1000);
        assert_eq!(z.len(), 1000);
        assert!(z.iter().all(|&x| x == 0));
        let v: Vec<u32> = (0..5000).collect();
        let b = HugeBuf::from_slice(&v);
        assert_eq!(&b[..], &v[..]);
        let mut w: HugeBuf<u16> = HugeBuf::for_overwrite(300);
        w.fill(7);
        assert!(w.iter().all(|&x| x == 7));
        let e: HugeBuf<u8> = HugeBuf::zeroed(0);
        assert!(e.is_empty());
    }

    #[test]
    fn huge_buf_try_and_drop_front() {
        let mut b: HugeBuf<u64> = HugeBuf::try_for_overwrite(1000).unwrap();
        for (i, x) in b.iter_mut().enumerate() {
            *x = i as u64;
        }
        b.drop_front(3);
        assert_eq!(b.len(), 997);
        assert!(b.iter().enumerate().all(|(i, &x)| x == i as u64 + 3));
        b.drop_front(5000);
        assert!(b.is_empty());
        // (a size beyond the address space is refused, not aborted on)
        assert!(HugeBuf::<u64>::try_for_overwrite(usize::MAX / 4).is_none());
    }
}
