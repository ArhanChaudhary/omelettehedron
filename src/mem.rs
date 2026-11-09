//! Large, 2 MiB-aligned buffers backed by transparent huge pages.
//!
//! Random lookups into tables of hundreds of megabytes are dominated by cache and TLB misses;
//! huge pages remove most of the TLB misses. `madvise` only works on page-aligned ranges, so
//! the allocation itself must be 2 MiB aligned (a plain `Vec` is not).

use rayon::prelude::*;
use std::alloc::{alloc, dealloc, Layout};
use std::ops::{Deref, DerefMut};
use std::ptr::NonNull;
use std::sync::atomic::{AtomicBool, Ordering};

const HUGE: usize = 2 << 20;

/// Fixed-size buffer of plain-old-data elements, zero-initialised, huge-page aligned.
pub struct HugeBuf<T: Copy> {
    ptr: NonNull<T>,
    len: usize,
    layout: Option<Layout>,
}

// SAFETY: HugeBuf owns its memory exclusively, like Vec<T>.
unsafe impl<T: Copy + Send> Send for HugeBuf<T> {}
// SAFETY: shared access only hands out &[T].
unsafe impl<T: Copy + Sync> Sync for HugeBuf<T> {}

static WARNED: AtomicBool = AtomicBool::new(false);

impl<T: Copy> HugeBuf<T> {
    /// Zero-filled buffer of `len` elements. `T` must be valid when all bits are zero.
    pub fn zeroed(len: usize) -> HugeBuf<T> {
        let bytes = (len * std::mem::size_of::<T>())
            .max(1)
            .next_multiple_of(HUGE);
        let align = HUGE.max(std::mem::align_of::<T>());
        let layout = Layout::from_size_align(bytes, align).expect("layout");
        // Allocate without touching the memory (an aligned zeroed allocation would memset it,
        // faulting in small pages before madvise), then advise, then zero in parallel.
        // SAFETY: layout has non-zero size
        let p = unsafe { alloc(layout) };
        let Some(ptr) = NonNull::new(p as *mut T) else {
            std::alloc::handle_alloc_error(layout)
        };
        #[cfg(target_os = "linux")]
        {
            // SAFETY: advisory call on memory we own; the range is 2 MiB aligned
            let r = unsafe { libc::madvise(p as *mut libc::c_void, bytes, libc::MADV_HUGEPAGE) };
            if r != 0 && !WARNED.swap(true, Ordering::Relaxed) {
                eprintln!("warning: madvise(MADV_HUGEPAGE) failed; tables use normal pages");
            }
        }
        // SAFETY: the allocation is `bytes` long; zero bytes are valid for T (documented)
        let raw = unsafe { std::slice::from_raw_parts_mut(p, bytes) };
        raw.par_chunks_mut(HUGE).for_each(|c| c.fill(0));
        HugeBuf {
            ptr,
            len,
            layout: Some(layout),
        }
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
        if let Some(l) = self.layout {
            // SAFETY: allocated with exactly this layout in `zeroed`
            unsafe { dealloc(self.ptr.as_ptr() as *mut u8, l) };
        }
    }
}
