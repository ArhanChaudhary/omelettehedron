//! Checks of the oriented phase-1 tables ([`super::p1disk`]) that do not trust them: the states
//! of a level and their classes, explicit solutions read off the table, a meet-in-the-middle
//! lower bound with no table, an exact strict table made by the two-phase code, and slot by
//! slot comparisons between the tables.

use super::oriented::{coset_parts, decode_ep, ep_of, P1Space, Target, N_EPERM, N_TURN, ROW_BYTES};
use super::p1disk::{nib, row_buf, Nibbles, P1Error, Plane, NIB_STRIDE};
use crate::coord::{move_tables, mset, Coords, N_ECOSET, N_MSET};
use crate::cube::{Cube, N_MOVES};
use crate::mem::HugeBuf;
use crate::sym::syms;
use crate::tables::{ae_targets, AeIndex, AeSym};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::sync::atomic::{AtomicU8, Ordering};

type Res<T> = Result<T, P1Error>;

/// A state: A-centre colouring and edge permutation (fixed frame).
pub type State = ([u8; 12], [u8; 12]);

/// Lehmer rank of a permutation of 12.
pub fn perm_rank(p: &[u8; 12]) -> u64 {
    // the later entries smaller than p[i]: those below p[i] less the earlier ones below it
    let mut seen = 0u32;
    let mut r = 0u64;
    for i in 0..12 {
        let x = u32::from(p[i] & 15);
        let smaller = x - (seen & ((1 << x) - 1)).count_ones();
        r = r * (12 - i as u64) + u64::from(smaller);
        seen |= 1 << x;
    }
    r
}

/// A key of a state, made without the index: the colouring as 2-bit digits, the edges by their
/// Lehmer rank.
pub fn state_key(xa: &[u8; 12], ep: &[u8; 12]) -> u64 {
    let a = xa.iter().fold(0u64, |a, &x| a << 2 | u64::from(x));
    a << 29 | perm_rank(ep)
}

/// A pair `(g, t)` applied to a state.
pub fn apply_pair(g: usize, t: usize, xa: &[u8; 12], ep: &[u8; 12]) -> State {
    let s = syms();
    let a = s.apply_a(g, t, xa);
    let mut n = [0u8; 12];
    for p in 0..12 {
        n[s.e_pos[t][p] as usize] = s.e_piece[g][ep[p] as usize];
    }
    (a, n)
}

/// The smallest key of a state's images under the 288 pairs (its class).
pub fn class_key(xa: &[u8; 12], ep: &[u8; 12]) -> u64 {
    syms()
        .pairs
        .iter()
        .map(|&(g, t)| {
            let (a, e) = apply_pair(g as usize, t as usize, xa, ep);
            state_key(&a, &e)
        })
        .min()
        .unwrap_or(u64::MAX)
}

/// The slots of level `d` of a search directory.
///
/// # Errors
///
/// On a file error.
pub fn level_slots(sp: &P1Space, dir: &Path, d: usize) -> Res<Vec<(usize, u32)>> {
    let nrows = sp.rows();
    let counts = super::p1disk::level_rows(dir, d, nrows)?
        .ok_or_else(|| P1Error::Bad(format!("level {d} is not complete")))?;
    let plane = Plane::open(&dir.join(format!("L{d:02}.bits")), nrows, false)?;
    let mut out = vec![];
    let mut b = row_buf();
    for row in 0..nrows {
        if counts[row] == 0 {
            continue;
        }
        plane.read(row, &mut b)?;
        for (i, &x) in b[..ROW_BYTES].iter().enumerate() {
            if x != 0 {
                for k in 0..8 {
                    if x >> k & 1 == 1 {
                        out.push((row, (8 * i + k) as u32));
                    }
                }
            }
        }
        if out.len() as u64 > counts.iter().take(row + 1).sum::<u64>() {
            return Err(P1Error::Bad(format!("row {row}: more slots than counted")));
        }
    }
    Ok(out)
}

/// Every state (fixed frame) of the given slots: for each colouring of the slot's class, the
/// state that maps to the slot.
pub fn expand(sp: &P1Space, slots: &[(usize, u32)]) -> Vec<State> {
    let s = syms();
    let ms = mset();
    let mut members: HashMap<usize, Vec<u32>> = HashMap::new();
    let rows: HashSet<usize> = slots.iter().map(|x| x.0).collect();
    for r in 0..N_MSET as u32 {
        let row = (sp.ae.cls.sym[r as usize] & 0xffff) as usize;
        if rows.contains(&row) {
            members.entry(row).or_default().push(r);
        }
    }
    // inverse piece relabellings
    let inv: Vec<[u8; 12]> = (0..s.e_piece.len())
        .map(|g| {
            let mut v = [0u8; 12];
            for k in 0..12 {
                v[s.e_piece[g][k] as usize] = k as u8;
            }
            v
        })
        .collect();
    let mut out = vec![];
    for &(row, col) in slots {
        let n = decode_ep(col);
        for &c in members.get(&row).map_or(&[][..], |v| &v[..]) {
            let v = sp.ae.cls.sym[c as usize];
            let (g, t) = ((v >> 16 & 0xff) as usize, (v >> 24) as usize);
            // n[e_pos[t][p]] = e_piece[g][ep[p]]
            let ep: [u8; 12] = std::array::from_fn(|p| inv[g][n[s.e_pos[t][p] as usize] as usize]);
            out.push((ms.unrank(c), ep));
        }
    }
    out
}

/// A shortest path to the target read off the table: from a state at `level`, a move to a
/// state one lower, until 0. `None` if the table has no such move somewhere (it would be
/// inconsistent).
///
/// # Errors
///
/// On a file error.
pub fn descend(sp: &P1Space, tab: &Nibbles, st: &State, level: u8) -> Res<Option<Vec<usize>>> {
    let mut c = Cube {
        ep: st.1,
        xa: st.0,
        ..Cube::solved()
    };
    let (r0, c0) = sp.slot_of(&c.xa, &c.ep);
    if tab.get(r0, c0)? != level {
        return Ok(None);
    }
    let mut sol = vec![];
    for d in (0..level).rev() {
        let mut found = None;
        for m in 0..N_MOVES {
            let n = c.apply(m);
            let (r, col) = sp.slot_of(&n.xa, &n.ep);
            if tab.get(r, col)? == d {
                found = Some((m, n));
                break;
            }
        }
        let Some((m, n)) = found else { return Ok(None) };
        sol.push(m);
        c = n;
    }
    Ok(Some(sol))
}

/// The class keys of every state within `radius` moves of the targets, by distance: a
/// breadth-first search over classes with cube moves and [`class_key`] (no table).
pub fn class_ball(sp: &P1Space, target: Target, radius: usize) -> Vec<HashSet<u64>> {
    let mut seen: HashSet<u64> = HashSet::new();
    let mut front: Vec<State> = vec![];
    let mut keys = HashSet::new();
    for (xa, ep) in sp.target_states(target) {
        let k = class_key(&xa, &ep);
        if seen.insert(k) {
            front.push((xa, ep));
            keys.insert(k);
        }
    }
    let mut out = vec![keys];
    for _ in 1..=radius {
        let mut cand: Vec<(u64, State)> = front
            .par_iter()
            .flat_map_iter(|&(xa, ep)| {
                let c = Cube {
                    ep,
                    xa,
                    ..Cube::solved()
                };
                (0..N_MOVES).map(move |m| {
                    let n = c.apply(m);
                    (class_key(&n.xa, &n.ep), (n.xa, n.ep))
                })
            })
            .collect();
        cand.par_sort_unstable_by_key(|x| x.0);
        cand.dedup_by_key(|x| x.0);
        let mut keys = HashSet::new();
        front = cand
            .into_iter()
            .filter(|(k, _)| seen.insert(*k))
            .map(|(k, s)| {
                keys.insert(k);
                s
            })
            .collect();
        out.push(keys);
    }
    out
}

/// The smallest `i + j` over states `j <= depth` moves from `st` whose class is at distance
/// `i` in `ball` ([`class_ball`]); `None` if none is. Every move sequence of length at most
/// `depth` is tried (no same face twice in a row).
pub fn mitm<S: std::hash::BuildHasher + Sync>(
    st: &State,
    ball: &[HashSet<u64, S>],
    depth: usize,
) -> Option<usize> {
    let c = Cube {
        ep: st.1,
        xa: st.0,
        ..Cube::solved()
    };
    // the first move in parallel
    let probe = |x: &Cube, j: usize| -> Option<usize> {
        let k = class_key(&x.xa, &x.ep);
        ball.iter().position(|b| b.contains(&k)).map(|i| i + j)
    };
    let mut best = probe(&c, 0);
    let rest: Option<usize> = (0..N_MOVES)
        .into_par_iter()
        .filter_map(|m| {
            let mut best: Option<usize> = None;
            let mut stack = vec![(c.apply(m), m, 1usize)];
            while let Some((x, last, j)) = stack.pop() {
                if let Some(v) = probe(&x, j) {
                    best = Some(best.map_or(v, |b: usize| b.min(v)));
                }
                if j < depth {
                    for m2 in 0..N_MOVES {
                        if m2 >> 1 != last >> 1 {
                            stack.push((x.apply(m2), m2, j + 1));
                        }
                    }
                }
            }
            best
        })
        .min();
    if let Some(r) = rest {
        best = Some(best.map_or(r, |b| b.min(r)));
    }
    best
}

/// The exact distance to the strict target (`H'`: every hexagon one colour, turned any way) of
/// every entry of the two-phase chiral A-centre/edge index ([`AeSym`]), by a breadth-first search
/// with that index's own moves and equivalents. Two entries a byte; 15 = not reached.
pub fn strict_exact() -> (AeSym, Vec<AtomicU8>, Vec<u64>) {
    let ix = AeSym::new(true);
    let mt = move_tables();
    let n = ix.size();
    let tab: Vec<AtomicU8> = (0..n.div_ceil(2)).map(|_| AtomicU8::new(0xff)).collect();
    let get = |i: usize| tab[i >> 1].load(Ordering::Relaxed) >> ((i & 1) * 4) & 15;
    let set = |i: usize, d: u8| -> bool {
        let sh = (i & 1) * 4;
        let cell = &tab[i >> 1];
        let mut old = cell.load(Ordering::Relaxed);
        loop {
            if old >> sh & 15 != 15 {
                return false;
            }
            let new = (old & !(15 << sh)) | d << sh;
            match cell.compare_exchange_weak(old, new, Ordering::Relaxed, Ordering::Relaxed) {
                Ok(_) => return true,
                Err(x) => old = x,
            }
        }
    };
    let mut hist = vec![];
    let mut cnt = 0u64;
    for i in ae_targets(&ix) {
        if set(i, 0) {
            cnt += 1;
        }
    }
    hist.push(cnt);
    for d in 1..15u8 {
        let cnt: u64 = (0..n)
            .into_par_iter()
            .with_min_len(1 << 16)
            .filter(|&i| get(i) == d - 1)
            .map(|i| {
                let (a, e, ebits) = ix.decode(i);
                let x = Coords {
                    c: 0,
                    a,
                    b: 0,
                    e,
                    ebits,
                };
                let mut c = 0u64;
                let mut eq = vec![];
                for m in 0..N_MOVES {
                    let j = ix.index(&mt.apply(&x, m));
                    eq.clear();
                    eq.push(j);
                    ix.equivalents(j, &mut eq);
                    for &k in &eq {
                        if set(k, d) {
                            c += 1;
                        }
                    }
                }
                c
            })
            .sum();
        if cnt == 0 {
            break;
        }
        hist.push(cnt);
    }
    (ix, tab, hist)
}

/// The value of entry `i` of a table from [`strict_exact`].
#[inline]
pub fn strict_get(tab: &[AtomicU8], i: usize) -> u8 {
    tab[i >> 1].load(Ordering::Relaxed) >> ((i & 1) * 4) & 15
}

/// Slot-by-slot comparisons ([`cross`]).
#[derive(Clone, Debug, Default)]
pub struct Cross {
    /// slots by `O - S` (0..=4; index 5: anything else)
    pub o_minus_s: [u64; 6],
    /// slots by `O - T1` (0, 1; index 2: anything else)
    pub o_minus_t1: [u64; 3],
    /// cosets whose smallest `O` over the 81 turns differs from `S`
    pub min_o_ne_s: u64,
    /// cosets whose smallest `T1` differs from `S`
    pub min_t1_ne_s: u64,
    /// slots where `T1 < S`
    pub t1_below_s: u64,
    /// slots compared with the giant tier's table, and mismatches mod 3
    pub giant: (u64, u64),
}

/// Compares the oriented table `O` and the `T1` table slot by slot, and both with the exact
/// strict table `S` ([`strict_exact`], same rows: entry `row * N_ECOSET + coset`), and `O` mod 3
/// with the giant tier's packed table if given (entry `row * N_EPERM + eperm rank`).
///
/// # Errors
///
/// On a file error.
pub fn cross(
    sp: &P1Space,
    o_dir: &Path,
    t1_dir: &Path,
    strict: &[AtomicU8],
    giant: Option<&Path>,
) -> Res<Cross> {
    let nrows = sp.rows();
    let o = Nibbles::open(o_dir, nrows, true)?;
    let t1 = Nibbles::open(t1_dir, nrows, true)?;
    let giant = match giant {
        Some(p) => Some(Giant::open(p, sp.slots())?),
        None => None,
    };
    let mut out = Cross::default();
    let mut bo = HugeBuf::<u8>::for_overwrite(NIB_STRIDE);
    let mut bt = HugeBuf::<u8>::for_overwrite(NIB_STRIDE);
    let mut gb: Vec<u8> = vec![];
    for row in 0..nrows {
        o.read(row, &mut bo)?;
        t1.read(row, &mut bt)?;
        if let Some(g) = &giant {
            g.read_row(row, &mut gb)?;
        }
        let part = (0..N_ECOSET)
            .into_par_iter()
            .with_min_len(4096)
            .fold(Cross::default, |mut a, k| {
                let s = strict_get(strict, row * N_ECOSET + k);
                let (mut mo, mut mt) = (15u8, 15u8);
                let (l, bits) = coset_parts(k as u32);
                for r in 0..N_TURN {
                    let col = k * N_TURN + r;
                    let (vo, vt) = (nib(&bo, col), nib(&bt, col));
                    mo = mo.min(vo);
                    mt = mt.min(vt);
                    let ds = usize::from(vo.wrapping_sub(s));
                    a.o_minus_s[ds.min(5)] += 1;
                    let dt = usize::from(vo.wrapping_sub(vt));
                    a.o_minus_t1[dt.min(2)] += 1;
                    if vt < s {
                        a.t1_below_s += 1;
                    }
                    if !gb.is_empty() {
                        let ep = ep_of(&l, bits, r as u8);
                        let mut e16 = [0u8; 16];
                        e16[..12].copy_from_slice(&ep);
                        let i = eperm_rank16(&e16) as usize;
                        let gv = trit(&gb, i);
                        a.giant.0 += 1;
                        if gv != vo % 3 {
                            a.giant.1 += 1;
                        }
                    }
                }
                if mo != s {
                    a.min_o_ne_s += 1;
                }
                if mt != s {
                    a.min_t1_ne_s += 1;
                }
                a
            })
            .reduce(Cross::default, |mut a, b| {
                for i in 0..6 {
                    a.o_minus_s[i] += b.o_minus_s[i];
                }
                for i in 0..3 {
                    a.o_minus_t1[i] += b.o_minus_t1[i];
                }
                a.min_o_ne_s += b.min_o_ne_s;
                a.min_t1_ne_s += b.min_t1_ne_s;
                a.t1_below_s += b.t1_below_s;
                a.giant.0 += b.giant.0;
                a.giant.1 += b.giant.1;
                a
            });
        for i in 0..6 {
            out.o_minus_s[i] += part.o_minus_s[i];
        }
        for i in 0..3 {
            out.o_minus_t1[i] += part.o_minus_t1[i];
        }
        out.min_o_ne_s += part.min_o_ne_s;
        out.min_t1_ne_s += part.min_t1_ne_s;
        out.t1_below_s += part.t1_below_s;
        out.giant.0 += part.giant.0;
        out.giant.1 += part.giant.1;
        if (row + 1) % 100 == 0 {
            println!("cross: {}/{nrows} rows; {out:?}", row + 1);
        }
    }
    Ok(out)
}

/// Rank of an even permutation of 12 as the optimal solver's giant table numbers it: the
/// mixed-radix Lehmer code of its first 10 entries.
pub fn eperm_rank16(ep: &[u8; 16]) -> u32 {
    let mut seen = 0u32;
    let mut r = 0u32;
    for i in 0..10 {
        let p = u32::from(ep[i] & 15);
        let c = p - (seen & ((1 << p) - 1)).count_ones();
        r = r * (12 - i as u32) + c;
        seen |= 1 << p;
    }
    r
}

/// Digit `i % 5` of byte `i / 5` (base 3).
#[inline]
fn trit(b: &[u8], i: usize) -> u8 {
    let mut x = b[i / 5];
    for _ in 0..i % 5 {
        x /= 3;
    }
    x % 3
}

/// The optimal solver's giant-tier table (`opt-eafull-sym288-p5.prun`: distance mod 3 to the
/// solved configuration in any frame, packed five entries to a byte), read a row at a time.
pub struct Giant {
    file: std::fs::File,
    offset: u64,
    /// entries per distance, from the header
    pub hist: Vec<u64>,
}

impl Giant {
    /// Opens the file and checks its header.
    ///
    /// # Errors
    ///
    /// If it cannot be read or is not a table of `n` entries.
    pub fn open(path: &Path, n: u64) -> Res<Giant> {
        use std::os::unix::fs::FileExt;
        let err = |e: std::io::Error| P1Error::Io(format!("read {}", path.display()), e);
        let file = std::fs::File::open(path).map_err(err)?;
        let mut h = [0u8; 24];
        file.read_exact_at(&mut h, 0).map_err(err)?;
        let word = |b: &[u8]| u64::from_le_bytes(b.try_into().unwrap_or([0; 8]));
        if word(&h[8..16]) != n {
            return Err(P1Error::Bad(format!(
                "{}: {} entries, expected {n}",
                path.display(),
                word(&h[8..16])
            )));
        }
        let nh = word(&h[16..24]) as usize;
        let mut hb = vec![0u8; 8 * nh];
        file.read_exact_at(&mut hb, 24).map_err(err)?;
        let hist = hb
            .as_chunks::<8>()
            .0
            .iter()
            .map(|c| u64::from_le_bytes(*c))
            .collect();
        Ok(Giant {
            file,
            offset: (32 + 8 * nh) as u64,
            hist,
        })
    }

    /// The bytes of a row (`N_EPERM / 5`).
    ///
    /// # Errors
    ///
    /// On a read error.
    pub fn read_row(&self, row: usize, buf: &mut Vec<u8>) -> Res<()> {
        use std::os::unix::fs::FileExt;
        let len = N_EPERM / 5;
        buf.resize(len, 0);
        self.file
            .read_exact_at(buf, self.offset + (row * len) as u64)
            .map_err(|e| P1Error::Io(format!("read giant row {row}"), e))
    }
}
