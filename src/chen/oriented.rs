//! The oriented phase-1 targets: every class-A hexagon one colour *and unturned*.
//!
//! The block split's phase 1 ([`super::relaxed`], the strict `H'` table) takes the hexagons up to
//! turning each in place, so its states are edge cosets (the 81 turns of the four triples merged).
//! Here a hexagon must also be unturned: the target is the solved configuration (A-centres and
//! edges) in each of the 12 frames ([`Target::Oriented`]), or that and the 8 configurations one
//! class-A turn from it ([`Target::T1`], 9 per frame). Neither is closed under turning a hexagon,
//! so a state is the whole edge permutation: 369,600 A-centre colourings x 239,500,800 even edge
//! permutations = 8.85e13 states.
//!
//! **Symmetry.** Both targets are closed under the 288 pairs of [`crate::sym`] (relabel by a
//! class-preserving symmetry `g`, reposition by `t`, both proper or both improper): `t g^-1` is a
//! rotation, and a pair is that rotation times the conjugation by `g`, which fixes the solved
//! configuration and maps face turns to face turns (class-A turns to class-A turns). So the
//! distance is a class function, and a table needs one row per A-colouring class (1,396, the
//! classes of the two-phase table) of all 239,500,800 edge permutations: 3.34e11 slots. A slot
//! `(row, col)` stands for the state (representative colouring of `row`, edges `col`); every
//! state of a colouring in the class maps to one slot of the row (through the pair that takes the
//! colouring to the representative), so a row's count times its class size counts states.
//!
//! **Edge column.** `col = (e * 8 + b) * 81 + r`: `e` the triple labelling (which A-triple's
//! piece is at each position, [`N_MSET`]), `b` three chirality bits (the fourth follows), as in the
//! two-phase table, and `r` (base 3, digit `t` for triple `t`) the index in [`triple_pieces`] of
//! the piece at the lowest position of triple `t`. Relabelling a triple's pieces cyclically
//! adds one to its digit; a move adds to each digit a number that depends only on `(e, b)` and
//! the move; a pair maps the digits by a signed permutation (from `g`) plus a number that
//! depends on `(e, b)`. So the 81 slots of a coset `(e, b)` move together as a block of 81 bits,
//! by an affine map of `Z_3^4` ([`P1Space::step`], [`P1Space::permute`]).

use crate::coord::{edge_chirality, move_tables, mset, triple_pieces, MoveTables};
use crate::coord::{N_ECOSET, N_MSET};
use crate::cube::{defs, is_a_move, Cube, N_MOVES};
use crate::sym::{edge_syms, syms, EdgeSyms, N_GEO};
use crate::tables::AeSym;
use std::collections::BTreeSet;

/// Turns of the four triples (relabellings) per edge coset: `3^4`.
pub const N_TURN: usize = 81;
/// Edge permutations: `12! / 2`.
pub const N_EPERM: usize = N_ECOSET * N_TURN;
/// Bytes of the 8 cosets (81 bits each) of one triple labelling.
pub const GROUP_BYTES: usize = 81;
/// Bytes of a row of bits (one per edge permutation).
pub const ROW_BYTES: usize = N_MSET * GROUP_BYTES;
/// A row on disk and in memory: whole 4 KiB pages (for direct I/O), with room to read 16 bytes
/// at any bit of the row.
pub const ROW_STRIDE: usize = ROW_BYTES.div_ceil(4096) * 4096;
const _: () = assert!(ROW_STRIDE >= ROW_BYTES + 16);

/// The phase-1 target set.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Target {
    /// the solved configuration in each frame (12 states)
    Oriented,
    /// the oriented configurations and those one class-A turn from them (108 states)
    T1,
}

impl Target {
    pub fn name(self) -> &'static str {
        match self {
            Target::Oriented => "oriented",
            Target::T1 => "t1",
        }
    }
}

/// A construction check failed (a bug).
#[derive(Debug, Clone, thiserror::Error)]
#[error("oriented phase-1 space: {0}")]
pub struct OrientedError(pub String);

/// The neighbour of a row by a move: the row of the moved representative and the pair that
/// takes the moved colouring to that row's representative.
#[derive(Clone, Copy, Debug, Default)]
pub struct Nbr {
    pub row: u32,
    pub g: u8,
    pub t: u8,
}

/// The labelling part of a pull step ([`P1Space::step_e`]).
#[derive(Clone, Copy, Debug)]
pub struct StepE {
    /// the source triple labelling
    pub e3: u32,
    f1: u8,
    j1: u8,
    f3: u8,
    j3: u8,
}

/// Rows, move and symmetry tables of the oriented phase-1 space.
pub struct P1Space {
    /// the A-colouring classes (rows) and their stabilisers
    pub ae: AeSym,
    /// colourings per class
    pub weight: Vec<u64>,
    /// per row and move
    pub nbr: Vec<[Nbr; N_MOVES]>,
    mt: &'static MoveTables,
    es: &'static EdgeSyms,
    /// per triple labelling and move: the old slot (0..3, ascending positions) of each triple
    /// that lands lowest, two bits per triple
    jm: Vec<[u8; N_MOVES]>,
    /// the same for the position map of each symmetry `t`: `js[t * N_MSET + e]`
    js: Vec<u8>,
    /// the linear part of each `g` on the turn digits
    pub lin: [[u8; N_TURN]; N_GEO],
    /// the constant part of each `g`
    cst: [u8; N_GEO],
    /// whether `g` maps every digit to itself
    lin_id: [bool; N_GEO],
    /// the turn added by a move or repositioning: `delta[bits][j]`
    delta: Box<[[u8; 256]; 16]>,
    /// addition in `Z_3^4`
    pub add: Box<[[u8; N_TURN]; N_TURN]>,
    /// per `g`: the linear part on an 81-bit block, a byte at a time (`(g * 11 + byte) * 256 +
    /// value`)
    ptab: Vec<u128>,
    /// digit masks of an 81-bit block: `dmask[i][v]` the bits whose digit `i` is `v`
    dmask: [[u128; 3]; 4],
}

/// The digits of a turn index.
#[inline]
pub fn digits(r: u8) -> [u8; 4] {
    [r % 3, r / 3 % 3, r / 9 % 3, r / 27]
}

#[inline]
fn undigits(d: [u8; 4]) -> u8 {
    d[0] + 3 * d[1] + 9 * d[2] + 27 * d[3]
}

/// Bit mask of the low 81 bits.
pub const MASK81: u128 = (1u128 << 81) - 1;

/// The packed lowest-slot indices of a position map `dest` for triple labelling `l`.
fn low_slots(l: &[u8; 12], dest: &[u8; 12]) -> u8 {
    let mut best = [(u8::MAX, 0u8); 4];
    let mut n = [0u8; 4];
    for p in 0..12 {
        let t = l[p] as usize;
        if dest[p] < best[t].0 {
            best[t] = (dest[p], n[t]);
        }
        n[t] += 1;
    }
    (0..4).fold(0, |a, t| a | best[t].1 << (2 * t))
}

impl P1Space {
    /// Builds the tables (a few seconds) and checks them against the definitions.
    ///
    /// # Errors
    ///
    /// If a construction check fails.
    pub fn new() -> Result<P1Space, OrientedError> {
        let err = |s: String| OrientedError(s);
        let d = defs();
        let s = syms();
        let ms = mset();
        let mt = move_tables();
        let es = edge_syms();
        let ae = AeSym::new(true);
        let mut weight = vec![0u64; ae.cls.rep.len()];
        for &v in ae.cls.sym.iter() {
            weight[(v & 0xffff) as usize] += 1;
        }
        let nbr: Vec<[Nbr; N_MOVES]> = ae
            .cls
            .rep
            .iter()
            .map(|&r| {
                std::array::from_fn(|m| {
                    let v = ae.cls.sym[mt.a[r as usize][m] as usize];
                    Nbr {
                        row: v & 0xffff,
                        g: (v >> 16 & 0xff) as u8,
                        t: (v >> 24) as u8,
                    }
                })
            })
            .collect();
        // where each position goes under each move (`moves[m].ep` is a gather)
        let mdest: [[u8; 12]; N_MOVES] = std::array::from_fn(|m| {
            let mut dest = [0u8; 12];
            for p in 0..12 {
                dest[d.moves[m].ep[p] as usize] = p as u8;
            }
            dest
        });
        let jm: Vec<[u8; N_MOVES]> = (0..N_MSET as u32)
            .map(|e| {
                let l = ms.unrank(e);
                std::array::from_fn(|m| low_slots(&l, &mdest[m]))
            })
            .collect();
        let mut js = vec![0u8; N_GEO * N_MSET];
        for t in 0..N_GEO {
            for e in 0..N_MSET {
                js[t * N_MSET + e] = low_slots(&ms.unrank(e as u32), &s.e_pos[t]);
            }
        }
        // relabelling: triple t's piece i becomes triple a_col[g][t]'s piece psi(i)
        let tp = triple_pieces();
        let mut lin = [[0u8; N_TURN]; N_GEO];
        let mut cst = [0u8; N_GEO];
        let mut lin_id = [false; N_GEO];
        for g in 0..N_GEO {
            let mut sign = [1u8; 4];
            let mut to = [0usize; 4];
            let mut c = [0u8; 4];
            for t in 0..4 {
                let t2 = s.a_col[g][t] as usize;
                let psi: Vec<u8> = (0..3)
                    .map(|i| {
                        let k = s.e_piece[g][tp[t][i] as usize];
                        tp[t2]
                            .iter()
                            .position(|&x| x == k)
                            .map_or(u8::MAX, |j| j as u8)
                    })
                    .collect();
                if psi.contains(&u8::MAX) {
                    return Err(err(format!("g {g}: triple {t} is not mapped to a triple")));
                }
                let even = (psi[1] + 3 - psi[0]) % 3 == 1;
                if (psi[2] + 3 - psi[1]) % 3 != if even { 1 } else { 2 } {
                    return Err(err(format!("g {g}: triple {t}: {psi:?}")));
                }
                to[t] = t2;
                sign[t] = if even { 1 } else { 2 };
                c[t2] = psi[0];
            }
            for r in 0..N_TURN as u8 {
                let dg = digits(r);
                let mut o = [0u8; 4];
                for t in 0..4 {
                    o[to[t]] = dg[t] * sign[t] % 3;
                }
                lin[g][r as usize] = undigits(o);
            }
            cst[g] = undigits(c);
            lin_id[g] = (0..4).all(|t| to[t] == t && sign[t] == 1);
        }
        let mut delta = Box::new([[0u8; 256]; 16]);
        for bits in 0..16 {
            for j in 0..256usize {
                let dg: [u8; 4] = std::array::from_fn(|t| {
                    let jt = (j >> (2 * t) & 3) as u8 % 3;
                    if bits >> t & 1 == 1 {
                        (3 - jt) % 3
                    } else {
                        jt
                    }
                });
                delta[bits][j] = undigits(dg);
            }
        }
        let mut add = Box::new([[0u8; N_TURN]; N_TURN]);
        for a in 0..N_TURN as u8 {
            for b in 0..N_TURN as u8 {
                let (x, y) = (digits(a), digits(b));
                add[a as usize][b as usize] = undigits(std::array::from_fn(|i| (x[i] + y[i]) % 3));
            }
        }
        let mut dmask = [[0u128; 3]; 4];
        for r in 0..N_TURN as u8 {
            let dg = digits(r);
            for i in 0..4 {
                dmask[i][dg[i] as usize] |= 1u128 << r;
            }
        }
        // dest[r] = h[lin(r)]: bit u of h goes to every r with lin(r) = u (one, lin is a bijection)
        let mut ptab = vec![0u128; N_GEO * 11 * 256];
        for g in 0..N_GEO {
            let mut inv = [0u8; N_TURN];
            for r in 0..N_TURN {
                inv[lin[g][r] as usize] = r as u8;
            }
            for i in 0..11 {
                for v in 0..256usize {
                    let mut o = 0u128;
                    for k in 0..8 {
                        let u = 8 * i + k;
                        if v >> k & 1 == 1 && u < N_TURN {
                            o |= 1u128 << inv[u];
                        }
                    }
                    ptab[(g * 11 + i) * 256 + v] = o;
                }
            }
        }
        let sp = P1Space {
            ae,
            weight,
            nbr,
            mt,
            es,
            jm,
            js,
            lin,
            cst,
            lin_id,
            delta,
            add,
            ptab,
            dmask,
        };
        sp.check()?;
        Ok(sp)
    }

    /// Rows (A-colouring classes).
    pub fn rows(&self) -> usize {
        self.ae.cls.rep.len()
    }

    /// Slots: rows x edge permutations.
    pub fn slots(&self) -> u64 {
        self.rows() as u64 * N_EPERM as u64
    }

    /// The fourth chirality bit of a coset.
    #[inline(always)]
    pub fn bits4(&self, e: u32, b3: u8) -> u8 {
        let bit3 = self.mt.eref[e as usize] ^ (b3 ^ b3 >> 1 ^ b3 >> 2) & 1;
        b3 | bit3 << 3
    }

    /// The source coset and turn of a pull: the state (rep of a row, coset `(e, bits)`, turn
    /// `r`) moved by `m` and taken through the pair `(g, t)` is the coset `(e3, bits3)` with
    /// turn `lin[g][r] + tr`. Returns `(e3, bits3, tr)`.
    #[inline(always)]
    pub fn step(&self, m: usize, g: usize, t: usize, e: u32, bits: u8) -> (u32, u8, u8) {
        let pe = self.mt.e[e as usize][m];
        let e1 = pe & 0xfffff;
        let bits1 = bits ^ (pe >> 20) as u8;
        let d1 = self.delta[bits as usize][self.jm[e as usize][m] as usize];
        let e2 = self.es.relab_e[g * N_MSET + e1 as usize];
        let bits2 = self.es.relab_bits[g][bits1 as usize];
        let p = self.es.repos_e[t * N_MSET + e2 as usize];
        let e3 = p & 0xfffff;
        let bits3 = bits2 ^ (p >> 20) as u8;
        let d3 = self.delta[bits2 as usize][self.js[t * N_MSET + e2 as usize] as usize];
        let tr = self.add[self.lin[g][d1 as usize] as usize]
            [self.add[self.cst[g] as usize][d3 as usize] as usize];
        (e3, bits3, tr)
    }

    /// The part of [`Self::step`] that depends only on the triple labelling `e` (shared by its
    /// 8 cosets).
    #[inline(always)]
    pub fn step_e(&self, m: usize, g: usize, t: usize, e: u32) -> StepE {
        let pe = self.mt.e[e as usize][m];
        let e2 = self.es.relab_e[g * N_MSET + (pe & 0xfffff) as usize];
        let p = self.es.repos_e[t * N_MSET + e2 as usize];
        StepE {
            e3: p & 0xfffff,
            f1: (pe >> 20) as u8,
            j1: self.jm[e as usize][m],
            f3: (p >> 20) as u8,
            j3: self.js[t * N_MSET + e2 as usize],
        }
    }

    /// The rest of [`Self::step`] for one coset of `e`: `(bits3, tr)`.
    #[inline(always)]
    pub fn step_b(&self, se: &StepE, g: usize, bits: u8) -> (u8, u8) {
        let d1 = self.delta[bits as usize][se.j1 as usize];
        let bits2 = self.es.relab_bits[g][(bits ^ se.f1) as usize];
        let d3 = self.delta[bits2 as usize][se.j3 as usize];
        let tr = self.add[self.lin[g][d1 as usize] as usize]
            [self.add[self.cst[g] as usize][d3 as usize] as usize];
        (bits2 ^ se.f3, tr)
    }

    /// The block `dest` with `dest[r] = src[lin[g][r] + tr]` (bits 0..81).
    #[inline(always)]
    pub fn permute(&self, g: usize, tr: u8, src: u128) -> u128 {
        let mut h = src & MASK81;
        if tr != 0 {
            let dg = digits(tr);
            for i in 0..4 {
                let s = [1u32, 3, 9, 27][i];
                let dm = &self.dmask[i];
                h = match dg[i] {
                    1 => (h >> s) & (dm[0] | dm[1]) | (h << (2 * s)) & dm[2],
                    2 => (h >> (2 * s)) & dm[0] | (h << s) & (dm[1] | dm[2]),
                    _ => h,
                };
            }
        }
        if self.lin_id[g] {
            return h;
        }
        let pt = &self.ptab[g * 11 * 256..(g + 1) * 11 * 256];
        let mut o = 0u128;
        for i in 0..11 {
            o |= pt[i * 256 + (h >> (8 * i) & 0xff) as usize];
        }
        o
    }

    /// The slot of the neighbour of slot `(row, col)` by move `m`.
    #[inline]
    pub fn nbr_slot(&self, row: usize, col: u32, m: usize) -> (usize, u32) {
        let n = self.nbr[row][m];
        let (k, r) = (col / N_TURN as u32, (col % N_TURN as u32) as u8);
        let (e, b3) = (k / 8, (k % 8) as u8);
        let (e3, bits3, tr) = self.step(m, n.g as usize, n.t as usize, e, self.bits4(e, b3));
        let r3 = self.add[self.lin[n.g as usize][r as usize] as usize][tr as usize];
        (
            n.row as usize,
            (e3 * 8 + u32::from(bits3 & 7)) * N_TURN as u32 + u32::from(r3),
        )
    }

    /// The slot of a state (A-centres and edges; the rest is ignored).
    pub fn slot_of(&self, xa: &[u8; 12], ep: &[u8; 12]) -> (usize, u32) {
        let s = syms();
        let v = self.ae.cls.sym[mset().rank(xa) as usize];
        let (g, t) = ((v >> 16 & 0xff) as usize, (v >> 24) as usize);
        let mut n = [0u8; 12];
        for p in 0..12 {
            n[s.e_pos[t][p] as usize] = s.e_piece[g][ep[p] as usize];
        }
        ((v & 0xffff) as usize, encode_ep(&n))
    }

    /// The state a slot stands for: (representative colouring, edges).
    pub fn state_of(&self, row: usize, col: u32) -> ([u8; 12], [u8; 12]) {
        (mset().unrank(self.ae.cls.rep[row]), decode_ep(col))
    }

    /// The states of a target set in the fixed frame: the images of the solved configuration
    /// (and, for [`Target::T1`], of its class-A turns) under the 288 pairs.
    pub fn target_states(&self, target: Target) -> Vec<([u8; 12], [u8; 12])> {
        let s = syms();
        let s0 = Cube::solved();
        let mut base = vec![(s0.xa, s0.ep)];
        if target == Target::T1 {
            for m in (0..N_MOVES).filter(|&m| is_a_move(m)) {
                let c = s0.apply(m);
                base.push((c.xa, c.ep));
            }
        }
        let mut out = BTreeSet::new();
        for &(xa, ep) in &base {
            for &(g, t) in &s.pairs {
                let a = s.apply_a(g as usize, t as usize, &xa);
                let mut n = [0u8; 12];
                for p in 0..12 {
                    n[s.e_pos[t as usize][p] as usize] = s.e_piece[g as usize][ep[p] as usize];
                }
                out.insert((a, n));
            }
        }
        out.into_iter().collect()
    }

    /// The target slots: every slot whose state is a target state.
    pub fn target_slots(&self, target: Target) -> Vec<(usize, u32)> {
        let ms = mset();
        let mut v: Vec<(usize, u32)> = self
            .target_states(target)
            .iter()
            .filter_map(|(xa, ep)| {
                let r = ms.rank(xa);
                let row = (self.ae.cls.sym[r as usize] & 0xffff) as usize;
                (self.ae.cls.rep[row] == r).then(|| (row, encode_ep(ep)))
            })
            .collect();
        v.sort_unstable();
        v
    }

    /// Construction checks: the tables against the definitions on a sample of cosets.
    fn check(&self) -> Result<(), OrientedError> {
        let ms = mset();
        // every 997th coset, every move, a few turns: the table step equals the cube move
        for k in (0..N_ECOSET as u32).step_by(997) {
            let (e, b3) = (k / 8, (k % 8) as u8);
            for r in [0u8, 40, 80] {
                let col = k * N_TURN as u32 + u32::from(r);
                let ep = decode_ep(col);
                if encode_ep(&ep) != col {
                    return Err(OrientedError(format!("column {col} does not round-trip")));
                }
                for row in [0usize, self.rows() / 2, self.rows() - 1] {
                    let xa = ms.unrank(self.ae.cls.rep[row]);
                    let c = Cube {
                        ep,
                        xa,
                        ..Cube::solved()
                    };
                    for m in 0..N_MOVES {
                        let n = c.apply(m);
                        let want = self.slot_of(&n.xa, &n.ep);
                        let got = self.nbr_slot(row, col, m);
                        if want != got {
                            return Err(OrientedError(format!(
                                "row {row} coset {e}/{b3} turn {r} move {m}: {got:?} != {want:?}"
                            )));
                        }
                    }
                }
            }
        }
        Ok(())
    }
}

/// The column of an edge permutation.
pub fn encode_ep(ep: &[u8; 12]) -> u32 {
    let d = defs();
    let tp = triple_pieces();
    let l: [u8; 12] = std::array::from_fn(|p| d.edge_atri[ep[p] as usize]);
    let e = mset().rank(&l);
    let bits = edge_chirality(ep);
    let mut dg = [u8::MAX; 4];
    for p in 0..12 {
        let t = l[p] as usize;
        if dg[t] == u8::MAX {
            dg[t] = tp[t].iter().position(|&x| x == ep[p]).unwrap_or(0) as u8;
        }
    }
    (e * 8 + u32::from(bits & 7)) * N_TURN as u32 + u32::from(undigits(dg))
}

/// The edge permutation of a column.
pub fn decode_ep(col: u32) -> [u8; 12] {
    let (l, bits) = coset_parts(col / N_TURN as u32);
    ep_of(&l, bits, (col % N_TURN as u32) as u8)
}

/// The triple labelling and the four chirality bits of coset `k` (shared by its 81 columns).
pub fn coset_parts(k: u32) -> ([u8; 12], u8) {
    let mt = move_tables();
    let (e, b3) = (k / 8, (k % 8) as u8);
    let bit3 = mt.eref[e as usize] ^ (b3 ^ b3 >> 1 ^ b3 >> 2) & 1;
    (mset().unrank(e), b3 | bit3 << 3)
}

/// The edge permutation of a coset's column with turn `r` ([`coset_parts`]).
pub fn ep_of(l: &[u8; 12], bits: u8, r: u8) -> [u8; 12] {
    let tp = triple_pieces();
    let dg = digits(r);
    let mut ep = [0u8; 12];
    let mut n = [0u8; 4];
    for p in 0..12 {
        let t = l[p] as usize;
        let j = n[t];
        n[t] += 1;
        let i = if bits >> t & 1 == 1 {
            (dg[t] + 3 - j) % 3
        } else {
            (dg[t] + j) % 3
        };
        ep[p] = tp[t][i as usize];
    }
    ep
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cube::{inverse_move, perm_parity};
    use rand::{Rng, SeedableRng};

    fn space() -> &'static P1Space {
        static S: std::sync::OnceLock<P1Space> = std::sync::OnceLock::new();
        S.get_or_init(|| P1Space::new().unwrap())
    }

    #[test]
    fn sizes() {
        let sp = space();
        assert_eq!(sp.rows(), 1396);
        assert_eq!(sp.weight.iter().sum::<u64>(), N_MSET as u64);
        assert_eq!(N_EPERM, 239_500_800);
        assert_eq!(sp.target_states(Target::Oriented).len(), 12);
        assert_eq!(sp.target_states(Target::T1).len(), 108);
    }

    #[test]
    fn columns_round_trip_and_are_even_permutations() {
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(5);
        for _ in 0..200_000 {
            let col = rng.gen_range(0..N_EPERM as u32);
            let ep = decode_ep(col);
            assert_eq!(perm_parity(&ep), 0);
            let mut seen = [false; 12];
            for &k in &ep {
                seen[k as usize] = true;
            }
            assert!(seen.iter().all(|&x| x));
            assert_eq!(encode_ep(&ep), col);
        }
        assert_eq!(encode_ep(&Cube::solved().ep) % N_TURN as u32, 0);
    }

    #[test]
    fn steps_follow_cube_moves() {
        let sp = space();
        let ms = mset();
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(7);
        for _ in 0..20_000 {
            let row = rng.gen_range(0..sp.rows());
            let col = rng.gen_range(0..N_EPERM as u32);
            let (xa, ep) = sp.state_of(row, col);
            assert_eq!(ms.rank(&xa), sp.ae.cls.rep[row]);
            let c = Cube {
                ep,
                xa,
                ..Cube::solved()
            };
            for m in 0..N_MOVES {
                let n = c.apply(m);
                assert_eq!(sp.nbr_slot(row, col, m), sp.slot_of(&n.xa, &n.ep));
            }
        }
    }

    #[test]
    fn split_step_equals_step() {
        let sp = space();
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(13);
        for _ in 0..100_000 {
            let (m, g, t) = (
                rng.gen_range(0..N_MOVES),
                rng.gen_range(0..N_GEO),
                rng.gen_range(0..N_GEO),
            );
            let e = rng.gen_range(0..N_MSET as u32);
            let bits = sp.bits4(e, rng.gen_range(0..8));
            let se = sp.step_e(m, g, t, e);
            let (b3, tr) = sp.step_b(&se, g, bits);
            assert_eq!((se.e3, b3, tr), sp.step(m, g, t, e, bits));
        }
    }

    #[test]
    fn permute_matches_digit_map() {
        let sp = space();
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(9);
        for g in 0..N_GEO {
            for tr in 0..N_TURN as u8 {
                let src: u128 = rng.gen::<u128>() & MASK81;
                let got = sp.permute(g, tr, src);
                for r in 0..N_TURN {
                    let u = sp.add[sp.lin[g][r] as usize][tr as usize];
                    assert_eq!(got >> r & 1, src >> u & 1, "g {g} tr {tr} r {r}");
                }
                assert_eq!(got & !MASK81, 0);
            }
        }
    }

    #[test]
    fn pairs_map_moves_to_moves_and_fix_the_targets() {
        // a pair takes the neighbours of a state to the neighbours of its image, up to a
        // rotation (a pair (identity, rotation)), so it preserves the distance to a rotation-closed
        // target set; the targets are closed under the pairs
        let s = syms();
        let sp = space();
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(11);
        let apply = |g: usize, t: usize, xa: &[u8; 12], ep: &[u8; 12]| {
            let a = s.apply_a(g, t, xa);
            let mut n = [0u8; 12];
            for p in 0..12 {
                n[s.e_pos[t][p] as usize] = s.e_piece[g][ep[p] as usize];
            }
            (a, n)
        };
        // the rotations: pairs (identity colour map, proper t)
        let id = (0..N_GEO)
            .find(|&g| {
                s.a_col[g] == [0, 1, 2, 3] && s.e_piece[g] == std::array::from_fn(|i| i as u8)
            })
            .unwrap();
        let rots: Vec<usize> = s
            .pairs
            .iter()
            .filter(|p| p.0 as usize == id)
            .map(|p| p.1 as usize)
            .collect();
        assert_eq!(rots.len(), 12);
        for _ in 0..300 {
            let c = Cube::random(&mut rng);
            for &(g, t) in &s.pairs {
                let (g, t) = (g as usize, t as usize);
                let (ya, ye) = apply(g, t, &c.xa, &c.ep);
                let img: BTreeSet<_> = (0..N_MOVES)
                    .map(|m| {
                        let n = c.apply(m);
                        apply(g, t, &n.xa, &n.ep)
                    })
                    .collect();
                let y = Cube {
                    ep: ye,
                    xa: ya,
                    ..c
                };
                // every image neighbour is a neighbour of y up to a rotation
                for z in &img {
                    let ok = (0..N_MOVES).any(|m| {
                        let n = y.apply(m);
                        rots.iter().any(|&r| apply(id, r, &n.xa, &n.ep) == *z)
                    });
                    assert!(ok);
                }
            }
        }
        for target in [Target::Oriented, Target::T1] {
            let ts: BTreeSet<_> = sp.target_states(target).into_iter().collect();
            for &(xa, ep) in &ts {
                for &(g, t) in &s.pairs {
                    assert!(ts.contains(&apply(g as usize, t as usize, &xa, &ep)));
                }
                // T1 is the oriented states and their neighbours by class-A moves
                if target == Target::T1 {
                    let c = Cube {
                        ep,
                        xa,
                        ..Cube::solved()
                    };
                    let o: BTreeSet<_> = sp.target_states(Target::Oriented).into_iter().collect();
                    let near = o.contains(&(xa, ep))
                        || (0..N_MOVES).filter(|&m| is_a_move(m)).any(|m| {
                            let n = c.apply(inverse_move(m));
                            o.contains(&(n.xa, n.ep))
                        });
                    assert!(near);
                }
            }
        }
    }

    #[test]
    fn target_slots_count_the_target_states() {
        let sp = space();
        for (target, n) in [(Target::Oriented, 12u64), (Target::T1, 108)] {
            let slots = sp.target_slots(target);
            let w: u64 = slots.iter().map(|&(row, _)| sp.weight[row]).sum();
            // a row's slots that are targets, times its class size, counts the target states
            // (the stabiliser images of a target slot are target slots too)
            assert_eq!(w, n, "{target:?}");
        }
    }
}
