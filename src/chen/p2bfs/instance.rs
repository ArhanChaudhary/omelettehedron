//! The phase-2 state graph as the breadth-first search indexes it.
//!
//! A state is seen in the pattern view ([`super::super::space`]: the 8 face turns, nothing
//! rotates the puzzle) and reduced by the 24 symmetries that keep the face classes: its
//! configuration is conjugated to the representative of its symmetry class (1,611 classes),
//! and its B-centre colouring and corners with it. A *slot* is then (class `c`, B-centre
//! colouring `x`, corner state `p`); a *block* is the `n_p` slots of one `(c, x)`, kept as a
//! bitmap of `words` 64-bit words (11,520 corner states: 180 words, 1,440 bytes).
//!
//! The full instance has 1,611 x 369,600 x 11,520 = 6.86 x 10^12 slots for the 1.59 x 10^14
//! positions of phase 2 (37,260 configurations x 11,520 corner states x 369,600 colourings).
//! Smaller instances replace the corners or the colourings by a quotient that commutes with
//! the moves and with the symmetries ([`CornerView`], [`BView`]); they are the same graph
//! seen coarser, used to check the search against independent ones.
//!
//! A configuration fixed by some symmetries (its class's stabiliser, [`Instance::stab`]) has
//! several slots for one position: (c, x, p) and its images under the stabiliser. The search
//! sets all of them, so that the set of each class is closed under its stabiliser; a class
//! then holds `24 / (1 + stab.len())` configurations, each with the class's set of slots.

use super::super::space::{space, Space, N_PMOVES};
use crate::coord::{mset, N_MSET};
use crate::cube::{Cube, N_MOVES};
use crate::fast::shuffle;
use std::collections::HashMap;

/// Symmetries that keep the face classes.
pub const N_SYM: usize = 24;
/// No move: a move that is not allowed.
pub const NO_EDGE: u16 = u16::MAX;

/// What the search keeps of the corners.
#[derive(Clone, Copy, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum CornerView {
    /// all 11,520 pattern-view corner states
    Full,
    /// the corner permutation only (360)
    Perm,
    /// nothing (1)
    None,
}

/// What the search keeps of the B-centres.
#[derive(Clone, Copy, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum BView {
    /// all 369,600 colourings
    Full,
    /// the colouring up to a renaming of the 4 colours (15,400 partitions into 4 triples)
    Unlabelled,
    /// nothing (1)
    None,
}

/// Which instance: the full one or a quotient.
#[derive(Clone, Copy, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct Kind {
    /// corners kept
    pub corners: CornerView,
    /// B-centres kept
    pub b: BView,
}

impl Kind {
    /// The full phase-2 space.
    pub const FULL: Kind = Kind {
        corners: CornerView::Full,
        b: BView::Full,
    };
}

/// The instance's tables.
pub struct Instance {
    /// which instance
    pub kind: Kind,
    /// symmetry classes of pattern-view configurations
    pub n_class: usize,
    /// B-centre states (blocks per class)
    pub n_x: usize,
    /// corner states (bits per block)
    pub n_p: usize,
    /// 64-bit words per block
    pub words: usize,
    /// per class and move: the class of the child's configuration (or [`NO_EDGE`]) and the
    /// symmetry that takes the child to it
    pub edges: Vec<[(u16, u8); N_MOVES]>,
    /// per B-centre state: the state after each move
    pub xm: Vec<[u32; N_MOVES]>,
    /// per symmetry and B-centre state (`k * n_x + x`): the state conjugated
    pub xs: Vec<u32>,
    /// per move, symmetry and corner state (`(m * N_SYM + k) * n_p + p`): the corner state
    /// after the move, conjugated by the symmetry
    pub pperm: Vec<u16>,
    /// per corner state and symmetry (`p * N_SYM + k`): the state conjugated
    pub ps: Vec<u16>,
    /// per class: the symmetries other than the identity that fix its representative
    pub stab: Vec<Vec<u8>>,
    /// per class: whether its configurations are in `H'` (every hexagon one colour)
    pub hprime: Vec<bool>,
    /// per class: its representative (pattern-view configuration)
    pub rep: Vec<u16>,
    /// the solved slot (class, x, p)
    pub solved: (u32, u32, u32),
    /// pattern-view configurations, for the totals (37,260)
    pub n_config: usize,
    /// the classes in the order the search splits them into parts: the classes joined by
    /// U, F, BR and BL turns (always allowed) together, so that most moves stay in a part
    pub order: Vec<u32>,
    /// the position of each class in [`Instance::order`]
    pub pos: Vec<u32>,
    /// full corner view: the pattern-view corner state ([`super::super::corners::PCorners`])
    /// of each corner state of the instance
    pub p_to_pc: Vec<u16>,
    /// full corner view: the block permutations a word at a time
    pub word: Option<WordPerm>,
}

impl Instance {
    /// Build the tables of an instance (a second or two, on all threads).
    pub fn new(kind: Kind) -> Instance {
        let s = space();
        let pat = &s.pat;
        let n_class = pat.rep.len();
        let edges: Vec<[(u16, u8); N_MOVES]> = (0..n_class)
            .map(|c| {
                let q = pat.rep[c] as usize;
                std::array::from_fn(|m| {
                    let n = pat.next[q][m];
                    if n == u16::MAX {
                        (NO_EDGE, 0)
                    } else {
                        let (c2, k) = pat.class[n as usize];
                        (c2, k)
                    }
                })
            })
            .collect();
        let (n_x, xm, xs, x_of_solved) = b_tables(s, kind.b);
        let (n_p, pnext, ps, p_of_solved, p_to_pc) = corner_tables(s, kind.corners);
        let mut pperm = vec![0u16; N_MOVES * N_SYM * n_p];
        for m in 0..N_MOVES {
            for k in 0..N_SYM {
                let base = (m * N_SYM + k) * n_p;
                for p in 0..n_p {
                    let after = pnext[p * N_MOVES + m] as usize;
                    pperm[base + p] = ps[after * N_SYM + k];
                }
            }
        }
        let word = (kind.corners == CornerView::Full).then(|| WordPerm::new(&pperm, n_p));
        let (c0, _) = pat.class[0];
        let order = a_move_order(&edges, &pat.pmove);
        let mut pos = vec![0u32; n_class];
        for (i, &c) in order.iter().enumerate() {
            pos[c as usize] = i as u32;
        }
        let hprime = (0..n_class)
            .map(|c| s.pshape[pat.cfg[pat.rep[c] as usize] as usize] == 0)
            .collect();
        Instance {
            kind,
            n_class,
            n_x,
            n_p,
            words: n_p.div_ceil(64),
            edges,
            xm,
            xs,
            pperm,
            ps,
            stab: pat.stab.clone(),
            hprime,
            rep: pat.rep.clone(),
            solved: (u32::from(c0), x_of_solved, p_of_solved),
            n_config: pat.cfg.len(),
            order,
            pos,
            p_to_pc,
            word,
        }
    }

    /// The fiber permutation of move `m` followed by symmetry `k`.
    #[inline]
    pub fn fiber_perm(&self, m: usize, k: usize) -> &[u16] {
        let base = (m * N_SYM + k) * self.n_p;
        &self.pperm[base..base + self.n_p]
    }

    /// Configurations in class `c`.
    #[inline]
    pub fn class_size(&self, c: usize) -> u64 {
        (N_SYM / (1 + self.stab[c].len())) as u64
    }

    /// Slots of the instance.
    pub fn n_slots(&self) -> u64 {
        (self.n_class * self.n_x * self.n_p) as u64
    }

    /// Positions of the instance (configurations x B-centre states x corner states).
    pub fn n_positions(&self) -> u128 {
        self.n_config as u128 * self.n_x as u128 * self.n_p as u128
    }

    /// Positions up to the 24 symmetries (by fixed points, per class under its stabiliser).
    pub fn n_orbits(&self) -> u128 {
        let fix_x: Vec<u128> = (0..N_SYM)
            .map(|k| {
                (0..self.n_x)
                    .filter(|&x| self.xs[k * self.n_x + x] as usize == x)
                    .count() as u128
            })
            .collect();
        let fix_p: Vec<u128> = (0..N_SYM)
            .map(|k| {
                (0..self.n_p)
                    .filter(|&p| self.ps[p * N_SYM + k] as usize == p)
                    .count() as u128
            })
            .collect();
        (0..self.n_class)
            .map(|c| {
                let fixed: u128 = self.stab[c]
                    .iter()
                    .map(|&k| fix_x[k as usize] * fix_p[k as usize])
                    .sum();
                (self.n_x as u128 * self.n_p as u128 + fixed) / (1 + self.stab[c].len() as u128)
            })
            .sum()
    }

    /// The mask of the valid bits of the last word of a block.
    #[inline]
    pub fn last_mask(&self) -> u64 {
        let r = self.n_p % 64;
        if r == 0 {
            u64::MAX
        } else {
            (1u64 << r) - 1
        }
    }
}

/// The full corner view's numbering: `permutation id * 32 + orientation bits of places
/// 0..5` for each pattern-view state, and its inverse.
fn corner_numbering(states: &[([u8; 6], [u8; 6])]) -> (Vec<u16>, Vec<u16>) {
    let mut perms: HashMap<[u8; 6], u16> = HashMap::new();
    let mut of_state = vec![0u16; states.len()];
    let mut back = vec![0u16; states.len()];
    for (st, (cp, co)) in states.iter().enumerate() {
        let next = perms.len() as u16;
        let id = *perms.entry(*cp).or_insert(next);
        let code = (0..5).fold(0u16, |a, v| a | u16::from(co[v] & 1) << v);
        let p = id * 32 + code;
        of_state[st] = p;
        if let Some(b) = back.get_mut(p as usize) {
            *b = st as u16;
        }
    }
    (of_state, back)
}

/// The permutation of a block (11,520 bits: 360 32-bit halves, one per corner permutation)
/// by a move and a symmetry, a half at a time: each half goes to one half, its bits permuted by
/// an affine map of their index (orientations: a linear map, then a flip). Built from the
/// per-bit tables, and checked against them; a (move, symmetry) whose map is not of that form
/// keeps the per-bit path.
pub struct WordPerm {
    /// per (move, symmetry): whether the word path holds
    pub ok: Vec<bool>,
    /// per (move, symmetry): the linear part as byte lookups (`[byte][value]`)
    lut: Vec<[[u32; 256]; 4]>,
    /// per (move, symmetry) and permutation (`mk * 360 + perm`): the target half and the flip
    to: Vec<(u16, u8)>,
}

/// The halves (corner permutations) of a full block.
pub const HALVES: usize = 360;

impl WordPerm {
    fn new(pperm: &[u16], n_p: usize) -> WordPerm {
        let combos = N_MOVES * N_SYM;
        let mut ok = vec![false; combos];
        let mut lut = vec![[[0u32; 256]; 4]; combos];
        let mut to = vec![(0u16, 0u8); combos * HALVES];
        if n_p != HALVES * 32 {
            return WordPerm { ok, lut, to };
        }
        for mk in 0..combos {
            let perm = &pperm[mk * n_p..(mk + 1) * n_p];
            // the flip of each half (the image of orientation 0) and the linear part (the
            // images of the unit vectors, from half 0)
            let lin: [u8; 5] = std::array::from_fn(|j| {
                let base = perm[0] as usize % 32;
                ((perm[1 << j] as usize % 32) ^ base) as u8
            });
            let apply_lin = |o: usize| -> usize {
                (0..5)
                    .filter(|&j| o >> j & 1 == 1)
                    .fold(0, |a, j| a ^ lin[j] as usize)
            };
            let mut good = true;
            for h in 0..HALVES {
                let base = perm[h * 32] as usize;
                let (th, t) = (base / 32, base % 32);
                to[mk * HALVES + h] = (th as u16, t as u8);
                for o in 0..32 {
                    if perm[h * 32 + o] as usize != th * 32 + (apply_lin(o) ^ t) {
                        good = false;
                    }
                }
            }
            ok[mk] = good;
            for b in 0..4 {
                for v in 0..256usize {
                    lut[mk][b][v] = (0..8)
                        .filter(|&i| v >> i & 1 == 1)
                        .fold(0u32, |a, i| a | 1 << apply_lin(8 * b + i));
                }
            }
        }
        WordPerm { ok, lut, to }
    }

    /// OR the image of block `src` (180 words) under (move, symmetry) `mk` into `out`.
    #[inline]
    pub fn apply(&self, mk: usize, src: &[u64], out: &mut [u64]) {
        let lut = &self.lut[mk];
        let to = &self.to[mk * HALVES..(mk + 1) * HALVES];
        for (h, &(th, t)) in to.iter().enumerate() {
            let w = (src[h / 2] >> (32 * (h % 2))) as u32;
            if w == 0 {
                continue;
            }
            let mut v = lut[0][(w & 255) as usize]
                | lut[1][((w >> 8) & 255) as usize]
                | lut[2][((w >> 16) & 255) as usize]
                | lut[3][(w >> 24) as usize];
            // the flip: the bits move from index o to o ^ t
            for (j, mask) in [
                0x5555_5555u32,
                0x3333_3333,
                0x0f0f_0f0f,
                0x00ff_00ff,
                0x0000_ffff,
            ]
            .into_iter()
            .enumerate()
            {
                if t >> j & 1 == 1 {
                    let sh = 1 << j;
                    v = ((v >> sh) & mask) | ((v & mask) << sh);
                }
            }
            let th = th as usize;
            out[th / 2] |= u64::from(v) << (32 * (th % 2));
        }
    }
}

/// The classes ordered by their components under the always-allowed moves (U, F, BR, BL):
/// the components in the order of their smallest class, each in class order.
fn a_move_order(edges: &[[(u16, u8); N_MOVES]], pmove: &[u8; N_MOVES]) -> Vec<u32> {
    fn find(root: &mut [usize], mut i: usize) -> usize {
        while root[i] != i {
            root[i] = root[root[i]];
            i = root[i];
        }
        i
    }
    let n = edges.len();
    let mut root: Vec<usize> = (0..n).collect();
    for c in 0..n {
        for m in 0..N_MOVES {
            let (c2, _) = edges[c][m];
            let face = usize::from(pmove[m] / 2);
            if c2 != NO_EDGE && crate::geometry::A_FACES.contains(&face) {
                let (a, b) = (find(&mut root, c), find(&mut root, c2 as usize));
                root[a.max(b)] = a.min(b);
            }
        }
    }
    let mut order: Vec<(usize, u32)> = (0..n).map(|c| (find(&mut root, c), c as u32)).collect();
    order.sort_unstable();
    order.into_iter().map(|(_, c)| c).collect()
}

/// B-centre tables of a view: states, moves (`[x][m]`), symmetries (`[k * n + x]`), and the
/// solved state.
fn b_tables(s: &Space, view: BView) -> (usize, Vec<[u32; N_MOVES]>, Vec<u32>, u32) {
    let pat = &s.pat;
    let ms = mset();
    let solved = Cube::solved().xb;
    match view {
        BView::None => (1, vec![[0; N_MOVES]], vec![0; N_SYM], 0),
        BView::Full | BView::Unlabelled => {
            // per colouring rank: its state in the view
            let (n, of_rank): (usize, Vec<u32>) = if view == BView::Full {
                (N_MSET, (0..N_MSET as u32).collect())
            } else {
                let perms = colour_perms();
                let mut ids: HashMap<u32, u32> = HashMap::new();
                let mut of_rank = vec![0u32; N_MSET];
                for r in 0..N_MSET as u32 {
                    let l = ms.unrank(r);
                    let canon = perms
                        .iter()
                        .map(|p| ms.rank(&std::array::from_fn(|i| p[l[i] as usize & 3])))
                        .min()
                        .unwrap_or(r);
                    let next = ids.len() as u32;
                    of_rank[r as usize] = *ids.entry(canon).or_insert(next);
                }
                (ids.len(), of_rank)
            };
            // a representative rank of each state
            let mut rep_rank = vec![u32::MAX; n];
            for r in 0..N_MSET as u32 {
                let st = of_rank[r as usize] as usize;
                if rep_rank[st] == u32::MAX {
                    rep_rank[st] = r;
                }
            }
            let to16 = |l: &[u8; 12]| -> [u8; 16] {
                let mut x = [0u8; 16];
                x[..12].copy_from_slice(l);
                x
            };
            let gathers: [[u8; 16]; N_MOVES] = std::array::from_fn(|m| {
                let mut g = [0x80u8; 16];
                g[..12].copy_from_slice(&pat.xb_src[m]);
                g
            });
            let xm: Vec<[u32; N_MOVES]> = rep_rank
                .iter()
                .map(|&r| {
                    let l = to16(&ms.unrank(r));
                    std::array::from_fn(|m| {
                        let after = shuffle(&l, &gathers[m]);
                        of_rank[ms.rank(&std::array::from_fn(|i| after[i])) as usize]
                    })
                })
                .collect();
            let mut xs = vec![0u32; N_SYM * n];
            for k in 0..N_SYM {
                for (st, &r) in rep_rank.iter().enumerate() {
                    let after = pat.sym_xb(k, &to16(&ms.unrank(r)));
                    xs[k * n + st] = of_rank[ms.rank(&std::array::from_fn(|i| after[i])) as usize];
                }
            }
            let x0 = of_rank[ms.rank(&solved) as usize];
            (n, xm, xs, x0)
        }
    }
}

/// The 24 permutations of the 4 colours.
fn colour_perms() -> Vec<[u8; 4]> {
    let mut out = vec![];
    for a in 0..4u8 {
        for b in 0..4u8 {
            for c in 0..4u8 {
                for d in 0..4u8 {
                    let p = [a, b, c, d];
                    let mut seen = [false; 4];
                    for &v in &p {
                        seen[v as usize] = true;
                    }
                    if seen.iter().all(|&s| s) {
                        out.push(p);
                    }
                }
            }
        }
    }
    out
}

/// Corner tables of a view: states, moves (`[p * N_MOVES + m]`), symmetries
/// (`[p * N_SYM + k]`), the solved state, and (full view) the pattern-view corner state of each.
#[allow(clippy::type_complexity)]
fn corner_tables(s: &Space, view: CornerView) -> (usize, Vec<u16>, Vec<u16>, u32, Vec<u16>) {
    let pc = &s.corners;
    match view {
        CornerView::None => (1, vec![0; N_MOVES], vec![0; N_SYM], 0, vec![]),
        CornerView::Full => {
            // states renumbered `permutation * 32 + orientations of places 0..5` (the sixth
            // follows), so that a block's 32-bit halves are the orientations of one
            // permutation ([`WordPerm`]); `p_to_pc` maps back
            let (of_state, p_to_pc) = corner_numbering(&pc.states);
            let n = pc.next.len();
            let mut next = vec![0u16; n * N_MOVES];
            let mut sym = vec![0u16; n * N_SYM];
            for (st, &id) in of_state.iter().enumerate() {
                for m in 0..N_MOVES {
                    next[id as usize * N_MOVES + m] = of_state[pc.next[st][m] as usize];
                }
                for k in 0..N_SYM {
                    sym[id as usize * N_SYM + k] = of_state[pc.sym[st][k] as usize];
                }
            }
            (
                n,
                next,
                sym,
                u32::from(of_state[pc.solved as usize]),
                p_to_pc,
            )
        }
        CornerView::Perm => {
            let mut ids: HashMap<[u8; 6], u16> = HashMap::new();
            let of_state: Vec<u16> = pc
                .states
                .iter()
                .map(|(cp, _)| {
                    let next = ids.len() as u16;
                    *ids.entry(*cp).or_insert(next)
                })
                .collect();
            let n = ids.len();
            let mut next = vec![0u16; n * N_MOVES];
            let mut sym = vec![0u16; n * N_SYM];
            for (st, &id) in of_state.iter().enumerate() {
                for m in 0..N_MOVES {
                    next[id as usize * N_MOVES + m] = of_state[pc.next[st][m] as usize];
                }
                for k in 0..N_SYM {
                    sym[id as usize * N_SYM + k] = of_state[pc.sym[st][k] as usize];
                }
            }
            (
                n,
                next,
                sym,
                u32::from(of_state[pc.solved as usize]),
                vec![],
            )
        }
    }
}

/// The pattern-view moves, for a check that the internal moves give each face turn once.
pub fn pattern_moves_covered() -> bool {
    let pm = space().pat.pmove;
    let mut seen = [false; N_PMOVES];
    for &m in &pm {
        seen[m as usize % N_PMOVES] = true;
    }
    seen.iter().all(|&s| s)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The quotient views are well defined: the moves and the symmetries map a class of the
    /// finer view into one class of the coarser one (checked through the tables' build, which
    /// takes the first member's image: every member must give the same).
    #[test]
    fn views_commute_with_moves_and_symmetries() {
        let s = space();
        let pc = &s.corners;
        // corner permutation: the image's permutation depends on the permutation only
        let mut img: HashMap<([u8; 6], usize), [u8; 6]> = HashMap::new();
        for (st, (cp, _)) in pc.states.iter().enumerate() {
            for m in 0..N_MOVES {
                let after = pc.states[pc.next[st][m] as usize].0;
                assert_eq!(*img.entry((*cp, m)).or_insert(after), after);
            }
            for k in 0..N_SYM {
                let after = pc.states[pc.sym[st][k] as usize].0;
                assert_eq!(*img.entry((*cp, 100 + k)).or_insert(after), after);
            }
        }
        assert!(pattern_moves_covered());
        // the views' sizes
        let i = Instance::new(Kind {
            corners: CornerView::Perm,
            b: BView::Unlabelled,
        });
        assert_eq!((i.n_x, i.n_p, i.n_class), (15_400, 360, 1611));
        let f = Instance::new(Kind::FULL);
        assert_eq!((f.n_x, f.n_p, f.words), (N_MSET, 11_520, 180));
        assert_eq!(f.n_positions(), 37_260 * 369_600 * 11_520);
    }

    /// The full view's numbering is a bijection, and the word-at-a-time permutation of every
    /// (move, symmetry) is the per-bit one, on random blocks of several densities.
    #[test]
    fn word_permutation_is_the_bit_permutation() {
        use rand::{Rng, SeedableRng};
        let f = Instance::new(Kind::FULL);
        let mut seen = vec![false; f.n_p];
        for &pc in &f.p_to_pc {
            assert!(!seen[pc as usize]);
            seen[pc as usize] = true;
        }
        let wp = f.word.as_ref().unwrap();
        assert!(
            wp.ok.iter().all(|&o| o),
            "some (move, symmetry) is not affine per half"
        );
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(4);
        for density in [0.01, 0.3, 0.9] {
            let mut src = vec![0u64; f.words];
            for p in 0..f.n_p {
                if rng.gen_bool(density) {
                    src[p / 64] |= 1 << (p % 64);
                }
            }
            for m in 0..N_MOVES {
                for k in 0..N_SYM {
                    let mut by_word = vec![0u64; f.words];
                    wp.apply(m * N_SYM + k, &src, &mut by_word);
                    let mut by_bit = vec![0u64; f.words];
                    let perm = f.fiber_perm(m, k);
                    for p in 0..f.n_p {
                        if src[p / 64] >> (p % 64) & 1 == 1 {
                            let t = perm[p] as usize;
                            by_bit[t / 64] |= 1 << (t % 64);
                        }
                    }
                    assert_eq!(by_word, by_bit, "move {m} symmetry {k}");
                }
            }
        }
    }
}
