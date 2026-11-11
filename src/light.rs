//! The tiny tier's A-centre/edge bound (the demo page's "Light" tables): two small tables on
//! projections of the A-centre/edge coordinate, each read through several relabellings,
//! instead of the small tier's 129 MB table.
//!
//! * **P** (the pair table, 36.3 MB): the A-centres x the positions of the edges of two
//!   A-triples (labels 0 and 1; the other six edges are not told apart), reduced by the 48
//!   symmetries that keep that split of the labels: 7864 A-centre classes x 18,480 = 145.3M
//!   entries. Its distance to the projection of `H_A` (0-11, mean 8.87) is a lower bound on
//!   the A-centre/edge distance. Read for all six pairs of triples (each through a rotation
//!   that relabels the pair as 0, 1), the largest is the bound.
//! * A bound that is the largest of several lookups cannot be decoded from its parent's (that
//!   would need each lookup's own), so the entries are not distances mod 3 but codes: P holds
//!   0 = at most 7, 1 = 8, 2 = 9, 3 = 10 or more.
//! * **S** (the split table, 0.9 MB): the A-centres x which edges belong to triples 0 or 1
//!   rather than 2 or 3 (96 symmetries), read for the three ways to split the four triples
//!   into two pairs when every P lookup says "at most 7": 0 = at most 4, 1 = 5, 2 = 6, 3 = 7
//!   or more. Without it the bound there would be the corner/B-centre table's alone.
//!
//! The two tables are one file (`ae-pair48-split96.prun`), P first. Estimated from the
//! distributions of the bounds over 4M random positions (Korf's formula, with the
//! corner/B-centre table): 1.50x the nodes of the small tier (1.72x with four P lookups, 2.0x
//! with three, 81x without S).

use crate::coord::*;
use crate::cube::*;
use crate::fast::*;
use crate::mem::HugeBuf;
use crate::prune::{saturating_codes, Graph, Table2};
use crate::sym::*;
use crate::tables::*;

/// Which labels of the edge triples a projection tells apart.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EdgeProj {
    /// triples 0 and 1 told apart, 2 and 3 merged: 18,480 labellings
    Pair,
    /// triples 0 and 1 merged, and 2 and 3: 924 labellings
    Split,
}

impl EdgeProj {
    /// Labellings of the projection.
    pub const fn count(self) -> usize {
        match self {
            EdgeProj::Pair => 18_480,
            EdgeProj::Split => 924,
        }
    }
    /// The block of each label (labels in one block are merged).
    const fn blocks(self) -> [u8; 4] {
        match self {
            EdgeProj::Pair => [0, 1, 2, 2],
            EdgeProj::Split => [0, 0, 1, 1],
        }
    }
}

/// The full A-centres x a projection of the edge triple labelling, reduced by the symmetry
/// pairs (g, t) of the 288 (see [`crate::sym`]) whose relabelling `g` keeps the projection's
/// blocks of labels. Entries: A-centre class x projected labelling.
pub struct AeProj {
    pub proj: EdgeProj,
    /// rank -> class id | g << 16 | t << 24 (the pair mapping the rank to its class's
    /// representative)
    pub cls: AClasses,
    /// symmetry pairs (other than the identity) fixing each class representative
    stab: Vec<Vec<(u8, u8)>>,
    /// colex rank of each 6-position mask (`Split`: the positions of labels 0 and 1)
    split_rank: Box<[u16]>,
    /// a multiset rank of each projected labelling (its full labelling's rank), for `decode`
    erep: Box<[u32]>,
    /// per relabelling `g`: the labels that become 0 and 1 (padded to 32 symmetries so that
    /// `& 31` replaces bounds checks)
    from01: [[u8; 2]; 32],
    /// per repositioning `t`: the gather of the positions, `new[p] = old[gather[p]]`
    gather: [[u8; 16]; 32],
    /// the fast-coordinate and label tables, held here so that the per-entry methods read no
    /// `OnceLock` (see `AeGraph` in [`crate::tables`])
    f: &'static FastTables,
    labels: &'static [[u8; 16]],
}

impl AeProj {
    pub fn new(proj: EdgeProj) -> AeProj {
        let s = syms();
        let blocks = proj.blocks();
        // the pairs whose relabelling maps blocks to blocks (the identity first, as listed)
        let pairs: Vec<(u8, u8)> = s
            .pairs
            .iter()
            .copied()
            .filter(|&(g, _)| {
                let c = &s.a_col[g as usize];
                (0..4).all(|x| {
                    (0..4).all(|y| {
                        (blocks[x] == blocks[y]) == (blocks[c[x] as usize] == blocks[c[y] as usize])
                    })
                })
            })
            .collect();
        let mut cls = AClasses::new(&pairs);
        for v in cls.sym.iter_mut() {
            let (g, t) = pairs[(*v >> 16) as usize];
            *v = (*v & 0xffff) | u32::from(g) << 16 | u32::from(t) << 24;
        }
        let ms = mset();
        let stab = cls
            .rep
            .iter()
            .map(|&r| {
                let l = ms.unrank(r);
                pairs[1..]
                    .iter()
                    .copied()
                    .filter(|&(g, t)| s.apply_a(g as usize, t as usize, &l) == l)
                    .collect()
            })
            .collect();
        let mut split_rank = vec![0u16; 4096].into_boxed_slice();
        let mut erep = vec![0u32; proj.count()].into_boxed_slice();
        match proj {
            // (the pair rank is the multiset rank without its last term: the labelling whose
            // label 2 takes the first three of the other six positions has rank r itself)
            EdgeProj::Pair => {
                for (r, e) in erep.iter_mut().enumerate() {
                    *e = r as u32;
                }
            }
            EdgeProj::Split => {
                let mut k = 0usize;
                for m in 0..4096u32 {
                    if m.count_ones() != 6 {
                        continue;
                    }
                    split_rank[m as usize] = k as u16;
                    // labels 0, 1 on the positions of m, 2, 3 on the others, three each in
                    // position order
                    let (mut inside, mut outside) = (0, 0);
                    let l: [u8; 12] = std::array::from_fn(|p| {
                        if m >> p & 1 == 1 {
                            inside += 1;
                            u8::from(inside > 3)
                        } else {
                            outside += 1;
                            2 + u8::from(outside > 3)
                        }
                    });
                    erep[k] = ms.rank(&l);
                    k += 1;
                }
            }
        }
        let from01 = std::array::from_fn(|g| {
            let c = s.a_col.get(g).copied().unwrap_or([0, 1, 2, 3]);
            let from = |y: u8| (0..4).find(|&x| c[usize::from(x)] == y).unwrap_or(y);
            [from(0), from(1)]
        });
        let gather = std::array::from_fn(|t| {
            let mut gth = [0x80u8; 16];
            if let Some(pos) = s.e_pos.get(t) {
                for p in 0..12 {
                    gth[usize::from(pos[p])] = p as u8;
                }
            }
            gth
        });
        AeProj {
            proj,
            cls,
            stab,
            split_rank,
            erep,
            from01,
            gather,
            f: fast(),
            labels: mset_labels(),
        }
    }

    /// Rank of the projected labelling with label masks `m`.
    #[inline(always)]
    fn edge_rank(&self, m: [u16; 4]) -> usize {
        match self.proj {
            EdgeProj::Pair => self.f.rank_pair(&m) as usize,
            EdgeProj::Split => self.split_rank[usize::from((m[0] | m[1]) & 0xfff)] as usize,
        }
    }

    /// Rank of the edge labelling `e` gathered by `gather` and relabelled by `g`. The
    /// projection needs only the masks of the two labels that become 0 and 1, so this is one
    /// byte shuffle and two label masks instead of relabelling every edge (the generator's
    /// depth 9 of P in one WebAssembly thread 6.28 -> 6.14 s with [`ProjGraph`], which also
    /// moves only the edges; the rest of that time is cache misses in the table).
    #[inline(always)]
    fn rank_moved(&self, e: &[u8; 16], gather: &[u8; 16], g: usize) -> usize {
        let e2 = shuffle(e, gather);
        let [l0, l1] = self.from01[g & 31];
        self.edge_rank([label_mask(&e2, l0), label_mask(&e2, l1), 0, 0])
    }

    /// The entries of the projection of `H_A` (where the distance is 0): the A-move closure
    /// of the solved position, followed in the unreduced projected coordinate, with every
    /// index of each.
    pub fn targets(&self) -> Vec<usize> {
        let f = self.f;
        let key = |y: &Fc| (f.a_rank(y), self.edge_rank(masks(&y.e)));
        let start = f.of(&Cube::solved());
        let mut seen = std::collections::HashSet::from([key(&start)]);
        let mut queue = vec![start];
        let mut k = 0;
        while let Some(&y) = queue.get(k) {
            k += 1;
            let em = masks(&y.e);
            for m in (0..N_MOVES).filter(|&m| is_a_move(m)) {
                let z = f.apply(&y, &em, m);
                if seen.insert(key(&z)) {
                    queue.push(z);
                }
            }
        }
        let mut v: Vec<usize> = queue
            .iter()
            .map(|y| self.index_fast(f, y, f.a_rank(y)))
            .collect();
        for i in 0..v.len() {
            let j = v[i];
            self.equivalents(j, &mut v);
        }
        v.sort_unstable();
        v.dedup();
        v
    }
}

/// Moves between the entries of an [`AeProj`] table, for its generation: an entry's A-centres
/// are its class representative, so the class-table entry after each move is tabulated per
/// class (as `AeGraph` in [`crate::tables`] does), and the edges need one byte shuffle per
/// move, the move's gather composed with each repositioning.
struct ProjGraph<'a, const N: usize> {
    ix: &'a AeProj,
    /// class-table entry of each class representative after each move
    next: Vec<[u32; N]>,
    /// the gather of move `m` and then repositioning `t`, at `m * 32 + t`
    gathers: Vec<[u8; 16]>,
}

impl<'a, const N: usize> ProjGraph<'a, N> {
    fn new(ix: &'a AeProj) -> Self {
        let f = ix.f;
        let d = defs();
        let next = ix
            .cls
            .rep
            .iter()
            .map(|&r| {
                let x = Fc {
                    a: ix.labels[r as usize],
                    ..Fc::ZERO
                };
                std::array::from_fn(|m| *ix.class_entry(f.a_rank(&f.apply_any::<N>(&x, 0, 0, m))))
            })
            .collect();
        let mut gathers = vec![[0u8; 16]; N * 32];
        for m in 0..N {
            let mut mv = [0x80u8; 16];
            mv[..12].copy_from_slice(&d.moves[m].ep);
            for t in 0..32 {
                gathers[m * 32 + t] = shuffle(&mv, &ix.gather[t]);
            }
        }
        ProjGraph { ix, next, gathers }
    }
}

impl<const N: usize> Graph for ProjGraph<'_, N> {
    const MOVES: usize = N;
    /// the edge labels and the class
    type Node = ([u8; 16], usize);
    type Cache = ();
    #[inline(always)]
    fn node(&self, i: usize, (): &mut ()) -> ([u8; 16], usize) {
        let cs = self.ix.proj.count();
        (self.ix.labels[self.ix.erep[i % cs] as usize], i / cs)
    }
    #[inline(always)]
    fn step(&self, (e, c): &([u8; 16], usize), m: usize) -> usize {
        let v = self.next.get(*c).map_or(0, |n| n[m % N]);
        let (g, t) = ((v >> 16 & 0xff) as usize, (v >> 24) as usize & 31);
        let gather = self.gathers.get((m % N) * 32 + t).unwrap_or(&[0; 16]);
        (v & 0xffff) as usize * self.ix.proj.count() + self.ix.rank_moved(e, gather, g)
    }
    fn equivalents(&self, j: usize, out: &mut Vec<usize>) {
        self.ix.equivalents(j, out);
    }
}

impl AeProj {
    /// The table of saturating codes ([`saturating_codes`]: 0 = at most `base`, 1 and 2 the
    /// next two distances, 3 = more) in the metric of `moves`.
    pub fn codes(&self, base: u8, moves: MoveSet, verbose: bool) -> Table2 {
        let name = format!("{}{}", self.name(), moves.table_suffix());
        let targets = self.targets();
        if moves.slice {
            let g = ProjGraph::<MAX_MOVES>::new(self);
            saturating_codes(&g, self.size(), &targets, base, &name, verbose)
        } else {
            let g = ProjGraph::<N_MOVES>::new(self);
            saturating_codes(&g, self.size(), &targets, base, &name, verbose)
        }
    }
}

impl AeIndex for AeProj {
    fn size(&self) -> usize {
        self.cls.rep.len() * self.proj.count()
    }
    fn index(&self, x: &Coords) -> usize {
        self.index_fast(self.f, &self.f.of_coords(x), x.a)
    }
    fn decode(&self, i: usize) -> (u32, u32, u8) {
        let cs = self.proj.count();
        let c = (i / cs).min(self.cls.rep.len().saturating_sub(1));
        (self.cls.rep[c], self.erep[i % cs], 0)
    }
    fn name(&self) -> &'static str {
        match self.proj {
            EdgeProj::Pair => "ae-pair48",
            EdgeProj::Split => "ae-split96",
        }
    }
    fn equivalents(&self, i: usize, out: &mut Vec<usize>) {
        let cs = self.proj.count();
        let c = i / cs;
        let Some(stab) = self.stab.get(c) else {
            return;
        };
        if stab.is_empty() {
            return;
        }
        let e = &self.labels[self.erep[i % cs] as usize];
        for &(g, t) in stab {
            out.push(c * cs + self.rank_moved(e, &self.gather[usize::from(t) & 31], g.into()));
        }
    }
    #[inline(always)]
    fn index_entry(&self, _f: &FastTables, x: &Fc, v: u32) -> usize {
        let (g, t) = ((v >> 16 & 0xff) as usize, (v >> 24) as usize);
        (v & 0xffff) as usize * self.proj.count() + self.rank_moved(&x.e, &self.gather[t & 31], g)
    }
    fn class_size(&self) -> usize {
        self.proj.count()
    }
    fn entry_parts(&self, _f: &FastTables, _x: &Fc, _v: u32) -> Option<(usize, usize, u8)> {
        None
    }
    fn class_entry(&self, a_rank: u32) -> &u32 {
        &self.cls.sym[a_rank as usize]
    }
    #[inline(always)]
    fn prefetch_class(&self, a_rank: u32) {
        prefetch_at(&self.cls.sym, a_rank as usize);
    }
}

/// Codes of P: at most `P_BASE`, then the next two distances, then more.
const P_BASE: u8 = 7;
/// Codes of S.
const S_BASE: u8 = 4;
/// The bound of each P code (code 0: S decides).
const P_VALUE: [u8; 4] = [0, P_BASE + 1, P_BASE + 2, P_BASE + 3];
/// The bound of each S code.
const S_VALUE: [u8; 4] = [0, S_BASE + 1, S_BASE + 2, S_BASE + 3];

/// P lookups per position (the six pairs of triples; with four, estimated 1.72x the nodes of
/// the small tier instead of 1.50x, and measured on one thread: first solutions of random
/// positions in 1.02-1.06x the time, 6 instead of 8 of 10 within 10 s).
const P_LOOKUPS: usize = 6;
/// S lookups (the three splits).
const S_LOOKUPS: usize = 3;

/// Pending keys ([`AeLight::pending`]): a final bound (flag `DONE`), or the index of the first
/// P lookup (bits `INDEX`) with the S bound above it (`S_SHIFT`). (P has fewer than 2^28
/// entries, so a key fits 32 bits.)
const DONE: usize = 1 << 31;
const S_SHIFT: u32 = 28;
const INDEX: usize = (1 << S_SHIFT) - 1;

/// The 2-bit entry `i` of a table's words (0 outside the table).
#[inline(always)]
fn code(tab: &[u64], i: usize) -> u8 {
    tab.get(i >> 5)
        .map_or(0, |w| (w >> ((i & 31) << 1) & 3) as u8)
}

/// The position `x` with its A-centre colours and edge triples relabelled by the rotation `r`
/// (as far from `H_A`: the pair (r, identity) is one of the 288 symmetries).
fn relabelled(s: &Syms, x: &Fc, r: usize) -> Fc {
    let mut y = *x;
    for p in 0..12 {
        y.a[p] = s.a_col[r][usize::from(x.a[p] & 3)];
        y.e[p] = s.a_col[r][usize::from(x.e[p] & 3)];
    }
    y
}

/// The tiny tier's A-centre/edge bound (see the module documentation). As an [`AeIndex`] it is
/// the pair table's index (the first part of its table); the search reads it through
/// [`AeIndex::keys_batch_pf`] and [`AeIndex::keys_sel_pf`], whose keys, once
/// [`AeIndex::resolve`] has made the rest of the lookups, are the bounds themselves.
///
/// The lookups are the positions relabelled by a rotation `r_j` for each pair of triples
/// (taking it to triples 0, 1), and for each split. They need no class lookup each: if `σ` maps
/// the A-centres to their P class representative `R`, then `r_j σ⁻¹ = h c_k` with `h` in P's
/// symmetry group `H` and `c_k = r_k` (the `r_j σ⁻¹` lie in the six different cosets `H r_k` as
/// `j` varies), so the lookup of `r_j` is that of `c_k` applied to `R`, whose P class and
/// symmetry `θ` (to its own representative, after `c_k`) are tabulated per class and coset:
/// the entry is that class and the edges moved by `θ σ`. The same with S's group and the
/// three splits.
pub struct AeLight {
    pub pair: AeProj,
    pub split: AeProj,
    /// per P class and pair coset: the P class of `c_k` applied to the representative, and
    /// `θ` (class | g << 16 | t << 24)
    tp: Box<[[u32; P_LOOKUPS]]>,
    /// per P class and split coset: the same for S
    ts: Box<[[u32; S_LOOKUPS]]>,
    /// `mul[i][j]`: the symmetry `i` after `j` ([`Syms::mul`]; padded to 32 so that `& 31`
    /// replaces bounds checks)
    mul: [[u8; 32]; 32],
}

impl AeLight {
    pub fn new() -> AeLight {
        let s = syms();
        let rotations: Vec<usize> = (0..N_GEO)
            .filter(|&g| crate::geometry::det(&s.mats[g]) == 1)
            .collect();
        // a rotation relabelling triples i, j as 0, 1 (the rotations act on the four
        // A-faces as the alternating group, which is 2-transitive), and one relabelling
        // {0, x} as {0, 1}
        let find = |ok: &dyn Fn(&[u8; 4]) -> bool| {
            rotations
                .iter()
                .copied()
                .find(|&g| ok(&s.a_col[g]))
                .unwrap_or(0)
        };
        let ij = [(0, 1), (2, 3), (0, 2), (1, 3), (0, 3), (1, 2)];
        let pairs = ij.map(|(i, j)| find(&|c| c[i] == 0 && c[j] == 1));
        let splits = [1, 2, 3].map(|x| find(&|c| c[0] <= 1 && c[x] <= 1));
        let pair = AeProj::new(EdgeProj::Pair);
        let split = AeProj::new(EdgeProj::Split);
        let ms = mset();
        // the entry of `c_k` = (r, identity) applied to the representative `l`, in `table`
        let entry = |table: &AeProj, l: &[u8; 12], r: usize| -> u32 {
            let v = *table.class_entry(ms.rank(&s.apply_a(r, 0, l)));
            let (g, t) = ((v >> 16 & 0xff) as usize, (v >> 24) as usize);
            (v & 0xffff) | u32::from(s.mul[g][r]) << 16 | (t as u32) << 24
        };
        let reps: Vec<[u8; 12]> = pair.cls.rep.iter().map(|&r| ms.unrank(r)).collect();
        let tp = reps
            .iter()
            .map(|l| pairs.map(|r| entry(&pair, l, r)))
            .collect();
        let ts = reps
            .iter()
            .map(|l| splits.map(|r| entry(&split, l, r)))
            .collect();
        AeLight {
            pair,
            split,
            tp,
            ts,
            mul: std::array::from_fn(|i| {
                std::array::from_fn(|j| {
                    if i < N_GEO && j < N_GEO {
                        s.mul[i][j]
                    } else {
                        0
                    }
                })
            }),
        }
    }

    /// The entry `w` of [`AeLight::tp`] or [`AeLight::ts`] with the symmetry `σ` of the
    /// position's class entry `v` composed in (`θ σ`).
    #[inline(always)]
    fn compose(&self, w: u32, v: u32) -> u32 {
        let (g, t) = ((v >> 16 & 0xff) as usize, (v >> 24) as usize);
        let (tg, tt) = ((w >> 16 & 0xff) as usize, (w >> 24) as usize);
        (w & 0xffff)
            | u32::from(self.mul[tg & 31][g & 31]) << 16
            | u32::from(self.mul[tt & 31][t & 31]) << 24
    }

    /// The S bound of position `x` with P class entry `v` (S is 0.9 MB, so its words are
    /// usually cached and read at once).
    #[inline(always)]
    fn split_value(&self, f: &FastTables, x: &Fc, v: u32, tab: &[u64]) -> u8 {
        let base = self.pair.size();
        let row = self
            .ts
            .get((v & 0xffff) as usize)
            .unwrap_or(&[0; S_LOOKUPS]);
        let c = row
            .iter()
            .map(|&w| code(tab, base + self.split.index_entry(f, x, self.compose(w, v))))
            .max()
            .unwrap_or(0);
        S_VALUE[usize::from(c & 3)]
    }

    /// The pending keys of the positions `xs[i]` for `i` in `ids`, written to `out[i]`, for a
    /// search that compares their bounds with `cut` (`usize::MAX`: in full): where S can
    /// exceed `cut` it is read first, and a position it puts above `cut` gets its final key
    /// (`DONE` | bound); the others get the index of their first P lookup, whose table word is
    /// prefetched, with their S bound (`<< S_SHIFT`; 0 where S was not read).
    /// [`AeLight::resolve_keys`] makes them bounds. (The class entries of all positions are
    /// prefetched first.)
    #[inline(always)]
    fn pending<const N: usize>(
        &self,
        f: &FastTables,
        xs: &[Fc],
        ids: &[u8],
        out: &mut [usize],
        tab: &[u64],
        cut: usize,
    ) {
        let n = ids.len().min(N);
        let ids = &ids[..n];
        let x = |j: usize| xs.get(usize::from(ids[j])).unwrap_or(&Fc::ZERO);
        let mut ar = [0u32; N];
        for j in 0..n {
            ar[j] = f.a_rank(x(j));
            self.pair.prefetch_class(ar[j]);
        }
        let s_can = cut < usize::from(S_VALUE[3]) || cut == usize::MAX;
        // Below P's smallest bound P is not read: S decides (in time-limited searches P cut
        // 0.2% of the children S kept at cut 6 and none below; 1T nodes/s 30.1 -> 31.9 M
        // without it, when added).
        let p_can = cut >= usize::from(P_VALUE[1]) - 1;
        for j in 0..n {
            let v = *self.pair.class_entry(ar[j]);
            let sv = if s_can {
                self.split_value(f, x(j), v, tab)
            } else {
                0
            };
            let key = if usize::from(sv) > cut || !p_can {
                DONE | usize::from(sv)
            } else {
                let w = self.tp.get((v & 0xffff) as usize).map_or(0, |r| r[0]);
                let i = self.pair.index_entry(f, x(j), self.compose(w, v));
                prefetch_at(tab, i >> 5);
                i | usize::from(sv) << S_SHIFT
            };
            if let Some(o) = out.get_mut(usize::from(ids[j])) {
                *o = key;
            }
        }
    }

    /// Turn the pending keys `keys[i]` of the positions `xs[i]`, `i` in `ids` (made by
    /// [`AeLight::pending`] with the same `cut`), into their bounds. The first P lookup's word
    /// was prefetched; the other P lookups are made (in one round, their cache misses in
    /// flight together) only for the positions whose bound so far is at most `cut`, if a P
    /// lookup can still exceed it.
    #[inline(always)]
    fn resolve_keys<const N: usize>(
        &self,
        f: &FastTables,
        xs: &[Fc],
        ids: &[u8],
        keys: &mut [usize],
        tab: &[u64],
        cut: usize,
    ) {
        let n = ids.len().min(N);
        let ids = &ids[..n];
        let full = cut == usize::MAX;
        let x = |j: usize| xs.get(usize::from(ids[j])).unwrap_or(&Fc::ZERO);
        // (a P lookup can exceed `cut` only if it is below the largest P value)
        let p_can = full || cut < usize::from(P_VALUE[3]);
        let mut best = [0u8; N];
        let mut sv = [0u8; N];
        let mut done = [false; N];
        let mut more = [0u8; N];
        let mut m = 0;
        for j in 0..n {
            let key = keys.get(usize::from(ids[j])).copied().unwrap_or(DONE);
            done[j] = key & DONE != 0;
            sv[j] = (key >> S_SHIFT & 7) as u8;
            best[j] = if done[j] { 0 } else { code(tab, key & INDEX) };
            more[m % N] = j as u8;
            m += usize::from(
                !done[j] && p_can && usize::from(P_VALUE[usize::from(best[j] & 3)]) <= cut,
            );
        }
        if m > 0 {
            let mut cv = [0u32; N];
            for &j in &more[..m] {
                let j = usize::from(j) % N;
                cv[j] = f.a_rank(x(j));
                self.pair.prefetch_class(cv[j]);
            }
            let mut pi = [[0usize; P_LOOKUPS]; N];
            for &j in &more[..m] {
                let j = usize::from(j) % N;
                cv[j] = *self.pair.class_entry(cv[j]);
                let w = self
                    .tp
                    .get((cv[j] & 0xffff) as usize)
                    .unwrap_or(&[0; P_LOOKUPS]);
                for k in 1..P_LOOKUPS {
                    pi[j][k] = self.pair.index_entry(f, x(j), self.compose(w[k], cv[j]));
                    prefetch_at(tab, pi[j][k] >> 5);
                }
            }
            for &j in &more[..m] {
                let j = usize::from(j) % N;
                let c = pi[j][1..].iter().map(|&i| code(tab, i)).max().unwrap_or(0);
                best[j] = best[j].max(c);
            }
        }
        for j in 0..n {
            if let Some(o) = keys.get_mut(usize::from(ids[j])) {
                if done[j] {
                    *o &= INDEX;
                } else {
                    *o = usize::from(if best[j] > 0 {
                        P_VALUE[usize::from(best[j] & 3)]
                    } else {
                        sv[j]
                    });
                }
            }
        }
    }

    /// The bound of one position.
    pub fn bound_of(&self, f: &FastTables, x: &Fc, tab: &[u64]) -> u8 {
        let mut keys = [0usize; 1];
        let xs = std::slice::from_ref(x);
        self.pending::<1>(f, xs, &[0], &mut keys, tab, usize::MAX);
        self.resolve_keys::<1>(f, xs, &[0], &mut keys, tab, usize::MAX);
        keys[0] as u8
    }

    /// The codes of position `x` in its P lookups and its S lookups, the slow way: each
    /// relabelled position looked up in full (for checks).
    pub fn codes(&self, f: &FastTables, x: &Fc, tab: &[u64]) -> ([u8; P_LOOKUPS], [u8; S_LOOKUPS]) {
        let s = syms();
        let rotations: Vec<usize> = (0..N_GEO)
            .filter(|&g| crate::geometry::det(&s.mats[g]) == 1)
            .collect();
        let find = |ok: &dyn Fn(&[u8; 4]) -> bool| {
            rotations
                .iter()
                .copied()
                .find(|&g| ok(&s.a_col[g]))
                .unwrap_or(0)
        };
        let ij = [(0, 1), (2, 3), (0, 2), (1, 3), (0, 3), (1, 2)];
        let p = ij.map(|(i, j)| {
            let y = relabelled(s, x, find(&|c| c[i] == 0 && c[j] == 1));
            code(tab, self.pair.index_fast(f, &y, f.a_rank(&y)))
        });
        let base = self.pair.size();
        let sp = [1, 2, 3].map(|k| {
            let y = relabelled(s, x, find(&|c| c[0] <= 1 && c[k] <= 1));
            code(tab, base + self.split.index_fast(f, &y, f.a_rank(&y)))
        });
        (p, sp)
    }
}

impl AeIndex for AeLight {
    const DIRECT: bool = true;
    fn size(&self) -> usize {
        self.pair.size() + self.split.size()
    }
    fn index(&self, x: &Coords) -> usize {
        self.pair.index(x)
    }
    fn decode(&self, i: usize) -> (u32, u32, u8) {
        self.pair.decode(i)
    }
    fn name(&self) -> &'static str {
        "ae-pair48-split96"
    }
    fn equivalents(&self, i: usize, out: &mut Vec<usize>) {
        self.pair.equivalents(i, out);
    }
    fn index_entry(&self, f: &FastTables, x: &Fc, v: u32) -> usize {
        self.pair.index_entry(f, x, v)
    }
    fn class_size(&self) -> usize {
        self.pair.class_size()
    }
    fn entry_parts(&self, _f: &FastTables, _x: &Fc, _v: u32) -> Option<(usize, usize, u8)> {
        None
    }
    fn class_entry(&self, a_rank: u32) -> &u32 {
        self.pair.class_entry(a_rank)
    }
    fn prefetch_class(&self, a_rank: u32) {
        self.pair.prefetch_class(a_rank);
    }
    #[inline(always)]
    fn keys_batch_pf<const N: usize>(
        &self,
        f: &FastTables,
        xs: &[Fc],
        out: &mut [usize],
        tab: &[u64],
        cut: usize,
    ) {
        let n = xs.len().min(ALL_IDS.len());
        self.pending::<N>(f, xs, &ALL_IDS[..n], out, tab, cut);
    }
    #[inline(always)]
    fn keys_sel_pf<const N: usize>(
        &self,
        f: &FastTables,
        xs: &[Fc; N],
        sel: &[u8],
        out: &mut [usize; N],
        tab: &[u64],
        cut: usize,
    ) {
        self.pending::<N>(f, xs, sel, out, tab, cut);
    }
    #[inline(always)]
    fn resolve<const N: usize>(
        &self,
        f: &FastTables,
        xs: &[Fc; N],
        sel: &[u8],
        count: usize,
        keys: &mut [usize; N],
        tab: &[u64],
        cut: usize,
    ) {
        let sel = sel.get(..count).unwrap_or(sel);
        self.resolve_keys::<N>(f, xs, sel, keys, tab, cut);
    }
    #[inline(always)]
    fn value(&self, _ae: &Table2, _parent: u8, key: usize) -> u8 {
        key as u8
    }
    fn bound(&self, f: &FastTables, x: &Fc, ae: &Table2) -> Option<u8> {
        Some(self.bound_of(f, x, &ae.data))
    }
    fn build(&self, moves: MoveSet, verbose: bool) -> Table2 {
        let p = self.pair.codes(P_BASE, moves, verbose);
        let s = self.split.codes(S_BASE, moves, verbose);
        // (P has a whole number of words: 7864 classes x 18,480 entries, a multiple of 32)
        let mut data = HugeBuf::<u64>::for_overwrite(p.data.len() + s.data.len());
        data[..p.data.len()].copy_from_slice(&p.data);
        data[p.data.len()..].copy_from_slice(&s.data);
        if verbose {
            eprintln!(
                "{}: P codes (at most {P_BASE}, then single distances, then more) {:?}, S {:?}",
                self.name(),
                p.hist,
                s.hist
            );
        }
        Table2 {
            data,
            n: self.size(),
            hist: p.hist,
        }
    }
    /// The codes of every lookup are those of distances: at a position and its neighbours
    /// they differ by at most one, a code for a single distance (1 or 2) has a neighbour one
    /// lower, and the solved position has code 0 everywhere. Checked at a few random positions
    /// and short scrambles (a table of another layout fails at once).
    fn check_table(&self, ae: &Table2, moves: MoveSet) -> bool {
        use rand::{Rng, SeedableRng};
        let f = self.pair.f;
        let tab = &ae.data[..];
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(0x1167);
        let n_moves = moves.n_moves();
        let ok = |x: &Fc| {
            let (p, s) = self.codes(f, x, tab);
            let em = masks(&x.e);
            let ns: Vec<([u8; P_LOOKUPS], [u8; S_LOOKUPS])> = (0..n_moves)
                .map(|m| self.codes(f, &f.apply_m(x, &em, m), tab))
                .collect();
            let fine = |c: u8, nb: &mut dyn Iterator<Item = u8>| {
                let nb: Vec<u8> = nb.collect();
                nb.iter().all(|&d| c.abs_diff(d) <= 1)
                    && (!(1..=2).contains(&c) || nb.contains(&(c - 1)))
            };
            (0..P_LOOKUPS).all(|l| fine(p[l], &mut ns.iter().map(|n| n.0[l])))
                && (0..S_LOOKUPS).all(|l| fine(s[l], &mut ns.iter().map(|n| n.1[l])))
        };
        if self.codes(f, &f.of(&Cube::solved()), tab) != ([0; P_LOOKUPS], [0; S_LOOKUPS]) {
            return false;
        }
        (0..24).all(|i| {
            let mut c = Cube::solved();
            // random positions, and scrambles of 4-9 moves (low codes)
            if i < 8 {
                c = Cube::random(&mut rng);
            } else {
                for _ in 0..(i % 6 + 4) {
                    c = c.apply(rng.gen_range(0..n_moves));
                }
            }
            ok(&f.of(&c))
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{Rng, SeedableRng};

    /// The generator's neighbours are the indices of the moved positions, for both projections
    /// and both move sets.
    fn check_graph<const N: usize>(ix: &AeProj) {
        let g = ProjGraph::<N>::new(ix);
        let f = fast();
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(45);
        for _ in 0..3000 {
            let i = rng.gen_range(0..ix.size());
            let (a, e, _) = ix.decode(i);
            let x = f.of_coords(&Coords {
                c: 0,
                a,
                b: 0,
                e,
                ebits: 0,
            });
            // (the decoded position has index i, or an equivalent one)
            let mut eq = vec![ix.index_fast(f, &x, f.a_rank(&x))];
            ix.equivalents(eq[0], &mut eq);
            assert!(eq.contains(&i), "{i}");
            let node = g.node(i, &mut ());
            let em = masks(&x.e);
            for m in 0..N {
                let y = f.apply_m(&x, &em, m);
                let j = g.step(&node, m);
                let mut eq = vec![ix.index_fast(f, &y, f.a_rank(&y))];
                ix.equivalents(eq[0], &mut eq);
                assert!(eq.contains(&j), "{i} {m}");
            }
        }
    }

    #[test]
    fn graph_steps_match_index() {
        for proj in [EdgeProj::Pair, EdgeProj::Split] {
            let ix = AeProj::new(proj);
            check_graph::<N_MOVES>(&ix);
            check_graph::<MAX_MOVES>(&ix);
        }
    }
}
