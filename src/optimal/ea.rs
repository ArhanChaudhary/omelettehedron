//! The edge x A-centre table of the optimal solver: the distance to solved of the whole edge
//! permutation and the A-centre colouring, up to whole-puzzle rotation and up to cycling the
//! edges of every A-face together.
//!
//! The projection keeps no corners, so it lives up to the 12 whole-puzzle rotations (every
//! rotation of the solved state counts as solved), and the 288 pairs of the two-phase
//! A-centre/edge table ([`crate::sym`]: relabel by a class-preserving symmetry `g`, reposition
//! by `t`, both proper or both improper) preserve the distance: `t g^-1` is a rotation, and a
//! pair is that rotation times the conjugation by `g`, which fixes the solved state and maps
//! face turns to face turns. So the table needs 1/288 of the projection's states, against 1/48
//! for a table with the corners (corners x edges) and 1/6 for the corner x B-centre x
//! edge-triple table ([`super::cbe`]).
//!
//! The full projection would be 1,396 A-colouring classes x 239,500,800 edge permutations =
//! 3.34e11 entries (83.6 GB). It is quotiented by `C`, which relabels every edge piece as the
//! next piece of its A-triple (all four triples turned one step together): a relabelling
//! commutes with the moves, a proper symmetry commutes with `C` and an improper one inverts it,
//! so the quotient is a projection with the same symmetries. `C` of the solved state is four
//! A-face turns away (an A-face turn moves no A-centre off its face), so the quotient merges
//! the solved state with two states 4 moves away: 1.11e11 entries, 27.9 GB.
//!
//! Entry `class * N_EQ + qrank`: `class` is the symmetry class of the A-colouring (as in the
//! two-phase table), and `qrank` the rank of the edge permutation after the pair that maps the
//! colouring to its class representative, in the normal form of its `C`-orbit (the piece at
//! position 0 is the first of its triple): its triple (4) x the Lehmer code of positions 1..10
//! (11!/2).
//!
//! Two lookups per position: the A-centres with the edges, and through `SWAP_AB` the B-centres
//! with the edges (cycled by B-triple).

use super::coords::{eperm_unrank, opt_moves, OptMoves, Pos, N_EPERM};
use super::sym::sym48;
use crate::coord::{triple_pieces, N_MSET};
use crate::fast::{fast, masks, mset_labels, shuffle, FastTables};
use crate::geometry::det;
use crate::prune::{bfs_until, Graph, Table2};
use crate::sym::{syms, AClasses};
use crate::tables::prefetch_at;
use std::sync::atomic::AtomicBool;

/// Edge permutations up to `C`: `N_EPERM` / 3.
pub const N_EQ: usize = N_EPERM / 3;
/// Lehmer codes of positions 1..10 of an even permutation of 12: 11!/2.
const N_REST: usize = N_EPERM / 12;

/// Lehmer code of positions 1..10 of an even permutation (the digits of [`eperm_rank`] after
/// the first), `0..N_REST`.
///
/// [`eperm_rank`]: super::coords::eperm_rank
#[inline(always)]
fn rest_rank(e: &[u8; 16]) -> u32 {
    let mut seen = 1u32 << (e[0] & 15);
    let mut r = 0u32;
    for i in 1..10 {
        let p = u32::from(e[i] & 15);
        // (entries below 12: the shift stays in range)
        let c = p - (seen & ((1 << p) - 1)).count_ones();
        r = r * (12 - i as u32) + c;
        seen |= 1 << p;
    }
    r
}

/// Index of the edge x A-centre table.
pub struct EaIndex {
    /// A-colouring rank -> class id (low 16 bits) | g << 16 | t << 24
    sym: Vec<u32>,
    /// class id -> representative A-colouring rank
    pub rep: Vec<u32>,
    /// per class: the pairs (g, t) other than the identity that fix the representative
    stab: Vec<Vec<(u8, u8)>>,
    /// per g (padded to 32): edge piece relabelling
    relab: Box<[[u8; 16]; 32]>,
    /// per t (padded to 32): gather of the edge positions
    gath: Box<[[u8; 16]; 32]>,
    /// `C^j` as a piece relabelling, `j < 3`
    cpow: [[u8; 16]; 4],
    /// per piece: the `j` such that `C^j` of the piece is the first of its triple
    to_first: [u8; 16],
    /// per piece: its A-triple, if it is the first of its triple (else 0)
    first_tri: [u8; 16],
    /// the first piece of each triple
    firsts: [u8; 4],
    /// the `SWAP_AB` conjugate: its edges (relabel, gather) and its A-centres from the
    /// B-centres (relabel, gather)
    w_e_relab: [u8; 16],
    w_e_gath: [u8; 16],
    w_a_relab: [u8; 16],
    w_a_gath: [u8; 16],
    f: &'static FastTables,
}

/// The cyclic order of each A-triple that the proper symmetries keep and the improper ones
/// reverse (`None` if there is none, which the geometry rules out). Returns, per triple, its
/// pieces in that order.
fn cyclic_triples() -> Option<[[u8; 3]; 4]> {
    let s = syms();
    let tp = triple_pieces();
    (0..16u32).find_map(|mask| {
        let order: [[u8; 3]; 4] = std::array::from_fn(|t| {
            let [a, b, c] = tp[t];
            if mask >> t & 1 == 0 {
                [a, b, c]
            } else {
                [a, c, b]
            }
        });
        let mut next = [0u8; 12];
        for o in &order {
            for i in 0..3 {
                next[usize::from(o[i])] = o[(i + 1) % 3];
            }
        }
        let ok = (0..s.mats.len()).all(|g| {
            let e = &s.e_piece[g];
            // proper: e(next(x)) == next(e(x)); improper: e(next(x)) == prev(e(x)), i.e.
            // next(e(next(x))) == e(x)
            (0..12).all(|x| {
                let en = e[usize::from(next[x])];
                if det(&s.mats[g]) > 0 {
                    en == next[usize::from(e[x])]
                } else {
                    next[usize::from(en)] == e[x]
                }
            })
        });
        ok.then_some(order)
    })
}

impl EaIndex {
    /// # Errors
    ///
    /// If no cyclic order of the triples is kept by the symmetries (not for this geometry).
    pub fn new() -> std::io::Result<EaIndex> {
        let s = syms();
        let order = cyclic_triples()
            .ok_or_else(|| std::io::Error::other("no symmetric cyclic order of the A-triples"))?;
        let AClasses { sym: raw, rep } = AClasses::new(&s.pairs);
        let sym: Vec<u32> = raw
            .iter()
            .map(|&v| {
                let (g, t) = s.pairs[(v >> 16) as usize];
                (v & 0xffff) | u32::from(g) << 16 | u32::from(t) << 24
            })
            .collect();
        drop(raw);
        let ms = crate::coord::mset();
        let stab = rep
            .iter()
            .map(|&r| {
                let l = ms.unrank(r);
                s.pairs[1..]
                    .iter()
                    .copied()
                    .filter(|&(g, t)| s.apply_a(g as usize, t as usize, &l) == l)
                    .collect()
            })
            .collect();
        let mut relab = Box::new([[0u8; 16]; 32]);
        let mut gath = Box::new([[0x80u8; 16]; 32]);
        for g in 0..s.mats.len() {
            relab[g][..12].copy_from_slice(&s.e_piece[g]);
            for p in 0..12 {
                gath[g][usize::from(s.e_pos[g][p])] = p as u8;
            }
        }
        let mut cpow = [[0u8; 16]; 4];
        let mut to_first = [0u8; 16];
        let mut first_tri = [0u8; 16];
        let mut firsts = [0u8; 4];
        for (t, o) in order.iter().enumerate() {
            firsts[t] = o[0];
            first_tri[usize::from(o[0])] = t as u8;
            for i in 0..3 {
                let x = usize::from(o[i]);
                for j in 0..3 {
                    cpow[j][x] = o[(i + j) % 3];
                }
                to_first[x] = ((3 - i) % 3) as u8;
            }
        }
        let w48 = sym48();
        let w = w48.swap;
        let mut w_e_relab = [0u8; 16];
        let mut w_e_gath = [0x80u8; 16];
        let mut w_a_relab = [0u8; 16];
        let mut w_a_gath = [0x80u8; 16];
        w_e_relab[..12].copy_from_slice(&w48.epos[w]);
        w_a_relab[..4].copy_from_slice(&w48.ccol[w][1]);
        // (SWAP_AB fixes the front corner: positions move by it too)
        for p in 0..12 {
            w_e_gath[usize::from(w48.epos[w][p])] = p as u8;
            w_a_gath[usize::from(w48.cpos[w][1][p])] = p as u8;
        }
        Ok(EaIndex {
            sym,
            rep,
            stab,
            relab,
            gath,
            cpow,
            to_first,
            first_tri,
            firsts,
            w_e_relab,
            w_e_gath,
            w_a_relab,
            w_a_gath,
            f: fast(),
        })
    }

    /// Entries: A-colouring classes x edge permutations up to `C`.
    /// (More than 2^32: 64-bit targets only; saturates elsewhere.)
    pub fn size(&self) -> usize {
        self.rep.len().saturating_mul(N_EQ)
    }

    /// `qrank` of edge permutation `ep` after pair (g, t).
    #[inline(always)]
    fn qrank(&self, ep: &[u8; 16], g: usize, t: usize) -> usize {
        let e = shuffle(&shuffle(&self.relab[g & 31], ep), &self.gath[t & 31]);
        let j = usize::from(self.to_first[usize::from(e[0] & 15)] & 3);
        let e = shuffle(&self.cpow[j], &e);
        usize::from(self.first_tri[usize::from(e[0] & 15)]) * N_REST + rest_rank(&e) as usize
    }

    /// Index of the A-colouring `xa` (lanes 0..12) with the edges `ep`.
    #[inline(always)]
    pub fn index_raw(&self, xa: &[u8; 16], ep: &[u8; 16]) -> usize {
        self.index_ranked(self.f.rank(&masks(xa)) as usize, ep)
    }

    /// Index of the A-colouring of rank `r` with the edges `ep`.
    #[inline(always)]
    fn index_ranked(&self, r: usize, ep: &[u8; 16]) -> usize {
        let v = self.sym.get(r).copied().unwrap_or(0);
        let (g, t) = ((v >> 16 & 0xff) as usize, (v >> 24) as usize);
        (v & 0xffff) as usize * N_EQ + self.qrank(ep, g, t)
    }

    /// The A-centres and edges of lookup `k` of `p`: 0 the position's, 1 those of its
    /// `SWAP_AB` conjugate (the B-centres as A-centres, and the edges cycled by B-triple).
    #[inline(always)]
    fn view(&self, k: usize, p: &Pos) -> ([u8; 16], [u8; 16]) {
        if k == 0 {
            (p.xa, p.ep)
        } else {
            (
                shuffle(&shuffle(&self.w_a_relab, &p.xb), &self.w_a_gath),
                shuffle(&shuffle(&self.w_e_relab, &p.ep), &self.w_e_gath),
            )
        }
    }

    /// Index of lookup `k` of position `p` (see [`EaIndex::view`]).
    #[inline(always)]
    pub fn index(&self, k: usize, p: &Pos) -> usize {
        let (xa, ep) = self.view(k, p);
        self.index_raw(&xa, &ep)
    }

    /// The A-colouring rank of lookup `k` of `p`, its class entry prefetched (the first step of
    /// the index, as `Bounds::pre`).
    #[inline(always)]
    pub fn pre(&self, k: usize, p: &Pos) -> u32 {
        let xa = if k == 0 {
            p.xa
        } else {
            shuffle(&shuffle(&self.w_a_relab, &p.xb), &self.w_a_gath)
        };
        let r = self.f.rank(&masks(&xa));
        prefetch_at(&self.sym, r as usize);
        r
    }

    /// [`EaIndex::index`] given the A-colouring rank from [`EaIndex::pre`].
    #[inline(always)]
    pub fn index_with(&self, k: usize, p: &Pos, a_rank: u32) -> usize {
        let ep = if k == 0 {
            p.ep
        } else {
            shuffle(&shuffle(&self.w_e_relab, &p.ep), &self.w_e_gath)
        };
        self.index_ranked(a_rank as usize, &ep)
    }

    /// A projection (A-colouring, edge permutation) of index `i`.
    pub fn decode(&self, i: usize) -> ([u8; 16], [u8; 16]) {
        let cl = i / N_EQ;
        let q = i % N_EQ;
        let r = self.rep.get(cl).copied().unwrap_or(0) as usize;
        let first = self.firsts[(q / N_REST) & 3] as usize;
        let ep = eperm_unrank((first * N_REST + q % N_REST) as u32);
        (mset_labels()[r % N_MSET], ep)
    }

    /// The other indices of the projection of index `i` (its colouring representative has a
    /// non-trivial stabiliser).
    pub fn equivalents(&self, i: usize, out: &mut Vec<usize>) {
        let cl = i / N_EQ;
        let Some(stab) = self.stab.get(cl) else {
            return;
        };
        if stab.is_empty() {
            return;
        }
        let (_, ep) = self.decode(i);
        for &(g, t) in stab {
            out.push(cl * N_EQ + self.qrank(&ep, g.into(), t.into()));
        }
    }

    /// The index of the solved state and its equivalents.
    pub fn targets(&self) -> Vec<usize> {
        let i = self.index(0, &Pos::solved());
        let mut v = vec![i];
        self.equivalents(i, &mut v);
        v.sort_unstable();
        v.dedup();
        v
    }
}

/// Moves between entries of the table.
struct EaGraph<'a> {
    ix: &'a EaIndex,
    mt: &'static OptMoves,
}

impl Graph for EaGraph<'_> {
    const MOVES: usize = crate::cube::N_MOVES;
    type Node = ([u8; 16], [u8; 16]);
    type Cache = ();
    fn node(&self, i: usize, (): &mut ()) -> ([u8; 16], [u8; 16]) {
        self.ix.decode(i)
    }
    #[inline(always)]
    fn step(&self, (xa, ep): &([u8; 16], [u8; 16]), m: usize) -> usize {
        let (xa2, ep2) = self.mt.apply_ea(xa, ep, m);
        self.ix.index_raw(&xa2, &ep2)
    }
    fn equivalents(&self, j: usize, out: &mut Vec<usize>) {
        self.ix.equivalents(j, out);
    }
}

/// The table; `None` if `stop` was set before it was complete.
pub fn build_ea(ix: &EaIndex, verbose: bool, stop: &AtomicBool) -> Option<Table2> {
    let g = EaGraph {
        ix,
        mt: opt_moves(),
    };
    bfs_until(&g, ix.size(), &ix.targets(), "opt-ea", verbose, stop)
}

#[cfg(test)]
#[cfg(target_pointer_width = "64")]
mod tests {
    use super::*;
    use crate::cube::{Cube, SWAP_AB};
    use rand::{Rng, SeedableRng};

    fn ea_of(c: &Cube) -> ([u8; 16], [u8; 16]) {
        let p = Pos::of(c);
        (p.xa, p.ep)
    }

    /// The index is invariant (up to equivalent indices) under the 288 pairs and under `C`, a
    /// decoded entry has its own index, and the second lookup is the first of the `SWAP_AB`
    /// conjugate.
    #[test]
    fn index_is_symmetric() {
        let ix = EaIndex::new().unwrap();
        eprintln!("A-colouring classes: {}", ix.rep.len());
        let s = syms();
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(11);
        for _ in 0..200 {
            let c = Cube::random(&mut rng);
            let (xa, ep) = ea_of(&c);
            let i = ix.index_raw(&xa, &ep);
            let mut eq = vec![i];
            ix.equivalents(i, &mut eq);
            for &(g, t) in &s.pairs {
                let (g, t) = (usize::from(g), usize::from(t));
                let mut l = [0u8; 12];
                l.copy_from_slice(&xa[..12]);
                let la = s.apply_a(g, t, &l);
                let mut xa2 = [0u8; 16];
                xa2[..12].copy_from_slice(&la);
                let ep2 = shuffle(&shuffle(&ix.relab[g], &ep), &ix.gath[t]);
                assert!(eq.contains(&ix.index_raw(&xa2, &ep2)), "pair ({g}, {t})");
            }
            for j in 0..3 {
                let ep2 = shuffle(&ix.cpow[j], &ep);
                assert!(eq.contains(&ix.index_raw(&xa, &ep2)), "C^{j}");
            }
            let (dx, de) = ix.decode(i);
            assert_eq!(ix.index_raw(&dx, &de), i);
            let p = Pos::of(&c);
            let w = Pos::of(&c.conjugate(&SWAP_AB).unwrap());
            assert_eq!(ix.index(1, &p), ix.index(0, &w));
        }
        let _ = rng.gen::<u8>();
    }

    /// Brute force: the distances of the reduced graph (a sparse breadth-first search over the
    /// table's indices with the generator's own steps and equivalents) equal those of the
    /// projections of every position within 5 moves, the projection taken up to the 288 pairs
    /// and `C` by a key made independently of the index (the least of all 864 images as bytes).
    /// (A projection's distance is the least distance of a position with that projection; for a
    /// projection of a position in the ball, the ball holds a witness.)
    #[test]
    fn reduced_graph_matches_brute_force() {
        use std::collections::HashMap;
        const R: u8 = 5;
        let ix = EaIndex::new().unwrap();
        let g = EaGraph {
            ix: &ix,
            mt: opt_moves(),
        };
        let mt = opt_moves();
        let s = syms();
        let key = |p: &Pos| -> [u8; 24] {
            let mut l = [0u8; 12];
            l.copy_from_slice(&p.xa[..12]);
            let mut best = [u8::MAX; 24];
            for &(g, t) in &s.pairs {
                let (g, t) = (usize::from(g), usize::from(t));
                let a = s.apply_a(g, t, &l);
                let mut e = [0u8; 12];
                for q in 0..12 {
                    e[usize::from(s.e_pos[t][q])] = s.e_piece[g][usize::from(p.ep[q])];
                }
                for j in 0..3 {
                    let mut k = [0u8; 24];
                    k[..12].copy_from_slice(&a);
                    for q in 0..12 {
                        k[12 + q] = ix.cpow[j][usize::from(e[q])];
                    }
                    best = best.min(k);
                }
            }
            best
        };
        // positions within R moves
        let mut ball: HashMap<Pos, u8> = HashMap::from([(Pos::solved(), 0)]);
        let mut front = vec![Pos::solved()];
        for d in 1..=R {
            let mut next = vec![];
            for p in &front {
                for m in 0..crate::cube::N_MOVES {
                    let y = mt.apply(p, m);
                    ball.entry(y).or_insert_with(|| {
                        next.push(y);
                        d
                    });
                }
            }
            front = next;
        }
        // the reduced graph to depth R
        let mut red: HashMap<usize, u8> = HashMap::new();
        let mut front: Vec<usize> = ix.targets();
        for &t in &front {
            red.insert(t, 0);
        }
        for d in 1..=R {
            let mut next = vec![];
            let mut eq = vec![];
            for &i in &front {
                let x = g.node(i, &mut ());
                for m in 0..crate::cube::N_MOVES {
                    let j = g.step(&x, m);
                    eq.clear();
                    eq.push(j);
                    g.equivalents(j, &mut eq);
                    for &e in &eq {
                        red.entry(e).or_insert_with(|| {
                            next.push(e);
                            d
                        });
                    }
                }
            }
            front = next;
        }
        let keys: HashMap<Pos, [u8; 24]> = ball.keys().map(|p| (*p, key(p))).collect();
        let mut raw: HashMap<[u8; 24], u8> = HashMap::new();
        for (p, &d) in &ball {
            let e = raw.entry(keys[p]).or_insert(d);
            *e = (*e).min(d);
        }
        for k in 0..2 {
            let mut checked = 0;
            for p in ball.keys() {
                // (lookup 1 is lookup 0 of the conjugate, which is in the ball too)
                let q = if k == 0 {
                    *p
                } else {
                    Pos::of(&p.to_cube().conjugate(&SWAP_AB).unwrap())
                };
                let Some(kq) = keys.get(&q) else { continue };
                assert_eq!(
                    red.get(&ix.index(k, p)).copied(),
                    Some(raw[kq]),
                    "lookup {k}"
                );
                checked += 1;
            }
            assert!(checked > 100_000, "{checked}");
        }
        // (the ball reaches distance R in the reduced graph, and no further)
        assert_eq!(red.values().max(), Some(&R));
    }

    /// Every entry's neighbours under the moves are entries of the same symmetry class as the
    /// neighbours of any equivalent position (the graph is well defined on classes).
    #[test]
    fn moves_respect_classes() {
        let ix = EaIndex::new().unwrap();
        let mt = opt_moves();
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(12);
        for _ in 0..200 {
            let c = Cube::random(&mut rng);
            let (xa, ep) = ea_of(&c);
            let i = ix.index_raw(&xa, &ep);
            let (dx, de) = ix.decode(i);
            for m in 0..crate::cube::N_MOVES {
                let (a1, e1) = mt.apply_ea(&xa, &ep, m);
                let j1 = ix.index_raw(&a1, &e1);
                let mut eq = vec![j1];
                ix.equivalents(j1, &mut eq);
                // the decoded representative's neighbours are the neighbours of some image of
                // the position: their classes appear among the position's neighbours' classes
                let found = (0..crate::cube::N_MOVES).any(|m2| {
                    let (a2, e2) = mt.apply_ea(&dx, &de, m2);
                    eq.contains(&ix.index_raw(&a2, &e2))
                });
                assert!(found, "move {m}");
            }
        }
    }
}
