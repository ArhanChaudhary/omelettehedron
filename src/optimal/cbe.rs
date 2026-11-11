//! The corner x B-centre x edge-triple table of the optimal solver (large tier): the distance
//! to solved of the corners, the B-centre colouring and the positions of the three edges of
//! one A-face, in order.
//!
//! The two-phase solver has a table over the same pieces (`crate::cbe`), but it holds distances
//! to `H_A` and keeps only the set of the three positions and the chirality (440 per row, not
//! 1,320), saturated at 9. This one is a distance to the solved state, exact (mod 3) at every
//! depth, and keeps the order of the three edges.
//!
//! One table serves eight lookups: each A-triple `k` of the position is mapped onto triple 0 by
//! conjugating with a class-preserving symmetry that maps A-face `k` onto A-face 0, and the
//! B-triples through `SWAP_AB` composed with such a symmetry (which turns the A-centres into
//! B-centres as in [`super::coords::CcIndex`]). The distance to solved is invariant under the
//! conjugation (see [`super::sym`]). The corners are reduced by the 6 symmetries that fix
//! A-face 0: 172 corner classes.
//!
//! Entry `(class * N_MSET + colouring) * N_TRIPLE + positions`: 172 x 369,600 x 1,320 =
//! 8.39e10 entries, 21.0 GB.
//!
//! (Chosen by an exact estimate on the recorded medium-tier search of the 13-move set:
//! 18,582 -> 6,350 nodes (0.34x), against 0.62x with the set of positions and the chirality
//! (7.0 GB) and 0.68x for an A-centre/edge table to solved without chirality (8.6 GB).)

use super::coords::{opt_moves, OptMoves, Pos};
use super::sym::*;
use crate::coord::{triple_pieces, N_CORNER, N_MSET};
use crate::cube::N_MOVES;
use crate::fast::{fast, masks, mset_labels, shuffle, FastTables};
use crate::prune::{bfs_until, Graph, Table2};
use std::collections::HashMap;
use std::sync::atomic::AtomicBool;

/// Ordered positions of three distinct edges: 12 x 11 x 10.
pub const N_TRIPLE: usize = 1320;
/// Lookups per position: the four A-triples, then the four A-triples of the `SWAP_AB`
/// conjugate (the B-triples).
pub const N_LOOKUPS: usize = 8;

/// Rank of three distinct positions, in order.
#[inline(always)]
fn triple_rank(p: [u8; 3]) -> usize {
    let [a, b, c] = p.map(usize::from);
    // (distinct positions; saturating so that no input panics)
    let b2 = b.saturating_sub(usize::from(b > a));
    let c2 = c.saturating_sub(usize::from(c > a) + usize::from(c > b));
    (a * 11 + b2) * 10 + c2
}

/// Index of the large tier's corner x B-centre x edge-triple table.
pub struct CbeIndex {
    /// per corner id (padded to 1024) and lookup: (class, relabelling symmetry s,
    /// repositioning symmetry t)
    kt: Box<[[(u16, u8, u8); N_LOOKUPS]; 1024]>,
    /// class -> representative corner id
    pub rep: Vec<u16>,
    /// per class: the (s, t) pairs other than the identity that fix the representative (all
    /// fix A-face 0)
    stab: Vec<Vec<(u8, u8)>>,
    /// per symmetry s (padded to 64): the colour map from the orbit it maps onto the B-centres
    /// (B for a class-preserving s, A for a class-swapping one)
    relab: Box<[[u8; 16]; 64]>,
    /// per symmetry t (padded to 64): gather of the B positions from that orbit's positions
    gath: Box<[[u8; 16]; 64]>,
    /// per symmetry s: the edge pieces that s maps onto triple 0's pieces, in their order
    src: Box<[[u8; 3]; 64]>,
    /// per symmetry s fixing A-face 0: for each piece of triple 0, the index (in triple 0) of
    /// the piece that s maps onto it
    zsrc: Box<[[u8; 3]; 64]>,
    /// per symmetry t: the image of each edge position
    epos: Box<[[u8; 16]; 64]>,
    /// the positions of each triple rank
    unrank: Vec<[u8; 3]>,
    /// per move: the new position of the edge at each position
    pmove: [[u8; 16]; N_MOVES],
    f: &'static FastTables,
}

impl CbeIndex {
    pub fn new() -> CbeIndex {
        let s = sym48();
        let tp = triple_pieces();
        let z = tp[0];
        // the class-preserving symmetries that map A-triple k onto A-triple 0 (6 each)
        let onto0 = |k: usize| -> Vec<usize> {
            (0..N_SYM_CP)
                .filter(|&m| z.contains(&s.epos[m][usize::from(tp[k][0])]))
                .collect()
        };
        let maps: Vec<Vec<usize>> = (0..4).map(onto0).collect();
        // class representative: the smallest corner id among the images under those
        // symmetries, and the first symmetry that gives it (the same classes for every k: the
        // symmetries for k are those fixing A-face 0 times one of them)
        let best = |k: usize, c: usize| -> (u16, usize) {
            maps[k]
                .iter()
                .map(|&m| (s.conj[m][c].0, m))
                .min()
                .unwrap_or((c as u16, 0))
        };
        let mut rep: Vec<u16> = (0..N_CORNER).map(|c| best(0, c).0).collect();
        rep.sort_unstable();
        rep.dedup();
        let cls: HashMap<u16, u16> = rep
            .iter()
            .enumerate()
            .map(|(i, &r)| (r, i as u16))
            .collect();
        #[allow(clippy::large_stack_arrays)] // 40 KB, once
        let mut kt = Box::new([[(0u16, 0u8, 0u8); N_LOOKUPS]; 1024]);
        for c in 0..N_CORNER {
            for k in 0..4 {
                let (r, m) = best(k, c);
                kt[c][k] = (
                    cls.get(&r).copied().unwrap_or(0),
                    m as u8,
                    s.repos(m, c) as u8,
                );
            }
        }
        // the B-triples: conjugate by SWAP_AB (B-centres <-> A-centres), then as above
        let w = s.swap;
        for c in 0..N_CORNER {
            let (cw, _) = s.conj[w][c];
            let tw = s.repos(w, c);
            for k in 0..4 {
                let (cl, m, t) = kt[cw as usize % N_CORNER][k];
                kt[c][4 + k] = (cl, s.mul[m as usize][w], s.mul[t as usize][tw]);
            }
        }
        let stab = rep
            .iter()
            .map(|&r| {
                maps[0]
                    .iter()
                    .filter(|&&m| m != 0 && s.conj[m][r as usize].0 == r)
                    .map(|&m| (m as u8, s.repos(m, r as usize) as u8))
                    .collect()
            })
            .collect();
        let orbit = |m: usize| usize::from(!swaps_classes(m));
        let mut relab = Box::new([[0u8; 16]; 64]);
        let mut gath = Box::new([[0x80u8; 16]; 64]);
        let mut src = Box::new([[0u8; 3]; 64]);
        let mut zsrc = Box::new([[0u8; 3]; 64]);
        let mut epos = Box::new([[0u8; 16]; 64]);
        for m in 0..N_SYM {
            relab[m][..4].copy_from_slice(&s.ccol[m][orbit(m)]);
            for p in 0..12 {
                gath[m][s.cpos[m][orbit(m)][p] as usize] = p as u8;
                epos[m][p] = s.epos[m][p];
            }
            for j in 0..3 {
                let e = (0..12u8)
                    .find(|&e| s.epos[m][usize::from(e)] == z[j])
                    .unwrap_or(0);
                src[m][j] = e;
                zsrc[m][j] = z.iter().position(|&x| x == e).unwrap_or(0) as u8;
            }
        }
        let mut unrank = vec![[0u8; 3]; N_TRIPLE];
        for a in 0..12u8 {
            for b in (0..12u8).filter(|&b| b != a) {
                for c in (0..12u8).filter(|&c| c != a && c != b) {
                    unrank[triple_rank([a, b, c])] = [a, b, c];
                }
            }
        }
        let mt = opt_moves();
        let mut id = [0u8; 16];
        for (p, x) in id.iter_mut().enumerate().take(12) {
            *x = p as u8;
        }
        let pmove = std::array::from_fn(|m| {
            // (after the move, position q holds what was at e2[q])
            let e2 = mt.apply_ep(&id, m);
            let mut pm = [0u8; 16];
            for q in 0..12 {
                pm[usize::from(e2[q] & 15)] = q as u8;
            }
            pm
        });
        CbeIndex {
            kt,
            rep,
            stab,
            relab,
            gath,
            src,
            zsrc,
            epos,
            unrank,
            pmove,
            f: fast(),
        }
    }

    /// Entries: corner classes x B-colourings x ordered positions.
    /// (More than 2^32: 64-bit targets only; saturates elsewhere.)
    pub fn size(&self) -> usize {
        self.rep
            .len()
            .saturating_mul(N_MSET)
            .saturating_mul(N_TRIPLE)
    }

    /// Rank of the colouring `x` after symmetry pair (s, t).
    #[inline(always)]
    fn row(&self, x: &[u8; 16], s: usize, t: usize) -> usize {
        let y = shuffle(&shuffle(&self.relab[s & 63], x), &self.gath[t & 63]);
        self.f.rank(&masks(&y)) as usize
    }

    #[inline(always)]
    fn entry(cl: u16, row: usize, tri: [u8; 3]) -> usize {
        (usize::from(cl) * N_MSET + row) * N_TRIPLE + triple_rank(tri)
    }

    /// Index of lookup `k` (`k < N_LOOKUPS`) of position `p`.
    #[inline(always)]
    pub fn index(&self, k: usize, p: &Pos) -> usize {
        let (cl, s, t) = self.kt[p.c as usize & 1023][k & 7];
        let (s, t) = (usize::from(s), usize::from(t));
        let x = if k < 4 { &p.xb } else { &p.xa };
        let mut at = [0u8; 16];
        for q in 0..12 {
            at[usize::from(p.ep[q] & 15)] = q as u8;
        }
        let src = &self.src[s & 63];
        let ep = &self.epos[t & 63];
        let tri = [0, 1, 2].map(|j| ep[usize::from(at[usize::from(src[j] & 15)] & 15)]);
        Self::entry(cl, self.row(x, s, t), tri)
    }

    /// Index of the projection (corner id `c`, B-colouring `x`, positions `pz` of triple 0's
    /// pieces).
    #[inline(always)]
    fn index_node(&self, c: u16, x: &[u8; 16], pz: [u8; 3]) -> usize {
        let (cl, s, t) = self.kt[c as usize & 1023][0];
        self.entry_sym(cl, s.into(), t.into(), x, pz)
    }

    /// The entry of the projection after a symmetry pair (s, t) that fixes A-face 0.
    #[inline(always)]
    fn entry_sym(&self, cl: u16, s: usize, t: usize, x: &[u8; 16], pz: [u8; 3]) -> usize {
        let zs = &self.zsrc[s & 63];
        let ep = &self.epos[t & 63];
        let tri = [0, 1, 2].map(|j| ep[usize::from(pz[usize::from(zs[j]) % 3] & 15)]);
        Self::entry(cl, self.row(x, s, t), tri)
    }

    /// A projection (corner id, B-colouring, positions of triple 0's pieces) of index `i`.
    pub fn decode(&self, i: usize) -> (u16, [u8; 16], [u8; 3]) {
        let per = N_MSET * N_TRIPLE;
        let c = self.rep.get(i / per).copied().unwrap_or(0);
        let r = i % per;
        (c, mset_labels()[r / N_TRIPLE], self.unrank[r % N_TRIPLE])
    }

    /// The other indices of the projection of index `i` (its corner representative has a
    /// non-trivial stabiliser).
    pub fn equivalents(&self, i: usize, out: &mut Vec<usize>) {
        let cl = i / (N_MSET * N_TRIPLE);
        let Some(stab) = self.stab.get(cl) else {
            return;
        };
        if stab.is_empty() {
            return;
        }
        let (_, x, pz) = self.decode(i);
        for &(s, t) in stab {
            out.push(self.entry_sym(cl as u16, s.into(), t.into(), &x, pz));
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
struct CbeGraph<'a> {
    ix: &'a CbeIndex,
    mt: &'static OptMoves,
}

impl Graph for CbeGraph<'_> {
    const MOVES: usize = N_MOVES;
    type Node = (u16, [u8; 16], [u8; 3]);
    type Cache = ();
    fn node(&self, i: usize, (): &mut ()) -> (u16, [u8; 16], [u8; 3]) {
        self.ix.decode(i)
    }
    #[inline(always)]
    fn step(&self, &(c, x, pz): &(u16, [u8; 16], [u8; 3]), m: usize) -> usize {
        let (c2, x2) = self.mt.apply_cb(c, &x, m);
        let pm = &self.ix.pmove[m & 15];
        self.ix
            .index_node(c2, &x2, pz.map(|p| pm[usize::from(p & 15)]))
    }
    fn equivalents(&self, j: usize, out: &mut Vec<usize>) {
        self.ix.equivalents(j, out);
    }
}

/// The table; `None` if `stop` was set (by a memory guard) before it was complete.
pub fn build_cbe(ix: &CbeIndex, verbose: bool, stop: &AtomicBool) -> Option<Table2> {
    let g = CbeGraph {
        ix,
        mt: opt_moves(),
    };
    bfs_until(&g, ix.size(), &ix.targets(), "opt-cbe", verbose, stop)
}

#[cfg(test)]
#[cfg(target_pointer_width = "64")]
mod tests {
    use super::*;
    use crate::cube::{Cube, SWAP_AB};
    use rand::{Rng, SeedableRng};
    use std::collections::HashSet;

    #[test]
    fn triple_rank_is_a_bijection() {
        let ix = CbeIndex::new();
        let mut seen = HashSet::new();
        for (r, &t) in ix.unrank.iter().enumerate() {
            assert_eq!(triple_rank(t), r);
            assert!(t[0] != t[1] && t[0] != t[2] && t[1] != t[2] && t.iter().all(|&p| p < 12));
            assert!(seen.insert(t));
        }
        assert_eq!(seen.len(), N_TRIPLE);
    }

    /// Every lookup's symmetry maps the corner state onto its class representative, and the
    /// symmetries fixing a representative fix A-face 0.
    #[test]
    fn corner_classes() {
        let ix = CbeIndex::new();
        let s = sym48();
        assert_eq!(ix.rep.len(), 172);
        assert_eq!(ix.size(), 172 * N_MSET * N_TRIPLE);
        let z = triple_pieces()[0];
        for c in 0..N_CORNER {
            for k in 0..N_LOOKUPS {
                let (cl, m, t) = ix.kt[c][k];
                let m = usize::from(m);
                assert_eq!(
                    s.conj[m][c].0,
                    ix.rep[usize::from(cl)],
                    "corner {c} lookup {k}"
                );
                assert_eq!(usize::from(t), s.repos(m, c));
                assert_eq!(swaps_classes(m), k >= 4);
                // the triple of lookup k is mapped onto triple 0
                let mut img: Vec<u8> = ix.src[m]
                    .iter()
                    .map(|&e| s.epos[m][usize::from(e)])
                    .collect();
                img.sort_unstable();
                assert_eq!(img, z.to_vec());
            }
        }
        for st in &ix.stab {
            for &(m, _) in st {
                let mut img: Vec<u8> = z
                    .iter()
                    .map(|&e| s.epos[usize::from(m)][usize::from(e)])
                    .collect();
                img.sort_unstable();
                assert_eq!(img, z.to_vec());
            }
        }
    }

    /// Canonical index of a projection: the smallest of its index and its equivalents.
    fn canon(ix: &CbeIndex, i: usize) -> usize {
        let mut v = vec![i];
        ix.equivalents(i, &mut v);
        v.into_iter().min().unwrap_or(i)
    }

    /// A decoded entry has its own index, a conjugate (by any of the 48 symmetries) has the
    /// same eight lookups up to their order and equivalent indices, and a lookup depends only
    /// on its pieces.
    #[test]
    fn index_is_symmetric() {
        let ix = CbeIndex::new();
        let s = sym48();
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(31);
        for _ in 0..300 {
            let c = Cube::random(&mut rng);
            let p = Pos::of(&c);
            let mut mine: Vec<usize> = (0..N_LOOKUPS)
                .map(|k| {
                    let i = ix.index(k, &p);
                    let (dc, dx, dz) = ix.decode(i);
                    assert_eq!(ix.index_node(dc, &dx, dz), i, "lookup {k}");
                    canon(&ix, i)
                })
                .collect();
            mine.sort_unstable();
            for m in 0..N_SYM {
                let q = Pos::of(&s.conj_by_maps(&c, m));
                let mut theirs: Vec<usize> = (0..N_LOOKUPS)
                    .map(|k| canon(&ix, ix.index(k, &q)))
                    .collect();
                theirs.sort_unstable();
                assert_eq!(mine, theirs, "symmetry {m}");
            }
            // the B-triples are the A-triples of the SWAP_AB conjugate
            let w = Pos::of(&c.conjugate(&SWAP_AB).unwrap());
            for k in 0..4 {
                assert_eq!(ix.index(4 + k, &p), ix.index(k, &w));
            }
            // other A-centres and other edges do not matter for the A-triple lookups
            let mut p2 = p;
            p2.xa = Pos::of(&Cube::random(&mut rng)).xa;
            for k in 0..4 {
                assert_eq!(ix.index(k, &p2), ix.index(k, &p));
            }
        }
        // the solved state's projection is fixed by every symmetry of its corner class
        assert_eq!(ix.targets().len(), 1);
    }

    /// The raw projection of lookup k: corners, the centre colouring the lookup sees and the
    /// positions of the triple's pieces in their order (an independent definition).
    fn raw_key(p: &Pos, k: usize) -> (u16, [u8; 16], [u8; 3]) {
        let q = if k < 4 {
            *p
        } else {
            Pos::of(&p.to_cube().conjugate(&SWAP_AB).unwrap())
        };
        let tp = triple_pieces()[k % 4];
        let at = |e: u8| (0..12u8).find(|&i| q.ep[usize::from(i)] == e).unwrap();
        (q.c, q.xb, tp.map(at))
    }

    /// Brute force: the distances of the reduced graph (a sparse breadth-first search over the
    /// table's indices with the generator's own steps and equivalents) equal the distances of
    /// the raw projections of every position within 5 moves, for all eight lookups. (A
    /// projection's distance is the least distance of a position with that projection; for a
    /// projection of a position in the ball, the ball holds a witness.)
    #[test]
    fn reduced_graph_matches_brute_force() {
        const R: u8 = 5;
        let ix = CbeIndex::new();
        let g = CbeGraph {
            ix: &ix,
            mt: opt_moves(),
        };
        let mt = opt_moves();
        // positions within R moves
        let mut ball: HashMap<Pos, u8> = HashMap::from([(Pos::solved(), 0)]);
        let mut front = vec![Pos::solved()];
        for d in 1..=R {
            let mut next = vec![];
            for p in &front {
                for m in 0..N_MOVES {
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
                for m in 0..N_MOVES {
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
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(32);
        for k in 0..N_LOOKUPS {
            let mut raw: HashMap<(u16, [u8; 16], [u8; 3]), u8> = HashMap::new();
            for (p, &d) in &ball {
                let e = raw.entry(raw_key(p, k)).or_insert(d);
                *e = (*e).min(d);
            }
            let mut checked = 0;
            for p in ball.keys() {
                // (all of the ball for lookup 0, a sample for the others)
                if k > 0 && rng.gen_range(0..8) != 0 {
                    continue;
                }
                let want = raw[&raw_key(p, k)];
                let got = red.get(&ix.index(k, p)).copied();
                assert_eq!(got, Some(want), "lookup {k}");
                checked += 1;
            }
            assert!(checked > 10_000);
        }
    }
}

/// Tests of the generated table (21.0 GB, read from `$FTO_TABLES`).
#[cfg(test)]
#[cfg(target_pointer_width = "64")]
mod table_tests {
    use super::*;
    use rand::{Rng, SeedableRng};

    /// The table is a distance table: from sampled entries, a walk down (each step to a
    /// neighbour one closer, which the mod-3 values identify) reaches the solved entry within
    /// the largest distance, equivalent entries hold the same value, and every entry was
    /// reached.
    #[test]
    #[ignore = "needs the large tier's table (21.0 GB; `omelettehedron optimal gen --tier large`)"]
    fn table_walks_down_from_sampled_entries() {
        let ix = CbeIndex::new();
        let path = crate::tables::table_dir(None).join("opt-cbe-sym6.prun");
        let t = Table2::load(&path, ix.size()).unwrap();
        assert_eq!(t.hist.iter().sum::<u64>(), ix.size() as u64);
        let g = CbeGraph {
            ix: &ix,
            mt: opt_moves(),
        };
        let targets = ix.targets();
        let max = t.hist.len() - 1;
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(33);
        let mut eq = vec![];
        let mut depths = vec![0u64; max + 1];
        for _ in 0..3000 {
            let mut i = rng.gen_range(0..ix.size());
            let mut v = t.get(i);
            let v0 = v;
            assert!(v < 3, "entry {i} not reached");
            let mut d = 0;
            while !targets.contains(&i) {
                assert!(d < max, "no walk down from entry {i}");
                eq.clear();
                ix.equivalents(i, &mut eq);
                assert!(eq.iter().all(|&j| t.get(j) == v));
                let x = g.node(i, &mut ());
                let want = (v + 2) % 3;
                let next = (0..N_MOVES)
                    .map(|m| g.step(&x, m))
                    .find(|&j| t.get(j) == want);
                i = next.unwrap();
                v = want;
                d += 1;
            }
            assert_eq!(d % 3, usize::from(v0));
            depths[d] += 1;
        }
        eprintln!("walk lengths of the sampled entries: {depths:?}");
    }
}
