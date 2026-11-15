//! The joint check of the two phases: every position within 35 moves, from the relaxed phase 1
//! (at most 11 moves to the phase-2 space) and the complete phase-2 search (at most 25).
//!
//! A position `p` with relaxed phase-1 distance `d1(p) <= 10` is within `10 + 25 = 35`. One with
//! `d1(p) = 11` is within 35 as soon as one of its optimal phase-1 solutions ends outside the
//! phase-2 search's last layer `S25` (distance 25): the rest is at most 24. The positions where
//! no optimal phase-1 solution does are among the *candidates* `q = s v`: `s` in `S25`, `v` 11
//! moves along which `d1` rises by one at each step (the reverse of an optimal solution of `q`
//! that ends at `s`). For each candidate the check looks for an optimal phase-1 solution that
//! ends outside `S25`; the candidates without one are listed for an explicit solution.
//!
//! * [`ae_depths`]: exact relaxed phase-1 distances of every entry of the A-centre/edge index
//!   (288 symmetries, edge chirality), by breadth-first search from the relaxed targets;
//! * [`geodesic`]: the entries on a shortest path from a target to an entry at the largest
//!   distance (a backward search from those, one distance down at a time): every step of an
//!   upward path `v` and of an optimal solution of a candidate is in it;
//! * [`Phase2Slots`]: membership of a fixed-frame position in a phase-2 layer, through the
//!   search's slot index and the layer's files;
//! * [`Joint`]: the enumeration from one position of each symmetry class of `S25` (the 24
//!   symmetries that keep the face classes conjugate positions, moves, `d1` and `S25` alike, so
//!   candidates and their resolution are conjugate too).

use super::engine::{for_slots, BfsError, Layout, Manifest};
use super::instance::{Instance, Kind, N_SYM};
use super::solve::{PState, Views};
use crate::chen::relaxed::RelaxedAe;
use crate::chen::space::space;
use crate::chen::tables::P2;
use crate::coord::{corners, mset, Coords, N_CORNER};
use crate::cube::{Cube, N_MOVES};
use crate::prune::Graph;
use crate::tables::{AeGraph, AeIndex, AeSym};
use rayon::prelude::*;
use std::collections::HashSet;
use std::path::Path;
use std::sync::atomic::{AtomicU64, AtomicU8, Ordering};
use std::sync::Mutex;

/// An entry not reached (yet).
const UNSET: u8 = u8::MAX;
/// Entries per parallel task of the scans.
const CHUNK: usize = 1 << 16;

/// Exact relaxed phase-1 distance of every entry, and the entries per distance.
pub struct AeDepths {
    /// distance of each entry
    pub depth: Vec<AtomicU8>,
    /// entries per distance
    pub hist: Vec<u64>,
}

impl AeDepths {
    /// The distance of entry `i`.
    #[inline]
    pub fn get(&self, i: usize) -> u8 {
        self.depth
            .get(i)
            .map_or(UNSET, |d| d.load(Ordering::Relaxed))
    }
}

fn atomic_u8(n: usize, v: u8) -> Vec<AtomicU8> {
    (0..n).into_par_iter().map(|_| AtomicU8::new(v)).collect()
}

/// Exact distances of all entries of `ix` to `targets` (sorted, with their equivalents), in
/// the face-turn metric: a breadth-first search one distance at a time, forward from each
/// entry at the current distance.
pub fn ae_depths(ix: &AeSym, targets: &[usize], verbose: bool) -> AeDepths {
    let g = AeGraph::<AeSym, N_MOVES>::new(ix);
    let n = ix.size();
    let depth = atomic_u8(n, UNSET);
    for &t in targets {
        if let Some(d) = depth.get(t) {
            d.store(0, Ordering::Relaxed);
        }
    }
    let mut hist = vec![targets.len() as u64];
    for d in 0u8.. {
        let found = AtomicU64::new(0);
        (0..n.div_ceil(CHUNK)).into_par_iter().for_each(|k| {
            let mut cache = None;
            let mut eq = Vec::new();
            let mut set = 0u64;
            for i in k * CHUNK..((k + 1) * CHUNK).min(n) {
                if depth[i].load(Ordering::Relaxed) != d {
                    continue;
                }
                g.neighbours(i, &mut cache, |j| {
                    if depth[j].load(Ordering::Relaxed) == UNSET {
                        depth[j].store(d + 1, Ordering::Relaxed);
                        set += 1;
                        eq.clear();
                        g.equivalents(j, &mut eq);
                        for &e in &eq {
                            if let Some(x) = depth.get(e) {
                                x.store(d + 1, Ordering::Relaxed);
                            }
                        }
                    }
                });
            }
            found.fetch_add(set, Ordering::Relaxed);
        });
        if found.load(Ordering::Relaxed) == 0 {
            break;
        }
        // (the count above may count an entry twice when two threads reach it at once)
        let c = (0..n.div_ceil(CHUNK))
            .into_par_iter()
            .map(|k| {
                (k * CHUNK..((k + 1) * CHUNK).min(n))
                    .filter(|&i| depth[i].load(Ordering::Relaxed) == d + 1)
                    .count() as u64
            })
            .sum();
        hist.push(c);
        if verbose {
            eprintln!("relaxed phase 1: depth {} {c} entries", d + 1);
        }
    }
    AeDepths { depth, hist }
}

/// The entries on a shortest path from a target to an entry at the largest distance.
pub struct Geodesic {
    bits: Vec<AtomicU64>,
    /// the largest distance
    pub top: u8,
    /// entries per distance
    pub per_level: Vec<u64>,
}

impl Geodesic {
    /// Whether entry `i` is in the set.
    #[inline]
    pub fn contains(&self, i: usize) -> bool {
        self.bits
            .get(i / 64)
            .is_some_and(|w| w.load(Ordering::Relaxed) >> (i % 64) & 1 == 1)
    }

    fn insert(&self, i: usize) {
        if let Some(w) = self.bits.get(i / 64) {
            w.fetch_or(1 << (i % 64), Ordering::Relaxed);
        }
    }
}

/// The geodesic set to the entries at the largest distance of `dep`: those, then per distance
/// below the neighbours (and their equivalents) of the set one distance up.
pub fn geodesic(ix: &AeSym, dep: &AeDepths) -> Geodesic {
    let g = AeGraph::<AeSym, N_MOVES>::new(ix);
    let n = ix.size();
    let top = (dep.hist.len() - 1) as u8;
    let geo = Geodesic {
        bits: (0..n.div_ceil(64))
            .into_par_iter()
            .map(|_| AtomicU64::new(0))
            .collect(),
        top,
        per_level: vec![0; top as usize + 1],
    };
    (0..n).into_par_iter().for_each(|i| {
        if dep.get(i) == top {
            geo.insert(i);
        }
    });
    for k in (0..top).rev() {
        (0..n.div_ceil(CHUNK)).into_par_iter().for_each(|c| {
            let mut cache = None;
            let mut eq = Vec::new();
            for i in c * CHUNK..((c + 1) * CHUNK).min(n) {
                if dep.get(i) != k + 1 || !geo.contains(i) {
                    continue;
                }
                g.neighbours(i, &mut cache, |j| {
                    if dep.get(j) == k && !geo.contains(j) {
                        geo.insert(j);
                        eq.clear();
                        g.equivalents(j, &mut eq);
                        for &e in &eq {
                            geo.insert(e);
                        }
                    }
                });
            }
        });
    }
    let mut per_level = vec![0u64; top as usize + 1];
    let counts: Vec<Vec<u64>> = (0..n.div_ceil(CHUNK))
        .into_par_iter()
        .map(|c| {
            let mut h = vec![0u64; top as usize + 1];
            for i in c * CHUNK..((c + 1) * CHUNK).min(n) {
                if geo.contains(i) {
                    if let Some(x) = h.get_mut(dep.get(i) as usize) {
                        *x += 1;
                    }
                }
            }
            h
        })
        .collect();
    for h in counts {
        for (a, b) in per_level.iter_mut().zip(h) {
            *a += b;
        }
    }
    Geodesic { per_level, ..geo }
}

/// Membership of fixed-frame positions in a layer of the phase-2 search.
pub struct Phase2Slots {
    /// the full instance
    pub inst: Instance,
    views: Views,
    /// the instance's corner state of each pattern-view corner state
    pc_to_p: Vec<u16>,
    /// every slot of the layer (with the slots equivalent under the stabilisers)
    slots: HashSet<(u32, u32, u32)>,
    /// one slot per position up to the 24 symmetries
    pub reps: Vec<(u32, u32, u32)>,
}

impl Phase2Slots {
    /// Read layer `d` of the search in `dir`.
    ///
    /// # Errors
    ///
    /// If the manifest or a layer file cannot be read.
    pub fn load(dir: &Path, d: usize) -> Result<Phase2Slots, BfsError> {
        let m = Manifest::load(dir)?.ok_or_else(|| BfsError::Manifest {
            path: dir.display().to_string(),
            what: "no search here".into(),
        })?;
        let inst = Instance::new(Kind::FULL);
        let lay = Layout::of(dir, &m);
        let mut slots = HashSet::new();
        let mut reps = Vec::new();
        for_slots(&lay, &inst, d, 0..inst.n_class, |c, x, p| {
            slots.insert((c as u32, x as u32, p as u32));
            if super::engine::canonical(&inst, c, x, p) {
                reps.push((c as u32, x as u32, p as u32));
            }
        })?;
        let mut pc_to_p = vec![0u16; inst.p_to_pc.len()];
        for (p, &pc) in inst.p_to_pc.iter().enumerate() {
            if let Some(v) = pc_to_p.get_mut(pc as usize) {
                *v = p as u16;
            }
        }
        Ok(Phase2Slots {
            inst,
            views: Views::new(),
            pc_to_p,
            slots,
            reps,
        })
    }

    /// Slots of the layer.
    pub fn len(&self) -> usize {
        self.slots.len()
    }

    /// Whether the layer is empty.
    pub fn is_empty(&self) -> bool {
        self.slots.is_empty()
    }

    /// The slot of a phase-2 position (fixed frame): its configuration's class, and its
    /// B-centres and corners seen in the pattern view, conjugated to the class representative.
    pub fn slot_of(&self, p2: &P2) -> (u32, u32, u32) {
        let s = space();
        let inst = &self.inst;
        let cfg = p2.cfg as usize;
        let q = s.pat.of_cfg[cfg] as usize;
        let (cls, k) = s.pat.class[q];
        let k = k as usize;
        let r = s.rot[cfg] as usize;
        let pc = s.corners.of_fixed[r * N_CORNER + p2.c as usize] as usize;
        let p = inst.ps[self.pc_to_p[pc] as usize * N_SYM + k];
        let xb = s.pattern_xb(p2.cfg, &p2.xb);
        let x = mset().rank(&std::array::from_fn(|i| xb[i]));
        let x = inst.xs[k * inst.n_x + x as usize];
        (u32::from(cls), x, u32::from(p))
    }

    /// Whether a fixed-frame position is in the layer (`None` if it is no phase-2 position).
    pub fn contains(&self, cube: &Cube) -> Option<bool> {
        let p2 = P2::of(cube)?;
        Some(self.slots.contains(&self.slot_of(&p2)))
    }

    /// The fixed-frame position of slot `(c, x, p)` (the class's representative configuration).
    pub fn cube(&self, (c, x, p): (u32, u32, u32)) -> Option<Cube> {
        let st = PState::of_slot(&self.inst, c as usize, x as usize, p as usize)?;
        let p2 = self.views.to_p2(&st);
        let (xa, ep) = *space().cfgs.get(p2.cfg as usize)?;
        let (cp, co) = *corners().id_to_state.get(p2.c as usize)?;
        let mut xb = [0u8; 12];
        xb.copy_from_slice(&p2.xb[..12]);
        Some(Cube { cp, co, ep, xa, xb })
    }
}

/// The counts and the unresolved candidates of a joint check.
#[derive(Default)]
pub struct JointResult {
    /// representatives of `S25` processed
    pub starts: u64,
    /// those whose A-centres and edges are in the geodesic set (at distance 0)
    pub starts_on_geodesics: u64,
    /// candidates (an upward path from a start; a position reached by two paths counts twice)
    pub candidates: u64,
    /// candidates with an optimal phase-1 solution ending outside the layer
    pub resolved: u64,
    /// nodes of the downward searches
    pub down_nodes: u64,
    /// candidates without one (each once)
    pub unresolved: Vec<Cube>,
    /// optimal solutions found ending in a position that is no phase-2 position (must be 0)
    pub anomalies: u64,
}

/// The joint check's state: the A-centre/edge index, distances, geodesic set and the layer.
pub struct Joint<'a> {
    /// the index
    pub ix: &'a AeSym,
    /// exact distances
    pub dep: &'a AeDepths,
    /// the geodesic set
    pub geo: &'a Geodesic,
    /// the last layer
    pub layer: &'a Phase2Slots,
}

impl Joint<'_> {
    /// The entry of a position.
    #[inline]
    pub fn entry(&self, c: &Cube) -> usize {
        self.ix.index(&Coords::of(c))
    }

    /// Whether `c` is on a geodesic at distance `k`.
    #[inline]
    fn on(&self, c: &Cube, k: u8) -> bool {
        let i = self.entry(c);
        self.geo.contains(i) && self.dep.get(i) == k
    }

    /// Whether some optimal phase-1 solution of `c` (at distance `k`, on the geodesic set) ends
    /// outside the layer; counts the nodes. A solution ending in no phase-2 position counts as
    /// an anomaly (and not as outside).
    fn resolve(&self, c: &Cube, k: u8, nodes: &mut u64, anomalies: &mut u64) -> bool {
        *nodes += 1;
        if k == 0 {
            let Some(inside) = self.layer.contains(c) else {
                *anomalies += 1;
                return false;
            };
            return !inside;
        }
        (0..N_MOVES).any(|m| {
            let n = c.apply(m);
            self.on(&n, k - 1) && self.resolve(&n, k - 1, nodes, anomalies)
        })
    }

    /// The candidates above `c` (at distance `k`, on the geodesic set): each upward path to the
    /// top, the candidate resolved or kept.
    fn up(&self, c: &Cube, k: u8, out: &mut JointResult, seen: &mut HashSet<Cube>) {
        if k == self.geo.top {
            out.candidates += 1;
            let (mut nodes, mut anomalies) = (0, 0);
            if self.resolve(c, k, &mut nodes, &mut anomalies) {
                out.resolved += 1;
            } else if seen.insert(*c) {
                out.unresolved.push(*c);
            }
            out.down_nodes += nodes;
            out.anomalies += anomalies;
            return;
        }
        for m in 0..N_MOVES {
            let n = c.apply(m);
            if self.on(&n, k + 1) {
                self.up(&n, k + 1, out, seen);
            }
        }
    }

    /// Run the check from one position of each symmetry class of the layer (all threads).
    pub fn run(&self, verbose: bool) -> JointResult {
        let total = Mutex::new(JointResult::default());
        let done = AtomicU64::new(0);
        let reps = &self.layer.reps;
        reps.par_iter().for_each(|&slot| {
            let mut r = JointResult {
                starts: 1,
                ..JointResult::default()
            };
            let mut seen = HashSet::new();
            match self.layer.cube(slot) {
                Some(s) if self.on(&s, 0) => {
                    r.starts_on_geodesics = 1;
                    self.up(&s, 0, &mut r, &mut seen);
                }
                Some(_) => {}
                None => r.anomalies += 1,
            }
            let mut t = total
                .lock()
                .unwrap_or_else(std::sync::PoisonError::into_inner);
            t.starts += r.starts;
            t.starts_on_geodesics += r.starts_on_geodesics;
            t.candidates += r.candidates;
            t.resolved += r.resolved;
            t.down_nodes += r.down_nodes;
            t.anomalies += r.anomalies;
            t.unresolved.extend(r.unresolved);
            let k = done.fetch_add(1, Ordering::Relaxed) + 1;
            if verbose && k.is_multiple_of(5000) {
                eprintln!(
                    "joint: {k} of {} starts, {} candidates, {} unresolved",
                    reps.len(),
                    t.candidates,
                    t.unresolved.len()
                );
            }
        });
        let mut t = total
            .into_inner()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        // (a candidate reached from two starts of one class is listed once)
        let mut uniq = HashSet::new();
        t.unresolved.retain(|c| uniq.insert(*c));
        t
    }
}

/// The relaxed table's distance and the exact distance array agree on random positions, and
/// on their conjugates by the 24 symmetries that keep the face classes (which the index's 288
/// include): returns the positions checked.
///
/// # Errors
///
/// The first disagreement found, described.
pub fn check_d1(ix: &AeSym, relaxed: &RelaxedAe, dep: &AeDepths, n: u64) -> Result<u64, String> {
    use rand::SeedableRng;
    let f = crate::fast::fast();
    let s48 = crate::optimal::sym::sym48();
    (0..n)
        .into_par_iter()
        .map(|seed| {
            let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(0x1d1 + seed);
            let c = Cube::random(&mut rng);
            let d = relaxed
                .exact(ix, &f.of(&c))
                .ok_or_else(|| format!("seed {seed}: the table does not walk down"))?;
            let e = dep.get(ix.index(&Coords::of(&c)));
            if d != e {
                return Err(format!("seed {seed}: table {d}, exact array {e}"));
            }
            for s in 0..N_SYM {
                let (cc, _) = s48
                    .conj_cube(&c, s)
                    .ok_or_else(|| format!("seed {seed}: no conjugate by {s}"))?;
                let ds = dep.get(ix.index(&Coords::of(&cc)));
                if ds != d {
                    return Err(format!("seed {seed}: d1 {d}, conjugate by {s}: {ds}"));
                }
            }
            Ok(1)
        })
        .sum()
}

/// Whether the geodesic set is exactly the backward closure: every entry at the largest
/// distance is in it, and an entry at distance `k` below is in it exactly when a neighbour at
/// distance `k + 1` is (a check of [`geodesic`] over all entries). Returns the entries that
/// break the rule.
pub fn check_geodesic(ix: &AeSym, dep: &AeDepths, geo: &Geodesic) -> u64 {
    let g = AeGraph::<AeSym, N_MOVES>::new(ix);
    let n = ix.size();
    (0..n.div_ceil(CHUNK))
        .into_par_iter()
        .map(|c| {
            let mut cache = None;
            let mut bad = 0u64;
            for i in c * CHUNK..((c + 1) * CHUNK).min(n) {
                let d = dep.get(i);
                let want = if d == geo.top {
                    true
                } else {
                    let mut up = false;
                    g.neighbours(i, &mut cache, |j| {
                        up |= dep.get(j) == d + 1 && geo.contains(j);
                    });
                    up
                };
                bad += u64::from(want != geo.contains(i));
            }
            bad
        })
        .sum()
}

impl Joint<'_> {
    /// The candidates above one start (at most `max`).
    pub fn candidates_from(&self, s: &Cube, max: usize) -> Vec<Cube> {
        fn up(j: &Joint, c: &Cube, k: u8, out: &mut Vec<Cube>, max: usize) {
            if out.len() >= max {
                return;
            }
            if k == j.geo.top {
                out.push(*c);
                return;
            }
            for m in 0..N_MOVES {
                let n = c.apply(m);
                if j.on(&n, k + 1) {
                    up(j, &n, k + 1, out, max);
                }
            }
        }
        let mut out = Vec::new();
        if self.on(s, 0) {
            up(self, s, 0, &mut out, max);
        }
        out
    }

    /// An optimal phase-1 solution of `c` (at the largest distance) that ends outside the
    /// layer, and its end.
    pub fn alternative(&self, c: &Cube) -> Option<(Vec<usize>, Cube)> {
        fn down(j: &Joint, c: &Cube, k: u8, path: &mut Vec<usize>) -> Option<Cube> {
            if k == 0 {
                return (j.layer.contains(c) == Some(false)).then_some(*c);
            }
            for m in 0..N_MOVES {
                let n = c.apply(m);
                if j.on(&n, k - 1) {
                    path.push(m);
                    if let Some(e) = down(j, &n, k - 1, path) {
                        return Some(e);
                    }
                    path.pop();
                }
            }
            None
        }
        let mut path = Vec::new();
        let e = down(self, c, self.geo.top, &mut path)?;
        Some((path, e))
    }
}

impl Joint<'_> {
    /// All optimal phase-1 solutions of `c` (at the largest distance): how many end inside the
    /// layer and how many outside (paths, not distinct ends).
    pub fn count_ends(&self, c: &Cube) -> (u64, u64) {
        fn down(j: &Joint, c: &Cube, k: u8, n: &mut (u64, u64)) {
            if k == 0 {
                if j.layer.contains(c) == Some(true) {
                    n.0 += 1;
                } else {
                    n.1 += 1;
                }
                return;
            }
            for m in 0..N_MOVES {
                let x = c.apply(m);
                if j.on(&x, k - 1) {
                    down(j, &x, k - 1, n);
                }
            }
        }
        let mut n = (0, 0);
        down(self, c, self.geo.top, &mut n);
        n
    }
}
