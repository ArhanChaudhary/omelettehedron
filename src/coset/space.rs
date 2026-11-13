//! A coset's 933,120 elements as the states of `H_A` they are indexed by (the phase-2 index
//! `corner id * 972 + edge id` of [`Phase2`]), the depth at which each element is first known to
//! be solvable, and the prepass that extends those marks by one class-A move.
//!
//! The marks are one byte per element (912 KB) rather than one bitmap per depth: a byte holds
//! the whole distance distribution, and the coset is small enough that every pass over it takes
//! milliseconds (the 933,120 x 8 neighbours of a full prepass are about 10 ms on one thread).

use crate::coord::*;
use crate::cube::*;
use crate::tables::{Phase2, N_P2_EDGE};

/// Elements of a coset: `|H_A|` = 960 corner states x 972 edge states.
pub const N_COSET: usize = N_CORNER * N_P2_EDGE;
/// Class-A moves (internal BR, BL, Bw, Dw and their inverses): the moves of the prepass.
pub const N_A: usize = 8;
/// The mark of an element with no known solution.
pub const UNSET: u8 = u8::MAX;

/// The right action of the class-A moves on the phase-2 index. Corners and edges move
/// independently, so it is two small tables (960 and 972 rows).
pub struct Space {
    /// the internal class-A moves, in the order of the table columns
    pub moves: [usize; N_A],
    corner: Vec<[u16; N_A]>,
    edge: Vec<[u16; N_A]>,
}

impl Space {
    pub fn new(p2: &Phase2) -> Space {
        let mut moves = [0usize; N_A];
        for (k, m) in (0..N_MOVES).filter(|&m| is_a_move(m)).enumerate().take(N_A) {
            moves[k] = m;
        }
        let mt = move_tables();
        let corner =
            mt.c.iter()
                .map(|row| std::array::from_fn(|k| row[moves[k]]))
                .collect();
        let edge = p2
            .emoves
            .iter()
            .map(|row| std::array::from_fn(|k| row[moves[k]]))
            .collect();
        Space {
            moves,
            corner,
            edge,
        }
    }

    /// An optimal phase-2 solution (class-A moves) of the state with index `p`: a walk down
    /// the phase-2 distances on the index (`None` if `p` is out of range).
    pub fn solve(&self, p2: &Phase2, p: u32) -> Option<Vec<usize>> {
        let mut i = p;
        let mut d = *p2.dist.get(i as usize)?;
        let mut out = vec![];
        while d > 0 {
            let k = (0..N_A).find(|&k| p2.dist[self.step(i, k) as usize] + 1 == d)?;
            out.push(self.moves[k]);
            i = self.step(i, k);
            d -= 1;
        }
        Some(out)
    }

    /// The index of state `i` after class-A move number `k` (move `self.moves[k]`).
    #[inline(always)]
    pub fn step(&self, i: u32, k: usize) -> u32 {
        let (c, e) = (i as usize / N_P2_EDGE, i as usize % N_P2_EDGE);
        u32::from(self.corner[c][k & 7]) * N_P2_EDGE as u32 + u32::from(self.edge[e][k & 7])
    }
}

/// Upper bounds of a coset's elements: the depth at which each element was first marked
/// (`UNSET`: none yet), with the elements first marked at each depth (the next prepass's
/// sources).
#[derive(Clone)]
pub struct Marks {
    pub d: Vec<u8>,
    levels: Vec<Vec<u32>>,
}

impl Marks {
    pub fn new() -> Marks {
        Marks {
            d: vec![UNSET; N_COSET],
            levels: vec![],
        }
    }

    /// Mark element `i` at depth `n` (the end state of a phase-1 sequence of `n` moves) if it
    /// has no mark yet; whether it was new.
    #[inline]
    pub fn mark(&mut self, i: u32, n: usize) -> bool {
        let Some(v) = self.d.get_mut(i as usize) else {
            return false;
        };
        if *v != UNSET {
            return false;
        }
        *v = n as u8;
        if self.levels.len() <= n {
            self.levels.resize(n + 1, vec![]);
        }
        self.levels[n].push(i);
        true
    }

    /// The prepass to depth `n`: every element first marked at `n - 1` extended by each class-A
    /// move (elements marked earlier were extended by earlier prepasses). Returns the number of
    /// new marks.
    pub fn prepass(&mut self, sp: &Space, n: usize) -> u64 {
        let Some(src) = n.checked_sub(1).and_then(|p| self.levels.get(p)) else {
            return 0;
        };
        let mut new = vec![];
        for &i in src {
            for k in 0..N_A {
                let j = sp.step(i, k);
                let v = &mut self.d[j as usize];
                if *v == UNSET {
                    *v = n as u8;
                    new.push(j);
                }
            }
        }
        let count = new.len() as u64;
        if self.levels.len() <= n {
            self.levels.resize(n + 1, vec![]);
        }
        self.levels[n].extend(new);
        count
    }

    /// Number of marked elements.
    pub fn count(&self) -> u64 {
        self.levels.iter().map(|l| l.len() as u64).sum()
    }

    /// The bounds if no phase-1 sequence were added any more: every mark extended by class-A
    /// moves as far as it reaches (all later prepasses at once). The largest is the depth at
    /// which every element would be covered.
    pub fn project(&self, sp: &Space) -> Vec<u8> {
        let mut p = self.d.clone();
        let mut buckets = self.levels.clone();
        let mut n = 0;
        while n < buckets.len() {
            let cur = std::mem::take(&mut buckets[n]);
            let next = (n + 1).min(usize::from(UNSET - 1)) as u8;
            let mut grown = vec![];
            for &i in &cur {
                // (an endpoint marked here may have been reached earlier through another one)
                if usize::from(p[i as usize]) != n {
                    continue;
                }
                for k in 0..N_A {
                    let j = sp.step(i, k);
                    let v = &mut p[j as usize];
                    if *v > next {
                        *v = next;
                        grown.push(j);
                    }
                }
            }
            if !grown.is_empty() {
                if buckets.len() <= n + 1 {
                    buckets.resize(n + 2, vec![]);
                }
                buckets[n + 1].extend(grown);
            }
            n += 1;
        }
        p
    }
}

impl Default for Marks {
    fn default() -> Self {
        Self::new()
    }
}

/// The B-centres of `H_A` told apart: label `l` (0..12) is the centre that starts at B-centre
/// position `l`, of colour `l / 3`. In `H_A` the B-centres travel with the corners, so each
/// corner state has one arrangement of labels; this is it, by corner id (`None` if that were
/// not so: checked over the whole group as it is built).
pub fn labelled_attachment() -> Option<&'static [[u8; 12]]> {
    static A: std::sync::OnceLock<Option<Vec<[u8; 12]>>> = std::sync::OnceLock::new();
    A.get_or_init(|| {
        let s0 = Cube {
            xb: std::array::from_fn(|i| i as u8),
            ..Cube::solved()
        };
        let cs = corners();
        let mut att: Vec<Option<[u8; 12]>> = vec![None; N_CORNER];
        att[usize::from(cs.id(&s0))] = Some(s0.xb);
        let mut q = vec![s0];
        let mut k = 0;
        while k < q.len() {
            let c = q[k];
            k += 1;
            for m in (0..N_MOVES).filter(|&m| is_a_move(m)) {
                let n = c.apply(m);
                let slot = att.get_mut(usize::from(cs.id(&n)))?;
                match slot {
                    None => {
                        *slot = Some(n.xb);
                        q.push(n);
                    }
                    Some(xb) if *xb != n.xb => return None,
                    Some(_) => {}
                }
            }
        }
        att.into_iter().collect()
    })
    .as_deref()
}

/// The position with its B-centres labelled as [`labelled_attachment`] labels them: the centres
/// of each colour get that colour's three labels in the order of their positions. (Any order
/// gives one of the position's lifts; this is the one the coset solver uses.)
pub fn lift(cube: &Cube) -> Cube {
    let mut next = [0u8; 4];
    let mut xb = [0u8; 12];
    for (p, &c) in cube.xb.iter().enumerate() {
        let c = usize::from(c) & 3;
        xb[p] = 3 * c as u8 + next[c];
        next[c] += 1;
    }
    Cube { xb, ..*cube }
}

/// The position of a labelled lift (colours from the labels).
pub fn unlabel(cube: &Cube) -> Cube {
    Cube {
        xb: cube.xb.map(|l| l / 3),
        ..*cube
    }
}

/// Number of elements at each value of `d` (index `UNSET` excluded: the elements with no mark
/// are `N_COSET` minus the sum).
pub fn histogram(d: &[u8]) -> Vec<u64> {
    let mut h = vec![0u64; 64];
    for &v in d {
        if v != UNSET {
            if let Some(x) = h.get_mut(usize::from(v)) {
                *x += 1;
            }
        }
    }
    while h.last() == Some(&0) {
        h.pop();
    }
    h
}
