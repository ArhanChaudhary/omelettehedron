//! Phase-1 enumeration for the coset solver: every canonical sequence of exactly `len` moves
//! that takes the coset into `H_A` and ends in a class-B move, and the `H_A` state each one
//! ends in.
//!
//! The search is the two-phase solver's phase 1 on one side only (target `H_A`), stripped of
//! everything that is not a phase-1 lower bound: no class-B bound, no distance-to-solved bound
//! and no best length, since every sequence counts, not only the ones that improve a solution.
//! A child is cut only when its exact A-centre/edge distance, its exact corner/B-centre
//! distance or (large and huge tiers, with at most 8 moves left) the joint table's bound
//! exceeds the moves left: all of them are distances of projections of the position, so no
//! sequence that reaches `H_A` is ever cut. The last move must be a class-B move: a sequence
//! ending in a class-A move is a shorter sequence plus that move, which the prepass covers
//! ([`super::space::Marks::prepass`]). A leaf is checked in full ([`Phase2::index`]), since the
//! tables only bound the distance to `H_A`.
//!
//! The work is split into the subtrees below depth 3 or 4 (as the two-phase search does),
//! handed to rayon's workers.

use crate::cube::*;
use crate::fast::*;
use crate::prune::decode;
use crate::tables::*;
use rayon::prelude::*;
use std::sync::atomic::{AtomicBool, Ordering};

/// Longest phase-1 sequence the enumeration handles.
pub const MAX_LEN: usize = 32;

/// With at most this many moves left the joint table is read (as in the two-phase search,
/// where it adds information only near the leaves).
const CBE_BUDGET: usize = 8;

/// The class-A moves as a bit mask over the 16 internal moves.
const A_MASK: u32 = {
    let mut a = 0u32;
    let mut m = 0;
    while m < N_MOVES {
        let s = m >> 1;
        if s == A_SLOTS[0] || s == A_SLOTS[1] || s == A_SLOTS[2] || s == A_SLOTS[3] {
            a |= 1 << m;
        }
        m += 1;
    }
    a
};

/// The moves that can end a phase-1 sequence: the class-B moves.
const B_MASK: u32 = ((1 << N_MOVES) - 1) & !A_MASK;

/// Moves allowed after each move by [`allowed_after`] (index `N_MOVES`: no previous move).
static ALLOWED: [u32; N_MOVES + 1] = {
    let mut a = [0u32; N_MOVES + 1];
    let mut last = 0;
    while last <= N_MOVES {
        let prev = if last == N_MOVES { MAX_MOVES } else { last };
        let mut m = 0;
        while m < N_MOVES {
            if allowed_after(false, prev, m) {
                a[last] |= 1 << m;
            }
            m += 1;
        }
        last += 1;
    }
    a
};

/// The end states of the sequences found, and the sequences themselves if kept.
#[derive(Clone, Default)]
pub struct Endpoints {
    /// sequence length
    pub len: usize,
    /// phase-2 index of the `H_A` state each sequence ends in
    pub index: Vec<u32>,
    /// the sequences, `len` internal moves each, in the order of `index` (empty if not kept)
    pub paths: Vec<u8>,
}

impl Endpoints {
    /// The `k`-th sequence, if kept.
    pub fn path(&self, k: usize) -> Option<&[u8]> {
        self.paths.get(k * self.len..(k + 1) * self.len)
    }

    fn append(&mut self, o: &mut Endpoints) {
        self.index.append(&mut o.index);
        self.paths.append(&mut o.paths);
    }
}

/// Counters of one enumeration.
#[derive(Clone, Copy, Debug, Default)]
pub struct Stats {
    /// evaluated nodes (children generated)
    pub nodes: u64,
    /// leaves the tables put in `H_A` (each checked in full)
    pub leaves: u64,
    /// of those, the sequences that end in `H_A` as patterns
    pub in_h: u64,
    /// of those, the ones that end in `H_A` with the B-centres told apart (the lift's own)
    pub exact: u64,
}

impl Stats {
    pub fn add(&mut self, o: &Stats) {
        self.nodes += o.nodes;
        self.leaves += o.leaves;
        self.in_h += o.in_h;
        self.exact += o.exact;
    }
}

/// A node at the split depth: the root of one task.
#[derive(Clone)]
struct Task {
    x: Fc,
    ha: u8,
    hc: u8,
    g: usize,
    last: usize,
    path: [u8; MAX_LEN],
}

/// One worker's state during a depth-first search.
struct Walker {
    path: [u8; MAX_LEN],
    found: Endpoints,
    /// sequences that end in `H_A` only as patterns (kept with `keep`)
    other: Endpoints,
    keep: bool,
    stats: Stats,
    /// stop at the first sequence found
    first: bool,
    done: bool,
    /// collect the nodes at this depth as tasks instead of searching below them
    split: usize,
    tasks: Vec<Task>,
}

impl Walker {
    fn new(len: usize, keep: bool, first: bool, split: usize) -> Walker {
        Walker {
            path: [0; MAX_LEN],
            found: Endpoints {
                len,
                ..Endpoints::default()
            },
            other: Endpoints {
                len,
                ..Endpoints::default()
            },
            keep,
            stats: Stats::default(),
            first,
            done: false,
            split,
            tasks: vec![],
        }
    }
}

/// Phase-1 enumeration from one position (any element of the coset: they all have the same
/// phase-1 sequences).
pub struct Phase1<'a, I: AeIndex> {
    t: &'a Tables<I>,
    f: &'static FastTables,
    cube: Cube,
    /// the position with its B-centres told apart (a lift: see [`Phase1::new`])
    lift: Cube,
    attach: &'static [[u8; 12]],
    /// record only the sequences that end in `H_A` with the B-centres told apart (the coset of
    /// the lift; see the module documentation of [`super`]); otherwise every sequence that ends
    /// in `H_A` as a pattern
    labelled: bool,
    x: Fc,
    ha: u8,
    hc: u8,
}

impl<'a, I: AeIndex> Phase1<'a, I> {
    /// # Errors
    ///
    /// If the tables are of a kind the enumeration does not read: the tiny tier (whose
    /// A-centre/edge bound is not a walked distance) or the slice metric.
    /// With `lift` (the position with its B-centres told apart, labels as
    /// [`super::space::labelled_attachment`] numbers them) only the sequences that take the lift
    /// itself into `H_A` are recorded: the coset of the lift, whose elements share them. Without
    /// it, every sequence that takes the position into `H_A` as a pattern.
    pub fn new(t: &'a Tables<I>, cube: Cube, lift: Option<Cube>) -> Result<Phase1<'a, I>, String> {
        let attach = super::space::labelled_attachment()
            .ok_or("the B-centres of H_A do not follow the corners")?;
        if lift.is_some_and(|l| super::space::unlabel(&l) != cube) {
            return Err("the lift is not a lift of the position".into());
        }
        if I::DIRECT {
            return Err("the coset solver needs the small tier or larger, not tiny".into());
        }
        if t.moves.slice {
            return Err("the coset solver works in the face-turn metric only".into());
        }
        let f = fast();
        let x = f.of(&cube);
        let (ha, hc) = t.exact_fast(&x);
        Ok(Phase1 {
            t,
            f,
            cube,
            lift: lift.unwrap_or_else(|| super::space::lift(&cube)),
            attach,
            labelled: lift.is_some(),
            x,
            ha,
            hc,
        })
    }

    /// The position the sequences start from.
    pub fn cube(&self) -> &Cube {
        &self.cube
    }

    /// The position with its B-centres told apart ([`super::space::lift`]).
    pub fn lift(&self) -> &Cube {
        &self.lift
    }

    /// Whether only the lift's own sequences are recorded (see [`Phase1::new`]).
    pub fn labelled(&self) -> bool {
        self.labelled
    }

    /// A lower bound on the phase-1 depth: the larger exact table distance (and the joint
    /// table's bound, if loaded).
    pub fn lower_bound(&self) -> usize {
        let mut h = self.ha.max(self.hc);
        if let Some(cbe) = &self.t.cbe {
            let mut ci = [0usize; 4];
            cbe.prefetch_all(self.f, &self.x, &mut ci);
            h = h.max(cbe.bound(&ci));
        }
        usize::from(h)
    }

    /// Every sequence of `len` moves (`len < MAX_LEN`), in parallel on the current rayon pool.
    /// The tasks are searched in `chunks` groups; after each, `on_chunk` gets that group's
    /// sequences and counters and says whether to go on (`false` ends the enumeration, which
    /// is then incomplete). Returns whether it completed.
    pub fn level(
        &self,
        len: usize,
        keep: bool,
        chunks: usize,
        stop: &AtomicBool,
        mut on_chunk: impl FnMut(Endpoints, Endpoints, Stats) -> bool,
    ) -> bool {
        if len >= MAX_LEN {
            return false;
        }
        let split = if len >= 15 { 4 } else { 3 };
        if len <= split {
            let mut w = Walker::new(len, keep, false, usize::MAX);
            self.root(&mut w, len);
            on_chunk(w.found, w.other, w.stats);
            return true;
        }
        let mut w = Walker::new(len, keep, false, split);
        self.root(&mut w, len);
        let tasks = std::mem::take(&mut w.tasks);
        // (the prefix nodes' counts belong to the first chunk)
        let mut pre = Some(w.stats);
        let size = tasks.len().div_ceil(chunks.max(1)).max(1);
        let groups = tasks.len().div_ceil(size);
        for (gi, group) in tasks.chunks(size).enumerate() {
            let parts: Vec<(Endpoints, Endpoints, Stats)> = group
                .par_iter()
                .map(|task| {
                    let mut w = Walker::new(len, keep, false, usize::MAX);
                    if stop.load(Ordering::Relaxed) {
                        return (w.found, w.other, w.stats);
                    }
                    w.path = task.path;
                    self.dfs(&mut w, &task.x, task.ha, task.hc, task.g, len, task.last);
                    (w.found, w.other, w.stats)
                })
                .collect();
            let mut found = Endpoints {
                len,
                ..Endpoints::default()
            };
            let mut other = found.clone();
            let mut stats = pre.take().unwrap_or_default();
            for (mut e, mut o, s) in parts {
                found.append(&mut e);
                other.append(&mut o);
                stats.add(&s);
            }
            if stop.load(Ordering::Relaxed) {
                return false;
            }
            // (after the last group the length is complete whatever the answer)
            if !on_chunk(found, other, stats) && gi + 1 < groups {
                return false;
            }
        }
        true
    }

    /// The first sequence of `len` moves found (one thread), and the counters.
    pub fn first(&self, len: usize) -> (Option<Vec<usize>>, Stats) {
        if len >= MAX_LEN {
            return (None, Stats::default());
        }
        let mut w = Walker::new(len, true, true, usize::MAX);
        self.root(&mut w, len);
        let p = w
            .found
            .path(0)
            .map(|p| p.iter().map(|&m| usize::from(m)).collect());
        (p, w.stats)
    }

    /// The exact phase-1 depth (IDA* from the lower bound, one thread) and a shortest
    /// sequence, if it is at most `max_len`; and the counters.
    pub fn depth(&self, max_len: usize) -> (Option<(usize, Vec<usize>)>, Stats) {
        let mut stats = Stats::default();
        for len in self.lower_bound()..=max_len.min(MAX_LEN - 1) {
            let (p, s) = self.first(len);
            stats.add(&s);
            if let Some(p) = p {
                return (Some((len, p)), stats);
            }
        }
        (None, stats)
    }

    /// Search from the root for sequences of `len` moves.
    fn root(&self, w: &mut Walker, len: usize) {
        if len == 0 {
            // (the empty sequence: the position itself is in H_A)
            if self.ha == 0 && self.hc == 0 {
                self.leaf(w, 0);
            }
            return;
        }
        self.dfs(w, &self.x, self.ha, self.hc, 0, len, N_MOVES);
    }

    /// A leaf (depth `g`, AE and CB both solved): check membership in full and record it.
    fn leaf(&self, w: &mut Walker, g: usize) {
        w.stats.leaves += 1;
        let lab = w.path[..g]
            .iter()
            .fold(self.lift, |c, &m| c.apply(usize::from(m)));
        let Some(i) = self.t.p2.index(&super::space::unlabel(&lab)) else {
            return;
        };
        w.stats.in_h += 1;
        let exact = self
            .attach
            .get(usize::from(crate::coord::corners().id(&lab)))
            .is_some_and(|a| *a == lab.xb);
        w.stats.exact += u64::from(exact);
        let e = if exact || !self.labelled {
            &mut w.found
        } else {
            &mut w.other
        };
        e.index.push(i as u32);
        if w.keep {
            e.paths.extend_from_slice(&w.path[..g]);
        }
        if w.first && (exact || !self.labelled) {
            w.done = true;
        }
    }

    /// Depth-first search below a node at depth `g < len` with exact distances `ha`, `hc`
    /// (`last`: the move leading to it, `N_MOVES` at the root).
    #[allow(clippy::too_many_arguments)]
    fn dfs(&self, w: &mut Walker, x: &Fc, ha: u8, hc: u8, g: usize, len: usize, last: usize) {
        if g == w.split {
            let mut task = Task {
                x: *x,
                ha,
                hc,
                g,
                last,
                path: [0; MAX_LEN],
            };
            task.path[..g].copy_from_slice(&w.path[..g]);
            w.tasks.push(task);
            return;
        }
        let (t, f) = (self.t, self.f);
        let budget = len - g - 1;
        let mut moves = ALLOWED[last.min(N_MOVES)];
        if budget == 0 {
            moves &= B_MASK;
        }
        // the children, then all their table lookups issued before any is read
        let fw = f.flip_word(&masks(&x.e));
        let mut ys = [Fc::ZERO; N_MOVES];
        let mut ms = [0u8; N_MOVES];
        let mut n = 0;
        while moves != 0 {
            let m = moves.trailing_zeros() as usize;
            moves &= moves - 1;
            ys[n % N_MOVES] = f.apply_w(x, fw, m);
            ms[n % N_MOVES] = m as u8;
            n += 1;
        }
        let n = n.min(N_MOVES);
        w.stats.nodes += n as u64;
        let mut ai = [0usize; N_MOVES];
        t.ae_ix
            .index_batch_fast_pf::<N_MOVES>(f, &ys[..n], &mut ai[..n], &t.ae.data);
        let mut ci = [0usize; N_MOVES];
        for i in 0..n {
            ci[i] = CbTable::index_fast(f, &ys[i]);
            prefetch_at(&t.cb.tab.data, ci[i] >> 5);
        }
        let mut sel = [0usize; N_MOVES];
        let mut hs = [(0u8, 0u8); N_MOVES];
        let mut k = 0;
        for i in 0..n {
            let a = decode(ha, t.ae.get(ai[i]));
            let c = decode(hc, t.cb.tab.get(ci[i]));
            if usize::from(a) <= budget && usize::from(c) <= budget {
                hs[i] = (a, c);
                sel[k % N_MOVES] = i;
                k += 1;
            }
        }
        if budget <= CBE_BUDGET && k > 0 {
            if let Some(cbe) = &t.cbe {
                let mut qi = [[0usize; 4]; N_MOVES];
                for &i in &sel[..k] {
                    cbe.prefetch_all(f, &ys[i % N_MOVES], &mut qi[i % N_MOVES]);
                }
                let mut k2 = 0;
                for j in 0..k {
                    let i = sel[j];
                    if usize::from(cbe.bound(&qi[i % N_MOVES])) <= budget {
                        sel[k2 % N_MOVES] = i;
                        k2 += 1;
                    }
                }
                k = k2;
            }
        }
        for &i in &sel[..k] {
            let i = i % N_MOVES;
            w.path[g % MAX_LEN] = ms[i];
            if g + 1 == len {
                self.leaf(w, len);
            } else {
                let (a, c) = hs[i];
                self.dfs(w, &ys[i], a, c, g + 1, len, usize::from(ms[i]));
            }
            if w.done {
                return;
            }
        }
    }
}
