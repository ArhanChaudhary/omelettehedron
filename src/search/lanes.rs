//! Lane-interleaved phase-1 search, after egg-solver's `lanes.rs`.
//!
//! Almost every table lookup misses the cache. Instead of waiting, each worker thread runs
//! several independent "lanes": each lane is an explicit-stack DFS of one task that generates
//! one child at a time and stops whenever it needs a table value: it computes the index,
//! prefetches the entry and yields. Lanes are visited round-robin; a visit first resolves the
//! lookup issued on the previous visit, then advances the lane to its next lookup, so every
//! miss gets a whole round of slack while about one miss per lane is in flight.
//!
//! A child goes through up to four checks, each one lookup (or a few independent ones):
//! AE and CB (AE first when 7+ moves remain, where CB cannot prune), then the joint CBE
//! table (8 or fewer moves left, large tiers), then the class-B bound against the best
//! solution (once one is known). The searched tree is the same as the recursive search's.

use super::*;

/// Lanes per worker thread unless overridden with `FTO_LANES`.
pub(super) const DEFAULT_LANES: usize = 24;
const MAXD: usize = 40;

#[derive(Clone, Copy)]
struct Frame {
    x: Fc,
    xb: Fc,
    em: [u16; 4],
    emb: [u16; 4],
    h_ae: u8,
    h_cb: u8,
    hb_ae: u8,
    hb_cb: u8,
    /// moves still to try
    todo: u16,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Check {
    Ae,
    Cb,
    Cbe,
    B,
}

struct Lane {
    task: usize,
    v: usize,
    use_b: bool,
    /// length of the task prefix
    base: usize,
    /// stack top (0 = task root)
    depth: usize,
    stack: Box<[Frame; MAXD]>,
    path: [usize; 64],
    stat: f32,
    // pending child of the top frame
    m: usize,
    y: Fc,
    ha: u8,
    hc: u8,
    yb: Fc,
    hba: u8,
    hbc: u8,
    checks: [Check; 4],
    nchecks: usize,
    ci: usize,
    idx: [usize; 4],
}

const EMPTY_FC: Fc = Fc {
    a: [0; 16],
    b: [0; 16],
    e: [0; 16],
    c: 0,
    ebits: 0,
};
const EMPTY_FRAME: Frame = Frame {
    x: EMPTY_FC,
    xb: EMPTY_FC,
    em: [0; 4],
    emb: [0; 4],
    h_ae: 0,
    h_cb: 0,
    hb_ae: 0,
    hb_cb: 0,
    todo: 0,
};

/// Moves allowed after `last` (usize::MAX at the root), only non-A moves if no move is left
/// after this one (a phase-1 solution ends with a move outside H_A).
#[inline(always)]
fn allowed_mask(last: usize, budget_after: usize) -> u16 {
    let mut mask = 0u16;
    for m in 0..N_MOVES {
        if allowed_after(last, m) && !(budget_after == 0 && is_a_move(m)) {
            mask |= 1 << m;
        }
    }
    mask
}

impl Lane {
    fn new() -> Lane {
        Lane {
            task: 0,
            v: 0,
            use_b: false,
            base: 0,
            depth: 0,
            stack: Box::new([EMPTY_FRAME; MAXD]),
            path: [0; 64],
            stat: 0.0,
            m: 0,
            y: EMPTY_FC,
            ha: 0,
            hc: 0,
            yb: EMPTY_FC,
            hba: 0,
            hbc: 0,
            checks: [Check::Ae; 4],
            nchecks: 0,
            ci: 0,
            idx: [0; 4],
        }
    }
}

pub(super) struct Iter<'s, 'a, I: AeIndex> {
    pub sh: &'s Shared<'a, I>,
    pub tasks: &'s [Task],
    pub next: &'s std::sync::atomic::AtomicUsize,
    pub stats: &'s [std::sync::atomic::AtomicU32],
    pub bound: usize,
    pub count_mode: bool,
}

impl<'s, 'a, I: AeIndex> Iter<'s, 'a, I> {
    /// One worker thread: run `n_lanes` lanes until the iteration's tasks are exhausted.
    pub fn worker(&self, n_lanes: usize) {
        let mut lanes: Vec<Lane> = (0..n_lanes).map(|_| Lane::new()).collect();
        let mut nodes = 0u64;
        let mut reported = 0u64;
        let mut n_active = 0;
        for i in 0..lanes.len() {
            if self.start_lane(&mut lanes[i], &mut nodes) {
                lanes.swap(i, n_active);
                n_active += 1;
            }
        }
        let mut checked = 0u64;
        while n_active > 0 {
            let mut i = 0;
            while i < n_active {
                if self.step(&mut lanes[i], &mut nodes) {
                    i += 1;
                } else {
                    n_active -= 1;
                    lanes.swap(i, n_active);
                }
            }
            if nodes - checked >= CHECK_EVERY {
                checked = nodes;
                self.sh.nodes.fetch_add(nodes - reported, Ordering::Relaxed);
                reported = nodes;
                if self.sh.time_up() {
                    break;
                }
            }
            if self.sh.stop.load(Ordering::Relaxed) {
                break;
            }
        }
        self.sh.nodes.fetch_add(nodes - reported, Ordering::Relaxed);
    }

    /// Give `lane` its next task and advance it to its first lookup; false if none is left.
    fn start_lane(&self, lane: &mut Lane, nodes: &mut u64) -> bool {
        loop {
            let i = self.next.fetch_add(1, Ordering::Relaxed);
            if i >= self.tasks.len() || self.sh.stop.load(Ordering::Relaxed) {
                return false;
            }
            let t = &self.tasks[i];
            lane.task = i;
            lane.v = t.v;
            lane.use_b = !self.count_mode && self.sh.best_len.load(Ordering::Relaxed) < usize::MAX;
            lane.base = t.path.len();
            lane.depth = 0;
            lane.stat = 0.0;
            lane.path[..lane.base].copy_from_slice(&t.path);
            if lane.base == self.bound {
                // the task itself is a phase-1 leaf candidate
                if t.h_ae == 0
                    && t.h_cb == 0
                    && (lane.base == 0 || !is_a_move(t.path[lane.base - 1]))
                {
                    self.sh.leaf(lane.v, &lane.path[..lane.base]);
                }
                self.finish_task(lane);
                continue;
            }
            let last = t.path.last().copied().unwrap_or(usize::MAX);
            let budget_after = self.bound - lane.base - 1;
            lane.stack[0] = Frame {
                x: t.x,
                xb: t.xb,
                em: masks(&t.x.e),
                emb: masks(&t.xb.e),
                h_ae: t.h_ae,
                h_cb: t.h_cb,
                hb_ae: t.hb_ae,
                hb_cb: t.hb_cb,
                todo: allowed_mask(last, budget_after),
            };
            if self.next_child(lane, nodes) {
                return true;
            }
        }
    }

    fn finish_task(&self, lane: &mut Lane) {
        self.stats[lane.task].store(lane.stat.to_bits(), Ordering::Relaxed);
    }

    /// Generate the next child of the top frame (popping finished frames) and issue its first
    /// lookup. False if the lane's task is finished and no new task is available.
    #[inline(always)]
    fn next_child(&self, lane: &mut Lane, nodes: &mut u64) -> bool {
        loop {
            let todo = lane.stack[lane.depth].todo;
            if todo == 0 {
                if lane.depth == 0 {
                    self.finish_task(lane);
                    return self.start_lane_inner(lane, nodes);
                }
                lane.depth -= 1;
                continue;
            }
            let m = todo.trailing_zeros() as usize;
            let top = &mut lane.stack[lane.depth];
            top.todo = todo & (todo - 1);
            *nodes += 1;
            let f = fast();
            lane.m = m;
            lane.y = f.apply(&top.x, &top.em, m);
            let g = lane.base + lane.depth;
            let budget = self.bound - g - 1;
            let (first, second) = if budget >= 7 {
                (Check::Ae, Check::Cb)
            } else {
                (Check::Cb, Check::Ae)
            };
            lane.checks[0] = first;
            lane.checks[1] = second;
            let mut n = 2;
            if budget <= 8 && self.sh.t.cbe.is_some() {
                lane.checks[n] = Check::Cbe;
                n += 1;
            }
            if lane.use_b {
                lane.checks[n] = Check::B;
                n += 1;
            }
            lane.nchecks = n;
            lane.ci = 0;
            self.issue(lane);
            return true;
        }
    }

    /// `start_lane` without the recursion through `next_child` (tail position).
    #[inline(never)]
    fn start_lane_inner(&self, lane: &mut Lane, nodes: &mut u64) -> bool {
        self.start_lane(lane, nodes)
    }

    /// Compute the index of the current check and prefetch it.
    #[inline(always)]
    fn issue(&self, lane: &mut Lane) {
        let t = self.sh.t;
        let f = fast();
        match lane.checks[lane.ci] {
            Check::Ae => {
                lane.idx[0] = t.ae_ix.index_fast(f, &lane.y, f.a_rank(&lane.y));
                prefetch(&t.ae.data[lane.idx[0] >> 5]);
            }
            Check::Cb => {
                lane.idx[0] = CbTable::index_fast(f, &lane.y);
                prefetch(&t.cb.tab.data[lane.idx[0] >> 5]);
            }
            Check::Cbe => {
                if let Some(cbe) = &t.cbe {
                    let mut ci = [0usize; 4];
                    cbe.prefetch_all(f, &lane.y, &mut ci);
                    lane.idx = ci;
                }
            }
            Check::B => {
                let top = &lane.stack[lane.depth];
                lane.yb = f.apply(&top.xb, &top.emb, self.sh.mb[lane.m]);
                lane.idx[0] = t.ae_ix.index_fast(f, &lane.yb, f.a_rank(&lane.yb));
                lane.idx[1] = CbTable::index_fast(f, &lane.yb);
                prefetch(&t.ae.data[lane.idx[0] >> 5]);
                prefetch(&t.cb.tab.data[lane.idx[1] >> 5]);
            }
        }
    }

    /// Resolve the pending lookup and advance to the next one. False if the lane is done.
    #[inline(always)]
    fn step(&self, lane: &mut Lane, nodes: &mut u64) -> bool {
        let t = self.sh.t;
        let g = lane.base + lane.depth;
        let budget = self.bound - g - 1;
        let top = &lane.stack[lane.depth];
        let pass = match lane.checks[lane.ci] {
            Check::Ae => {
                lane.ha = decode(top.h_ae, t.ae.get(lane.idx[0]));
                lane.ha as usize <= budget
            }
            Check::Cb => {
                lane.hc = decode(top.h_cb, t.cb.tab.get(lane.idx[0]));
                lane.hc as usize <= budget
            }
            Check::Cbe => match &t.cbe {
                Some(cbe) => cbe.bound(&lane.idx) as usize <= budget,
                None => true,
            },
            Check::B => {
                lane.hba = decode(top.hb_ae, t.ae.get(lane.idx[0]));
                lane.hbc = decode(top.hb_cb, t.cb.tab.get(lane.idx[1]));
                let best = self.sh.best_len.load(Ordering::Relaxed);
                g + 1 + (lane.hba.max(lane.hbc) as usize) < best
            }
        };
        if !pass {
            return self.next_child(lane, nodes);
        }
        lane.ci += 1;
        if lane.ci < lane.nchecks {
            self.issue(lane);
            return true;
        }
        // the child passed every check
        lane.stat += STAT_W[lane.ha.max(lane.hc) as usize];
        lane.path[g] = lane.m;
        if g + 1 == self.bound {
            // (bounds are 0 here and the move is outside H_A)
            self.sh.leaf(lane.v, &lane.path[..g + 1]);
        } else {
            let budget_after = self.bound - g - 2;
            let child = Frame {
                x: lane.y,
                xb: lane.yb,
                em: masks(&lane.y.e),
                emb: if lane.use_b {
                    masks(&lane.yb.e)
                } else {
                    [0; 4]
                },
                h_ae: lane.ha,
                h_cb: lane.hc,
                hb_ae: lane.hba,
                hb_cb: lane.hbc,
                todo: allowed_mask(lane.m, budget_after),
            };
            lane.depth += 1;
            lane.stack[lane.depth] = child;
        }
        self.next_child(lane, nodes)
    }
}
