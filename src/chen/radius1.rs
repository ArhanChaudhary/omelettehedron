//! The radius-1 phase-2 space: the positions whose relaxed phase-1 distance ([`super::relaxed`])
//! is at most 1, that is, whose A-centres and edges form a configuration of the block split's
//! phase-2 space ([`super::space`]: distance 0) or are one face turn from one (distance 1).
//!
//! The distance depends on the A-centres and edges alone and is the same in every frame (the
//! configurations are closed under the 12 rotations that keep the face classes, and a rotation
//! takes moves to moves), so the space is a set of pattern-view configurations with every one of
//! the 11,520 corner states and 369,600 B-centre colourings, and a move is allowed exactly when
//! it stays in the set. The breadth-first search ([`super::p2bfs`]) indexes it like the block
//! split's space, with these configurations in place of the 37,260.
//!
//! In the fixed frame the set is the 447,120 configurations and their neighbours under the 16
//! moves: [`N_CFG1`] configurations in 12 frames of [`N_PAT1`] (37,260 at distance 0, 278,640 at
//! distance 1). A configuration's frame follows its parent's by the move's frame change (the
//! same rotation for every configuration of a frame, read off the block split's space), and the
//! build checks that every configuration gets one frame however it is reached, that frame 0 is
//! closed under the 24 symmetries, and that every move can be undone.

use super::space::{key, patterns, space, sym_maps, Patterns, QuickMap, SymMaps, NONE};
use crate::cube::{inverse_move, Cube, N_MOVES};
use crate::geometry::{class_preserving_group, rotation_group, transpose, Mat};
use std::sync::OnceLock;

/// Configurations of the radius-1 space in the fixed frame.
pub const N_CFG1: usize = 3_790_800;
/// Its pattern-view configurations (one frame).
pub const N_PAT1: usize = 315_900;
/// Frames: the rotations that keep the face classes.
const N_FRAMES: usize = 12;

/// The set cannot be built: a check of its construction failed (a bug).
#[derive(Debug, Clone, thiserror::Error)]
#[error("the radius-1 configuration set: {0}")]
pub struct Radius1Error(String);

/// The radius-1 configuration set with its move tables.
pub struct Radius1 {
    /// configurations in the fixed frame: the block split's 447,120 first (its numbering), then
    /// those at distance 1
    pub cfgs: Vec<([u8; 12], [u8; 12])>,
    /// configuration id of each key ([`key`])
    ids: QuickMap<u128, u32>,
    /// configuration after each move, or [`NONE`] where the move leaves the set
    pub next: Vec<[u32; N_MOVES]>,
    /// the frame of each configuration (index into [`super::space::Space::rot_xb`])
    pub rot: Vec<u8>,
    /// the pattern-view configurations and their classes under the 24 symmetries
    pub pat: Patterns,
    /// per pattern-view configuration: its relaxed phase-1 distance, 0 or 1
    pub d1: Vec<u8>,
    /// per pattern-view configuration: whether it is in `H'` (every hexagon one colour)
    pub hprime: Vec<bool>,
}

impl Radius1 {
    /// Configuration id of a position, if its configuration is in the set.
    pub fn id(&self, cube: &Cube) -> Option<u32> {
        self.ids.get(&key(&cube.xa, &cube.ep)).copied()
    }

    /// The relaxed phase-1 distance of a fixed-frame configuration: 0 for the block split's.
    #[inline]
    pub fn cfg_d1(&self, cfg: u32) -> u8 {
        u8::from(cfg as usize >= space().cfgs.len())
    }
}

fn build() -> Result<Radius1, Radius1Error> {
    let err = |what: String| Radius1Error(what);
    let s = space();
    let s0 = Cube::solved();
    let n0 = s.cfgs.len();
    let mut cfgs = s.cfgs.clone();
    let mut ids: QuickMap<u128, u32> = QuickMap::default();
    for (i, &(xa, ep)) in cfgs.iter().enumerate() {
        ids.insert(key(&xa, &ep), i as u32);
    }
    // the frame change of each move from each frame, from the block split's space
    let mut after = [[u8::MAX; N_MOVES]; N_FRAMES];
    for c in 0..n0 {
        let r = s.rot[c] as usize;
        for m in 0..N_MOVES {
            let n = s.next[c][m];
            if n == NONE {
                continue;
            }
            let r2 = s.rot[n as usize];
            let e = after
                .get_mut(r)
                .and_then(|row| row.get_mut(m))
                .ok_or_else(|| err(format!("frame {r} out of range")))?;
            if *e == u8::MAX {
                *e = r2;
            } else if *e != r2 {
                return Err(err(format!(
                    "move {m} from frame {r} leads to frames {} and {r2}",
                    *e
                )));
            }
        }
    }
    if after.iter().flatten().any(|&a| a == u8::MAX) {
        return Err(err("a frame change is not seen".into()));
    }
    let mut rot = s.rot.clone();
    let mut next: Vec<[u32; N_MOVES]> = vec![[NONE; N_MOVES]; n0];
    // the neighbours of the block split's configurations: at distance 1 where new
    for c in 0..n0 {
        let (xa, ep) = cfgs[c];
        let cube = Cube { ep, xa, ..s0 };
        for m in 0..N_MOVES {
            let n = cube.apply(m);
            let k = key(&n.xa, &n.ep);
            let want = after[rot[c] as usize % N_FRAMES][m];
            let id = if let Some(&id) = ids.get(&k) {
                id
            } else {
                let id = cfgs.len() as u32;
                ids.insert(k, id);
                cfgs.push((n.xa, n.ep));
                rot.push(want);
                id
            };
            if rot[id as usize] != want {
                return Err(err(format!(
                    "configuration {id} is in frame {} and in frame {want}",
                    rot[id as usize]
                )));
            }
            if s.next[c][m] != NONE && s.next[c][m] != id {
                return Err(err(format!("move {m} of configuration {c} disagrees")));
            }
            next[c][m] = id;
        }
    }
    // the moves of the configurations at distance 1: those that stay in the set
    for c in n0..cfgs.len() {
        let (xa, ep) = cfgs[c];
        let cube = Cube { ep, xa, ..s0 };
        let mut row = [NONE; N_MOVES];
        for m in 0..N_MOVES {
            let n = cube.apply(m);
            if let Some(&id) = ids.get(&key(&n.xa, &n.ep)) {
                let want = after[rot[c] as usize % N_FRAMES][m];
                if rot[id as usize] != want {
                    return Err(err(format!(
                        "configuration {id} is in frame {} and in frame {want}",
                        rot[id as usize]
                    )));
                }
                row[m] = id;
            }
        }
        next.push(row);
    }
    if cfgs.len() != N_CFG1 {
        return Err(err(format!("{} configurations, not {N_CFG1}", cfgs.len())));
    }
    let frame0: usize = rot.iter().map(|&r| usize::from(r == 0)).sum();
    if frame0 != N_PAT1 {
        return Err(err(format!(
            "{frame0} configurations in frame 0, not {N_PAT1}"
        )));
    }
    for (c, row) in next.iter().enumerate() {
        for (m, &n) in row.iter().enumerate() {
            if n != NONE && next[n as usize][inverse_move(m)] != c as u32 {
                return Err(err(format!(
                    "move {m} of configuration {c} cannot be undone"
                )));
            }
        }
    }
    // the pattern view
    let rots: Vec<Mat> = rotation_group().iter().map(transpose).collect();
    if rots.len() != N_FRAMES {
        return Err(err(format!("{} rotations, not {N_FRAMES}", rots.len())));
    }
    let pat = patterns(&cfgs, &ids, &next, &rot, &rots, &s.rot_xb);
    if pat.cfg.len() != N_PAT1 {
        return Err(err(format!(
            "{} pattern-view configurations, not {N_PAT1}",
            pat.cfg.len()
        )));
    }
    // every configuration turned by its frame is the frame-0 configuration recorded for it
    let turn: Vec<SymMaps> = rots.iter().map(sym_maps).collect();
    for (c, &(xa, ep)) in cfgs.iter().enumerate() {
        let img = turn[rot[c] as usize % N_FRAMES].apply(&xa, &ep, false);
        let p = pat.of_cfg[c] as usize;
        if pat.cfg.get(p).map(|&q| cfgs[q as usize]) != Some(img) {
            // (diagnostics: the frames that do turn it into the set, and its parents)
            let into: Vec<(usize, Option<u8>)> = (0..N_FRAMES)
                .map(|r| {
                    let (xa2, ep2) = turn[r].apply(&xa, &ep, false);
                    (r, ids.get(&key(&xa2, &ep2)).map(|&id| rot[id as usize]))
                })
                .collect();
            let frames = &rot;
            let parents: Vec<(usize, usize, u8)> = next
                .iter()
                .enumerate()
                .flat_map(|(c2, row)| {
                    row.iter()
                        .enumerate()
                        .filter(|&(_, &n)| n as usize == c)
                        .map(move |(m, _)| (c2, m, frames[c2]))
                })
                .collect();
            return Err(err(format!(
                "configuration {c} (frame {}, distance {}) turned is not in frame 0; turned by \
                 each frame's rotation it lands in the set with these frames: {into:?}; its \
                 parents (configuration, move, frame): {parents:?}",
                rot[c],
                u8::from(c >= n0)
            )));
        }
    }
    // frame 0 is closed under the 24 symmetries, and the classes map onto their
    // representatives
    let syms: Vec<SymMaps> = class_preserving_group().iter().map(sym_maps).collect();
    let frame0_id = |img: ([u8; 12], [u8; 12])| -> Option<u32> {
        ids.get(&key(&img.0, &img.1))
            .copied()
            .filter(|&id| rot[id as usize] == 0)
    };
    for (p, &c) in pat.cfg.iter().enumerate() {
        let (xa, ep) = cfgs[c as usize];
        for sym in &syms {
            if frame0_id(sym.apply(&xa, &ep, true)).is_none() {
                return Err(err(format!(
                    "frame 0 is not closed under the symmetries (configuration {c})"
                )));
            }
        }
        let (cl, k) = pat.class[p];
        let r = pat.rep.get(cl as usize).copied().unwrap_or(u32::MAX) as usize;
        let img = syms
            .get(k as usize)
            .map(|sym| sym.apply(&xa, &ep, true))
            .and_then(frame0_id);
        if img != pat.cfg.get(r).copied() {
            return Err(err(format!(
                "pattern-view configuration {p} is not mapped onto its representative"
            )));
        }
    }
    let d1: Vec<u8> = pat
        .cfg
        .iter()
        .map(|&c| u8::from(c as usize >= n0))
        .collect();
    let n_d0: usize = d1.iter().map(|&d| usize::from(d == 0)).sum();
    if n_d0 != s.pat.cfg.len() {
        return Err(err(format!(
            "{n_d0} pattern-view configurations at distance 0, not {}",
            s.pat.cfg.len()
        )));
    }
    let hprime = pat
        .cfg
        .iter()
        .map(|&c| (c as usize) < n0 && s.pshape[c as usize] == 0)
        .collect();
    Ok(Radius1 {
        cfgs,
        ids,
        next,
        rot,
        pat,
        d1,
        hprime,
    })
}

/// The radius-1 configuration set, built once (a few seconds).
///
/// # Errors
///
/// [`Radius1Error`] if a check of the construction fails.
pub fn radius1() -> Result<&'static Radius1, Radius1Error> {
    static S: OnceLock<Result<Radius1, Radius1Error>> = OnceLock::new();
    S.get_or_init(build).as_ref().map_err(Clone::clone)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chen::relaxed::RelaxedAe;
    use crate::fast::fast;
    use crate::tables::AeSym;
    use std::collections::HashMap;

    /// The set builds with its checks; the distances 0 and 1 are the relaxed table's (where
    /// a table is at hand, else by the target set), and every configuration two moves from
    /// the block split's set but not in this one is at distance 2.
    #[test]
    fn set_matches_the_relaxed_distances() {
        let r = radius1().unwrap();
        let s = space();
        assert_eq!(
            r.d1.iter().map(|&d| usize::from(d == 1)).sum::<usize>(),
            278_640
        );
        assert_eq!(r.hprime.iter().filter(|&&h| h).count(), 81);
        // per frame, 315,900 configurations
        let mut per = HashMap::new();
        for &f in &r.rot {
            *per.entry(f).or_insert(0usize) += 1;
        }
        assert_eq!(per.len(), 12);
        assert!(per.values().all(|&n| n == N_PAT1));
        let s0 = Cube::solved();
        let dir = std::env::var("FTO_CHEN_TABLES")
            .ok()
            .map(std::path::PathBuf::from);
        let ix = AeSym::new(true);
        let relaxed = dir.and_then(|d| RelaxedAe::load(&d, &ix, false).ok());
        let f = fast();
        let mut rng = <rand_xoshiro::Xoshiro256PlusPlus as rand::SeedableRng>::seed_from_u64(5);
        for _ in 0..3000 {
            let c = rand::Rng::gen_range(&mut rng, 0..r.cfgs.len());
            let (xa, ep) = r.cfgs[c];
            let cube = Cube { ep, xa, ..s0 };
            let d = r.cfg_d1(c as u32);
            assert_eq!(s.id(&cube).is_some(), d == 0);
            if let Some(t) = &relaxed {
                assert_eq!(t.exact(&ix, &f.of(&cube)), Some(d), "configuration {c}");
            }
            // a move out of the set lands at distance 2
            for m in 0..N_MOVES {
                if r.next[c][m] == NONE {
                    let n = cube.apply(m);
                    assert!(r.id(&n).is_none());
                    if let Some(t) = &relaxed {
                        assert_eq!(t.exact(&ix, &f.of(&n)), Some(2));
                    }
                }
            }
        }
    }
}
