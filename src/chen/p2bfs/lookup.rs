//! Positions looked up in the layers on disk: the slot of a fixed-frame position in an
//! instance of either configuration set, and whether a layer's file holds it (one chunk read),
//! for checks against sampled positions while a search runs and its layers come and go.

use super::codec::{ClassReader, CHUNK};
use super::engine::{canonical_slot, BfsError, Layout};
use super::instance::{CornerView, Instance, N_SYM};
use crate::chen::radius1::radius1;
use crate::chen::space::space;
use crate::coord::{corners, mset, N_CORNER};
use crate::cube::Cube;
use crate::fast::shuffle;
use crate::pattern;

/// The slot maps of an instance.
pub struct Lookup<'a> {
    inst: &'a Instance,
    /// the instance's corner state of each pattern-view corner state
    pc_to_p: Vec<u16>,
}

impl<'a> Lookup<'a> {
    /// The maps of `inst` (the full corner view).
    pub fn new(inst: &'a Instance) -> Option<Lookup<'a>> {
        if inst.kind.corners != CornerView::Full {
            return None;
        }
        let mut pc_to_p = vec![0u16; inst.n_p];
        for (p, &pc) in inst.p_to_pc.iter().enumerate() {
            *pc_to_p.get_mut(pc as usize)? = p as u16;
        }
        Some(Lookup { inst, pc_to_p })
    }

    /// The slot of a fixed-frame position: its configuration's class, and its B-centres and
    /// corners seen in the pattern view through the configuration's frame, conjugated to the
    /// class representative (one of the images under the stabiliser; a layer holds all of
    /// them). `None` if the configuration is not in the instance's set.
    pub fn slot_of(&self, cube: &Cube) -> Option<(u32, u32, u32)> {
        let s = space();
        let inst = self.inst;
        // the configuration's frame, and its pattern-view configuration's class and symmetry
        let (r, (cls, k)) = if inst.kind.radius == 0 {
            let cfg = s.id(cube)? as usize;
            (s.rot[cfg] as usize, s.pat.class[s.pat.of_cfg[cfg] as usize])
        } else {
            let r1 = radius1().ok()?;
            let cfg = r1.id(cube)? as usize;
            (
                r1.rot[cfg] as usize,
                r1.pat.class[r1.pat.of_cfg[cfg] as usize],
            )
        };
        let k = k as usize;
        let pc = *s
            .corners
            .of_fixed
            .get(r * N_CORNER + corners().id(cube) as usize)? as usize;
        let p = *inst.ps.get(*self.pc_to_p.get(pc)? as usize * N_SYM + k)?;
        let mut xb = [0u8; 16];
        xb[..12].copy_from_slice(&cube.xb);
        let xb = shuffle(&xb, s.rot_xb.get(r)?);
        let x = mset().rank(&std::array::from_fn(|i| xb[i]));
        let x = *inst.xs.get(k * inst.n_x + x as usize)?;
        Some((u32::from(cls), x, u32::from(p)))
    }

    /// The pattern of slot `(c, x, p)`: the position (the class's representative
    /// configuration, the corners and B-centres as the slot has them) as held in the pattern
    /// view; `None` for a slot outside the instance.
    pub fn pattern_of(&self, (c, x, p): (u32, u32, u32)) -> Option<String> {
        let s = space();
        let inst = self.inst;
        let q = *inst.rep.get(c as usize)? as usize;
        let (xa, ep) = if inst.kind.radius == 0 {
            *s.cfgs.get(*s.pat.cfg.get(q)? as usize)?
        } else {
            let r1 = radius1().ok()?;
            *r1.cfgs.get(*r1.pat.cfg.get(q)? as usize)?
        };
        let (cp, co) = *s
            .corners
            .states
            .get(*inst.p_to_pc.get(p as usize)? as usize)?;
        if x as usize >= inst.n_x {
            return None;
        }
        let xb = mset().unrank(x);
        let cube = Cube { cp, co, ep, xa, xb };
        Some(pattern::format(&pattern::of_cube(&cube)))
    }

    /// The canonical slot of a pattern's position: the smallest image under the class's
    /// stabiliser (`None` if the pattern does not parse or is outside the instance).
    pub fn canonical_slot_of(&self, pat: &str) -> Option<(u32, u32, u32)> {
        let cube = pattern::position(&pattern::parse(pat).ok()?).ok()?;
        let (c, x, p) = self.slot_of(&cube)?;
        let (x, p) = canonical_slot(self.inst, c as usize, x as usize, p as usize);
        Some((c, x as u32, p as u32))
    }

    /// Whether layer `d` holds slot `(c, x, p)`: reads the chunk of its block (a missing class
    /// file, of a layer being made or removed, holds nothing).
    ///
    /// # Errors
    ///
    /// If the class's file cannot be read.
    pub fn layer_has(
        &self,
        lay: &Layout,
        d: usize,
        (c, x, p): (u32, u32, u32),
    ) -> Result<bool, BfsError> {
        let inst = self.inst;
        let path = lay.file(d, c as usize);
        let err = |source: std::io::Error| BfsError::Io {
            path: path.display().to_string(),
            source,
        };
        let rd = match ClassReader::open(&path, inst.n_x, inst.n_p) {
            Ok(rd) => rd,
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(false),
            Err(e) => return Err(err(e)),
        };
        let mut buf = vec![];
        let mut bits: Vec<u16> = vec![];
        let mut found = false;
        rd.for_chunk(x as usize / CHUNK, &mut buf, false, |xx, b| {
            if xx == x as usize {
                bits.clear();
                b.bits(&mut bits, inst.n_p);
                found = bits.contains(&(p as u16));
            }
        })
        .map_err(err)?;
        Ok(found)
    }
}
