//! Geometric sticker model of the face-turning octahedron.
//!
//! The octahedron is `|x|+|y|+|z| <= 1`. Face `f` is the octant with outward normal
//! [`FACES`]`[f]`. Every face is cut into 9 triangles; sticker centroids are stored scaled by
//! 9 so that all coordinates are small integers. A face turn rotates every sticker whose
//! centroid satisfies `n·p > 1/3` (i.e. `> 3` after scaling).
//!
//! Notation (puzzle held with a corner pointing at you): the front corner is `+x`, the four
//! front faces are U (top), R (right), F (bottom) and L (left); the back faces are B (top),
//! BR, D (bottom) and BL.

pub type V3 = [i32; 3];
pub type Mat = [[i32; 3]; 3];

/// Outward normals of the 8 faces.
pub const FACES: [V3; 8] = [
    [1, 1, 1],
    [1, 1, -1],
    [1, -1, 1],
    [1, -1, -1],
    [-1, 1, 1],
    [-1, 1, -1],
    [-1, -1, 1],
    [-1, -1, -1],
];

/// Standard FTO names of the faces, indexed like [`FACES`].
pub const FACE_NAMES: [&str; 8] = ["U", "R", "L", "F", "B", "BR", "BL", "D"];

/// Class-A faces (U, F, BR, BL): product of normal coordinates is +1.
pub const A_FACES: [usize; 4] = [0, 3, 5, 6];
/// Class-B faces (R, L, B, D).
pub const B_FACES: [usize; 4] = [1, 2, 4, 7];

/// The six vertices; vertex 0 (`+x`, the front corner) is the internally fixed corner.
pub const VERTS: [V3; 6] = [
    [1, 0, 0],
    [-1, 0, 0],
    [0, 1, 0],
    [0, -1, 0],
    [0, 0, 1],
    [0, 0, -1],
];

#[inline]
pub fn opposite(f: usize) -> usize {
    7 - f
}

#[inline]
pub fn is_a_face(f: usize) -> bool {
    let n = FACES[f];
    n[0] * n[1] * n[2] == 1
}

pub fn dot(a: V3, b: V3) -> i32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub fn cross(a: V3, b: V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

pub fn mat_vec(m: &Mat, v: V3) -> V3 {
    [dot(m[0], v), dot(m[1], v), dot(m[2], v)]
}

pub fn mat_mul(a: &Mat, b: &Mat) -> Mat {
    let mut r = [[0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = (0..3).map(|k| a[i][k] * b[k][j]).sum();
        }
    }
    r
}

pub fn transpose(a: &Mat) -> Mat {
    let mut r = [[0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = a[j][i];
        }
    }
    r
}

pub fn det(m: &Mat) -> i32 {
    dot(m[0], cross(m[1], m[2]))
}

pub const IDENTITY: Mat = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

/// Index of the face whose normal is `n` (n must be a ±1 vector).
pub fn face_of_normal(n: V3) -> usize {
    FACES
        .iter()
        .position(|&f| f == n)
        .expect("not a face normal")
}

/// Index of the face containing surface point `p` (no zero coordinates).
pub fn face_of_point(p: V3) -> usize {
    face_of_normal([p[0].signum(), p[1].signum(), p[2].signum()])
}

/// Rotation by 120° about the normal of face `f`, clockwise when looking at the face from
/// outside (`cw == true`) or anticlockwise.
pub fn face_rotation(f: usize, cw: bool) -> Mat {
    // cw about (1,1,1): (x,y,z) -> (y,z,x)
    let p: Mat = [[0, 1, 0], [0, 0, 1], [1, 0, 0]];
    let pt = transpose(&p);
    let s = FACES[f];
    let d: Mat = [[s[0], 0, 0], [0, s[1], 0], [0, 0, s[2]]];
    // conjugating by the (possibly improper) diagonal matrix reverses the sense when det(d) = -1
    let base = if (s[0] * s[1] * s[2] == 1) == cw {
        p
    } else {
        pt
    };
    mat_mul(&mat_mul(&d, &base), &d)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    Corner,
    Edge,
    Centre,
}

#[derive(Clone, Copy, Debug)]
pub struct Sticker {
    pub face: usize,
    pub p: V3,
    pub kind: Kind,
}

/// The 72 stickers in a fixed order (face by face).
pub fn stickers() -> Vec<Sticker> {
    let pats: [(Kind, [V3; 3]); 3] = [
        (Kind::Corner, [[7, 1, 1], [1, 7, 1], [1, 1, 7]]),
        (Kind::Edge, [[1, 4, 4], [4, 1, 4], [4, 4, 1]]),
        (Kind::Centre, [[5, 2, 2], [2, 5, 2], [2, 2, 5]]),
    ];
    let mut v = vec![];
    for (f, s) in FACES.iter().enumerate() {
        for (kind, ps) in pats.iter() {
            for q in ps {
                v.push(Sticker {
                    face: f,
                    p: [q[0] * s[0], q[1] * s[1], q[2] * s[2]],
                    kind: *kind,
                });
            }
        }
    }
    v
}

pub struct StickerModel {
    pub st: Vec<Sticker>,
}

impl StickerModel {
    pub fn new() -> Self {
        StickerModel { st: stickers() }
    }
    pub fn index_of(&self, p: V3) -> usize {
        self.st
            .iter()
            .position(|s| s.p == p)
            .expect("no sticker at point")
    }
    /// Sticker permutation (`perm[i]` = destination of sticker i) of rotating every sticker
    /// that satisfies `pred` by `m`.
    pub fn perm_where(&self, m: &Mat, pred: impl Fn(V3) -> bool) -> Vec<usize> {
        self.st
            .iter()
            .enumerate()
            .map(|(i, s)| {
                if pred(s.p) {
                    self.index_of(mat_vec(m, s.p))
                } else {
                    i
                }
            })
            .collect()
    }
    /// Turn of the single layer of face `f`.
    pub fn face_turn(&self, f: usize, cw: bool) -> Vec<usize> {
        let n = FACES[f];
        self.perm_where(&face_rotation(f, cw), |p| dot(n, p) > 3)
    }
    /// Two-layer ("wide") turn of face `f`: everything except the layer of the opposite face.
    pub fn wide_turn(&self, f: usize, cw: bool) -> Vec<usize> {
        let o = FACES[opposite(f)];
        self.perm_where(&face_rotation(f, cw), |p| dot(o, p) <= 3)
    }
    /// Whole-puzzle symmetry (any signed permutation matrix).
    pub fn whole(&self, m: &Mat) -> Vec<usize> {
        self.perm_where(m, |_| true)
    }
}

/// All 48 symmetries of the octahedron (signed permutation matrices).
pub fn octahedral_group() -> Vec<Mat> {
    let perms = [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ];
    let mut v = vec![];
    for p in perms {
        for s in 0..8 {
            let mut m = [[0; 3]; 3];
            for i in 0..3 {
                m[i][p[i]] = if s >> i & 1 == 1 { -1 } else { 1 };
            }
            v.push(m);
        }
    }
    v
}

/// Symmetries mapping class-A faces to class-A faces (T_d, 24 elements), identity first.
pub fn class_preserving_group() -> Vec<Mat> {
    let mut v: Vec<Mat> = octahedral_group()
        .into_iter()
        .filter(|m| {
            A_FACES
                .iter()
                .all(|&f| is_a_face(face_of_normal(mat_vec(m, FACES[f]))))
        })
        .collect();
    v.sort_by_key(|m| (*m != IDENTITY, det(m) < 0));
    v
}

/// The 12 proper rotations of [`class_preserving_group`] (these are reachable by face turns).
pub fn rotation_group() -> Vec<Mat> {
    class_preserving_group()
        .into_iter()
        .filter(|m| det(m) == 1)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rotations_are_clockwise_from_outside() {
        for f in 0..8 {
            let r = face_rotation(f, true);
            assert_eq!(det(&r), 1);
            assert_eq!(mat_vec(&r, FACES[f]), FACES[f]);
            for v in VERTS {
                if dot(v, FACES[f]) == 1 {
                    let w = mat_vec(&r, v);
                    assert!(dot(cross(v, w), FACES[f]) < 0, "face {f} not clockwise");
                }
            }
            let rr = mat_mul(&r, &face_rotation(f, false));
            assert_eq!(rr, IDENTITY);
        }
    }

    #[test]
    fn group_sizes() {
        assert_eq!(octahedral_group().len(), 48);
        assert_eq!(class_preserving_group().len(), 24);
        assert_eq!(rotation_group().len(), 12);
        let sm = StickerModel::new();
        assert_eq!(sm.st.len(), 72);
    }

    #[test]
    fn opposite_faces() {
        for f in 0..8 {
            let o = opposite(f);
            for i in 0..3 {
                assert_eq!(FACES[o][i], -FACES[f][i]);
            }
            assert_ne!(is_a_face(f), is_a_face(o));
        }
    }
}
