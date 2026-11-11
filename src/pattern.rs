//! Positions entered as sticker colours ("patterns") instead of move sequences: what the
//! puzzle looks like, for example a scrambled puzzle whose scramble is not known.
//!
//! A pattern lists the 72 sticker colours as the puzzle is held with a corner pointing at
//! you, in picture order: the four front faces U, R, F, L (clockwise from the top), then the
//! four back faces B, BL, D, BR as seen after turning the whole puzzle around (left or right,
//! so that the top stays at the top: B is then at the top, BL on the right). Each face is read
//! with its outer edge at the top: the 5 stickers along that edge from left to right, the 3 in
//! the middle row, then the one at the corner in the middle of the picture.
//!
//! Colours are any 8 different characters, 9 stickers each. The colour scheme is not needed:
//! which colour belongs on which face is worked out from the pieces (every assignment of the 8
//! colours to the 8 faces is tried, and the one that makes the most pieces real FTO pieces
//! wins), so the puzzle may be held in any of its 24 orientations and the scheme may be any
//! arrangement of colours. A solution found for the relabelled position solves the pattern as
//! entered: renaming colours changes neither which stickers match nor what a turn does, so a
//! sequence that makes every face of the relabelled position one colour does the same for the
//! pattern. The assignment is also chosen so that the front corner is the solver's fixed
//! corner, in place (any assignment that is a whole-puzzle rotation of a valid one is valid,
//! and exactly one of the 12 that keep the face classes puts the front corner piece there), so
//! the position is solved in the frame of the puzzle as held and the solution's moves are the
//! turns to make, as held.

use crate::cube::{defs, perm_parity, Cube, Stickers};
use crate::geometry::{
    det, face_of_normal, is_a_face, mat_vec, octahedral_group, transpose, Kind, Mat, FACES,
    FACE_NAMES, V3,
};
use crate::notation::{apply_physical, PhysMove};
use std::fmt;
use std::sync::OnceLock;

/// Stickers in a pattern.
pub const STICKERS: usize = 72;

const X: V3 = [1, 0, 0];
const Y: V3 = [0, 1, 0];
const Z: V3 = [0, 0, 1];
const NX: V3 = [-1, 0, 0];
const NY: V3 = [0, -1, 0];
const NZ: V3 = [0, 0, -1];

/// The faces in picture order, each with its three corners: the corner in the middle of the
/// picture (front `+x`, back `-x`), then the ends of the face's outer edge, left then right as
/// seen with that edge at the top. Seen from the front, `+z` is the top-left corner of the
/// picture and `+y` the top-right; from the back (turned around), `+y` is the top-left.
const PICTURE: [(usize, [V3; 3]); 8] = [
    (0, [X, Z, Y]),    // U, the top of the front picture
    (1, [X, Y, NZ]),   // R
    (3, [X, NZ, NY]),  // F
    (2, [X, NY, Z]),   // L
    (4, [NX, Y, Z]),   // B, the top of the back picture
    (6, [NX, Z, NY]),  // BL
    (7, [NX, NY, NZ]), // D
    (5, [NX, NZ, Y]),  // BR
];

/// Faces in picture order (indices into [`FACE_NAMES`]).
pub const PICTURE_FACES: [usize; 8] = [0, 1, 3, 2, 4, 6, 7, 5];

/// The colour letter of each face (indexed like [`FACE_NAMES`]) when a pattern is printed: the
/// standard colour scheme (that of the Diansheng, Dayan and most other FTOs) held the standard
/// way, white on U and green on F (the scrambling orientation of FTO competitions): white U,
/// red R, purple L, green F, blue B, grey BR, orange BL, yellow D. Opposite faces are
/// white/yellow, green/blue, red/orange and purple/grey (checked on a physical puzzle, which
/// also settles left and right: red is on R). (Lanlan FTOs have red and green swapped.)
/// Patterns that are read may use any 8 characters.
pub const DEFAULT_SYMBOLS: [char; 8] = ['W', 'R', 'P', 'G', 'B', 'A', 'O', 'Y'];

/// `a * p + b * q + c * r`: a sticker centroid (scaled by 9, as in [`crate::geometry`]).
fn comb(a: i32, p: V3, b: i32, q: V3, c: i32, r: V3) -> V3 {
    std::array::from_fn(|i| a * p[i] + b * q[i] + c * r[i])
}

/// Sticker index (the numbering of [`crate::geometry::stickers`]) at each picture position.
pub fn picture_stickers() -> &'static [usize; STICKERS] {
    static P: OnceLock<[usize; STICKERS]> = OnceLock::new();
    P.get_or_init(|| {
        let sm = &defs().sm;
        let mut out = [0; STICKERS];
        for (i, &(_, [a, l, r])) in PICTURE.iter().enumerate() {
            // reading order: the outer row (corner, centre, edge, centre, corner), the middle
            // row (edge, centre, edge), the corner in the middle of the picture
            let at = [
                comb(7, l, 1, a, 1, r),
                comb(5, l, 2, a, 2, r),
                comb(1, a, 4, l, 4, r),
                comb(5, r, 2, a, 2, l),
                comb(7, r, 1, a, 1, l),
                comb(1, r, 4, l, 4, a),
                comb(5, a, 2, l, 2, r),
                comb(1, l, 4, r, 4, a),
                comb(7, a, 1, l, 1, r),
            ];
            for (j, &p) in at.iter().enumerate() {
                out[9 * i + j] = sm.index_of(p);
            }
        }
        out
    })
}

/// Picture position of each sticker index (the inverse of [`picture_stickers`]).
fn sticker_pictures() -> &'static [usize; STICKERS] {
    static P: OnceLock<[usize; STICKERS]> = OnceLock::new();
    P.get_or_init(|| {
        let mut inv = [0; STICKERS];
        for (i, &s) in picture_stickers().iter().enumerate() {
            inv[s] = i;
        }
        inv
    })
}

/// A piece named in an error: its kind, the picture positions of its stickers and the faces it
/// is on (indices into [`FACE_NAMES`], in picture order).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Piece {
    pub kind: Kind,
    pub stickers: Vec<usize>,
    pub faces: Vec<usize>,
}

impl Piece {
    fn new(kind: Kind, sticker_indices: &[usize]) -> Piece {
        let pics = sticker_pictures();
        let mut stickers: Vec<usize> = sticker_indices.iter().map(|&s| pics[s]).collect();
        stickers.sort_unstable();
        let faces = stickers.iter().map(|&p| PICTURE_FACES[p / 9]).collect();
        Piece {
            kind,
            stickers,
            faces,
        }
    }
}

impl fmt::Display for Piece {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let names: Vec<&str> = self.faces.iter().map(|&x| FACE_NAMES[x]).collect();
        match (self.kind, names.as_slice()) {
            (Kind::Corner, [a, b, c, d]) => write!(f, "the corner between {a}, {b}, {c} and {d}"),
            (Kind::Edge, [a, b]) => write!(f, "the edge between {a} and {b}"),
            (_, [a, ..]) => write!(f, "a centre on {a}"),
            _ => write!(f, "a piece"),
        }
    }
}

/// Why a pattern cannot be made by turning (several can apply at once).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Unreachable {
    /// an odd permutation of the corners
    CornersSwapped,
    /// an odd number of corners turned over
    CornerFlipped,
    /// an odd permutation of the edges
    EdgesSwapped,
}

/// What is wrong with a pattern.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum PatternError {
    /// not 72 stickers (the count given)
    Length(usize),
    /// not 8 different colours (the count given)
    Colours(usize),
    /// colours that are not on exactly 9 stickers, with their counts
    Counts(Vec<(char, usize)>),
    /// pieces that are no FTO piece under the best assignment of colours to faces
    Pieces(Vec<Piece>),
    /// pairs of positions holding the same piece
    Duplicates(Vec<[Piece; 2]>),
    /// a pattern of real pieces that no sequence of turns makes
    Unreachable(Vec<Unreachable>),
}

impl PatternError {
    /// Picture positions of the stickers the error is about (empty if it is about the whole
    /// pattern).
    pub fn stickers(&self) -> Vec<usize> {
        let mut v: Vec<usize> = match self {
            PatternError::Pieces(p) => p.iter().flat_map(|p| p.stickers.clone()).collect(),
            PatternError::Duplicates(d) => d
                .iter()
                .flat_map(|[a, b]| a.stickers.iter().chain(&b.stickers).copied())
                .collect(),
            _ => vec![],
        };
        v.sort_unstable();
        v.dedup();
        v
    }
}

fn list(items: &[String]) -> String {
    match items {
        [] => String::new(),
        [a] => a.clone(),
        [init @ .., last] => format!("{} and {last}", init.join(", ")),
    }
}

impl fmt::Display for PatternError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PatternError::Length(n) => write!(
                f,
                "a pattern has 72 stickers (9 on each of the 8 faces), but this one has {n}"
            ),
            PatternError::Colours(n) => {
                write!(f, "an FTO has 8 colours, but this pattern has {n}")
            }
            PatternError::Counts(c) => {
                let parts: Vec<String> = c.iter().map(|(s, n)| format!("{s} is on {n}")).collect();
                write!(
                    f,
                    "each colour should be on exactly 9 stickers, but {}",
                    list(&parts)
                )
            }
            PatternError::Pieces(p) => {
                let names: Vec<String> = p.iter().map(ToString::to_string).collect();
                write!(
                    f,
                    "{} on an FTO with these colours (a sticker is probably entered wrong): {}",
                    if p.len() == 1 {
                        "this piece cannot exist"
                    } else {
                        "these pieces cannot exist"
                    },
                    list(&names)
                )
            }
            PatternError::Duplicates(d) => {
                let parts: Vec<String> = d.iter().map(|[a, b]| format!("{a} and {b}")).collect();
                write!(
                    f,
                    "two pieces have the same colours (a sticker is probably entered wrong): {}",
                    parts.join("; ")
                )
            }
            PatternError::Unreachable(u) => {
                let parts: Vec<String> = u
                    .iter()
                    .map(|x| {
                        match x {
                            Unreachable::CornersSwapped => "two corners would have to be swapped",
                            Unreachable::CornerFlipped => {
                                "a corner would have to be turned over in place"
                            }
                            Unreachable::EdgesSwapped => "two edges would have to be swapped",
                        }
                        .to_string()
                    })
                    .collect();
                write!(
                    f,
                    "this pattern cannot be made by turning the puzzle: {} (usually two stickers \
                     are entered the wrong way round)",
                    list(&parts)
                )
            }
        }
    }
}

impl std::error::Error for PatternError {}

/// Whether a character of a pattern's text is a sticker (not a separator).
fn is_sticker(c: char) -> bool {
    !c.is_whitespace() && !matches!(c, '/' | ',' | '|')
}

/// Read a pattern: 72 colour characters in picture order. Whitespace, `/`, `,` and `|` are
/// ignored, so faces and the two halves can be set apart.
///
/// # Errors
/// [`PatternError::Length`] if there are not 72 characters.
pub fn parse(text: &str) -> Result<[char; STICKERS], PatternError> {
    let v: Vec<char> = text.chars().filter(|&c| is_sticker(c)).collect();
    v.as_slice()
        .try_into()
        .map_err(|_| PatternError::Length(v.len()))
}

/// One pattern of a list ([`records`]): the line it starts on (from 1) and its text.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Record {
    pub line: usize,
    pub text: String,
}

/// Split a list of patterns (`solve --patterns`) into records, for [`parse`]. A `#` starts a
/// comment, to the end of its line. A pattern takes one line or several: lines are joined until
/// they hold 72 stickers. A blank line ends a pattern early, and so does a line that would take
/// it past 72 stickers (that line starts the next pattern), so a pattern with a sticker too
/// few or too many is reported on its own and does not run into the next one. A line of more
/// than 72 stickers is a pattern on its own.
pub fn records(text: &str) -> Vec<Record> {
    let mut out = vec![];
    // the pattern being read: its first line, its text and its stickers so far
    let mut cur: Option<(usize, String, usize)> = None;
    let mut flush = |cur: &mut Option<(usize, String, usize)>| {
        if let Some((line, text, _)) = cur.take() {
            out.push(Record { line, text });
        }
    };
    for (i, raw) in text.lines().enumerate() {
        let line = raw.split('#').next().unwrap_or("");
        let n = line.chars().filter(|&c| is_sticker(c)).count();
        if raw.trim().is_empty() {
            flush(&mut cur);
            continue;
        }
        if n == 0 {
            // a comment or separators only: neither part of a pattern nor the end of one
            continue;
        }
        if cur.as_ref().is_some_and(|c| c.2 + n > STICKERS) {
            flush(&mut cur);
        }
        let c = cur.get_or_insert_with(|| (i + 1, String::new(), 0));
        if !c.1.is_empty() {
            c.1.push(' ');
        }
        c.1.push_str(line.trim());
        c.2 += n;
        if c.2 >= STICKERS {
            flush(&mut cur);
        }
    }
    flush(&mut cur);
    out
}

/// The pattern as text: faces of 9 in picture order, the front and back halves set apart.
pub fn format(p: &[char; STICKERS]) -> String {
    let face = |i: usize| p[9 * i..9 * i + 9].iter().collect::<String>();
    let half = |r: std::ops::Range<usize>| r.map(face).collect::<Vec<_>>().join(" ");
    format!("{} / {}", half(0..4), half(4..8))
}

/// The pattern of the puzzle after `moves` from solved, held as at the start, with the colour
/// letters of [`DEFAULT_SYMBOLS`].
pub fn of_moves(moves: &[PhysMove]) -> [char; STICKERS] {
    let faces: Stickers = std::array::from_fn(|i| defs().sm.st[i].face as u8);
    let after = apply_physical(&faces, moves);
    let pics = picture_stickers();
    std::array::from_fn(|i| DEFAULT_SYMBOLS[after[pics[i]] as usize])
}

/// The pattern of a position (in the frame of [`crate::cube`], held as the solver holds it),
/// with the colour letters of [`DEFAULT_SYMBOLS`].
pub fn of_cube(c: &Cube) -> [char; STICKERS] {
    let st = c.to_stickers();
    let pics = picture_stickers();
    std::array::from_fn(|i| DEFAULT_SYMBOLS[st[pics[i]] as usize])
}

/// The distinct colours of a pattern in order of first appearance, and the pattern as their
/// indices in sticker order (not picture order).
fn colour_indices(p: &[char; STICKERS]) -> (Vec<char>, [u8; STICKERS]) {
    let mut colours: Vec<char> = vec![];
    let mut s = [0u8; STICKERS];
    for (i, &c) in p.iter().enumerate() {
        let k = colours.iter().position(|&x| x == c).unwrap_or_else(|| {
            colours.push(c);
            colours.len() - 1
        });
        s[picture_stickers()[i]] = k as u8;
    }
    (colours, s)
}

/// The pattern after `moves`.
pub fn after(p: &[char; STICKERS], moves: &[PhysMove]) -> [char; STICKERS] {
    let (colours, s) = colour_indices(p);
    let s = apply_physical(&s, moves);
    let pics = picture_stickers();
    std::array::from_fn(|i| colours.get(s[pics[i]] as usize).copied().unwrap_or('?'))
}

/// The colour of each face (picture order) if every face is one colour.
pub fn face_colours(p: &[char; STICKERS]) -> Option<[char; 8]> {
    let mut out = ['?'; 8];
    for (i, o) in out.iter_mut().enumerate() {
        let face = &p[9 * i..9 * i + 9];
        if face.iter().any(|&c| c != face[0]) {
            return None;
        }
        *o = face[0];
    }
    Some(out)
}

/// True if `moves` turn the pattern into one with every face a single colour.
pub fn solves(p: &[char; STICKERS], moves: &[PhysMove]) -> bool {
    face_colours(&after(p, moves)).is_some()
}

/// A solved pattern with the face colours `colours` (picture order).
pub fn solved(colours: &[char; 8]) -> [char; STICKERS] {
    std::array::from_fn(|i| colours[i / 9])
}

/// A scramble that makes a pattern: the solution reversed, applied to a solved puzzle.
///
/// The puzzle has no fixed centres, so a solved puzzle can be held 24 ways, and the reversed
/// solution makes the pattern exactly only from the hold the solution ends in (`start`): from
/// another hold the same turns make a different puzzle. A usual hold is turned into `start` by
/// whole-puzzle rotations of 120° about a face (`v` turns such as `Lv`, at most two) for 12 of
/// the 24 holds: those keep each colour in its half of the faces (U F BR BL or R L B D, which
/// alternate around every corner). For the other 12 (`start` is the usual hold turned so the
/// halves swap, 90° about a corner) no turns from the usual hold make the pattern as entered;
/// the scramble with each turn renamed for the turned puzzle makes the same puzzle, which then
/// matches the pattern once the whole puzzle is turned.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Scramble {
    /// colour of each face (picture order) of the solved puzzle to start from
    pub start: [char; 8],
    /// the turns
    pub moves: Vec<PhysMove>,
    /// from the usual hold, if one was given and it is a hold of this puzzle's colours
    pub from_usual: Option<FromUsual>,
}

/// A scramble from the usual hold.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum FromUsual {
    /// the rotations to `start` (empty if the usual hold is `start`), then the turns: the
    /// pattern exactly as entered
    Exact(Vec<String>),
    /// turns that make the same puzzle, held another way than when the pattern was entered
    SamePuzzle(Vec<PhysMove>),
}

impl Scramble {
    /// The scramble from the usual hold as text, and whether it makes the pattern exactly.
    pub fn from_usual_text(&self) -> Option<(String, bool)> {
        match self.from_usual.as_ref()? {
            FromUsual::Exact(r) => {
                let m = crate::notation::format_moves(&self.moves);
                let t = r
                    .iter()
                    .map(String::as_str)
                    .chain((!m.is_empty()).then_some(m.as_str()))
                    .collect::<Vec<_>>()
                    .join(" ");
                Some((t, true))
            }
            FromUsual::SamePuzzle(m) => Some((crate::notation::format_moves(m), false)),
        }
    }
}

/// The scramble that `solution` (which must solve the pattern) gives, and the one from
/// `usual` (the face colours of the solved puzzle as usually held, picture order) if given.
/// `None` if `solution` does not solve the pattern.
pub fn scramble(
    p: &[char; STICKERS],
    solution: &[PhysMove],
    usual: Option<&[char; 8]>,
) -> Option<Scramble> {
    let start = face_colours(&after(p, solution))?;
    let moves = crate::notation::invert_moves(solution);
    let from_usual = usual.and_then(|u| {
        if let Some(r) = rotations_between(u, &start) {
            return Some(FromUsual::Exact(r));
        }
        // start is the usual hold turned by some rotation s: the turns renamed by s^-1 make,
        // from the usual hold, the pattern held turned back by s^-1
        let s = rotations()
            .into_iter()
            .find(|r| face_colours(&held(r, &solved(u))).as_ref() == Some(&start))?;
        Some(FromUsual::SamePuzzle(turned(&transpose(&s), &moves)))
    });
    Some(Scramble {
        start,
        moves,
        from_usual,
    })
}

/// The 24 rotations of the whole puzzle.
fn rotations() -> Vec<Mat> {
    octahedral_group()
        .into_iter()
        .filter(|m| det(m) == 1)
        .collect()
}

/// `moves` as turned on the puzzle held rotated by `rot`.
fn turned(rot: &Mat, moves: &[PhysMove]) -> Vec<PhysMove> {
    moves
        .iter()
        .map(|&pm| PhysMove {
            face: face_of_normal(mat_vec(rot, FACES[pm.face])),
            ..pm
        })
        .collect()
}

/// The pattern `p` with the whole puzzle turned by `rot`.
fn held(rot: &Mat, p: &[char; STICKERS]) -> [char; STICKERS] {
    let whole = defs().sm.whole(rot);
    let pics = picture_stickers();
    let inv = sticker_pictures();
    let mut q = ['?'; STICKERS];
    for i in 0..STICKERS {
        q[inv[whole[pics[i]]]] = p[i];
    }
    q
}

/// The shortest `v` turns (at most two) that turn the solved puzzle with face colours `from`
/// into the one with `to`.
fn rotations_between(from: &[char; 8], to: &[char; 8]) -> Option<Vec<String>> {
    if from == to {
        return Some(vec![]);
    }
    let singles: Vec<(String, Vec<PhysMove>)> = PICTURE_FACES
        .iter()
        .flat_map(|&f| {
            [
                format!("{}v", FACE_NAMES[f]),
                format!("{}v'", FACE_NAMES[f]),
            ]
        })
        .filter_map(|n| crate::notation::parse_moves(&n).ok().map(|m| (n, m)))
        .collect();
    let start = solved(from);
    let reaches = |m: &[PhysMove]| face_colours(&after(&start, m)).as_ref() == Some(to);
    if let Some((n, _)) = singles.iter().find(|(_, m)| reaches(m)) {
        return Some(vec![n.clone()]);
    }
    for (a, ma) in &singles {
        for (b, mb) in &singles {
            let m: Vec<PhysMove> = ma.iter().chain(mb).copied().collect();
            if reaches(&m) {
                return Some(vec![a.clone(), b.clone()]);
            }
        }
    }
    None
}

/// Lookup tables for scoring an assignment of colours to faces.
struct Lookup {
    /// corner key (4 faces of 3 bits, in slot order) -> home corner * 2 + orientation + 1, or 0
    corner: Vec<u8>,
    /// edge key (A-slot face, B-slot face) -> home edge + 1, or 0
    edge: [u8; 64],
}

fn lookup() -> &'static Lookup {
    static L: OnceLock<Lookup> = OnceLock::new();
    L.get_or_init(|| {
        let d = defs();
        let mut corner = vec![0u8; 1 << 12];
        for k in 0..6 {
            for o in 0..2 {
                let key = (0..4).fold(0, |a, j| {
                    a | (d.corner_cols[k][(j + 2 * o) % 4] as usize) << (3 * j)
                });
                corner[key] = (2 * k + o + 1) as u8;
            }
        }
        let mut edge = [0u8; 64];
        for k in 0..12 {
            edge[d.edge_cols[k][0] as usize | (d.edge_cols[k][1] as usize) << 3] = k as u8 + 1;
        }
        Lookup { corner, edge }
    })
}

/// The pieces of `s` (colour indices in sticker order) under the assignment `f` (colour ->
/// face): the (home corner * 2 + orientation + 1, or 0 if no corner) at each vertex, the home
/// edge + 1 (or 0) at each edge position, and for each centre sticker whether its face is of
/// the position's class.
struct Decoded {
    corners: [u8; 6],
    edges: [u8; 12],
    centres_ok: [bool; 24],
}

fn decode(s: &[u8; STICKERS], f: [u8; 8]) -> Decoded {
    let d = defs();
    let l = lookup();
    let face = |i: usize| f[s[i] as usize] as usize;
    let corners = std::array::from_fn(|v| {
        let key = (0..4).fold(0, |a, j| a | face(d.corner_slots[v][j]) << (3 * j));
        l.corner[key]
    });
    let edges =
        std::array::from_fn(|e| l.edge[face(d.edge_slots[e][0]) | face(d.edge_slots[e][1]) << 3]);
    let centres_ok = std::array::from_fn(|i| {
        if i < 12 {
            is_a_face(face(d.a_pos[i]))
        } else {
            !is_a_face(face(d.b_pos[i - 12]))
        }
    });
    Decoded {
        corners,
        edges,
        centres_ok,
    }
}

fn score(s: &[u8; STICKERS], f: [u8; 8]) -> usize {
    let x = decode(s, f);
    x.corners.iter().filter(|&&c| c != 0).count()
        + x.edges.iter().filter(|&&e| e != 0).count()
        + x.centres_ok.iter().filter(|&&ok| ok).count()
}

/// Calls `visit` with every permutation of `0..8` (Heap's algorithm).
#[cfg(test)]
fn permutations(mut visit: impl FnMut(&[u8; 8])) {
    let mut a: [u8; 8] = std::array::from_fn(|i| i as u8);
    let mut c = [0usize; 8];
    visit(&a);
    let mut i = 1;
    while i < 8 {
        if c[i] < i {
            if i % 2 == 0 {
                a.swap(0, i);
            } else {
                a.swap(c[i], i);
            }
            visit(&a);
            c[i] += 1;
            i = 1;
        } else {
            c[i] = 0;
            i += 1;
        }
    }
}

/// The score of an assignment under which every piece is real.
const PERFECT: usize = 6 + 12 + 24;

/// The assignment read straight off a right pattern: the front corner's four colours go to
/// the faces of corner 0 in place, and each of the other four to the face opposite the one
/// front colour it never shares a corner with. `None` if the corners do not give a full
/// assignment (the pattern has a wrong sticker, and [`best_assignment`] finds the closest).
fn direct(s: &[u8; STICKERS]) -> Option<[u8; 8]> {
    const UNSET: u8 = 8;
    let d = defs();
    let mut f = [UNSET; 8];
    let front: [usize; 4] = std::array::from_fn(|j| s[d.corner_slots[0][j]] as usize);
    for (j, &c) in front.iter().enumerate() {
        if f[c] != UNSET {
            return None;
        }
        f[c] = d.corner_cols[0][j];
    }
    let mut meet = [[false; 8]; 8];
    for slots in &d.corner_slots {
        for &a in slots {
            for &b in slots {
                meet[s[a] as usize][s[b] as usize] = true;
            }
        }
    }
    for &c in &front {
        let mut far = (0..8).filter(|&x| x != c && !meet[c][x]);
        let x = far.next()?;
        if far.next().is_some() || f[x] != UNSET {
            return None;
        }
        f[x] = 7 - f[c];
    }
    Some(f)
}

/// A piece for [`best_assignment`]: its stickers' colours and what makes it real.
struct Part {
    /// colour index of each sticker (corners in slot order, edges A then B)
    colours: [u8; 4],
    kind: PartKind,
}

#[derive(Clone, Copy)]
enum PartKind {
    Corner,
    Edge,
    /// a centre sticker at a position of class A (true) or B (false)
    Centre(bool),
}

impl Part {
    fn real(&self, f: [u8; 8]) -> bool {
        let face = |i: usize| f[self.colours[i] as usize] as usize;
        match self.kind {
            PartKind::Corner => lookup().corner[(0..4).fold(0, |a, j| a | face(j) << (3 * j))] != 0,
            PartKind::Edge => lookup().edge[face(0) | face(1) << 3] != 0,
            PartKind::Centre(a) => is_a_face(face(0)) == a,
        }
    }
}

/// Branch and bound over the assignments of faces to colours 0, 1, 2, ...: `score` counts
/// the real pieces among those already decided, and a branch that cannot reach the best score
/// even if every undecided piece turns out real is cut.
struct Search<'a> {
    s: &'a [u8; STICKERS],
    /// pieces by the colour that decides them
    by_last: [Vec<Part>; 8],
    f: [u8; 8],
    best: (usize, bool, [u8; 8]),
}

impl Search<'_> {
    fn go(&mut self, k: usize, used: u8, score: usize, undecided: usize) {
        if k == 8 {
            let home = decode(self.s, self.f).corners[0] == 1;
            if (score, home) > (self.best.0, self.best.1) {
                self.best = (score, home, self.f);
            }
            return;
        }
        let decided = self.by_last[k].len();
        for face in 0..8u8 {
            if used >> face & 1 != 0 {
                continue;
            }
            self.f[k] = face;
            let gain = self.by_last[k].iter().filter(|p| p.real(self.f)).count();
            // (`<`, not `<=`: an equal score may still win on the tie-break)
            if score + gain + undecided - decided < self.best.0 {
                continue;
            }
            self.go(k + 1, used | 1 << face, score + gain, undecided - decided);
        }
    }
}

/// The assignment that makes the most pieces real; among those, one that puts corner 0 in
/// place at vertex 0 (one exists whenever every piece is real, see the module notes). Used only
/// for a pattern with a mistake, to name the pieces that are wrong. Branch and bound: a
/// pattern with a few wrong stickers visits a small part of the 8! assignments (the brute
/// force over all of them took about 10 ms).
fn best_assignment(s: &[u8; STICKERS]) -> [u8; 8] {
    let d = defs();
    let mut by_last: [Vec<Part>; 8] = Default::default();
    let mut add = |stickers: &[usize], kind: PartKind| {
        let mut colours = [0u8; 4];
        for (c, &i) in colours.iter_mut().zip(stickers) {
            *c = s[i];
        }
        // (a piece is decided once colours 0..=last have faces)
        let last = stickers.iter().map(|&i| s[i]).max().unwrap_or(0);
        by_last[usize::from(last)].push(Part { colours, kind });
    };
    for v in 0..6 {
        add(&d.corner_slots[v], PartKind::Corner);
    }
    for e in 0..12 {
        add(&d.edge_slots[e], PartKind::Edge);
    }
    for i in 0..12 {
        add(&[d.a_pos[i]], PartKind::Centre(true));
        add(&[d.b_pos[i]], PartKind::Centre(false));
    }
    let undecided = by_last.iter().map(Vec::len).sum();
    let mut search = Search {
        s,
        by_last,
        f: [0; 8],
        best: (0, false, [0; 8]),
    };
    search.go(0, 0, 0, undecided);
    search.best.2
}

/// The position of a pattern, in the frame of the puzzle as held (solve it with
/// [`crate::notation::Frame::identity`]), or what is wrong with the pattern.
///
/// # Errors
/// What is wrong with the pattern: the first of wrong sticker or colour counts, pieces that
/// do not exist, pieces that appear twice, and a pattern that turning cannot make.
pub fn position(p: &[char; STICKERS]) -> Result<Cube, PatternError> {
    let (colours, s) = colour_indices(p);
    if colours.len() != 8 {
        return Err(PatternError::Colours(colours.len()));
    }
    let counts: Vec<(char, usize)> = colours
        .iter()
        .enumerate()
        .map(|(k, &c)| (c, s.iter().filter(|&&x| x as usize == k).count()))
        .filter(|&(_, n)| n != 9)
        .collect();
    if !counts.is_empty() {
        return Err(PatternError::Counts(counts));
    }
    let f = match direct(&s).filter(|f| score(&s, *f) == PERFECT) {
        Some(f) => f,
        None => best_assignment(&s),
    };
    let x = decode(&s, f);
    let d = defs();
    let mut bad = vec![];
    for v in 0..6 {
        if x.corners[v] == 0 {
            bad.push(Piece::new(Kind::Corner, &d.corner_slots[v]));
        }
    }
    for e in 0..12 {
        if x.edges[e] == 0 {
            bad.push(Piece::new(Kind::Edge, &d.edge_slots[e]));
        }
    }
    for i in 0..24 {
        if !x.centres_ok[i] {
            let st = if i < 12 { d.a_pos[i] } else { d.b_pos[i - 12] };
            bad.push(Piece::new(Kind::Centre, &[st]));
        }
    }
    if !bad.is_empty() {
        return Err(PatternError::Pieces(bad));
    }
    let mut dup = vec![];
    for a in 0..6 {
        for b in a + 1..6 {
            if (x.corners[a] - 1) / 2 == (x.corners[b] - 1) / 2 {
                dup.push([
                    Piece::new(Kind::Corner, &d.corner_slots[a]),
                    Piece::new(Kind::Corner, &d.corner_slots[b]),
                ]);
            }
        }
    }
    for a in 0..12 {
        for b in a + 1..12 {
            if x.edges[a] == x.edges[b] {
                dup.push([
                    Piece::new(Kind::Edge, &d.edge_slots[a]),
                    Piece::new(Kind::Edge, &d.edge_slots[b]),
                ]);
            }
        }
    }
    if !dup.is_empty() {
        return Err(PatternError::Duplicates(dup));
    }
    let faces: Stickers = std::array::from_fn(|i| f[s[i] as usize]);
    // every piece is real and appears once, so this succeeds (with corner 0 in place)
    let c = Cube::from_stickers(&faces).ok_or(PatternError::Unreachable(vec![]))?;
    let mut why = vec![];
    if perm_parity(&c.cp) != 0 {
        why.push(Unreachable::CornersSwapped);
    }
    if c.co.iter().fold(0, |a, &o| a ^ o) != 0 {
        why.push(Unreachable::CornerFlipped);
    }
    if perm_parity(&c.ep) != 0 {
        why.push(Unreachable::EdgesSwapped);
    }
    if !why.is_empty() || c.cp[0] != 0 || c.co[0] != 0 {
        return Err(PatternError::Unreachable(why));
    }
    Ok(c)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::notation::{format_moves, invert_moves, parse_moves};
    use rand::seq::SliceRandom;
    use rand::{Rng, SeedableRng};

    /// True if the position's stickers are the pattern's with the colours renamed.
    fn same_up_to_colours(p: &[char; STICKERS], c: &Cube) -> bool {
        let st = c.to_stickers();
        let pics = picture_stickers();
        let mut map: Vec<(char, u8)> = vec![];
        for i in 0..STICKERS {
            let f = st[pics[i]];
            match map.iter().find(|(c, _)| *c == p[i]) {
                Some(&(_, g)) if g != f => return false,
                Some(_) => {}
                None => {
                    if map.iter().any(|&(_, g)| g == f) {
                        return false;
                    }
                    map.push((p[i], f));
                }
            }
        }
        true
    }

    #[test]
    fn scrambles_from_the_usual_hold() {
        // a pattern entered holding the puzzle any of the 24 ways: the scramble from the hold
        // the solution ends in makes it exactly, and so does the one from the usual hold (with
        // rotations) for the 12 holds that keep the face halves; for the other 12 the one from
        // the usual hold makes the same puzzle, held another way
        let usual: [char; 8] = PICTURE_FACES.map(|f| DEFAULT_SYMBOLS[f]);
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(11);
        let mut from_usual = 0;
        for rot in rotations() {
            for _ in 0..4 {
                let n = rng.gen_range(0..20);
                let m: Vec<PhysMove> = (0..n)
                    .map(|_| PhysMove::face(rng.gen_range(0..8), rng.gen_bool(0.5)))
                    .collect();
                let q = held(&rot, &of_moves(&m));
                let sol = turned(&rot, &invert_moves(&m));
                let sc = scramble(&q, &sol, Some(&usual)).unwrap();
                assert_eq!(after(&solved(&sc.start), &sc.moves), q);
                let keeps = is_a_face(face_of_normal(mat_vec(&rot, FACES[0])));
                let (text, exact) = sc.from_usual_text().unwrap();
                assert_eq!(exact, keeps);
                let made = after(&solved(&usual), &parse_moves(&text).unwrap());
                if let Some(FromUsual::Exact(r)) = &sc.from_usual {
                    assert!(r.len() <= 2);
                    assert_eq!(made, q, "{text}");
                    from_usual += 1;
                } else {
                    // the same puzzle: the pattern once the whole puzzle is turned
                    assert!(rotations().iter().any(|r| held(r, &made) == q), "{text}");
                }
            }
        }
        assert_eq!(from_usual, 12 * 4);
    }

    #[test]
    fn picture_order_covers_every_sticker_on_its_face() {
        let p = picture_stickers();
        let mut seen = [false; STICKERS];
        for (i, &s) in p.iter().enumerate() {
            assert!(!seen[s], "sticker {s} twice");
            seen[s] = true;
            assert_eq!(defs().sm.st[s].face, PICTURE[i / 9].0);
            assert_eq!(PICTURE_FACES[i / 9], PICTURE[i / 9].0);
            // the corner, edge and centre stickers are where the reading order puts them
            let kind = match i % 9 {
                0 | 4 | 8 => Kind::Corner,
                2 | 5 | 7 => Kind::Edge,
                _ => Kind::Centre,
            };
            assert_eq!(defs().sm.st[s].kind, kind, "picture position {i}");
        }
        for (f, &(face, [a, l, r])) in PICTURE.iter().enumerate() {
            // the corners are the face's, and l -> r runs clockwise around the picture as seen
            // from outside the face: (l - a) x (r - a) points out of the face
            let n = crate::geometry::cross(comb(1, l, -1, a, 0, r), comb(1, r, -1, a, 0, l));
            assert!(
                crate::geometry::dot(n, FACES[face]) < 0,
                "face {f}: l -> r is anticlockwise seen from outside"
            );
        }
    }

    #[test]
    fn neighbouring_pictures_share_their_edges() {
        // stickers 0..4 of a face run along its outer edge; the two halves' outer edges are the
        // same octahedron edges, so their edge pieces are the same
        let d = defs();
        let pics = sticker_pictures();
        for e in 0..12 {
            let [a, b] = d.edge_slots[e];
            let (pa, pb) = (pics[a], pics[b]);
            assert_ne!(pa / 9, pb / 9);
        }
    }

    #[test]
    fn solved_pattern() {
        let p = of_moves(&[]);
        assert_eq!(
            format(&p),
            "WWWWWWWWW RRRRRRRRR GGGGGGGGG PPPPPPPPP / BBBBBBBBB OOOOOOOOO YYYYYYYYY AAAAAAAAA"
        );
        assert_eq!(position(&p), Ok(Cube::solved()));
    }

    #[test]
    fn single_turns_read_back_as_the_same_position() {
        // a turn of a face not touching the front corner keeps the frame, so the position is
        // the scramble's internal position exactly
        for name in ["B", "BR", "BL", "D", "B'", "D'"] {
            let m = parse_moves(name).unwrap();
            let (c, _) = crate::notation::scramble_to_cube(&m);
            assert_eq!(position(&of_moves(&m)), Ok(c), "{name}");
        }
    }

    #[test]
    fn scrambles_read_back_and_their_solutions_solve_the_pattern() {
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(5);
        let faces = ["U", "R", "L", "F", "B", "BR", "BL", "D"];
        for _ in 0..200 {
            let n = rng.gen_range(0..30);
            let text: Vec<String> = (0..n)
                .map(|_| {
                    let f = faces[rng.gen_range(0..8)];
                    let k = ["", "w", "s"][rng.gen_range(0..3)];
                    let p = if rng.gen_bool(0.5) { "'" } else { "" };
                    format!("{f}{k}{p}")
                })
                .collect();
            let m = parse_moves(&text.join(" ")).unwrap();
            let p = of_moves(&m);
            let c = position(&p).unwrap_or_else(|e| panic!("{}: {e}", text.join(" ")));
            // the position shows the pattern (in the frame of the puzzle as held), and the
            // scramble's inverse solves the pattern
            assert!(same_up_to_colours(&p, &c), "{}", text.join(" "));
            // a right pattern is read without the search over every assignment
            assert!(direct(&colour_indices(&p).1)
                .is_some_and(|f| score(&colour_indices(&p).1, f) == PERFECT));
            let inv = invert_moves(&m);
            assert!(solves(&p, &inv), "{}", format_moves(&inv));
        }
    }

    #[test]
    fn any_colours_and_orientation() {
        // recolour a scramble's pattern with an arbitrary bijection and hold it in any of the
        // 24 orientations: it still reads, and the solution of the read position solves it
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(9);
        let rotations = rotations();
        assert_eq!(rotations.len(), 24);
        for t in 0..60 {
            let n = rng.gen_range(5..25);
            let m: Vec<PhysMove> = (0..n)
                .map(|_| PhysMove::face(rng.gen_range(0..8), rng.gen_bool(0.5)))
                .collect();
            let p = of_moves(&m);
            // hold the puzzle differently
            let rot = rotations[t % rotations.len()];
            let q = held(&rot, &p);
            // recolour
            let letters = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'];
            let mut shuffled = letters;
            shuffled.shuffle(&mut rng);
            let q: [char; STICKERS] = std::array::from_fn(|i| {
                shuffled[DEFAULT_SYMBOLS.iter().position(|&c| c == q[i]).unwrap()]
            });
            let c = position(&q).unwrap_or_else(|e| panic!("{e}"));
            assert!(same_up_to_colours(&q, &c));
            // the scramble's inverse, as turned on the rotated puzzle, solves it
            assert!(solves(&q, &turned(&rot, &invert_moves(&m))));
        }
    }

    #[test]
    fn errors() {
        let solved = of_moves(&[]);
        assert_eq!(parse("W"), Err(PatternError::Length(1)));
        // one sticker of the wrong colour
        let mut p = solved;
        p[2] = 'R';
        assert!(matches!(position(&p), Err(PatternError::Counts(_))));
        // two stickers of one edge exchanged: the counts are right, the edge is no piece
        let mut p = solved;
        let d = defs();
        let inv = sticker_pictures();
        let [a, b] = d.edge_slots[0];
        p.swap(inv[a], inv[b]);
        let e = position(&p).unwrap_err();
        assert!(
            matches!(e, PatternError::Pieces(ref v) if v.len() == 1),
            "{e:?}"
        );
        assert_eq!(e.stickers().len(), 2);
        // two edges exchanged (both stickers): real pieces, but odd
        let mut p = solved;
        let [c, dd] = d.edge_slots[1];
        p.swap(inv[a], inv[c]);
        p.swap(inv[b], inv[dd]);
        assert_eq!(
            position(&p),
            Err(PatternError::Unreachable(vec![Unreachable::EdgesSwapped]))
        );
        // a corner turned over (its four stickers rotated by half a turn)
        let mut p = solved;
        let s = d.corner_slots[3].map(|x| inv[x]);
        let cols = s.map(|x| p[x]);
        for j in 0..4 {
            p[s[j]] = cols[(j + 2) % 4];
        }
        assert_eq!(
            position(&p),
            Err(PatternError::Unreachable(vec![Unreachable::CornerFlipped]))
        );
        assert!(position(&p)
            .unwrap_err()
            .to_string()
            .contains("turned over"));
    }

    #[test]
    fn any_input_gives_a_position_or_an_error() {
        // no input panics: random letters (few, so the counts sometimes pass), and right
        // patterns with stickers swapped or overwritten. A pattern that reads is shown by its
        // position, every error names stickers that exist, and its message is not empty.
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(17);
        let (mut ok, mut errors) = (0, 0);
        for t in 0..3000 {
            let p: [char; STICKERS] = if t % 3 == 0 {
                let k = rng.gen_range(1..10);
                std::array::from_fn(|_| (b'a' + rng.gen_range(0..k)) as char)
            } else {
                let m: Vec<PhysMove> = (0..rng.gen_range(0..20))
                    .map(|_| PhysMove::face(rng.gen_range(0..8), rng.gen_bool(0.5)))
                    .collect();
                let mut p = of_moves(&m);
                for _ in 0..rng.gen_range(0..4) {
                    let (i, j) = (rng.gen_range(0..STICKERS), rng.gen_range(0..STICKERS));
                    if t % 3 == 1 {
                        p.swap(i, j);
                    } else {
                        p[i] = p[j];
                    }
                }
                p
            };
            match position(&p) {
                Ok(c) => {
                    assert!(same_up_to_colours(&p, &c), "{}", format(&p));
                    ok += 1;
                }
                Err(e) => {
                    assert!(e.stickers().iter().all(|&i| i < STICKERS), "{e:?}");
                    assert_ne!(e.to_string(), "");
                    errors += 1;
                }
            }
            let m: Vec<PhysMove> = (0..rng.gen_range(0..6))
                .map(|_| PhysMove::face(rng.gen_range(0..8), rng.gen_bool(0.5)))
                .collect();
            let _ = scramble(&p, &m, Some(&['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']));
        }
        assert!(ok > 500 && errors > 500, "{ok} read, {errors} errors");
    }

    #[test]
    fn branch_and_bound_finds_the_best_score() {
        // against the brute force over all 8! assignments, on corrupted patterns
        let mut rng = rand_xoshiro::Xoshiro256PlusPlus::seed_from_u64(23);
        let mut checked = 0;
        for t in 0..30 {
            let m: Vec<PhysMove> = (0..rng.gen_range(0..20))
                .map(|_| PhysMove::face(rng.gen_range(0..8), rng.gen_bool(0.5)))
                .collect();
            let mut p = of_moves(&m);
            for _ in 0..=(t % 6) {
                let (i, j) = (rng.gen_range(0..STICKERS), rng.gen_range(0..STICKERS));
                p.swap(i, j);
            }
            let (colours, s) = colour_indices(&p);
            if colours.len() != 8 {
                continue;
            }
            let mut brute = (0, false);
            permutations(|f| {
                let x = (score(&s, *f), decode(&s, *f).corners[0] == 1);
                brute = brute.max(x);
            });
            let f = best_assignment(&s);
            assert_eq!(
                (score(&s, f), decode(&s, f).corners[0] == 1),
                brute,
                "{}",
                format(&p)
            );
            checked += 1;
        }
        assert!(checked > 25);
    }

    #[test]
    fn records_of_a_list() {
        let p = format(&of_moves(&parse_moves("R U' BL F").unwrap()));
        let q = format(&of_moves(&parse_moves("D B").unwrap()));
        let (front, back) = p.split_once(" / ").unwrap();
        let short = &q[1..];
        // one per line, a comment, a pattern over two lines (with a comment between), a blank
        // line ending a pattern a sticker short, and a pattern a sticker short followed at once
        // by a whole one
        let text = format!(
            "# my patterns\n{p}\n\n  {q}  # D B\n{front}\n# the back\n{back}\n{short}\n\n\
             {short}\n{p}\n{p}X\n"
        );
        let r = records(&text);
        let texts: Vec<&str> = r.iter().map(|x| x.text.as_str()).collect();
        let joined = format!("{front} {back}");
        let lines: Vec<usize> = r.iter().map(|x| x.line).collect();
        let longer = format!("{p}X");
        assert_eq!(
            texts,
            [
                p.as_str(),
                q.as_str(),
                joined.as_str(),
                short,
                short,
                p.as_str(),
                longer.as_str()
            ]
        );
        assert_eq!(lines, [2, 4, 5, 8, 10, 11, 12]);
        assert_eq!(
            parse(&r[2].text),
            Ok(of_moves(&parse_moves("R U' BL F").unwrap()))
        );
        assert_eq!(parse(&r[3].text), Err(PatternError::Length(71)));
        assert_eq!(parse(&r[6].text), Err(PatternError::Length(73)));
        assert_eq!(records("\n# nothing\n  \n"), Vec::<Record>::new());
    }

    #[test]
    fn pattern_of_a_position() {
        // a scramble whose turns keep the fixed corner: the same pattern as from the moves
        let m = parse_moves("B BR' D BL").unwrap();
        let (c, _) = crate::notation::scramble_to_cube(&m);
        assert_eq!(of_cube(&c), of_moves(&m));
        assert_eq!(position(&of_cube(&c)), Ok(c));
    }

    #[test]
    fn piece_names() {
        let d = defs();
        let p = Piece::new(Kind::Corner, &d.corner_slots[0]);
        assert_eq!(p.to_string(), "the corner between U, R, F and L");
    }
}
