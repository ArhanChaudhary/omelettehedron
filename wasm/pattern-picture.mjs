// The picture of the puzzle for entering a pattern (src/pattern.rs): the front and the back as
// two squares, each cut along its diagonals into four faces of 9 stickers.
//
// Hold the puzzle with a corner pointing at you: the front square shows U (top), R (right),
// F (bottom) and L (left) around that corner. Turn the whole puzzle around, keeping the top at
// the top, and the back square shows B (top), BL (right), D (bottom) and BR (left). Sticker i
// of the picture (0..71) is pattern position i: the faces in the order U R F L B BL D BR, and
// each face read with its outer edge at the top (5 stickers along that edge from left to right,
// the 3 of the middle row, then the one at the corner in the middle of the square).

/** Face names in picture order. */
export const PICTURE_FACES = ['U', 'R', 'F', 'L', 'B', 'BL', 'D', 'BR'];

/** Pattern letters of the solved faces in picture order (`pattern::DEFAULT_SYMBOLS`). */
export const SOLVED_LETTERS = ['W', 'R', 'G', 'P', 'B', 'O', 'Y', 'A'];

// The 9 stickers of a face in reading order, as the grid points of their corners: weights
// (out of 3) of the face's left outer corner l, the middle corner a and the right outer
// corner r. Their sums are the sticker centroids of src/pattern.rs (`picture_stickers`).
const STICKER_POINTS = [
  [[3, 0, 0], [2, 1, 0], [2, 0, 1]], // corner at l
  [[2, 1, 0], [2, 0, 1], [1, 1, 1]], // centre at l
  [[1, 1, 1], [2, 0, 1], [1, 0, 2]], // edge on the outer side
  [[0, 1, 2], [1, 0, 2], [1, 1, 1]], // centre at r
  [[0, 0, 3], [0, 1, 2], [1, 0, 2]], // corner at r
  [[1, 1, 1], [2, 1, 0], [1, 2, 0]], // edge on the side l-a
  [[1, 2, 0], [0, 2, 1], [1, 1, 1]], // centre at a
  [[1, 1, 1], [0, 1, 2], [0, 2, 1]], // edge on the side r-a
  [[0, 3, 0], [1, 2, 0], [0, 2, 1]], // corner at a
];

/**
 * The 72 sticker triangles for squares of side `size` whose top-left corners are `front` and
 * `back` ([x, y]), each shrunk by `gap` toward its centre. Returns, in picture order,
 * { points: [[x, y] x 3], face: picture face index }.
 */
export function stickerTriangles({ size = 200, front = [0, 0], back = [220, 0], gap = 0.12 } = {}) {
  const out = [];
  for (let f = 0; f < 8; f++) {
    const [ox, oy] = f < 4 ? front : back;
    // the square's corners clockwise from the top-left, and its centre
    const corners = [[ox, oy], [ox + size, oy], [ox + size, oy + size], [ox, oy + size]];
    const a = [ox + size / 2, oy + size / 2];
    const l = corners[f % 4];
    const r = corners[(f + 1) % 4];
    for (const tri of STICKER_POINTS) {
      const pts = tri.map(([wl, wa, wr]) => [
        (wl * l[0] + wa * a[0] + wr * r[0]) / 3,
        (wl * l[1] + wa * a[1] + wr * r[1]) / 3,
      ]);
      const cx = (pts[0][0] + pts[1][0] + pts[2][0]) / 3;
      const cy = (pts[0][1] + pts[1][1] + pts[2][1]) / 3;
      out.push({
        points: pts.map(([x, y]) => [x + (cx - x) * gap, y + (cy - y) * gap]),
        face: f,
      });
    }
  }
  return out;
}

/** Where to write each face's name: on the face's sticker 6 (the middle of the face). */
export function faceLabels({ size = 200, front = [0, 0], back = [220, 0] } = {}) {
  const tris = stickerTriangles({ size, front, back, gap: 0 });
  return PICTURE_FACES.map((name, f) => {
    const pts = tris[9 * f + 6].points;
    return {
      name,
      x: (pts[0][0] + pts[1][0] + pts[2][0]) / 3,
      y: (pts[0][1] + pts[1][1] + pts[2][1]) / 3,
    };
  });
}
