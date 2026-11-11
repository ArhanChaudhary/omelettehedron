// The pattern entry of the demo page: a palette of the eight colours and the picture of the
// puzzle (pattern-picture.mjs), painted by clicking or dragging over the stickers, plus a
// pattern code that can be copied, edited or pasted. Letters are those of the solved faces
// (SOLVED_LETTERS), so the palette is also the visitor's colour scheme: the colour of each face
// of their solved puzzle as they usually hold it.

import { stickerTriangles, faceLabels, PICTURE_FACES, SOLVED_LETTERS } from './pattern-picture.mjs';

const DEFAULTS = {
  W: ['#f4f4f0', 'white'],
  R: ['#d42a2f', 'red'],
  G: ['#1f9a3b', 'green'],
  P: ['#7b41bf', 'purple'],
  B: ['#2260d0', 'blue'],
  O: ['#f47c0c', 'orange'],
  Y: ['#f4d10b', 'yellow'],
  A: ['#8b8b87', 'grey'],
};
const STORE_COLOURS = 'omelettehedron-pattern-colours';
const STORE_CELLS = 'omelettehedron-pattern';
const SVG = 'http://www.w3.org/2000/svg';
const SIZE = 200;
const FRONT = [10, 30];
const BACK = [240, 30];

function load(key) {
  try {
    return JSON.parse(localStorage.getItem(key) ?? 'null');
  } catch {
    return null;
  }
}

function save(key, value) {
  try {
    localStorage.setItem(key, JSON.stringify(value));
  } catch {
    // private windows and blocked storage: the page works without it
  }
}

function el(tag, attrs = {}, ns = null) {
  const e = ns ? document.createElementNS(ns, tag) : document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
  return e;
}

/** Black or white text, whichever reads better on `hex`. */
function ink(hex) {
  const n = parseInt(hex.slice(1), 16);
  const lum = 0.299 * (n >> 16) + 0.587 * ((n >> 8) & 255) + 0.114 * (n & 255);
  return lum > 140 ? '#111' : '#fff';
}

const STYLE = `
.pi-palette { display: flex; flex-wrap: wrap; gap: 0.35rem; align-items: center; margin: 0.4rem 0; }
.pi-swatch { width: 2.6rem; height: 2.1rem; padding: 0; border: 2px solid var(--line); border-radius: 6px; font-size: 0.8rem; font-weight: 600; }
.pi-swatch[aria-pressed="true"] { outline: 3px solid var(--accent); outline-offset: 2px; }
.pi-swatch input { display: block; width: 100%; height: 0.9rem; border: 0; padding: 0; background: none; }
.pi-picture { width: 100%; max-width: 34rem; display: block; touch-action: none; user-select: none; -webkit-user-select: none; }
.pi-picture polygon { stroke: var(--bg); stroke-width: 1.2; cursor: pointer; }
.pi-picture polygon.pi-blank { fill: color-mix(in srgb, var(--line) 55%, var(--bg)); }
.pi-picture polygon.pi-bad { stroke: #e0119a; stroke-width: 3.5; }
.pi-picture text { font: 600 12px system-ui, sans-serif; pointer-events: none; paint-order: stroke; stroke-width: 3px; }
.pi-picture text.pi-title { fill: var(--muted); stroke: none; font-weight: 400; }
.pi-code { width: 100%; box-sizing: border-box; font-family: ui-monospace, monospace; }
`;

export class PatternInput {
  /**
   * Build the pattern entry inside `root`. `onChange()` is called after every edit.
   */
  constructor(root, { onChange = () => {} } = {}) {
    this.onChange = onChange;
    const stored = load(STORE_COLOURS);
    this.colours = Object.fromEntries(
      SOLVED_LETTERS.map((l) => [l, typeof stored?.[l] === 'string' ? stored[l] : DEFAULTS[l][0]]),
    );
    const cells = load(STORE_CELLS);
    this.cells = Array.isArray(cells) && cells.length === 72 && cells.every((c) => c === null || c in DEFAULTS)
      ? cells
      : Array(72).fill(null);
    this.brush = 'W';
    this.editing = false;

    if (!document.getElementById('pi-style')) {
      const style = el('style', { id: 'pi-style' });
      style.textContent = STYLE;
      document.head.append(style);
    }

    // palette
    this.palette = el('div', { class: 'pi-palette', role: 'group', 'aria-label': 'Colours' });
    this.swatches = {};
    SOLVED_LETTERS.forEach((letter, f) => {
      const b = el('button', { type: 'button', class: 'pi-swatch' });
      b.addEventListener('click', (e) => {
        if (e.target.tagName === 'INPUT') return;
        this.brush = letter;
        this.renderPalette();
      });
      const picker = el('input', { type: 'color', 'aria-label': `Colour of ${PICTURE_FACES[f]}` });
      picker.hidden = true;
      picker.addEventListener('input', () => {
        this.colours[letter] = picker.value;
        save(STORE_COLOURS, this.colours);
        this.render();
      });
      b.append(document.createTextNode(PICTURE_FACES[f]), picker);
      this.swatches[letter] = { button: b, picker };
      this.palette.append(b);
    });
    const blank = el('button', { type: 'button', class: 'pi-swatch', title: 'Eraser' });
    blank.textContent = '×';
    blank.addEventListener('click', () => {
      this.brush = null;
      this.renderPalette();
    });
    this.swatches[''] = { button: blank };
    this.palette.append(blank);
    this.editButton = el('button', { type: 'button' });
    this.editButton.addEventListener('click', () => {
      this.editing = !this.editing;
      this.renderPalette();
    });
    this.palette.append(this.editButton);

    // picture
    this.svg = el('svg', { class: 'pi-picture', viewBox: `0 0 450 ${FRONT[1] + SIZE + 6}`, role: 'img' }, SVG);
    this.svg.setAttribute('aria-label', 'Picture of the puzzle: front and back');
    const title = (x, text) => {
      const t = el('text', { x, y: 18, 'text-anchor': 'middle', class: 'pi-title' }, SVG);
      t.textContent = text;
      this.svg.append(t);
    };
    title(FRONT[0] + SIZE / 2, 'Front');
    title(BACK[0] + SIZE / 2, 'Back (puzzle turned around)');
    this.polys = stickerTriangles({ size: SIZE, front: FRONT, back: BACK, gap: 0.06 }).map((t, i) => {
      const p = el('polygon', { points: t.points.map(([x, y]) => `${x.toFixed(2)},${y.toFixed(2)}`).join(' ') }, SVG);
      p.dataset.i = i;
      this.svg.append(p);
      return p;
    });
    this.labels = faceLabels({ size: SIZE, front: FRONT, back: BACK }).map(({ name, x, y }) => {
      const t = el('text', { x, y: y + 4, 'text-anchor': 'middle' }, SVG);
      t.textContent = name;
      this.svg.append(t);
      return t;
    });
    let painting = false;
    const paintAt = (e) => {
      const target = document.elementFromPoint(e.clientX, e.clientY);
      const i = target?.dataset?.i;
      if (i === undefined || this.cells[i] === this.brush) return;
      this.cells[i] = this.brush;
      this.edited();
    };
    this.svg.addEventListener('pointerdown', (e) => {
      painting = true;
      this.svg.setPointerCapture?.(e.pointerId);
      paintAt(e);
      e.preventDefault();
    });
    this.svg.addEventListener('pointermove', (e) => painting && paintAt(e));
    for (const end of ['pointerup', 'pointercancel', 'lostpointercapture']) {
      this.svg.addEventListener(end, () => { painting = false; });
    }

    // pattern code
    this.code = el('input', {
      class: 'pi-code',
      'aria-label': 'Pattern code',
      autocomplete: 'off',
      autocapitalize: 'characters',
      spellcheck: 'false',
    });
    this.code.addEventListener('change', () => this.setCode(this.code.value));

    root.append(this.palette, this.svg, this.code);
    this.render();
  }

  /** The name of a colour for messages: its colour name, or the face it is on when solved. */
  colourName(letter) {
    if (this.colours[letter]?.toLowerCase() === DEFAULTS[letter]?.[0]) return DEFAULTS[letter][1];
    return `your ${PICTURE_FACES[SOLVED_LETTERS.indexOf(letter)]} colour`;
  }

  /** Number of stickers without a colour. */
  get blanks() {
    return this.cells.filter((c) => c === null).length;
  }

  /** The pattern as the solver reads it (src/pattern.rs), `.` for a blank sticker. */
  get pattern() {
    const faces = [];
    for (let f = 0; f < 8; f++) faces.push(this.cells.slice(9 * f, 9 * f + 9).map((c) => c ?? '.').join(''));
    return `${faces.slice(0, 4).join(' ')} / ${faces.slice(4).join(' ')}`;
  }

  /** The face colours of the usual hold, for the solver's `usual` (always the solved letters). */
  get usual() {
    return SOLVED_LETTERS.join('');
  }

  /** Set the pattern from a code; returns an error message, or null. */
  setCode(text) {
    const letters = [...text.toUpperCase()].filter((c) => !/[\s/,|]/.test(c));
    if (letters.length !== 72) return this.codeError(`A pattern code has 72 letters, but this one has ${letters.length}.`);
    const bad = letters.find((c) => c !== '.' && !(c in DEFAULTS));
    if (bad) return this.codeError(`The letters of a pattern code are ${SOLVED_LETTERS.join(' ')} (and . for a blank sticker), not ${bad}.`);
    this.cells = letters.map((c) => (c === '.' ? null : c));
    this.edited();
    return null;
  }

  codeError(message) {
    this.code.setCustomValidity(message);
    this.code.reportValidity();
    return message;
  }

  /** Every sticker the colour of its solved face. */
  fillSolved() {
    this.cells = Array.from({ length: 72 }, (_, i) => SOLVED_LETTERS[Math.floor(i / 9)]);
    this.edited();
  }

  /** Every sticker blank. */
  clear() {
    this.cells = Array(72).fill(null);
    this.edited();
  }

  /** Outline the stickers at these picture positions (until the next edit). */
  highlight(positions) {
    const set = new Set(positions);
    this.polys.forEach((p, i) => p.classList.toggle('pi-bad', set.has(i)));
  }

  edited() {
    this.code.setCustomValidity('');
    this.highlight([]);
    save(STORE_CELLS, this.cells);
    this.render();
    this.onChange();
  }

  renderPalette() {
    for (const [letter, { button, picker }] of Object.entries(this.swatches)) {
      button.setAttribute('aria-pressed', String((letter || null) === this.brush));
      if (!picker) continue;
      const colour = this.colours[letter];
      button.style.background = colour;
      button.style.color = ink(colour);
      picker.value = colour;
      picker.hidden = !this.editing;
      const face = PICTURE_FACES[SOLVED_LETTERS.indexOf(letter)];
      button.title = `Paint with ${this.colourName(letter)} (the colour of ${face} on your solved puzzle)`;
    }
    this.editButton.textContent = this.editing ? 'Done' : 'Change colours';
  }

  render() {
    this.renderPalette();
    this.polys.forEach((p, i) => {
      const c = this.cells[i];
      p.classList.toggle('pi-blank', c === null);
      if (c === null) p.removeAttribute('fill');
      else p.setAttribute('fill', this.colours[c]);
    });
    this.labels.forEach((t, f) => {
      // the label takes the ink of the sticker under it
      const c = this.cells[9 * f + 6];
      const fill = c ? ink(this.colours[c]) : 'var(--fg)';
      t.setAttribute('fill', fill);
      t.setAttribute('stroke', fill === '#fff' ? 'rgba(0,0,0,0.45)' : 'rgba(255,255,255,0.55)');
    });
    if (document.activeElement !== this.code) this.code.value = this.pattern;
  }
}
