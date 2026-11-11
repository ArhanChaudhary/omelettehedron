// The browser demo's solver (index.html; also driven by test.mjs in Bun): runs the
// WebAssembly module in a worker, so that the page stays responsive and the module may block.
// Messages in: { id, cmd, ...args }; out: { id, result } or { id, error }. The module's
// output lines (progress while generating tables, each better solution while solving) go to
// the BroadcastChannel named at `init`, from this worker and from the pool's threads.

import { FtoSolver, wasmFile } from './omelettehedron.mjs';

let solver = null;

const commands = {
  /**
   * Start the module: the threads build if `threads` > 1, the one with room for the tables
   * `tier` ('tiny', 'small' or 'medium'). Returns `{ threads, control }`.
   */
  async init({ threads, channel, tier = 'small' }) {
    solver?.close();
    solver = null;
    const wasm = new URL(`dist/${wasmFile(tier, threads)}`, import.meta.url);
    solver = await FtoSolver.create({ wasm, threads: threads > 1 ? threads : undefined, channel });
    // (the control array goes to the page for `cancel`, if it can be shared)
    const shared =
      typeof SharedArrayBuffer !== 'undefined' && solver.control.buffer instanceof SharedArrayBuffer;
    return { threads: solver.threads, control: shared ? solver.control : null };
  },
  /** Generate the tables `tier` (progress on the channel), then warm the engine up. */
  async generate({ tier = 'small' } = {}) {
    const info = solver.generate({ tier });
    await solver.warmUp();
    return info;
  },
  /** Load the two table files (Blobs, URLs or bytes), then warm the engine up. */
  async load({ ae, cb }) {
    const info = await solver.loadTables(ae, cb);
    await solver.warmUp();
    return info;
  },
  /** The tables as files: `[{ name, blob }, { name, blob }]`. */
  export: () => solver.exportTables(),
  /** Solve, printing improvements on the channel. */
  solve: ({ scramble, seconds = 0, maxLength = 0 }) =>
    solver.solve(scramble, { seconds, maxLength, improvements: true }),
  /** Solve a pattern (72 sticker colours), printing improvements on the channel. */
  solvePattern: ({ pattern, usual = '', seconds = 0, maxLength = 0 }) =>
    solver.solvePattern(pattern, { usual, seconds, maxLength, improvements: true }),
  pattern: ({ scramble }) => solver.pattern(scramble),
  random: ({ seed }) => solver.randomScramble({ seed }),
  isSolved: ({ moves }) => solver.isSolved(moves),
  bench: (args) => solver.bench(args),
  /** The module's memory: `{ bytes, fixed }` (fixed: the threads build, whose memory never grows). */
  memory: () => ({ bytes: solver.memory.buffer.byteLength, fixed: solver.hasThreads }),
};

onmessage = async ({ data: { id, cmd, ...args } }) => {
  try {
    if (cmd !== 'init' && !solver) throw new Error('not started');
    postMessage({ id, result: await commands[cmd](args) });
  } catch (e) {
    postMessage({ id, error: e.message });
  }
};
