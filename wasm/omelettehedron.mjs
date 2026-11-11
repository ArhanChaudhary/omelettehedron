// omelettehedron in WebAssembly: a JavaScript wrapper for browsers and Bun.
//
//   import { FtoSolver } from './omelettehedron.mjs';
//   const solver = await FtoSolver.create({ wasm: 'dist/omelettehedron.wasm' });
//   solver.generate();                          // the small tier's tables, in memory
//   solver.solve("U R' BL F D B' L BR");        // optimal: { solution, length, optimal, ... }
//   solver.solve(scramble, { seconds: 5 });     // the best solution found in 5 s
//
// Three sets of tables (`TIERS`, the native tiers of the same names): tiny (42 MB; quick to
// make, slower solves), small (134 MB) and medium (1.04 GB; fewer positions to search, so
// faster solves). The single-thread module takes any; the threads build comes in two modules
// (`wasmFile`: tiny and small share one), since its memory is fixed (see wasm/build.rs).
//
// The module is built for WASI (wasm32-wasip1, or wasm32-wasip1-threads for several
// threads): this file provides the few WASI functions it imports (clock, console output,
// random numbers, sleep), and for the threads build the `thread-spawn` of wasi-threads with a
// pool of workers (`omelettehedron-thread.mjs`) on a shared memory. In browsers the threads
// build needs a cross-origin isolated page (COOP/COEP headers, see `serve.mjs`), and either
// build should run in a worker: the solver blocks until it is done. See README.md,
// "WebAssembly demo".
//
// Output lines of the module (stdout: each better solution while solving, if asked for;
// stderr: progress while generating tables) go to a BroadcastChannel if `channel` is given,
// from every thread (the pool's threads have their own workers, whose event loops are blocked
// while they work; the channel delivers at once), else to the console. `onLine(fd, line)`
// takes the lines of the calling thread instead (those of the pool's threads still go to the
// channel or the console).

const ERRNO_SUCCESS = 0;
const ERRNO_AGAIN = 6;
const ERRNO_BADF = 8;
const ERRNO_NOSYS = 52;

// The control array (Int32Array, shared with the pool's workers): word 0 is the last thread
// id, word 1 the number of cancellations so far (see `cancel`); thread pool slot i is words
// 4 + 4i .. 8 + 4i: state, thread id, argument.
const CANCEL_WORD = 1;
const SLOT_STARTING = 0;
const SLOT_IDLE = 1;
const SLOT_CLAIMED = 2;
const SLOT_RUN = 3;
const slotBase = (i) => 4 + 4 * i;

// How far each cancellation moves the module's monotonic clock: 10^7 s, more than any time
// limit (the module's optimal solves have 10^6 s), and few enough to keep the clock in u64
// nanoseconds for a thousand cancellations.
const CANCEL_JUMP_NS = 10n ** 16n;

/** A copy of wasm memory bytes (TextDecoder, Blob and others refuse shared views). */
function copyOut(memory, ptr, len) {
  const copy = new Uint8Array(len);
  copy.set(new Uint8Array(memory.buffer, ptr, len));
  return copy;
}

/** Where output lines go: `onLine`, else the BroadcastChannel `channel`, else the console. */
export function lineSink({ onLine, channel } = {}) {
  if (onLine) return onLine;
  if (channel && typeof BroadcastChannel !== 'undefined') {
    const bc = new BroadcastChannel(channel);
    return (fd, line) => bc.postMessage({ fd, line });
  }
  return (fd, line) => (fd === 1 ? console.log : console.error)(line);
}

/**
 * The WASI preview 1 functions the module imports (others answer ENOSYS): no files, no
 * arguments or environment; stdout and stderr go to `sink(fd, line)`, line by line. The
 * monotonic clock runs `CANCEL_JUMP_NS` ahead per cancellation counted in `control` (see
 * `cancel`). (Addresses arrive as signed i32: `>>> 0` makes those from 2 GiB up positive.)
 */
export function wasiImports(getMemory, module, { sink, control } = {}) {
  const lines = { 1: '', 2: '' };
  const out = sink ?? lineSink();
  const view = () => new DataView(getMemory().buffer);
  const implemented = {
    args_sizes_get(countPtr, sizePtr) {
      [countPtr, sizePtr] = [countPtr >>> 0, sizePtr >>> 0];
      view().setUint32(countPtr, 0, true);
      view().setUint32(sizePtr, 0, true);
      return ERRNO_SUCCESS;
    },
    args_get: () => ERRNO_SUCCESS,
    environ_sizes_get(countPtr, sizePtr) {
      [countPtr, sizePtr] = [countPtr >>> 0, sizePtr >>> 0];
      view().setUint32(countPtr, 0, true);
      view().setUint32(sizePtr, 0, true);
      return ERRNO_SUCCESS;
    },
    environ_get: () => ERRNO_SUCCESS,
    clock_time_get(id, _precision, outPtr) {
      outPtr >>>= 0;
      // Realtime (0) in ns since 1970, else monotonic: the same epoch-based clock in every
      // thread (each worker has its own `performance.now()` origin), moved on by
      // cancellations.
      let ns;
      if (id === 0) {
        ns = BigInt(Date.now()) * 1000000n;
      } else {
        ns = BigInt(Math.round((performance.timeOrigin + performance.now()) * 1e6));
        const cancels = control ? Atomics.load(control, CANCEL_WORD) : 0;
        ns += BigInt(cancels) * CANCEL_JUMP_NS;
      }
      view().setBigUint64(outPtr, ns, true);
      return ERRNO_SUCCESS;
    },
    clock_res_get(_id, outPtr) {
      outPtr >>>= 0;
      view().setBigUint64(outPtr, 1000n, true);
      return ERRNO_SUCCESS;
    },
    fd_write(fd, iovs, iovsLen, writtenPtr) {
      [iovs, writtenPtr] = [iovs >>> 0, writtenPtr >>> 0];
      if (fd !== 1 && fd !== 2) return ERRNO_BADF;
      const memory = getMemory();
      let written = 0;
      let text = '';
      const decoder = new TextDecoder();
      for (let i = 0; i < iovsLen; i++) {
        const ptr = view().getUint32(iovs + 8 * i, true);
        const len = view().getUint32(iovs + 8 * i + 4, true);
        text += decoder.decode(copyOut(memory, ptr, len));
        written += len;
      }
      const all = lines[fd] + text;
      const parts = all.split('\n');
      lines[fd] = parts.pop();
      for (const line of parts) out(fd, line);
      view().setUint32(writtenPtr, written, true);
      return ERRNO_SUCCESS;
    },
    fd_fdstat_get(fd, statPtr) {
      statPtr >>>= 0;
      if (fd > 2) return ERRNO_BADF;
      const v = view();
      v.setUint8(statPtr, 2); // character device
      v.setUint16(statPtr + 2, 0, true);
      v.setBigUint64(statPtr + 8, 0n, true);
      v.setBigUint64(statPtr + 16, 0n, true);
      return ERRNO_SUCCESS;
    },
    fd_prestat_get: () => ERRNO_BADF,
    fd_close: () => ERRNO_SUCCESS,
    // Sleeps (std::thread::sleep): clock subscriptions only.
    poll_oneoff(inPtr, outPtr, count, eventsPtr) {
      [inPtr, outPtr, eventsPtr] = [inPtr >>> 0, outPtr >>> 0, eventsPtr >>> 0];
      const v = view();
      let timeoutMs = 0;
      for (let i = 0; i < count; i++) {
        const s = inPtr + 48 * i;
        if (v.getUint8(s + 8) === 0) {
          const ns = v.getBigUint64(s + 24, true);
          timeoutMs = Math.max(timeoutMs, Number(ns) / 1e6);
        }
      }
      sleep(timeoutMs);
      for (let i = 0; i < count; i++) {
        const s = inPtr + 48 * i;
        const e = outPtr + 32 * i;
        v.setBigUint64(e, v.getBigUint64(s, true), true);
        v.setUint16(e + 8, ERRNO_SUCCESS, true);
        v.setUint8(e + 10, v.getUint8(s + 8));
      }
      v.setUint32(eventsPtr, count, true);
      return ERRNO_SUCCESS;
    },
    sched_yield: () => ERRNO_SUCCESS,
    random_get(ptr, len) {
      [ptr, len] = [ptr >>> 0, len >>> 0];
      const bytes = new Uint8Array(len);
      for (let i = 0; i < len; i += 65536) {
        crypto.getRandomValues(bytes.subarray(i, Math.min(len, i + 65536)));
      }
      new Uint8Array(getMemory().buffer, ptr, len).set(bytes);
      return ERRNO_SUCCESS;
    },
    proc_exit(code) {
      throw new Error(`omelettehedron: exit ${code}`);
    },
  };
  const imports = {};
  for (const { module: from, name } of WebAssembly.Module.imports(module)) {
    if (from === 'wasi_snapshot_preview1') {
      imports[name] = implemented[name] ?? (() => ERRNO_NOSYS);
    }
  }
  return imports;
}

/** Block for `ms` (Atomics.wait where allowed, else a busy wait). */
function sleep(ms) {
  if (ms <= 0) return;
  try {
    Atomics.wait(new Int32Array(new SharedArrayBuffer(4)), 0, 0, ms);
  } catch {
    const end = performance.now() + ms;
    while (performance.now() < end);
  }
}

/**
 * The `thread-spawn` import of wasi-threads: hands the thread to an idle pool worker. The
 * module starts its threads once and never ends them (rayon's global pool), so the pool has
 * room for all of them. It must not wait for a worker: the spawning thread holds wasi-libc's
 * thread list lock.
 */
export function threadSpawn(control) {
  const slots = (control.length - 4) / 4;
  return (startArg) => {
    for (let i = 0; i < slots; i++) {
      const base = slotBase(i);
      if (Atomics.compareExchange(control, base, SLOT_IDLE, SLOT_CLAIMED) === SLOT_IDLE) {
        const tid = Atomics.add(control, 0, 1) + 1;
        Atomics.store(control, base + 1, tid);
        Atomics.store(control, base + 2, startArg);
        Atomics.store(control, base, SLOT_RUN);
        Atomics.notify(control, base);
        return tid;
      }
    }
    console.error('omelettehedron: no idle worker for a new thread');
    return -ERRNO_AGAIN;
  };
}

/**
 * Cancel the solve in progress: moves the module's monotonic clock past every time limit, so
 * the search stops at its next time check (every 16384 nodes per thread) and returns the best
 * solution so far. (Solves started later measure their time from the moved clock.) Call it from any thread that
 * has the solver's `control` array, e.g. the page while a worker solves; it works only if
 * `control` is shared (a SharedArrayBuffer: Bun, or a cross-origin isolated page). Returns
 * whether it could.
 */
export function cancel(control) {
  const shared = typeof SharedArrayBuffer !== 'undefined' && control?.buffer instanceof SharedArrayBuffer;
  if (!shared) return false;
  Atomics.add(control, CANCEL_WORD, 1);
  return true;
}

/** The imported memory's limits (pages), read from the module's import section. */
function importedMemoryLimits(bytes) {
  let at = 8;
  const leb = () => {
    let result = 0;
    let shift = 0;
    for (;;) {
      const b = bytes[at++];
      result += (b & 0x7f) * 2 ** shift;
      shift += 7;
      if (b < 0x80) return result;
    }
  };
  const skipName = () => {
    const len = leb();
    const name = new TextDecoder().decode(bytes.subarray(at, at + len));
    at += len;
    return name;
  };
  while (at < bytes.length) {
    const id = bytes[at++];
    const size = leb();
    const end = at + size;
    if (id === 2) {
      for (let n = leb(); n > 0; n--) {
        skipName();
        skipName();
        const kind = bytes[at++];
        if (kind === 0) leb();
        else if (kind === 1) {
          at++;
          const flags = bytes[at++];
          leb();
          if (flags & 1) leb();
        } else if (kind === 2) {
          const flags = bytes[at++];
          const initial = leb();
          const maximum = flags & 1 ? leb() : undefined;
          return { initial, maximum, shared: (flags & 2) !== 0 };
        } else if (kind === 3) at += 2;
        else if (kind === 4) {
          at++;
          leb();
        }
      }
    }
    at = end;
  }
  return undefined;
}

/**
 * A local file as a `Blob` (Bun: a path, or a `file:` URL), else undefined (browsers, and
 * http(s), data and blob URLs, which are fetched).
 */
function localFile(source) {
  if (typeof Bun === 'undefined' || !(typeof source === 'string' || source instanceof URL)) {
    return undefined;
  }
  return /^(https?|data|blob):/.test(String(source)) ? undefined : Bun.file(source);
}

/** The bytes of `source`: a path (Bun), a URL, a Blob, or bytes already. */
async function readSource(source) {
  if (source instanceof Uint8Array) return source;
  if (source instanceof ArrayBuffer) return new Uint8Array(source);
  if (typeof Blob !== 'undefined' && source instanceof Blob) {
    return new Uint8Array(await source.arrayBuffer());
  }
  const file = localFile(source);
  if (file) return new Uint8Array(await file.arrayBuffer());
  const response = await fetch(source);
  if (!response.ok) throw new Error(`${source}: ${response.status} ${response.statusText}`);
  return new Uint8Array(await response.arrayBuffer());
}

/** Start a pool worker (`omelettehedron-thread.mjs`) and wait until it is idle. */
async function startWorker(module, memory, control, slot, channel) {
  const url = new URL('./omelettehedron-thread.mjs', import.meta.url);
  const worker = new Worker(url, { type: 'module' });
  await new Promise((resolve, reject) => {
    worker.onmessage = (e) => (e.data === 'ready' ? resolve() : reject(new Error(e.data)));
    worker.onerror = (e) => reject(new Error(e.message));
    worker.postMessage({ module, memory, control, slot, channel });
  });
  // (Bun: idle pool workers must not keep the process alive.)
  worker.unref?.();
  return worker;
}

/**
 * The sets of tables: the module's number for each, its two files (the native `gen --tier`
 * files) and the threads module with room for it.
 */
export const TIERS = {
  tiny: { code: 2, ae: 'ae-pair48-split96.prun', cb: 'cb-sym24.prun', threadsWasm: 'omelettehedron-threads.wasm' },
  small: { code: 0, ae: 'ae-sym288-triples.prun', cb: 'cb-sym24.prun', threadsWasm: 'omelettehedron-threads.wasm' },
  medium: { code: 1, ae: 'ae-sym288-chiral.prun', cb: 'cb-sym24.prun', threadsWasm: 'omelettehedron-threads-medium.wasm' },
};

/** The module file (in `dist/`) for the tables `tier` on `threads` threads (> 1: a threads build). */
export function wasmFile(tier = 'small', threads = 1) {
  const t = TIERS[tier];
  if (!t) throw new Error(`omelettehedron: no tables called ${tier} (tiny, small or medium)`);
  return threads > 1 ? t.threadsWasm : 'omelettehedron.wasm';
}

/** `source` as a Blob, if it is one or (Bun) a local file; else undefined. */
async function blobOf(source) {
  if (typeof Blob !== 'undefined' && source instanceof Blob) return source;
  const file = localFile(source);
  if (file && !(await file.exists())) throw new Error(`${source}: no such file`);
  return file;
}

/** A random 32-bit seed. */
const randomSeed = () => crypto.getRandomValues(new Uint32Array(1))[0];

/** The solver: one WebAssembly instance with its tables. */
export class FtoSolver {
  /**
   * Instantiate the module and start its threads. `wasm`: a path (Bun), URL, bytes or a
   * compiled `WebAssembly.Module`; `threads`: threads for the threads build (default: all
   * cores; the single-thread build uses 1). `onLine`, `channel`: where the module's output
   * lines go (see the top of this file).
   */
  static async create({ wasm, threads, onLine, channel } = {}) {
    const bytes = wasm instanceof WebAssembly.Module ? null : await readSource(wasm);
    const module = bytes ? await WebAssembly.compile(bytes) : wasm;
    const needsThreads = WebAssembly.Module.imports(module).some((i) => i.module === 'wasi');
    const shared = typeof SharedArrayBuffer !== 'undefined';
    let memory;
    const imports = {};
    if (needsThreads) {
      if (!shared) {
        throw new Error('the threads build needs SharedArrayBuffer (a cross-origin isolated page)');
      }
      // (The threads build imports all its memory at once, see wasm/build.rs.)
      const limits = bytes ? importedMemoryLimits(bytes) : { initial: 16384, maximum: 16384 };
      memory = new WebAssembly.Memory({ ...limits, shared: true });
      threads = Math.max(1, threads ?? globalThis.navigator?.hardwareConcurrency ?? 1);
    } else {
      threads = 1;
    }
    // One pool worker per thread the module starts (threads - 1), plus one.
    const slots = needsThreads ? threads : 0;
    const control = new Int32Array(shared ? new SharedArrayBuffer(4 * slotBase(slots)) : new ArrayBuffer(4 * slotBase(slots)));
    if (needsThreads) {
      imports.env = { memory };
      imports.wasi = { 'thread-spawn': threadSpawn(control) };
    }
    let instanceMemory = memory;
    const sink = lineSink({ onLine, channel });
    imports.wasi_snapshot_preview1 = wasiImports(() => instanceMemory, module, { sink, control });
    const instance = await WebAssembly.instantiate(module, imports);
    instanceMemory = memory ?? instance.exports.memory;
    instance.exports._initialize?.();
    const workers = [];
    try {
      for (let slot = 0; slot < slots; slot++) {
        workers.push(await startWorker(module, memory, control, slot, channel));
      }
    } catch (e) {
      for (const worker of workers) worker.terminate();
      throw e;
    }
    const solver = new FtoSolver(instance, instanceMemory, control, workers);
    solver.threads = solver.call('fto_init', threads).threads;
    return solver;
  }

  constructor(instance, memory, control, workers) {
    this.exports = instance.exports;
    this.memory = memory;
    /** The control array: pass it to `cancel` (from another thread) to stop a solve. */
    this.control = control;
    this.workers = workers;
    /** Threads the solver uses (1 for the single-thread build). */
    this.threads = 1;
    /** `{ tier, ae_bytes, cb_bytes, seconds }` of the tables, or null before there are any. */
    this.tables = null;
  }

  /** Whether this is the threads build. */
  get hasThreads() {
    return this.exports.fto_has_threads() === 1;
  }

  /** The output of the last call (JSON or an error message). */
  output() {
    const ptr = this.exports.fto_output() >>> 0;
    const len = this.exports.fto_output_len() >>> 0;
    return new TextDecoder().decode(copyOut(this.memory, ptr, len));
  }

  /** Run an export; returns its parsed JSON output or throws its error. */
  call(name, ...args) {
    let status;
    try {
      status = this.exports[name](...args);
    } catch (e) {
      // A Rust panic or an allocation failure aborts the instance.
      this.broken = true;
      throw new Error(`omelettehedron: ${name} failed (${e.message}); create a new instance`);
    }
    const text = this.output();
    if (status !== 0) throw new Error(`omelettehedron: ${text}`);
    return JSON.parse(text);
  }

  /** Pass a string to `f(ptr, len)` in a module buffer. */
  withText(text, f) {
    const bytes = new TextEncoder().encode(text);
    const ptr = this.exports.fto_alloc(bytes.length) >>> 0;
    if (ptr === 0) throw new Error('omelettehedron: out of memory');
    try {
      new Uint8Array(this.memory.buffer, ptr, bytes.length).set(bytes);
      return f(ptr, bytes.length);
    } finally {
      this.exports.fto_free(ptr, bytes.length);
    }
  }

  /**
   * Generate the tables of `tier` ('tiny', 'small' or 'medium', see `TIERS`) in memory:
   * natively on one core ~5 s for tiny, ~30 s for small and ~3 min for medium; more threads:
   * faster. The threads module
   * must be the one for the tier (`wasmFile`). Progress lines go to stderr (`onLine(2,
   * line)`). Returns `{ tier, ae_bytes, cb_bytes, seconds }`.
   */
  generate({ tier = 'small' } = {}) {
    const t = TIERS[tier];
    if (!t) throw new Error(`omelettehedron: no tables called ${tier} (tiny, small or medium)`);
    this.tables = null;
    this.tables = this.call('fto_generate', t.code);
    return this.tables;
  }

  /**
   * Load the two table files of a tier (`TIERS[tier].ae` and `cb-sym24.prun`, from the native
   * `gen --tier tiny`, `small` or `medium`, or from `exportTables`): each a path (Bun), URL,
   * Blob or bytes. The tier is that of the first file; the checksums are checked. Blobs and
   * Bun's files are copied into the module in pieces (no second whole copy in JavaScript's
   * memory). Returns `{ tier, ae_bytes, cb_bytes, seconds }`.
   */
  async loadTables(ae, cb) {
    this.tables = null;
    for (const [which, source] of [ae, cb].entries()) {
      const blob = await blobOf(source);
      const bytes = blob ? null : await readSource(source);
      const len = blob ? blob.size : bytes.length;
      const ptr = this.exports.fto_table_buffer(which, len) >>> 0;
      if (ptr === 0) throw new Error(`omelettehedron: no room for a ${len}-byte table file`);
      if (bytes) {
        new Uint8Array(this.memory.buffer, ptr, len).set(bytes);
        continue;
      }
      const piece = 1 << 24;
      for (let at = 0; at < len; at += piece) {
        const part = new Uint8Array(await blob.slice(at, Math.min(len, at + piece)).arrayBuffer());
        new Uint8Array(this.memory.buffer, ptr + at, part.length).set(part);
      }
    }
    this.tables = this.call('fto_load');
    return this.tables;
  }

  /**
   * The tables as their two files (the format of `loadTables` and the native `gen`), as
   * `[{ name, blob }, { name, blob }]`: e.g. to keep them in IndexedDB.
   */
  exportTables() {
    return [0, 1].map((which) => {
      const e = this.call('fto_export', which);
      const parts = [copyOut(this.memory, e.header_ptr, e.header_len)];
      // (copied in pieces: no second copy of the whole table at once)
      const piece = 1 << 24;
      for (let at = 0; at < e.data_len; at += piece) {
        parts.push(copyOut(this.memory, e.data_ptr + at, Math.min(piece, e.data_len - at)));
      }
      return { name: e.name, blob: new Blob(parts) };
    });
  }

  /**
   * Solve a scramble. `seconds` > 0: return the best solution found in that time; 0 (the
   * default): run until it is proven optimal. `maxLength` > 0: stop at the first solution of
   * at most that many moves. `improvements`: print each better solution as it is found
   * (stdout: `<seconds>s  <length>  <moves>`, see `parseImprovement`). Returns `{ solution,
   * length, optimal, verified, seconds, nodes, solutions, events, cancelled }` (`solution`
   * null if none was found in time; `verified`: it solves the scramble on the sticker model;
   * `cancelled`: stopped by `cancel`).
   */
  solve(scramble, { seconds = 0, maxLength = 0, improvements = false } = {}) {
    return this.timed(() => this.withText(scramble, (ptr, len) =>
      this.call('fto_solve', ptr, len, seconds, maxLength, improvements ? 1 : 0)));
  }

  /**
   * Run a solve (`fto_solve` or `fto_solve_pattern`) and set `cancelled`. A cancellation moves
   * the module's clock on by `CANCEL_JUMP_NS` (that is how the search sees its time run out),
   * so the times of a cancelled solve are measured here instead: without this a stopped solve
   * reported about 10,000,000 s.
   */
  timed(run) {
    const cancels = Atomics.load(this.control, CANCEL_WORD);
    const start = performance.now();
    const r = run();
    r.cancelled = Atomics.load(this.control, CANCEL_WORD) !== cancels;
    if (r.cancelled && typeof r.seconds === 'number') {
      r.seconds = (performance.now() - start) / 1000;
      r.events = (r.events ?? []).map(([t, length]) => [Math.min(t, r.seconds), length]);
    }
    return r;
  }

  /**
   * A random-state scramble: a random position (from `seed`), solved up to its first
   * solution, which is inverted. Returns `{ scramble, length }`.
   */
  randomScramble({ seed = randomSeed() } = {}) {
    return this.call('fto_random_scramble', seed >>> 0);
  }

  /** Whether `moves` (e.g. a scramble and its solution) leave the puzzle solved: `{ solved, moves }`. */
  isSolved(moves) {
    return this.withText(moves, (ptr, len) => this.call('fto_is_solved', ptr, len));
  }

  /** The pattern (72 sticker colours, src/pattern.rs) that a scramble makes: `{ pattern }`. */
  pattern(scramble) {
    return this.withText(scramble, (ptr, len) => this.call('fto_pattern', ptr, len));
  }

  /**
   * Solve a pattern (72 sticker colours) as `solve` solves a scramble; `usual` is the colours
   * of U R F L B BL D BR on the solved puzzle as usually held (8 letters), for a scramble from
   * that hold. Returns `solve`'s result plus `scramble: { start, moves, usual, exact }`, or
   * `{ error: { message, kind, stickers, counts } }` for a wrong pattern.
   */
  solvePattern(pattern, { usual = '', seconds = 0, maxLength = 0, improvements = false } = {}) {
    return this.timed(() => this.withText(pattern, (ptr, len) =>
      this.withText(usual, (uptr, ulen) =>
        this.call('fto_solve_pattern', ptr, len, uptr, ulen, seconds, maxLength, improvements ? 1 : 0))));
  }

  /**
   * The native `bench`: `count` seeded random positions (or random-move scrambles of
   * `scrambleLength` moves) solved one at a time for `seconds` each (0: optimal). Returns the
   * native `--json` summary.
   */
  bench({ count = 10, seed = 2026, seconds = 1, scrambleLength = 0 } = {}) {
    return this.call('fto_bench', count, seed, seconds, scrambleLength);
  }

  /**
   * Have the engine compile the search with its optimizing compiler now, so that the first
   * real solve gets it: V8 (Chrome, Node) runs new code with its baseline compiler and
   * switches a function to optimized code only on its next call, so a long first solve would
   * run slowly to its end. Solves short scrambles optimally for about `seconds`, then gives
   * the background compiler time to finish. Not in Bun: JavaScriptCore moves a running loop
   * to optimized code.
   */
  async warmUp({ seconds = 0.3 } = {}) {
    if (typeof Bun !== 'undefined') return;
    const start = performance.now();
    let seed = 1;
    while (performance.now() - start < seconds * 1000) {
      this.solve(randomMoves(8, seed++));
    }
    await new Promise((resolve) => setTimeout(resolve, 300));
    this.solve(randomMoves(8, seed));
  }

  /** Stop the thread pool's workers (the instance cannot be used afterwards). */
  close() {
    for (const worker of this.workers) worker.terminate();
    this.workers = [];
  }
}

/** The 8 faces (in the order of the crate's FACE_NAMES). */
export const FACES = ['U', 'R', 'L', 'F', 'B', 'BR', 'BL', 'D'];

/** A seeded random-move scramble of `length` face turns, no face twice in a row. */
export function randomMoves(length, seed) {
  let x = (seed >>> 0) || 1;
  const next = () => {
    // xorshift32
    x ^= x << 13;
    x >>>= 0;
    x ^= x >>> 17;
    x ^= x << 5;
    x >>>= 0;
    return x;
  };
  const moves = [];
  let last = -1;
  while (moves.length < length) {
    const face = next() % 8;
    if (face === last) continue;
    last = face;
    moves.push(FACES[face] + (next() & 1 ? "'" : ''));
  }
  return moves.join(' ');
}

/**
 * The inverse of a move sequence, as text: the moves in reverse order, each turned the other
 * way. A solution inverted is a scramble of the position it solves. It follows the notation
 * of `parse_moves` (src/notation.rs): a face name, an optional `w` / `s` / `v`, then nothing,
 * `1`, `'`, `1'`, `3`, `2` or `2'`; a token it does not know is kept as it is (the solver then
 * reports it). Runs without the module, so it works while tables are being made.
 */
export function invertMoves(text) {
  const inverse = { '': "'", 1: "'", "'": '', "1'": '', 3: '', 2: "2'", "2'": '2' };
  return text
    .replace(/[‘’‛′´`ʼʹ]/g, "'")
    .split(/[\s,]+/)
    .filter((t) => t)
    .reverse()
    .map((t) => {
      const m = /^([A-Za-z]+)((?:[123]'?)|'|)$/.exec(t);
      return m && m[2] in inverse ? m[1] + inverse[m[2]] : t;
    })
    .join(' ');
}

/** An improvement line printed while solving (`{ seconds, length, solution }`), or null. */
export function parseImprovement(line) {
  const m = /^\s*([\d.]+)s\s+(\d+)\s+(.*)$/.exec(line);
  return m ? { seconds: Number(m[1]), length: Number(m[2]), solution: m[3].trim() } : null;
}

export { slotBase, SLOT_STARTING, SLOT_IDLE, SLOT_CLAIMED, SLOT_RUN, CANCEL_WORD };
