// A thread of the threads build (wasi-threads): a pool worker of `omelettehedron.mjs`. It
// instantiates the module on the shared memory, then waits in its control slot;
// `thread-spawn` (in any thread) hands it a thread, which it runs with `wasi_thread_start`.
// The module's threads (rayon's global pool) never end, so in practice each worker runs one.
// Its output lines go to the BroadcastChannel `channel` if given (this worker's own event
// loop is blocked while the thread runs), else to the console.
//
// Adapted from egg-solver's wasm/egg-solver-thread.mjs (same author).

import { wasiImports, threadSpawn, lineSink, slotBase, SLOT_IDLE, SLOT_RUN } from './omelettehedron.mjs';

const post = (message) => self.postMessage(message);

/**
 * A thread that traps (a Rust panic aborts) never tells the threads waiting for it: they
 * would wait forever. Say so, and in Bun end the process (a page stays stuck: its solver
 * worker may be the one waiting).
 */
function fatal(message) {
  console.error(message);
  globalThis.process?.kill(globalThis.process.pid, 'SIGTERM');
}

async function serve({ module, memory, control, slot, channel }) {
  let instance;
  try {
    instance = await WebAssembly.instantiate(module, {
      env: { memory },
      wasi: { 'thread-spawn': threadSpawn(control) },
      wasi_snapshot_preview1: wasiImports(() => memory, module, {
        sink: lineSink({ channel }),
        control,
      }),
    });
  } catch (e) {
    post(`omelettehedron thread: ${e.message}`);
    return;
  }
  const base = slotBase(slot);
  Atomics.store(control, base, SLOT_IDLE);
  post('ready');
  for (;;) {
    let state;
    while ((state = Atomics.load(control, base)) !== SLOT_RUN) {
      Atomics.wait(control, base, state);
    }
    const tid = Atomics.load(control, base + 1);
    const startArg = Atomics.load(control, base + 2);
    try {
      instance.exports.wasi_thread_start(tid, startArg);
    } catch (e) {
      fatal(`omelettehedron: thread ${tid} failed (${e.message}); the module cannot go on`);
      return;
    }
    Atomics.store(control, base, SLOT_IDLE);
  }
}

self.onmessage = (e) => serve(e.data);
