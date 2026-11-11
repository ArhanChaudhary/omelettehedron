#!/usr/bin/env bun
// Headless test of the built modules (scripts/build-wasm.sh) in Bun, driven through the
// demo's worker (demo-worker.mjs) as the page drives it:
//  - the threads build (FTO_TEST_THREADS threads, default 4) generates the small tier while
//    its progress lines arrive on the BroadcastChannel, and exports the tables; with
//    FTO_TABLES set (a directory with the native small tier) the files must equal the
//    native ones byte for byte;
//  - the single-thread build loads the exported files;
//  - both solve short scrambles optimally: every solution must solve its scramble on the
//    sticker model (`isSolved`, apart from the solver's own check), have the native optimal
//    length and, on one thread, take exactly the native number of nodes (the search is
//    deterministic); each better solution must arrive on the channel while solving. The
//    native values are recorded below; FTO_NATIVE=<binary> with FTO_TABLES recomputes them;
//  - a random-state scramble is well formed;
//  - an optimal solve of a random position, cancelled after a second, returns at once.
//
//   bun wasm/test.mjs        (exit status 0 if all passed)

import { cancel, parseImprovement } from './omelettehedron.mjs';

const threads = Number(process.env.FTO_TEST_THREADS ?? 4);
const failures = [];
const check = (ok, message) => {
  if (!ok) {
    failures.push(message);
    console.error(`FAIL: ${message}`);
  }
};
const seconds = (t0) => ((performance.now() - t0) / 1000).toFixed(2);
const sleep = (ms) => new Promise((resolve) => setTimeout(resolve, ms));

// Short random-move scrambles (`randomMoves` of omelettehedron.mjs: lengths 6..14, seeds
// 11..17) with the native optimal length and one-thread node count (FTO_DEBUG=1
// omelettehedron solve --tier small -j 1 -v, after fix/cancel: the module's optimal solves have
// no time limit, only a cancellation deadline, so they keep the one-thread target-length passes).
const CASES = [
  ["F' BL' D' F' L' BR'", 6, 62],
  ["B BR U' BR L B' D F'", 8, 56],
  ["BR' U' R' D B' BL' U B U F'", 10, 58],
  ["BL D L B U' BL L' F' D R L'", 11, 486],
  ["D' L' F BL BR' F B L' D' U' F L", 12, 214],
  ["R' B BR BL' U F' R' L U' BR' B' D F'", 13, 2000],
  ["U R' B R' U BL' BR' R B BL' U L BL' U", 14, 169652],
];

/** Recompute the native values with FTO_NATIVE and FTO_TABLES, if set. */
function nativeCases() {
  const native = process.env.FTO_NATIVE;
  const tables = process.env.FTO_TABLES;
  if (!native || !tables) return CASES;
  return CASES.map(([scramble]) => {
    const p = Bun.spawnSync([native, 'solve', '--tier', 'small', '--tables', tables, '-j', '1', '-v', scramble], {
      env: { ...process.env, FTO_DEBUG: '1' },
    });
    const text = p.stdout.toString() + p.stderr.toString();
    const length = Number(/^length (\d+) in/m.exec(text)?.[1]);
    const nodes = [...text.matchAll(/depth \d+: \d+ tasks, (\d+) nodes/g)].map((m) => Number(m[1])).pop();
    return [scramble, length, nodes];
  });
}

/** The demo worker, as the page starts it; `call(cmd, args)` and the channel's lines. */
function demoWorker(name) {
  const worker = new Worker(new URL('./demo-worker.mjs', import.meta.url), { type: 'module' });
  const pending = new Map();
  let nextId = 0;
  worker.onmessage = ({ data }) => {
    const p = pending.get(data.id);
    pending.delete(data.id);
    'error' in data ? p.reject(new Error(data.error)) : p.resolve(data.result);
  };
  const channelName = `omelettehedron-test-${name}-${Math.random().toString(36).slice(2)}`;
  const lines = [];
  const channel = new BroadcastChannel(channelName);
  channel.onmessage = ({ data }) => lines.push(data);
  const call = (cmd, args = {}) => new Promise((resolve, reject) => {
    const id = nextId++;
    pending.set(id, { resolve, reject });
    worker.postMessage({ id, cmd, ...args });
  });
  const close = () => {
    worker.terminate();
    channel.close();
  };
  return { call, lines, channelName, close };
}

/** Solve every case optimally on `w` and check the results. */
async function solveCases(w, label, cases, exactNodes) {
  for (const [scramble, length, nodes] of cases) {
    w.lines.length = 0;
    const r = await w.call('solve', { scramble });
    await sleep(50); // (the channel's last lines)
    const ok = r.solution !== null && (await w.call('isSolved', { moves: `${scramble} ${r.solution}` })).solved;
    check(r.solution !== null && ok, `${label}: ${r.solution} does not solve ${scramble}`);
    check(r.verified, `${label}: the module's own check failed for ${scramble}`);
    check(r.optimal, `${label}: ${scramble} not proven optimal`);
    if (length !== null) check(r.length === length, `${label}: ${scramble}: length ${r.length}, native ${length}`);
    if (exactNodes && nodes !== null) check(r.nodes === nodes, `${label}: ${scramble}: ${r.nodes} nodes, native ${nodes}`);
    const improvements = w.lines.filter((l) => l.fd === 1).map((l) => parseImprovement(l.line));
    check(improvements.length > 0 && improvements.at(-1)?.solution === r.solution,
      `${label}: ${scramble}: improvements ${JSON.stringify(improvements)} do not end in ${r.solution}`);
    console.log(`${label}: ${scramble}: ${r.length} moves (${r.optimal ? 'optimal' : 'not optimal'}), `
      + `${r.nodes} nodes, ${(r.seconds * 1000).toFixed(1)} ms: ${r.solution}`);
  }
}

/** An optimal solve of a random position, cancelled after `after` ms from this thread. */
async function cancelled(w, control, label, after = 1000) {
  const { scramble } = await w.call('random', { seed: 7 });
  const t0 = performance.now();
  const solving = w.call('solve', { scramble });
  await sleep(after);
  check(cancel(control), `${label}: cancel refused`);
  const t1 = performance.now();
  const r = await solving;
  const late = (performance.now() - t1) / 1000;
  check(r.cancelled && !r.optimal, `${label}: cancelled solve: ${JSON.stringify(r)}`);
  check(late < 2, `${label}: the cancelled solve took ${late.toFixed(2)} s to return`);
  check(r.seconds > after / 1000 - 0.1 && r.seconds < (performance.now() - t0) / 1000 + 0.1,
    `${label}: cancelled solve reports ${r.seconds} s`);
  console.log(`${label}: cancelled an optimal solve after ${after} ms: returned ${late.toFixed(3)} s later, `
    + `best ${r.length ?? 'none'}, ${r.nodes} nodes`);
}

const cases = nativeCases();

// The threads build: generate, export, compare with the native files.
let t0 = performance.now();
const mt = demoWorker('threads');
const init = await mt.call('init', { threads, channel: mt.channelName });
check(init.threads === threads && init.control, `threads build: init ${JSON.stringify(init)}`);
console.log(`threads build: ${init.threads} threads, started in ${seconds(t0)} s`);
t0 = performance.now();
const gen = await mt.call('generate');
await sleep(50);
const progress = mt.lines.filter((l) => l.fd === 2 && /depth \d+/.test(l.line));
check(progress.length >= 15, `threads build: ${progress.length} progress lines`);
console.log(`threads build: generated in ${gen.seconds.toFixed(2)} s (call ${seconds(t0)} s), `
  + `AE ${gen.ae_bytes} + CB ${gen.cb_bytes} bytes, ${progress.length} progress lines`);
const files = await mt.call('export');
const bytes = await Promise.all(files.map((f) => f.blob.bytes()));
if (process.env.FTO_TABLES) {
  for (const [i, f] of files.entries()) {
    const native = await Bun.file(`${process.env.FTO_TABLES}/${f.name}`).bytes();
    check(Buffer.compare(native, bytes[i]) === 0, `exported ${f.name} differs from the native file`);
    console.log(`exported ${f.name} (${bytes[i].length} bytes) equals the native file`);
  }
}
await solveCases(mt, `threads build, ${threads} threads`, cases, false);
await cancelled(mt, init.control, `threads build, ${threads} threads`);

// The single-thread build: load the exported files.
const st = demoWorker('single');
const init1 = await st.call('init', { threads: 1, channel: st.channelName });
check(init1.threads === 1, `single build: init ${JSON.stringify(init1)}`);
t0 = performance.now();
const loaded = await st.call('load', { ae: files[0].blob, cb: files[1].blob });
console.log(`single build: loaded the exported tables in ${loaded.seconds.toFixed(2)} s (call ${seconds(t0)} s)`);
let refused = false;
try {
  const damaged = bytes[1].slice();
  damaged[damaged.length - 1] ^= 1;
  await st.call('load', { ae: files[0].blob, cb: damaged });
} catch (e) {
  refused = /checksum/.test(e.message);
}
check(refused, 'single build: a damaged table file was not refused');
await st.call('load', { ae: files[0].blob, cb: files[1].blob });
await solveCases(st, 'single build', cases, true);
const random = await st.call('random', { seed: 1 });
const moves = random.scramble.split(' ');
check(moves.length === random.length && random.length >= 15 && random.length <= 40,
  `single build: random scramble ${JSON.stringify(random)}`);
check(!(await st.call('isSolved', { moves: random.scramble })).solved, 'single build: random scramble is solved');
console.log(`single build: random-state scramble (${random.length} moves): ${random.scramble}`);
let bad = false;
try {
  await st.call('solve', { scramble: 'U R X' });
} catch (e) {
  bad = /unknown move/.test(e.message);
}
check(bad, 'single build: a bad scramble was accepted');
// Patterns: a scramble's pattern solved as a pattern (the same optimal length), with the
// scramble from the usual hold making the pattern again; and a pattern with a wrong sticker.
for (const [scramble, length] of cases.slice(0, 3)) {
  const { pattern } = await st.call('pattern', { scramble });
  const r = await st.call('solvePattern', { pattern, usual: 'WRGPBOYA' });
  check(r.verified && r.optimal && (length === null || r.length === length),
    `single build: pattern of ${scramble}: ${JSON.stringify(r)}`);
  check(r.scramble?.exact, `single build: no exact scramble from the usual hold: ${JSON.stringify(r.scramble)}`);
  const again = await st.call('pattern', { scramble: r.scramble?.usual ?? '' });
  check(again.pattern === pattern, `single build: ${r.scramble?.usual} makes ${again.pattern}, not ${pattern}`);
  console.log(`single build: pattern of ${scramble}: ${r.solution} (${r.length}), scramble ${r.scramble?.usual}`);
}
const solvedPattern = (await st.call('pattern', { scramble: '' })).pattern;
const wrong = await st.call('solvePattern', { pattern: solvedPattern.replace('W', 'R') });
check(wrong.error?.kind === 'counts' && wrong.error.counts.length === 2,
  `single build: a wrong sticker gave ${JSON.stringify(wrong)}`);
await cancelled(st, init1.control, 'single build');

mt.close();
st.close();
if (failures.length) {
  console.error(`wasm test: ${failures.length} failure(s)`);
  process.exit(1);
}
console.log('wasm test: all passed');
process.exit(0);
