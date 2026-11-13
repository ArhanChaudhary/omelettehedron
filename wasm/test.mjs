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
//  - an optimal solve of a random position, cancelled after a second, returns at once;
//  - the threads build for the small tables refuses to generate the medium ones (no room in
//    its fixed 1 GiB), with an error instead of a trap;
//  - the tiny tables: the threads build generates them (with FTO_TABLES holding the native
//    tiny tier, byte for byte its files), the single-thread build loads them, and both solve
//    the same scrambles with the native optimal lengths (and on one thread node counts).
// The medium tables (1.04 GB), if FTO_TABLES has the native medium tier
// (`gen two-phase --tier medium`):
//  - the single-thread build loads the native files and solves the same scrambles with the
//    native optimal lengths and one-thread node counts (the medium index on 32-bit wasm);
//  - the threads build for the medium tables (2 GiB, `omelettehedron-threads-medium.wasm`)
//    loads them too and solves them; a solve is cancelled.
// With FTO_TEST_MEDIUM_GENERATE=1 the threads build for the medium tables also generates them
// (slow: minutes) and exports them; with FTO_TABLES the files must equal the native ones.
//
//   bun wasm/test.mjs        (exit status 0 if all passed)

import { cancel, invertMoves, parseImprovement, TIERS } from './omelettehedron.mjs';

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
// The fourth value is the node count with the medium tables (`--tier medium`, same lengths),
// the fifth with the tiny ones (`--tier tiny`).
const CASES = [
  ["F' BL' D' F' L' BR'", 6, 62, 62, 170],
  ["B BR U' BR L B' D F'", 8, 56, 56, 110],
  ["BR' U' R' D B' BL' U B U F'", 10, 58, 58, 58],
  ["BL D L B U' BL L' F' D R L'", 11, 486, 74, 616],
  ["D' L' F BL BR' F B L' D' U' F L", 12, 214, 82, 278],
  ["R' B BR BL' U F' R' L U' BR' B' D F'", 13, 2000, 316, 4806],
  ["U R' B R' U BL' BR' R B BL' U L BL' U", 14, 169652, 19846, 315594],
];

/** The cases as [scramble, length, nodes] for the tables `tier`. */
const casesFor = (cases, tier) =>
  cases.map(([scramble, length, small, medium, tiny]) =>
    [scramble, length, { small, medium, tiny }[tier]]);

/** Recompute the native values with FTO_NATIVE and FTO_TABLES, if set. */
function nativeCases() {
  const native = process.env.FTO_NATIVE;
  const tables = process.env.FTO_TABLES;
  if (!native || !tables) return CASES;
  const run = (tier, scramble) => {
    const p = Bun.spawnSync([native, 'solve', '--tier', tier, '--tables', tables, '-j', '1', '-v', scramble], {
      env: { ...process.env, FTO_DEBUG: '1' },
    });
    const text = p.stdout.toString() + p.stderr.toString();
    const length = Number(/^length (\d+) in/m.exec(text)?.[1]);
    const nodes = [...text.matchAll(/depth \d+: \d+ tasks, (\d+) nodes/g)].map((m) => Number(m[1])).pop();
    return [length, nodes];
  };
  const medium = Bun.file(`${tables}/${TIERS.medium.ae}`).size > 0;
  const tiny = Bun.file(`${tables}/${TIERS.tiny.ae}`).size > 0;
  return CASES.map(([scramble]) => {
    const [length, nodes] = run('small', scramble);
    return [scramble, length, nodes, medium ? run('medium', scramble)[1] : null,
      tiny ? run('tiny', scramble)[1] : null];
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
    const out = w.lines.filter((l) => l.fd === 1).map((l) => l.line);
    const improvements = out.map(parseImprovement).filter(Boolean);
    check(improvements.length > 0 && improvements.at(-1)?.solution === r.solution,
      `${label}: ${scramble}: improvements ${JSON.stringify(improvements)} do not end in ${r.solution}`);
    // proven lower bounds (`<s>s  lower <n>`) never exceed the optimum
    const lowers = out.map((l) => /^\s*[\d.]+s\s+lower\s+(\d+)\s*$/.exec(l)).filter(Boolean).map((m) => Number(m[1]));
    check(lowers.every((n) => n <= r.length), `${label}: ${scramble}: lower bounds ${lowers} above ${r.length}`);
    // keep improving (no time limit): better and better solutions, ending proven optimal
    const k = await w.call('solve', { scramble, seconds: -1 });
    check(k.optimal && k.length === r.length && k.verified,
      `${label}: ${scramble}: keep improving gave ${JSON.stringify(k)}, optimal mode ${r.length}`);
    console.log(`${label}: ${scramble}: ${r.length} moves (${r.optimal ? 'optimal' : 'not optimal'}), `
      + `${r.nodes} nodes, ${(r.seconds * 1000).toFixed(1)} ms: ${r.solution}`);
  }
}

/** Whether two Blobs have the same bytes (compared in pieces). */
async function sameBytes(a, b) {
  if (a.size !== b.size) return false;
  const piece = 1 << 24;
  for (let at = 0; at < a.size; at += piece) {
    const [x, y] = await Promise.all([a, b].map((f) => f.slice(at, at + piece).bytes()));
    if (Buffer.compare(x, y) !== 0) return false;
  }
  return true;
}

/**
 * An optimal solve of a random position, cancelled after `after` ms from this thread; as a
 * scramble, then as a pattern (a stopped pattern solve once reported about 10,000,000 s).
 */
async function cancelled(w, control, label, after = 1000) {
  const { scramble } = await w.call('random', { seed: 7 });
  const { pattern } = await w.call('pattern', { scramble });
  for (const [kind, call, args] of [['solve', 'solve', { scramble }], ['pattern', 'solvePattern', { pattern }]]) {
    await cancelledOne(w, control, `${label} (${kind})`, after, call, args);
  }
}

async function cancelledOne(w, control, label, after, call, args) {
  const t0 = performance.now();
  const solving = w.call(call, args);
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
await solveCases(mt, `threads build, ${threads} threads`, casesFor(cases, 'small'), false);
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
await solveCases(st, 'single build', casesFor(cases, 'small'), true);
const random = await st.call('random', { seed: 1 });
const moves = random.scramble.split(' ');
check(moves.length === random.length && random.length >= 15 && random.length <= 40,
  `single build: random scramble ${JSON.stringify(random)}`);
check(!(await st.call('isSolved', { moves: random.scramble })).solved, 'single build: random scramble is solved');
// The page's Invert button: a sequence followed by its inverse leaves the puzzle solved, for
// every suffix and layer of the notation (and a solution inverted is a scramble it solves).
for (const seq of [random.scramble, "U R' BL2 F2' D3 B1 L1' BR", "Uw Rs' BLv F2 Dw2' Us3 BRv'", 'F’ U']) {
  const inv = invertMoves(seq);
  check((await st.call('isSolved', { moves: `${seq} ${inv}` })).solved,
    `single build: ${seq} followed by its inverse ${inv} is not solved`);
  check(invertMoves(inv).split(' ').length === inv.split(' ').length, `single build: inverse of ${inv}`);
}
check(invertMoves("U R' BL") === "BL' R U'", `single build: invertMoves gave ${invertMoves("U R' BL")}`);
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

// The threads build for the small tables has no room for the medium ones: an error, no trap.
let noRoom = '';
try {
  await mt.call('generate', { tier: 'medium' });
} catch (e) {
  noRoom = e.message;
}
check(/not enough memory/.test(noRoom), `threads build: generating the medium tables: ${noRoom || 'no error'}`);
console.log(`threads build: generating the medium tables is refused: ${noRoom}`);

// --- The tiny tables (the same threads module as the small ones) ---
const MiB = 2 ** 20;
t0 = performance.now();
mt.lines.length = 0;
const genT = await mt.call('generate', { tier: 'tiny' });
await sleep(50);
const progressT = mt.lines.filter((l) => l.fd === 2 && /depth \d+/.test(l.line));
check(genT.tier === 'tiny' && genT.ae_bytes === 37260304 && progressT.length >= 10,
  `threads build: generated the tiny tables ${JSON.stringify(genT)}, ${progressT.length} progress lines`);
console.log(`threads build: generated the tiny tables in ${genT.seconds.toFixed(2)} s (call ${seconds(t0)} s), `
  + `AE ${genT.ae_bytes} + CB ${genT.cb_bytes} bytes, ${progressT.length} progress lines`);
const filesT = await mt.call('export');
check(filesT[0].name === TIERS.tiny.ae, `threads build: exported ${filesT[0].name}`);
if (process.env.FTO_TABLES && (await Bun.file(`${process.env.FTO_TABLES}/${TIERS.tiny.ae}`).exists())) {
  for (const f of filesT) {
    const same = await sameBytes(Bun.file(`${process.env.FTO_TABLES}/${f.name}`), f.blob);
    check(same, `exported ${f.name} (tiny) differs from the native file`);
    console.log(`exported ${f.name} (${f.blob.size} bytes) ${same ? 'equals' : 'DIFFERS FROM'} the native file`);
  }
}
await solveCases(mt, `threads build, ${threads} threads, tiny`, casesFor(cases, 'tiny'), false);
mt.close();
st.close();
const stT = demoWorker('single-tiny');
await stT.call('init', { threads: 1, channel: stT.channelName, tier: 'tiny' });
t0 = performance.now();
const loadedT = await stT.call('load', { ae: filesT[0].blob, cb: filesT[1].blob });
const grownT = await stT.call('memory');
check(loadedT.tier === 'tiny', `single build: loaded ${JSON.stringify(loadedT)}`);
console.log(`single build: loaded the tiny tables in ${loadedT.seconds.toFixed(2)} s (call ${seconds(t0)} s); `
  + `memory grown to ${grownT.bytes / MiB} MiB`);
await solveCases(stT, 'single build, tiny', casesFor(cases, 'tiny'), true);
stT.close();

// --- The medium tables ---
const tablesDir = process.env.FTO_TABLES;
const nativeMedium = tablesDir && (await Bun.file(`${tablesDir}/${TIERS.medium.ae}`).exists());
const generateMedium = process.env.FTO_TEST_MEDIUM_GENERATE === '1';
if (!nativeMedium && !generateMedium) {
  console.log('medium tables: not tested (FTO_TABLES without the medium tier, and FTO_TEST_MEDIUM_GENERATE unset)');
} else {
  const mediumCases = casesFor(cases, 'medium');
  const nativeFiles = nativeMedium
    ? { ae: `${tablesDir}/${TIERS.medium.ae}`, cb: `${tablesDir}/${TIERS.medium.cb}` }
    : null;

  // The threads build for the medium tables: its memory, generation (if asked for), loading.
  t0 = performance.now();
  const mm = demoWorker('threads-medium');
  const initM = await mm.call('init', { threads, channel: mm.channelName, tier: 'medium' });
  check(initM.threads === threads && initM.control, `medium threads build: init ${JSON.stringify(initM)}`);
  const memory = await mm.call('memory');
  check(memory.bytes === 2 * 1024 * MiB && memory.fixed, `medium threads build: memory ${JSON.stringify(memory)}`);
  console.log(`medium threads build: ${initM.threads} threads, ${memory.bytes / MiB} MiB fixed memory, `
    + `started in ${seconds(t0)} s`);
  let medium = nativeFiles;
  if (generateMedium) {
    t0 = performance.now();
    mm.lines.length = 0;
    const g = await mm.call('generate', { tier: 'medium' });
    await sleep(50);
    const lines = mm.lines.filter((l) => l.fd === 2).map((l) => l.line);
    check(g.tier === 'medium' && g.ae_bytes === 1031923200, `medium threads build: generated ${JSON.stringify(g)}`);
    console.log(`medium threads build: generated in ${g.seconds.toFixed(2)} s (call ${seconds(t0)} s), `
      + `AE ${g.ae_bytes} + CB ${g.cb_bytes} bytes; progress:\n  ${lines.join('\n  ')}`);
    const exported = await mm.call('export');
    check(exported[0].name === TIERS.medium.ae, `medium threads build: exported ${exported[0].name}`);
    if (nativeFiles) {
      for (const [i, f] of exported.entries()) {
        const same = await sameBytes(Bun.file(`${tablesDir}/${f.name}`), f.blob);
        check(same, `exported ${f.name} (medium) differs from the native file`);
        console.log(`exported ${f.name} (${f.blob.size} bytes) ${same ? 'equals' : 'DIFFERS FROM'} the native file`);
      }
    }
    medium = { ae: exported[0].blob, cb: exported[1].blob };
  }
  t0 = performance.now();
  const loadedM = await mm.call('load', medium);
  check(loadedM.tier === 'medium', `medium threads build: loaded ${JSON.stringify(loadedM)}`);
  console.log(`medium threads build: loaded the ${generateMedium ? 'exported' : 'native'} medium tables in `
    + `${loadedM.seconds.toFixed(2)} s (call ${seconds(t0)} s)`);
  await solveCases(mm, `medium threads build, ${threads} threads`, mediumCases, false);
  await cancelled(mm, initM.control, `medium threads build, ${threads} threads`);
  mm.close();

  // The single-thread build: the same files; one-thread node counts as natively.
  const sm = demoWorker('single-medium');
  await sm.call('init', { threads: 1, channel: sm.channelName, tier: 'medium' });
  t0 = performance.now();
  const loaded1 = await sm.call('load', medium);
  const grown = await sm.call('memory');
  check(loaded1.tier === 'medium' && loaded1.ae_bytes === 1031923200, `single build: loaded ${JSON.stringify(loaded1)}`);
  console.log(`single build: loaded the medium tables in ${loaded1.seconds.toFixed(2)} s (call ${seconds(t0)} s); `
    + `memory grown to ${grown.bytes / MiB} MiB`);
  await solveCases(sm, 'single build, medium', mediumCases, true);
  sm.close();
}

if (failures.length) {
  console.error(`wasm test: ${failures.length} failure(s)`);
  process.exit(1);
}
console.log('wasm test: all passed');
process.exit(0);
