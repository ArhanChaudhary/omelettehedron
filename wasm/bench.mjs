#!/usr/bin/env bun
// Timings of the WebAssembly build in Bun (a browser-like engine: JavaScriptCore, as Safari),
// comparable with the native `bench two-phase`:
//
//   bun wasm/bench.mjs [--tier tiny|small|medium] [--threads N] [--tables DIR] [--count 10]
//                      [--seed 2026] [--seconds 5]
//
// --tier T      the tables: tiny, small (default) or medium
// --threads N   the threads build with N threads (default: the single-thread build)
// --tables DIR  load the tier's files from DIR (ae-pair48-split96.prun, ae-sym288-triples.prun
//               or ae-sym288-chiral.prun, and cb-sym24.prun) instead of generating
//
// Prints the table generation (or load) time and the module's memory, the optimal solves of
// short scrambles (those of test.mjs), then the native `bench two-phase` summary for `count` seeded
// random positions with a time limit of `seconds` each (the same positions as the native
// `bench --tier T -N count --seed seed -t seconds`).

import { FtoSolver, randomMoves, TIERS, wasmFile } from './omelettehedron.mjs';

const args = process.argv.slice(2);
const option = (name, fallback) => {
  const i = args.indexOf(name);
  return i >= 0 ? args[i + 1] : fallback;
};
const threads = Number(option('--threads', 0));
const tier = option('--tier', 'small');
const dist = new URL('./dist/', import.meta.url);
const wasm = new URL(wasmFile(tier, threads ? Math.max(threads, 2) : 1), dist);
const MiB = 2 ** 20;

let t0 = performance.now();
const solver = await FtoSolver.create({ wasm, threads: threads || undefined, onLine: () => {} });
console.log(`module: ${wasm.pathname.split('/').pop()}, ${solver.threads} thread(s), started in `
  + `${((performance.now() - t0) / 1000).toFixed(2)} s`);
t0 = performance.now();
const tables = option('--tables');
const info = tables
  ? await solver.loadTables(`${tables}/${TIERS[tier].ae}`, `${tables}/${TIERS[tier].cb}`)
  : solver.generate({ tier });
console.log(`tables (${info.tier}): ${tables ? 'loaded' : 'generated'} in ${info.seconds.toFixed(2)} s `
  + `(call ${((performance.now() - t0) / 1000).toFixed(2)} s); memory ${solver.memory.buffer.byteLength / MiB} MiB`
  + `${solver.hasThreads ? ' (fixed)' : ' (grown to)'}`);

// Optimal solves of short scrambles (test.mjs's).
let total = 0;
for (const [length, seed] of [[6, 11], [8, 12], [10, 13], [11, 14], [12, 15], [13, 16], [14, 17]]) {
  const scramble = randomMoves(length, seed);
  const r = solver.solve(scramble);
  total += r.seconds;
  console.log(`optimal ${String(r.length).padStart(2)} moves in ${r.seconds.toFixed(3).padStart(7)} s, `
    + `${String(r.nodes).padStart(10)} nodes (${(r.nodes / Math.max(r.seconds, 1e-9) / 1e6).toFixed(1)} M/s): ${scramble}`);
}
console.log(`optimal solves: ${total.toFixed(3)} s in all`);

// The native bench.
const count = Number(option('--count', 10));
const seed = Number(option('--seed', 2026));
const seconds = Number(option('--seconds', 5));
if (count === 0) {
  solver.close();
  process.exit(0);
}
const s = solver.bench({ count, seed, seconds });
console.log(`bench: ${count} random positions (seed ${seed}), ${seconds} s each, ${s.threads} thread(s)`);
for (const [t, mean, solved] of s.checkpoints) {
  console.log(`  best after ${String(t).padStart(5)} s: ${Number.isFinite(mean) && mean !== null ? mean.toFixed(2) : '-'} (${solved}/${count} solved)`);
}
console.log(`  first solution: ${s.mean_first_time?.toFixed(3)} s mean, ${s.mean_first_length?.toFixed(2)} moves mean; `
  + `max ${Math.max(...s.results.map((r) => r.first_time ?? Infinity)).toFixed(3)} s`);
console.log(`  mean length ${s.mean_length?.toFixed(2)}, ${(s.nodes_per_sec / 1e6).toFixed(1)} M nodes/s, `
  + `${s.solved_optimally} proven optimal`);
console.log(`memory after: ${solver.memory.buffer.byteLength / MiB} MiB`);
solver.close();
process.exit(0);
