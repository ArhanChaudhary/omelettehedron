// Run a WASI preview-1 module under bun: bun wasi-run.mjs module.wasm [args...]
import { WASI } from "node:wasi";
import { readFileSync } from "node:fs";
const [file, ...args] = process.argv.slice(2);
const wasi = new WASI({ version: "preview1", args: [file, ...args], env: process.env, preopens: { "/": "/" } });
const mod = await WebAssembly.compile(readFileSync(file));
const inst = await WebAssembly.instantiate(mod, { wasi_snapshot_preview1: wasi.wasiImport });
const code = wasi.start(inst);
process.exit(code ?? 0);
