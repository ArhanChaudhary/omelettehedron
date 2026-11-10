#!/usr/bin/env bash
# Build the WebAssembly modules of the browser demo into wasm/dist/:
#   omelettehedron.wasm          one thread (wasm32-wasip1): runs anywhere
#   omelettehedron-threads.wasm  threads (wasm32-wasip1-threads): Bun, or a cross-origin
#                                isolated page
# Needs the targets: rustup target add wasm32-wasip1 wasm32-wasip1-threads
# The build is reproducible from Cargo.lock (--locked) with the same toolchain; SIMD128 is
# always on (every current browser, Bun and Node have it). The threads build's fixed memory
# size is set in wasm/build.rs.
# usage: scripts/build-wasm.sh   (CARGO_TARGET_DIR is honoured)
# (Adapted from egg-solver's wasm/build.sh, same author.)
set -euo pipefail
cd "$(dirname "$0")/.."
target_dir=${CARGO_TARGET_DIR:-target}
export RUSTFLAGS="-C target-feature=+simd128"
for target in wasm32-wasip1 wasm32-wasip1-threads; do
    if ! rustup target list --installed | grep -qx "$target"; then
        echo "missing Rust target $target: rustup target add $target" >&2
        exit 1
    fi
done
mkdir -p wasm/dist
for target in wasm32-wasip1 wasm32-wasip1-threads; do
    cargo build --locked -p omelettehedron-wasm --profile wasm-release --target "$target"
done
cp "$target_dir/wasm32-wasip1/wasm-release/omelettehedron_wasm.wasm" wasm/dist/omelettehedron.wasm
cp "$target_dir/wasm32-wasip1-threads/wasm-release/omelettehedron_wasm.wasm" \
    wasm/dist/omelettehedron-threads.wasm
ls -l wasm/dist/*.wasm
