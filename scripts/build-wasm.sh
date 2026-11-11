#!/usr/bin/env bash
# Build the WebAssembly modules of the browser demo into wasm/dist/:
#   omelettehedron.wasm                 one thread (wasm32-wasip1): runs anywhere, either set of
#                                       tables (its memory grows)
#   omelettehedron-threads.wasm         threads (wasm32-wasip1-threads): Bun, or a cross-origin
#                                       isolated page; 1 GiB of memory, the small tables
#   omelettehedron-threads-medium.wasm  the same with 2 GiB, for the medium tables
# Needs the targets: rustup target add wasm32-wasip1 wasm32-wasip1-threads
# The build is reproducible from Cargo.lock (--locked) with the same toolchain; SIMD128 is
# always on (every current browser, Bun and Node have it). The threads builds' fixed memory
# sizes are set in wasm/build.rs.
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
# (each module is copied as soon as it is built: the two threads builds have the same output)
build() {
    local target=$1 output=$2
    shift 2
    cargo build --locked -p omelettehedron-wasm --profile wasm-release --target "$target" "$@"
    cp "$target_dir/$target/wasm-release/omelettehedron_wasm.wasm" "wasm/dist/$output"
}
build wasm32-wasip1 omelettehedron.wasm
build wasm32-wasip1-threads omelettehedron-threads.wasm
build wasm32-wasip1-threads omelettehedron-threads-medium.wasm --features medium
ls -l wasm/dist/*.wasm
