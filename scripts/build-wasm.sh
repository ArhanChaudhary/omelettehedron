#!/usr/bin/env bash
# Build the WebAssembly modules of the browser demo into wasm/dist/:
#   omelettehedron.wasm                 one thread (wasm32-wasip1): runs anywhere, either set of
#                                       tables (its memory grows)
#   omelettehedron-threads.wasm         threads (wasm32-wasip1-threads): Bun, or a cross-origin
#                                       isolated page; 1 GiB of memory, the small tables
#   omelettehedron-threads-medium.wasm  the same with 2 GiB, for the medium tables
# Built with the Rust version in wasm/rust-version, the one the GitHub Pages site is built with
# (.github/workflows/pages.yml), so that a local test checks the modules the site ships;
# RUSTUP_TOOLCHAIN overrides it. Needs that toolchain with the targets:
#   rustup toolchain install "$(cat wasm/rust-version)" --profile minimal \
#     --target wasm32-wasip1 --target wasm32-wasip1-threads
# The build is reproducible from Cargo.lock (--locked) with the same toolchain; SIMD128 is
# always on (every current browser, Bun and Node have it). The threads builds' fixed memory
# sizes are set in wasm/build.rs.
# usage: scripts/build-wasm.sh   (CARGO_TARGET_DIR is honoured)
set -euo pipefail
cd "$(dirname "$0")/.."
target_dir=${CARGO_TARGET_DIR:-target}
export RUSTFLAGS="-C target-feature=+simd128"
export RUSTUP_TOOLCHAIN=${RUSTUP_TOOLCHAIN:-$(cat wasm/rust-version)}
if ! rustc --version > /dev/null 2>&1; then
    echo "missing Rust $RUSTUP_TOOLCHAIN (wasm/rust-version): rustup toolchain install" \
        "$RUSTUP_TOOLCHAIN --profile minimal --target wasm32-wasip1 --target wasm32-wasip1-threads" >&2
    exit 1
fi
for target in wasm32-wasip1 wasm32-wasip1-threads; do
    if ! rustup target list --installed | grep -qx "$target"; then
        echo "missing Rust target $target: rustup target add --toolchain $RUSTUP_TOOLCHAIN $target" >&2
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
