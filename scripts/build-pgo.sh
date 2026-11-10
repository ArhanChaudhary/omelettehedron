#!/bin/bash
# Profile-guided build (x86-64): an instrumented build solves a training set, and the release
# build then uses the recorded profile. The training positions use seeds other than those of the
# standard benchmarks (scripts/bench.sh, scripts/bench-opt.sh).
# usage: scripts/build-pgo.sh <tier>   (the tier's tables must exist or are generated first)
#   -> ${CARGO_TARGET_DIR:-target}/pgo/release/omelettehedron
# Needs llvm-profdata of the active toolchain: rustup component add llvm-tools
set -euo pipefail
tier=${1:?tier (small|medium|large|huge)}
cd "$(dirname "$0")/.."
target=${CARGO_TARGET_DIR:-target}
host=$(rustc -vV | sed -n 's/^host: //p')
profdata="$(rustc --print sysroot)/lib/rustlib/$host/bin/llvm-profdata"
[ -x "$profdata" ] || { echo "llvm-profdata not found: rustup component add llvm-tools" >&2; exit 1; }
base="-C target-cpu=x86-64-v3"
prof=$target/pgo-data
rm -rf "$prof"
RUSTFLAGS="$base -Cprofile-generate=$prof/raw" CARGO_TARGET_DIR=$target/pgo-gen cargo build --release -q
gen=$target/pgo-gen/release/omelettehedron
$gen bench --tier "$tier" -N 10 -t 2 --seed 11 > /dev/null
$gen bench --tier "$tier" -N 4 -t 2 --seed 12 -j 1 > /dev/null
$gen bench --tier "$tier" -N 6 --scramble-len 18 --seed 13 > /dev/null
$gen bench --tier "$tier" -N 6 --scramble-len 16 --seed 14 -j 1 > /dev/null
"$profdata" merge -o "$prof/merged.profdata" "$prof/raw"
RUSTFLAGS="$base -Cprofile-use=$(realpath "$prof/merged.profdata")" CARGO_TARGET_DIR=$target/pgo \
  cargo build --release -q
echo "$target/pgo/release/omelettehedron"
