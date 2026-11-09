#!/usr/bin/env bash
# Standard benchmark set used to compare optimisations.
# usage: scripts/bench.sh <label> [tier]
# Writes bench-results/<label>-{mt,st}.json. Positions are seeded random states.
set -euo pipefail
label=${1:?label}
tier=${2:-medium}
cd "$(dirname "$0")/.."
cargo build --release -q
bin=target/release/fto
mkdir -p bench-results
load=$(cut -d' ' -f1 /proc/loadavg)
echo "load before benchmark: $load" >&2
# multi-threaded: 40 positions x 5 s
$bin bench -N 40 -t 5 --seed 2026 --tier "$tier" --json "bench-results/$label-mt.json" > "bench-results/$label-mt.txt"
# single-threaded: 20 positions x 5 s
$bin bench -N 20 -t 5 --seed 7 -j 1 --tier "$tier" --json "bench-results/$label-st.json" > "bench-results/$label-st.txt"
cat "bench-results/$label-mt.txt" "bench-results/$label-st.txt"
