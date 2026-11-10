#!/usr/bin/env bash
# Optimal-mode benchmark: time to prove optimality on seeded random-move scrambles.
# usage: scripts/bench-opt.sh <label> [tier]   -> bench-results/<label>-opt{mt,st}.json
set -euo pipefail
label=${1:?label}
tier=${2:-large}
cd "$(dirname "$0")/.."
cargo build --release -q
bin=target/release/fto
mkdir -p bench-results
echo "load before benchmark: $(cut -d' ' -f1 /proc/loadavg)" >&2
$bin bench -N 10 --scramble-len 19 --seed 1919 --tier "$tier" --json "bench-results/$label-optmt.json" > "bench-results/$label-optmt.txt"
$bin bench -N 10 --scramble-len 17 --seed 1717 -j 1 --tier "$tier" --json "bench-results/$label-optst.json" > "bench-results/$label-optst.txt"
grep -E "elapsed|optimal|final" "bench-results/$label-optmt.txt" "bench-results/$label-optst.txt"
