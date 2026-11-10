#!/usr/bin/env bash
# Optimal-mode benchmark: time to prove optimality on seeded random-move scrambles.
# usage: scripts/bench-opt.sh <label> [tier]   -> bench-results/<label>-opt{mt,st}.json
set -euo pipefail
label=${1:?label}
tier=${2:?tier (small|medium|large|huge)}
cd "$(dirname "$0")/.."
# Serialize builds and runs with other benchmarks on this machine (lock shared with the
# egg-solver sessions), so no build or other benchmark overlaps a measurement.
exec 9>/tmp/cube-bench.lock
flock 9
cargo build --release -q
bin=${CARGO_TARGET_DIR:-target}/release/fto
mkdir -p bench-results
echo "commit $(git rev-parse --short HEAD), FTO_DIRECT=${FTO_DIRECT:-}, load before benchmark: $(cut -d' ' -f1 /proc/loadavg)" >&2
$bin bench -N 10 --scramble-len 19 --seed 1919 --tier "$tier" --json "bench-results/$label-optmt.json" > "bench-results/$label-optmt.txt"
$bin bench -N 10 --scramble-len 17 --seed 1717 -j 1 --tier "$tier" --json "bench-results/$label-optst.json" > "bench-results/$label-optst.txt"
grep -E "elapsed|optimal|final" "bench-results/$label-optmt.txt" "bench-results/$label-optst.txt"
