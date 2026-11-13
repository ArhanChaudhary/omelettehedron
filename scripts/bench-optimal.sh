#!/usr/bin/env bash
# Optimal solver benchmark: time (and nodes, lookups) to solve seeded random-move scrambles
# optimally with `omelettehedron bench optimal`: 10 x 17-move scrambles on 1 thread and 10 x
# 19-move scrambles on all threads. The two-phase solver's optimal-mode benchmark is
# scripts/bench-opt.sh.
# usage: scripts/bench-optimal.sh <label> <tier> [moves]
#   -> bench-results/<label>-optimal{1t,mt}.{json,txt}; compare with
#      scripts/compare-optimal.py <before> <after>
set -euo pipefail
label=${1:?label}
tier=${2:?tier (small|medium)}
moves=${3:-}
cd "$(dirname "$0")/.."
# Serialize builds and runs with other benchmarks on this machine (the lock is shared with
# other projects' benchmarks), so no build or other benchmark overlaps a measurement.
exec 9>/tmp/cube-bench.lock
flock 9
cargo build --release -q
bin=${CARGO_TARGET_DIR:-target}/release/omelettehedron
mv=()
[ -n "$moves" ] && mv=(--moves "$moves")
mkdir -p bench-results
echo "commit $(git rev-parse --short HEAD), load before benchmark: $(cut -d' ' -f1 /proc/loadavg)" >&2
$bin bench optimal --tier "$tier" "${mv[@]}" -N 10 --scramble-len 17 --seed 1717 -j 1 \
    --json "bench-results/$label-optimal1t.json" > "bench-results/$label-optimal1t.txt"
$bin bench optimal --tier "$tier" "${mv[@]}" -N 10 --scramble-len 19 --seed 1919 \
    --json "bench-results/$label-optimalmt.json" > "bench-results/$label-optimalmt.txt"
grep -E "solved|mean length|total nodes|total lookups|total time|nodes/s" \
    "bench-results/$label-optimal1t.txt" "bench-results/$label-optimalmt.txt"
