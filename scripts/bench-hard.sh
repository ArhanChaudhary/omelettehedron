#!/bin/bash
# Solve the hard positions of bench/hard.tsv: scripts/bench-hard.sh <tier> [solve options...]
# e.g. scripts/bench-hard.sh large -t 10        (best solution within 10 s each)
#      scripts/bench-hard.sh large               (optimal; can take hours for 23-move positions)
# Prints name, known length, the solver's result line. Takes the shared benchmark lock per position.
set -euo pipefail
tier=${1:?tier (small|medium|large|huge)}
shift
cd "$(dirname "$0")/.."
bin=${CARGO_TARGET_DIR:-target}/release/omelettehedron
grep -v '^#' bench/hard.tsv | while IFS=$'\t' read -r name known config scramble; do
  echo "== $name (known $known, $config)"
  flock /tmp/cube-bench.lock "$bin" solve --tier "$tier" -v "$@" "$scramble" 2>&1 | grep -vE "^loaded|^  depth .*split done" || true
done
