#!/bin/bash
# Solve the hard positions of bench/hard.tsv with the two-phase solver:
#   scripts/bench-hard.sh <tier> [solve two-phase options...]
# e.g. scripts/bench-hard.sh large -t 10        (best solution within 10 s each)
#      scripts/bench-hard.sh large               (optimal, with target-length passes; can take
#                                                 hours for 23-move positions)
# Prints name, known length, the solver's output. Takes the shared benchmark lock per position.
set -euo pipefail
tier=${1:?tier (small|medium|large|huge)}
shift
cd "$(dirname "$0")/.."
bin=${CARGO_TARGET_DIR:-target}/release/omelettehedron
# without a limit, the proof with target-length passes (with one, the best solution by then)
mode=(--passes)
for a in "$@"; do
  case $a in -t* | --time* | -l* | --max-length* | -n* | --max-solutions*) mode=() ;; esac
done
grep -v '^#' bench/hard.tsv | while IFS=$'\t' read -r name known config scramble solution; do
  echo "== $name (known $known, $config)"
  flock /tmp/cube-bench.lock "$bin" solve two-phase "${mode[@]}" --tier "$tier" -v "$@" "$scramble" 2>&1 |
    grep -vE "^loaded|^found|^  depth .*split done" || true
done
