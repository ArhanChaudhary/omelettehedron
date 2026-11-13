#!/usr/bin/env bash
# Optimal mode on random positions (the representative case): one position per seed, 48 threads,
# each capped at a time limit. Slow (minutes per position), so it is for headline numbers and
# optimal-mode changes rather than every A/B comparison.
# usage: scripts/bench-randopt.sh <label> <tier> [seeds="1 2 3"] [cap_seconds=1200]
#   -> bench-results/<label>-randopt-seed<S>.json and a summary on stdout
set -euo pipefail
label=${1:?label}
tier=${2:?tier (small|medium|large|huge)}
seeds=${3:-1 2 3}
cap=${4:-1200}
cd "$(dirname "$0")/.."
mkdir -p bench-results
bin=${CARGO_TARGET_DIR:-target}/release/omelettehedron
for seed in $seeds; do
  # one lock turn per position, so other jobs can run in between
  flock /tmp/cube-bench.lock sh -c "cargo build --release -q && \
    echo \"commit \$(git rev-parse --short HEAD), seed $seed, load \$(cut -d' ' -f1 /proc/loadavg)\" >&2 && \
    FTO_DEBUG=1 $bin bench two-phase -N 1 --seed $seed -t $cap --tier $tier -v \
      --json bench-results/$label-randopt-seed$seed.json" 2>&1 | grep -E "^#|commit|depth" || true
done
python3 - "$label" $seeds <<'PY'
import json, sys
label, seeds = sys.argv[1], sys.argv[2:]
print(f"| seed | best | proven | time (s) | nodes |")
print("|---|---|---|---|---|")
for s in seeds:
    r = json.load(open(f"bench-results/{label}-randopt-seed{s}.json"))["results"][0]
    print(f"| {s} | {r['length']} | {'yes' if r['optimal'] else 'no'} | {r['elapsed']:.1f} | {r['nodes']:.3e} |")
PY
