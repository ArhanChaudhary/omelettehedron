#!/usr/bin/env bash
# Cross-check of optimal mode on deep positions: every position is solved with the early end
# of the optimality proof (the default) and with the full proof (FTO_FULL_PROOF=1); the proven
# lengths must be identical. Run it for every change that touches optimal mode: the early end
# is only valid while both variants complete every iteration (see search.rs, `solve`), and the
# brute-force tests only reach about 10 moves.
# usage: scripts/check-optimal.sh <tier>   (about 2 minutes at 48 threads with the large tier)
set -euo pipefail
tier=${1:?tier}
cd "$(dirname "$0")/.."
exec 9>/tmp/cube-bench.lock
flock 9
cargo build --release -q
bin=${CARGO_TARGET_DIR:-target}/release/fto
out=$(mktemp -d)
run() { # name, extra bench arguments...
    local name=$1; shift
    "$bin" bench --tier "$tier" "$@" --json "$out/$name-early.json" > /dev/null
    FTO_FULL_PROOF=1 "$bin" bench --tier "$tier" "$@" --json "$out/$name-full.json" > /dev/null
}
run scr19 -N 10 --scramble-len 19 --seed 1919   # 19-move scrambles (optimal about 18)
run rand3 -N 1 --seed 3                         # a random position (optimal 20)
python3 - "$out" <<'PY'
import json, sys, os
out = sys.argv[1]
bad = 0
for name in ("scr19", "rand3"):
    early = json.load(open(os.path.join(out, f"{name}-early.json")))["results"]
    full = json.load(open(os.path.join(out, f"{name}-full.json")))["results"]
    for e, f in zip(early, full):
        ok = e["optimal"] and f["optimal"] and e["length"] == f["length"]
        bad += not ok
        print(f"{name} #{e['index']}: early {e['length']} ({e['elapsed']:.2f}s)  "
              f"full {f['length']} ({f['elapsed']:.2f}s)  {'ok' if ok else 'MISMATCH'}")
print("all proven lengths agree" if not bad else f"{bad} MISMATCHES")
sys.exit(1 if bad else 0)
PY
rm -rf "$out"
