#!/usr/bin/env bash
# Cross-check of optimal mode on deep positions: every position is solved with the early end
# of the optimality proof (the default) and with the full proof (FTO_FULL_PROOF=1); the proven
# lengths must be identical. Run it for every change that touches optimal mode: the early end
# is only valid while both variants complete every iteration (see search.rs, `solve`), and the
# brute-force tests only reach about 10 moves. A one-thread set checks the target-length passes
# (on by default with one thread, off with FTO_FULL_PROOF) the same way.
#
# With a move set that has slices (`--moves slice`) the early end is off (a solving path can end
# in a slice, which neither variant finds earlier), so without target passes (FTO_TARGET_GAP=0)
# both runs must do the same search: on one thread their node counts must be identical too. The
# proven lengths are also checked against the meet-in-the-middle brute force in the slice metric
# (the slice tests of tests/solver.rs, including the ignored 400-position one).
# usage: scripts/check-optimal.sh <tier> [moves]   (moves as for --moves, e.g. slice)
#   (about 2 minutes at 48 threads with the large tier, with or without slices)
set -euo pipefail
tier=${1:?tier}
moves=${2:-}
cd "$(dirname "$0")/.."
exec 9>/tmp/cube-bench.lock
flock 9
cargo build --release -q
bin=${CARGO_TARGET_DIR:-target}/release/omelettehedron
out=$(mktemp -d)
mv=()
[ -n "$moves" ] && mv=(--moves "$moves")
gap=()
run() { # name, extra bench arguments... (with the environment in `gap` for both runs)
    local name=$1; shift
    env "${gap[@]}" "$bin" bench two-phase --tier "$tier" "${mv[@]}" "$@" \
        --json "$out/$name-early.json" > /dev/null
    env "${gap[@]}" FTO_FULL_PROOF=1 "$bin" bench two-phase --tier "$tier" "${mv[@]}" "$@" \
        --json "$out/$name-full.json" > /dev/null
}
same_nodes=()
if [[ "$moves" == *slice* ]]; then
    run scr17 -N 10 --scramble-len 17 --seed 1717   # 17-move scrambles (optimal about 15)
    run scr15st -N 10 --scramble-len 15 --seed 1515 -j 1   # one thread: target passes
    gap=(FTO_TARGET_GAP=0)
    run scr13st -N 10 --scramble-len 13 --seed 1313 -j 1   # one thread: identical searches
    sets=(scr17 scr15st scr13st)
    same_nodes=(scr13st)
else
    run scr19 -N 10 --scramble-len 19 --seed 1919   # 19-move scrambles (optimal about 18)
    run rand3 -N 1 --seed 3                         # a random position (optimal 20)
    run scr17st -N 10 --scramble-len 17 --seed 1717 -j 1   # one thread: target passes
    sets=(scr19 rand3 scr17st)
fi
status=0
python3 - "$out" "${sets[*]}" "${same_nodes[*]}" <<'PY' || status=1
import json, sys, os
out, sets, same = sys.argv[1], sys.argv[2].split(), sys.argv[3].split()
bad = 0
for name in sets:
    early = json.load(open(os.path.join(out, f"{name}-early.json")))["results"]
    full = json.load(open(os.path.join(out, f"{name}-full.json")))["results"]
    for e, f in zip(early, full):
        ok = e["optimal"] and f["optimal"] and e["length"] == f["length"]
        if name in same:
            ok = ok and e["nodes"] == f["nodes"]
        bad += not ok
        print(f"{name} #{e['index']}: early {e['length']} ({e['elapsed']:.2f}s, {e['nodes']} nodes)  "
              f"full {f['length']} ({f['elapsed']:.2f}s, {f['nodes']} nodes)  {'ok' if ok else 'MISMATCH'}")
print("all proven lengths agree" if not bad else f"{bad} MISMATCHES")
sys.exit(1 if bad else 0)
PY
rm -rf "$out"
if [[ "$moves" == *slice* ]]; then
    cargo test -q --test solver with_slices -- --include-ignored || status=1
fi
exit $status
