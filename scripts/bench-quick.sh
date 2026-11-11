#!/usr/bin/env bash
# Quick benchmark of both solvers, the default for A/B comparisons: one round takes about a
# minute with the large tiers on 48 threads (plus table loading), so a comparison of two builds
# in two interleaved rounds stays under 5 minutes. The long suites (scripts/bench.sh,
# bench-opt.sh, bench-optimal.sh, bench-randopt.sh, bench-hard.sh) stay for occasional runs.
#
# usage:
#   scripts/bench-quick.sh [options] <label>               build this tree and run one round
#   scripts/bench-quick.sh [options] --compare <A> <B>     two binaries, interleaved rounds
# options:
#   --rounds N          rounds (default 1; 2 with --compare)
#   --bin PATH          (single run) use this binary instead of building the tree
#   --tier T            two-phase tier (default large; "none" skips the two-phase solver)
#   --opt-tier T[,T]    optimal solver tier(s) (default large; "none" skips it)
# Results: bench-results/quick/<label>/<item>-r<round>.json (with --compare, the directory is
# bench-results/quick/compare-<time>/ and the files start with A- and B-), and a markdown table.
#
# The set (each item a fixed set of seeded positions):
#   two-phase, optimal mode: 10 x 17-move scrambles (seed 1717) on 1 thread, 10 x 19-move
#     scrambles (seed 1919) on all threads, with FTO_CBS_EAGER=1 (the distance-to-solved table
#     read at start-up: node counts on 1 thread repeat exactly, and the time leaves its reading
#     out);
#   two-phase, time limit 1 s: 10 random positions (seed 2026) on all threads, 10 (seed 7) on
#     1 thread (mean length after 1 s, first solution, nodes/s);
#   optimal solver: the same 17-move set on 1 thread and 19-move set on all threads.
# Every run is capped at 150 s. The runs hold the bench lock (/tmp/cube-bench.lock)
# exclusively; building takes it shared, so builds elsewhere wait for the timing and not
# the other way round.
set -euo pipefail
cd "$(dirname "$0")/.."

rounds=""
bin=""
tier=large
opt_tiers=large
compare=()
label=""
while [ $# -gt 0 ]; do
    case $1 in
        --rounds) rounds=$2; shift 2 ;;
        --bin) bin=$2; shift 2 ;;
        --tier) tier=$2; shift 2 ;;
        --opt-tier) opt_tiers=$2; shift 2 ;;
        --compare) compare=("$2" "$3"); shift 3 ;;
        -h|--help) sed -n '2,/^set -euo/p' "$0" | sed '$d'; exit 0 ;;
        -*) echo "unknown option $1" >&2; exit 2 ;;
        *) label=$1; shift ;;
    esac
done

exec 9>/tmp/cube-bench.lock
if [ ${#compare[@]} -eq 2 ]; then
    names=(A B)
    bins=("${compare[0]}" "${compare[1]}")
    rounds=${rounds:-2}
    out=bench-results/quick/compare-$(date -u +%Y%m%d-%H%M%S)
else
    [ -n "$label" ] || { echo "usage: $0 [options] <label> | --compare <A> <B>" >&2; exit 2; }
    if [ -z "$bin" ]; then
        flock -s 9
        cargo build --release -q
        bin=${CARGO_TARGET_DIR:-target}/release/omelettehedron
    fi
    names=("$label")
    bins=("$bin")
    rounds=${rounds:-1}
    out=bench-results/quick/$label
fi
for b in "${bins[@]}"; do
    [ -x "$b" ] || { echo "not an executable: $b" >&2; exit 2; }
done
mkdir -p "$out"
flock 9
{
    echo "commit $(git rev-parse --short HEAD 2>/dev/null || echo -), load before benchmark: $(cut -d' ' -f1 /proc/loadavg)"
    for i in "${!bins[@]}"; do echo "${names[$i]}: ${bins[$i]}"; done
} | tee "$out/info.txt" >&2

# item name -> command line (without the binary and --json)
items=()
declare -A cmd envs
if [ "$tier" != none ]; then
    items+=(2p-opt17-1t 2p-opt19-mt 2p-t1s-mt 2p-t1s-1t)
    cmd[2p-opt17-1t]="bench --tier $tier -N 10 --scramble-len 17 --seed 1717 -j 1"
    cmd[2p-opt19-mt]="bench --tier $tier -N 10 --scramble-len 19 --seed 1919"
    envs[2p-opt17-1t]=FTO_CBS_EAGER=1
    envs[2p-opt19-mt]=FTO_CBS_EAGER=1
    cmd[2p-t1s-mt]="bench --tier $tier -N 10 -t 1 --seed 2026"
    cmd[2p-t1s-1t]="bench --tier $tier -N 10 -t 1 --seed 7 -j 1"
fi
if [ "$opt_tiers" != none ]; then
    for t in ${opt_tiers//,/ }; do
        items+=("opt-$t-17-1t" "opt-$t-19-mt")
        cmd[opt-$t-17-1t]="optimal bench --tier $t -N 10 --scramble-len 17 --seed 1717 -j 1"
        cmd[opt-$t-19-mt]="optimal bench --tier $t -N 10 --scramble-len 19 --seed 1919"
    done
fi

n=${#bins[@]}
for r in $(seq 1 "$rounds"); do
    for item in "${items[@]}"; do
        for k in $(seq 0 $((n - 1))); do
            # rotate the order every round
            i=$(( (k + r - 1) % n ))
            f="$out/${names[$i]}-$item-r$r"
            [ $n -eq 1 ] && f="$out/$item-r$r"
            echo "round $r: ${names[$i]} $item" >&2
            # shellcheck disable=SC2086 # the command lines are split on purpose
            env ${envs[$item]:-} timeout 150 "${bins[$i]}" ${cmd[$item]} --json "$f.json" \
                > "$f.txt" 2>&1 || echo "  failed or over 150 s (see $f.txt)" >&2
        done
    done
done
flock -u 9

python3 - "$out" "$rounds" "${names[@]}" -- "${items[@]}" <<'PY'
import json, os, statistics, sys

out, rounds = sys.argv[1], int(sys.argv[2])
sep = sys.argv.index("--")
names, items = sys.argv[3:sep], sys.argv[sep + 1:]
single = len(names) == 1
med = statistics.median


def load(name, item, r):
    p = f"{out}/{item}-r{r}.json" if single else f"{out}/{name}-{item}-r{r}.json"
    try:
        return json.load(open(p))
    except (OSError, ValueError):
        return None


def metrics(item, d):
    """(time, nodes, text) of one run: time is the total time of an optimal set or 1 / (nodes
    per second) of a time-limited set, so that lower is better for both."""
    if item.startswith("opt-"):
        return d["total_time"], d["total_nodes"], (
            f"{d['total_time']:.3f} s, {d['total_nodes']:.4g} nodes, {d['solved']}/{len(d['results'])} "
            f"solved, mean {d['mean_length']:.2f}")
    if "-t1s-" in item:
        return 1 / d["nodes_per_sec"], None, (
            f"mean {d['mean_length']:.3f} after 1 s, first {d['mean_first_length']:.2f} @ "
            f"{d['mean_first_time'] * 1e3:.2f} ms, {d['nodes_per_sec']:.3e} nodes/s")
    t = sum(x["elapsed"] for x in d["results"])
    nodes = sum(x["nodes"] for x in d["results"])
    return t, nodes, (f"{t:.3f} s, {nodes:.4g} nodes, {d['solved_optimally']}/{len(d['results'])} "
                      f"proven, mean {d['mean_length']:.2f}")


if single:
    print(f"| item | {' | '.join(f'round {r}' for r in range(1, rounds + 1))} |")
    print("|---|" + "---|" * rounds)
    for item in items:
        cells = []
        for r in range(1, rounds + 1):
            d = load(names[0], item, r)
            cells.append(metrics(item, d)[2] if d else "-")
        print(f"| {item} | {' | '.join(cells)} |")
    sys.exit(0)

a, b = names
print(f"| item | {a} (round 1) | {b} (round 1) | {b} / {a} (medians) | per round | nodes |")
print("|---|---|---|---|---|---|")
for item in items:
    runs = [(load(a, item, r), load(b, item, r)) for r in range(1, rounds + 1)]
    runs = [(x, y) for x, y in runs if x and y]
    if not runs:
        print(f"| {item} | - | - | - | - | - |")
        continue
    ma = [metrics(item, x) for x, _ in runs]
    mb = [metrics(item, y) for _, y in runs]
    ta, tb = med(m[0] for m in ma), med(m[0] for m in mb)
    kind = "nodes/s" if "-t1s-" in item else "time"
    ratio = ta / tb if kind == "nodes/s" else tb / ta
    per = ", ".join(f"{(x[0] / y[0]) if kind == 'nodes/s' else (y[0] / x[0]):.3f}" for x, y in zip(ma, mb))
    if ma[0][1] is None:
        nodes = "-"
    elif "-1t" in item:
        nodes = "same" if all(x[1] == y[1] for x, y in zip(ma, mb)) else f"{mb[0][1] / ma[0][1]:.4f}x"
    else:
        nodes = f"{med(m[1] for m in mb) / med(m[1] for m in ma):.4f}x"
    print(f"| {item} | {ma[0][2]} | {mb[0][2]} | {kind} {ratio:.3f}x | {per} | {nodes} |")
print(f"\n(the ratio is of the medians over {rounds} rounds: time for optimal sets, nodes/s for "
      "time-limited ones; nodes on 1 thread must be the same for a change that keeps the search; "
      "single items vary by a few percent between rounds)")
PY
