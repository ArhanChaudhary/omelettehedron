#!/usr/bin/env python3
"""Compare two benchmark labels: scripts/compare.py <before> <after> -> markdown table."""
import json, sys
def f3(v):
    return "-" if v is None else f"{v:.3f}"
def d3(a, b):
    return "-" if a is None or b is None else f"{b-a:+.3f}"
def load(label, kind):
    return json.load(open(f"bench-results/{label}-{kind}.json"))
b, a = sys.argv[1], sys.argv[2]
rows = []
for kind, name in (("mt", "48 threads"), ("st", "1 thread")):
    x, y = load(b, kind), load(a, kind)
    n = x["scrambles"]
    for (t, lx, kx), (_, ly, ky) in zip(x["checkpoints"], y["checkpoints"]):
        rows.append((f"{name}: mean length after {t:g}s (solved/{n})", f"{f3(lx)} ({kx})", f"{f3(ly)} ({ky})", d3(lx, ly)))
    rows.append((f"{name}: mean first-solution time (s)", f"{x['mean_first_time']:.4f}", f"{y['mean_first_time']:.4f}", f"{y['mean_first_time']/x['mean_first_time']:.2f}x"))
    rows.append((f"{name}: mean first-solution length", f"{x['mean_first_length']:.3f}", f"{y['mean_first_length']:.3f}", f"{y['mean_first_length']-x['mean_first_length']:+.3f}"))
    rows.append((f"{name}: nodes/s", f"{x['nodes_per_sec']:.3e}", f"{y['nodes_per_sec']:.3e}", f"{y['nodes_per_sec']/x['nodes_per_sec']:.2f}x"))
print(f"| metric | before ({b}) | after ({a}) | change |")
print("|---|---|---|---|")
for r in rows:
    print("| " + " | ".join(r) + " |")
