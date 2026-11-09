#!/usr/bin/env python3
"""Compare two benchmark labels: scripts/compare.py <before> <after> -> markdown table."""
import json, sys
def load(label, kind):
    return json.load(open(f"bench-results/{label}-{kind}.json"))
b, a = sys.argv[1], sys.argv[2]
rows = []
for kind, name in (("mt", "48 threads"), ("st", "1 thread")):
    x, y = load(b, kind), load(a, kind)
    for (t, lx), (_, ly) in zip(x["checkpoints"], y["checkpoints"]):
        rows.append((f"{name}: mean length after {t:g}s", f"{lx:.3f}", f"{ly:.3f}", f"{ly-lx:+.3f}"))
    rows.append((f"{name}: mean first-solution time (s)", f"{x['mean_first_time']:.4f}", f"{y['mean_first_time']:.4f}", f"{y['mean_first_time']/x['mean_first_time']:.2f}x"))
    rows.append((f"{name}: mean first-solution length", f"{x['mean_first_length']:.3f}", f"{y['mean_first_length']:.3f}", f"{y['mean_first_length']-x['mean_first_length']:+.3f}"))
    rows.append((f"{name}: nodes/s", f"{x['nodes_per_sec']:.3e}", f"{y['nodes_per_sec']:.3e}", f"{y['nodes_per_sec']/x['nodes_per_sec']:.2f}x"))
print(f"| metric | before ({b}) | after ({a}) | change |")
print("|---|---|---|---|")
for r in rows:
    print("| " + " | ".join(r) + " |")
