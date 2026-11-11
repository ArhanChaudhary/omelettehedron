#!/usr/bin/env python3
"""Compare two optimal-solver benchmark labels (scripts/bench-optimal.sh):
scripts/compare-optimal.py <before> <after> -> markdown table."""
import json, os, sys
b, a = sys.argv[1], sys.argv[2]
rows = []
for kind, name in (("optimal1t", "1 thread, 17-move scrambles"), ("optimalmt", "all threads, 19-move scrambles")):
    pb, pa = f"bench-results/{b}-{kind}.json", f"bench-results/{a}-{kind}.json"
    if not (os.path.exists(pb) and os.path.exists(pa)):
        continue
    x, y = json.load(open(pb)), json.load(open(pa))
    same = [r["length"] for r in x["results"]] == [r["length"] for r in y["results"]]
    rows.append((f"{name}: solved / mean length", f"{x['solved']} / {x['mean_length']:.2f}", f"{y['solved']} / {y['mean_length']:.2f}", "same lengths" if same else "LENGTHS DIFFER"))
    rows.append((f"{name}: total time (s)", f"{x['total_time']:.2f}", f"{y['total_time']:.2f}", f"{y['total_time']/x['total_time']:.2f}x"))
    rows.append((f"{name}: nodes", f"{x['total_nodes']}", f"{y['total_nodes']}", f"{y['total_nodes']/x['total_nodes']:.3f}x"))
    rows.append((f"{name}: lookups", f"{x['total_lookups']}", f"{y['total_lookups']}", f"{y['total_lookups']/max(x['total_lookups'], 1):.3f}x"))
    rows.append((f"{name}: nodes/s", f"{x['nodes_per_sec']:.3e}", f"{y['nodes_per_sec']:.3e}", f"{y['nodes_per_sec']/x['nodes_per_sec']:.2f}x"))
print(f"| metric (optimal solver) | before ({b}) | after ({a}) | change |")
print("|---|---|---|---|")
for r in rows:
    print("| " + " | ".join(r) + " |")
