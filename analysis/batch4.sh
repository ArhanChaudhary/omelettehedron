#!/usr/bin/env bash
set -u
D=/mnt/egg/super/fto-analysis-data
B=/mnt/egg/super/fto-target-analysis/bin/fto-a3
export FTO_TABLES=/mnt/egg/super/fto-tables
echo "load $(cut -d' ' -f1 /proc/loadavg)"
FTO_CBA=1 FTO_SKIPLAST=1 FTO_TARGET=1 FTO_DEBUG=1 $B bench -N 10 --scramble-len 19 --seed 1919 --tier large -v --json $D/b4-opt-cba.json > $D/b4-opt-cba.txt 2>&1
FTO_SKIPLAST=1 FTO_TARGET=1 FTO_DEBUG=1 $B bench -N 10 --scramble-len 19 --seed 1919 --tier large -v --json $D/b4-opt-nocba.json > $D/b4-opt-nocba.txt 2>&1
FTO_CBA=1 FTO_DEBUG=1 $B bench -N 2 -t 5 --seed 2026 --tier large -v --json $D/b4-t5-cba.json > $D/b4-t5-cba.txt 2>&1
FTO_DEBUG=1 $B bench -N 2 -t 5 --seed 2026 --tier large -v --json $D/b4-t5-nocba.json > $D/b4-t5-nocba.txt 2>&1
cd /mnt/egg/super/fto-wt-analysis && CARGO_TARGET_DIR=/mnt/egg/super/fto-target-analysis timeout 300 cargo test --release --test growth -- --ignored --nocapture > $D/growth.txt 2>&1
echo "load $(cut -d' ' -f1 /proc/loadavg)"
