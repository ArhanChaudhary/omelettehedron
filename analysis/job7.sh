#!/usr/bin/env bash
set -u
D=/mnt/egg/super/fto-analysis-data
cd /mnt/egg/super/fto-wt-analysis
export CARGO_TARGET_DIR=/mnt/egg/super/fto-target-analysis
cargo build --release 2>&1 | grep -E "^(warning|error)|Finished|-->" -A7
cp /mnt/egg/super/fto-target-analysis/release/fto /mnt/egg/super/fto-target-analysis/bin/fto-a4 || exit 1
B=/mnt/egg/super/fto-target-analysis/bin/fto-a4
export FTO_TABLES=/mnt/egg/super/fto-tables
echo "load $(cut -d' ' -f1 /proc/loadavg)"
for k in 0 64 256 0; do
  FTO_DELAY=$k $B bench -N 6 -t 3 --seed 2026 --tier large > $D/b7-mt-d$k.txt 2>&1
  echo "mt delay $k: $(grep nodes/s $D/b7-mt-d$k.txt)"
done
for k in 0 64 256; do
  FTO_DELAY=$k $B bench -N 3 -t 3 -j 1 --seed 7 --tier large > $D/b7-st-d$k.txt 2>&1
  echo "st delay $k: $(grep nodes/s $D/b7-st-d$k.txt)"
done
echo "load $(cut -d' ' -f1 /proc/loadavg)"
