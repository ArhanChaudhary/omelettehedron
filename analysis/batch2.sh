#!/usr/bin/env bash
set -u
D=/mnt/egg/super/fto-analysis-data
BIN=/mnt/egg/super/fto-target-analysis/bin
export FTO_TABLES=/mnt/egg/super/fto-tables
echo "load $(cut -d' ' -f1 /proc/loadavg)"
FTO_SKIPLAST=1 timeout 300 $BIN/fto-stats bench -N 10 --scramble-len 19 --seed 1919 --tier large -v > $D/stats-opt19.txt 2>&1
timeout 120 $BIN/fto-stats bench -N 2 -t 5 --seed 2026 --tier large -v > $D/stats-t5.txt 2>&1
FTO_DEBUG=1 FTO_SKIPLAST=1 $BIN/fto-a1 bench -N 10 --scramble-len 19 --seed 1919 --tier large -v > $D/debug-opt19.txt 2>&1
FTO_DEBUG=1 $BIN/fto-a1 bench -N 3 -t 5 --seed 2026 --tier large -v > $D/debug-t5.txt 2>&1
PERF_SECS=15 $D/perf-timed.sh $BIN/fto-a1 $D/perf-t5 -N 4 -t 5 --seed 2026 --tier large
echo "load $(cut -d' ' -f1 /proc/loadavg)"
