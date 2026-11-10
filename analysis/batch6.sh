#!/usr/bin/env bash
set -u
D=/mnt/egg/super/fto-analysis-data
B=/mnt/egg/super/fto-target-analysis/bin/fto-a3
export FTO_TABLES=/mnt/egg/super/fto-tables
echo "load $(cut -d' ' -f1 /proc/loadavg)"
FTO_DEBUG=1 FTO_SKIPLAST=1 timeout 300 $B solve --random 1 --seed 1 --tier large -v -p improvements > $D/b6-rand1-skip.txt 2>&1
FTO_DEBUG=1 FTO_SKIPLAST=1 FTO_TARGET=1 timeout 300 $B solve --random 1 --seed 1 --tier large -v -p improvements > $D/b6-rand1-target.txt 2>&1
echo "load $(cut -d' ' -f1 /proc/loadavg)"
FTO_SKIPLAST=1 PERF_SECS=10 $D/perf-timed.sh $B $D/perf-opt -N 10 --scramble-len 19 --seed 1919 --tier large
