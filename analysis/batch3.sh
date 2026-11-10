#!/usr/bin/env bash
set -u
D=/mnt/egg/super/fto-analysis-data
B=/mnt/egg/super/fto-target-analysis/bin/fto-a3
export FTO_TABLES=/mnt/egg/super/fto-tables
echo "load $(cut -d' ' -f1 /proc/loadavg)"
$B bench -N 10 --scramble-len 19 --seed 1919 --tier large --json $D/b3-opt-base.json > $D/b3-opt-base.txt
FTO_SKIPLAST=1 $B bench -N 10 --scramble-len 19 --seed 1919 --tier large --json $D/b3-opt-skip.json > $D/b3-opt-skip.txt
FTO_SKIPLAST=1 FTO_TARGET=1 $B bench -N 10 --scramble-len 19 --seed 1919 --tier large --json $D/b3-opt-target.json > $D/b3-opt-target.txt
$B bench -N 10 --scramble-len 17 --seed 1717 -j 1 --tier large --json $D/b3-optst-base.json > $D/b3-optst-base.txt
FTO_SKIPLAST=1 $B bench -N 10 --scramble-len 17 --seed 1717 -j 1 --tier large --json $D/b3-optst-skip.json > $D/b3-optst-skip.txt
FTO_SKIPLAST=1 FTO_TARGET=1 $B bench -N 10 --scramble-len 17 --seed 1717 -j 1 --tier large --json $D/b3-optst-target.json > $D/b3-optst-target.txt
PERF_SECS=12 $D/perf-timed.sh $B $D/perf-t5-1t -N 4 -t 5 -j 1 --seed 7 --tier large
echo "load $(cut -d' ' -f1 /proc/loadavg)"
