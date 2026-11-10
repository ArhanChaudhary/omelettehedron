#!/usr/bin/env bash
set -u
B=/mnt/egg/super/fto-target-analysis/bin/fto-a1
D=/mnt/egg/super/fto-analysis-data
export FTO_TABLES=/mnt/egg/super/fto-tables
echo "load $(cut -d' ' -f1 /proc/loadavg)"
$B bench -N 10 --scramble-len 19 --seed 1919 --tier large --json $D/opt-base.json > $D/opt-base.txt
FTO_SKIPLAST=1 $B bench -N 10 --scramble-len 19 --seed 1919 --tier large --json $D/opt-skip.json > $D/opt-skip.txt
FTO_CBSYM=1 $B bench -N 10 -t 3 --seed 2026 --tier large -v --json $D/t3-cbsym.json > $D/t3-cbsym.txt 2>&1
$B bench -N 10 -t 3 --seed 2026 --tier large --json $D/t3-base.json > $D/t3-base.txt
FTO_CBSYM=1 $B bench -N 5 -t 3 -j 1 --seed 7 --tier large --json $D/t3st-cbsym.json > $D/t3st-cbsym.txt 2>&1
$B bench -N 5 -t 3 -j 1 --seed 7 --tier large --json $D/t3st-base.json > $D/t3st-base.txt
echo "load $(cut -d' ' -f1 /proc/loadavg)"
