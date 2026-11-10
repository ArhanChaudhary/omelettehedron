#!/usr/bin/env bash
B=/mnt/egg/super/fto-target-analysis/bin/fto-a4
export FTO_TABLES=/mnt/egg/super/fto-tables
for k in 0 256 1024 4096; do echo "k=$k $(FTO_DELAY=$k $B bench -N 2 -t 2 -j 1 --seed 7 --tier large | grep nodes/s)"; done
