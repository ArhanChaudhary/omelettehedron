#!/usr/bin/env bash
# usage: perf-timed.sh <bin> <out-prefix> <bench args...>   (run under flock)
set -u
bin=$1; out=$2; shift 2
export FTO_TABLES=/mnt/egg/super/fto-tables
$bin bench "$@" > $out.txt 2>&1 &
pid=$!
sleep 4
chrt -i 0 perf record -q -e cpu-clock -F 499 -p $pid -o $out.perf -- sleep ${PERF_SECS:-15} 2>/dev/null
wait $pid
