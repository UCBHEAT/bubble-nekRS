#!/usr/bin/env bash
# Progress of the cases of submit-bundle-frontier.sh jobs (Slurm), e.g. from
# the directory holding the runs:
#   <gpu/common>/bundle-status.sh . 5580267 5580268
# prints each job's state and, for every run directory with a log from one of
# them (<dir>/logfile-<jobid>), its last step, time and CFL, the seconds per
# step over its last 1000 steps, its last Sh and any error lines.
#
# With --wait <seconds> (first argument) it first blocks until one of the jobs
# changes state, a case log shows a new error or exit status line, or the
# given time passes, and says which, e.g. to be told when a queued job starts:
#   <gpu/common>/bundle-status.sh --wait 7200 . 5580267
set -euo pipefail

wait_s=0
if [ "${1:-}" = --wait ]; then
    wait_s=$2
    shift 2
fi
if [ $# -lt 2 ]; then
    echo "Usage: $0 [--wait <seconds>] <study dir> <jobid> [<jobid> ...]"
    exit 1
fi
study=$1
shift
jobs=("$@")
errpat='ERROR|Abort|exit status|Segmentation fault|hipError|out of memory|srun: error|Unreasonable|[^a-z]nan[^a-z]|NaN'

jobstate() {
    local s
    s=$(squeue -h -j "$1" -o %T 2>/dev/null || true)
    [ -n "$s" ] || s=$(sacct -X -n -j "$1" -o State%20 2>/dev/null | head -1 | awk '{print $1}')
    echo "${s:-UNKNOWN}"
}
nerrors() {
    local n=0 j log
    for j in "${jobs[@]}"; do
        for log in "$study"/*/logfile-"$j"; do
            [ -f "$log" ] && n=$((n + $(grep -c -E "$errpat" "$log" || true)))
        done
    done
    echo $n
}

if [ "$wait_s" -gt 0 ]; then
    declare -A state0
    for j in "${jobs[@]}"; do state0[$j]=$(jobstate "$j"); done
    err0=$(nerrors)
    start=$(date +%s)
    event=""
    while [ -z "$event" ]; do
        sleep 60
        for j in "${jobs[@]}"; do
            s=$(jobstate "$j")
            [ "$s" = "${state0[$j]}" ] || event+="job $j: ${state0[$j]} -> $s. "
        done
        [ "$(nerrors)" -le "$err0" ] || event+="new error or exit status lines in the logs. "
        if [ -z "$event" ] && [ $(($(date +%s) - start)) -ge "$wait_s" ]; then
            event="nothing after $((wait_s / 60)) min. "
        fi
    done
    echo "$(date '+%m-%d %H:%M') $event"
fi

for j in "${jobs[@]}"; do
    echo "job $j: $(jobstate "$j")"
done
for j in "${jobs[@]}"; do
    for log in "$study"/*/logfile-"$j"; do
        [ -f "$log" ] || continue
        d=$(basename "$(dirname "$log")")
        st=$(grep -E '^step= *[0-9]+ +t= ' "$log" | tail -1 |
             awk '{gsub("step=", ""); printf "step=%s t=%s CFL=%s", $1, $3, $6}' || true)
        rate=$(grep -E '^step= *[0-9]+ +elapsedStep=' "$log" | tail -1001 |
               awk '{gsub("step=", ""); n[NR] = $1; s[NR] = $5 + 0}
                    END {if (NR > 1) printf "%.4f s/step", (s[NR] - s[1])/(n[NR] - n[1])}' || true)
        sh=$(grep -E '^mdot=' "$log" | tail -1 | awk '{print $3}' || true)
        echo "  $j $d: $st $rate $sh"
        grep -E "$errpat" "$log" | tail -2 | cut -c1-250 | sed 's/^/      /' || true
    done
done
