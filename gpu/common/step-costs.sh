#!/usr/bin/env bash
# Mean wall time per time step of a nekRS-LS run, from its log, split into
# ordinary steps and steps that also ran a TLSR or CLSR reinitialization
# (pseudo-step lines before the step's elapsedStep line), from a given step on
# (default 200, which skips the start-up transient), e.g.
#   <gpu/common>/step-costs.sh 0.444x/logfile-5580267 500
# The overall mean sets the run time; the split shows how much of it the
# level set reinitialization takes, which is latency-bound on small meshes.
set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 <nekRS log> [first step]"
    exit 1
fi

awk -v first="${2:-200}" '
/^Pseudo-step=.*TLSR/ {tlsr = 1}
/^Pseudo-step=.*CLSR/ {clsr = 1}
/^step= *[0-9]+ +elapsedStep=/ {
    s = $0; gsub("step=", "", s); split(s, f, " ")
    e = f[3]; sub("s$", "", e)
    if (f[1] >= first) {
        k = tlsr ? "TLSR" : (clsr ? "CLSR" : "ordinary")
        sum[k] += e; n[k]++; total += e; all++
    }
    tlsr = 0; clsr = 0
}
END {
    if (!all) { print "no steps from step " first " on"; exit 1 }
    for (k in n) printf "  %-8s %8d steps  mean %.4f s\n", k, n[k], sum[k]/n[k]
    printf "  overall  %8d steps  mean %.4f s/step\n", all, total/all
}' "$1"
