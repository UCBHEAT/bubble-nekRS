#!/usr/bin/env bash
# Prepare run directories to continue from their last checkpoint after a job
# ended before endTime, e.g. ./prepare-restart.sh 1x 1.5x
#
# nekRS restarts checkpoint numbering at 0, so the finished part's field
# files and logs move to part<N>/ first. data.csv keeps accumulating in
# place; rows written after the restart checkpoint are dropped so the
# restarted run doesn't duplicate them.
set -euo pipefail

for dir in "$@"; do
    (
    cd "$dir"
    files=($(ls bubble3d0.f[0-9]* 2>/dev/null || true))
    if [ ${#files[@]} -eq 0 ]; then
        echo "$dir: no checkpoints to restart from"
        exit 1
    fi

    # A job killed while writing leaves a short last file; don't restart from it.
    last=${files[-1]}
    if [ ${#files[@]} -gt 2 ] && \
       [ $(stat -c %s "$last") -ne $(stat -c %s "${files[-2]}") ]; then
        echo "$dir: $last is incomplete, moving it to $last.incomplete"
        mv "$last" "$last.incomplete"
        last=${files[-2]}
    fi
    t_restart=$(head -c 132 "$last" | awk '{print $8}')

    part=1
    while [ -e part$part ]; do part=$((part + 1)); done
    mkdir part$part
    mv bubble3d0.f[0-9]* part$part/
    mv logfile-* nodes-* nekRS_*.e* part$part/ 2>/dev/null || true
    if [ -e bubble3d.nek5000 ]; then mv bubble3d.nek5000 part$part/; fi

    ln -sfn part$part/$last restart.fld
    sed -i 's/^#\?startFrom *=.*/startFrom = restart.fld/' bubble3d.par

    if [ -e data.csv ]; then
        # data.csv times are rounded to 4 decimals.
        awk -F, -v t="$t_restart" 'NR == 1 || $1 <= t + 6e-5' data.csv > data.csv.tmp
        mv data.csv.tmp data.csv
    fi
    echo "$dir: restarting from part$part/$last (t = $t_restart)"
    )
done
