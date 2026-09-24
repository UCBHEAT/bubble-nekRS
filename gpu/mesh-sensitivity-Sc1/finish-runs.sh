#!/usr/bin/env bash
# After the last job of each run, e.g. ./finish-runs.sh 1x 1.5x 2.25x
# moves the final part's field files to part<N>/ like prepare-restart.sh,
# then links every checkpoint of every part, in time order, as
# bubble3d0.f00000, bubble3d0.f00001, ... and writes bubble3d.nek5000 so the
# whole run opens as one time series in ParaView or VisIt.
set -euo pipefail

for dir in "$@"; do
    (
    cd "$dir"
    if compgen -G "bubble3d0.f[0-9]*" > /dev/null && [ ! -L "$(ls bubble3d0.f[0-9]* | head -n 1)" ]; then
        part=1
        while [ -e part$part ]; do part=$((part + 1)); done
        mkdir part$part
        mv bubble3d0.f[0-9]* part$part/
        mv logfile-* nodes-* part$part/ 2>/dev/null || true
        if [ -e bubble3d.nek5000 ]; then mv bubble3d.nek5000 part$part/; fi
    fi
    find . -maxdepth 1 -type l -name "bubble3d0.f[0-9]*" -delete

    # Order all checkpoints by their simulation time (header field 8).
    i=0
    for f in $(for g in part*/bubble3d0.f[0-9]*; do
                   echo "$(head -c 132 "$g" | awk '{print $8}') $g"
               done | sort -g | awk '{print $2}'); do
        ln -s "$f" "$(printf "bubble3d0.f%05d" $i)"
        i=$((i + 1))
    done
    cat > bubble3d.nek5000 <<EOF
 filetemplate: bubble3d%01d.f%05d
 firsttimestep: 0
 numtimesteps: $i
EOF
    echo "$dir: linked $i checkpoints"
    )
done
