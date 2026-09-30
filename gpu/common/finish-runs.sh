#!/usr/bin/env bash
# After the last job of each run, e.g. ./finish-runs.sh 1x 1.5x 2.25x
# moves the final part's field files to part<N>/ like prepare-restart.sh,
# then links every checkpoint of every part, in time order, as
# <case>0.f00000, <case>0.f00001, ... and writes <case>.nek5000 so the
# whole run opens as one time series in ParaView or VisIt (<case>: the name
# of the run directory's .par file, e.g. bubble3d).
set -euo pipefail

for dir in "$@"; do
    (
    cd "$dir"
    pars=(*.par)
    if [ ${#pars[@]} -ne 1 ] || [ ! -f "${pars[0]}" ]; then
        echo "$dir: must hold exactly one .par file"
        exit 1
    fi
    casename=${pars[0]%.par}
    if compgen -G "${casename}0.f[0-9]*" > /dev/null && [ ! -L "$(ls ${casename}0.f[0-9]* | head -n 1)" ]; then
        part=1
        while [ -e part$part ]; do part=$((part + 1)); done
        mkdir part$part
        mv ${casename}0.f[0-9]* part$part/
        mv logfile-* nodes-* nekRS_*.e* part$part/ 2>/dev/null || true
        if [ -e $casename.nek5000 ]; then mv $casename.nek5000 part$part/; fi
    fi
    find . -maxdepth 1 -type l -name "${casename}0.f[0-9]*" -delete

    # Order all checkpoints by their simulation time (header field 8).
    i=0
    for f in $(for g in part*/${casename}0.f[0-9]*; do
                   echo "$(head -c 132 "$g" | awk '{print $8}') $g"
               done | sort -g | awk '{print $2}'); do
        ln -s "$f" "$(printf "${casename}0.f%05d" $i)"
        i=$((i + 1))
    done
    cat > $casename.nek5000 <<EOF
 filetemplate: $casename%01d.f%05d
 firsttimestep: 0
 numtimesteps: $i
EOF
    echo "$dir: linked $i checkpoints"
    )
done
