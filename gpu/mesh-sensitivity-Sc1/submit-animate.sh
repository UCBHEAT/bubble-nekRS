#!/usr/bin/env bash
# Postprocess and animate finished runs on one Polaris GPU node, e.g.
#   PROJ_ID=nek-vf QUEUE=debug ./submit-animate.sh 01:00 \
#       "1x:20:1x Kolmogorov (27x54x27)" "1.5x:8:1.5x Kolmogorov (18x36x18)" ...
# Each <run dir>:<chunks>:<label> splits the run's checkpoints into <chunks>
# pvbatch processes (one CPU core each) for the bubble statistics
# (post.csv) and again for the frames, which ffmpeg encodes to
# <run dir>/animation.mp4. Finished chunks and frames are skipped, so a job
# that runs out of time can simply be resubmitted.
set -euo pipefail

: ${PROJ_ID:?PROJ_ID must be set}
: ${QUEUE:?QUEUE must be set}
: ${FFMPEG:=/lus/eagle/projects/nek-vf/benl/answinter26/tools/ffmpeg}
: ${PARAVIEW_MODULE:=visualization/paraview/paraview-5.13.1-EGL}

if [ $# -lt 2 ]; then
    echo "Usage: PROJ_ID=<project> QUEUE=<queue> $0 <hh:mm> <run dir>:<chunks>:<label> ..."
    exit 1
fi
time=$1
shift
script=$(cd "$(dirname "$0")" && pwd)/animate.py

SFILE=animate.batch
cat > $SFILE <<EOF
#!/bin/bash
#PBS -A $PROJ_ID
#PBS -N animate
#PBS -q $QUEUE
#PBS -l walltime=${time}:00
#PBS -l filesystems=home:eagle:grand
#PBS -l select=1:system=polaris
#PBS -k doe
#PBS -j oe

cd \$PBS_O_WORKDIR
module use /soft/modulefiles
module load $PARAVIEW_MODULE
script=$script
ffmpeg=$FFMPEG
EOF
cat >> $SFILE <<'EOF'
ngpu=4

# run_chunks <mode> <run dir> <chunks> [label]
run_chunks() {
    local mode=$1 dir=$2 chunks=$3 label=${4:-}
    local n=$(ls $dir/bubble3d0.f[0-9]* | wc -l)
    local per=$(( (n + chunks - 1)/chunks ))
    mkdir -p $dir/animate-logs
    for ((k = 0; k < chunks; k++)); do
        local first=$((k*per)) last=$(((k+1)*per))
        [ $first -lt $n ] || break
        local dev=$((gpu % ngpu))
        gpu=$((gpu + 1))
        VTK_DEFAULT_EGL_DEVICE_INDEX=$dev pvbatch $script $mode $dir $first $last "$label" \
            > $dir/animate-logs/$mode-$first.log 2>&1 &
    done
}

EOF
dirs=()
for spec in "$@"; do
    IFS=: read -r d chunks label <<< "$spec"
    d=$(cd "$d" && pwd)
    dirs+=("$d")
    echo "specs+=(\"$d:$chunks:$label\")" >> $SFILE
done
cat >> $SFILE <<'EOF'

echo "$(date) bubble statistics"
gpu=0
for spec in "${specs[@]}"; do
    IFS=: read -r d chunks label <<< "$spec"
    run_chunks post $d $chunks
done
wait
for spec in "${specs[@]}"; do
    IFS=: read -r d chunks label <<< "$spec"
    pvpython $script merge $d
done

echo "$(date) frames"
gpu=0
for spec in "${specs[@]}"; do
    IFS=: read -r d chunks label <<< "$spec"
    run_chunks frames $d $chunks "$label"
done
wait

echo "$(date) encoding"
for spec in "${specs[@]}"; do
    IFS=: read -r d chunks label <<< "$spec"
    n=$(ls $d/bubble3d0.f[0-9]* | wc -l)
    if [ $(ls $d/frames/frame_*.png | grep -v tmp | wc -l) -eq $n ]; then
        $ffmpeg -hide_banner -loglevel error -y -framerate 10 -i $d/frames/frame_%05d.png \
            -c:v libx264 -pix_fmt yuv420p -crf 18 $d/animation.mp4
        echo "$d/animation.mp4"
    else
        echo "$d: frames incomplete, resubmit to finish"
    fi
done
echo "$(date) done"
EOF

qsub -q $QUEUE $SFILE
