#!/usr/bin/env bash
# Frontier counterpart of submit-animate.sh: postprocess and animate finished
# runs, each on its own CPU node, e.g. from the directory holding them:
#   PROJ_ID=fus167 FFMPEG=<ffmpeg> PV_SITE=<dir> ./submit-animate-frontier.sh 01:00 \
#       "1x:2:Sc = 4, 1x Kolmogorov (23805 elements)" \
#       "0.296x:3:Sc = 4, 0.296x Kolmogorov (441864 elements)" ...
# As there, each <run dir>:<chunks>:<label> splits the run's checkpoints
# between <chunks> pvbatch processes for the bubble statistics (post.csv) and
# again for the frames, which ffmpeg encodes to <run dir>/animation.mp4, and
# finished chunks and frames are skipped, so a job that runs out of time can
# simply be resubmitted. Frontier's ParaView renders with OSMesa on the CPU,
# so a run's processes share its node's cores for VTK's threads (the
# STDThread backend: the default, Sequential, is about 10x slower) and its
# memory. They run without MPI, one task each of one srun step: animate.py
# reads whole checkpoints, which parallel pvbatch would split between ranks.
#   PARTITION, QOS  as for ../common/submit-bundle-frontier.sh
#   PARAVIEW_BIN    ParaView's bin directory (default: ParaView 5.13.1 with
#                   OSMesa from the Frontier software stack)
#   FFMPEG          ffmpeg executable with libx264 (default: ffmpeg from PATH;
#                   Frontier has none, so set it to a static build)
#   PV_SITE         directory with numpy for ParaView's Python 3.11, which has
#                   none (pip install --target <dir> 'numpy<2' with
#                   cray-python/3.11.7); required
set -euo pipefail

: ${PROJ_ID:?PROJ_ID must be set}
: ${PARTITION:=batch}
: ${QOS:=}
: ${PARAVIEW_BIN:=/sw/frontier/spack-envs/cpe24.11-cpu/opt/gcc-13.2/paraview-5.13.1-5tu7m74h2ny7rd24x5utrgqkseltz32r/bin}
# Frontier has no ffmpeg: set FFMPEG to a static build with libx264.
: ${FFMPEG:=ffmpeg}
: ${PV_SITE:?PV_SITE must be set: a directory with numpy for the ParaView Python}

if [ $# -lt 2 ]; then
    echo "Usage: PROJ_ID=<project> $0 <hh:mm> <run dir>:<chunks>:<label> ..."
    exit 1
fi
time=$1
shift
script=$(cd "$(dirname "$0")" && pwd -P)/animate.py
[ -x "$PARAVIEW_BIN/pvbatch" ] || { echo "Cannot find $PARAVIEW_BIN/pvbatch"; exit 1; }
command -v "$FFMPEG" > /dev/null || { echo "Cannot find $FFMPEG"; exit 1; }
[ -d "$PV_SITE/numpy" ] || { echo "Cannot find numpy in $PV_SITE"; exit 1; }

SFILE=animate.sbatch
cat > $SFILE <<EOF
#!/bin/bash
#SBATCH -A $PROJ_ID
#SBATCH -J animate
#SBATCH -o animate-%j.out
#SBATCH -t ${time}:00
#SBATCH -N $#
#SBATCH -p $PARTITION
#SBATCH --exclusive
EOF
if [ -n "$QOS" ]; then
    echo "#SBATCH -q $QOS" >> $SFILE
fi
cat >> $SFILE <<EOF

module load PrgEnv-gnu gcc-native/13.2 cray-mpich/8.1.31
export LD_LIBRARY_PATH=$(dirname $PARAVIEW_BIN)/lib64:\$LD_LIBRARY_PATH
export PYTHONPATH=$PV_SITE
export VTK_SMP_BACKEND_IN_USE=STDThread
pvbin=$PARAVIEW_BIN
script=$script
ffmpeg=$FFMPEG
EOF
for spec in "$@"; do
    IFS=: read -r d chunks label <<< "$spec"
    d=$(cd "$d" && pwd -P)
    [ -f "$d/data.csv" ] || { echo "Cannot find $d/data.csv"; exit 1; }
    echo "specs+=(\"$d:$chunks:$label\")" >> $SFILE
done
cat >> $SFILE <<'EOF'
mapfile -t hosts < <(scontrol show hostnames "$SLURM_JOB_NODELIST")

# run_chunks <mode> <run dir> <chunks> <host> [label]: the chunks as the tasks
# of one srun step on the run's node, each with an equal share of its CPUs.
run_chunks() {
    local mode=$1 dir=$2 chunks=$3 host=$4 label=${5:-}
    local case=$(basename $dir/*.par .par)
    local n=$(ls $dir/${case}0.f[0-9]* | wc -l)
    local per=$(( (n + chunks - 1)/chunks ))
    local ntasks=$(( (n + per - 1)/per ))
    local c=$((SLURM_CPUS_ON_NODE/ntasks))
    mkdir -p $dir/animate-logs
    # pvbatch cannot take an empty argument, so the label goes to frames only.
    VTK_SMP_MAX_THREADS=$c LP_NUM_THREADS=$c srun -N 1 -n $ntasks -c $c -w $host \
        --output=$dir/animate-logs/$mode-%t.log bash -c \
        'first=$((SLURM_PROCID*$3)); exec $0 --no-mpi $1 '$mode' $2 $first $((first + $3)) "${@:4}"' \
        $pvbin/pvbatch $script $dir $per ${label:+"$label"}
}

# animate <run dir> <chunks> <label> <host>
animate() {
    local d=$1 chunks=$2 label=$3 host=$4
    echo "$(date) $d: bubble statistics on $host"
    run_chunks post $d $chunks $host
    $pvbin/pvpython $script merge $d
    echo "$(date) $d: frames"
    run_chunks frames $d $chunks $host "$label"
    local case=$(basename $d/*.par .par)
    local n=$(ls $d/${case}0.f[0-9]* | wc -l)
    if [ $(ls $d/frames/frame_*.png | grep -v tmp | wc -l) -eq $n ]; then
        $ffmpeg -hide_banner -loglevel error -y -framerate 10 -i $d/frames/frame_%05d.png \
            -c:v libx264 -pix_fmt yuv420p -crf 18 $d/animation.mp4
        echo "$(date) $d/animation.mp4"
    else
        echo "$(date) $d: frames incomplete, resubmit to finish"
    fi
}

for k in "${!specs[@]}"; do
    IFS=: read -r d chunks label <<< "${specs[k]}"
    animate $d $chunks "$label" ${hosts[k]} &
done
wait
echo "$(date) done"
EOF

sbatch --parsable $SFILE
