#!/usr/bin/env bash
# Frontier (Slurm) counterpart of submit-bundle.sh: submit several run
# directories of a mesh study as one job, each case on its own set of nodes,
# e.g. from the directory holding the runs:
#   PROJ_ID=fus167 <gpu/common>/submit-bundle-frontier.sh 02:00 0.444x:8 0.667x:3 1x:1 1.5x:1 2.25x:1
# The case name is that of the run directory's .par file (bubble.par: bubble).
# Batch jobs of 1-91 nodes may run for at most 2 hours (92-183 nodes: 6 h);
# PARTITION=extended allows up to 24 h on at most 64 nodes, with one running
# job per user.
#
# The job environment follows nekRS_HPCsupport/Frontier/nrsqsub with
# PrgEnv-gnu and gcc-native/13.2, which nekRS must have been built with.
# The batch script is written to bundle.sbatch, the Slurm log to
# nekRS_bundle-<jobid>.out, and each case's nekRS output to
# <dir>/logfile-<jobid>. As in submit-bundle.sh, each case stops with a
# checkpoint STOP_MARGIN seconds (default 300) before the walltime runs out
# (BUBBLE_STOP_AT, see the case's .udf), so prepare-restart.sh continues from
# where it stopped, and nodes that fail a HIP context check (gpucheck-hip.cpp)
# at the start of the job are left out, the case with the most nodes running
# on that many fewer.
#
#   PARTITION=<p>      batch (default) or extended
#   QOS=<q>            e.g. debug (at most 2 h, one queued job per user)
#   PREPARE_RESTART=1  first run prepare-restart.sh on each case that has new
#                      checkpoints from an earlier job, and skip any case whose
#                      data.csv has reached endTime
#   DEPEND=<jobid>     start only after job <jobid> ends successfully
#   DRYRUN=1           write bundle.sbatch but don't submit it
#   NEKRS_GPU_MPI=1    GPU-aware MPI (default 0, as in nrsqsub)
# Prints the job id. The job exits non-zero if any case fails, which ends a
# chain of DEPEND jobs there.
set -euo pipefail

: ${PROJ_ID:?PROJ_ID must be set}
: ${NEKRS_HOME:?NEKRS_HOME must be set}
: ${PARTITION:=batch}
: ${QOS:=}
: ${NEKRS_JITC_NTHREADS:=7}
: ${PREPARE_RESTART:=0}
: ${DEPEND:=}
: ${STOP_MARGIN:=300}
: ${DRYRUN:=0}
: ${NEKRS_GPU_MPI:=0}

if [ $# -lt 2 ]; then
    echo "Usage: PROJ_ID=<project> [PARTITION=batch] [QOS=] $0 <hh:mm> <dir>:<nodes> [<dir>:<nodes> ...]"
    exit 1
fi

time=$1
shift
IFS=: read -r hh mm <<< "$time"
walltime_s=$((10#$hh*3600 + 10#$mm*60))
gpu_per_node=8
cores_per_task=7

dirs=()
cases=()
nodes=()
total_nodes=0
for spec in "$@"; do
    dir=$(cd "${spec%:*}" && pwd -P)
    n=${spec##*:}
    pars=("$dir"/*.par)
    if [ ${#pars[@]} -ne 1 ] || [ ! -f "${pars[0]}" ]; then
        echo "$dir must hold exactly one .par file"
        exit 1
    fi
    casename=$(basename "${pars[0]}" .par)
    for f in $casename.udf $casename.re2; do
        if [ ! -f "$dir/$f" ]; then
            echo "Cannot find $dir/$f"
            exit 1
        fi
    done
    dirs+=("$dir")
    cases+=("$casename")
    nodes+=("$n")
    total_nodes=$((total_nodes + n))
done

striping_factor=$((total_nodes / 2))
if [ $striping_factor -lt 1 ]; then striping_factor=1; fi
if [ $striping_factor -gt 400 ]; then striping_factor=400; fi

SFILE=bundle.sbatch
cat > $SFILE <<EOF
#!/bin/bash
#SBATCH -A $PROJ_ID
#SBATCH -J nekRS_bundle
#SBATCH -o nekRS_bundle-%j.out
#SBATCH -t ${time}:00
#SBATCH -N $total_nodes
#SBATCH -p $PARTITION
#SBATCH --exclusive
EOF
if [ -n "$QOS" ]; then
    echo "#SBATCH -q $QOS" >> $SFILE
fi
cat >> $SFILE <<EOF

# The case's .udf ends each case cleanly, with a checkpoint, once the job is
# within $STOP_MARGIN s of its walltime.
export BUBBLE_STOP_AT=\$((\$(date +%s) + $walltime_s - $STOP_MARGIN))
EOF
cat >> $SFILE <<'EOF'

cd $SLURM_SUBMIT_DIR
echo Jobid: $SLURM_JOB_ID
echo Running on nodes $SLURM_JOB_NODELIST

module reset
module load Core/24.07
module load PrgEnv-gnu
module load gcc-native/13.2
module load craype-accel-amd-gfx90a
module load cray-mpich
module load rocm
module load cmake
module unload cray-libsci
module list
ulimit -s unlimited

# As in nrsqsub (Frontier). The GTL variables must be set before the UDF is
# compiled so it links the GPU transport layer.
export PE_MPICH_GTL_DIR_amd_gfx90a="-L${CRAY_MPICH_ROOTDIR}/gtl/lib"
export PE_MPICH_GTL_LIBS_amd_gfx90a="-lmpi_gtl_hsa"
export MPICH_GPU_SUPPORT_ENABLED=1
export OOGS_ENABLE_NBC_DEVICE=1
export MPICH_OFI_NIC_POLICY=NUMA
export FI_CXI_RX_MATCH_MODE=hybrid
export PMI_MMAP_SYNC_WAIT_TIME=600
export MPICH_MPIIO_STATS=1
export NEKRS_CACHE_BCAST=0
EOF
cat >> $SFILE <<EOF
export NEKRS_HOME=$NEKRS_HOME
export NEKRS_GPU_MPI=$NEKRS_GPU_MPI
export MPICH_MPIIO_HINTS="*:cray_cb_write_lock_mode=2:cray_cb_nodes_multiplier=4:striping_unit=16777216:striping_factor=${striping_factor}:romio_cb_write=enable:romio_ds_write=disable:romio_no_indep_rw=true"

prepare_restart=$PREPARE_RESTART
scripts=$(cd "$(dirname "$0")" && pwd -P)
gpu_per_node=$gpu_per_node
cores_per_task=$cores_per_task
jitc_nthreads=$NEKRS_JITC_NTHREADS
dirs=(${dirs[*]})
cases=(${cases[*]})
nodes=(${nodes[*]})
EOF
cat >> $SFILE <<'EOF'
jobid=$SLURM_JOB_ID
bin=$NEKRS_HOME/bin/nekrs
mapfile -t all_hosts < <(scontrol show hostnames "$SLURM_JOB_NODELIST")

# Leave out nodes whose GPUs are missing or fail HIP context creation; the
# case with the most nodes then runs on fewer nodes.
srun -N ${#all_hosts[@]} --ntasks-per-node=1 --gpus-per-task=$gpu_per_node \
    ./.gpucheck $gpu_per_node > gpucheck-$jobid 2>&1
cat gpucheck-$jobid
hosts=()
for h in "${all_hosts[@]}"; do
    if grep -q "^${h%%.*} OK" gpucheck-$jobid; then
        hosts+=("$h")
    fi
done
missing=$((${#all_hosts[@]} - ${#hosts[@]}))
if [ $missing -gt 0 ]; then
    big=0
    for k in "${!nodes[@]}"; do
        if [ ${nodes[k]} -gt ${nodes[big]} ]; then big=$k; fi
    done
    if [ ${nodes[big]} -le $missing ]; then
        echo "$(date) $missing bad nodes, too few left to run"
        exit 1
    fi
    nodes[big]=$((nodes[big] - missing))
    echo "$(date) left out $missing bad nodes, ${dirs[big]} runs on ${nodes[big]} nodes"
fi

if [ $prepare_restart = 1 ]; then
    run_dirs=() run_cases=() run_nodes=()
    for k in "${!dirs[@]}"; do
        d=${dirs[k]} casename=${cases[k]}
        end=$(awk -F= '/^endTime/ {print $2+0}' $d/$casename.par)
        last=$(tail -1 $d/data.csv 2>/dev/null | cut -d, -f1)
        if [ -n "$last" ] && awk -v a="$last" -v b="$end" 'BEGIN {exit !(a >= b - 1e-3)}'; then
            echo "$(date) $d reached endTime $end, skipping"
            continue
        fi
        if compgen -G "$d/${casename}0.f[0-9]*" > /dev/null && [ ! -L "$(ls $d/${casename}0.f[0-9]* | head -n 1)" ]; then
            $scripts/prepare-restart.sh $d || exit 1
        fi
        run_dirs+=("$d") run_cases+=("$casename") run_nodes+=("${nodes[k]}")
    done
    if [ ${#run_dirs[@]} -eq 0 ]; then
        echo "$(date) nothing left to run"
        exit 0
    fi
    dirs=("${run_dirs[@]}") cases=("${run_cases[@]}") nodes=("${run_nodes[@]}")
fi

# run_case <dir> <case name> <index of first node> <number of nodes>
run_case() {
    local dir=$1 casename=$2 first=$3 n=$4
    local ntasks=$((n * gpu_per_node))
    local nodelist
    nodelist=$(IFS=,; echo "${hosts[*]:first:n}")
    cd $dir
    echo "$nodelist" > nodes-$jobid
    {
        echo "Running on nodes $nodelist"
        echo "$(date) precompilation"
        NEKRS_JITC_NTHREADS=$jitc_nthreads srun -N 1 -n $gpu_per_node -w ${hosts[first]} \
            -c $cores_per_task --gpus-per-task=1 --gpu-bind=closest \
            $bin --setup $casename --backend HIP --device-id 0 --build-only $ntasks
        status=$?
        if [ $status -eq 0 ]; then
            echo "$(date) actual run"
            srun -N $n -n $ntasks -w $nodelist \
                -c $cores_per_task --gpus-per-task=1 --gpu-bind=closest \
                $bin --setup $casename --backend HIP --device-id 0
            status=$?
        fi
        echo "$(date) exit status $status"
    } > logfile-$jobid 2>&1
    echo $status > .status-$jobid
}

first=0
for k in "${!dirs[@]}"; do
    run_case ${dirs[k]} ${cases[k]} $first ${nodes[k]} &
    first=$((first + nodes[k]))
done
wait
echo "$(date) all cases done"
failed=0
for d in "${dirs[@]}"; do
    s=$(cat $d/.status-$jobid 2>/dev/null || echo 1)
    rm -f $d/.status-$jobid
    [ "$s" = 0 ] || { echo "$d failed (exit status $s)"; failed=1; }
done
exit $failed
EOF

hipcc -O2 -o .gpucheck.new "$(dirname "$0")/gpucheck-hip.cpp"
mv -f .gpucheck.new .gpucheck
if [ $DRYRUN = 1 ]; then
    echo "DRYRUN: wrote $SFILE"
elif [ -n "$DEPEND" ]; then
    sbatch --parsable --dependency=afterok:$DEPEND $SFILE
else
    sbatch --parsable $SFILE
fi
