#!/usr/bin/env bash
# Submit several run directories of this study as one Polaris PBS job, each
# case on its own set of nodes, e.g. from the directory holding the runs:
#   PROJ_ID=nek-vf QUEUE=prod ./submit-bundle.sh 03:00 1x:7 1.5x:2 2.25x:1
# The prod queue needs at least 10 nodes, which is far more than any single
# case of this study needs, so all meshes share one job.
#
# The job environment follows $NEKRS_HOME/bin/nrsqsub_polaris. The batch
# script is written to bundle.batch, the PBS log to nekRS_bundle.e<jobid>,
# and each case's nekRS output to <dir>/logfile-<jobid>.
set -euo pipefail

: ${PROJ_ID:?PROJ_ID must be set}
: ${QUEUE:?QUEUE must be set}
: ${NEKRS_HOME:?NEKRS_HOME must be set}
: ${NEKRS_JITC_NTHREADS:=7}

if [ $# -lt 2 ]; then
    echo "Usage: PROJ_ID=<project> QUEUE=<queue> $0 <hh:mm> <dir>:<nodes> [<dir>:<nodes> ...]"
    exit 1
fi

casename=bubble3d
time=$1
shift
gpu_per_node=4
cores_per_numa=8

dirs=()
nodes=()
total_nodes=0
for spec in "$@"; do
    dir=$(cd "${spec%:*}" && pwd)
    n=${spec##*:}
    for f in $casename.par $casename.udf $casename.re2; do
        if [ ! -f "$dir/$f" ]; then
            echo "Cannot find $dir/$f"
            exit 1
        fi
    done
    dirs+=("$dir")
    nodes+=("$n")
    total_nodes=$((total_nodes + n))

    # Same GPU binding helper as nrsqsub_polaris.
    cat > "$dir/.lhelper" <<EOF
#!/bin/bash
gpu_id=\$(($gpu_per_node - 1 - \${PMI_LOCAL_RANK} % $gpu_per_node))
export CUDA_VISIBLE_DEVICES=\$gpu_id
\$*
EOF
    chmod 755 "$dir/.lhelper"
done

striping_factor=$((total_nodes / 2))
if [ $striping_factor -lt 1 ]; then striping_factor=1; fi
if [ $striping_factor -gt 128 ]; then striping_factor=128; fi

SFILE=bundle.batch
cat > $SFILE <<EOF
#!/bin/bash
#PBS -A $PROJ_ID
#PBS -N nekRS_bundle
#PBS -q $QUEUE
#PBS -l walltime=${time}:00
#PBS -l filesystems=home:eagle:grand
#PBS -l select=$total_nodes:system=polaris
#PBS -l place=scatter
#PBS -k doe
#PBS -j eo
EOF
cat >> $SFILE <<'EOF'

cd $PBS_O_WORKDIR
echo Jobid: $PBS_JOBID
echo Running on host `hostname`
echo Running on nodes `cat $PBS_NODEFILE`

module restore
module use /soft/modulefiles
module swap PrgEnv-nvidia PrgEnv-gnu
module load cudatoolkit-standalone/13.0.1
module load cuda/13.0
module load gcc-native/14
module load craype-x86-milan craype-accel-nvidia80
module load spack-pe-base cmake
module unload darshan
module list
nvidia-smi
ulimit -s unlimited

EOF
cat >> $SFILE <<EOF
export NEKRS_HOME=$NEKRS_HOME
export NEKRS_GPU_MPI=1
export MPICH_MPIIO_HINTS=*:striping_unit=16777216:striping_factor=${striping_factor}:romio_cb_write=enable:romio_ds_write=disable:romio_no_indep_rw=true
export MPICH_MPIIO_STATS=1
export NEKRS_CACHE_BCAST=0
export NEKRS_LOCAL_TMP_DIR=/local/scratch
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_OFI_NIC_POLICY=NUMA
export MPIR_CVAR_CH4_OFI_ENABLE_RMA=0

casename=$casename
gpu_per_node=$gpu_per_node
cores_per_numa=$cores_per_numa
jitc_nthreads=$NEKRS_JITC_NTHREADS
EOF
cat >> $SFILE <<'EOF'
jobid=${PBS_JOBID%%.*}
bin=$NEKRS_HOME/bin/nekrs
mapfile -t hosts < <(awk '!seen[$0]++' $PBS_NODEFILE)

# run_case <dir> <index of first node> <number of nodes>
run_case() {
    local dir=$1 first=$2 n=$3
    local ntasks=$((n * gpu_per_node))
    cd $dir
    printf "%s\n" "${hosts[@]:first:n}" > nodes-$jobid
    {
        echo "Running on nodes $(tr '\n' ' ' < nodes-$jobid)"
        echo "$(date) precompilation"
        NEKRS_JITC_NTHREADS=$jitc_nthreads mpiexec --hostfile nodes-$jobid \
            -n $gpu_per_node -ppn $gpu_per_node -d $cores_per_numa --cpu-bind depth \
            ./.lhelper $bin --setup $casename --backend CUDA --device-id 0 --build-only $ntasks
        status=$?
        if [ $status -eq 0 ]; then
            echo "$(date) actual run"
            mpiexec --hostfile nodes-$jobid \
                -n $ntasks -ppn $gpu_per_node -d $cores_per_numa --cpu-bind depth \
                ./.lhelper $bin --setup $casename --backend CUDA --device-id 0
            status=$?
        fi
        echo "$(date) exit status $status"
    } > logfile-$jobid 2>&1
}

EOF
first=0
for k in "${!dirs[@]}"; do
    echo "run_case ${dirs[k]} $first ${nodes[k]} &" >> $SFILE
    first=$((first + nodes[k]))
done
cat >> $SFILE <<'EOF'
wait
echo "$(date) all cases done"
EOF

qsub -q $QUEUE $SFILE
