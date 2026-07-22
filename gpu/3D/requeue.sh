#!/bin/bash
set -euo pipefail

gen_batch_file_header()
{
#--------------------------------------
: ${QUEUE:="EarlyAppAccess"}
: ${NEKRS_GPU_MPI:=0}
: ${NEKRS_BACKEND:="dpcpp"}
: ${RANKS_PER_NODE:=12}
: ${RANKS_FOR_BUILD:=12}
: ${CPU_BIND_LIST:="0-7,104-111:8-15,112-119:16-23,120-127:24-31,128-135:32-39,136-143:40-47,144-151:52-59,156-163:60-67,164-171:68-75,172-179:76-83,180-187:84-91,188-195:92-99,196-203"}
: ${OCCA_DPCPP_COMPILER_FLAGS:="-O3 -fsycl -fsycl-targets=intel_gpu_pvc -ftarget-register-alloc-mode=pvc:auto -fma"}
: ${ONEAPI_SDK:=""}
#--------------------------------------

source $NEKRS_HOME/bin/nrsqsub_utils
setup $# 1

TOTAL_RANKS=$(( nodes * RANKS_PER_NODE ))
gpus_per_node=6
tiles_per_gpu=2

chk_case $TOTAL_RANKS

#--------------------------------------
# Generate the submission script
SFILE="$1"
echo "#!/bin/bash" > $SFILE
echo "#PBS -A $PROJ_ID" >>$SFILE
# We set a custom jobname.
#echo "#PBS -N $jobname" >>$SFILE
echo "#PBS -l walltime=${time}:00" >>$SFILE
echo "#PBS -l select=$qnodes" >>$SFILE
echo "#PBS -l place=scatter" >>$SFILE
echo "#PBS -k doe" >>$SFILE
echo "#PBS -j oe" >>$SFILE
# This just doesn't work. Don't know why.
#echo "#PBS -M ben_li@berkeley.edu" >>$SFILE
#echo "#PBS -m bea" >>$SFILE
}

gen_batch_file_nekrs_invocation()
{
SFILE="$1"
echo "export TZ='/usr/share/zoneinfo/US/Central'" >> $SFILE

# job to "run" from your submission directory
echo "cd \$PBS_O_WORKDIR" >> $SFILE

echo "echo Jobid: \$PBS_JOBID" >>$SFILE
echo "echo Running on host \`hostname\`" >>$SFILE
echo "echo Running on nodes \`cat \$PBS_NODEFILE\`" >>$SFILE

if [ -n "$ONEAPI_SDK" ]; then
  echo "module load ${ONEAPI_SDK}" >> $SFILE
fi
echo "module load cmake" >> $SFILE
echo "module list" >> $SFILE

echo "export NEKRS_HOME=$NEKRS_HOME" >>$SFILE
echo "export NEKRS_GPU_MPI=$NEKRS_GPU_MPI" >>$SFILE
echo "export MPICH_GPU_SUPPORT_ENABLED=$NEKRS_GPU_MPI" >> $SFILE

echo "export OCCA_DPCPP_COMPILER_FLAGS=\"$OCCA_DPCPP_COMPILER_FLAGS\"" >> $SFILE

# Workaround for MPICH 52.2 see https://docs.alcf.anl.gov/aurora/known-issues/
#echo "unset MPIR_CVAR_CH4_COLL_SELECTION_TUNING_JSON_FILE" >> $SFILE
#echo "unset MPIR_CVAR_COLL_SELECTION_TUNING_JSON_FILE" >> $SFILE
#echo "unset MPIR_CVAR_CH4_POSIX_COLL_SELECTION_TUNING_JSON_FILE" >> $SFILE

# Cray MPI libfabric defaults
echo "export FI_CXI_RDZV_THRESHOLD=16384" >> $SFILE
echo "export FI_CXI_RDZV_EAGER_SIZE=2048" >> $SFILE
echo "export FI_CXI_DEFAULT_CQ_SIZE=131072" >> $SFILE
echo "export FI_CXI_DEFAULT_TX_SIZE=1024" >> $SFILE
echo "export FI_CXI_OFLOW_BUF_SIZE=12582912" >> $SFILE
echo "export FI_CXI_OFLOW_BUF_COUNT=3" >> $SFILE
echo "export FI_CXI_REQ_BUF_MIN_POSTED=6" >> $SFILE
echo "export FI_CXI_REQ_BUF_SIZE=12582912" >> $SFILE
echo "export FI_MR_CACHE_MAX_SIZE=-1" >> $SFILE
echo "export FI_MR_CACHE_MAX_COUNT=524288" >> $SFILE
echo "export FI_CXI_REQ_BUF_MAX_CACHED=0" >> $SFILE
echo "export FI_CXI_REQ_BUF_MIN_POSTED=6" >> $SFILE
# echo "export FI_CXI_RX_MATCH_MODE=hardware" >> $SFILE 
echo "export FI_CXI_RX_MATCH_MODE=hybrid" >> $SFILE # required by parRSB

CMD=.lhelper
if [ ! -e $CMD ]; then
  # There is no reason to recreate this file, tiles_per_gpu and gpus_per_node
  # are static.
  echo "#!/bin/bash" > $CMD
  echo "gpu_id=\$(((PALS_LOCAL_RANKID / ${tiles_per_gpu}) % ${gpus_per_node}))" >> $CMD
  echo "tile_id=\$((PALS_LOCAL_RANKID % ${tiles_per_gpu}))" >> $CMD
  echo "export ZE_AFFINITY_MASK=\$gpu_id.\$tile_id" >> $CMD
  echo "\"\$@\"" >> $CMD
  chmod u+x $CMD
fi

if [ $RUN_ONLY -eq 0 ]; then
  echo -e "\n# precompilation" >>$SFILE
  CMD_build="mpiexec --no-vni -n ${RANKS_FOR_BUILD} -ppn ${RANKS_FOR_BUILD} --cpu-bind=list:${CPU_BIND_LIST} -- ./${CMD} $bin --setup \${case_tmp} --backend ${NEKRS_BACKEND} --device-id 0 $extra_args --build-only \${ntasks_tmp}"
  add_build_CMD "$SFILE" "$CMD_build" "$TOTAL_RANKS"
fi

if [ $BUILD_ONLY -eq 0 ]; then
  link_neknek_logfile "$SFILE"
  echo -e "\n# actual run" >>$SFILE
  echo "mpiexec --no-vni -n ${TOTAL_RANKS} -ppn ${RANKS_PER_NODE} --cpu-bind=list:${CPU_BIND_LIST} -- ./${CMD} $bin --setup ${case} --backend ${NEKRS_BACKEND} --device-id 0 $extra_args" >> $SFILE
fi
}

qsub_aurora()
{
qsub -l "filesystems=home" -q $QUEUE "$@"
}

requeue="false"
if [[ "$1" == "-r" ]]; then
    requeue="true"
    shift
    # Don't clobber the original .batch file.
    mv $1.batch $1.pre-requeue.batch
    echo "Requeuing..."
fi

desc=$(cat JOB_DESC || echo "job")

gen_batch_file_header $1.batch $2 $3

echo "#PBS -N $desc" >> $1.batch
echo "cd $PWD" >> $1.batch
if [[ $requeue == "true" ]]; then
    echo "backup_folder=\$(date +%Y%m%d-%H%M)" >> $1.batch
    echo "last_checkpoint=\$(ls *.f0* | tail -n 1)" >> $1.batch
    echo "mkdir \$backup_folder && mv *.f0* data_*.csv \$backup_folder/" >> $1.batch
    echo "rm -f restart.fld && ln -s \$backup_folder/\$last_checkpoint restart.fld" >> $1.batch
    echo "sed -i 's/#\\?startFrom \\?=.*/startFrom = restart.fld/' $1.par" >> $1.batch
fi
echo "RUN_ONLY=1 NEKRS_HOME=$NEKRS_HOME QUEUE=$QUEUE ./requeue.sh -r $1 $2 $3" >> $1.batch

gen_batch_file_nekrs_invocation $1.batch

if [[ $requeue == "true" ]]; then
    if [[ -v PBS_JOBID ]]; then
        qsub_aurora -W depend=afternotok:"$PBS_JOBID" $1.batch
    else
        qsub_aurora $1.batch
        echo "Requeued outside of PBS context (manual requeue), job should start immediately."
    fi
else
    qsub_aurora $1.batch
fi
