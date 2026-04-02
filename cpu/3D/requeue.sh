#!/bin/bash
set -euo pipefail

requeue="false"
if [[ "$1" == "-r" ]]; then
    requeue="true"
    shift
    # Don't clobber the original .batch file.
    mv $1.batch $1.pre-requeue.batch
    echo "Requeuing..."
fi

desc=$(cat JOB_DESC || echo "job")

echo -e "$1\n$PWD" > SESSION.NAME
echo "#PBS -A CFDAPM" > $1.batch
echo "#PBS -l select=$2:mpiprocs=128" >> $1.batch
echo "#PBS -l walltime=$3" >> $1.batch
echo "#PBS -N $desc" >> $1.batch
echo "#PBS -q backfill" >> $1.batch
echo "set -euo pipefail" >> $1.batch
echo "module load intel-oneapi-mpi/2021.15" >> $1.batch
echo "cd $PWD" >> $1.batch
if [[ $requeue == "true" ]]; then
    echo "backup_folder=\$(date +%Y%m%d-%H%M)" >> $1.batch
    echo "last_checkpoint=\$(ls *.f0* | tail -n 1)" >> $1.batch
    echo "mkdir \$backup_folder && mv *.f0* *.csv \$backup_folder/" >> $1.batch
    echo "rm -f restart.fld && ln -s \$backup_folder/\$last_checkpoint restart.fld" >> $1.batch
fi
echo "./requeue.sh -r $1 $2 $3" >> $1.batch
echo "mpirun ./nek5000" >> $1.batch
if [[ $requeue == "true" ]]; then
    if [[ -v PBS_JOBID ]]; then
        qsub -W depend=afternotok:"$PBS_JOBID" $1.batch
    else
        qsub $1.batch
        echo "Requeued outside of PBS context (manual requeue), job should start immediately."
    fi
else
    qsub $1.batch
    echo "Wait for nek5000 to start, then manually edit your .par file to set startFrom=restart.fld."
fi
