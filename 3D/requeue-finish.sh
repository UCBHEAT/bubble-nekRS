#!/bin/bash
set -euo pipefail

# Backup last run's checkpoints if needed
if compgen -G "*.f0*"; then
    backup_folder=$(date +%Y%m%d-%H%M)
    last_checkpoint=$(ls *.f0* | tail -n 1)
    echo "Backup last set of checkpoints to $backup_folder"
    mkdir $backup_folder && mv *.f0* *.csv $backup_folder/
fi

echo "Making symlinks..."
# The CSVs are nicely numbered in sequence.
for file in *-*/*.csv; do
    ln -s $file $(basename $file)
done
# The .f##### files are not, annoyingly.
i=0
for file in $(find . -maxdepth 2 -regextype posix-extended -regex "\\./[0-9]{8}-[0-9]{4}/bubble3d0\.f[0-9]{5}" | sort); do
    i=$((i+1))
    ipad=$(printf "%05d" $i)
    ln -s $file bubble3d0.f$ipad
done

echo "Writing bubble3d.nek5000..."
cat >bubble3d.nek5000 <<END
 filetemplate: bubble3d%01d.f%05d
 firsttimestep: 1
 numtimesteps:            $i
END
