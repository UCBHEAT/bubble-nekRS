#!/usr/bin/env bash
# Create one run directory per mesh resolution under <dest>, e.g.
#   ./setup-runs.sh ~/answinter26/Sc1
# makes ~/answinter26/Sc1/{2.25x,1.5x,1x} with the case files and .re2 mesh.
# Other meshes are given as name:Nelx:dt:checkpointInterval, e.g.
#   GENBOX=~/answinter26/tools/genbox-maxnel1.5M ./setup-runs.sh ~/answinter26/Sc1 \
#       0.667x:40:5e-5:0.5 0.444x:60:2.5e-5:0.5 0.296x:90:1.25e-5:0.25
# (genbox must be built with MAXNEL >= Nelx*2*Nelx*Nelx).
#
# Resolutions are multiples of the Kolmogorov scale lambda_k = 0.0671 mm
# (lambda_k/d = 0.021536), using the mean unique GLL spacing dx = h/N with
# N = 7 and element size h = 4/Nelx:
#   2.25x  -> 12x24x12   (2.21x)     0.667x -> 40x80x40    (0.663x)
#   1.5x   -> 18x36x18   (1.47x)     0.444x -> 60x120x60   (0.442x)
#   1x     -> 27x54x27   (0.98x)     0.296x -> 90x180x90   (0.295x)
set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 <dest> [name:Nelx:dt:checkpointInterval ...]"
    exit 1
fi

src=$(cd "$(dirname "$0")" && pwd)
dest=$1
shift
specs=("$@")
if [ ${#specs[@]} -eq 0 ]; then
    specs=(2.25x:12:1e-4:0.1 1.5x:18:1e-4:0.1 1x:27:1e-4:0.1)
fi

for spec in "${specs[@]}"; do
    IFS=: read -r name nx dt ckpt <<< "$spec"
    dir=$dest/$name
    if [ -e "$dir/bubble3d.par" ]; then
        echo "$dir already set up, skipping"
        continue
    fi
    mkdir -p "$dir"
    cp "$src"/{bubble3d.udf,bubble3d.oudf,case.hpp,customhooks.hpp,util.hpp,mesh} "$dir"/
    sed -e "s/^dt = .*/dt = $dt/" -e "s/^checkpointInterval = .*/checkpointInterval = $ckpt/" \
        "$src/bubble3d.par" > "$dir/bubble3d.par"
    sed -e "s/^-27 -54 -27 .*/-$nx -$((2*nx)) -$nx            Nelx Nely Nelz ($name Kolmogorov)/" \
        "$src/bubble3d.box" > "$dir/bubble3d.box"
    (cd "$dir" && ./mesh > mesh.log)
    echo "$dir: ${nx}x$((2*nx))x${nx} elements, dt = $dt, checkpointInterval = $ckpt"
done
