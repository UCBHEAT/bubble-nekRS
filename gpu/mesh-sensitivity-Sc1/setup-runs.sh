#!/usr/bin/env bash
# Create one run directory per mesh resolution under <dest>, e.g.
#   ./setup-runs.sh ~/answinter26/Sc1
# makes ~/answinter26/Sc1/{2.25x,1.5x,1x} with the case files and .re2 mesh.
#
# Resolutions are multiples of the Kolmogorov scale lambda_k = 0.0671 mm
# (lambda_k/d = 0.021536), using the mean unique GLL spacing dx = h/N with
# N = 7 and element size h = 4/Nelx:
#   2.25x -> 12x24x12 (2.21x)
#   1.5x  -> 18x36x18 (1.47x)
#   1x    -> 27x54x27 (0.98x)
set -euo pipefail

if [ $# -ne 1 ]; then
    echo "Usage: $0 <dest>"
    exit 1
fi

src=$(cd "$(dirname "$0")" && pwd)
dest=$1

for spec in 2.25x:12 1.5x:18 1x:27; do
    name=${spec%%:*}
    nx=${spec##*:}
    dir=$dest/$name
    if [ -e "$dir/bubble3d.par" ]; then
        echo "$dir already set up, skipping"
        continue
    fi
    mkdir -p "$dir"
    cp "$src"/{bubble3d.par,bubble3d.udf,bubble3d.oudf,case.hpp,customhooks.hpp,util.hpp,mesh} "$dir"/
    sed -e "s/^-27 -54 -27 .*/-$nx -$((2*nx)) -$nx            Nelx Nely Nelz ($name Kolmogorov)/" \
        "$src/bubble3d.box" > "$dir/bubble3d.box"
    (cd "$dir" && ./mesh > mesh.log)
    echo "$dir: ${nx}x$((2*nx))x${nx} elements"
done
