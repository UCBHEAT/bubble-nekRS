#!/usr/bin/env bash
# Create one run directory per mesh under <dest>, e.g.
#   PYTHON=<python with gmsh> GMSH2NEK=<gmsh2nek> ./setup-runs.sh ~/answinter26/Sc1-pid
# makes <dest>/{2.25x,1.5x,1x,0.667x,0.444x} with the case files and the
# bubble.re2 mesh, and <dest>/common/pid.hpp for the udf. Meshes are given as
# name:Nk:dt:checkpointInterval, e.g.
#   ./setup-runs.sh ~/answinter26/Sc1-pid 1x:27:1e-4:0.5
# where the finest (interface) elements are cubes of h = 4/Nk, the element
# size of gpu/mesh-sensitivity's Nk x 2Nk x Nk mesh of the same name. Every
# size target of ../pid-centering-3D/generate_mesh.py's default mesh (made for
# h_iface = 0.2) is scaled by h/0.2, so the whole mesh is refined together.
# The fine cube around the bubble has an odd number of cells per side (the
# generator's even count plus one), so the bubble centre is inside an element
# rather than on a vertex with element edges along all three axes through it,
# where the gas core develops strong spurious currents at the physical
# density ratio (u_max 7 against 1.4 on the 1.5x mesh, and the 2.25x bubble
# breaks up).
# The interface width is eps = 1.5 h/N (N = 7), as interfaceWidthFactor = 1.5
# gives on the uniform meshes.
# SC sets the Schmidt number in each run's case.hpp (default 1), END_TIME the
# endTime in bubble.par (default 30).
#
# Resolutions are multiples of the Kolmogorov scale lambda_k = 0.0671 mm
# (lambda_k/d = 0.021536), using the mean unique GLL spacing h/N of the finest
# elements:
#   2.25x  -> Nk 12 (2.21x)    1x     -> Nk 27 (0.98x)    0.444x -> Nk 60 (0.442x)
#   1.5x   -> Nk 18 (1.47x)    0.667x -> Nk 40 (0.663x)
set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: [SC=1] [END_TIME=30] [PYTHON=python3] [GMSH2NEK=gmsh2nek] $0 <dest> [name:Nk:dt:checkpointInterval ...]"
    exit 1
fi

: ${SC:=1}
: ${END_TIME:=30}
src=$(cd "$(dirname "$0")" && pwd)
dest=$1
shift
specs=("$@")
if [ ${#specs[@]} -eq 0 ]; then
    specs=(2.25x:12:1e-4:0.5 1.5x:18:1e-4:0.5 1x:27:1e-4:0.5 0.667x:40:5e-5:1 0.444x:60:2.5e-5:1)
fi

mkdir -p "$dest/common"
cp "$src/../common/pid.hpp" "$dest/common/"

for spec in "${specs[@]}"; do
    IFS=: read -r name nk dt ckpt <<< "$spec"
    dir=$dest/$name
    if [ -e "$dir/bubble.par" ]; then
        echo "$dir already set up, skipping"
        continue
    fi
    mkdir -p "$dir"
    # h = 4/Nk, the mesh size targets (generate_mesh.py defaults times h/0.2)
    # and eps = 1.5 h/N.
    # The fine cube has n = 2 ceil(0.76/h) + 1 cells (0.76 = the generator's
    # r_iface_out + iface_margin), half-width n h/2.
    read -r h h_wake_cross h_wake_near h_wake_far h_far eps fine_half <<< "$(awk -v nk="$nk" 'BEGIN {
        h = 4/nk; s = h/0.2
        c = 0.76/h - 1e-9; n = int(c); if (c > n) n++; n = 2*n + 1
        printf "%.17g %.17g %.17g %.17g %.17g %.8g %.17g\n", h, 0.3*s, 0.25*s, 0.5*s, 1.25*s, 1.5*h/7, n*h/2 }')"
    MESH_DIR=$dir "$src/../pid-centering-3D/mesh" --name bubble --h_iface "$h" \
        --h_wake_cross "$h_wake_cross" --h_wake_near "$h_wake_near" --h_wake_far "$h_wake_far" \
        --h_far "$h_far" --fine_half "$fine_half" --rho_ratio 3691.4 --max_elements 1000000 \
        > "$dir/mesh.log" 2>&1 \
        || { tail -20 "$dir/mesh.log"; echo "$dir: meshing failed (see $dir/mesh.log)"; exit 1; }
    cp "$src"/{bubble.udf,bubble.oudf,bubble.usr,customhooks.hpp,util.hpp} "$dir"/
    sed -e "s/^static double Sc = .*/static double Sc = $SC;/" -e "s/(Sc=[0-9.]*)/(Sc=$SC)/" \
        "$src/case.hpp" > "$dir/case.hpp"
    sed -e "s/^dt = .*/dt = $dt/" -e "s/^checkpointInterval = .*/checkpointInterval = $ckpt/" \
        -e "s/^endTime = .*/endTime = $END_TIME/" \
        -e "s/^interfaceWidthValue = .*/interfaceWidthValue = $eps/" \
        -e "s|This value is for the 1x mesh (h = 4/27).|This value is for the $name mesh (h = 4/$nk).|" \
        "$src/bubble.par" > "$dir/bubble.par"
    nel=$("${PYTHON:-python3}" -c "import json, sys; print(json.load(open(sys.argv[1]))['n_elements'])" "$dir/bubble.plan.json")
    echo "$dir: h = 4/$nk, $nel elements, eps = $eps, dt = $dt, checkpointInterval = $ckpt, Sc = $SC, endTime = $END_TIME"
done
