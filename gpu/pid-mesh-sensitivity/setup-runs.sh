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
# A fifth field, name:Nk:dt:checkpointInterval:Nk_outer, refines only the fine
# cube to h = 4/Nk and sizes the rest of the mesh (wake column and far field)
# as for Nk_outer, e.g. 0.088x:304:1.5625e-6:1:90 (4.5M elements against
# 13.9M for the whole mesh at Nk = 304).
# SC sets the Schmidt number in each run's case.hpp (default 1), END_TIME the
# endTime in bubble.par (default 30). MAX_ELEMENTS (default 1000000) and
# MAX_MEMORY_GB (default 8) are generate_mesh.py's caps on the element count
# and its memory; 0.198x (1.37M elements) needs both raised.
# SPECIES="20 1150" adds species at these Schmidt numbers to every run, as
# scalars c20 and c1150 next to c, on the same flow (see bubble.udf; output
# in data-c20.csv and data-c1150.csv).
# MESH_ONLY=1 only makes the meshes (e.g. while FROM's checkpoint is still to
# come): a later run of the script reuses a run directory's mesh.
# FROM=<run dir> and START=<time> start every run at t = START from FROM's
# checkpoint at that time instead of from rest: the flow and level set
# (interpolated onto the new mesh, startFrom = <file>+int, unless it is
# FROM's) and the PID state (FROM's data.csv row at START), e.g. to run a
# finer mesh or another Schmidt number from a developed flow. Runs with FROM's
# species (same SC and SPECIES) continue them. Otherwise the species start
# afresh from FROM's c, their deficits scaled from FROM's Schmidt number to
# theirs ([CASEDATA] speciesStart and speciesStartSc, see bubble.udf), and
# Sh then reads high for a long time while the slowly renewed liquid behind
# the bubble adjusts (Sc = 20 from Sc = 4 on 0.444x: 1.2% after 1-1.5 time
# units, 0.6% after 2-2.5), so it is better to carry every species from
# t = 0 on the coarsest mesh and start the finer ones from it.
#
# Resolutions are multiples of the Kolmogorov scale lambda_k = 0.0671 mm
# (lambda_k/d = 0.021536), using the mean unique GLL spacing h/N of the finest
# elements:
#   2.25x  -> Nk 12 (2.21x)    1x     -> Nk 27 (0.98x)    0.444x -> Nk 60 (0.442x)
#   1.5x   -> Nk 18 (1.47x)    0.667x -> Nk 40 (0.663x)
set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: [SC=1] [SPECIES=\"<Sc> ...\"] [END_TIME=30] [MAX_ELEMENTS=1000000] [MAX_MEMORY_GB=8] [FROM=<run dir> START=<time>] [PYTHON=python3] [GMSH2NEK=gmsh2nek] $0 <dest> [name:Nk:dt:checkpointInterval[:Nk_outer] ...]"
    exit 1
fi

: ${SC:=1}
: ${END_TIME:=30}
: ${MAX_ELEMENTS:=1000000}
: ${MAX_MEMORY_GB:=8}
: ${SPECIES:=}
: ${FROM:=}
if [ -n "$FROM" ]; then
    : ${START:?START must be set with FROM}
    FROM=$(cd "$FROM" && pwd)
fi
src=$(cd "$(dirname "$0")" && pwd)
dest=$1
shift
specs=("$@")
if [ ${#specs[@]} -eq 0 ]; then
    specs=(2.25x:12:1e-4:0.5 1.5x:18:1e-4:0.5 1x:27:1e-4:0.5 0.667x:40:5e-5:1 0.444x:60:2.5e-5:1)
fi

mkdir -p "$dest/common"
cp "$src/../common/pid.hpp" "$dest/common/"

# add_species <par file>: the SPECIES scalars, with sections like [SCALAR C].
add_species() {
    awk -v species="$SPECIES" '
        /^scalars = / { n = split(species, sc, " "); for (i = 1; i <= n; i++) $0 = $0 ", c" sc[i] }
        /^\[SCALAR C\]$/ { insec = 1 }
        insec && /^$/ {
            insec = 0
            print
            for (i = 1; i <= n; i++) {
                print "[SCALAR C" sc[i] "]"
                print "# As c, at Sc = " sc[i] " (the udf takes it from the name)"
                for (j = 1; j <= nk; j++) print keys[j]
                print ""
            }
            next
        }
        insec && !/^#/ && !/^\[/ { keys[++nk] = $0 }
        { print }' "$1" > "$1.tmp"
    mv "$1.tmp" "$1"
}

# start_from <run dir>: start the run at t = START from FROM (see above).
start_from() {
    local dir=$1 f t ckpt="" first opt=""
    for f in "$FROM"/bubble0.f[0-9]* "$FROM"/part*/bubble0.f[0-9]*; do
        [ -f "$f" ] || continue
        t=$(head -c 132 "$f" | awk '{print $8}')
        if awk -v a="$t" -v b="$START" 'BEGIN {exit !(a > b - 1e-6 && a < b + 1e-6)}'; then
            ckpt=$(readlink -f "$f")
            break
        fi
    done
    [ -n "$ckpt" ] || { echo "$FROM has no checkpoint at t = $START"; exit 1; }
    # One copy in <dest> for all the runs, with the mesh coordinates, which an
    # interpolated restart needs and nekRS writes only into each job's first
    # checkpoint.
    fld=start-$(basename "$(dirname "$FROM")")-$(basename "$FROM")-t$START.fld
    if [ ! -e "$dest/$fld" ]; then
        if head -c 132 "$ckpt" | awk '{exit !($12 ~ /^X/)}'; then
            cp "$ckpt" "$dest/$fld"
        else
            first=$(ls "$(dirname "$ckpt")"/bubble0.f[0-9]* | head -n 1)
            "${PYTHON:-python3}" "$src/../../scripts/add-coords.py" "$ckpt" "$first" "$dest/$fld"
        fi
    fi
    ln -sfn "../$fld" "$dir/start.fld"
    cmp -s "$FROM/bubble.re2" "$dir/bubble.re2" || opt=+int
    # +u+p+s: read the fields but not the (single precision) coordinates.
    fromsc=$(sed -n 's/^static double Sc = \(.*\);/\1/p' "$FROM/case.hpp")
    sed -i "s/^#\?startFrom *=.*/startFrom = start.fld$opt+u+p+s/" "$dir/bubble.par"
    if [ "$fromsc" = "$SC" ] && [ "$(grep '^scalars' "$FROM/bubble.par")" = "$(grep '^scalars' "$dir/bubble.par")" ]; then
        species="its species continue"
    else
        species="species afresh from its c at Sc = $fromsc"
        sed -i -e "s/^#\?speciesStart *=.*/speciesStart = $START/" \
            -e "s/^#\?speciesStartSc *=.*/speciesStartSc = $fromsc/" "$dir/bubble.par"
    fi
    # data.csv times are rounded to 4 decimals.
    awk -F, -v t="$START" 'NR == 1 || ($1 > t - 6e-5 && $1 < t + 6e-5)' "$FROM/data.csv" > "$dir/data.csv"
    [ $(wc -l < "$dir/data.csv") -eq 2 ] || { echo "$FROM/data.csv has no row at t = $START"; exit 1; }
    echo "$dir: starts at t = $START from $ckpt${opt:+ (interpolated)}, $species"
}

for spec in "${specs[@]}"; do
    IFS=: read -r name nk dt ckpt outer <<< "$spec"
    outer=${outer:-$nk}
    dir=$dest/$name
    if [ -e "$dir/bubble.par" ]; then
        echo "$dir already set up, skipping"
        continue
    fi
    mkdir -p "$dir"
    # h = 4/Nk, the mesh size targets (generate_mesh.py defaults times h/0.2,
    # or (4/Nk_outer)/0.2 outside the fine cube) and eps = 1.5 h/N.
    # The fine cube has n = 2 ceil(0.76/h) + 1 cells (0.76 = the generator's
    # r_iface_out + iface_margin), half-width n h/2.
    read -r h h_wake_cross h_wake_near h_wake_far h_far eps fine_half <<< "$(awk -v nk="$nk" -v outer="$outer" 'BEGIN {
        h = 4/nk; s = (4/outer)/0.2
        c = 0.76/h - 1e-9; n = int(c); if (c > n) n++; n = 2*n + 1
        printf "%.17g %.17g %.17g %.17g %.17g %.8g %.17g\n", h, 0.3*s, 0.25*s, 0.5*s, 1.25*s, 1.5*h/7, n*h/2 }')"
    if [ -e "$dir/bubble.re2" ] && [ -e "$dir/bubble.plan.json" ]; then
        echo "$dir: reusing its mesh"
    else
        MESH_DIR=$dir "$src/../pid-centering-3D/mesh" --name bubble --h_iface "$h" \
            --h_wake_cross "$h_wake_cross" --h_wake_near "$h_wake_near" --h_wake_far "$h_wake_far" \
            --h_far "$h_far" --fine_half "$fine_half" --rho_ratio 3691.4 --max_elements "$MAX_ELEMENTS" --max_memory_gb "$MAX_MEMORY_GB" \
            > "$dir/mesh.log" 2>&1 \
            || { tail -20 "$dir/mesh.log"; echo "$dir: meshing failed (see $dir/mesh.log)"; exit 1; }
    fi
    if [ "${MESH_ONLY:-0}" = 1 ]; then
        echo "$dir: mesh only, $("${PYTHON:-python3}" -c "import json, sys; print(json.load(open(sys.argv[1]))['n_elements'])" "$dir/bubble.plan.json") elements"
        continue
    fi
    cp "$src"/{bubble.udf,bubble.oudf,bubble.usr,customhooks.hpp,util.hpp} "$dir"/
    sed -e "s/^static double Sc = .*/static double Sc = $SC;/" -e "s/(Sc=[0-9.]*)/(Sc=$SC)/" \
        "$src/case.hpp" > "$dir/case.hpp"
    sed -e "s/^dt = .*/dt = $dt/" -e "s/^checkpointInterval = .*/checkpointInterval = $ckpt/" \
        -e "s/^endTime = .*/endTime = $END_TIME/" \
        -e "s/^interfaceWidthValue = .*/interfaceWidthValue = $eps/" \
        -e "s|This value is for the 1x mesh (h = 4/27).|This value is for the $name mesh (h = 4/$nk).|" \
        "$src/bubble.par" > "$dir/bubble.par"
    if [ -n "$SPECIES" ]; then
        add_species "$dir/bubble.par"
    fi
    if [ -n "$FROM" ]; then
        start_from "$dir"
    fi
    nel=$("${PYTHON:-python3}" -c "import json, sys; print(json.load(open(sys.argv[1]))['n_elements'])" "$dir/bubble.plan.json")
    rest=""; [ "$outer" = "$nk" ] || rest=" (rest of the mesh as for Nk = $outer)"
    echo "$dir: h = 4/$nk$rest, $nel elements, eps = $eps, dt = $dt, checkpointInterval = $ckpt, Sc = $SC${SPECIES:+ (and $SPECIES)}, endTime = $END_TIME"
done
