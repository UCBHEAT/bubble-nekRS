# pid-mesh-sensitivity

Sc = 1 mesh sensitivity study of the Sherwood number, like
`gpu/mesh-sensitivity`, but with the bubble held in place by the PID
controller of `gpu/pid-centering-3D` and fresh liquid flowing in, on graded
gmsh meshes. It uses the same five refinement levels, 2.25x, 1.5x, 1x,
0.667x and 0.444x the Kolmogorov scale, where the level of a mesh is set by
its finest elements (around the bubble).

In `gpu/mesh-sensitivity` the bubble rose through a fully periodic 4x8x4 box
and kept meeting its own species-depleted wake, while a soft source
restored c in the liquid. It then drifted sideways at a time set by
numerical noise, left its depletion trail and Sh rose, so the meshes could
only be compared in windows where all of them were in the same phase. Here
the bubble never meets depleted liquid: the depleted wake leaves through the
outflow, and the liquid enters saturated.

What comes from where:
- Physics, species model and bubble centering from `gpu/pid-centering-3D`
  (see its README and `gpu/pid-centering`'s): a spherical bubble (D = 1)
  held at the origin by a PID controller that accelerates the reference
  frame (`../common/pid.hpp`), in a 6x15x6 domain that is periodic in x and
  z, with the liquid entering at the top at the frame velocity and leaving
  at the bottom. c is a passive scalar with the liquid diffusivity 1/Pe
  everywhere, c = 1 at the inflow and a hard sink (c = 0 where psi < 0.5)
  in the gas; no soft source. mdot comes from the species budget, and
  `MTC = mdot/(area * c_inflow)`, `Sh = MTC Pe`.
- Nondimensional numbers, numerics and job handling from
  `gpu/mesh-sensitivity`: Re = 231.4, Sc = 1 (Pe = 231.4), Fr = 1.633,
  We = 2.667, muratio = 148.1 and the physical density ratio 1/0.0002709 =
  3691; N = 7; the same dt for the same finest element size; TLSR and CLSR
  every 0.01 and 0.001 time units; data.csv every 0.1; the BUBBLE_STOP_AT
  clean stop and the job scripts in `../common`.
- Meshes from `../pid-centering-3D/generate_mesh.py` (the Cartesian core):
  a uniform cube of fine elements around the bubble, a wake column below it,
  and elements that grow by at most 1.5x per element towards the boundaries.

## Meshes

The Kolmogorov scale is lambda_K = 0.0671 mm (lambda_K/d = 0.021536). As in
`gpu/mesh-sensitivity`, the resolution is the mean unique GLL spacing
dx = h/N at N = 7, here of the finest elements, which are cubes of the same
size h = 4/Nk as the uniform elements of the gpu/mesh-sensitivity run of the
same name (Nk x 2Nk x Nk elements). Every size target of the PID case's
default mesh (made for h = 0.2: wake 1.25-2.5 h, far field 6.25 h) is scaled
by h/0.2, so the whole mesh is refined together:

| Run    | h (finest) | dx/lambda_K | fine cube | elements (x y z) | GLL points (E N^3) | uniform mesh of the same h | largest edge | Hmax/Hmin |
|--------|------------|-------------|-----------|------------------|--------------------|----------------------------|--------------|-----------|
| 2.25x  | 4/12       | 2.21        | 6^3, half-width 1.0   | 12x22x12 = 3168   | 1.09M | 3456 (1.2M)    | 1.63 | 2.91 |
| 1.5x   | 4/18       | 1.47        | 8^3, half-width 0.889 | 16x32x16 = 8192   | 2.81M | 11664 (4.0M)   | 0.89 | 3.20 |
| 1x     | 4/27       | 0.98        | 12^3, half-width 0.889 | 22x44x22 = 21296 | 7.30M | 39366 (13.5M)  | 0.81 | 4.11 |
| 0.667x | 4/40       | 0.663       | 16^3, half-width 0.8  | 28x62x28 = 48608  | 16.7M | 128000 (43.9M) | 0.60 | 5.27 |
| 0.444x | 4/60       | 0.442       | 24^3, half-width 0.8  | 40x92x40 = 147200 | 50.5M | 432000 (148M)  | 0.38 | 5.35 |

The fine cube is the smallest multiple of h that covers the interface zone
r <= 0.76, so its half-width varies a little with h. The interface width is
eps = 1.5 h/N ([LVLSET] interfaceWidthValue), which is what
interfaceWidthFactor = 1.5 gives on the uniform meshes; the lvlSet default
would take h from the largest element.

## Numerics

As in `gpu/mesh-sensitivity`, except:
- `customProperties` calls `scalar->mueSVV()`, so the scalar SVV in
  `[SCALAR *]` is active (gpu/mesh-sensitivity's scalars ran without it; see
  `../pid-centering/README.md`).
- `[TLSR]`/`[CLSR] maximumSteps = 200/80`. lvlSet takes the pseudo-time step
  from the smallest element and the distance to cover from the largest, so
  on a uniform mesh its default runs 25 (N + 1) = 200 TLSR and 10 (N + 1) = 80
  CLSR pseudo-steps. On these meshes (Hmax/Hmin 2.9-5.4) that would be
  580-1000 and 230-430, 3-5x the cost, only to redistance the far field.
  The caps keep the uniform-mesh counts, so near the bubble the reinit is
  the same as in gpu/mesh-sensitivity.
- The buoyancy reference density is the liquid's (open domain), not the
  domain average.

| Run    | dt      | TLSR/CLSR every | checkpoints |
|--------|---------|-----------------|-------------|
| 2.25x  | 1e-4    | 100/10 steps    | 0.5         |
| 1.5x   | 1e-4    | 100/10 steps    | 0.5         |
| 1x     | 1e-4    | 100/10 steps    | 0.5         |
| 0.667x | 5e-5    | 200/20 steps    | 1           |
| 0.444x | 2.5e-5  | 400/40 steps    | 1           |

## Running on Polaris

The meshes need a python with gmsh and numpy, and gmsh2nek:

```
python3 -m venv --system-site-packages ~/answinter26/tools/gmsh-venv   # conda python
~/answinter26/tools/gmsh-venv/bin/pip install --pre -i https://gmsh.info/python-packages-dev-nox "gmsh==4.15.0.dev1+nox"
cd ~/Nek5000/tools && bin_nek_tools=~/answinter26/tools ./maketools gmsh2nek
```

(The PyPI `gmsh` wheel needs libGLU, which the Polaris login nodes lack; the
`nox` build of gmsh 4.15, the version `../pid-centering-3D` was tested
with, does not.) Then

```
PYTHON=~/answinter26/tools/gmsh-venv/bin/python GMSH2NEK=~/answinter26/tools/gmsh2nek \
    ./setup-runs.sh ~/answinter26/Sc1-pid
```

makes the five run directories (name:Nk:dt:checkpointInterval specs select
others; SC sets the Schmidt number), each with its mesh (`bubble.re2`,
`mesh.log`, `bubble.plan.json`), and `common/pid.hpp` next to them for the
udf. Meshing takes 5-40 s per run on a login node.

Runs are submitted and restarted with the job scripts in `../common` (see
`../mesh-sensitivity/README.md`), e.g.

```
cd ~/answinter26/Sc1-pid
PROJ_ID=nek-vf QUEUE=prod PREPARE_RESTART=1 <common>/submit-bundle.sh 03:00 \
    0.444x:12 0.667x:4 1x:2 1.5x:1 2.25x:1
```

A restart takes the PID integral and frame velocity from the data.csv row at
the restart time, which the clean stop always writes.

## Output

data.csv has one row every 0.1 time units and at each job's last step:
- `mdot`, `MTC`, `Sh`: interface transfer from the species budget over the
  row's interval (see `../pid-centering/README.md`);
- `mdot_sink`, `Sh_sink`: the same from what the hard sink zeroed, which is
  how gpu/mesh-sensitivity measured mdot; with BDF2 it reads low;
- `total_area`, `c_bulk`, `c_min`, `c_max`, `gas_volume`;
- `bubble_dx/dy/dz`: gas centroid offset from the setpoint;
- `rise_u/v/w`: lab-frame velocity of the bubble over the interval;
- `u_gas_*`, `u_liq_*`: gas and liquid mean velocities in the frame;
- `F_pid_*`, `pid_int_*`, `frame_u/v/w`: controller state.

`sherwood.py --tmin 10 --tmax 30 <runs>` averages Sh (with the standard error
of 1 time unit batch means), Sh_sink, the rise velocity, the lateral speed
and the rise Reynolds number over a window. `path.py <runs>` integrates the
lab-frame lateral velocity into the bubble's sideways path, the counterpart
of gpu/mesh-sensitivity's drift.py: when it exceeds 0.01, 0.05 and 0.25 D,
and the growth rate of the lateral speed.
