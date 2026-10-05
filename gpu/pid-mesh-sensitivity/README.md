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
| 2.25x  | 4/12       | 2.21        | 7^3, half-width 1.17   | 13x23x13 = 3887   | 1.33M | 3456 (1.19M)   | 1.46 | 2.60 |
| 1.5x   | 4/18       | 1.47        | 9^3, half-width 1.0    | 17x31x17 = 8959   | 3.07M | 11664 (4.0M)   | 1.29 | 3.42 |
| 1x     | 4/27       | 0.98        | 13^3, half-width 0.963 | 23x45x23 = 23805  | 8.17M | 39366 (13.5M)  | 0.77 | 3.87 |
| 0.667x | 4/40       | 0.663       | 17^3, half-width 0.85  | 29x62x29 = 52142  | 17.9M | 128000 (43.9M) | 0.61 | 5.19 |
| 0.444x | 4/60       | 0.442       | 25^3, half-width 0.833 | 41x92x41 = 154652 | 53.0M | 432000 (148M)  | 0.41 | 5.42 |

![The z = 0 plane of the five meshes, and the 1x mesh around the bubble](doc/mesh.png)

The z = 0 plane of the five meshes and, on the right, the 1x mesh around the
bubble at t = 30 with its GLL points (`figures.py`).

The fine cube has the generator's cell count (the smallest even number of
cells that covers the interface zone r <= 0.76) plus one, so the bubble centre
lies inside an element. With an even count the centre is an element vertex
and all three axes through it are element edges; there the gas core carries
strong spurious currents at the physical density ratio (u_max up to 7 on the
1.5x mesh at t = 0.5, against 1.4 with the odd count), which then corrupt the
level set along the edges and seed spurious gas in the liquid, and the 2.25x
bubble breaks up at t = 1.7. (In gpu/mesh-sensitivity the bubble rose through
the mesh, so it only passed such points; here the PID holds it in place.)

The interface width is eps = 1.5 h/N ([LVLSET] interfaceWidthValue), which is
what interfaceWidthFactor = 1.5 gives on the uniform meshes; the lvlSet
default would take h from the largest element.

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

Besides the odd fine cube (see Meshes), three changes to the PID case were
needed at the physical density ratio (it had only run at rho ratio 40):
- **Fixed outflow.** The outflow imposes the frame velocity, like the
  inflow, so both fluxes balance exactly and the pressure has no Dirichlet
  boundary. The pressure outflow of `gpu/pid-centering-3D` (zeroNeumann
  velocity, p = 0 with Dong's backflow term) blew up at an outflow node at
  t = 0.073 on the 1x and 0.667x meshes alike, the velocity there growing
  ~500x per step, with either pressure extrapolation order of the rho
  splitting. The wake is still 10 D long before it meets the outflow.
- **No noise in the liquid (psi snap).** After every level set step psi is
  set to 1 wherever 1 - psi < 1e-6 ([CASEDATA] psiSnap). CLSR takes its
  normals from the TLS, which TLSR makes a distance function only within
  ~5 finest elements of the interface (200 pseudo-steps at CFL 0.4, the
  same reach as the default on gpu/mesh-sensitivity's uniform meshes).
  Farther out the normals are arbitrary, and CLSR amplifies the ~1e-7
  deficits the CLS solve leaves in the liquid wherever they diverge, until
  spurious gas pockets form at element edges and vertices and their
  buoyancy blows the run up. Without the snap this happened 1-2 D from the
  bubble on the 0.667x and 0.444x meshes (0.667x: deficit 5e-7 at t = 0.1,
  5e-3 at t = 0.9, blow-up at t = 1.26), and, before the far-field reset
  below, anywhere in the coarse far field (t = 0.35 on 1x, 0.78 on 2.25x;
  with the default reinit caps, which redistance farther, already from
  t = 0.035 near the inflow, where the zero-Neumann TLSR condition bends the
  distance function). The uniform periodic meshes of gpu/mesh-sensitivity
  are fine everywhere and have no boundaries. The real interface tail only
  falls below 1e-6 about 14 eps (3 finest elements) from the interface, so
  the snap removes noise and the outermost tail, which CLSR partly restores
  at every call: data.csv records the volume removed (`snap_removed`,
  ~1.6e-5 per 0.1 time units on 0.667x, i.e. ~1% of the bubble by t = 30,
  in proportion to eps).
- **Pure far-field liquid.** psi is also reset to 1 farther than 2 D from
  the PID setpoint after every level set step, which catches anything
  larger that gets there (the bubble stays within ~0.04 D of the setpoint
  and its interface band within ~1 D). data.csv records the largest deficit
  found there before the reset (`far_deficit_max`, ~1e-7) and the gas volume
  removed (`far_removed`).

`[CASEDATA]` switches for tests: `traceVelocity` and `traceFarField` print
where the largest velocity or the far-field deficit is while they exceed the
given value; `psiSnap = 0`, `farFieldClean`, `scalarSVV`,
`rhoSplittingFilter` and `pressureExtOrder` turn the corresponding features
off or on.

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
- `far_deficit_max`, `far_removed`, `snap_removed`: the far-field reset and
  the psi snap (see Numerics), and `u_max`, the largest velocity;
- `bubble_dx/dy/dz`: gas centroid offset from the setpoint;
- `rise_u/v/w`: lab-frame velocity of the bubble over the interval;
- `u_gas_*`, `u_liq_*`: gas and liquid mean velocities in the frame;
- `F_pid_*`, `pid_int_*`, `frame_u/v/w`: controller state.

`sherwood.py --tmin 10 --tmax 30 <runs>` averages Sh (with the standard error
of 1 time unit batch means), Sh_sink, the rise velocity, the lateral speed
and the rise Reynolds number over a window, and reports when Sh settled
(every later 1 time unit batch mean within 0.5% of the window mean) and the
gas volume lost. `path.py <runs>` integrates the lab-frame lateral velocity
into the bubble's sideways path, the counterpart of gpu/mesh-sensitivity's
drift.py: when it exceeds 0.01, 0.05 and 0.25 D, and the growth rate of the
lateral speed. `wake.py <checkpoint> [<first checkpoint of its job>]` reads a
checkpoint directly (no ParaView) and prints the bubble's aspect ratio and
the liquid velocity on the axis behind it, which shows whether a standing
eddy has formed. `figures.py` draws the two figures above (doc/), and
`submit-animate.sh` with `animate.py` (gpu/mesh-sensitivity's, adapted)
makes each run's animation.mp4 and post.csv, the bubble's volume, centroid,
velocities, interface area and extents at every checkpoint, on a Polaris GPU
node.

## Sc = 1 results (Polaris, October 2026)

### Runs

All five runs went from t = 0 to 30 in eight bundled jobs of 20-24 nodes on
the small queue (October 1-4, 2026), 450 node-hours in all, 80% of them for
0.444x. Measured seconds per step, and what a run to t = 30 takes at that
size:

| Run    | nodes | s/step | steps     | hours | node-hours |
|--------|-------|--------|-----------|-------|------------|
| 2.25x  | 1     | 0.048  | 300,000   | 4.0   | 4          |
| 1.5x   | 1     | 0.074  | 300,000   | 6.2   | 6          |
| 1x     | 2     | 0.096  | 300,000   | 8.0   | 16         |
| 0.667x | 4     | 0.070  | 600,000   | 11.6  | 46         |
| 0.444x | 16    | 0.0625 | 1,200,000 | 20.8  | 333        |

0.444x ran its last four jobs on 20 nodes at 0.061 s per step, only 2.5%
faster: 16 nodes is the economical size (4.2-4.3 time units per 3 hour job
on either). After the animations (`submit-animate.sh`), the checkpoints were
deleted except two self-contained ones per run, at t = 10 and 30
(`bubble_t10.fld`, `bubble_t30.fld`, which can also seed runs on other
meshes with `startFrom = <file>+int`). answinter26/Sc1-pid/<run>/ keeps the
inputs, the mesh, data.csv, post.csv, animation.mp4, logs/ (the jobs' logs,
gzipped) and those two checkpoints, 10 GB in all.

### Sherwood number

Sh (± the standard error of 1 time unit batch means) and Sh_sink, with
`sherwood.py`; settled is when Sh reached its plateau (every later batch
mean within 0.5% of the t = 10-30 mean, - if it still drifts at t = 30):

| Run    | t = 5-10     | t = 10-20    | t = 20-30    | t = 10-30    | Sh_sink, t = 10-30 | settled |
|--------|--------------|--------------|--------------|--------------|--------------------|---------|
| 2.25x  | 10.15 ± 0.09 | 9.65 ± 0.07  | 9.57 ± 0.03  | 9.61 ± 0.04  | 6.63               | -       |
| 1.5x   | 14.26 ± 0.01 | 14.36 ± 0.03 | 14.68 ± 0.01 | 14.52 ± 0.04 | 10.28              | -       |
| 1x     | 15.91 ± 0.01 | 15.92 ± 0.01 | 15.97 ± 0.01 | 15.94 ± 0.01 | 11.17              | 5       |
| 0.667x | 16.88 ± 0.04 | 16.94 ± 0.00 | 16.98 ± 0.01 | 16.96 ± 0.01 | 11.59              | 7       |
| 0.444x | 17.15 ± 0.02 | 17.20 ± 0.00 | 17.20 ± 0.00 | 17.20 ± 0.00 | 11.61              | 6       |

Sh dips to 11-12 at t = 1, while the bubble accelerates from rest, and is
steady from t = 5-7 on every resolved mesh. Afterwards it only creeps up with
the gas loss (below): from t = 10-15 to 25-30 by 2.7% on 1.5x, 0.5% on 1x and
0.3% on 0.667x and 0.1% on 0.444x. Unlike in gpu/mesh-sensitivity,
nothing changes later on: the bubble rises straight on every mesh (lateral
speed at most 1.1e-4, `path.py`).

Over t = 10-30, the three finest meshes converge at an observed order of 3.7,
with GCI = 1.4% (p = 2, F_s = 1.25) between 0.667x and 0.444x, and Richardson
extrapolation (p = 2) gives Sh = 17.39. With 1.5x, the observed order of the
coarser triplet is 0.8, but it depends on the window (1.0 over t = 10-15),
because the 1.5x Sh drifts with its gas loss. The figures and the GCI table
are made by 2026-11-nek-cst/data/mesh_sensitivity_sc1_pid.py in
UCBHEAT/papers.

At the Reynolds number of the simulated rise (Re = 263, from the 0.444x rise
speed and d_eq), the potential flow (Boussinesq) value is Sh = 18.3 and Feng
and Michaelides (2001) give 19.0, so the extrapolated Sh is 5% and 9% below
them; the rigid-sphere Frossling correlation gives 11.0.

**Comparison with gpu/mesh-sensitivity.** Before its bubbles drifted
(t = 5-10), gpu/mesh-sensitivity's Sh, which came from the hard-sink count,
was 10.01, 11.02, 11.41 and 11.47 on 1.5x to 0.444x; Sh_sink here is 10.06,
11.14, 11.54 and 11.58 over the same window, 0.5-1.1% higher. The budget Sh is
1.41 (1.5x) to 1.48 (0.444x) times Sh_sink, approaching 1/g0 = 1.5, the BDF2
undercount of the sink count (see `../pid-centering/README.md`). So the two
setups agree where they can be compared, and gpu/mesh-sensitivity's Sh (and
the figure and GCI table of mesh_sensitivity_sc1.py) read low by that factor.

### Rise, shape and wake

Averages over t = 10-30, with d_eq from the gas volume and Re = 231.4 x rise
speed x d_eq; u_max is the largest velocity in the domain after t = 1, and
the aspect ratio (equatorial radius over half height, where psi = 0.5,
`wake.py`) is at t = 30:

| Run    | rise velocity | d_eq  | Re  | area  | u_max / rise | aspect ratio |
|--------|---------------|-------|-----|-------|--------------|--------------|
| 2.25x  | 0.842         | 1.033 | 201 | 3.171 | 4.33         | 1.20         |
| 1.5x   | 1.019         | 1.013 | 239 | 3.274 | 2.81         | 1.65         |
| 1x     | 1.066         | 1.004 | 248 | 3.325 | 2.19         | 1.74         |
| 0.667x | 1.123         | 1.001 | 260 | 3.397 | 1.78         | 1.93         |
| 0.444x | 1.139         | 0.999 | 263 | 3.389 | 1.71         | 1.88         |

The rise velocity settles more slowly than Sh: on 0.444x it is 1.121 at t = 8,
1.131 at 10, 1.139 at 14 and 1.141 from t = 18 on. d_eq > 1 on the coarse
meshes because the gas volume, the integral of 1 - psi, exceeds the sharp
volume by (4/3) pi^3 R eps^2 for the CLS profile (0.105, 20%, on 2.25x at
t = 0, and 0.4% on 0.444x).

No standing eddy forms behind the bubble on the resolved meshes: on 1.5x to
0.444x the liquid on the axis below it moves away from it everywhere (on
2.25x, whose interface is 0.2 D thick, the gas circulation carries liquid
up to 0.27 D below the rear). The liquid right behind the rear moves more
slowly the finer the mesh, though: at less than 0.1 times the inflow speed
up to 0.13, 0.28 and 0.38 D below the rear on 1x, 0.667x and 0.444x, so the
converged flow may be close to separating.

![c on the z = 0 plane of the 0.444x run at t = 30](doc/wake.png)

c on the z = 0 plane of the 0.444x run at t = 30 (`figures.py`): the thin
boundary layer over the front of the bubble, and the depleted liquid that
leaves its rear in a wake that still has c = 0.42 on the axis 2.5 D behind
it and recovers slowly towards the outflow.

### Gas volume

| Run    | gas volume change, t = 0.1-30 | removed by the psi snap | by the far-field reset |
|--------|-------------------------------|-------------------------|------------------------|
| 2.25x  | -12.3%                        | 9.4%                    | 1.8e-5                 |
| 1.5x   | -6.5%                         | 4.3%                    | 1.3e-5                 |
| 1x     | -3.6%                         | 2.0%                    | 1.2e-5                 |
| 0.667x | -2.3%                         | 1.0%                    | 2.2e-5                 |
| 0.444x | -1.5%                         | 0.5%                    | 1.6e-5                 |

(the far-field column in volume, the others in percent of the gas volume at
t = 0.1). The snap's share halves with each refinement, about as h^2; the
rest, 2.9% on 2.25x down to 0.9% on 0.444x, is the conservative level set's
own drift, and the far-field reset removes nothing that matters.

### Run length

A run to t = 15, averaged over t = 10-15, would have given the same result
for half the cost: Sh over t = 10-15 is within 0.3% of the t = 10-30 mean on
1x and finer (0.03% on 0.444x), and the convergence study barely moves (GCI
1.51% instead of 1.38% between 0.667x and 0.444x, observed order 3.6 instead
of 3.7, extrapolated Sh 17.40 instead of 17.39). The rise velocity over
t = 10-15 is 0.35% lower than over t = 10-30 on 0.444x. Later studies (higher
Sc) can likely stop at t = 15: c is passive, so the flow and its start-up do
not depend on Sc, and the concentration boundary layer is renewed by the
flow past the bubble, on the time it takes to pass it, not by diffusion;
`settled` from `sherwood.py`, and `wake.py` for the nearly stagnant liquid
behind the bubble, tell whether that still holds.
