# pid-centering

Based on `sanitycheck/gpu` (benl/sanitycheck eb12acd): a quasi-2D rising
bubble (z-invariant cylinder of diameter D = 1, one periodic z element) held in
place by a PID controller that accelerates the reference frame, so the liquid
flows down past a stationary bubble. The liquid enters at the top, at rest in
the lab frame and saturated with species (c = 1), and leaves at the bottom.

- Domain: x in [-3, 3] (periodic, 6 D wide), y from -6 (outflow, 6 D below
  the bubble centre) to 3 (inflow), bubble at the origin. 16x19 elements at
  polynomial order 7, graded from an 8x8 block of h = 0.25 cubes around the
  bubble (the baseline element size), and one 0.25 cube element in z.
- dt = 0.0025, with TLSR/CLSR reinit every 100/10 steps.
- Species: c is a passive scalar (CST off, liquid diffusivity 1/Pe
  everywhere) with a hard sink (c = 0 wherever psi < 0.5); no soft source.
- Same FLiBe/Ar nondimensional properties as the baseline.

## PID centering force

The controller applies a spatially uniform acceleration (force per unit mass,
same units as the `1/Fr^2 = 0.375` gravity term in `buoyancySource`) to the
whole domain, i.e. the fictitious force of an accelerating frame:

    F_pid = -(kp*e + ki*integral(e dt) + kd*u_gas)

- `e` is the gas-phase centroid offset from the setpoint `[PID] x0, y0` (the
  initial bubble position), `integral((1-psi) x dV)/integral((1-psi) dV) - x_0`.
- `u_gas` is the gas-phase mean velocity, used as de/dt (derivative on
  measurement, which avoids differencing the centroid across level set
  reinit jumps). Each CLSR reinit also shifts the centroid up by ~3e-4,
  which `u_gas` does not see; the integral term cancels that drift.
- Gains are in the `[PID]` section of `bubble.par`. The plant (`F_pid` to
  centroid) is a double integrator, so the loop is stable for `kd > 0`,
  `ki > 0` and `kd*kp > ki` (`ki = 0`: `kd > 0`, `kp > 0`); the setup prints
  a warning otherwise. The defaults `kp=12, ki=8, kd=6` make the closed-loop
  characteristic polynomial `s^3 + kd s^2 + kp s + ki` equal to `(s + 2)^3`
  (s is the Laplace variable: solutions go like `e^(s t)`). So all three
  closed-loop poles are at `s = -2`, and an offset decays like `e^(-2t)`, with
  a time constant of 0.5, without oscillating. For poles at `s = -a`, use
  `kp = 3a^2, ki = a^3, kd = 3a`. The offset peaks at ~0.02 D near t = 1
  while the bubble accelerates.

The liquid far from the bubble is at rest in the lab frame, so in the
accelerated frame it moves with the frame velocity
`U_frame = integral(F_pid dt)`, which is imposed as the inflow velocity
(passed to `udfDirichlet` through `bc->o_usrwrk`). Together, the uniform acceleration and the matching
inflow are a change of reference frame: they do not change the bubble/liquid
relative motion, up to O(dt) differences between the forward-Euler inflow
update and the BDF/EXT integration of the forcing, which a uniform pressure
gradient absorbs (~0.1% of gravity on average here). The buoyancy reference
density is the liquid's, since the far-field pressure gradient is the
liquid's hydrostatic one.

Implementation: `../common/pid.hpp`. The controller is updated once per
step in `UDF_ExecuteStep` (after `lvlSet::solve`); `customSource` adds the
stored force after `applySurfaceTensionAcc` (which overwrites the explicit
term buffer) and sets that step's inflow velocity. On restart, the integral term
and `U_frame` are restored from the `data.csv` row at the restart time.

## Boundary conditions

| field | inflow (top, ID 1) | outflow (bottom, ID 2) |
|---|---|---|
| velocity | `udfDirichlet`: `U_frame` | `zeroNeumann`, with Dong's stabilized outflow pressure |
| tls | `zeroNeumann` | `zeroNeumann` |
| cls (and CLS reinit) | `udfDirichlet`: psi = 1 | `zeroNeumann` |
| c | `udfDirichlet`: c = 1 | `zeroNeumann` |

x and z are periodic. The boundary IDs come from the `1`/`2` placeholders in
`bubble.box`, mapped by `usrdat2` in `bubble.usr`.

## Level set and species settings

- `[LVLSET] interfaceWidthValue = 0.0536` (1.5 h/N on the fine block). By
  default lvlSet takes eps from the largest element, which on this graded
  mesh would give 0.080 (0.12 with interfaceWidthFactor = 1.5). Keep the
  bubble radius large compared with eps: at eps = 0.085 (R/eps = 5.9) the
  gas core goes unstable, with spurious velocities ~1 in a bubble rising at
  0.26, while R/eps = 9.3 is clean.
- `[TLSR]`/`[CLSR] maximumSteps = 200/80`: the pseudo-step counts scale with
  the largest/smallest element size, so they are capped at the uniform-mesh
  values.
- `customProperties` calls `scalar->mueSVV()`: nekRS-LS only computes the
  scalar SVV viscosity when no `userProperties` hook is set, so without it the
  `svv` regularization in `[SCALAR *]` does nothing, and c (Pe = 2.7e5)
  overshoots to c ~ 1.5 within t = 1.
- c starts as the smeared c = psi, so the liquid side of the interface band
  starts partly depleted, and early on the depleted wake (the c = 0.9 isoline
  in the animation) is mostly that band swept to the rear and trapped in the
  wake eddies. A saturated start (c = 1 wherever psi >= 0.5) is worse
  numerically: the hard sink then cuts a full 1 -> 0 jump at psi = 0.5, which
  the mesh cannot resolve at Pe = 2.7e5, and c rings to ~1.7 by t = 0.5.

## Running

```
nekls               # NEKRS_HOME=~/.local/nekrs-ls
nek5000             # for genbox
./mesh
rm -f data.csv      # data.csv is appended to by every run
nrsmpi bubble 1
```

To restart, copy a checkpoint (e.g. `bubble0.f00040`) to `restart.fld`,
uncomment `startFrom` in `bubble.par`, and keep `data.csv` (written with the
current column layout). If `data.csv` has no row at the restart time, the
integral starts from 0 and `U_frame` is taken from the restart field's
inflow velocity. nekRS numbers the new output from `bubble0.f00000` again
and rewrites `bubble.nek5000` to list only the new files, so move the
existing `bubble0.f*` aside first, or restart in a copy of the case
directory. The animation then shows the fields of the restarted segment only.

## Output

- `data.csv`, one row per checkpoint, appended to by every run. In a single
  uninterrupted run, row k is `bubble0.f<k-1>` (same `Time`); after a restart,
  match rows to field files by `Time`.
  - `mdot`: species removed per unit time since the previous row, from the
    budget
    `integral(net boundary inflow of c) dt - change in integral(c dV)`.
    It uses only the stored c and the boundary fluxes, not the time
    integrator. With the sink disabled it reads |mdot| <= 1.4e-5 per 0.05
    window (<= 0.4% of the signal, ~0.1% on average) with both BDF1 and
    BDF2 (t = 0-1.5 tests). Adding up the c the hard sink zeroes instead
    reads 1.34-1.41x low with BDF2: zeroing a gas node also zeroes its
    history, so it only holds 1/g0 = 2/3 of what flowed in during the step.
    With BDF1 (g0 = 1) that count agrees with the budget to 0.3%. The same
    undercount affects any hard-sink bookkeeping that sums the zeroed c
    under BDF2, including `gpu/match-nek5000-202605`.
  - `MTC = mdot/(area * (c_inflow - 0))`, `Sh = MTC Pe`, with the inflow
    concentration c_inflow = 1 and c = 0 at the interface; `total_area` is
    the integral of the interface delta function.
  - `c_bulk, c_min, c_max`: liquid-volume average and range of c.
  - `bubble_dx, bubble_dy, bubble_dz`: gas centroid offset from the setpoint (`e`).
  - `rise_u, rise_v, rise_w`: lab-frame velocity of the (level set) bubble
    centroid, i.e. its drift in the frame since the previous row minus
    `U_frame`; `rise_v` is the rise velocity. (The gas material velocity
    `u_gas - U_frame` reads ~2% lower, since it misses the reinit shifts.)
  - `F_pid_x, F_pid_y, F_pid_z`: PID force per unit mass applied in the next
    step (the current value, which swings by ~+-0.015 with the reinit cycle;
    the mean over a window is the change of `frame_*` divided by its length).
  - `pid_int_*`, `frame_u, frame_v, frame_w`: controller state (integral, `U_frame`).
- A `PID:` line in the log every step with `e`, `u_gas`, `F_pid`, the gas
  volume, and `U_frame`, for full time resolution.

## Animation

```
~/.local/paraview-5.13.2/bin/pvbatch ../common/animate.py [CASE_DIR]
```

`gpu/common/animate.py` is shared with `gpu/pid-centering-3D`; CASE_DIR
defaults to the current directory. It builds the pipeline:
- a 3D render of the interface over a liquid z-vorticity mid-plane, with
  velocity glyphs and the c = 0.9 isoline of species-depleted liquid (drawn
  only where psi > 0.95, because c rings inside the interface band next to
  the hard-sink cut);
- charts of `F_pid`, as its mean over the last 0.5 time units from `frame_*`,
  and of the centroid offset from `data.csv`, drawn up to the current time.

It writes `frames/`, encodes `bubble-pid.mp4` with ffmpeg, and saves this
case's `animate.pvsm`. To edit in the GUI, load `animate.pvsm` (File > Load
State, "Search files under specified directory" to point it at a case
directory with the same domain; the camera and glyph grid keep the saved
domain bounds, so for a different domain rerun `pvbatch ../common/animate.py
CASE_DIR`), or run `animate.py` from the Python shell, which only builds the
pipeline.

## Domain size

**Outlet.** The outlet 6 D below the bubble centre is far enough for these
results. A run with the outlet at 12 D (y from -12 to 3, otherwise identical)
agrees to 0.1-0.7% through t = 20, and neither outlet sees backflow. (These
two runs predate the budget diagnostics, so the table uses `U_frame` and the
hard-sink count, which both files have.)

| t = 0-20 | outlet at -6 | outlet at -12 |
|---|---|---|
| rise velocity (-`U_frame`) at t = 20 | 0.6584 | 0.6576 |
| Sh from the hard-sink count, t in 8-12 / 12-16 / 16-20 | 1371 / 1332 / 1304 | 1371 / 1333 / 1313 |
| centreline lab-frame liquid velocity / rise velocity at y = -1, -2, -3 (t = 20) | 1.333, 1.616, 1.382 | 1.335, 1.618, 1.384 |
| same at y = -5 (0.5 D above the short outlet) | 0.392 | 0.400 |

The wake is long, though. A pair of standing eddies about 2.5 D long forms
behind the (slightly oblate) bubble, where the liquid on the centreline
rises up to 1.6x faster than the bubble, and at t = 20 the wake still carries
0.2-0.4 of the rise velocity through y = -6 and is still lengthening. For much
longer runs, check the outlet again, or use the 12 D outlet: 400 instead of
304 elements, but it ran only ~4% slower (418 vs 404 s per unit when the two
shared the GPU), since the GPU is latency-bound at this size.

**Width and inflow distance.** Both matter for the rise velocity. Moving
the inflow from 3 D to 5 D above the bubble centre (W = 6) makes the bubble
rise 2.4% faster, and widening the periodic domain from 6 D to 10 D (inflow
at 5 D) adds another 5.2%:

| mean over t in | W = 6, inflow at 3 (committed) | W = 6, inflow at 5 | W = 10, inflow at 5 |
|---|---|---|---|
| rise velocity, 8-12 / 12-16 / 16-20 | 0.634 / 0.648 / 0.656 | 0.647 / 0.663 / 0.672 | 0.673 / 0.694 / 0.707 |
| budget Sh, 8-12 / 12-16 / 16-20 | 1837 / 1847 / 1834 | 1843 / 1846 / 1845 | 1872 / 1857 / 1822 |
| lateral return flow beside the bubble / V (t = 20) | 0.20-0.28 | 0.21-0.29 | 0.12-0.16 |

The return flow is the lab-frame downflow at |x| >= 1.5, y = 0 to -2: the
upward flux the wake carries, coming back down beside the bubble. It falls
roughly as 1/W, so W = 10 is still confined. Extrapolating the two widths
with a deficit that scales as 1/W or 1/W^2 puts the unconfined rise velocity
3-7% above the W = 10 value (two widths cannot fix the exponent). A longer
inflow section (e.g. 7 D) would show whether 5 D is enough. Sh differs by
at most 2% between the three runs while the rise velocity differs by 7.8%
(see "What Sh means here").

For rise velocity or wake studies, use at least the W = 10, inflow-at-5 mesh
(20x21 elements; it ran ~4% slower than the committed mesh when the two
shared the GPU) and check a wider one. The `bubble.box` node lines are

```
20 21 -1               Nelx Nely Nelz (positive: explicit node coordinates follow)
-5 -4 -3 -2.3 -1.75 -1.325 -1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1 1.325 1.75 2.3 3 4 5
-6 -5 -4 -3.05 -2.3 -1.75 -1.325 -1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1 1.325 1.75 2.3 3 4 5
```

The committed 6 D domain is kept as the fast testbed for the PID and species
numerics; the animation script fits its camera and glyph grid to whichever
domain it is run on.

The 2x2 D periodic box used earlier was far too small: a dense bubble
array, and with an outlet 0.55 D behind the bubble its wake would give
backflow over 25-29% of the outlet.

In the committed 6 D domain, by t = 20 the rise velocity is 0.66 and still
creeping up (0.648 and 0.656 over t = 12-16 and 16-20), the budget Sh has
settled to 1834-1847 over t = 8-20 (2048 over t = 4-8), and the centroid
stays within 2.4e-4 D of the setpoint for t > 12.

## What Sh means here

The `Sh` in `data.csv` is a well-defined budget of the species the sink
removes, but it is set by the numerics, not by resolved interfacial
transport, so do not read it as a physical transfer rate:

- For a clean (mobile) interface, penetration theory with the
  potential-flow surface speed 2V sin(theta) around a cylinder gives
  `Sh = (8/pi) sqrt(Pe_V/(2 pi)) = 1.02 sqrt(Pe_V)`. With the rise velocity
  V = 0.66, Pe_V = V Pe = 1.75e5 and Sh ~ 425, an upper estimate for a clean
  bubble (an oblate shape adds ~1.5%). The run gives Sh = 1834 (budget)
  for t in 16-20, 4.3x that.
- The physical concentration boundary layer, D/Sh ~ 2.4e-3, is 7-22x
  thinner than the GLL spacing on the fine block (0.016-0.052) and 23x
  thinner than eps = 0.0536. What sets the flux is numerical transport across
  the interface band: the SVV diffusion, the sink cut at psi = 0.5, and the
  CLSR reinit, which shifts the bubble centroid up by ~3e-4 per call
  (0.012-0.017 D per unit time), so the front of the psi = 0.5 contour
  advances into the liquid. If the liquid swept that way were at c = 1, that
  sweep alone would amount to Sh ~ 1250.
- It depends on the regularization: over t = 0-1 on the periodic box,
  turning on the scalar SVV raised the species loss by 48% and the counted Sh
  by 68%. The SVV viscosity scales with the local |u|, which here includes
  the frame velocity, so the regularization is not frame-invariant either.
- c is not bounded next to the sink cut: it rings to [-0.19, 1.23] for t in
  16-20, at the rear of the bubble where psi = 0.55-0.76, and the part of
  the liquid (psi >= 0.5) with c > 1.01 grows from 1% at t = 5 to 4% at
  t = 20. With the sink disabled, c stays in [0, 1].
- It hardly responds to the flow: across the three domains in "Domain
  size" the rise velocity varies by 7.8%, which would change a resolved Sh
  by ~4% (Sh ~ sqrt(V)), but the budget Sh varies by at most 2% (1.3% for
  t in 16-20, where the fastest bubble has the lowest).

A physically meaningful Sh needs the concentration boundary layer resolved
by several GLL points, e.g. Pe lowered to where D/Sh is comparable to the
fine-block spacing, and a mesh/eps convergence check.

## Time step and performance

On the WSL laptop's GTX 1050, this 304-element case runs at ~180-220 s of
wall time per unit of simulated time when it has the GPU to itself, so
t = 0-20 takes about 1-1.2 h (the 336-element inflow-at-5 variant took 75 min
alone). Runs sharing the GPU are ~1.8-2x slower each: the t = 0-20 runs above
took 404-426 s per unit in pairs.

- dt is limited by the explicit surface tension (capillary) constraint on the
  fine block around the bubble, not by advection. Stable dt is ~0.8x the
  estimate `sqrt((rho_l+rho_g)*(h/N)^3/(4*pi*sigma))`; above it the CFL
  develops a step-to-step sawtooth, which can blow up the run. From a sweep
  on a periodic 2x2 D box with the same fine-block elements:

  | mesh (h/N)          | dt     | t = 0-3 unless noted                                  |
  |---------------------|--------|-------------------------------------------------------|
  | 10x10, N=7 (0.0286) | 0.001  | stable (stopped at t = 1.75)                          |
  |                     | 0.0015 | stable, matches dt = 0.001 to t = 1.7                 |
  |                     | 0.002  | sawtooth from t = 1.5, blew up at t = 1.7             |
  | 8x8, N=7 (0.0357)   | 0.002  | stable                                                |
  |                     | 0.0025 | stable (used)                                         |
  |                     | 0.003  | sawtooth from t = 0.76, CFL up to 0.93 (stopped at t = 1.95) |

- Level set reinitialization dominates the run time: each CLSR call runs 80
  pseudo-steps and each TLSR call 200, while a normal step is cheap. The
  intervals are kept at 10 and 100 steps (the lvlSet defaults), so a larger
  dt means fewer calls per unit time.
- Reinit settings change the answer at the few percent level: on the
  periodic 2x2 box at t = 3, slip was 0.450 at dt = 0.002, 0.441 at
  dt = 0.0025, and 0.482 with `[CLSR] maximumSteps = 20` (~1.8x faster).
