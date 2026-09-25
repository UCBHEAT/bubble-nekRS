# pid-centering

Based on `sanitycheck/gpu` (benl/sanitycheck eb12acd): a quasi-2D rising
bubble (z-invariant cylinder of diameter 1, one periodic z element) in a fully
periodic box, reduced for fast iteration, with a PID controller that holds the
bubble at the domain center.

- Domain: 2x2 bubble diameters in x-y, centered on the bubble, 8x8 elements
  (the baseline element size h = 0.25) at polynomial order 7, and one 0.25
  cube element in z. The level set interface width uses the element volume,
  so the z element must stay a cube.
- dt = 0.0025, with TLSR/CLSR reinit every 100/10 steps.
- Same FLiBe/Ar nondimensional properties, level set settings
  (`interfaceWidthFactor = 1.5`), and `[SCALAR C] solver = none` as the
  baseline. Because the baseline's single z element is 1 thick, its interface
  width is eps = 1.5*(0.25*0.25*1)^(1/3)/9 = 0.066, against 0.054 here.

The 2x2 periodic box is a dense periodic array of bubbles (19.6% gas area
fraction, 1 D gap to the periodic images, and the bubble rises into its own
wake every 2 D of travel), and that confinement affects the rise velocity and
wake. The comparison with the baseline 4x8 box does not isolate it, since the
interface width differs too: on a 10x10 mesh (about the baseline in-plane
resolution, eps = 0.043) slip_v is +14% above the baseline at t = 0.5, +9% at
t = 1, crosses over near t = 1.7, and is -6% at t = 3 (the baseline data stop
there). The 8x8 mesh is another ~8% lower at t = 3. Centering does not remove
the confinement; it only keeps the bubble fixed in the mesh.

## PID centering force

The controller applies a spatially uniform acceleration (force per unit mass,
same units as the `1/Fr^2 = 0.375` gravity term in `buoyancySource`) to the
whole domain:

    F_pid = -(kp*e + ki*integral(e dt) + kd*u_gas)

- `e` is the gas-phase centroid offset from the domain center,
  `integral((1-psi) x dV)/integral((1-psi) dV) - x_0`.
- `u_gas` is the gas-phase mean velocity, used as de/dt (derivative on
  measurement, which avoids differencing the centroid across level set
  reinit jumps).
- Gains are in the `[PID]` section of `bubble.par`. The plant (`F_pid` to
  centroid) is a double integrator, so the loop is stable for `kd > 0`,
  `ki > 0` and `kd*kp > ki` (`ki = 0`: `kd > 0`, `kp > 0`); the setup prints
  a warning otherwise. The
  defaults `kp=12, ki=8, kd=6` put a triple closed-loop pole at `s = -2`. In
  the t = 0-10 run the offset peaks at 0.021 D near t = 1 while the bubble
  accelerates, then stays within 0.008 D while the rise velocity keeps
  changing (it steps up from ~0.45 to ~0.65 as the bubble rises into its
  periodic wake); x stays symmetric.

In the periodic box the gravity and pressure gradient terms exert zero net
force, so a uniform acceleration is exactly an acceleration of the reference
frame. It shifts the whole velocity field uniformly and does not change the
bubble/liquid relative motion: the frame rises with the bubble and the liquid
flows down past it. A constant rise velocity needs no steady force; the
integral term cancels slow drifts, mainly the small upward shift of the
centroid at each CLSR reinit that `u_gas` does not see (so at steady state
`ki*integral = kd*drift`, and `slip_v` is a few percent below the rise rate
of the level set centroid).

Implementation: `pid.hpp`. The controller is updated once per step in
`UDF_ExecuteStep` (after `lvlSet::solve`), and the stored force is added in
`customSource` after `applySurfaceTensionAcc` (which overwrites the explicit
term buffer). On restart the integral term is restored from the `data.csv`
row at the restart time.

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
current column layout, i.e. with `pid_int_*`, or the integral starts from 0).
nekRS numbers the new output from `bubble0.f00000` again and rewrites
`bubble.nek5000` to list only the new files, so move the existing `bubble0.f*`
aside first, or restart in a copy of the case directory. `animate.py` then
shows the fields of the restarted segment only.

## Output

- `data.csv`, one row per checkpoint, appended to by every run. In a single
  uninterrupted run, row k is `bubble0.f<k-1>` (same `Time`); after a restart,
  nekRS numbers the new field files from `bubble0.f00000` again and `tstep`
  restarts, so match rows to field files by `Time`. In addition to the
  baseline columns:
  - `bubble_dx, bubble_dy, bubble_dz`: gas centroid offset from the domain center (`e`).
  - `slip_u, slip_v, slip_w`: gas minus liquid mean velocity; `slip_v` is the rise velocity.
  - `F_pid_x, F_pid_y, F_pid_z`: PID force per unit mass applied in the next step.
  - `pid_int_x, pid_int_y, pid_int_z`: the integral term state.
- A `PID:` line in the log every step with `e`, `u_gas`, `F_pid`, and the gas
  volume, for full time resolution.

## Animation

```
~/.local/paraview-5.13.2/bin/pvbatch animate.py [CASE_DIR]
```

Builds the pipeline (3D render of the interface over a liquid z-vorticity
mid-plane with velocity glyphs, plus charts of `F_pid` and the centroid
offset from `data.csv` drawn up to the current time), writes `frames/`,
encodes `bubble-pid.mp4` with ffmpeg, and saves `animate.pvsm`. To edit in the
GUI, load `animate.pvsm` (File > Load State, "Search files under specified
directory" to point it at a case directory) or run `animate.py` from the
Python shell, which only builds the pipeline.

## Time step and performance

On the WSL laptop's GTX 1050, the run takes ~110-140 s of wall time per unit
of simulated time when it has the GPU to itself (~131 s on average, i.e. ~22
min of solve for t = 0-10, plus a minute or two of setup once the kernels are
cached). For comparison, `sanitycheck/gpu` (4x8 D, N = 9,
dt = 0.001) takes ~2200 s per unit on the same GPU, and this 2x2 case at
dt = 0.001 on a 10x10 mesh took ~315 s per unit.

- dt is limited by the explicit surface tension (capillary) constraint, not
  advection (CFL <= 0.225 over t = 0-10). Stable dt is ~0.8x the estimate
  `sqrt((rho_l+rho_g)*(h/N)^3/(4*pi*sigma))`; above it the CFL develops a
  step-to-step sawtooth, which can blow up the run:

  | mesh (h/N)          | dt     | t = 0-3 unless noted                                  |
  |---------------------|--------|-------------------------------------------------------|
  | 10x10, N=7 (0.0286) | 0.001  | stable (stopped at t = 1.75)                          |
  |                     | 0.0015 | stable, matches dt = 0.001 to t = 1.7                 |
  |                     | 0.002  | sawtooth from t = 1.5, blew up at t = 1.7             |
  | 8x8, N=7 (0.0357)   | 0.002  | stable                                                |
  |                     | 0.0025 | stable (used; also stable over t = 0-10)              |
  |                     | 0.003  | sawtooth from t = 0.76, CFL up to 0.93 (stopped at t = 1.95) |

- Level set reinitialization is ~80% of the run time: each CLSR call runs a
  fixed 10(N+1) = 80 pseudo-steps (~2 s) and each TLSR call 25(N+1) = 200
  (also ~2 s, as TLSR pseudo-steps are cheaper; a step that runs both takes
  ~4-4.5 s), while a normal step takes ~0.05 s. The intervals are kept at 10 and
  100 steps (the lvlSet defaults), so a larger dt means fewer calls per unit
  time.
- Reinit settings change the answer at the few percent level. At t = 3 on
  8x8, slip_v is 0.450 at dt = 0.002, 0.441 at dt = 0.0025, and 0.482 with
  `[CLSR] maximumSteps = 20` at dt = 0.002 (which runs ~1.8x faster); 10x10
  at dt = 0.0015 gives 0.478.
