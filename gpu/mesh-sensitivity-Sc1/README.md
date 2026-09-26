# mesh-sensitivity-Sc1

Mesh sensitivity study of the Sherwood number for the gpu/3D bubble case at
Sc = 1 with the physical FLiBe/Ar density ratio, at 2.25x, 1.5x and 1x the
Kolmogorov scale, followed by 0.667x, 0.444x and 0.296x ((1.5)^-1 to
(1.5)^-3).

Physics follows gpu/3D (4x8x4 fully periodic domain, bubble of diameter 1
starting at y = 1, Re = 231.4, Fr = 1.633, We = 2.667, muratio = 148.1) with:
* Sc = 1, so Pe = Re = 231.4 (flibe1/ar in cpu/common/calc_dimensionless_numbers.py).
* Physical density ratio rhoratio = 1/0.0002709 = 3691 instead of 40.

Species transport and the Sherwood number follow gpu/match-nek5000-202605
rather than gpu/3D's continuous species transfer (CST) flux integral:
* Uniform diffusivity D = 1/Pe in both phases, no CST terms.
* Hard sink: every step, c is removed where psi < 0.5 and the removed amount
  is accumulated; mdot is its average over each sample interval.
* Soft source: c is driven towards 1 in the liquid bulk (psi > 0.9).
* MTC = mdot/(total_area*c_bulk) and Sh = MTC*Pe, written to data.csv every
  [CASEDATA] sampleInterval (0.1 time units; the checkpoint interval if
  unset), together with the gas volume, the bubble centroid (a circular mean,
  as the domain is periodic), the gas and liquid mean velocities and the rise
  velocity (their difference in y). The 2.25x, 1.5x and 1x runs predate the
  centroid and velocity columns; their post.csv has the same quantities from
  the checkpoints (animate.py).

Numerics follow gpu/3D except:
* Polynomial order 7 and TLSR/CLSR at the nekRS-LS defaults (no targetCFL,
  maximumSteps or SVV), as recommended by the level set developer; TLSR runs
  every 0.01 and CLSR every 0.001 time units (every 100 and 10 steps at
  dt = 1e-4).
* Each mesh is generated directly by genbox (no hrefine, which leaves the
  partition imbalanced when the base mesh has few elements per rank).
* The 2.25x, 1.5x and 1x meshes all use dt = 1e-4 so only the spatial
  resolution changes between them; each finer mesh halves dt. With the
  physical density ratio, larger dt gives spurious velocity spikes in the gas
  core (on element edges along the bubble axis) near t = 0.32 on the coarse
  meshes: the 2.25x run blew up at t = 0.72 with dt = 4e-4, and with
  dt = 2e-4 it survived the spike but its Sh differed from dt = 1e-4 by 10% at
  t = 0.4. On the 1x mesh, dt = 2e-4 and 1e-4 agree within 1%.

## Meshes

The Kolmogorov scale is lambda_k = 0.0671 mm (lambda_k/d = 0.021536). Mesh
resolution is measured by the mean unique GLL spacing dx = h/N at N = 7, with
element size h = 4/Nelx.

| Run    | Elements  | dx/lambda_k | GLL points (E*N^3) | dt      | checkpoints |
|--------|-----------|-------------|--------------------|---------|-------------|
| 2.25x  | 12x24x12  | 2.21        | 1.19M              | 1e-4    | 0.1         |
| 1.5x   | 18x36x18  | 1.47        | 4.00M              | 1e-4    | 0.1         |
| 1x     | 27x54x27  | 0.98        | 13.5M              | 1e-4    | 0.1         |
| 0.667x | 40x80x40  | 0.663       | 43.9M              | 5e-5    | 0.5         |
| 0.444x | 60x120x60 | 0.442       | 148M               | 2.5e-5  | 0.5         |
| 0.296x | 90x180x90 | 0.295       | 500M               | 1.25e-5 | 0.25        |

## Running on Polaris

Create the run directories and meshes (genbox must be on PATH; the finer
meshes need genbox built with MAXNEL >= 2*Nelx^3, e.g. answinter26/tools/
genbox-maxnel1.5M):

```
./setup-runs.sh ~/answinter26/Sc1
GENBOX=~/answinter26/tools/genbox-maxnel1.5M ./setup-runs.sh ~/answinter26/Sc1 \
    0.667x:40:5e-5:0.5 0.444x:60:2.5e-5:0.5 0.296x:90:1.25e-5:0.25
```

The first three meshes are small compared to the 10-node minimum of the prod
queue, so they ran concurrently in one job, each on its own nodes:

```
cd ~/answinter26/Sc1
PROJ_ID=nek-vf QUEUE=prod <this dir>/submit-bundle.sh 03:00 1x:7 1.5x:2 2.25x:1
```

submit-bundle.sh exports BUBBLE_STOP_AT, 5 minutes before the walltime runs
out; bubble3d.udf then ends the run after the next step, which writes a
checkpoint. If a job ends before endTime = 30, prepare the unfinished runs for
a restart and resubmit them:

```
<this dir>/prepare-restart.sh 1x
PROJ_ID=nek-vf QUEUE=prod <this dir>/submit-bundle.sh 03:00 1x:10
```

prepare-restart.sh moves the finished part's field files and logs to part<N>/
(nekRS restarts checkpoint numbering at 0), trims data.csv back to the
restart time and sets startFrom. After the last part, finish-runs.sh links
all parts' checkpoints in time order for ParaView.

A run can also be queued as a chain of jobs that restart themselves: with
PREPARE_RESTART=1 each job first runs prepare-restart.sh on its cases (and
skips those that have reached endTime), and DEPEND=<jobid> holds a job until
the previous one ends successfully, e.g.

```
A=$(PROJ_ID=nek-vf QUEUE=prod PREPARE_RESTART=1 <this dir>/submit-bundle.sh 03:50 0.444x:25)
B=$(PROJ_ID=nek-vf QUEUE=prod PREPARE_RESTART=1 DEPEND=$A <this dir>/submit-bundle.sh 03:50 0.444x:25)
```

A job whose case fails exits non-zero, so the rest of the chain is not run.
Short walltimes that end before the next large job in `qstat -T` often start
at once as backfill (`pbsnodes -aS` shows the free nodes).

A single run can also be submitted with nrsqsub_polaris, e.g. for bringup in
the debug queue (without BUBBLE_STOP_AT):

```
cd ~/answinter26/Sc1-debug/2.25x
QUEUE=debug PROJ_ID=nek-vf nrsqsub_polaris bubble3d.par 1 1:00
```

## Postprocessing

* `sherwood.py --tmin 20 --tmax 30 <runs>`: window-averaged Sh from data.csv.
* `drift.py <runs>`: when the bubble leaves its initial rise axis, and the
  exponential growth rate of the lateral offset.
* `submit-animate.sh` and `animate.py`: animation.mp4 and post.csv (centroid,
  velocities, interface area, shape) from the checkpoints, with ParaView on a
  compute node.
* `add-coords.py`: makes a checkpoint self-contained (nekRS writes the mesh
  coordinates only into each job's first checkpoint), so it can seed a run on
  another mesh with `startFrom = <file>+int`.

## Results (Polaris, September 2026)

### 2.25x, 1.5x and 1x

All three used dt = 1e-4, except 2.25x, which blew up at t = 24.87 and was
continued from t = 24.8 at dt = 5e-5. After the animations, the field files
were deleted to make room for the finer meshes; answinter26/Sc1/<run>/ keeps
the inputs, data.csv, post.csv, animation.mp4, logs/ and bubble3d_t15.fld, a
self-contained checkpoint at t = 15.0001.

Sh averaged over t = 20-30 (100 checkpoint intervals, `sherwood.py`):

| Run   | Sh    | batch SE | std  | area | c_bulk |
|-------|-------|----------|------|------|--------|
| 2.25x | 7.84  | 0.19     | 0.62 | 4.04 | 0.967  |
| 1.5x  | 8.40  | 0.21     | 0.65 | 3.34 | 0.978  |
| 1x    | 11.31 | 0.09     | 0.32 | 3.42 | 0.974  |

This comparison is confounded by path instability. The bubble rises straight
at first, then drifts sideways out of its own periodic wake, and Sh rises
when it does. The lateral offset of the 1.5x and 1x bubbles grows
exponentially at nearly the same rate, but from seeds 300 times apart, so the
drift starts about 10 time units apart (`drift.py`):

| Run   | offset > 0.01 | > 0.05 | > 0.25 | growth rate | fit at t = 0 |
|-------|---------------|--------|--------|-------------|--------------|
| 2.25x | 1.1           | 1.5    | 2.4    | 1.76        | 2e-3         |
| 1.5x  | 22.2          | 24.7   | 27.3   | 0.67        | 3e-9         |
| 1x    | 12.5          | 14.8   | 17.5   | 0.74        | 1e-6         |

(times in time units, rates per time unit, offsets in bubble diameters; the
last column is the exponential fit over 1e-3 < offset < 0.1 extrapolated
back to t = 0).
The seed is numerical noise, so the onset is not a property of the mesh
resolution alone. Over 20-30 the 1x bubble has drifted but the 1.5x bubble
has only started to. Comparing the two while both still rise straight, Sh
(± batch standard error) differs by about 10% rather than 35%:

| Window  | 2.25x      | 1.5x         | 1x           |
|---------|------------|--------------|--------------|
| 5-10    | 16.8 ± 4.1 | 10.01 ± 0.12 | 11.02 ± 0.23 |
| 10-15   | 12.8 ± 2.9 | 8.64 ± 0.07  | 9.46 ± 0.05  |
| 15-20   | 7.01 ± 0.23| 8.30 ± 0.04  | 9.82 ± 0.30  |
| 20-25   | 7.34 ± 0.10| 8.06 ± 0.03  | 11.39 ± 0.02 |
| 25-30   | 8.34 ± 0.14| 8.75 ± 0.37  | 11.24 ± 0.19 |

On the 2.25x mesh (3 elements across the bubble) the bubble drifts almost at
once and is far more deformed (area 4.04 vs 3.34-3.42, aspect ratio 0.9 vs
1.6-1.8), rises at 0.64 instead of about 1.0, and before t = 1 its gas core
carried velocities 20-40 times the rise velocity along the bubble axis (about
5 times on 1x), so its values are not reliable. In the t = 20-30 window the
maximum gas velocity is 3.4-4.7 times the rise velocity on the 1.5x and 1x
meshes (5.8-8.8 on 2.25x).

### Finer meshes

0.667x runs from t = 0 in answinter26/Sc1/0.667x on 10 nodes. Level set
reinitialization dominates the cost: CLSR runs every 20 steps and TLSR every
200 at dt = 5e-5, and on 10 nodes together they take 40% of the time. Seconds
per step over steps 200-2000 (t = 0.01-0.1) of short test runs:

| Nodes (GPUs) | GLL points/GPU | mean   | ordinary step | CLSR step | TLSR step |
|--------------|----------------|--------|---------------|-----------|-----------|
| 2 (8)        | 5.5M           | 0.226  | 0.117         | 1.78      | 6.08      |
| 10 (40)      | 1.1M           | 0.0738 | 0.0445        | 0.486     | 1.75      |
| 24 (96)      | 0.46M          | 0.0586 | 0.0374        | 0.354     | 1.31      |

A run to t = 30 is 600,000 steps: at these rates 38 hours (75 node-hours) on
2 nodes, 12 hours (123) on 10 and 9.8 hours (244, with 25 nodes charged) on
24-25. On the 1x mesh the cost per step after the first 0.1 time units was
about 25% higher, so these are lower bounds.

0.667x checkpoints are 1.84 GB (2.6 GB for each job's first, which also holds
the coordinates), about 110 GB for t = 0-30 at 0.5. The 0.444x and 0.296x
checkpoints are 3.4 and 11 times larger, so 0.444x at 0.5 needs about 370 GB
and 0.296x at 0.25 about 2.5 TB, against the 1 TB nek-vf quota on eagle; their
checkpoint intervals (or which checkpoints are kept) have to change before
they run.
