# mesh-sensitivity-Sc1

Mesh sensitivity study of the Sherwood number for the gpu/3D bubble case at
Sc = 1 with the physical FLiBe/Ar density ratio, at 2.25x, 1.5x and 1x the
Kolmogorov scale.

Physics follows gpu/3D (4x8x4 fully periodic domain, bubble of diameter 1
starting at y = 1, Re = 231.4, Fr = 1.633, We = 2.667, muratio = 148.1) with:
* Sc = 1, so Pe = Re = 231.4 (flibe1/ar in cpu/common/calc_dimensionless_numbers.py).
* Physical density ratio rhoratio = 1/0.0002709 = 3691 instead of 40.

Species transport and the Sherwood number follow gpu/match-nek5000-202605
rather than gpu/3D's continuous species transfer (CST) flux integral:
* Uniform diffusivity D = 1/Pe in both phases, no CST terms.
* Hard sink: every step, c is removed where psi < 0.5 and the removed amount
  is accumulated; mdot is its average over each checkpoint interval.
* Soft source: c is driven towards 1 in the liquid bulk (psi > 0.9).
* MTC = mdot/(total_area*c_bulk) and Sh = MTC*Pe, written to data.csv every
  checkpoint (0.1 time units).

Numerics follow gpu/3D except:
* Polynomial order 7 and TLSR/CLSR at the nekRS-LS defaults (no targetCFL,
  maximumSteps or SVV), as recommended by the level set developer; TLSR runs
  every 100 and CLSR every 10 timesteps.
* Each mesh is generated directly by genbox (no hrefine, which leaves the
  partition imbalanced when the base mesh has few elements per rank).
* The same dt = 2e-4 is used for every mesh so only the spatial resolution
  changes. With the physical density ratio, dt = 4e-4 gave spurious velocity
  spikes in the gas core (on element edges along the bubble axis) near
  t = 0.35 on the 2.25x and 1.5x meshes, and the 2.25x run blew up at t = 0.72.

## Meshes

The Kolmogorov scale is lambda_k = 0.0671 mm (lambda_k/d = 0.021536). Mesh
resolution is measured by the mean unique GLL spacing dx = h/N at N = 7, with
element size h = 4/Nelx.

| Run   | Elements | dx/lambda_k | GLL points (E*N^3) |
|-------|----------|-------------|--------------------|
| 2.25x | 12x24x12 | 2.21        | 1.19M              |
| 1.5x  | 18x36x18 | 1.47        | 4.00M              |
| 1x    | 27x54x27 | 0.98        | 13.5M              |

## Running on Polaris

Create the run directories and meshes (genbox must be on PATH):

```
./setup-runs.sh ~/answinter26/Sc1
```

The meshes are small compared to the 10-node minimum of the prod queue, so
all three run concurrently in one job, each on its own nodes:

```
cd ~/answinter26/Sc1
PROJ_ID=nek-vf QUEUE=prod <this dir>/submit-bundle.sh 03:00 1x:7 1.5x:2 2.25x:1
```

Level set reinitialization dominates the cost. Measured on 4 GPUs, a step
averages 0.065 s (2.25x), 0.098 s (1.5x) and 0.24 s (1x), so each run needs
roughly 2.5-3 hours on the layout above. The small prod queue allows at most
3 hours. If a job ends before endTime = 30, prepare the unfinished runs for a
restart and resubmit them:

```
<this dir>/prepare-restart.sh 1x
PROJ_ID=nek-vf QUEUE=prod <this dir>/submit-bundle.sh 03:00 1x:10
```

prepare-restart.sh moves the finished part's field files to part<N>/ (nekRS
restarts checkpoint numbering at 0), trims data.csv back to the restart time
and sets startFrom.

A single run can also be submitted with nrsqsub_polaris, e.g. for bringup in
the debug queue:

```
cd ~/answinter26/Sc1-debug/2.25x
QUEUE=debug PROJ_ID=nek-vf nrsqsub_polaris bubble3d.par 1 1:00
```
