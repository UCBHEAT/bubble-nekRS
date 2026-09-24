# Ready-to-submit: 2D DNS production run

Command (from this directory, on Aurora):
    nekls && QUEUE=prod ./requeue.sh 2d-dns 256 47:00

Constraints & sizing
- prod queue minimum is 256 nodes (auto-routes to large/small); 47:00 is
  the max wall per allocation, requeue.sh chains restarts on
  afternotok using checkpointInterval = 0.5 sim-time.
- 256 nodes = 1536 tiles -> 1.76M GLL/tile for this 12.5M-elem N=5 mesh.
  That is under the ~8M GLL/tile scaling optimum (64 nodes would be the
  optimum) but the floor is enforced; total cost is GLL-seconds, roughly
  unchanged: ~500-700 node-days for t=30 (measured 0.474 s/step at 32
  nodes incl. TLSR/CLSR redistancing; 3M steps).
- Do NOT lower dt below... dt = 1e-5 is the ONLY demonstrated stable dt at
  full mesh (3e-5 diverged step 98, 5e-5 step 26, 1e-4 step 3-4).
- Validated config: N=5 (N=7 pressure GMRES stagnates), interfaceWidthFactor
  = 15, TLSR/CLSR maximumSteps 400/250, all-P box + boundaryTypeMap=none,
  autotuning off via UDF_Setup0, genmap not needed (mesh script removed it).
- Qualification: 3300 steps at dt=1e-5 clean (CFL flat 0.022 after
  startup transient), job 8856397.

Known open issue
- The CFL-triggered spurious acceleration from the level-set coupling on
  the 1-cell periodic-z slab (vacuum twin clean; interface twin diverges
  at target CFL). dt=1e-5 keeps us below it, but revisit if long-run CFL
  drifts upward near the bubble wake.
