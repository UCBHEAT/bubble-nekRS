# Comparison with Fan, Fang & Bolotnov (2020) "bubble PID"

Reference: Y. Fan, J. Fang, I. Bolotnov, *Complex bubble deformation and
break-up dynamics studies using interface capturing approach*, Experimental
and Computational Multiphase Flow 3(3), 139–151 (2020),
doi:10.1007/s42757-020-0073-3 (open access, CC-BY). PHASTA implementation
following Thomas et al. (2015), the original PHASTA PID bubble controller.

This note compares their controller with ours (`gpu/common/pid.hpp`, used by
`gpu/pid-centering` and `gpu/pid-centering-3D`) and records the changes that
comparison suggests.

## What they do

Same physical idea: switch from the container frame to a bubble-based frame so
the bubble sits at a fixed point in a static mesh while uniform liquid flows
past at the terminal velocity. The control force is a spatially uniform body
force (their Eqs. 7, from Thomas et al. 2015):

    F_i(n+1) = c1*F̄_i(n) + c2*[F̄_i(n) + c3*x_i(n) + c4*v_i(n) + c5*ẋ_i(n)]   (their Eq. 7)

- x_i is the centroid displacement from the setpoint, v_i and ẋ_i the bubble
  velocity and acceleration; F̄_i is a history-averaged control force playing
  the role of the integral term. The five coefficients c1..c5 are tuned, not
  derived from a loop analysis; no stability criterion is given.
- Two derivative channels (velocity and acceleration) instead of one.
- Domain: 10R cube, bubble at x = 3.3R from the inlet (streamwise x), uniform
  inflow, four cylindrical refinement shells, 52 elements/D chosen by GCI
  (~3.8% grid uncertainty on the film radius). ρl:ρg = 1000, μl:μg = 100.
- The controller output is a *measurement*: at steady state F_x balances
  drag, so C_D = 2|F_x|/(ρl A u²) is reported (their Eqs. 8–9), with V and A
  taken as the initial-sphere volume/cross-section. Terminal velocity becomes
  an input, removing its experimental uncertainty from Re.
- Verification: Popinet's DNS bubble topologies (Tripathi et al. 2015),
  Bhaga–Weber shapes (Eo/Mo match ±1%/±8%, Re deviates more at high
  deformation), and Sharaf et al. film/break-up experiments.

## What we do

Same frame-shift idea, different controller realisation (and different
discretisation: NekRS spectral elements vs PHASTA unstructured FEM):

    F_pid = -(kp*e + ki*∫e dt + kd*u_gas),  e = x_gas − x_0

- e is the gas centroid offset (CLS-weighted), ∫e dt a true trapezoid-free
  time integral, u_gas the gas mean velocity used as ė (derivative on
  measurement, so CLSR reinit centroid jumps don't enter the D term; the I
  term absorbs the residual drift instead).
- Gains are placed by closed-loop pole assignment: the F_pid→centroid plant is
  a double integrator, characteristic polynomial s³ + kd s² + kp s + ki;
  kd > 0, ki ≥ 0, kd·kp > ki is checked at setup (Routh–Hurwitz), and the
  defaults kp = 12, ki = 8, kd = 6 put a triple pole at s = −2.
- The frame shift is explicit and kinematically consistent: the uniform
  acceleration is added to the momentum source and the inflow gets
  U_frame(t) = ∫F_pid dt (forward-Euler updated one step ahead), so far-field
  liquid moves at exactly the frame velocity; O(dt) mismatches are absorbed by
  the uniform pressure gradient (~0.1% of gravity measured).
- Controller state (integral, U_frame) is checkpointed to data.csv and
  restored on restart, so restarts do not perturb the closed loop.
- 3D (gpu/pid-centering-3D): same controller, one setpoint per coordinate;
  our setup differs physically in that buoyancy (−1/Fr²) stays on and the
  controller *centers* rather than replaces it; PHASTA turns container gravity
  off and the control force is the entire substitute for buoyancy.

## Assessment against our stated goals

The two implementations are equivalent in purpose and both reach the same
steady state (bubble fixed, uniform counterflow at U_term). Our version is in
better shape analytically (explicit plant model, verifiable stability
condition, restart consistency) and their version is in better shape as a
metrology instrument (the control force is *the* reported drag, with
verification against experiments and another DNS code). Three changes follow.

## Changes adopted in this branch

1. **Steady-state force-balance check + C_D report** (their Eqs. 8–9 are the
   controller's payoff once it converges). Note the convention difference:
   PHASTA runs without container gravity, so their steady state is
   F_pid = buoyancy and the reported drag is F_x directly. We keep −1/Fr²
   explicit and the controller centers, so the frame co-moves at constant
   velocity once converged and F_pid → 0; the drag on the bubble is then the
   residual buoyancy of the gas volume. `pidSteadyReport()` in
   `gpu/common/pid.hpp` computes, per checkpoint window,
   - `F_drag = (1 − 1/rhoratio)·V_gas/Fr² + mean(F_pid)·V_gas/rhoratio`, and
   - `C_D = 2·F_drag/(A·U_term²)` with U_term = |rise| (the measured lab-frame
     rise speed, already in data.csv) and A = `[PID] frontalArea` if set, else
     the frontal area of the equal-volume sphere.
   and prints both; `F_drag,C_D` are appended to data.csv (column layout
   change: delete the old data.csv before restarting a case). The window-mean
   of F_pid doubles as the convergence indicator and is printed as a % of the
   full buoyancy. In 2D the area (and hence drag and C_D) are per unit z
   depth; `gpu/pid-centering/bubble.par` sets `frontalArea = 0.25` (D = 1
   cylinder × 0.25 periodic depth) and warns if it is missing. Diagnostic
   only; nothing about the dynamics changed.
2. **Document their alternative derivative/filter choice** (pid.hpp header
   comment): their two derivative channels (v and dv/dt, Eq. 7) damp the loop
   without a large kd on velocity alone. We keep the single-channel form
   (clean pole placement, no acceleration-measurement noise) but record the
   fallback recipe.
3. **Domain-length sanity note** for mesh work: their 10R cube with the bubble
   3.3R from the inlet is much tighter than our 6D-wide/9D-tall domain;
   streamwise maps to our −y (flow enters the top). Read `mesh-sensitivity/`
   and `sc1` results against that reference.

## Things we deliberately do NOT copy

- Tuned c1..c5 recursion (their Eq. 7): it hides the integral in a force
  average, has no stated stability map, and our pole-placement form already
  meets the offset target (peak ~0.02 D near t = 1).
- Density/viscosity ratios 1000/100: our FLiBe/Ar case physics fixes 40/148.
- Treating gravity as absent: our case reports Sh against a buoyancy-driven
  rise, so keeping −1/Fr² explicit and centering around it is required.
