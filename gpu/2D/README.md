Based on match-nek5000-202605.

- Moved to 2D, as Feng and Michaelides 2001 used 50 cells in the boundary layer, and that's not feasible in 3D.
  Note this still incurs (polynomial order + 1) compute overhead, as there are still GLL points in the z dir.
  For poly order 9, the 10x slowdown might defeat the GPU speedup... it certainly doesn't matter when doing
  small 2D cases but for production 2D cases nekRS is basically broken without real 2D support.

To do:
- Add implicit source term after Nadish pointed out how to do it in nekRS.
  https://github.com/nandu90/nekRS-LS/blob/nekLS\_dev/src/solver/scalar/scalarSolver.cpp#L865-L870
  Example implementation: https://github.com/nandu90/nekRS-LS/blob/nekLS\_dev/src/app/nrs/plugins/RANSktau.cpp#L411
  However it's nontrivial to implement. Need to
  1. Determine the convention, is it $a$ in the equation dc/dt = ac, or is it $-a$?
  2. Allocate and manage my own memory for this field (if needs to be -psi so I can't just return a slice of the psi field).
  3. Plumb through the userImplicitLinearTerm as it's called once per scalar ID. Function needs the ID of the scalar, not just the name.
     Completely different convention from all of the other callbacks like the explicit source term.
