#include <cmath>

// Must be set to the average element length.
// TODO: calculate from mesh on-demand, as in
// https://github.com/nandu90/Nek5000/blob/nekLS/core/experimental/lvlSet.f#L587
const double deltael = 0.08333333333;
/*
    dfloat deltael = 0.0;
    for (dlong i = 0; i < mesh->Nelements; i++) {
       deltael += tmp[i];
    }
    MPI_Allreduce(MPI_IN_PLACE, &deltael, 1, MPI_DFLOAT, MPI_SUM, platform->comm.mpiComm);
    dlong N = mesh->Nelements;
    MPI_Allreduce(MPI_IN_PLACE, &N, 1, MPI_DLONG, MPI_SUM, platform->comm.mpiComm);
    deltael /= N;
    if(platform->comm.mpiRank == 0)
      printf("Average element length = %.4e\n",deltael);
  }
*/

// Recommended value is epsilon = 1.5/(lx1-1.0), e.g.
//   polynomial order 5: 0.3
//   polynomial order 7: 0.2143
// etc.
// TODO: calculate from actual polynomial order (mesh->Np) rather than hardcoding.
const double eps_cls = 0.2143;

// FLiBe(625C,1bar)/Ar(625C,1bar)/3.12mm/0.29 m/s
static double Re = 231.4;
static double Sc = 1150;
static double Pe = Re*Sc;
static double Fr = 1.633;
static double rhoratio = 0.002709;
static double muratio = 24.92;
static double nuratio = muratio/rhoratio;
static double diffratio = 1.755e5;
static double solubilityratio = 356.7;

/**
 * Apply the smoothed Heaviside step function to input phi, mapping
 * (-infty, infty) to (0, 1), with thickness epsilon.
 *
 * @param phi input to produce a step function at the zero level set
 * @param epsilon step function width
 * @returns output between 0 and 1
 */
inline double heaviside(double phi, double epsilon) {
    return 0.5*(std::tanh(phi/(2.0*deltael*epsilon))+1.0);
}

/**
 * Reimplement std::clamp(value, min, max) as an inline function to
 * work around disallowed call to constexpr host functions from device
 * code.
 *
 * @param value input
 * @param min minimum
 * @param max maximum
 * @returns output clamped to [min, max]
 */
inline double clip(double value, double min, double max) {
    if (value < min) return min;
    if (value > max) return max;
    return value;
}
