#include <cmath>

// Must be set to the average element length.
// TODO: calculate from mesh on-demand, as in
// https://github.com/nandu90/Nek5000/blob/nekLS/core/experimental/lvlSet.f#L587
const double deltael = 0.05;

// Recommended value is epsilon = 1.5/(lx1-1.0), e.g.
//   polynomial order 5: 0.3
//   polynomial order 7: 0.2143
// etc.
// TODO: calculate from actual polynomial order (mesh->Np) rather than hardcoding.
const double eps_cls = 0.2143;

// FLiBe(625C,1bar)/Ar(625C,1bar)/3.12mm/0.29 m/s
const double Re = 231.4;
const double Fr = 1.633;
const double We = 2.667;
const double Sc = 1150;
const double Pe = 231.4*1150; // = 2.66e+05
const double Mo = 2.479e-09;
const double rhoratio = 0.0002709;
const double nuratio = 24.92;
const double muratio = 24.92*0.0002709;
const double diffratio = 1.755e+05;
const double solubilityratio = 356.7;

/**
 * Apply the smoothed Heaviside step function to input phi, mapping
 * (-infty, infty) to (0, 1), with thickness epsilon.
 */
inline double heaviside(double phi, double epsilon) {
    return 0.5*(std::tanh(phi/(2.0*deltael*epsilon))+1.0);
}
