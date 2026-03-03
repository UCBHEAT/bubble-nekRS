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

// Water/air/artificial solute
const double charl = 0.01;
const double gdim = 9.81;
const double rhol = 998.2;
const double rhog = 1.122;
const double mul = 1e-3;
const double mug = 1.824e-5;
const double diffl = 0.2e-4;
const double diffg = 1.0e-4;
static double udim = diffl/charl; // diffusive scaling
static double Re = rhol*udim*charl/mul;
static double Sc = mul/rhol/diffl; // = 1/Re in diffusive scaling
static double Pe = Re*Sc; // = 1 in diffusive scaling
static double rhoratio = rhog/rhol;
static double muratio = mug/mul;
static double nuratio = muratio/rhoratio;
static double diffratio = diffg/diffl;
static double solubilityratio = 3.0;

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
