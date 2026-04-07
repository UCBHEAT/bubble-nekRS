#include <cmath>

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

// Common case information struct, to make device kernel arguments
// more readable.
typedef struct caseinfo {
    // Number of elements local to this CPU node
    dlong Nelements;
    // Offset between vector field components
    dlong fieldOffset;
    // Element length
    dfloat deltael;
    // Target conservative level set (CLS) interface width (epsilon)
    dfloat eps_cls;
} caseinfo_t;

/**
 * Apply the smoothed Heaviside step function to input phi, mapping
 * (-infty, infty) to (0, 1), with thickness epsilon.
 *
 * @param phi input to produce a step function at the zero level set
 * @param epsilon step function width
 * @returns output between 0 and 1
 */
inline double heaviside(double phi, double deltael, double epsilon) {
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
