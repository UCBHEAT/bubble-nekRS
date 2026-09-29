// FLiBe(625C,1bar)(Sc=1)/Ar(625C,1bar)/3.12mm/0.29 m/s, physical density ratio
// (flibe1/ar in cpu/common/calc_dimensionless_numbers.py)
static double Re = 231.4;
static double Sc = 1;
static double Pe = Re*Sc;
static double Fr = 1.633;
static double We = 2.667;
static double rhoratio = 1.0/0.0002709;
static double muratio = 148.1;
static double nuratio = muratio/rhoratio;

// Static case info struct to be passed to device kernels.
static caseinfo_t info;
