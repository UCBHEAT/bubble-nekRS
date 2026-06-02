// FLiBe(625C,1bar)/Ar(625C,1bar)/3.12mm/0.29 m/s
static double Re = 231.4;
static double Sc = 1150;
static double Pe = Re*Sc;
static double Fr = 1.633;
static double We = 2.667;
static double rhoratio = 0.025;//0.0002709;
static double muratio = 24.92;
static double nuratio = muratio/rhoratio;
static double diffratio = 1.755e5;
static double solubilityratio = 356.7;

// Static case info struct to be passed to device kernels.
static caseinfo_t info;
