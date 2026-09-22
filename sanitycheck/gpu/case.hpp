// FLiBe(625C,1bar)/Ar(625C,1bar)/3.12mm/0.29 m/s + density ratio set to 40
static double Re = 231.4;
static double Sc = 1150;
static double Pe = Re*Sc;
static double Fr = 1.633;
static double We = 2.667;
static double rhoratio = 40;
static double muratio = 148.1;
static double nuratio = muratio/rhoratio;
static double diffratio = 5.698e-6;
static double solubilityratio = 0.002804;

// Static case info struct to be passed to device kernels.
static caseinfo_t info;
