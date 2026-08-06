import argparse
import sys
import numpy as np
from dataclasses import dataclass
from typing import Optional

# Generate dimensionless number parameter files for the 6 cases in our test matrix.

g = 9.80665  # CGPM 1901 conventional value
R = 8.31446261815324  # ideal gas constant

@dataclass
class Phase:
    name: str
    rho: float
    D: float
    kH: float
    nu: float
    sigma: Optional[float] = None
    @property
    def mu(self):
        return self.rho*self.nu

flibe = Phase("FLiBe(625C,1bar)", rho=1975., D=3.343e-9, kH=3.75464e-7, nu=3.843e-6)
flinak = Phase("FLiNaK(625C,1bar)", rho=2019, D=1.07566e-8, kH=3.7966e-5, nu=4.1537e-3/2019)
water = Phase("water(20C,1bar,oxygen)", rho=998.2, D=2.01e-9, kH=1.243e-5, nu=1.e-3/998.2)
ar = Phase("Ar(625C,1bar)", rho=0.535, D=5.8677e-4, kH=1./(R*898.15), nu=9.578e-5)
he = Phase("He(625C,1bar)", rho=0.05359, D=1.5547e-4, kH=1./(R*898.15), nu=4.2785e-5/0.05359)
air = Phase("air(20C,1bar,oxygen)", rho=1.122, D=1.916e-5, kH=1./(R*293.15), nu=1.824e-5/1.122)
flibe1 = Phase("FLiBe(625C,1bar)(Sc=1)", rho=1975., D=3.843e-6/1, kH=3.75464e-7, nu=3.843e-6)
flibe4 = Phase("FLiBe(625C,1bar)(Sc=4)", rho=1975., D=3.843e-6/4, kH=3.75464e-7, nu=3.843e-6)
flibe20 = Phase("FLiBe(625C,1bar)(Sc=20)", rho=1975., D=3.843e-6/20, kH=3.75464e-7, nu=3.843e-6)
flibe100 = Phase("FLiBe(625C,1bar)(Sc=100)", rho=1975., D=3.843e-6/100, kH=3.75464e-7, nu=3.843e-6)

combos = [
    (flibe, ar, 0.188, None, None),
    (flibe, he, 0.188, None, None),
    (flinak, ar, 0.182, None, None),
    (water, air, 7.2e-2, 2e-3, 0.232),
    (water, air, 7.2e-2, 4e-3, 0.221),
    (water, air, 7.2e-2, 6e-3, 0.208),
    (flibe1, ar, 0.188, None, None),
    (flibe4, ar, 0.188, None, None),
    (flibe20, ar, 0.188, None, None),
    (flibe100, ar, 0.188, None, None),
]

case_lookup = {}
for combo in combos:
    liquid, gas, sigma, d, u = combo
    if liquid.name.startswith("FLiBe") and not "Sc=" in liquid.name and gas.name.startswith("Ar"):
        key = "FLiBe/Ar"
    elif liquid.name.startswith("FLiBe") and gas.name.startswith("He"):
        key = "FLiBe/He"
    elif liquid.name.startswith("FLiNaK") and gas.name.startswith("Ar"):
        key = "FLiNaK/Ar"
    elif liquid.name.startswith("water") and gas.name.startswith("air"):
        if d == 2e-3:
            key = "water/air (2mm)"
        elif d == 4e-3:
            key = "water/air (4mm)"
        elif d == 6e-3:
            key = "water/air (6mm)"
        else:
            key = "water/air"
    elif "Sc=1" in liquid.name:
        key = "FLiBe/Ar (Sc=1)"
    elif "Sc=4" in liquid.name:
        key = "FLiBe/Ar (Sc=4)"
    elif "Sc=20" in liquid.name:
        key = "FLiBe/Ar (Sc=20)"
    elif "Sc=100" in liquid.name:
        key = "FLiBe/Ar (Sc=100)"
    else:
        key = f"{liquid.name.split('(')[0]}/{gas.name.split('(')[0]}"
    case_lookup[key] = combo

def get_combo_by_case(name):
    if name not in case_lookup:
        print(f"Unknown case '{name}'. Available:")
        for k in sorted(case_lookup):
            print(f"  {k}")
        sys.exit(1)
    return case_lookup[name]

def calc_dimensionless_numbers(liquid, gas, sigma, d, u):
    if d is None:
        d = np.sqrt(sigma/((liquid.rho - gas.rho)*g))
    if u is None:
        V = (4./3)*np.pi*((d/2)**3)
        A = np.pi*((d/2)**2)
        Cd = 0.5
        u = np.sqrt(2*(liquid.rho-gas.rho)*V*g / (liquid.rho*A*Cd))
    Re = u*d/liquid.nu
    Sc = liquid.nu/liquid.D
    Mo = g*(liquid.mu**4)*(liquid.rho - gas.rho)/((liquid.rho**2)*(sigma**3))
    correlation_validity = (Re >= 1.) and (Re <= 5000.) and \
            (Re >= 3.73*(Mo*(-0.209))) and (Re <= 3.1*(Mo**(-0.25)))
    Brauer71 = 2.0 + 9.45e-4*(Re**1.07)*(Sc**0.888)
    HongBrauer84 = 2.0 + 1.5e-2*(Re**0.89)*(Sc**0.7)
    print(f"""
! {liquid.name}/{gas.name}/{d*1000:.2f}mm/{u:.2f} m/s
      real Re, Fr, We, Sc, Pe
      parameter (Re = {Re:.4g})
      parameter (Fr = {u/np.sqrt(g*d):.4g})
      parameter (We = {liquid.rho*(u**2)*d/sigma:.4g})
      parameter (Sc = {Sc:.4g})
      parameter (Pe = Re*Sc) ! = {Re*Sc:.4g}
      ! Mo = {Mo:.4g}
      real rhoratio, nuratio, muratio, diffratio, solubilityratio
      parameter (rhoratio = {gas.rho/liquid.rho:.4g})
      parameter (nuratio = {gas.nu/liquid.nu:.4g})
      parameter (muratio = nuratio*rhoratio)
      parameter (diffratio = {gas.D/liquid.D:.4g})
      parameter (solubilityratio = {gas.kH/liquid.kH:.4g})

      ! Correlations (valid={correlation_validity}):
      ! Brauer71 = {Brauer71:.4g}
      ! HongBrauer84 = {HongBrauer84:.4g}""")
    V2 = (4./3)*np.pi*((d/2)**3)
    A2 = np.pi*((d/2)**2)
    Cd2 = 0.5
    u2 = np.sqrt(2*(liquid.rho-gas.rho)*V2*g / (liquid.rho*A2*Cd2))
    epsilon = u*g
    lambda_k = ((liquid.nu**3)/epsilon)**0.25
    lambda_kd = ((liquid.D**3)/epsilon)**0.25
    print(f"""
      ! u from simple calc = {u2:.4g}
      ! Bubble rise specific turbulent KE dissipation rate = {epsilon:.4g} W/kg
      ! Kolmogorov scale lambda_k = {lambda_k*10**3:.4g} mm
      ! Mass transfer Kolmogorov scale lambda_kd = {lambda_kd*10**3:.4g} mm
      ! At polynomial order 7, {d/lambda_kd/8:.4g} elements per non-dim length unit
      ! 2x4x2 -> {d/lambda_kd/4:.0f}x{d/lambda_kd/2:.0f}x{d/lambda_kd/4:.0f}""")

def write_nek5000_case(filepath, liquid, gas, sigma, d, u):
    if d is None:
        d = np.sqrt(sigma/((liquid.rho - gas.rho)*g))
    if u is None:
        V = (4./3)*np.pi*((d/2)**3)
        A = np.pi*((d/2)**2)
        Cd = 0.5
        u = np.sqrt(2*(liquid.rho-gas.rho)*V*g / (liquid.rho*A*Cd))
    Re = u*d/liquid.nu
    Sc = liquid.nu/liquid.D
    Fr = u/np.sqrt(g*d)
    We = liquid.rho*(u**2)*d/sigma
    Pe = Re*Sc
    rhoratio = gas.rho/liquid.rho
    nuratio = gas.nu/liquid.nu
    muratio = nuratio*rhoratio
    diffratio = gas.D/liquid.D
    solubilityratio = gas.kH/liquid.kH
    content = f"""      ! Case-specific variables.

      ! Fluid parameters. As they are all constant we will just set
      ! them in this include file as parameters.

      ! Standard dimensionless numbers.
      ! {liquid.name}/{gas.name}/{d*1000:.2f}mm/{u:.2f} m/s
      real Re, Fr, We, Sc, Pe
      parameter (Re = {Re:.4g})
      parameter (Fr = {Fr:.4g})
      parameter (We = {We:.4g})
      parameter (Sc = {Sc:.4g})
      parameter (Pe = Re*Sc) ! = {Pe:.4g}

      ! These values are given as gas over liquid.
      real rhoratio, nuratio, muratio, diffratio, solubilityratio
      parameter (rhoratio = {rhoratio:.4g})
      parameter (nuratio = {nuratio:.4g})
      parameter (muratio = nuratio*rhoratio)
      parameter (diffratio = {diffratio:.4g})
      parameter (solubilityratio = {solubilityratio:.4g})

      ! Field mappings.
      integer ifld_v, ifld_cls, ifld_clsr, ifld_tls, ifld_tlsr, ifld_c
      parameter (ifld_v = 1)    ! velocity
      parameter (ifld_cls = 2)  ! temperature = conservative level set
      parameter (ifld_tls = 3)  ! scalar01 = traditional level set
      parameter (ifld_clsr = 4) ! scalar02 = CLS internal re-distancing
      parameter (ifld_tlsr = 5) ! scalar03 = TLS internal re-distancing
      parameter (ifld_c = 6)    ! scalar04 = concentration

      ! User parameter mappings
      integer iprm_tlsr_freq, iprm_clsr_freq, iprm_pord, iprm_cst_ver,
     $    iprm_sinkmode, iprm_sourcemode
      parameter (iprm_tlsr_freq = 1) ! uparam01 = TLS redistancing freq
      parameter (iprm_clsr_freq = 2) ! uparam02 = CLS redistancing freq
      parameter (iprm_pord = 3)      ! uparam03 = p extrapolation order
      parameter (iprm_cst_ver = 4)   ! uparam04 = CST version
      parameter (iprm_sinkmode = 5)  ! uparam04 = sink term mode
      parameter (iprm_sourcemode = 6)! uparam04 = source term mode

      ! Set a reasonable sink strength. This can be estimated from a
      ! target c_bubble:
      !   bubble_sink_out = bubble_interface_in
      !   sink_str * V * c_bubble_target = A * MTC *
      !                                    (c_bulk - H*c_bubble_target)
      !   sink_str = (A/V) * MTC * (c_bulk/c_bubble_target - H)
      ! For 2D, area/volume = pi*d/(pi*d^2/4) = 4/d^2 = 4.
      !   sink_str = 4 * (~200) * (1/c_bubble_target - 0.0028)
      !            = ~80 if targetting c_bubble ~ 10
      real sink_str
      parameter (sink_str = 0.05)

      ! Weaken liquid bulk source by this arbitrary ratio, or else the bulk
      ! returns to c=1 too fast and we don't visualize the wake of depleted
      ! concentration that the bubble leaves behind.
      real source_str
      parameter (source_str = 0.05)
"""
    with open(filepath, "w") as f:
        f.write(content)

def write_nekRS_case(filepath, liquid, gas, sigma, d, u):
    if d is None:
        d = np.sqrt(sigma/((liquid.rho - gas.rho)*g))
    if u is None:
        V = (4./3)*np.pi*((d/2)**3)
        A = np.pi*((d/2)**2)
        Cd = 0.5
        u = np.sqrt(2*(liquid.rho-gas.rho)*V*g / (liquid.rho*A*Cd))
    Re = u*d/liquid.nu
    Sc = liquid.nu/liquid.D
    Fr = u/np.sqrt(g*d)
    We = liquid.rho*(u**2)*d/sigma
    Pe = Re*Sc
    rhoratio = liquid.rho/gas.rho
    nuratio = liquid.nu/gas.nu
    muratio = (liquid.rho*liquid.nu)/(gas.rho*gas.nu)
    diffratio = liquid.D/gas.D
    solubilityratio = liquid.kH/gas.kH
    content = f"""// {liquid.name}/{gas.name}/{d*1000:.2f}mm/{u:.2f} m/s
static double Re = {Re:.4g};
static double Sc = {Sc:.4g};
static double Pe = Re*Sc;
static double Fr = {Fr:.4g}; // correct definition = u/sqrt(gL)
static double We = {We:.4g};
static double rhoratio = {rhoratio:.4g}; // liquid/gas density ratio
static double muratio = {muratio:.4g}; // liquid/gas viscosity ratio
static double nuratio = muratio/rhoratio; // derived
static double diffratio = {diffratio:.4g}; // liquid/gas diffusivity ratio
static double solubilityratio = {solubilityratio:.4g}; // liquid/gas solubility ratio

// Static case info struct to be passed to device kernels.
static caseinfo_t info;
"""
    with open(filepath, "w") as f:
        f.write(content)

def main():
    parser = argparse.ArgumentParser(
        description="Generate dimensionless number parameter files for bubble-nekRS.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Available combos: " + ", ".join(sorted(case_lookup))
    )
    parser.add_argument("--case", default="FLiBe/Ar", choices=list(case_lookup.keys()),
                        help="Select combo (default: FLiBe/Ar)")

    subparsers = parser.add_subparsers(dest="command", help="Subcommands")

    # Default behavior when no subcommand given: print reference values
    # We handle this by checking args.command after parse_args
    # But argparse subparsers make the subcommand optional only if we don't set required=True
    # Actually we want subcommand optional. Let's not set required.

    subparsers.add_parser("list", help="List available case combos")
    subparsers.add_parser("help", help="Print usage")

    gen_5000 = subparsers.add_parser("gen_nek5000_case", help="Generate Nek5000 CASE file (Fortran, gas/liquid)")
    gen_5000.add_argument("file", help="Output CASE file path")

    gen_nekRS = subparsers.add_parser("gen_nekRS_case", help="Generate NekRS case.hpp file (C++, liquid/gas)")
    gen_nekRS.add_argument("file", help="Output case.hpp file path")

    args = parser.parse_args()

    combo = get_combo_by_case(args.case)

    if args.command is None:
        calc_dimensionless_numbers(*combo)
    elif args.command == "list":
        print("Available combos:")
        for k in sorted(case_lookup):
            print(f"  {k}")
    elif args.command == "help":
        parser.print_help()
    elif args.command == "gen_nek5000_case":
        write_nek5000_case(args.file, *combo)
    elif args.command == "gen_nekRS_case":
        write_nekRS_case(args.file, *combo)

if __name__ == "__main__":
    main()
