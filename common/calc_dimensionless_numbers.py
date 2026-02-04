from dataclasses import dataclass, field
from typing import Optional

import numpy as np

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
flinak = Phase("FLiNaK(625C,1bar)", rho=2019, D=1.07566e-8, kH=3.7966e-5, nu=4.1537e-3/2019) # see footnotes
water = Phase("water(20C,1bar,oxygen)", rho=998.2, D=2.01e-9, kH=1.243e-5, nu=1.e-3/998.2) # from [Marschall12], see footnotes for kH validation
ar = Phase("Ar(625C,1bar)", rho=0.535, D=5.8677e-4, kH=1./(R*898.15), nu=9.578e-5) # from NIST WebBook, H2-Ar diffusion coefficient from [MarreroMason72]
he = Phase("He(625C,1bar)", rho=0.05359, D=1.5547e-4, kH=1./(R*898.15), nu=4.2785e-5/0.05359) # from NIST WebBook, H2-He diffusion coefficient from [MarreroMason72], see footnotes
air = Phase("air(20C,1bar,oxygen)", rho=1.122, D=1.916e-5, kH=1./(R*293.15), nu=1.824e-5/1.122) # from [Marschall12]
#water_25C = Phase("water(25C,1bar,CO2)", rho=997., D=1.97e-9, kH=3.49e-4, nu=8.927e-7)
#co2_25C = Phase("CO2(25C,1bar,CO2)", rho=1.784, D=1.381e-5, kH=1./(R*298.15), nu=8.369e-6)

combos = [
    # liquid, gas, surface tension, bubble diameter, bubble velocity
    (flibe, ar, 0.188, None, None),
    (flibe, he, 0.188, None, None),
    (flinak, ar, 0.182, None, None), # see footnotes
    (water, air, 7.2e-2, 2e-3, 0.232),
    (water, air, 7.2e-2, 4e-3, 0.221),
    (water, air, 7.2e-2, 6e-3, 0.208),
]

def calc_dimensionless_numbers(liquid, gas, sigma, d, u):
    if d is None:
        # Use capillary length.
        d = np.sqrt(sigma/((liquid.rho - gas.rho)*g))

    if u is None:
        if False:
            """
            These formulae all seem to overestimate significantly...

            Hadamard-Rybczynski is supposed to consider dissipation
            of energy in both liquid and gas.
            
            Stokes' law is supposed to consider dissipation of energy
            in the liquid.
            """
            # Calculate terminal velocity (Hadamard-Rybczynski eq'n)
            u = (2./3)*((d/2)**2)*g*(liquid.rho - gas.rho)/liquid.mu
            u *= (liquid.mu + gas.mu)/(2*liquid.mu + 3*gas.mu)

            # Calculate terminal velocity (Stokes' law)
            u = 2*((d/2)**2)*g*(liquid.rho-gas.rho)/(9*liquid.mu)
            
        # Calculate terminal velocity (simple)
        V = (4./3)*np.pi*((d/2)**3)
        A = np.pi*((d/2)**2)
        Cd = 0.5 # sphere drag coefficient
        u = np.sqrt(2*(liquid.rho-gas.rho)*V*g / (liquid.rho*A*Cd))

    # Calculate correlations
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
    if True:
        # Calculate terminal velocity (simple)
        V = (4./3)*np.pi*((d/2)**3)
        A = np.pi*((d/2)**2)
        Cd = 0.5 # sphere drag coefficient
        u2 = np.sqrt(2*(liquid.rho-gas.rho)*V*g / (liquid.rho*A*Cd))
        print(f"! u from simple calc = {u2}")
        # Calculate Kolmogorov scale
        epsilon = u*g
        lambda_k = ((liquid.nu**3)/epsilon)**0.25
        lambda_kd = ((liquid.D**3)/epsilon)**0.25
        print(f"! Bubble rise specific turbulent KE dissipation rate = {epsilon} J/kg")
        print(f"! Kolmogorov scale lambda_k = {lambda_k*10**3} mm")
        print(f"! Mass transfer Kolmogorov scale lambda_kd = {lambda_kd*10**3} mm")

for combo in combos:
    calc_dimensionless_numbers(*combo)

# Footnotes

# Water-air-O2 at 20C kH
# [Marschall12] uses "He=33" (H=1/33); we confirm that here.
# NIST: kH(T) = kH(25C) exp(C1 * (1/T - 1/298.15))
#   kH(25C) = 1.3e-5 mol/(m^3 Pa)
#   C1 = 1500 K
# => kH(20C) = (1.3e-5 mol/(m^3 Pa)) exp(-1500*(1/293.15 - 1/298.15))
#            = 1.193e-5 mol/(m^3 Pa)
# which results in solubilityratio=34.39, not 33.
# To match Marschall12's values, we use kH=1.35e-5 mol/(m^3 Pa).

# FLiNaK properties from [RomatoskiHu17].
# rho = 2579-0.624T
#     = 2019
# mu = 0.04exp(4170/T) in cP (centipoise) = mPa*s
#    = 4.1537e-3 Pa*s
# FLiNaK/H2 properties from [HumrickhouseFuerst20].
# kH = 4.0e-7*exp(34000/(RT))
#    = 3.7966e-5 mol/(m^3 Pa)
# D = 8.7e-10*exp(-50000/(RT))
#    = 1.07566e-12
# (this appears to be misspecified in the original Fukada and Morisaki paper
# as m^2/s when it is really cm^2/s, and Humrickhouse and Fuerst did not fix
# it. It is correct on everyone's plots.)
#    = 1.07566e-8 m^2/s
# FLiNaK surface tension from Janz (1988)
# sigma = 0.2726 - 1.014e-4*T
#       = 0.1815

# H2-He diffusion coefficient, [MarreroMason72] page 54 lower plot.
# p*D = 2.70e-2 T^(1.510)/((ln(T/5.34e6))^2) atm cm^2/s
#     = 1.53432178 atm cm^2/s
# D = 1.5547e-4 m^2/s
