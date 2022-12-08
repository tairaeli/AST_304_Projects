import numpy as np
from structure import central_thermal, integrate
from eos import mean_molecular_weight
from zams import Teff, surface_luminosity
from astro_const import Rsun, Msun
from scipy.optimize import brentq

r0 = 0.000000001

delta_m = 1e-30

eta = 9e-20

xi = 0.05

XH = 0.706

Z = np.array([1,2,7])
A = np.array([1,4, 14])
X = np.array([XH, 0.275, 0.019])

mu = mean_molecular_weight(Z,A,X)

# S_mass = 0.3
R = 0.33333


def min_this(R_guess, S_mass, delta_m, eta, xi, mu):

    Pc, rhoc, Tc = central_thermal(S_mass, R_guess, mu)

    m, r, p, l = integrate(Pc,rhoc,Tc,delta_m,eta,xi,mu)

    # calculating desired radius
    L_surface = surface_luminosity(Teff(S_mass*Msun),R_guess*Rsun)

    return l[-1] - L_surface

mass_list = np.linspace(0.1,0.31,10)
# mass_list = [0.1]
L_list = np.zeros(len(mass_list))
Teff_list = np.zeros(len(mass_list))
rhoc_list = np.zeros(len(mass_list))
Tc_list = np.zeros(len(mass_list))

for i,S_mass in enumerate(mass_list):
    R_final = brentq(min_this,a = 0.01, b = 1, xtol = 1e-10,args = (S_mass, delta_m, eta, xi, mu))

    T_surf = Teff(S_mass)

    Teff_list[i] = T_surf

    L_list[i] = surface_luminosity(T_surf,R_final)

    Pc, rhoc, Tc = central_thermal(S_mass, R_final, mu)

    rhoc_list[i] = rhoc

    Tc_list[i] = Tc

print(R_final)

Pc, rhoc, Tc = central_thermal(S_mass, R_final, mu)

m, r, p, l = integrate(Pc,rhoc,Tc,delta_m,eta,xi,mu)

print(R*Rsun - r[-1])