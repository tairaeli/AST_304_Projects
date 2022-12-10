import sys
sys.path.append('.\Project_2')
from p2_structure import integrate, pressure_guess
from p2_eos import density
import matplotlib.pyplot as plt
import scipy.optimize as op
from astro_const import Msun, Rsun, G, Ke
import numpy as np

def mass_finder(Pc, wanted_mass, delta_m, eta, xi, mue):
    '''
    Function to be included into our rootfinding method.
    Runs the integrate function, then returns the difference
    between the final mass found by the integration and the
    specified desired mass
    '''
    
    m_step, r_step, p_step = integrate(Pc,delta_m,eta,xi,mue)

    return wanted_mass - m_step[-1]

# setting up array of desired white dwarf masses
mass_array = np.arange(0.1, 1.1, 0.1)

# Setting some initial conditions
mue = 2

eta = 1e-10

xi = 0.05

Smass = Msun * 0.1
    
# setting delta_m here as it depends on desired mass
delta_m = Smass*1e-6

K_mod_list = np.linspace(0.5,2,3)

for K_mod in K_mod_list:

    print(K_mod, "* K:")

    P_guess = pressure_guess(Smass,Ke*K_mod, mue = 2)

    # running the rootfinding method to get ideal central pressure
    desired_pressure = op.bisect(mass_finder, a = P_guess*1e-5, b = P_guess*1e5, xtol = 1e-10, args = (Smass, delta_m, eta, xi, mue))

    m, r, p = integrate(desired_pressure,delta_m,eta,xi,mue)

    # assigning data to arrays
    alpha = desired_pressure/(G*Smass**2*r[-1]**(-4))

    beta = density(desired_pressure, mue)/(3*Smass/(4*np.pi*r[-1]**3))

    print(r"\alpha =",alpha, r"\beta =",beta)