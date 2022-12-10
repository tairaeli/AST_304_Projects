import numpy as np
from structure import integrate, pressure_guess
from astro_const import Msun
import matplotlib.pyplot as plt
import scipy.optimize as op
from eos import pressure, density

mue = 2

S_mass = Msun*0.1

Pc = pressure_guess(S_mass, mue)

delta_m = S_mass*1e-6

eta = 1e-10

xi = 0.003

def min_this(Pc, S_mass, delta_m, eta, xi, mue):

    m_step, r_step, p_step = integrate(Pc,delta_m,eta,xi,mue)

    return S_mass - m_step[-1]

solve_pressure = op.bisect(min_this,a = 1e19, b = 1e23, xtol = 1e-10,args = (S_mass, delta_m, eta, xi, mue))

print(solve_pressure)

print()

new_m_step, new_r_step, new_p_step = integrate(solve_pressure,delta_m,eta,xi,mue)

print(S_mass - new_m_step[-1])