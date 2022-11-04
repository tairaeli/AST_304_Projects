import numpy as np
from structure import integrate, pressure_guess
from astro_const import Msun

mue = 2

S_mass = Msun

Pc = pressure_guess(S_mass, mue)

delta_m = 1e6

eta = 1e-10

xi = 0.0001

m_step = [0]

tolerance = 0.01

while S_mass-m_step[-1] > tolerance:
    m_step, r_step, p_step = integrate(Pc,delta_m,eta,xi,mue)

