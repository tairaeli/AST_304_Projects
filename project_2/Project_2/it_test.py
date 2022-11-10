import numpy as np
from structure import integrate, pressure_guess
from astro_const import Msun, Rsun, G
import matplotlib.pyplot as plt
import scipy.optimize as op

mue = 2

S_mass = Msun

Pc = pressure_guess(S_mass, mue)

# print(Pc)

delta_m = 1e-4

eta = 1e-10

xi = 0.05

m_step, r_step, p_step = integrate(Pc,delta_m,eta,xi,mue)

print(m_step[-1])