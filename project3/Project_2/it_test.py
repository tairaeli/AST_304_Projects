import numpy as np
from structure import integrate, pressure_guess
from astro_const import Msun, Rsun, G
import matplotlib.pyplot as plt
import scipy.optimize as op
from eos import density

mue = 2

S_mass = Msun

Pc = pressure_guess(S_mass, mue)

# print(Pc)

delta_m = 1e-100

eta = 1

xi = 1

m_step, r_step, p_step = integrate(Pc,delta_m,eta,xi,mue)

print(p_step[0])

print(density(p_step[0], mue))

print(m_step[-1])
# if eta and xi are at optimal, changing delta_m doesn't do much
# past xi ~ 1.75, w/ optimal eta and delta_m, start getting negative pressure
# lower than xi ~ 0.003, get too many iterations
