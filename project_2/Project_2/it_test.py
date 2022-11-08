import numpy as np
from structure import integrate, pressure_guess
from astro_const import Msun
import matplotlib.pyplot as plt
import scipy.optimize as op

mue = 2

S_mass = Msun

Pc = pressure_guess(S_mass, mue)

delta_m = 1e3

eta = 1e-24

xi = 0.1

tolerance = 0.01

m_step, r_step, p_step = integrate(Pc,delta_m,eta,xi,mue)

print(m_step[-1])

plt.plot(m_step)
# plt.xscale("log")
# plt.yscale("log")
plt.show()
# print(S_mass)
# print(len(m_step))

# print(p_step)

# def min_this(Pc):

#     m_step, r_step, p_step = integrate(Pc,delta_m,eta,xi,mue)

#     m_result = m_step[-1]

#     m_desired = 1.98e30

#     return m_desired - m_result