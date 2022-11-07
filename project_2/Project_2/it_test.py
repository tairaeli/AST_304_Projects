import numpy as np
from structure import integrate, pressure_guess
from astro_const import Msun
import matplotlib.pyplot as plt

mue = 2

S_mass = Msun

Pc = pressure_guess(S_mass, mue)

delta_m = 1e5

eta = 1e-24

xi = 0.1

m_step = [0]

tolerance = 0.01

m_step, r_step, p_step = integrate(Pc,delta_m,eta,xi,mue)

print(m_step[-1])

plt.plot(p_step)
# plt.xscale("log")
# plt.yscale("log")
plt.show()
# print(S_mass)
# print(len(m_step))

# print(p_step)

