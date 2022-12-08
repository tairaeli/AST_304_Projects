import numpy as np
from structure import central_thermal, integrate
from eos import mean_molecular_weight

# seems like changing the r0 significantly alters the resulting arrays
# not really sure what to do with this information
r0 = 0.000000001

delta_m = 1e-35

eta = 1e-30

xi = 0.05

XH = 0.706

Z = np.array([1,2,7])
A = np.array([1,4, 14])
X = np.array([XH, 0.275, 0.019])

mu = mean_molecular_weight(Z,A,X)

M = 0.3
R = 0.3

Pc, rhoc, Tc = central_thermal(M, R, mu)

m, r, p, l = integrate(Pc,rhoc,Tc,delta_m,eta,xi,mu)

print(r)
print(l)