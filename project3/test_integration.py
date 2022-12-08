import numpy as np
from structure import central_values, integrate
from eos import mean_molecular_weight

# seems like changing the r0 significantly alters the resulting arrays
# not really sure what to do with this information
r0 = 0.000000001

delta_m = 1e-4

eta = 1e-30

xi = 0.05

XH = 0.706

Z = np.array([1,2,7])
A = np.array([1,4, 14])
X = np.array([XH, 0.275, 0.019])

mu = mean_molecular_weight(Z,A,X)

m, r, p, l = integrate(r0, delta_m, eta, xi, mu, XH, max_steps=10000)

print(r)