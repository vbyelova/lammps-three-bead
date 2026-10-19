import matplotlib.pyplot as plt
import numpy as np

def ljYukawa(r, epsilon, sigma, kappa):
    return 4 * epsilon * ((sigma/r)**12 - (sigma/r)**6) + epsilon - (epsilon * sigma ) * (np.exp(- kappa * r) / r)

epsilon = 1
sigma = 1
kappa = 1
vals = np.arange(0.001, 3, 0.001)
ax = plt.gca()
ax.set_ylim(-5, 5)
ax.plot(vals, [ljYukawa(val, epsilon, sigma, kappa) for val in vals])
plt.show()
