import matplotlib.pyplot as plt
import numpy as np

def ljYukawa(r, a, b, k):
    return 4 * a * ((b / r)**12 - (b / r)**6) + a - (a * b) * (np.exp(- k * r) / r) - a * b

def ljYukawaFirstDeriv(r, a, b, k):
    return -48 * (a / b) * ((b / r)**13 - (b / r)**7) + (a * b) * (k**2 * np.exp(- k * r) / r**2)

def ljYukawaSecondDeriv(r, a, b, k):
    return 48 *(a / (b**2)) * (13 * (b / r)**14 - 3.5 * (b / r)**8) - (a * b) * (k**3 * np.exp(-k * r) / r**3)

def objective(params, r, a, b, k):
    """objective function to use scipy minimize and find values of epsilon, sigma and kappa
        such that the LJ well depth and width can be controlled.
        Conditions for successful values:
        - """

epsilon = 1
sigma = 1
kappa = 1
vals = np.arange(0.001, 3, 0.001)
ax = plt.gca()
ax.set_ylim(-5, 5)
ax.plot(vals, [ljYukawa(val, epsilon, sigma, kappa) for val in vals])
plt.show()
