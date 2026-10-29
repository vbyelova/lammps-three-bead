import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter
import matplotlib.colors as mcolors

from scipy.optimize import minimize, minimize_scalar, brentq
from itertools import cycle

def ljYukawa(r, epsilon, sigma, prefactor, kappa):
    lj = 4 * epsilon * ((sigma/r)**12 - (sigma/r)**6) + epsilon
    yukawa = prefactor * np.exp(- kappa * r) / r
    return lj - yukawa


def wellVals(epsilon, sigma, kappa, prefactor, depth_target):
    """minimizes to find values at energetic minimum and also well width, r2 - r1"""
    res = minimize_scalar(lambda r: ljYukawa(r, epsilon, sigma, prefactor, kappa), bounds = (0.5 * sigma, 3 * sigma),
                          method = "bounded")

    r_min, u_min = res.x, res.fun
    if u_min >= 0 or r_min > 2.99 * sigma:
        return None
    half = u_min / 2
    f = lambda r: ljYukawa(r, epsilon, sigma, prefactor, kappa) - half
    try:
        r1 = brentq(f, 0.5 * sigma, r_min)
        r2 = brentq(f, r_min, 10 * sigma)
    except ValueError:
        return None
    return - u_min, r_min, r2 - r1

def objective(params, targets, sigma):
    """objective function to use scipy minimize and find values of epsilon, sigma and kappa
        such that the LJ well depth and width can be controlled.
        Conditions for successful values:
        - epsilon has to be at dE/dr = 0, d^2E/dr^2 > 0
        - r_1 = r_2 at E = 0.5 epsilon
        - r_1 - r_2 = sigma
        
        a = epsilon, b = sigma"""

    depth_target, width_target = targets
    epsilon, kappa = np.exp(params[0]), np.exp(params[2])
    prefactor = params[1]
    vals = wellVals(epsilon, sigma, kappa, prefactor, depth_target)
    if vals is None:
        return 1e6
    depth, rmin, width = vals
    return((depth - depth_target) / depth_target)**2 + ((width - width_target)/ width_target)**2

def fit(objective, params0, args, tol = 1e-12, max_restarts = 50, n_random = 20, seed = 0):
    rng = np.random.default_rng(seed)
    best = None

    starts = [np.array(params0)]
    for _ in range(n_random):
        starts.append(np.array(params0) + rng.normal(scale=1.0, size=len(params0)))

    for p0 in starts:
        sol = minimize(objective, p0, args=args, method="Nelder-Mead",
                       options=dict(xatol=1e-10, fatol=1e-14, maxiter=5000))
        if best is None or sol.fun < best.fun:
            best = sol

    for i in range(max_restarts):
        if best.fun < tol:
            break
        sol = minimize(objective, best.x, args=args, method="Nelder-Mead",
                       options=dict(xatol=1e-12, fatol=1e-16, maxiter=5000))
        improved = best.fun - sol.fun
        if sol.fun < best.fun:
            best = sol
        if improved < 1e-3 * best.fun:   
            break


    return best



params0 = [np.log(0.5), np.log(1.0), np.log(2.0)]
sigma = 1

parDiameter = 2 * 2**(1/6)
parRadius = 2 ** (1/6)
depth_target = 1.4
width_target = 1 * parDiameter

targets = (depth_target, width_target)
solution = fit(objective, params0, args = (targets, sigma))
print(solution.fun)
#solution = minimize(objective, params0, args = (targets, sigma), method = "Nelder-Mead")
epsilon, prefactor, kappa = np.exp(solution.x[0]), solution.x[1], np.exp(solution.x[2])
vals = np.arange(0.001, 5, 0.001)
# ax = plt.gca()
# ax.set_ylim(-5, 5)
# ax.plot(vals, [ljYukawa(val, epsilon, sigma, prefactor, kappa) for val in vals])
# plt.show()

depth_targets = [0.1, 0.1, 1.5, 1.5]
width_targets = [1 * parRadius, 4 * parRadius, 1 * parRadius, 4 * parRadius]

ax = plt.gca()
ax.set_ylim(-5, 5)
colours = ["pink", "green", "blue", "purple"]
for i, (w, d) in enumerate(zip(width_targets, depth_targets)):
    targets = (d, w)
    solution = fit(objective, params0, args = (targets, sigma))
    sigma = 1
    eps, pref, k = np.exp(solution.x[0]), solution.x[1], np.exp(solution.x[2])
    vals = np.arange(0.001, 5, 0.001)
    ax.plot(vals, [ljYukawa(val, eps, sigma, pref, k) for val in vals], color = colours[i], label = f"{targets}")
plt.show()
    