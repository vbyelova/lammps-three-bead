# a script to find the optimal coefficient values based on desired energy barrier height.

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize, root

def quartic(x, x0, k2, k3, k4, y_max):
    return k2 * (x - x0)**2 + k3 * (x - x0)**3 + k4 * (x - x0)**4 + y_max

def quarticFirstDeriv(x, x0, k2, k3, k4):
    return 2 * k2 * (x - x0) + 3 * k3 * (x - x0)**2 + 4 * k4 * (x - x0)**3

def quarticSecondDeriv(x, x0, k2, k3, k4):
    return 2 * k2 + 6 * k3 * (x - x0) + 12 * k4 * (x - x0)**2

def objective(params, x_min1, y_min1, x_min2, y_min2, y_max):
    """An objective function to pass to scipy's minimizing function that optimizes
        quartic function coefficients.
        Checks if suggested minima and maxima are correct and optimizes value
        of equilibrium angle as well."""
    
    error = 0
    x0, k2, k3, k4 = params

    # checking first minimum values
    yAtMin1 = quartic(x_min1, x0, k2, k3, k4, y_max)
    min1firstDeriv = quarticFirstDeriv(x_min1, x0, k2, k3, k4)
    min1secondDeriv = quarticSecondDeriv(x_min1, x0, k2, k3, k4)
    if min1secondDeriv < 0:
        error += np.abs(min1secondDeriv) * 1000
    error += min1firstDeriv**2 * 1000
    error += (yAtMin1 - y_min1)**2 * 1000

    # checking second minimum values
    yAtMin2 = quartic(x_min2, x0, k2, k3, k4, y_max)
    min2firstDeriv = quarticFirstDeriv(x_min2, x0, k2, k3, k4)
    min2secondDeriv = quarticSecondDeriv(x_min2, x0, k2, k3, k4)
    if min2secondDeriv < 0:
        error += np.abs(min2secondDeriv) * 1000
    error += min2firstDeriv**2 * 1000
    error += (yAtMin2 - y_min2)**2 * 1000

    # checking maximum value / energy barrier height
    yAtMax = quartic(x0, x0, k2, k3, k4, y_max)
    maxFirstDeriv = quarticFirstDeriv(x0, x0, k2, k3, k4)
    maxSecondDeriv = quarticSecondDeriv(x0, x0, k2, k3, k4)
    if maxSecondDeriv > 0:
        error += np.abs(maxSecondDeriv) * 1000
    error += maxFirstDeriv**2 * 1000
    error += (yAtMax - y_max)**2 * 1000    

    return error

# parameters

# first minimum
x_min1 = 60
y_min1 = 0

# second minimum 
x_min2 = 180
y_min2 = 2

# first maximum i.e. barrier height

y_max = 10

# initial guesses

k2 = -0.0055
k3 = 5E-6
k4 = 1.5E-6
angle0 = 130

initialGuess = [angle0, k2, k3, k4]
counter = 0

results = minimize(objective, initialGuess, args = (x_min1, y_min1, x_min2, y_min2, y_max), method = "Nelder-Mead", options={'maxiter':1000})
angle0, k2, k3, k4 = results.x[0], results.x[1], results.x[2], results.x[3]


energy = []
angles = list(np.arange(40, 190, 1))

for angle in angles:
    energy.append(quartic(angle, angle0, k2, k3, k4, y_max))

plt.plot(angles, energy, "x")
plt.show()
