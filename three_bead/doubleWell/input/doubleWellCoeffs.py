# a script to find the optimal coefficient values based on desired energy barrier height.

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize, root, bisect

def quartic(theta, theta0, k2, k3, k4):
    return k2 * (theta - theta0)**2 + k3 * (theta - theta0)**3 + k4 * (theta - theta0)**4

def quarticFirstDeriv(theta, theta0, k2, k3, k4):
    return 2 * k2 * (theta - theta0) + 3 * k3 * (theta - theta0)**2 + 4 * k4 * (theta - theta0)**3

def quarticSecondDeriv(theta, theta0, k2, k3, k4):
    return 2 * k2 + 6 * k3 * (theta - theta0) + 12 * k4 * (theta - theta0)**2


def objective(params, theta_min1, theta_min2, barrierToUnfold, barrierToRefold):
    """An objective function to pass to scipy's minimizing function that optimizes
        quartic function coefficients.
        Checks if suggested minima and maxima are correct and optimizes value
        of equilibrium theta as well."""
    
    error = 0
    theta0, k2, k3, k4 = params

    # checking first minimum values
    yAtMin1 = quartic(theta_min1, theta0, k2, k3, k4)
    min1firstDeriv = quarticFirstDeriv(theta_min1, theta0, k2, k3, k4)
    min1secondDeriv = quarticSecondDeriv(theta_min1, theta0, k2, k3, k4)
    if min1secondDeriv < 0:
        error += np.abs(min1secondDeriv) * 1000
    error += min1firstDeriv**2 


    # checking second minimum values
    yAtMin2 = quartic(theta_min2, theta0, k2, k3, k4)
    min2firstDeriv = quarticFirstDeriv(theta_min2, theta0, k2, k3, k4)
    min2secondDeriv = quarticSecondDeriv(theta_min2, theta0, k2, k3, k4)
    if min2secondDeriv < 0:
        error += np.abs(min2secondDeriv) * 1000
    error += min2firstDeriv**2 


    # from folded to unfolded (left to right minima)
    yAtTheta0 = quartic(theta0, theta0, k2, k3, k4)
    error += (np.abs(yAtMin1 - yAtTheta0) - barrierToUnfold)**2

    # from unfolded back to folded (right to left minima)
    error += (np.abs(yAtMin2 - yAtTheta0) - barrierToRefold)**2

    print("error: ", error)
    print("args,", params)

    return error

def calcDoubleWellCoeffs(barrierToUnfold, barrierToRefold):

    # parameters
    # first minimum
    theta_min1 = np.pi / 3

    # second minimum 
    theta_min2 = np.pi


    # initial guesses from simultaneous equations + tinkering
    #k2 = -0.0055
    #k3 = 5E-6
    #k4 = 1.5E-6

    k2 = -10.3
    k3 = 2.14
    k4 = 4.8
    
    theta0 = 13/18 * np.pi

    initialGuess = [theta0, k2, k3, k4]
    counter = 0

    print("Finding optimal coefficients for quartic function...")
    results = minimize(objective, initialGuess, args = (theta_min1, theta_min2, barrierToUnfold, barrierToRefold), method = "Nelder-Mead", options={'maxiter':1000})
    theta0, k2, k3, k4 = results.x[0], results.x[1], results.x[2], results.x[3]
    print
    energy = []
    thetas = list(np.arange(2/9 * np.pi, 10/9 * np.pi, 0.01))

    for theta in thetas:
        energy.append(quartic(theta, theta0, k2, k3, k4))
    
    print("Found optimal values!")
    print(f"Theta_0 = {theta0}\n")
    print(f"k2 = {k2}\n")
    print(f"k3 = {k3}\n")
    print(f"k4 = {k4}\n")


    plt.plot(thetas, energy, color = "purple")
    plt.xlabel("theta")
    plt.ylabel("Potential energy")
    plt.title("Asymmetric double well for two-state protein")
    plt.show()

    return theta0, k2, k3, k4

def writeAngleFile(theta0, k2, k3, k4):
    with open("angleInfo.in", "w") as f:
        f.write(f"# angle coefficents as calculated \n")
        f.write(f"angle_style lepton no_offset \n")
        f.write(f"angle_coeff 1 {theta0} 'k2*(theta - {theta0})^2 + k3*(theta - {theta0})^3 + k4*(theta - {theta0})^4; k2={k2}; k3={k3}; k4={k4}' \n")
        print("written angle coefficients to file!")
    return
