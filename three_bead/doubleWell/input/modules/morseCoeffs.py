import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, minimize_scalar

def morsePotential(r, r0, alpha, d):
    """morse potential to represent lj but has controllable well width and depth
    - d in energy units
    - alpha in 1/distance units
    - r0 in distance units"""
    return d * (np.exp(-2 * alpha * (r - r0)) - 2 * np.exp(- alpha * (r - r0)))



def morseCoeffs(depth, width, r0):
    """
    calculates alpha for morse potential
    depth = maximum well depth below zero energy
    width = well width measured at E = -depth/2
    r0 = position of the minimum 
    """
    alpha = np.log(3 + 2*np.sqrt(2)) / width
    print(f"pair_coeff * * {depth:.6g} {alpha:.6g} {r0:.6g}")

    fig, ax = plt.subplots()
    vals = np.arange(0, 15, 0.01)
    ax.plot(vals, [morsePotential(r, r0, alpha, depth) for r in vals])
    plt.axhline(y = -0.5 * depth, linestyle = "--")
    plt.show()

    return depth, alpha, r0


def writeCoeffFile(depth, width, conditions):
    parDiameter = 2 * 2**(1/6)
    width = width * parDiameter
    if width < 2.5 * parDiameter:
        r0 = width
    else:
        r0 = 0.5 * width
    d, alpha, r0 = morseCoeffs(depth, width, r0)
    with open(f"../runs/{conditions}/morseInfo.in", "w") as f:
        f.write(f"# pair style morse coeffs \n")
        f.write(f"pair_style morse 14\n")
        f.write(f"pair_coeff * * {d:.6g} {alpha:.6g} {r0:.6g}\n")
