import numpy as np

from collections import defaultdict
from calcAngles import *
from calcFracDim import * 
from calcPercolation import *

from parseBonds import *

#finalFrame(f"../runs/{conditions}/{filename}/analysis/dillPercolation.pkl", totalBonds, Npar)
#systemData = parseSystemData(f"../runs/{conditions}/{filename}/output/systemData.txt")
#boxLength, Npar, Nsteps, Nwrite, equilTime = systemData[0], systemData[1], systemData[2], systemData[3], systemData[4]
#totalBonds = totalBonds(f"../runs/{conditions}/{filename}/output/nbonds.dat")

#calcAngles(totalBonds, f"../runs/{conditions}/{filename}/output/moleculeangles.dat", Npar)

unfoldBarriers = [3]
refoldBarrier = 2
numRuns = 5
Vf = 0.07
numMol = 3938

#checkBondLength(unfoldBarriers, refoldBarrier, numRuns, Vf, numMol)

checkClumping(unfoldBarriers, refoldBarrier, numRuns, Vf, numMol)

# fractalDims = findAllFractalDims(unfoldBarriers, refoldBarrier, numRuns, Vf, numMol)
# finalFractalDims, finalFractalDimsError = fractalDims[0], fractalDims[1]

# plt.bar(unfoldBarriers, finalFractalDims, yerr = finalFractalDimsError, color = "purple")
# plt.xlabel("Unfolding barrier (kT)")
# plt.ylabel("Fractal dimension")
# plt.title("Fractal dimension through box counting method")
# plt.show()

# angles = calcAngles(3, refoldBarrier, Vf, numMol)
# print(len(angles))
# hist, bins = np.histogram(angles, bins = 18)
# logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
# plt.hist(angles, bins = logbins, color = "pink")
# plt.legend(["Vf = 0.07\nunfolding barrier  = 3 kT\n n. molecules = 851"])
# plt.xlabel("molecule angle \u03B8")
# plt.ylabel("log(number of molecules)")
# plt.title("Bimodal distribution of folded and unfolded molecules")
# plt.show()

