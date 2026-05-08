import matplotlib.pyplot as plt
import numpy as np
from numpy import random
from scipy.optimize import minimize, root
import os

from modules.doubleWellCoeffs import *
from modules.generateMonomer import *
from modules.generateDimer import *
from modules.generateThreeBead import *
from modules.generateLammpsInput import *
from modules.generateRDF import *

#generate optimal values for quartic function

# first input value is unfolding barrier height
# second input value is minimum energy for unfolded state

barrierToUnfold = [1, 2, 5]
barrierToRefold = 2

numRuns = 1
vf = [0.04]

boxLength = 100
prob = 1 # probability of unfolding
bondsPerAtom = 2 

#generate molecule input file

# for volfrac in vf:
#     for barrier in barrierToUnfold:
#         numPar = (boxLength**3 * volfrac) / ((1 * 2 **(1 / 6))**3 * (np.pi / 6))
#         numMol = int(numPar / 3) + 1
#         conditions = f"unfold{barrier}_refold{barrierToRefold}_Vf{volfrac}_mol{numMol}"

#         files = []
#         for n in range (0, numRuns):
#             filename = f"Run{n}_{conditions}"
#             quarticVals = calcDoubleWellCoeffs(barrier, barrierToRefold)
#             angle0, k2, k3, k4 = quarticVals[0], quarticVals[1], quarticVals[2], quarticVals[3]

threeExamples([4], 2)
