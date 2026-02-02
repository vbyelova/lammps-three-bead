import matplotlib.pyplot as plt
import numpy as np
from numpy import random
from scipy.optimize import minimize, root
import os

from doubleWellCoeffs import *
from generateMonomer import *
from generateDimer import *
from generateThreeBead import *
from generateLammpsInput import *

#generate optimal values for quartic function

# first input value is unfolding barrier height
# second input value is minimum energy for unfolded state

barrierToUnfold = 3
barrierToRefold = 2

numRuns = 1
Vf = 0.07

boxLength = 20
numPar = (boxLength**3 * Vf) / ((1 * 2 **(1 / 6))**3 * (np.pi / 6))
numMol = int(numPar / 3) + 1
prob = 1 # probability of unfolding
bondsPerAtom = 2 

conditions = f"unfold{barrierToUnfold}_refold{barrierToRefold}_Vf{Vf}_mol{numMol}"

quarticVals = calcDoubleWellCoeffs(barrierToUnfold,barrierToRefold)
angle0, k2, k3, k4 = quarticVals[0], quarticVals[1], quarticVals[2], quarticVals[3]

#generate molecule input file


files = []
for n in range (0, numRuns):
    filename = f"Run{n}_{conditions}"
    files.append(filename)    

if not os.path.exists("../runs"):
        os.makedirs("../runs")

if not os.path.exists(f"../runs/{conditions}"):
     os.makedirs(f"../runs/{conditions}")
     writeAngleFile(angle0, k2, k3, k4, conditions)

for filename in files:
    if not os.path.exists(f"../runs/{conditions}/{filename}/input"):
        os.makedirs(f"../runs/{conditions}/{filename}/input")
        generateThreeBead(f"../runs/{conditions}/{filename}/input/{filename}.in", numMol, boxLength)
    if not os.path.exists(f"../runs/{conditions}/{filename}/output"):
        os.makedirs(f"../runs/{conditions}/{filename}/output")
    if not os.path.exists(f"../runs/{conditions}/{filename}/analysis"):
        os.makedirs(f"../runs/{conditions}/{filename}/analysis")
        generateLammpsInput(conditions, filename, boxLength, numMol, prob, bondsPerAtom)

with open(f"../runs/{conditions}/runScripts.py", "w") as f:
    f.write(f"# run this script from {conditions} to run the simulations.\n")
    f.write(f"import subprocess\n"
            f"import os\n\n")
    for filename in files:
        f.write(f"if not os.path.exists(f'../runs/{conditions}/{filename}/input/log.lammps'):\n")
        f.write(f"\tprint('running simulation {filename}..')\n")
        f.write(f"\tsubprocess.run(['lmp','-in','in.lammps'], cwd = '{filename}/input')\n\n")

