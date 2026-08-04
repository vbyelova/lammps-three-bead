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

barrierToUnfold = [1, 5]
barrierToRefold = 2

numRuns = 10
vf = [0.07]

boxLength = 55
prob = 1 # probability of unfolding
bondsPerAtom = 5

#generate molecule input file

for volfrac in vf:
    for barrier in barrierToUnfold:
        numPar = (boxLength**3 * volfrac) / ((1 * 2 **(1 / 6))**3 * (np.pi / 6))
        numMol = int(numPar / 3) + 1
        conditions = f"unfold{barrier}_refold{barrierToRefold}_Vf{volfrac}_mol{numMol}_bondsPerAtom{bondsPerAtom}"

        files = []
        for n in range (0, numRuns):
            filename = f"Run{n}_{conditions}"
            files.append(filename)    

        if not os.path.exists("../runs"):
                os.makedirs("../runs")

        if not os.path.exists(f"../runs/{conditions}"):
            os.makedirs(f"../runs/{conditions}")
            quarticVals = calcDoubleWellCoeffs(barrier, barrierToRefold)
            angle0, k2, k3, k4 = quarticVals[0], quarticVals[1], quarticVals[2], quarticVals[3]
            writeAngleFile(angle0, k2, k3, k4, conditions)

        for filename in files:
            if not os.path.exists(f"../runs/{conditions}/{filename}/input"):
                os.makedirs(f"../runs/{conditions}/{filename}/input")
                generateThreeBead(f"../runs/{conditions}/{filename}/input/{filename}.in", numMol, boxLength)
            if not os.path.exists(f"../runs/{conditions}/{filename}/output"):
                os.makedirs(f"../runs/{conditions}/{filename}/output")
                generateLammpsInput(conditions, filename, boxLength, numMol, prob, bondsPerAtom)
            if not os.path.exists(f"../runs/{conditions}/{filename}/analysis"):
                os.makedirs(f"../runs/{conditions}/{filename}/analysis")


        with open(f"../runs/{conditions}/runScripts.py", "w") as f:
            f.write(f"# run this script from {conditions} to run the simulations.\n")
            f.write(f"import subprocess\n"
                    f"import os\n\n")
            for filename in files:
                f.write(f"if not os.path.exists(f'../runs/{conditions}/{filename}/input/log.lammps'):\n")
                f.write(f"\tprint('running simulation {filename}..')\n")
                f.write(f"\tsubprocess.run(['lmp','-in','in.lammps'], cwd = '{filename}/input')\n\n")

