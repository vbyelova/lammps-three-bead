import matplotlib.pyplot as plt
import numpy as np
from numpy import random
from scipy.optimize import minimize, root
import os

from modules.doubleWellCoeffs import *
from modules.generateMonomer import *
from modules.generateThreeBead import *
from modules.generateLammpsInput import *

#generate optimal values for quartic function

# first input value is unfolding barrier height
# second input value is minimum energy for unfolded state

barrierToUnfold = [10]
barrierToRefold = 2

numRuns = 2
vf = [0.07]

boxLength = 30
prob = 1 # probability of unfolding
bondsPerAtom = 2

simType = "step strain"
strainFractions = [0.05, 0.07]

#generate molecule input file

for volfrac in vf:
    for barrier in barrierToUnfold:
        numPar = (boxLength**3 * volfrac) / ((1 * 2 **(1 / 6))**3 * (np.pi / 6))
        numMol = int(numPar / 3) + 1
        conditions = f"unfold{barrier}_refold{barrierToRefold}_Vf{volfrac}_mol{numMol}"

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
            for strainFraction in strainFractions:
                if simType == "step strain":
                    prefix = f"stepStrain{strainFraction}_"
                elif simType == "gelation":
                    prefix = ""
                if not os.path.exists(f"../runs/{conditions}/{filename}/input"):
                    os.makedirs(f"../runs/{conditions}/{filename}/input")
                    generateThreeBead(f"../runs/{conditions}/{filename}/input/{filename}.in", numMol, boxLength)
                if not os.path.exists(f"../runs/{conditions}/{filename}/output"):
                    os.makedirs(f"../runs/{conditions}/{filename}/output")
                generateLammpsInput(conditions, filename, boxLength, numMol, prob, bondsPerAtom)

                if not os.path.exists(f"../runs/{conditions}/{filename}/{prefix}output"):
                    os.makedirs(f"../runs/{conditions}/{filename}/{prefix}output")

                if not os.path.exists(f"../runs/{conditions}/{filename}/{prefix}analysis"):
                    os.makedirs(f"../runs/{conditions}/{filename}/{prefix}analysis")

                if simType == "step strain":
                    stepStrainInput(conditions, filename, boxLength, numMol, prob, bondsPerAtom, strainFraction)

        with open(f"../runs/{conditions}/runScripts.py", "w") as f:
            f.write(f"# run this script from {conditions} to run the simulations.\n")
            f.write(f"import subprocess\n"
                    f"import os\n\n")
            for filename in files:
                f.write(f"if not os.path.exists(f'../runs/{conditions}/{filename}/input/log.lammps'):\n")
                f.write(f"\tprint('running gelation simulation {filename}..')\n")
                f.write(f"\tsubprocess.run(['lmp','-in','in.lammps'], cwd = '{filename}/input')\n\n")
                for strainFraction in strainFractions:
                    f.write(f"\tprint('running strain fraction = {strainFraction} simulation {filename}..')\n")
                    f.write(f"\tsubprocess.run(['lmp','-in','stepStrain{strainFraction}_in.lammps'], cwd = '{filename}/input')\n\n")
