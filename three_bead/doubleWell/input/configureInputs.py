import matplotlib.pyplot as plt
import numpy as np
from numpy import random
from scipy.optimize import minimize, root
import os

from doubleWellCoeffs import *
from generateMonomer import *
from generateDimer import *
from generateThreeBead import *

#generate optimal values for quartic function

# first input value is unfolding barrier height
# second input value is minimum energy for unfolded state

barrierToUnfold = 2
barrierToRefold = 1.5

numRuns = 3
Vf = 0.07

quarticVals = calcDoubleWellCoeffs(barrierToUnfold,barrierToRefold)

angle0, k2, k3, k4 = quarticVals[0], quarticVals[1], quarticVals[2], quarticVals[3]

writeAngleFile(angle0, k2, k3, k4)

#generate molecule input file

#first value is number of molecules
#second value is length of box

numMol = 3938
boxLength = 50

conditions = f"langevin_10_unfold{barrierToUnfold}_refold{barrierToRefold}_Vf{Vf}_mol{numMol}"

files = []
for n in range (0, numRuns):
    filename = f"Run{n}_{conditions}"
    files.append(filename)    

if not os.path.exists("../runs"):
        os.makedirs("../runs")

if not os.path.exists(f"../runs/{conditions}"):
     os.makedirs(f"../runs/{conditions}")

for filename in files:
    if not os.path.exists(f"../runs/{conditions}/{filename}/input"):
        os.makedirs(f"../runs/{conditions}/{filename}/input")
        generateThreeBead(f"../runs/{conditions}/{filename}/input/{filename}.in", numMol, boxLength)
    if not os.path.exists(f"../runs/{conditions}/{filename}/output"):
        os.makedirs(f"../runs/{conditions}/{filename}/output")
    if not os.path.exists(f"../runs/{conditions}/{filename}/analysis"):
        os.makedirs(f"../runs/{conditions}/{filename}/analysis")

   

with open(f"{conditions}.txt", "w") as f:
    for filename in files:
        f.write(f"{filename}\n")

