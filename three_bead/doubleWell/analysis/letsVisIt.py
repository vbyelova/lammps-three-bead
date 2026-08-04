import numpy as np
import pickle
import re
import os
import subprocess
from collections import defaultdict

from modules.parseDump import *
from modules.calcAngles import *
from modules.calcFracDim import * 
from modules.calcPercolation import *
from modules.calcStressTensor import *
from modules.threeBeadClasses import *
from modules.parseRDF import *
from modules.loadForVis import *

# let's get the system data first
unfoldBarriers = [1, 5]
refoldBarrier = 2
numRuns = 1
vf = 0.07
numMol = 5243
boxLength = 55
suffixes = ["", ""]
bondsPerAtom = 5

for barrier, suffix in zip(unfoldBarriers, suffixes):
    conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}"
    for runNum in range(numRuns):
        filename = f"Run{runNum}_{conditions}"
        systemData = parseSystemData(f"../runs/{conditions}/{filename}/output/systemData.txt")

        boxLength, Npar, Nsteps, Nwrite, equilTime = systemData[0], systemData[1], systemData[2], systemData[3], systemData[4]

        # # let's parse the dump file first and save it for later

        # make directories for data and figures
        if not os.path.exists(f"../runs/{conditions}/{filename}/analysis/figs"):
            os.makedirs(f"../runs/{conditions}/{filename}/analysis/figs")
        if not os.path.exists(f"../runs/{conditions}/averagedFigs"):
            os.makedirs(f"../runs/{conditions}/averagedFigs")
        if not os.path.exists(f"../runs/{conditions}/data"):
            os.makedirs(f"../runs/{conditions}/data")
        if not os.path.exists(f"../runs/boxLength{boxLength}"):
            os.makedirs(f"../runs/boxLength{boxLength}")


        # check if particle trajectories have been parsed
        if not os.path.exists(f"../runs/{conditions}/{filename}/analysis/particleTraj.pkl"):
            data = readData(f"../runs/{conditions}/{filename}/output/vis.lammpstrj", Nsteps, Nwrite, Npar, equilTime)
            particles, timesteps = data[0], data[1]
            with open(f"../runs/{conditions}/{filename}/analysis/particleTraj.pkl", "wb") as f:
                pickle.dump(particles, f)
            with open(f"../runs/{conditions}/{filename}/analysis/timesteps.pkl", "wb") as f:
                pickle.dump(timesteps, f)
        if os.path.exists(f"../runs/{conditions}/{filename}/analysis/particleTraj.pkl"):
            with open(f"../runs/{conditions}/{filename}/analysis/particleTraj.pkl", "rb") as f:
                particles = pickle.load(f)
            with open(f"../runs/{conditions}/timesteps.pkl", "rb") as f:
                timesteps = pickle.load(f)




        angles = calcAngles(barrier, refoldBarrier, runNum, vf, numMol, bondsPerAtom, suffix)
        unfoldedMols = unfoldedMolecules(barrier, refoldBarrier, runNum, vf, numMol, angles, bondsPerAtom, suffix)     
        loadForceVis(barrier, refoldBarrier, runNum, vf, numMol,boxLength, particles, unfoldedMols, angles, bondsPerAtom, suffix)

#parCoordination = particleCoordination(barrier, refoldBarrier, runNum, vf, numMol, nBonds, bondInfo, timesteps)
#loadCoordVis(barrier, refoldBarrier, runNum, vf, numMol, particles, timesteps, parCoordination)

print("done :D")
