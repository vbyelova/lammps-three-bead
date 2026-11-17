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
unfoldBarriers = [1]
refoldBarrier = 2
numRuns = 1
vf = 0.07
numMol = 3939
boxLength = 50
suffixes = [""]
bondsPerAtom = 2
prob = 1

simType = "step strain"
strainFraction = 0.2
prefix = f"stepStrain{strainFraction}_"
for barrier, suffix in zip(unfoldBarriers, suffixes):
    conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
    for runNum in range(numRuns):
        filename = f"Run{runNum}_{conditions}"
        systemData = parseSystemData(f"../runs/{conditions}/{filename}/{prefix}output/systemData.txt", simType)

        boxLength, Npar, Nsteps, Nwrite, equilTime = systemData[0], systemData[1], systemData[2], systemData[3], systemData[4]

        # # let's parse the dump file first and save it for later
        
        # make directories for data and figures
        if not os.path.exists(f"../runs/{conditions}/{filename}/{prefix}analysis/figs"):
            os.makedirs(f"../runs/{conditions}/{filename}/{prefix}analysis/figs")
        if not os.path.exists(f"../runs/{conditions}/averagedFigs"):
            os.makedirs(f"../runs/{conditions}/averagedFigs")
        if not os.path.exists(f"../runs/{conditions}/{prefix}data"):
            os.makedirs(f"../runs/{conditions}/{prefix}data")
        if not os.path.exists(f"../runs/boxLength{boxLength}"):
            os.makedirs(f"../runs/boxLength{boxLength}")


        # check if particle trajectories have been parsed
        if not os.path.exists(f"../runs/{conditions}/{filename}/{prefix}analysis/particleTraj.pkl"):
            data = readData(f"../runs/{conditions}/{filename}/{prefix}output/dump.lammpstrj", Nsteps, Nwrite, Npar, equilTime)
            particles, timesteps = data[0], data[1]
            with open(f"../runs/{conditions}/{filename}/{prefix}analysis/particleTraj.pkl", "wb") as f:
                pickle.dump(particles, f)
            with open(f"../runs/{conditions}/{filename}/{prefix}analysis/timesteps.pkl", "wb") as f:
                pickle.dump(timesteps, f)
        if os.path.exists(f"../runs/{conditions}/{filename}/{prefix}analysis/particleTraj.pkl"):
            with open(f"../runs/{conditions}/{filename}/{prefix}analysis/particleTraj.pkl", "rb") as f:
                particles = pickle.load(f)
            with open(f"../runs/{conditions}/{prefix}timesteps.pkl", "rb") as f:
                timesteps = pickle.load(f)
                print(timesteps)




        angles = calcAngles(barrier, refoldBarrier, runNum, vf, numMol, bondsPerAtom, prob, suffix, simType, strainFraction)
        unfoldedMols = unfoldedMolecules(barrier, refoldBarrier, runNum, vf, numMol, angles, bondsPerAtom, prob, suffix, simType, strainFraction)     
        loadForceVis(barrier, refoldBarrier, runNum, vf, numMol,boxLength, particles, unfoldedMols, angles, bondsPerAtom, prob, suffix, simType, strainFraction)
        #loadBondForceVis(barrier, refoldBarrier, runNum ,vf, numMol, boxLength, particles, bondsPerAtom, prob, suffix, simType, strainFraction)
#parCoordination = particleCoordination(barrier, refoldBarrier, runNum, vf, numMol, nBonds, bondInfo, timesteps)
#loadCoordVis(barrier, refoldBarrier, runNum, vf, numMol, particles, timesteps, parCoordination)

print("done :D")
