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
from modules.calcPorosity import *

# let's get the system data first
unfoldBarriers = [1, 2, 3, 4]
refoldBarrier = 2
numRuns = 5
vf = 0.04
numMol = 18007
boxLength = 100

for barrier in unfoldBarriers:
    conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}"
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
        # if not os.path.exists(f"../runs/{conditions}/{filename}/analysis/particleTraj.pkl"):
        #     data = readData(f"../runs/{conditions}/{filename}/output/dump.lammpstrj", Nsteps, Nwrite, Npar, equilTime)
        #     particles, timesteps = data[0], data[1]
        #     with open(f"../runs/{conditions}/{filename}/analysis/particleTraj.pkl", "wb") as f:
        #         pickle.dump(particles, f)
        #     with open(f"../runs/{conditions}/{filename}/analysis/timesteps.pkl", "wb") as f:
        #         pickle.dump(timesteps, f)
        # if os.path.exists(f"../runs/{conditions}/{filename}/analysis/particleTraj.pkl"):
        #     with open(f"../runs/{conditions}/{filename}/analysis/particleTraj.pkl", "rb") as f:
        #         particles = pickle.load(f)
        #     with open(f"../runs/{conditions}/{filename}/analysis/timesteps.pkl", "rb") as f:
        #         timesteps = pickle.load(f)

        # # check if bonds have been counted yet
        # if not os.path.exists(f"../runs/{conditions}/{filename}/analysis/nBonds.pkl"):
        #     nBonds = totalBonds(f"../runs/{conditions}/{filename}/output/nbonds.dat")
        #     with open(f"../runs/{conditions}/{filename}/analysis/nBonds.pkl", "wb") as f:
        #         pickle.dump(nBonds, f)
        # if os.path.exists(f"../runs/{conditions}/{filename}/analysis/nBonds.pkl"):
        #     with open(f"../runs/{conditions}/{filename}/analysis/nBonds.pkl", "rb") as f:
        #         nBonds = pickle.load(f)
        
        # # extract bond information e.g. bonded pairs, forces etc
        # if not os.path.exists(f"../runs/{conditions}/{filename}/analysis/bondInfo.pkl"):
        #     bondInfo = parseBondInfo(barrier, refoldBarrier, runNum, vf, numMol, nBonds, boxLength, timesteps)
        #     with open(f"../runs/{conditions}/{filename}/analysis/bondInfo.pkl", "wb") as f:
        #         pickle.dump(bondInfo, f)
        # if os.path.exists(f"../runs/{conditions}/{filename}/analysis/bondInfo.pkl"):
        #     with open(f"../runs/{conditions}/{filename}/analysis/bondInfo.pkl", "rb") as f:
        #         bondInfo = pickle.load(f)


        with open(f"../runs/{conditions}/{filename}/analysis/bondInfo.pkl", "rb") as f:
                bondInfo = pickle.load(f)
        with open(f"../runs/{conditions}/{filename}/analysis/nBonds.pkl", "rb") as f:
                nBonds = pickle.load(f)
        with open(f"../runs/{conditions}/{filename}/analysis/timesteps.pkl", "rb") as f:
                timesteps = pickle.load(f)
        with open(f"../runs/{conditions}/{filename}/analysis/particleTraj.pkl", "rb") as f:
                particles = pickle.load(f)


        # extract bonded pair information
        # if not os.path.exists(f"../runs/{conditions}/{filename}/analysis/percBonds.pkl"):    
        #     bondedPairs = parseForPercolation(particles,f"../runs/{conditions}/{filename}/output/bondinfo.dat",
        #                                        nBonds, boxLength, timesteps)
        #     percolatedBonds, bondedAtoms = bondedPairs[0], bondedPairs[1]
        #     with open(f"../runs/{conditions}/{filename}/analysis/percBonds.pkl", "wb") as f:
        #         pickle.dump(percolatedBonds, f)
        #     with open(f"../runs/{conditions}/{filename}/analysis/bondedPairs.pkl",  "wb") as f:
        #         pickle.dump(bondedPairs, f)
        # if os.path.exists(f"../runs/{conditions}/{filename}/analysis/percBonds.pkl"):
        #     with open(f"../runs/{conditions}/{filename}/analysis/percBonds.pkl", "rb") as f:
        #         percolatedBonds = pickle.load(f)
        #     with open(f"../runs/{conditions}/{filename}/analysis/bondedPairs.pkl", "rb") as f:
        #         bondedPairs = pickle.load(f)

        # if not os.path.exists(f"../runs/{conditions}/{filename}/analysis/percDims.pkl"):
        #     frameByFramePerc(particles, f"../runs/{conditions}/{filename}/output/bondinfo.dat",
        #                       f"../runs/{conditions}/{filename}/analysis/percinfo.txt",nBonds, boxLength, timesteps)
        #     systemDataFile = f"../runs/{conditions}/{filename}/output/systemData.txt"
        #     percInfoFile = f"../runs/{conditions}/{filename}/analysis/percinfo.txt"
        #     outputFile = f"../runs/{conditions}/{filename}/analysis/percDimsPerFrame.txt"
        #     subprocess.run(["./addingData/addingData", systemDataFile, percInfoFile, outputFile])
        #     percDims = getPercDims(barrier, refoldBarrier, runNum, vf, numMol, timesteps)
        #     with open(f"../runs/{conditions}/{filename}/analysis/percDims.pkl", "wb") as f:
        #         pickle.dump(percDims, f)
        
        # if os.path.exists(f"../runs/{conditions}/{filename}/analysis/percDims.pkl"):
        #     with open(f"../runs/{conditions}/{filename}/analysis/percDims.pkl", "rb") as f:
        #         percDims = pickle.load(f)
        

#avUnfoldPerFrame(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, timesteps)
#plotAvPressure(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, timesteps)
#plotAverageFractalDim(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, timesteps, percDims)
#anglePopulation(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, timesteps)
#avCoordination(unfoldBarriers, refoldBarrier, numRuns, numMol, vf, boxLength, timesteps)
avRDFs, avRDFsErr, shellBoundsavPosRDF = (unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, timesteps)
#calcPorosity(particles, barrier, refoldBarrier, runNum, vf, numMol, boxLength)
print("done :D")
