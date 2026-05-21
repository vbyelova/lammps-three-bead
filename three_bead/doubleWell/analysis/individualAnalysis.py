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
unfoldBarriers = [1, 2, 3, 4, 5]
refoldBarrier = 2
numRuns = 10
vf = 0.04
numMol = 18007
boxLength = 100

# config = configparser.ConfigParser()
# config.read("analysisParams.txt")

# barrier = int(config.get("analysisParams", "barrier"))
# refoldBarrier = int(config.get("analysisParams", "refoldBarrier"))
# runNum = int(config.get("analysisParams", "runNum"))
# vf = float(config.get("analysisParams", "vf"))
# numMol = int(config.get("analysisParams", "numMol"))
# boxLength = int(config.get("analysisParams", "boxLength"))

conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}"
filename = f"Run{runNum}_{conditions}"
systemData = parseSystemData(f"../runs/{conditions}/{filename}/output/systemData.txt")

boxLength, nPar, nSteps, nWrite, equilTime = systemData[0], systemData[1], systemData[2], systemData[3], systemData[4]

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
if not os.path.exists(f"../runs/boxLength{boxLength}/vf{vf}"):
    os.makedirs(f"../runs/boxLength{boxLength}/vf{vf}")
if not os.path.exists(f"../runs/boxLength{boxLength}/vf{vf}/data"):
    os.makedirs(f"../runs/boxLength{boxLength}/vf{vf}/data")
    with open(f"../runs/{conditions}/analysisNotes.txt", "w") as f:
        f.write(f"# any extra info e.g. runs skipped will be noted here.\n")
        f.close()


# check if particle trajectories have been parsed
if not os.path.exists(f"../runs/{conditions}/{filename}/analysis/particleTraj.pkl"):
    data = readData(f"../runs/{conditions}/{filename}/output/dump.lammpstrj", nSteps, nWrite, nPar, equilTime)
    particles, timesteps = data[0], data[1]
    with open(f"../runs/{conditions}/{filename}/analysis/particleTraj.pkl", "wb") as f:
        pickle.dump(particles, f)

    if not os.path.exists(f"../runs/{conditions}/timesteps.pkl"):
        with open(f"../runs/{conditions}/timesteps.pkl", "wb") as f:
            pickle.dump(timesteps, f)


# check if bonds have been counted yet
if not os.path.exists(f"../runs/{conditions}/{filename}/analysis/nBonds.pkl"):
    nBonds = totalBonds(f"../runs/{conditions}/{filename}/output/nbonds.dat")
    with open(f"../runs/{conditions}/{filename}/analysis/nBonds.pkl", "wb") as f:
        pickle.dump(nBonds, f)

# extract bond information e.g. bonded pairs, forces etc
if not os.path.exists(f"../runs/{conditions}/{filename}/analysis/bondInfo.pkl"):
    bondInfo = parseBondInfo(barrier, refoldBarrier, runNum, vf, numMol, nBonds, boxLength)
    with open(f"../runs/{conditions}/{filename}/analysis/bondInfo.pkl", "wb") as f:
        pickle.dump(bondInfo, f)

if not os.path.exists(f"../runs/{conditions}/{filename}/analysis/percDimsPerFrame.txt"):
    frameByFramePerc(barrier, refoldBarrier, runNum, vf, numMol, boxLength)
    systemDataFile = f"../runs/{conditions}/{filename}/output/systemData.txt"
    percInfoFile = f"../runs/{conditions}/{filename}/analysis/percinfo.txt"
    outputFile = f"../runs/{conditions}/{filename}/analysis/percDimsPerFrame.txt"
    subprocess.run(["./addingData/addingData", systemDataFile, percInfoFile, outputFile])
    percDims = getPercDims(barrier, refoldBarrier, runNum, vf, numMol)
    with open(f"../runs/{conditions}/{filename}/analysis/percDims.pkl", "wb") as f:
        pickle.dump(percDims, f)
        
print("done :D")
