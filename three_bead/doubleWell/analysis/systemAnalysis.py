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
unfoldBarriers = [5]
refoldBarrier = 2
numRuns = 10
vf = 0.04
numMol = 18007
boxLength = 100


for barrier in unfoldBarriers:
    conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}"
    for runNum in range(numRuns):
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


        # extract bonded pair information
        # if not os.path.exists(f"../runs/{conditions}/{filename}/analysis/percBonds.pkl"):    
        #     bondedPairs = parseForPercolation(particles,f"../runs/{conditions}/{filename}/output/bondinfo.dat",
        #                                        nBonds, boxLength)
        #     percolatedBonds, bondedAtoms = bondedPairs[0], bondedPairs[1]
        #     with open(f"../runs/{conditions}/{filename}/analysis/percBonds.pkl", "wb") as f:
        #         pickle.dump(percolatedBonds, f)
        #     with open(f"../runs/{conditions}/{filename}/analysis/bondedPairs.pkl",  "wb") as f:
        #         pickle.dump(bondedPairs, f)


        if not os.path.exists(f"../runs/{conditions}/{filename}/analysis/percDims.pkl"):
            frameByFramePerc(barrier, refoldBarrier, runNum, vf, numMol, boxLength)
            systemDataFile = f"../runs/{conditions}/{filename}/output/systemData.txt"
            percInfoFile = f"../runs/{conditions}/{filename}/analysis/percinfo.txt"
            outputFile = f"../runs/{conditions}/{filename}/analysis/percDimsPerFrame.txt"
            if not os.path.exists(f"../runs/{conditions}/{filename}/analysis/percDimsPerFrame.txt"):
                subprocess.run(["./addingData/addingData", systemDataFile, percInfoFile, outputFile])
            percDims = getPercDims(barrier, refoldBarrier, runNum, vf, numMol)
            with open(f"../runs/{conditions}/{filename}/analysis/percDims.pkl", "wb") as f:
                pickle.dump(percDims, f)
        
        
        # if not os.path.exists(f"../runs/{conditions}/{filename}/analysis/mainPores.pkl"):
        #      mainPores = calcPorosity(barrier, refoldBarrier, runNum, vf, numMol, boxLength)
        #      with open(f"../runs/{conditions}/{filename}/analysis/mainPores.pkl", "wb") as f:
        #           pickle.dump(mainPores, f)
        
        
# if not os.path.exists(f"../runs/boxLength{boxLength}/data/atPercFractalDims.pkl"):
#     atPercFinalFractalDims, atPercFinalFractalDimsError, atPercCorrLength, atPercCorrLengthErr = (
#                                                                             findAllFractalDims(unfoldBarriers, refoldBarrier,
#                                                                             numRuns, vf, numMol, boxLength,
#                                                                             True))
#     with open(f"../runs/boxLength{boxLength}/data/atPercFractalDims.pkl", "wb") as f:
#          pickle.dump((atPercFinalFractalDims, atPercFinalFractalDimsError, atPercCorrLength, atPercCorrLengthErr), f)

#     simEndFinalFractalDims, simEndFinalFractalDimsError, simEndCorrLength, simEndCorrLengthErr = (
#                                                                             findAllFractalDims(unfoldBarriers, refoldBarrier,
#                                                                             numRuns, vf, numMol, boxLength,
#                                                                             False))
#     with open(f"../runs/boxLength{boxLength}/data/simEndFractalDims.pkl", "wb") as f:
#          pickle.dump((simEndFinalFractalDims, simEndFinalFractalDimsError, simEndCorrLength, simEndCorrLengthErr), f)

# if not os.path.exists(f"../runs/boxLength{boxLength}/vf{vf}/data/avPercolation.pkl"):
#     avPercolation, avPercolationErr = plotAvPercolation(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength)
#     with open(f"../runs/boxLength{boxLength}/vf{vf}/data/avPercolation.pkl", "wb") as f:
#             pickle.dump((avPercolation, avPercolationErr), f)

# if not os.path.exists(f"../runs/boxLength{boxLength}/vf{vf}/data/avRDF.pkl"):
#     avRDFs, avRDFsErr, shellBounds = avPosRDF(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength)
#     with open(f"../runs/boxLength{boxLength}/vf{vf}/data/avRDF.pkl", "wb") as f:
#         pickle.dump((avRDFs, avRDFsErr, shellBounds), f)

# if not os.path.exists(f"../runs/boxLength{boxLength}/vf{vf}/data/avNewUnfoldedMol.pkl"):
#     avNewUnfoldedMol, avNewUnfoldedMolErr = avUnfoldPerFrame(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength)
#     with open(f"../runs/boxLength{boxLength}/vf{vf}/data/avNewUnfoldedMol.pkl", "wb") as f:
#         pickle.dump((avNewUnfoldedMol, avNewUnfoldedMolErr), f)

# if not os.path.exists(f"../runs/boxLength{boxLength}/vf{vf}/data/unfoldedOverTime.pkl"):
#     avUnfoldOverTime, avUnfoldOverTimeErr = unfoldedOverTime(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength)
#     with open(f"../runs/boxLength{boxLength}/vf{vf}/data/unfoldedOverTime.pkl", "wb") as f:
#         pickle.dump((avUnfoldOverTime, avUnfoldOverTimeErr), f)

# if not os.path.exists(f"../runs/boxLength{boxLength}/vf{vf}/data/avPressure.pkl"):
#     avPressure, avPressureErr = plotAvPressure(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength)
#     with open(f"../runs/boxLength{boxLength}/vf{vf}/data/avPressure.pkl", "wb") as f:
#         pickle.dump((avPressure, avPressureErr), f)

# if not os.path.exists(f"../runs/boxLength{boxLength}/vf{vf}/data/avCoordination/pkl"):
#     avCoordination, avCoordinationErr = plotAvCoordination(unfoldBarriers, refoldBarrier, numRuns, numMol, vf, boxLength)
#     with open(f"../runs/boxLength{boxLength}/vf{vf}/data/avCoordination/pkl", "wb") as f:
#         pickle.dump((avCoordination, avCoordinationErr), f)

#plotTogetherShifted(unfoldBarriers, refoldBarrier, numRuns, numMol, vf, boxLength)
#plotTogether(unfoldBarriers, refoldBarrier, numRuns, numMol ,vf, boxLength)
#anglePopulation(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength)
#plotAvCoordination(unfoldBarriers, refoldBarrier, numRuns, numMol, vf, boxLength)
print("done :D")
