import numpy as np
import pickle
import re
import os
import subprocess
from collections import defaultdict
import multiprocessing as mp

from modules.parseDump import *
from modules.calcAngles import *
from modules.calcFracDim import * 
from modules.calcPercolation import *
from modules.calcStressTensor import *
from modules.calcCorrelation import *
from modules.threeBeadClasses import *
from modules.parseRDF import *
from modules.calcPorosity import *
from modules.parallelPercolation import *

# let's get the system data first
unfoldBarriers = [1, 2, 3, 4, 5]#[1, 1, 2, 2, 3, 3, 4, 4, 5, 5]
suffixes = ["", "", "", "", ""]#["", "NOPERCOLATIONNOLJ", "", "NOPERCOLATIONNOLJ", "", "NOPERCOLATIONNOLJ", "", "NOPERCOLATIONNOLJ", "", "NOPERCOLATIONNOLJ"]
refoldBarrier = 2
numRuns = 10
vf = 0.07
numMol = 31512
boxLength = 100
bondsPerAtom = 2
prob = 1

for barrier, suffix in zip(unfoldBarriers, suffixes):
    if bondsPerAtom == 2:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
    else:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
        if prob < 1:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}_prob{prob}{suffix}"

    for runNum in range(numRuns):
        filename = f"Run{runNum}_{conditions}"
        systemData = parseSystemData(f"../runs/{conditions}/{filename}/output/systemData.txt")

        boxLength, nPar, nSteps, nWrite, equilTime = systemData[0], systemData[1], systemData[2], systemData[3], systemData[4]
        #densityUnfoldingDegreeCorrelation(barrier, refoldBarrier, vf, runNum, numMol, boxLength)
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
        if not os.path.exists(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}"):
            os.makedirs(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}")
        if not os.path.exists(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}"):
            os.makedirs(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}")
        if not os.path.exists(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/data"):
            os.makedirs(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/data")
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

        if not os.path.exists(f"../runs/{conditions}/{filename}/analysis/coordination.pkl"):
            coord = coordination(barrier, refoldBarrier, runNum, vf, numMol, bondsPerAtom, prob, suffix)
            with open(f"../runs/{conditions}/{filename}/analysis/coordination.pkl", "wb") as f:
                pickle.dump(coord, f)

        # check if bonds have been counted yet
        if not os.path.exists(f"../runs/{conditions}/{filename}/analysis/nBonds.pkl"):
            nBonds = totalBonds(f"../runs/{conditions}/{filename}/output/nbonds.dat")
            with open(f"../runs/{conditions}/{filename}/analysis/nBonds.pkl", "wb") as f:
                pickle.dump(nBonds, f)
        
        # extract bond information e.g. bonded pairs, forces etc
        if not os.path.exists(f"../runs/{conditions}/{filename}/analysis/bondInfo.pkl"):
            bondInfo = parseBondInfo(barrier, refoldBarrier, runNum, vf, numMol, boxLength, bondsPerAtom, prob, suffix)
            with open(f"../runs/{conditions}/{filename}/analysis/bondInfo.pkl", "wb") as f:
                pickle.dump(bondInfo, f)

        # if not os.path.exists(f"../runs/{conditions}/{filename}/analysis/densityUnfoldedCorr.pkl"):
        #     meanCorr = densityUnfoldedCorrelation(barrier, refoldBarrier, runNum, vf, numMol, boxLength, bondsPerAtom, prob, suffix)
        #     with open(f"../runs/{conditions}/{filename}/analysis/densityUnfoldedCorr.pkl", "wb") as f:
        #         pickle.dump(meanCorr, f)

        if not os.path.exists(f"../runs/{conditions}/{filename}/analysis/angleDist.pkl"):
            hist, bins, finalFrameAngles = bimodalAngle(barrier, refoldBarrier, runNum, vf, numMol, bondsPerAtom, prob, suffix)
            with open(f"../runs/{conditions}/{filename}/analysis/angleDist.pkl", "wb") as f:
                pickle.dump((hist, bins, finalFrameAngles), f)

        # percDims = getPercDims(barrier, refoldBarrier, runNum, vf, numMol, bondsPerAtom, prob, suffix)
        # if not os.path.exists(f"../runs/{conditions}/{filename}/analysis/percDims.pkl"):
        #     with open(f"../runs/{conditions}/{filename}/analysis/percDims.pkl", "wb") as f:
        #         pickle.dump(percDims, f)
        
        
# if not os.path.exists(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/data/atPercFractalDims.pkl"):
#     atPercFinalFractalDims, atPercFinalFractalDimsError, atPercCorrLength, atPercCorrLengthErr = (
#                                                                             findAllFractalDims(unfoldBarriers, refoldBarrier,
#                                                                             numRuns, vf, numMol, boxLength,
#                                                                             True, bondsPerAtom, prob, suffixes))
#     with open(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/data/atPercFractalDims.pkl", "wb") as f:
#          pickle.dump((atPercFinalFractalDims, atPercFinalFractalDimsError, atPercCorrLength, atPercCorrLengthErr), f)

#     simEndFinalFractalDims, simEndFinalFractalDimsError, simEndCorrLength, simEndCorrLengthErr = (
#                                                                             findAllFractalDims(unfoldBarriers, refoldBarrier,
#                                                                             numRuns, vf, numMol, boxLength,
#                                                                             False, bondsPerAtom, prob, suffixes))
#     with open(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/data/simEndFractalDims.pkl", "wb") as f:
#          pickle.dump((simEndFinalFractalDims, simEndFinalFractalDimsError, simEndCorrLength, simEndCorrLengthErr), f)

# if not os.path.exists(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/data/avPercolation.pkl"):
#     avPercolation, avPercolationErr = calcAvPercolation(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, prob, suffix)
#     with open(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/data/avPercolation.pkl", "wb") as f:
#             pickle.dump((avPercolation, avPercolationErr), f)

# perc_args = [
#     (barrier, suffix, refoldBarrier, vf, numMol, boxLength, bondsPerAtom, prob, runNum)
#     for barrier, suffix in zip(unfoldBarriers, suffixes)
#     for runNum in range(numRuns)
#     if suffix not in ("NOPERCOLATIONNOLJ", "NOPERCOLATION")   # skip non-percolation runs early
# ]

# # Use all available cores, or cap it (e.g. mp.cpu_count() - 1)
# with mp.Pool(processes=mp.cpu_count()) as pool:
#     pool.map(process_single_run, perc_args)

# if not os.path.exists(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/data/avRDF.pkl"):
#     avRDFs, avRDFsErr, shellBounds = calcPosRDF(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, prob, suffixes)
#     with open(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/data/avRDF.pkl", "wb") as f:
#         pickle.dump((avRDFs, avRDFsErr, shellBounds), f)
#     plotAvRDF(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, prob, suffixes)

# if not os.path.exists(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/data/avNewUnfoldedMol.pkl"):
#     avNewUnfoldedMol, avNewUnfoldedMolErr = avUnfoldPerFrame(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, prob, suffixes)
#     with open(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/data/avNewUnfoldedMol.pkl", "wb") as f:
#         pickle.dump((avNewUnfoldedMol, avNewUnfoldedMolErr), f)

# if not os.path.exists(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/data/unfoldedOverTime.pkl"):
#     avUnfoldOverTime, avUnfoldOverTimeErr = unfoldedOverTime(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, prob, suffixes)
#     with open(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/data/unfoldedOverTime.pkl", "wb") as f:
#         pickle.dump((avUnfoldOverTime, avUnfoldOverTimeErr), f)

# if not os.path.exists(f"../runs/boxLength{boxLength}/vf{vf}/data/avPressure.pkl"):
#     avPressure, avPressureErr = calcAvPressure(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, prob, suffixes)
#     with open(f"../runs/boxLength{boxLength}/vf{vf}/data/avPressure.pkl", "wb") as f:
#         pickle.dump((avPressure, avPressureErr), f)

# if not os.path.exists(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/data/avCoordination.pkl"):
#     avCoordination, avCoordinationErr = plotAvCoordination(unfoldBarriers, refoldBarrier, numRuns, numMol, vf, boxLength, bondsPerAtom, prob, suffixes)
#     with open(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/data/avCoordination.pkl", "wb") as f:
#         pickle.dump((avCoordination, avCoordinationErr), f)



#plotAvBimodalAngle(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, prob, suffixes)
#collapseUnfoldedOverTime(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, prob, suffixes)
# avAllCounts, avAllCountsErr, avAllMeanAngles, avAllMeanAnglesErr, avCorrFuncs, avCorrFuncsErr = calcAvDensityUnfoldedCorr(unfoldBarriers, refoldBarrier,
#                                                                                                                            numRuns, vf, numMol, boxLength, bondsPerAtom, prob, suffixes)
# with open(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/data/densityUnfoldedCorr.pkl", "wb") as f:
#    pickle.dump((avAllCounts, avAllCountsErr, avAllMeanAngles, avAllMeanAnglesErr, avCorrFuncs, avCorrFuncsErr), f)
#plotAvDensityUnfoldedCorr(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, prob, suffixes)
#plotAverageFractalDim(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, prob)
#avBoxCounting(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, prob, suffixes)
#plotTogetherShifted(unfoldBarriers, refoldBarrier, numRuns, numMol, vf, boxLength)
#plotTogether(unfoldBarriers, refoldBarrier, numRuns, numMol ,vf, boxLength)
#anglePopulation(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, prob, suffixes)
#plotAvCoordination(unfoldBarriers, refoldBarrier, numRuns, numMol, vf, boxLength, bondsPerAtom, prob, suffixes)
#plotBimodalAndDensityUnfoldedCorr(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, prob, suffixes)
#plotUnfoldedOverTime(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, prob, suffixes)
#plotAvPercolation(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, suffixes)
#plotTotalInterMol(unfoldBarriers, refoldBarrier, numRuns, numMol, vf, boxLength, bondsPerAtom)
#plotAvPressure(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, prob, suffixes)
#plotFolded(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, prob, suffixes)
#plotTwoVfFractalDim(unfoldBarriers, refoldBarrier, numRuns, [0.07, 0.07], [31512, 3939], [100, 50], bondsPerAtom, prob)
plotCentralCoordination(unfoldBarriers, refoldBarrier, numRuns, numMol, vf, boxLength, bondsPerAtom, prob, suffixes)
print("done :D")
