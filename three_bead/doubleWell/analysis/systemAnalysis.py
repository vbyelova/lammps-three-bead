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
numRuns = 5
vf = 0.04
numMol = 18007
boxLength = 100

with open(f"../runs/timesteps.pkl", "rb") as f:
    timesteps = pickle.load(f)

for barrier in unfoldBarriers:
    conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}"
    for runNum in range(numRuns):
        filename = f"Run{runNum}_{conditions}"
        systemData = parseSystemData(f"../runs/{conditions}/{filename}/output/systemData.txt")

        boxLength, Npar, Nsteps, Nwrite, equilTime = systemData[0], systemData[1], systemData[2], systemData[3], systemData[4]

        # # let's parse the dump file first and save it for later

        # make directories for data and figures
        # if not os.path.exists(f"../runs/{conditions}/{filename}/analysis/figs"):
        #     os.makedirs(f"../runs/{conditions}/{filename}/analysis/figs")
        # if not os.path.exists(f"../runs/{conditions}/averagedFigs"):
        #     os.makedirs(f"../runs/{conditions}/averagedFigs")
        # if not os.path.exists(f"../runs/{conditions}/data"):
        #     os.makedirs(f"../runs/{conditions}/data")
        # if not os.path.exists(f"../runs/boxLength{boxLength}"):
        #     os.makedirs(f"../runs/boxLength{boxLength}")
        # if not os.path.exists(f"../runs/boxLength{boxLength}/vf{vf}"):
        #     os.makedirs(f"../runs/boxLength{boxLength}/vf{vf}")
        # if not os.path.exists(f"../runs/boxLength{boxLength}/vf{vf}/data"):
        #     os.makedirs(f"../runs/boxLength{boxLength}/vf{vf}/data")
        #     with open(f"../runs/{conditions}/analysisNotes.txt", "w") as f:
        #         f.write(f"# any extra info e.g. runs skipped will be noted here.\n")
        #         f.close()


        # check if particle trajectories have been parsed
        # if not os.path.exists(f"../runs/{conditions}/{filename}/analysis/particleTraj.pkl"):
        #     data = readData(f"../runs/{conditions}/{filename}/output/dump.lammpstrj", Nsteps, Nwrite, Npar, equilTime)
        #     particles = data[0], data[1]
        #     with open(f"../runs/{conditions}/{filename}/analysis/particleTraj.pkl", "wb") as f:
        #         pickle.dump(particles, f)

        # if not os.path.exists(f"../runs/{conditions}/timesteps.pkl"):
        #     with open(f"../runs/{conditions}/timesteps.pkl", "wb") as f:
        #         pickle.dump(timesteps, f)
        # if os.path.exists(f"../runs/{conditions}/{filename}/analysis/particleTraj.pkl"):
        #     with open(f"../runs/{conditions}/{filename}/analysis/particleTraj.pkl", "rb") as f:
        #         particles = pickle.load(f)
        #     with open(f"../runs/{conditions}/{filename}/analysis/timesteps.pkl", "rb") as f:
        #         timesteps = pickle.load(f)

        # check if bonds have been counted yet
        # if not os.path.exists(f"../runs/{conditions}/{filename}/analysis/nBonds.pkl"):
        #     nBonds = totalBonds(f"../runs/{conditions}/{filename}/output/nbonds.dat")
        #     with open(f"../runs/{conditions}/{filename}/analysis/nBonds.pkl", "wb") as f:
        #         pickle.dump(nBonds, f)
        # if os.path.exists(f"../runs/{conditions}/{filename}/analysis/nBonds.pkl"):
        #     with open(f"../runs/{conditions}/{filename}/analysis/nBonds.pkl", "rb") as f:
        #         nBonds = pickle.load(f)
        
        # # extract bond information e.g. bonded pairs, forces etc
        # if not os.path.exists(f"../runs/{conditions}/{filename}/analysis/bondInfo.pkl"):
        #     bondInfo = parseBondInfo(barrier, refoldBarrier, runNum, vf, numMol, nBonds, boxLength)
        #     with open(f"../runs/{conditions}/{filename}/analysis/bondInfo.pkl", "wb") as f:
        #         pickle.dump(bondInfo, f)
        # if os.path.exists(f"../runs/{conditions}/{filename}/analysis/bondInfo.pkl"):
        #     with open(f"../runs/{conditions}/{filename}/analysis/bondInfo.pkl", "rb") as f:
        #         bondInfo = pickle.load(f)


        # with open(f"../runs/{conditions}/{filename}/analysis/bondInfo.pkl", "rb") as f:
        #         bondInfo = pickle.load(f)
        # with open(f"../runs/{conditions}/{filename}/analysis/nBonds.pkl", "rb") as f:
        #         nBonds = pickle.load(f)
        # with open(f"../runs/{conditions}/{filename}/analysis/timesteps.pkl", "rb") as f:
        #         timesteps = pickle.load(f)
        # with open(f"../runs/{conditions}/{filename}/analysis/particleTraj.pkl", "rb") as f:
        #         particles = pickle.load(f)


        # extract bonded pair information
        # if not os.path.exists(f"../runs/{conditions}/{filename}/analysis/percBonds.pkl"):    
        #     bondedPairs = parseForPercolation(particles,f"../runs/{conditions}/{filename}/output/bondinfo.dat",
        #                                        nBonds, boxLength)
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
        #     frameByFramePerc(barrier, refoldBarrier, runNum, vf, numMol, boxLength)
        #     systemDataFile = f"../runs/{conditions}/{filename}/output/systemData.txt"
        #     percInfoFile = f"../runs/{conditions}/{filename}/analysis/percinfo.txt"
        #     outputFile = f"../runs/{conditions}/{filename}/analysis/percDimsPerFrame.txt"
        #     if not os.path.exists(f"../runs/{conditions}/{filename}/analysis/percDimsPerFrame.txt"):
        #         subprocess.run(["./addingData/addingData", systemDataFile, percInfoFile, outputFile])
        #     percDims = getPercDims(barrier, refoldBarrier, runNum, vf, numMol)
        #     with open(f"../runs/{conditions}/{filename}/analysis/percDims.pkl", "wb") as f:
        #         pickle.dump(percDims, f)
        
#         if os.path.exists(f"../runs/{conditions}/{filename}/analysis/percDims.pkl"):
#             with open(f"../runs/{conditions}/{filename}/analysis/percDims.pkl", "rb") as f:
#                 percDims = pickle.load(f)
        
        # if not os.path.exists(f"../runs/{conditions}/{filename}/analysis/mainPores.pkl"):
        #      mainPores = calcPorosity(barrier, refoldBarrier, runNum, vf, numMol, boxLength)
        #      with open(f"../runs/{conditions}/{filename}/analysis/mainPores.pkl", "wb") as f:
        #           pickle.dump(mainPores, f)
        
        
# if not os.path.exists(f"../runs/boxLength{boxLength}/data/atPercFractalDims.pkl"):
#     atPercFinalFractalDims, atPercFinalFractalDimsError, atPercCorrLength, atPercCorrLengthErr = (
#                                                                             findAllFractalDims(unfoldBarriers, refoldBarrier,
#                                                                             numRuns, vf, numMol, boxLength,
#                                                                             timesteps, True))
#     with open(f"../runs/boxLength{boxLength}/data/atPercFractalDims.pkl", "wb") as f:
#          pickle.dump((atPercFinalFractalDims, atPercFinalFractalDimsError, atPercCorrLength, atPercCorrLengthErr), f)

#     simEndFinalFractalDims, simEndFinalFractalDimsError, simEndCorrLength, simEndCorrLengthErr = (
#                                                                             findAllFractalDims(unfoldBarriers, refoldBarrier,
#                                                                             numRuns, vf, numMol, boxLength,
#                                                                             timesteps, False))
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
#         pickle.dump((avPressure, avPressureErr))

# if not os.path.exists(f"../runs/boxLength{boxLength}/vf{vf}/data/avCoordination/pkl"):
#     avCoordination, avCoordinationErr = plotAvCoordination(unfoldBarriers, refoldBarrier, numRuns, numMol, vf, boxLength)
#     with open(f"../runs/boxLength{boxLength}/vf{vf}/data/avCoordination/pkl", "wb") as f:
#         pickle.dump((avCoordination, avCoordinationErr), f)

def plotTogether(unfoldBarriers, refoldBarrier, numRuns, numMol, vf, boxLength):
    with open(f"../runs/boxLength{boxLength}/vf{vf}/data/avNewUnfoldedMol.pkl", "rb") as f:
        avNewUnfoldedMol, avNewUnfoldedMolErr = pickle.load(f)
    with open(f"../runs/boxLength{boxLength}/vf{vf}/data/avPercolation.pkl", "rb") as f:
        avPercolation, avPercolationErr = pickle.load(f)
    # with open(f"../runs/boxLength{boxLength}/vf{vf}/data/avPressure.pkl", "rb") as f:
    #     avPressure, avPressureErr = pickle.load(f)
    with open(f"../runs/boxLength{boxLength}/vf{vf}/data/unfoldedOverTime.pkl", "rb") as f:
        avUnfoldOverTime, avUnfoldOverTimeErr = pickle.load(f)
    
    fig, (ax1, ax2, ax3) = plt.subplots(3, sharex = True, figsize = (15, 15))
    fig.subplots_adjust(wspace = 0, hspace = 0)
    for barrier in unfoldBarriers:
        ax1.errorbar(timesteps, avUnfoldOverTime[barrier], yerr = avUnfoldOverTimeErr[barrier],
                     label = f"barrier = {barrier}kT")
        ax2.errorbar(timesteps, avNewUnfoldedMol[barrier], yerr = avNewUnfoldedMolErr[barrier],
                label = f"unfolding barrier = {barrier}kT")
        ax3.errorbar(timesteps, avPercolation[barrier], yerr = avPercolationErr[barrier],
                    label = f"barrier = {barrier}kT")
        # ax4.errorbar(timesteps, avPressure[barrier], yerr = avPressureErr[barrier],
        #             label = f"barrier = {barrier}kT")
        
    ax1.set_ylabel("number of unfolded molecules")
    ax2.set_ylabel("new unfolded molecules")
    ax3.set_ylabel("percolation dimension")
    ax1.legend()
    ax2.legend()
    ax3.legend()
#    ax4.set_ylabel("average pressure")
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/sharedaxis_vf{vf}.png")
    ax1.semilogx()
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/sharedaxis_vf{vf}_semilog.png")
    plt.close()

plotTogether(unfoldBarriers, refoldBarrier, numRuns, numMol ,vf, boxLength)
#anglePopulation(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength)
#plotAvCoordination(unfoldBarriers, refoldBarrier, numRuns, numMol, vf, boxLength)
print("done :D")
