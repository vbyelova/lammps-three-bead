import numpy as np
import matplotlib.pyplot as plt
import pickle
from collections import defaultdict

from .calcAngles import *

def densityUnfoldedCorrelation(barrier, refoldBarrier, runNum , vf, numMol, boxLength, suffix):
    """correlation function between local density and degree of unfolding """
    voxelSize = 0.07 * boxLength
    voxelUnfoldedCount = defaultdict(int)
    voxelParticleCount = defaultdict(int)
    moleculeAngles = defaultdict(list)
    conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
    filename = f"Run{runNum}_{conditions}"

    with open(f"../runs/{conditions}/{filename}/analysis/particleTraj.pkl", "rb") as f:
        particles = pickle.load(f)
    with open(f"../runs/{conditions}/timesteps.pkl", "rb") as f:
        timesteps = pickle.load(f)
    # first find local density
    # split box into voxels
    bins = np.arange(-0.5 * boxLength, 0.5 * boxLength, voxelSize) 
    
    gx, gy, gz = np.meshgrid(bins, bins, bins, indexing = 'ij')
    coords = np.column_stack([gx.ravel(), gy.ravel(), gz.ravel()])

    
    angles = calcAngles(barrier, refoldBarrier, runNum, vf, numMol, suffix)
    frame = max(angles.keys())

    # count how many particles per voxel (discrete)
    for pNum, p in enumerate(particles):
        x = p.properties[frame, 1]
        y = p.properties[frame, 2]
        z = p.properties[frame, 3]

        xVox = np.digitize(x, bins) - 1   #digitize returns 1-based, so -1
        yVox = np.digitize(y, bins) - 1
        zVox = np.digitize(z, bins) - 1

        voxelIndex = ((xVox, yVox, zVox))
        voxelParticleCount[voxelIndex] += 1
    # next let's find the degree of unfolding

        # check for the central particle
        if pNum % 3 == 1:
            molNum = pNum // 3
            # for each molecule in the box, get the angle of unfolding
            centralAngle = angles[frame][molNum]
            if centralAngle >= 120:
                voxelUnfoldedCount[voxelIndex] += 1
            moleculeAngles[voxelIndex].append(centralAngle)

    allCounts = []
    allMeanAngles = []
    voxelData = {}

    for voxelID in voxelParticleCount.keys():
        if voxelID in moleculeAngles and len(moleculeAngles[voxelID]) > 0:
            localMeanAngle = np.mean(moleculeAngles[voxelID])

            allCounts.append(voxelParticleCount[voxelID])
            allMeanAngles.append(localMeanAngle)

            voxelData[voxelID]= { "meanAngle": localMeanAngle,
                                 "totalPar": voxelParticleCount[voxelID],
                                 "foldedPar": voxelParticleCount[voxelID] - voxelUnfoldedCount[voxelID],
                                 "unfoldedPar": voxelUnfoldedCount[voxelID]}
            
    meanCount = np.mean(allCounts)
    meanAngle = np.mean(allMeanAngles)

    corrFunc = []
    for voxelID, data in voxelData.items():
        corr = ((data["totalPar"] - meanCount) / numMol) * (data["meanAngle"] - meanAngle)
        corrFunc.append(corr)

    #print(allCounts, allMeanAngles)
    plt.scatter(allCounts, allMeanAngles)
    plt.ylabel("molecule angle")
    plt.xlabel("particle density")
    meanCorr = np.mean(corrFunc)
    print(meanCorr)
    plt.savefig(f"../runs/{conditions}/{filename}/analysis/figs/densityUnfoldingCorr.png")
    plt.clf()
    #plt.show()
    return allCounts, allMeanAngles, meanCorr


def calcAvDensityUnfoldedCorr(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, suffixes):
    """takes an average of all the correlation functions calculated"""
    allCorrFuncs = defaultdict(list)
    barrierAllCounts = defaultdict(list)
    barrierAllMeanAngles = defaultdict(list)
    avAllCounts = {}
    avAllCountsErr = {}
    avAllMeanAngles = {}
    avAllMeanAnglesErr = {}
    avCorrFuncs = {}
    avCorrFuncsErr = {}

    for barrier, suffix in zip(unfoldBarriers, suffixes):
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}"
        with open(f"../runs/{conditions}/timesteps.pkl", "rb") as f:
            timesteps = pickle.load(f)

        for runNum in range(numRuns):
            filename = f"Run{runNum}_{conditions}"
            with open(f"../runs/{conditions}/{filename}/analysis/densityUnfoldedCorr.pkl", "rb") as f:
                allCounts, allMeanAngles, meanCorr = pickle.load(f)
            barrierAllCounts[barrier].append(allCounts)
            barrierAllMeanAngles[barrier].append(allMeanAngles)
            allCorrFuncs[barrier].append(meanCorr)
        
        avCorrFuncsArray = np.array(allCorrFuncs[barrier])
        avCorrFuncs[barrier] = np.mean(avCorrFuncsArray, axis = 0)
        avCorrFuncsErr[barrier] = np.std(avCorrFuncsArray, axis = 0)
    
        avAllCountsArray = np.array(barrierAllCounts[barrier])
        avAllCounts[barrier] = np.mean(avAllCountsArray, axis = 0)
        avAllCounts[barrier] = np.std(avAllCountsArray, axis = 0)

        avAllMeanAnglesArray = np.array(barrierAllMeanAngles[barrier])
        avAllMeanAngles[barrier] = np.mean(barrierAllMeanAngles, axis = 0)
        avAllMeanAnglesErr[barrier] = np.std(barrierAllMeanAngles, axis = 0)

    return avAllCounts, avAllCountsErr, avAllMeanAngles, avAllMeanAnglesErr, avCorrFuncs, avCorrFuncsErr

def plotAvDensityUnfoldedCorr(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, suffixes):
    with open(f"../runs/boxLength{boxLength}/vf{vf}/data/densityUnfoldedCorr.pkl", "rb") as f:
        avAllCounts, avAllCountsErr, avAllMeanAngles, avAllMeanAnglesErr, avCorrFuncs, avCorrFuncsErr = pickle.load(f)

    fig, ax = plt.subplots(2, 3, sharex = True, sharey = True)

    ax[0, 0].scatter(avAllCounts[unfoldBarriers[0]], avAllMeanAngles[unfoldBarriers[0]], label = f"unfolding barrier = {unfoldBarriers[0]}kT")
    ax[1, 0].scatter(avAllCounts[unfoldBarriers[1]], avAllMeanAngles[unfoldBarriers[1]], label = f"unfolding barrier = {unfoldBarriers[1]}kT")
    ax[0, 1].scatter(avAllCounts[unfoldBarriers[2]], avAllMeanAngles[unfoldBarriers[2]], label = f"unfolding barrier = {unfoldBarriers[2]}kT")
    ax[0, 2].scatter(avAllCounts[unfoldBarriers[3]], avAllMeanAngles[unfoldBarriers[3]], label = f"unfolding barrier = {unfoldBarriers[3]}kT")
    ax[1, 1].scatter(avAllCounts[unfoldBarriers[4]], avAllMeanAngles[unfoldBarriers[4]], label = f"unfolding barrier = {unfoldBarriers[4]}kT")
#    ax6.scatter(avCorrFuncs[5][0], avCorrFuncs[5][1], label = f"unfolding barrier = {unfoldBarriers[5]}kT")

    ax.supylabel("unfolding degree")
    ax.supxlabel("density of voxel")

    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/avDensityUnfoldedCorr.png")
    plt.show()
    plt.close()
    return
