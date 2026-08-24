import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pickle
from collections import defaultdict
import pandas as pd

from .calcAngles import *

def densityUnfoldedCorrelation(barrier, refoldBarrier, runNum , vf, numMol, boxLength, bondsPerAtom, prob, suffix):
    """correlation function between local density and degree of unfolding """
    voxelSize = 0.07 * boxLength
    voxelUnfoldedCount = defaultdict(int)
    voxelParticleCount = defaultdict(int)
    moleculeAngles = defaultdict(list)
    if bondsPerAtom == 2:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
    else:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
        if prob < 1:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}_prob{prob}{suffix}"
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

    
    angles = calcAngles(barrier, refoldBarrier, runNum, vf, numMol, bondsPerAtom, prob, suffix)
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
    allCountsIncEmpty = []
    allMeanAnglesIncEmpty = []
    voxelData = {}

    allVoxelIDs = [(x, y, z) for x in range(len(bins)) for y in range(len(bins)) for z in range(len(bins))]
    for voxelID in allVoxelIDs:
        if voxelID in moleculeAngles and len(moleculeAngles[voxelID]) > 0:
            localMeanAngle = np.mean(moleculeAngles[voxelID])

            allCounts.append(voxelParticleCount[voxelID])
            allMeanAngles.append(localMeanAngle)
            allCountsIncEmpty.append(voxelParticleCount[voxelID])
            allMeanAnglesIncEmpty.append(localMeanAngle)
            voxelData[voxelID]= { "meanAngle": localMeanAngle,
                                 "totalPar": voxelParticleCount[voxelID],
                                 "foldedPar": voxelParticleCount[voxelID] - voxelUnfoldedCount[voxelID],
                                 "unfoldedPar": voxelUnfoldedCount[voxelID]}
        
        else:
            allCountsIncEmpty.append(0)
            allMeanAnglesIncEmpty.append(0)

    meanCount = np.mean(allCounts)
    meanAngle = np.mean(allMeanAngles)

    corrFunc = []
    for voxelID, data in voxelData.items():
        corr = ((data["totalPar"] - meanCount) / numMol) * (data["meanAngle"] - meanAngle)
        corrFunc.append(corr)

    #print(allCounts, allMeanAngles)
    plt.scatter(allMeanAngles, allCounts)
    plt.ylabel("count of particles in voxel")
    plt.xlabel("mean molecule angle")
    meanCorr = np.mean(corrFunc)
    print(meanCorr)
    plt.savefig(f"../runs/{conditions}/{filename}/analysis/figs/densityUnfoldingCorr.png")
    plt.clf()
    #plt.show()
    return allCounts, allMeanAngles, meanCorr, allCountsIncEmpty, allMeanAnglesIncEmpty


def calcAvDensityUnfoldedCorr(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, prob, suffixes):
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
        if bondsPerAtom == 2:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
        else:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
            if prob < 1:
                conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}_prob{prob}{suffix}"
        with open(f"../runs/{conditions}/timesteps.pkl", "rb") as f:
            timesteps = pickle.load(f)

        for runNum in range(numRuns):
            filename = f"Run{runNum}_{conditions}"
            with open(f"../runs/{conditions}/{filename}/analysis/densityUnfoldedCorr.pkl", "rb") as f:
                allCounts, allMeanAngles, meanCorr, allCountsIncEmpty, allMeanAnglesIncEmpty = pickle.load(f)
                print(len(allCountsIncEmpty))
            barrierAllCounts[barrier].append(allCountsIncEmpty)
            barrierAllMeanAngles[barrier].append(allMeanAnglesIncEmpty)
            allCorrFuncs[barrier].append(meanCorr)
        
        avCorrFuncsArray = np.array(allCorrFuncs[barrier])
        avCorrFuncs[barrier] = np.mean(avCorrFuncsArray, axis = 0)
        avCorrFuncsErr[barrier] = np.std(avCorrFuncsArray, axis = 0)
    
        avAllCountsArray = np.array(barrierAllCounts[barrier])
        avAllCounts[barrier] = np.mean(avAllCountsArray, axis = 0)
        avAllCountsErr[barrier] = np.std(avAllCountsArray, axis = 0)

        avAllMeanAnglesArray = np.array(barrierAllMeanAngles[barrier])
        avAllMeanAngles[barrier] = np.mean(barrierAllMeanAngles[barrier], axis = 0)
        avAllMeanAnglesErr[barrier] = np.std(barrierAllMeanAngles[barrier], axis = 0)

    return avAllCounts, avAllCountsErr, avAllMeanAngles, avAllMeanAnglesErr, avCorrFuncs, avCorrFuncsErr

def plotAvDensityUnfoldedCorr(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, prob, suffixes):
    with open(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/data/densityUnfoldedCorr.pkl", "rb") as f:
        avAllCounts, avAllCountsErr, avAllMeanAngles, avAllMeanAnglesErr, avCorrFuncs, avCorrFuncsErr = pickle.load(f)

    axes = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)]
    fig, ax = plt.subplots(2, 3, sharex = True, sharey = True)
    for barrier in unfoldBarriers:
        if bondsPerAtom == 2:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}"
        else:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}"
            if prob < 1:
                conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}_prob{prob}"
        key = axes[barrier - 1]
        with open(f"../runs/{conditions}/timesteps.pkl", "rb") as f:
            timesteps = pickle.load(f)
        for runNum in range(numRuns):
            filename = f"Run{runNum}_{conditions}"
            with open(f"../runs/{conditions}/{filename}/analysis/densityUnfoldedCorr.pkl", "rb") as f:
                allCounts, allMeanAngles, meanCorr, allCountsIncEmpty, allMeanAnglesIncEmpty = pickle.load(f)
            filtered = [(v, n) for v, n in zip(allCountsIncEmpty, allMeanAnglesIncEmpty)
                        if not (v == 0 and n == 0)]
            filteredCounts, filteredAngles = zip(*filtered)
        ax[key].hist2d(filteredAngles, filteredCounts, bins = 50, norm = mcolors.PowerNorm(0.5), cmap = plt.cm.summer)
        ax[key].set_title(f"E_U = {barrier}kT")
    fig.supxlabel(f"molecule angle (\u03B8)")
    fig.supylabel("number of particles in voxel")
    
    plt.tight_layout()
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/avDensityUnfoldedCorr.png")
    plt.show()
    plt.close()
    return
