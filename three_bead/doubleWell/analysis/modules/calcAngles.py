# all functions that require processing of the protein angles are here.

import numpy as np
import re
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as mcolors
import matplotlib.patches as patches
from itertools import cycle
from collections import defaultdict
from .parseDump import *
from scipy.optimize import curve_fit
import pickle
import pandas as pd

def calcAngles(barrier, refoldBarrier, runNum, vf, numMol, bondsPerAtom, prob, suffix, simType = "gelation"):
    if simType == "step strain":
        prefix = "stepStrain_"
    elif simType == "gelation":
        prefix = ""
    """saves angles of three bead molecules in each simulation frame."""
    if bondsPerAtom == 2:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
    else:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
        if prob < 1:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}_prob{prob}{suffix}"
    filename = f"Run{runNum}_{conditions}"
    systemData = parseSystemData(f"../runs/{conditions}/{filename}/{prefix}output/systemData.txt")
    boxLength, Npar, Nsteps, Nwrite, equilTime = systemData[0], systemData[1], systemData[2], systemData[3], systemData[4]
    
    angles = defaultdict(list)
    frame = 0
    counter = 0
    with open(f"../runs/{conditions}/{filename}/{prefix}output/moleculeangles.dat", "r") as f:
        print("opened file")
        for line in f.readlines():
            if re.search(r"\d+\s[\d.]+\s-?[\d.]+", line):
                angles[frame].append(float(line.rsplit()[1]))
                counter += 1
                if counter == numMol:
                    frame += 1
                    counter = 0
    
    print("parsed angles..")
    return angles

def unfoldedMolecules(barrier, refoldBarrier, runNum, vf, numMol, angles, bondsPerAtom, prob, suffix, simType = "gelation"):
    """makes a list of molecules that are unfolded (have an angle of 120-180). Calculates the unfolded angle
        based on the bimodal distribution"""
    if bondsPerAtom == 2:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
    else:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
        if prob < 1:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}_prob{prob}{suffix}"

    with open(f"../runs/{conditions}/angleInfo.in", "r") as f:
        lines = f.readlines()
        info = lines[2].rsplit()
        maximum = float(info[2]) * 180 / np.pi
    unfoldedMols = defaultdict(list)
    for frame in range(len(angles)):
        for molID, angle in enumerate(angles[frame]):
             if angle > maximum:
                 unfoldedMols[frame].append(molID)

    print("generated list of unfolded molecules..")
    return unfoldedMols

def unfoldedOverTime(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, prob, suffixes, simType = "gelation"):
    """calculates number of molecules unfolded over course of simulation."""
    allUnfold = defaultdict(list)
    avUnfold = {}
    avUnfoldErr = {}

    with open(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/{prefix}data/avPercolation.pkl", "rb") as f:
        avPercolation, avPercolationErr = pickle.load(f)
    if simType == "step strain":
        prefix = "stepStrain_"
    elif simType == "gelation":
        prefix = ""
    for barrier, suffix in zip(unfoldBarriers, suffixes):
        key = (barrier, suffix)
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
            with open(f"../runs/{conditions}/{filename}/{prefix}analysis/nBonds.pkl", "rb") as f:
                nBonds = pickle.load(f)
            angles = calcAngles(barrier, refoldBarrier, runNum, vf, numMol, bondsPerAtom, prob, suffix, simType = "gelation")
            unfoldedMols = unfoldedMolecules(barrier, refoldBarrier, runNum, vf, numMol, angles, bondsPerAtom, prob, suffix, simType = "gelation")
            totalMols = [len(unfoldedMols[frame]) for frame in range(len(timesteps))]
            allUnfold[key].append(totalMols)

        avUnfoldArray = np.array(allUnfold[key])
        avUnfold[key] = avUnfoldArray.mean(axis = 0)
        avUnfoldErr[key] = avUnfoldArray.std(axis = 0)
    
    return avUnfold, avUnfoldErr

def plotUnfoldedOverTime(unfoldBarrier, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, prob, suffixes, simType = "gelation"):
    """produces a plot of the percentage of molecules in the systems that have unfolded with time. """
    percAt = []
    perc = []
    noperc = []
    if simType == "step strain":
        prefix = "stepStrain_"
    elif simType == "gelation":
        prefix = ""
    #colours = ["maroon", "firebrick", "indianred", "lightcoral", "lightpink"]
    colours = (["black", "dimgrey", "grey", "darkgrey", "silver"])

    fig, ax = plt.subplots()
    with open(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/{prefix}data/avPercolation.pkl", "rb") as f:
        avPercolation, avPercolationErr = pickle.load(f)
    with open(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/{prefix}data/unfoldedOverTime.pkl", "rb") as f:
        avUnfold, avUnfoldErr = pickle.load(f)

        
    for i, (barrier, suffix) in enumerate(zip(unfoldBarrier, suffixes)):
        if bondsPerAtom == 2:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
        else:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
            if prob < 1:
                conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}_prob{prob}{suffix}"
        with open(f"../runs/{conditions}/timesteps.pkl", "rb") as f:
            timesteps = pickle.load(f) 
        key = (barrier, suffix)
        if suffix != "":
            percAt.append(None)
            noperc.append((key))
        if suffix == "":
            perc.append(key)
            percAt.append(np.where(avPercolation[barrier] == 3)[0][0])

            ax.scatter(timesteps[percAt[i]], avUnfold[key][percAt[i]] / numMol * 100, color = "deeppink", marker = "x", zorder = 10)
            print(f"barrier = {barrier}kT, percolates at {timesteps[percAt[i]]}")
    for i, barrier in enumerate(sorted(set(unfoldBarrier))):
       
        ax.errorbar(timesteps, avUnfold[perc[i]] / numMol * 100, yerr = avUnfoldErr[perc[i]] / numMol * 100,
                    label = f"$E_U$ = {barrier}kT", color = colours[i])
        if len(noperc) > 0:
            ax.errorbar(timesteps, avUnfold[noperc[i]] / numMol * 100, yerr = avUnfoldErr[noperc[i]] / numMol * 100,
                        color = colours[i], linestyle = "dashdot")
        else:
            pass

    ax.set_xlabel("simulation frame")
    ax.set_ylabel(f"number of unfolded molecules (%)")
    ax.legend()
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/unfoldedMol_vf{vf}.jpg")
    ax.semilogx()
    ax.set_xlabel("timesteps", fontsize = 12)
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/unfoldedMol_vf{vf}_semilog.jpg")
    ax.semilogy()
    ax.set_ylabel("number of unfolded molecules", fontsize = 12)
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/unfoldedMol_vf{vf}_log.jpg")

    plt.close()

    print("plotted graphs of total unfolded molecules over time..")
    return 


def bimodalAngle(barrier, refoldBarrier, runNum, vf, numMol, bondsPerAtom, prob, suffix, simType = "gelation"):
    """plots a histogram of the molecule angles in the system."""
    if simType == "step strain":
        prefix = "stepStrain_"
    elif simType == "gelation":
        prefix = ""
    if bondsPerAtom == 2:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
    else:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
        if prob < 1:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}_prob{prob}{suffix}"
    filename = f"Run{runNum}_{conditions}"
    angles = calcAngles(barrier, refoldBarrier, runNum, vf, numMol, bondsPerAtom, prob, suffix, simType = "gelation")
    finalFrameAngles = angles[int(len(angles)-1)]
    hist, bins = np.histogram(finalFrameAngles, bins = np.linspace(0,200, 31))
    plt.hist(finalFrameAngles, bins = bins, color = "pink")
    plt.legend([f"vf = {vf}\nunfolding barrier  = {barrier} kT\n num. mol. = {numMol}"])
    plt.xlabel("molecule angle \u03B8")
    plt.ylabel("number of molecules")
    plt.title("Bimodal distribution of folded and unfolded molecules")
    plt.savefig(f"../runs/{conditions}/{filename}/{prefix}analysis/figs/bimodalAngleDist.png")
    plt.close()
    print("plotted histogram of angle distribution in system..")
    return hist, bins, finalFrameAngles


def densityUnfoldedCorrelation(barrier, refoldBarrier, runNum , vf, numMol, boxLength, bondsPerAtom, prob, suffix, simType = "gelation"):
    """correlation function between local density and degree of unfolding """
    if simType == "step strain":
        prefix = "stepStrain_"
    elif simType == "gelation":
        prefix = ""
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
    with open(f"../runs/{conditions}/{filename}/{prefix}analysis/particleTraj.pkl", "rb") as f:
        particles = pickle.load(f)
    with open(f"../runs/{conditions}/timesteps.pkl", "rb") as f:
        timesteps = pickle.load(f)
    # first find local density
    # split box into voxels
    bins = np.arange(-0.5 * boxLength, 0.5 * boxLength, voxelSize) 
    
    gx, gy, gz = np.meshgrid(bins, bins, bins, indexing = 'ij')
    coords = np.column_stack([gx.ravel(), gy.ravel(), gz.ravel()])

    
    angles = calcAngles(barrier, refoldBarrier, runNum, vf, numMol, bondsPerAtom, prob, suffix, simType = "gelation")
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
    plt.savefig(f"../runs/{conditions}/{filename}/{prefix}analysis/figs/densityUnfoldingCorr.jpg")
    plt.clf()
    #plt.show()
    return allCounts, allMeanAngles, meanCorr, allCountsIncEmpty, allMeanAnglesIncEmpty


def calcAvDensityUnfoldedCorr(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, prob, suffixes, simType = "gelation"):
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
    if simType == "step strain":
        prefix = "stepStrain_"
    elif simType == "gelation":
        prefix = ""
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
            with open(f"../runs/{conditions}/{filename}/{prefix}analysis/densityUnfoldedCorr.pkl", "rb") as f:
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

def plotAvDensityUnfoldedCorr(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, prob, suffixes, simType = "gelation"):
    """plots heat maps of local densities and local average degree of unfolding."""
    with open(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/{prefix}data/densityUnfoldedCorr.pkl", "rb") as f:
        avAllCounts, avAllCountsErr, avAllMeanAngles, avAllMeanAnglesErr, avCorrFuncs, avCorrFuncsErr = pickle.load(f)
    if simType == "step strain":
        prefix = "stepStrain_"
    elif simType == "gelation":
        prefix = ""
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
            with open(f"../runs/{conditions}/{filename}/{prefix}analysis/densityUnfoldedCorr.pkl", "rb") as f:
                allCounts, allMeanAngles, meanCorr, allCountsIncEmpty, allMeanAnglesIncEmpty = pickle.load(f)
            filtered = [(v, n) for v, n in zip(allCountsIncEmpty, allMeanAnglesIncEmpty)
                        if not (v == 0 and n == 0)]
            filteredCounts, filteredAngles = zip(*filtered)
        ax[key].hist2d(filteredAngles, filteredCounts, bins = 50, norm = mcolors.PowerNorm(0.5), cmap = plt.cm.summer)
        ax[key].set_title(f"E_U = {barrier}kT")
    fig.supxlabel(f"molecule angle (\u03B8)")
    fig.supylabel("number of particles in voxel")
    
    plt.tight_layout()
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/avDensityUnfoldedCorr.jpg")
    plt.show()
    plt.close()
    return


def plotBimodalAndDensityUnfoldedCorr(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, prob, suffixes, simType = "gelation"):
    """plots the average distribution of molecule angles in the system on the top and the heat map distribution on the bottom.
        The heat maps depict average degree of unfolding in a voxel against density of a voxel, essentially describing how
        low density and high density regions look."""
    if simType == "step strain":
        prefix = "stepStrain_"
    elif simType == "gelation":
        prefix = ""    
    allHistograms = defaultdict(list)
    histogramsMean = {}
    histogramsErr = {}
    fig, axs = plt.subplots(2, 5, sharex = "row", sharey = "row", figsize = (20, 8))
    labels1 = [(0, 0), (0, 1), (0, 2), (0, 3), (0, 4)]
    labels2 = [(1, 0), (1, 1), (1, 2), (1, 3), (1, 4)]


    allFilteredCounts = []
    allFilteredAngles = []

    colours = (["black", "dimgrey", "grey", "darkgrey", "silver"])
    for i, (barrier, suffix) in enumerate(zip(unfoldBarriers, suffixes)):
        ax = axs[labels1[barrier - 1]]
        if bondsPerAtom == 2:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
        else:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
            if prob < 1:
                conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}_prob{prob}{suffix}"
        for runNum in range(numRuns):
            filename = f"Run{runNum}_{conditions}"
            with open(f"../runs/{conditions}/{filename}/{prefix}analysis/angleDist.pkl", "rb") as f:
                hist, bins, finalFrameAngles = pickle.load(f)
            allHistograms[barrier].append(hist)
        key = (barrier, suffix)
        allHistogramsArray = np.array(allHistograms[barrier])
        histogramsMean[key] = np.mean(allHistogramsArray, axis = 0)
        histogramsErr[key] = np.std(allHistogramsArray, axis = 0)
        ax.set_xticks([60, 120, 180])
        if suffix != "":
            simtype = f"\nno intermol."
        else:
            simtype = ""

        bin_centers = (bins[:-1] + bins[1:]) / 2
        ax.bar(bin_centers, histogramsMean[key], width=np.diff(bins),
            yerr=histogramsErr[key], color = colours[barrier - 1],
            label=f"$E_U$ = {barrier}kT")
    
        
        ax.set_xlim([20, 200])
        ax.set_ylim([0, np.max(histogramsMean[key] + 500)])
        ax.axvline(60, linestyle = "--", color = "mediumblue")
        ax.axvline(180, linestyle = "--", color = "red")
        ax.set_title(f"$E_U$ = {barrier}kT", fontsize = 15)

        # check the difference between mode and mean of the raw distribution
        mode_angle = bin_centers[np.argmax(histogramsMean[key])]
        mean_angle = np.average(bin_centers, weights = histogramsMean[key])
        print(mode_angle, mean_angle)

    for barrier in unfoldBarriers:
        if bondsPerAtom == 2:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
        else:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
            if prob < 1:
                conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}_prob{prob}{suffix}"
        ax = axs[labels2[barrier - 1]]
        with open(f"../runs/{conditions}/timesteps.pkl", "rb") as f:
            timesteps = pickle.load(f)
        for runNum in range(numRuns):
            filename = f"Run{runNum}_{conditions}"
            with open(f"../runs/{conditions}/{filename}/{prefix}analysis/densityUnfoldedCorr.pkl", "rb") as f:
                allCounts, allMeanAngles, meanCorr, allCountsIncEmpty, allMeanAnglesIncEmpty = pickle.load(f)
            filtered = [(v, n) for v, n in zip(allCountsIncEmpty, allMeanAnglesIncEmpty)
                        if not (v == 0 and n == 0)]
            filteredCounts, filteredAngles = zip(*filtered)
            allFilteredCounts.extend(filteredCounts)
            allFilteredAngles.extend(filteredAngles)

        h = ax.hist2d(allFilteredAngles, allFilteredCounts, bins = (40, 40), range = [[20, 200], [0, 200]], norm = mcolors.PowerNorm(0.4), cmap = plt.cm.plasma)
        if barrier == 5: 
            ax.add_patch(patches.Rectangle(xy = (120, 0), width = 65, height = 25, linewidth = 1.5, color = "white", linestyle = "--", fill = False))
    ax.set_xticks([60, 120, 180])
    fig.text(0.08, 0.65, "count of molecules", va='center', rotation='vertical', fontsize=15)
    fig.text(0.08, 0.28, "count of molecules in voxel", va='center', rotation='vertical', fontsize=15)

    #fig.text(0.5, 0.52, "molecule angle", ha='center', fontsize=12)
    fig.text(0.5, 0.04, f"molecule angle (\u03B8)", ha='center', fontsize=15)

    #plt.tight_layout()
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/bimodalAndDensityUnfoldedCorr.jpg")
    return



def anglePopulation(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, prob, suffixes, simType = "gelation"):
    """plots a histogram of average distribution of folded vs unfolded molecules with
        energy barrier."""
    if simType == "step strain":
        prefix = "stepStrain_"
    elif simType == "gelation":
        prefix = ""
    totalFolded = defaultdict(list)
    totalUnfolded = defaultdict(list)
    foldedMean = {}
    unfoldedMean = {}
    foldedError = {}
    unfoldedError = {}
    plotLabels = []
    fig, ax = plt.subplots()
    barWidth = 0.4 
    perc = []
    noperc = []
    x = np.arange(len(unfoldBarriers))
    for barrier, suffix in zip(unfoldBarriers, suffixes):
        if suffix != "":
            plotLabels.append(f"{barrier} no intermol.")
            groupID = 0
            noperc.append((barrier, suffix))
        if suffix == "":
            plotLabels.append(f"{barrier}")
            groupID = 1
            perc.append((barrier, suffix))
        if bondsPerAtom == 2:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
        else:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
            if prob < 1:
                conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}_prob{prob}{suffix}"

        with open(f"../runs/{conditions}/timesteps.pkl", "rb") as f:
            timesteps = pickle.load(f)
        finalframe = len(timesteps) - 1

        for runNum in range(numRuns):
            filename = f"Run{runNum}_{conditions}"
            with open(f"../runs/{conditions}/{filename}/{prefix}analysis/nBonds.pkl", "rb") as f:
                nBonds = pickle.load(f)
            angles = calcAngles(barrier, refoldBarrier, runNum, vf, numMol, bondsPerAtom, prob, suffix, simType = "gelation")
            unfoldedMol = unfoldedMolecules(barrier, refoldBarrier, runNum, vf, numMol, angles, bondsPerAtom, prob, suffix, simType = "gelation")
            print(f"barrier {barrier}, run {runNum}")
        
            numUnfolded = len(unfoldedMol[finalframe])
            numFolded = numMol - numUnfolded
            
            totalFolded[(barrier, suffix)].append(numFolded)
            totalUnfolded[(barrier, suffix)].append(numUnfolded)
            
        foldedArray = np.array(totalFolded[(barrier, suffix)])
        unfoldedArray = np.array(totalUnfolded[(barrier, suffix)])
            
        foldedMean[(barrier, suffix)] = foldedArray.mean()
        unfoldedMean[(barrier, suffix)] = unfoldedArray.mean()
        foldedError[(barrier, suffix)] = (foldedArray.std())
        unfoldedError[(barrier, suffix)] =(unfoldedArray.std())

        print(f"unfolded percent is {unfoldedMean[barrier, suffix] / numMol * 100}, error is {unfoldedError[(barrier, suffix)] / numMol * 100}")

    x = np.arange(len(sorted(set(unfoldBarriers))))

    offset = barWidth * 0.5

    ax.bar(x - offset, [foldedMean[key] / numMol * 100 for key in perc],
          width = barWidth, yerr =  [foldedError[key] / numMol * 100 for key in perc], color = "darkblue", label = "folded")
    ax.bar(x - offset, [unfoldedMean[key] / numMol * 100 for key in perc], width = barWidth, yerr = [unfoldedError[key] / numMol * 100 for key in perc], color = "red", bottom = [foldedMean[key] / numMol * 100 for key in perc], label = "unfolded")
    if len(noperc) > 0:
        ax.bar(x + offset, [foldedMean[key] / numMol * 100 for key in noperc],
            width = barWidth, yerr =  [foldedError[key] / numMol * 100 for key in noperc], color = "dimgrey", label = "folded (no intermol.)")
        ax.bar(x + offset, [unfoldedMean[key] / numMol * 100 for key in noperc], width = barWidth, yerr = [unfoldedError[key] / numMol * 100 for key in noperc], color = "lightgrey", bottom = [foldedMean[key] / numMol * 100 for key in noperc], label = "unfolded (no intermol.)")    
    else:
        pass
    ax.set_xticks(x)
    ax.set_xticklabels(sorted(set(unfoldBarriers)))
    ax.set_xlabel(f"$E_U$ (kT)", fontsize = 12)
    ax.set_ylabel("molecule count (%)", fontsize = 12)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
        fancybox=True, shadow=True, ncol=2)
    plt.tight_layout()
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/barrierAnglePop_vf{vf}.jpg")
    print("plotted folded vs unfolded population by simulation end..")
    return

def avUnfoldPerFrame(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, prob, suffixes, simType = "gelation"):
    """plots a scatter graph of average unfolding events per timestep,
        compares against multiple unfolding barriers."""
    if simType == "step strain":
        prefix = "stepStrain_"
    elif simType == "gelation":
        prefix = ""
    totalNewUnfoldedMol = defaultdict(list)
    avNewUnfoldedMol = {}
    avNewUnfoldedMolErr = {}
    fig, ax = plt.subplots()

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
            nBonds = totalBonds(f"../runs/{conditions}/{filename}/{prefix}output/nbonds.dat")
            #print(f"len nBonds {len(nBonds)}, run{runNum} barrier{barrier}")
            if len(nBonds) < 100: 
                print(f"skipping run{runNum} barrier{barrier}")
                continue
            angles = calcAngles(barrier, refoldBarrier, runNum, vf, numMol, bondsPerAtom, prob, suffix, simType = "gelation")
            unfoldedMols = unfoldedMolecules(barrier, refoldBarrier, runNum, vf, numMol, angles, bondsPerAtom, prob, suffix, simType = "gelation")
            
            newUnfoldedCount = [0]
            #print(nBonds)
            for frame in range(1, len(nBonds)):
                newUnfoldedCount.append(len(unfoldedMols[frame]) - len(unfoldedMols[frame - 1]))
            
            totalNewUnfoldedMol[barrier].append(newUnfoldedCount)
            #print(totalNewUnfoldedMol[barrier])
        
        newUnfoldedArray = np.array(totalNewUnfoldedMol[barrier])
        avNewUnfoldedMol[barrier] = newUnfoldedArray.mean(axis = 0)
        avNewUnfoldedMolErr[barrier] = newUnfoldedArray.std(axis = 0)
        ax.errorbar(timesteps, avNewUnfoldedMol[barrier], yerr = avNewUnfoldedMolErr[barrier],
                label = f"unfolding barrier = {barrier}kT {suffix}")

    ax.set_xlabel("simulation frame")
    ax.set_ylabel("number of new unfolding events")
    ax.legend()
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/unfoldingEvents_vf{vf}.jpg")

    print("plotted new unfolding events..")

    ax.semilogx()
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/unfoldingEvents_vf{vf}_semilog.jpg")
    plt.close()

    print("plotted unfolding events on semilog axis..")
    return avNewUnfoldedMol, avNewUnfoldedMolErr


def plotCrosslinkPercentAtPercAndSimEnd(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, prob, suffixes, simType = "gelation"):
    """plots the percentage of possible crosslinks made at percolation and at the end of the simulation. Also plots
        percentage of molecules in the system that have unfolded at percolation and at the end of the simulation."""
    percAt = []
    colours = (["black", "dimgrey", "grey", "darkgrey", "silver"])
    allnBonds = defaultdict(list)
    avnBonds = {}
    avnBondsErr = {}
    if simType == "step strain":
        prefix = "stepStrain_"
    elif simType == "gelation":
        prefix = ""
    fig, ax1 = plt.subplots()
    for i, (barrier, suffix) in enumerate(zip(unfoldBarriers, suffixes)):
        if bondsPerAtom == 2:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
        else:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
            if prob < 1:
                conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}_prob{prob}{suffix}"
        with open(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/{prefix}data/avPercolation.pkl", "rb") as f:
            avPercolation, avPercolationErr = pickle.load(f)
        with open(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/{prefix}data/unfoldedOverTime.pkl", "rb") as f:
            avUnfold, avUnfoldErr = pickle.load(f)
        for runNum in range(numRuns):
            filename = f"Run{runNum}_{conditions}"
            with open(f"../runs/{conditions}/{filename}/{prefix}analysis/nBonds.pkl", "rb") as f:
                nBonds = pickle.load(f)
            allnBonds[barrier].append(nBonds)
        avnBondsArray = np.array(allnBonds[barrier])
        percentage = ((avnBondsArray - numMol * 2) / (numMol * 3) )* 100
        avnBonds[barrier] = percentage.mean(axis = 0)
        avnBondsErr[barrier] = percentage.std(axis = 0)
        percAt.append(np.where(avPercolation[barrier] == 3)[0][0])
    atPercBonds = ax1.errorbar(unfoldBarriers, [avnBonds[b][percAt[i]]for b in unfoldBarriers],
                                yerr = [avnBondsErr[b][percAt[i]] for b in unfoldBarriers], color = "cornflowerblue", fmt = "D--", mfc = "none", label = "at gelation")
    simEndBonds = ax1.errorbar(unfoldBarriers, [avnBonds[b][-1] for b in unfoldBarriers],
                                yerr = [avnBondsErr[b][-1] for b in unfoldBarriers], color = "cornflowerblue", fmt = "D--", label = "at simulation end")
    ax1.set_xlabel("$E_U$ (kT)", fontsize = 12)
    ax1.set_xticks([1, 2, 3, 4, 5])
    ax1.set_xticklabels([1, 2, 3, 4, 5])
    ax1.set_ylabel("possible crosslinks made (%)", fontsize = 12, color = "cornflowerblue")
    ax2 = ax1.twinx()

    atPercUnfold = ax2.errorbar(unfoldBarriers, [(avUnfold[b, s][percAt[i]] / numMol) * 100 for b, s in zip(unfoldBarriers, suffixes) ],
                  yerr = [(avUnfoldErr[b, s][percAt[i]] / numMol) * 100 for b, s in zip(unfoldBarriers, suffixes)],
                color = "red", fmt = "D--", mfc = "none", label = "at gelation")
    simEndUnfold = ax2.errorbar(unfoldBarriers, [(avUnfold[b, s][-1] / numMol) * 100 for b, s in zip(unfoldBarriers, suffixes)],
                  yerr = [(avUnfoldErr[b, s][-1] / numMol) * 100 for b, s in zip(unfoldBarriers, suffixes)], color = "red", fmt = "D--",
                 label = "at simulation end")
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc = "right")
    ax2.set_ylabel("unfolded molecules (%)", color = "red")
    plt.tight_layout()
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/crosslinkPercent.jpg")

    return
