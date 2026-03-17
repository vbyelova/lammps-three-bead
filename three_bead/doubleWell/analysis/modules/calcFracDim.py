# a python script to calculate the fractal dimension of a network consisting
# of three-bead molecules.
# Victoria Byelova

import matplotlib.pyplot as plt
import numpy as np
import pickle
import math
import pwlf


from collections import defaultdict
from .parseDump import *

def boxCounting(filename, boxLength, timesteps, percDims, atPerc):
    """A function to count the number of particles in given boxLength intervals. First the number of voxels and
        their sizes is decided. Then bins are made to represent the voxels and digitize decides which voxel
        each particle belongs to. The unique voxels are counted and added to a list."""
    
    with open(filename, "rb") as f:
        particles = pickle.load(f)

    particleSize = 2**(1/6)
    halfLength = 0.5 * boxLength
    
    totalUniqueVoxels = []

    voxelSizes = np.logspace(np.log10(particleSize), np.log10(halfLength), num = 10)
    numVoxels = [int(boxLength/vSize) for vSize in voxelSizes]
    voxelSizes = [boxLength/ numDivs for numDivs in numVoxels]
    print("initialised arrays for box counting..")

    if atPerc == True:
        for dim in percDims:
            if dim == 3:
                percAtFrame = percDims.index(dim)
                continue
    elif atPerc == False:
        percAtFrame = len(timesteps) - 1

    for vNum, v in enumerate(numVoxels):
        bins = np.linspace(- halfLength, halfLength, v + 1)
      
        uniqueVoxels = set()
        print("made bins for boxLength size", vNum + 1 )
        
        for pNum, p in enumerate(particles):
            x = p.properties[percAtFrame, 1]
            y = p.properties[percAtFrame, 2]
            z = p.properties[percAtFrame, 3]

            xVox = np.digitize(x, bins) - 1   #digitize returns 1-based, so -1
            yVox = np.digitize(y, bins) - 1
            zVox = np.digitize(z, bins) - 1

            uniqueVoxels.add((xVox, yVox, zVox))
        
        #print("unique boxLength counted")
        #print(len(uniqueVoxels), v)
        totalUniqueVoxels.append([voxelSizes[vNum], len(uniqueVoxels)])
    print(totalUniqueVoxels)
    return np.array(totalUniqueVoxels)
        

def calcFractalDimension(unfoldBarrier, refoldBarrier, runNum, vf, numMol, totalUniqueVoxels, boxLength, timesteps, percDims, atPerc):
    """calculates the fractal dimension for a single run"""
    conditions = f"unfold{unfoldBarrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}"
    filename = f"Run{runNum}_{conditions}"
    
    if atPerc == True:
        atPercString = "at percolation"
    elif atPerc == False:
        atPercString = "at end of sim"
    logNumUniqueVoxels = np.array([np.log(i) for i in totalUniqueVoxels[:,1]])
    logInverseVoxelSizes = np.array([np.log(1 / i) for i in totalUniqueVoxels[:,0]])

    plt.plot(logInverseVoxelSizes, logNumUniqueVoxels, '*', color = "hotpink")
    plt.legend([f"vf = {vf}\nunfolding barrier  = {unfoldBarrier}kT\n n. molecules = {numMol}"])
    plt.title("Finding fractal dimension by boxLength counting method")
    plt.xlabel("log(1 / R)")
    plt.ylabel("log( N )")
    plt.title(f"box counting {atPercString}")


    # if boxLength > 30:
    model = pwlf.PiecewiseLinFit(logInverseVoxelSizes, logNumUniqueVoxels)
    breaks = model.fit_guess([-1.5])

    xHat = np.linspace(min(logInverseVoxelSizes), max(logInverseVoxelSizes), 100)
    yHat = model.predict(xHat)
    plt.plot(xHat, yHat, "-", color = "purple")
    plt.xlim(math.floor(min(logInverseVoxelSizes)), int(max(logInverseVoxelSizes)))
    plt.tight_layout()
    plt.savefig(f"../runs/{conditions}/{filename}/analysis/figs/boxCounting")
    plt.close()
    return 2, np.abs(model.slopes[1]), breaks[1]

    # if boxLength < 30:
    #     trend = np.polyfit(logInverseVoxelSizes, logNumUniqueVoxels, 1)
    #     trendpoly = np.poly1d(trend)
    #     plt.plot(logInverseVoxelSizes, trendpoly(logInverseVoxelSizes), color = "purple")
    #     plt.tight_layout()
    #     plt.savefig(f"../runs/{conditions}/{filename}/analysis/figs/boxCounting")
    #     plt.close()
    #     return 1, trend[0]

def findAllFractalDims(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, timesteps, percDims, atPerc):
    """finds the fractal dimension for each run for each unfolding barrier and also calculates the average
        correlation length. Can find this for either at percolation point or in final frame of sim."""

    fractalDim = defaultdict(list)
    fractalDimError = defaultdict(list)

    fractalDimMean = {}
    finalFractalDims = []
    finalFractalDimsError = []

    corrLength = defaultdict(list)
    avCorrLength = {}
    avCorrLengthErr = {}

    for barrier in unfoldBarriers:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}"
        if atPerc == True:
            atPercString = "AtPerc"
        elif atPerc == False:
            atPercString = "SimEnd"
        with open(f"../runs/{conditions}/data/corrLengths{atPercString}.txt", "w") as f:
            f.write(f"# correlation length\t# standard deviation\n")

            for runNum in range(numRuns):
                filename = f"Run{runNum}_{conditions}"

                totalUniqueVoxels = boxCounting(f"../runs/{conditions}/{filename}/analysis/particleTraj.pkl",
                                                boxLength, timesteps, percDims, atPerc)
                dimsAndBreak = calcFractalDimension(barrier, refoldBarrier, runNum, vf, numMol, 
                                                    totalUniqueVoxels, boxLength, timesteps, percDims, atPerc)
                if dimsAndBreak[0] == 2:
                    df = dimsAndBreak[1]
                    breakpoint = dimsAndBreak[2]
                    fractalDim[barrier].append(df)
                    print("breakpoint = ", breakpoint)
                    
                    corrLength[barrier].append(10 ** -breakpoint)

                if dimsAndBreak[0] == 1:
                    fractalDim[barrier].append(dimsAndBreak[1])

            fractalDimMean[barrier] = sum(fractalDim[barrier]) / numRuns
            for val in fractalDim[barrier]:
                fractalDimError[barrier].append((val - fractalDimMean[barrier])**2)
        
            finalFractalDimsError.append(np.sqrt(sum(fractalDimError[barrier]) / numRuns))  
            finalFractalDims.append(sum(fractalDim[barrier]) / numRuns)

            print(f"fractal dimensions for unfolding barrier = {barrier}: {fractalDim[barrier]}")
            
            corrLengthArray = np.array(corrLength[barrier])
            avCorrLength[barrier] = corrLengthArray.mean(axis = 0)
            avCorrLengthErr[barrier] = corrLengthArray.std(axis = 0)
                
            f.write(f"{avCorrLength[barrier]}\t{avCorrLengthErr[barrier]}\n")

    return finalFractalDims, finalFractalDimsError, avCorrLength, avCorrLengthErr

def plotAverageFractalDim(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, timesteps, percDims):
    """plots a scatter graph of what the average fractal dimension is with each unfolding barrier
        at percolation and at the end of the simulation."""
    atPercFractalDims = findAllFractalDims(unfoldBarriers, refoldBarrier, numRuns, 
                                           vf, numMol, boxLength, timesteps, percDims, True)
    atPercFinalFractalDims, atPercFinalFractalDimsError = atPercFractalDims[0], atPercFractalDims[1]

    simEndFractalDims = findAllFractalDims(unfoldBarriers, refoldBarrier, numRuns, 
                                           vf, numMol, boxLength, timesteps, percDims, False)
    simEndFinalFractalDims, simEndFinalFractalDimsError = simEndFractalDims[0], simEndFractalDims[1]    

    barWidth = 0.35
    x = np.arange(len(unfoldBarriers))
    fig, ax = plt.subplots()
    ax.bar(x - barWidth / 2, atPercFinalFractalDims,
            yerr = atPercFinalFractalDimsError,
            width = barWidth, label = "at percolation", color = "hotpink")
    ax.bar(x + barWidth / 2, simEndFinalFractalDims,
            yerr = simEndFinalFractalDimsError,
            width = barWidth, label = "at end of sim", color = "royalblue")
    ax.set_xticks(x)
    ax.set_xticklabels(unfoldBarriers)
    ax.set_xlabel("Unfolding barrier (kT)")
    ax.set_ylabel("Fractal dimension")
    ax.set_title("Fractal dimension through box counting method")
    plt.savefig(f"../runs/boxLength{boxLength}/fractalDim_vf{vf}.png")
    plt.close()

    atPercAvCorrLength, atPercAvCorrLengthErr = atPercFractalDims[2], atPercFractalDims[3]
    simEndAvCorrLength, simEndAvCorrLengthErr = simEndFractalDims[2], simEndFractalDims[3]

    fig, ax = plt.subplots()
    ax.bar(x - barWidth / 2, [atPercAvCorrLength[b] for b in unfoldBarriers],
            yerr = [atPercAvCorrLengthErr[b] for b in unfoldBarriers],
            width = barWidth, color = "hotpink", label = "at percolation")
    ax.bar(x + barWidth / 2, [simEndAvCorrLength[b] for b in unfoldBarriers],
            yerr = [simEndAvCorrLengthErr[b] for b in unfoldBarriers],
            width = barWidth, color = "royalblue", label = "at end of sim")
    
    ax.set_xticks(x)
    ax.set_xticklabels(unfoldBarriers)
    ax.set_xlabel("Unfolding barrier (kT)")
    ax.set_ylabel(r"correlation length $\xi$")
    ax.set_title("Cluster sizes at different simulation points")
    plt.savefig(f"../runs/boxLength{boxLength}/corrLength_vf{vf}.png")
    plt.close()

    return

