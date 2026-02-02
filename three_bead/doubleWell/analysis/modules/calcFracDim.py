# a python script to calculate the fractal dimension of a network consisting
# of three-bead molecules.
# Victoria Byelova

import matplotlib.pyplot as plt
import numpy as np
import pickle
import pwlf

from collections import defaultdict
from .parseDump import *

def boxCounting(filename, boxLength):
    """A function to count the number of particles in given box intervals. First the number of voxels and
        their sizes is decided. Then bins are made to represent the voxels and digitize decides which voxel
        each particle belongs to. The unique voxels are counted and added to a list."""
    
    with open(filename, "rb") as f:
        particles = pickle.load(f)

    smallestBox = 0

    particleSize = 2**(1/6)
    halfLength = 0.5 * boxLength
    
    while 2**smallestBox <= boxLength:
        smallestBox += 1

    numVoxels = [2**_ for _ in range(0, smallestBox + 1)]
    voxelSizes = [boxLength / numDivs for numDivs in numVoxels] 
    totalUniqueVoxels = []
    print("initialised arrays")

    for vNum, v in enumerate(numVoxels):
        bins = np.linspace(- halfLength, halfLength, v + 1)
      
        uniqueVoxels = set()
        print("made bins for box size", vNum + 1, )
        print(bins)
        
        for pNum, p in enumerate(particles):
            x = p.properties[-1, 1]
            y = p.properties[-1, 2]
            z = p.properties[-1, 3]

            xVox = np.digitize(x, bins) - 1   #digitize returns 1-based, so -1
            yVox = np.digitize(y, bins) - 1
            zVox = np.digitize(z, bins) - 1

            uniqueVoxels.add((xVox, yVox, zVox))
        
        #print("unique box counted")
        #print(len(uniqueVoxels), v)
        totalUniqueVoxels.append([voxelSizes[vNum], len(uniqueVoxels)])

    return np.array(totalUniqueVoxels)

def calcFractalDimension(unfoldBarrier, refoldBarrier, runNum, Vf, numMol, totalUniqueVoxels):
    conditions = f"unfold{unfoldBarrier}_refold{refoldBarrier}_Vf{Vf}_mol{numMol}"
    filename = f"Run{runNum}_{conditions}"
    logNumUniqueVoxels = np.array([np.log(i) for i in totalUniqueVoxels[:,1]])
    logInverseVoxelSizes = np.array([np.log(1 / i) for i in totalUniqueVoxels[:,0]])
    print(logNumUniqueVoxels)
    print(logInverseVoxelSizes)
    plt.plot(logInverseVoxelSizes, logNumUniqueVoxels, '*', color = "hotpink")
    plt.legend([f"Vf = {Vf}\nunfolding barrier  = {unfoldBarrier}kT\n n. molecules = {numMol}"])
    plt.title("Finding fractal dimension by box counting method")
    plt.xlabel("log(1 / R)")
    plt.ylabel("log( N )")

    trend = np.polyfit(logInverseVoxelSizes, logNumUniqueVoxels, 1)
    trendpoly = np.poly1d(trend)
    #plt.plot(logInverseVoxelSizes, trendpoly(logInverseVoxelSizes), color = "purple")

    model = pwlf.PiecewiseLinFit(logInverseVoxelSizes, logNumUniqueVoxels)
    breaks = model.fitfast(2, pop = 3)
    xHat = np.linspace(min(logInverseVoxelSizes), max(logNumUniqueVoxels), 100)
    yHat = model.predict(xHat)
    plt.plot(xHat, yHat, "-", color = "purple")
    plt.title(f"Run {runNum}")
    plt.savefig(f"../runs/{conditions}/{filename}/analysis/figs/boxCounting")
    plt.close()
    #print(f"fractal dims: {trend[0]}")
    return trend[0]


def findAllFractalDims(UnfoldBarriers, refoldBarrier, numRuns, Vf, numMol, boxLength):
    fractalDim = defaultdict(list)
    fractalDimError = defaultdict(list)
    fractalDimMean = {}
    finalFractalDims = []
    finalFractalDimsError = []
    for barrier in UnfoldBarriers:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{Vf}_mol{numMol}"

        for runNum in range(0, numRuns):
            filename = f"Run{runNum}_{conditions}"

            totalUniqueVoxels = boxCounting(f"../runs/{conditions}/{filename}/analysis/dillParticles.pkl", boxLength)
            fractalDim[barrier].append(calcFractalDimension(barrier, refoldBarrier, runNum, Vf, numMol, totalUniqueVoxels))

        fractalDimMean[barrier] = sum(fractalDim[barrier]) / numRuns
        for val in fractalDim[barrier]:
            fractalDimError[barrier].append((val - fractalDimMean[barrier])**2)
    

        finalFractalDimsError.append(np.sqrt(sum(fractalDimError[barrier]) / numRuns))  
        finalFractalDims.append(sum(fractalDim[barrier]) / numRuns)

    for barrier in UnfoldBarriers:
        print(f"fractal dimensions for unfolding barrier = {barrier}: {fractalDim[barrier]}")
        
    return finalFractalDims, finalFractalDimsError

def plotAverageFractalDim(unfoldBarriers, refoldBarrier, numRuns, Vf, numMol, boxLength):
    """plots a scatter graph of what the average fractal dimension is with each unfolding barrier"""
    fractalDims = findAllFractalDims(unfoldBarriers, refoldBarrier, numRuns, Vf, numMol, boxLength)
    finalFractalDims, finalFractalDimsError = fractalDims[0], fractalDims[1]

    plt.bar(unfoldBarriers, finalFractalDims, yerr = finalFractalDimsError, color = "purple")
    plt.xlabel("Unfolding barrier (kT)")
    plt.ylabel("Fractal dimension")
    plt.title("Fractal dimension through box counting method")
    plt.show()

