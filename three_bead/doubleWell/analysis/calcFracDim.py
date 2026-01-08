# a python script to calculate the fractal dimension of a network consisting
# of three-bead molecules.
# Victoria Byelova

import matplotlib.pyplot as plt
import numpy as np
import pickle

from collections import defaultdict
from parseDump import *

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

    numVoxels = [2**_ for _ in range(1, smallestBox)]
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

def calcFractalDimension(totalUniqueVoxels):
    logNumUniqueVoxels = np.array([np.log(i) for i in totalUniqueVoxels[:,1]])
    logInverseVoxelSizes = np.array([np.log(1 / i) for i in totalUniqueVoxels[:,0]])
    print(logNumUniqueVoxels)
    print(logInverseVoxelSizes)
    plt.plot(logInverseVoxelSizes, logNumUniqueVoxels, '*', color = "hotpink")
    plt.legend(["Vf = 0.07\nunfolding barrier  = 3kT\n n. molecules = 851"])
    plt.title("Finding fractal dimension by box counting method")
    plt.xlabel("log(1 / R)")
    plt.ylabel("log( N )")

    trend = np.polyfit(logInverseVoxelSizes, logNumUniqueVoxels, 1)
    trendpoly = np.poly1d(trend)
    plt.plot(logInverseVoxelSizes, trendpoly(logInverseVoxelSizes), color = "purple")
    plt.legend
    plt.show()
    #print(f"fractal dims: {trend[0]}")
    return trend[0]

def plotAllFractalDims(unfoldBarriers, refoldBarrier, numRuns, Vf, numMol):
    allUniqueVoxels = defaultdict(list)
    allVoxelSizes = defaultdict(list)

    for barrier in unfoldBarriers:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{Vf}_mol{numMol}"

        for n in range(0, numRuns):
            filename = f"threeBead_Run{n}_{conditions}"
            systemData = parseSystemData(f"../runs/{conditions}/{filename}/output/systemData.txt")
            boxLength, Npar, Nsteps, Nwrite, equilTime = systemData[0], systemData[1], systemData[2], systemData[3], systemData[4]        
            totalUniqueVoxels = boxCounting(f"../runs/{conditions}/{filename}/analysis/dillParticles.pkl", boxLength)
            
            logNumUniqueVoxels = np.array([np.log(i) for i in totalUniqueVoxels[:,1]])
            logInverseVoxelSizes = np.array([np.log(1 / i) for i in totalUniqueVoxels[:,0]])
            
            allUniqueVoxels[barrier].append(logNumUniqueVoxels)
            allVoxelSizes[barrier].append(logInverseVoxelSizes)

            trend = np.polyfit(logInverseVoxelSizes, logNumUniqueVoxels, 1)
            trendpoly = np.poly1d(trend)


def findAllFractalDims(UnfoldBarriers, refoldBarrier, numRuns, Vf, numMol):
    fractalDim = defaultdict(list)
    fractalDimError = defaultdict(list)
    fractalDimMean = {}
    finalFractalDims = []
    finalFractalDimsError = []
    for barrier in UnfoldBarriers:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{Vf}_mol{numMol}"

        for n in range(0, numRuns):
            filename = f"threeBead_Run{n}_{conditions}"
            systemData = parseSystemData(f"../runs/{conditions}/{filename}/output/systemData.txt")
            boxLength, Npar, Nsteps, Nwrite, equilTime = systemData[0], systemData[1], systemData[2], systemData[3], systemData[4]

            totalUniqueVoxels = boxCounting(f"../runs/{conditions}/{filename}/analysis/dillParticles.pkl", boxLength)
            fractalDim[barrier].append(calcFractalDimension(totalUniqueVoxels))

        fractalDimMean[barrier] = sum(fractalDim[barrier]) / numRuns
        for val in fractalDim[barrier]:
            fractalDimError[barrier].append((val - fractalDimMean[barrier])**2)
    

        finalFractalDimsError.append(np.sqrt(sum(fractalDimError[barrier]) / numRuns))  
        finalFractalDims.append(sum(fractalDim[barrier]) / numRuns)

    for barrier in UnfoldBarriers:
        print(f"fractal dimensions for unfolding barrier = {barrier}: {fractalDim[barrier]}")
        
    return finalFractalDims, finalFractalDimsError

    

def sierpinskiTetrahedron():
    """An example of how to calculate fractal dimension for a system, in this case our
        system is a sierpinski tetrahedron.
        Adapted from https://www.hellotriangle.io/post/fractal-geometry-sierpinski-tetrahedron
    """

    # Example of a 2-level sierpinski triangle      
#                           
#                           For a 3D tetrahedron:
#         /\                there are 4 faces (including the base)
#        /__\               that have this pattern.
#       /\  /\
#      /__\/__\
#     /\      /\   
#    /__\    /__\
#   /\  /\  /\  /\
#  /__\/__\/__\/__\

    levels = 4

    h = 1
    p1 = [- 0.5 * h, 0.0, - np.sqrt(3)/2 * h]
    p2 = [0.0, 0.0, 0.0]
    p3 = [0.5 * h, 0.0, - np.sqrt(3)/2 * h]
    p4 = [0.0, np.sqrt(3)/2, 0.0]
    particles = [p1, p2, p3, p4]

def sierpinskiCube():
    """An example of how to calculate fractal dimension for a system, in this case our
        system is a sierpinski cube, or a menger sponge.
    """

# Example of a 2D menger square
#
#  ___________________
# |  _      _      _  |
# | |_|    |_|    |_| |
# |       _____       |
# |  _   |     |   _  |
# | |_|  |     |  |_| |
# |      |_____|      |
# |  _      _      _  |
# | |_|    |_|    |_| |
# |___________________|




#sierpinskiTetrahedron()
