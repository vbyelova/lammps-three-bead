import numpy as np
import matplotlib.pyplot as plt
import pickle
from collections import defaultdict

from .calcAngles import *

def densityUnfoldingDegreeCorrelation(barrier, refoldBarrier, vf, runNum, numMol, boxLength, voxelSize):
    """correlation function between local density and degree of unfolding """

    voxelUnfoldedCount = defaultdict(int)
    voxelParticleCount = defaultdict(int)
    moleculeAngles = defaultdict(list)
    frame = -1
    conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}"
    filename = f"Run{runNum}_{conditions}"

    with open(f"../runs/{conditions}/{filename}/analysis/particleTraj.pkl", "rb") as f:
        particles = pickle.load(f)
    with open(f"../runs/{conditions}/{filename}/timesteps.pkl", "rb") as f:
        timesteps = pickle.load(f)
    # first find local density
    # split box into voxels
    bins = np.arange(-0.5 * boxLength, 0.5 * boxLength, voxelSize) 
    
    gx, gy, gz = np.meshgrid(bins, bins, bins, indexing = 'ij')
    coords = np.column_stack([gx.ravel(), gy.ravel(), gz.ravel()])

    
    angles = calcAngles(barrier, refoldBarrier, runNum, vf, numMol)
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
        if pNum // 3 == (pNum + 1) // 3 and pNum // 3 == (pNum - 1) // 3:
            molNum = pNum // 3
            # for each molecule in the box, get the angle of unfolding
            centralAngle = angles[frame, molNum]
            if centralAngle >= 120:
                voxelUnfoldedCount += 1
            else:
                continue
            moleculeAngles[voxelIndex].append(centralAngle)
        else:
            continue

    for voxelID in voxelParticleCount.keys():
        localMeanAngle = np.mean(moleculeAngles[voxelID])
        localNumberDensity = (voxelParticleCount[voxelID] - voxelUnfoldedCount[voxelID]) / voxelParticleCount[voxelID]

