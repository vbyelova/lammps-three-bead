import numpy as np
import matplotlib.pyplot as plt
import pickle
from collections import defaultdict

def densityUnfoldingDegreeCorrelation(barrier, refoldBarrier, vf, runNum, numMol, boxLength, voxelSize):

    voxels = defaultdict(list)
    localDensities = {}
    moleculeAngles = {}

    """correlation function between local density and degree of unfolding """
    conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}"
    filename = f"Run{runNum}_{conditions}"

    with open(f"../runs/{conditions}/{filename}/analysis/particleTraj.pkl", "rb") as f:
        particles = pickle.load(f)
    with open(f"../runs/{conditions}/{filename}/timesteps.pkl", "rb") as f:
        timesteps = pickle.load(f)
    # first find local density
    # split box into voxels
    gridX = np.arange(- 0.5 * boxLength, 0.5 * boxLength, voxelSize)
    gridY = np.arange(- 0.5 * boxLength, 0.5 * boxLength, voxelSize)
    gridZ = np.arange(- 0.5 * boxLength, 0.5 * boxLength, voxelSize)
    
    gx, gy, gz = np.meshgrid(gridX, gridY, gridZ, indexing = 'ij')
    voxelCoords = np.column_stack([gx.ravel(), gy.ravel(), gz.ravel()])

    # count how many particles per voxel (discrete) and calculate volume of these particles  
    for pNum, p in enumerate(particles):
        x = p.properties[-1, 1]
        y = p.properties[-1, 2]
        z = p.properties[-1, 3]

    # divide total particle volume by voxel volume

    # next let's find the degree of unfolding

    # for each molecule in the box, get the angle of unfolding
    # find the mean angle

    # apply these to a correlation function
    # loop over each voxel and r is how many voxels away
