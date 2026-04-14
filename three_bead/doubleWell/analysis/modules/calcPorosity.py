# maximal ball algorithm
# https://www.sciencedirect.com/science/article/pii/S0098300416305180?pes=vor&entityID=https%3A%2F%2Fpassport01.leeds.ac.uk%2Fidp%2Fshibboleth&utm_source=acs&getft_integrator=acs
# https://www.sciencedirect.com/science/article/pii/S037843710600464X?pes=vor&entityID=https%3A%2F%2Fpassport01.leeds.ac.uk%2Fidp%2Fshibboleth&utm_source=acs&getft_integrator=acs

import pickle
import numpy as np
from collections import defaultdict
from scipy.spatial import KDTree
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt

class Voxel():
    def __init__(self):
        self.pos = np.array([0, 0, 0])
        self.void = False
        self.filled = False
        self.poreSize = 0
        self.child = False
        self.parent = False

def calcPorosity(particles, barrier, refoldBarrier, runNum, vf, numMol, boxLength):
    parRadius = 2 ** (1/6)
    voxelSize = parRadius * 0.1

    # extract particle coordinates
    sphereCentres = np.zeros((len(particles), 3))
    for num, p in enumerate(particles):
        x = particles[num].properties[-1, 1]
        y = particles[num].properties[-1, 2]
        z = particles[num].properties[-1, 3]
        
        sphereCentres[num] = np.array([x, y, z])

    # find smallest and largest voxels a sphere can occupy
    boxMin = np.min(sphereCentres, axis = 0) - parRadius
    boxMax = np.max(sphereCentres, axis = 0) + parRadius
    
    # generate voxel grid
    gridX = np.arange(boxMin[0], boxMax[0], voxelSize)
    gridY = np.arange(boxMin[1], boxMax[1], voxelSize)
    gridZ = np.arange(boxMin[2], boxMax[2], voxelSize)
    
    gx, gy, gz = np.meshgrid(gridX, gridY, gridZ, indexing = 'ij')
    coords = np.column_stack([gx.ravel(), gy.ravel(), gz.ravel()])

    allVoxels = []

    for c in coords:
        v = Voxel()
        v.pos = c
        allVoxels.append(v)
    
    voxelCoords = np.array([v.pos for v in allVoxels])

    # find distance between pore centre and particles
    tree = KDTree(sphereCentres)
    dists, n = tree.query(voxelCoords, k = 1)

    voidVoxels = []

    for num, v in enumerate(allVoxels):
        if dists[num] < parRadius:
            v.filled = True
        else:
            v.void = True
            v.poreSize = dists[num] - parRadius
            voidVoxels.append(v)

    # maximum pore size at each voxel
    maxBallRadius = dists[dists >= parRadius] - parRadius
    print(maxBallRadius)

    # filter through parent and child voxels
    parents = []
    children = []
    for num1, i in enumerate(voidVoxels):
        for num2, j in enumerate(voidVoxels):
            if i.poreSize == j.poreSize:
                pass
            elif i.poreSize <= j.poreSize:
                dist = np.linalg.norm(i.pos - j.pos)

                if dist + j.poreSize <= i.poreSize:
                    i.parent = True
                    j.child = True

                    parents.append(i)
                    children.append(j)
    # decide which pores are main pores and which are throats

    mainPores = [v for v in voidVoxels if not v.child]

    return mainPores

def plotAvPorosity(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength):
    allPores = defaultdict(list)
    avPoreSize = {}
    avPoreSizeErr = {}

    parRadius = 2 ** (1/6)
    voxelSize = parRadius * 0.1


    for barrier in unfoldBarriers:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}"
        for runNum in range(numRuns):
            filename = f"Run{runNum}_{conditions}"
            fig, ax = plt.subplots()
            with open(f"../runs/{conditions}/{filename}/mainPores.pkl", "rb") as f:
                mainPores = pickle.load(f)
            allPores[barrier].append(v.poreSize for v in mainPores)
            allPoresArray = np.array(allPores[barrier])
            avPoreSize[barrier] = allPoresArray.mean(axis = 0)
            avPoreSizeErr[barrier] = allPoresArray.std(axis = 0)
            
            bins = np.arange(min(avPoreSize[barrier]), max(avPoreSize[barrier]), 0.1 * voxelSize)
            ax.hist(avPoreSize[barrier], bins, yerr = avPoreSizeErr[barrier],
                    label = f"unfolding barrier = {barrier}kT\nrun {runNum}")
            ax.set_xlabel("pore radius")
            ax.set_ylabel("number of pores")
            ax.set_title("pore size distribution")
            plt.savefig(f"../runs/{conditions}//averagedfigs/poreDistribution")
            plt.close()

    return
