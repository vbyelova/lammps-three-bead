
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
    __slots__ = ['pos', 'void', 'filled', 'poreSize', 'child', 'parent', 'throat']
    def __init__(self, pos = None):
        self.pos = np.array([0, 0, 0])
        self.void = False
        self.filled = False
        self.poreSize = 0
        self.child = False
        self.parent = False
        self.throat = False

def calcPorosity(barrier, refoldBarrier, runNum, vf, numMol, boxLength):
    parRadius = 2 ** (1/6)
    voxelSize = parRadius * 0.2

    conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}"
    filename = f"Run{runNum}_{conditions}"
    with open(f"../runs/{conditions}/{filename}/analysis/particleTraj.pkl", "rb") as f:
        particles = pickle.load(f)

    # extract particle coordinates
    sphereCentres = np.array([p.properties[-1, 1:4] for p in particles])

    print("extracted particle coordinates for porosity..")

    # find smallest and largest voxels a sphere can occupy
    boxMin = np.min(sphereCentres, axis = 0) - parRadius
    boxMax = np.max(sphereCentres, axis = 0) + parRadius
    
    # generate voxel grid
    gridX = np.arange(boxMin[0], boxMax[0], voxelSize) + voxelSize / 2
    gridY = np.arange(boxMin[1], boxMax[1], voxelSize) + voxelSize / 2
    gridZ = np.arange(boxMin[2], boxMax[2], voxelSize) + voxelSize / 2
    
    gx, gy, gz = np.meshgrid(gridX, gridY, gridZ, indexing = 'ij')
    coords = np.column_stack([gx.ravel(), gy.ravel(), gz.ravel()])

    print("generated coords..")
    # find distance between pore centre and particles
    tree = KDTree(sphereCentres)
    filledMask = tree.query_ball_point(coords, r = parRadius, return_length = True) > 0

    print("generated mask..")
 
    voidCoords = coords[~filledMask]
    dists, _ = tree.query(voidCoords, k = 1)
    poreSizes = dists - parRadius

    voidTree = KDTree(voidCoords)
    maxPore = poreSizes.max()

    isChild = np.zeros(len(voidCoords), dtype = bool)


    # filter through parent and child voxels
    group = 2000
    for start in range(0, len(voidCoords), group):
        end = min(start + group, len(voidCoords))
        
        groupCoords = voidCoords[start:end]
        groupSizes = poreSizes[start:end]

        neighbours = voidTree.query_ball_point(groupCoords, r = maxPore)

        for i, (j, jSize, jNeigh) in enumerate(zip(range(start, end), groupSizes, neighbours)):
            if len(jNeigh) == 0:
                continue
            jNeighArray = np.array(jNeigh)
            iSize = poreSizes[jNeighArray]
            candidate = jNeighArray[iSize > jSize]
            diffs = voidCoords[candidate] - voidCoords[j]
            dist2 = np.einsum("ij, ij->i", diffs, diffs)
            contained = dist2 + jSize**2 <= poreSizes[candidate]**2
            if contained.any():
                isChild[j] = True
    
    print("filtered parents and children..")

    # decide which pores are main pores and which are throats

    parentMask = ~isChild
    parentCoords = voidCoords[parentMask]
    parentSize = poreSizes[parentMask]
    print("found parent and child voxels..")

    maxParent = parentSize.max()
    parentTree = KDTree(parentCoords)

    isThroat = np.zeros(len(parentCoords), dtype = bool)
    for i in range(len(parentCoords)):
        neighbours = parentTree.query_ball_point(parentCoords[i], r = (np.sqrt(parentSize[i]**2 + maxParent**2)))
        nb = np.array([n for n in neighbours if n != i])
        diffs = parentCoords[nb] - parentCoords[i]
        dist2 = np.einsum("ij, ij->i", diffs, diffs)
        overlaps = dist2 < parentSize[i]**2 + parentSize[nb]**2
        if overlaps.any():
            isThroat[i] = True
            isThroat[nb[overlaps]] = True
    
    mainMask = ~isThroat

    mainPoreCoords = parentCoords[mainMask]
    mainPoreSize = parentSize[mainMask]

    mainPores = [{"pos" : mainPoreCoords[k], "size": mainPoreSize[k]} for k in range(len(mainPoreCoords))]
    print("found pores and throats..")
    print(mainPores)
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
            with open(f"../runs/{conditions}/{filename}/analysis/mainPores.pkl", "rb") as f:
                mainPores = pickle.load(f)
            allPores[barrier].append([v.poreSize for v in mainPores])
            allPoresArray = np.array(allPores[barrier])
            avPoreSize[barrier] = allPoresArray.mean(axis = 0)
            avPoreSizeErr[barrier] = allPoresArray.std(axis = 0)
            
            bins = np.arange(min(avPoreSize[barrier]), max(avPoreSize[barrier]), 0.1 * voxelSize)
            ax.hist(avPoreSize[barrier], bins, yerr = avPoreSizeErr[barrier],
                    label = f"unfolding barrier = {barrier}kT\nrun {runNum}")
            ax.set_xlabel("pore radius")
            ax.set_ylabel("number of pores")
            ax.set_title("pore size distribution")
            plt.savefig(f"../runs/{conditions}/averagedfigs/poreDistribution")
            plt.close()

    return
