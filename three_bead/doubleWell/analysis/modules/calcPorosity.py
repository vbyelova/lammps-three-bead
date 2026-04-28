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
    sphereCentres = np.zeros((len(particles), 3))
    for num, p in enumerate(particles):
        x = particles[num].properties[-1, 1]
        y = particles[num].properties[-1, 2]
        z = particles[num].properties[-1, 3]
        
        sphereCentres[num] = np.array([x, y, z])

    print("extracted particle coordinates for porosity..")

    # find smallest and largest voxels a sphere can occupy
    boxMin = np.min(sphereCentres, axis = 0) - parRadius
    boxMax = np.max(sphereCentres, axis = 0) + parRadius
    
    # generate voxel grid
    gridX = np.arange(boxMin[0], boxMax[0], voxelSize)
    gridY = np.arange(boxMin[1], boxMax[1], voxelSize)
    gridZ = np.arange(boxMin[2], boxMax[2], voxelSize)
    
    gx, gy, gz = np.meshgrid(gridX, gridY, gridZ, indexing = 'ij')
    coords = np.column_stack([gx.ravel(), gy.ravel(), gz.ravel()])


    allVoxels = [Voxel(pos = c) for c in coords]
    
    print("generated all voxels..")

    voxelCoords = np.array([v.pos for v in allVoxels])

    # find distance between pore centre and particles
    tree = KDTree(sphereCentres)
    filledVoxels = tree.query_ball_point(voxelCoords, r = parRadius, return_sorted = False)
    filledMask = np.array([len(id) > 0 for id in filledVoxels])

    print("generated mask..")

    for num, v in enumerate(allVoxels):
        if filledMask[num]:
            v.filled = True
        else:
            v.void = True
            dist, n = tree.query(v.pos, k = 1)
            v.poreSize = dist - parRadius

    voidVoxels = [v for v in allVoxels if v.void]

    print("found void voxels..")

    # filter through parent and child voxels
    for num1, i in enumerate(voidVoxels):
        for num2, j in enumerate(voidVoxels):
            if i == j:
                pass
            # is i a potential parent?
            elif i.poreSize > j.poreSize:
                dist = (i.pos[0] - j.pos[0])**2 + (i.pos[1] - j.pos[1])**2 +(i.pos[2] - j.pos[2])**2
                
                # does i fully encompass j?
                if dist + j.poreSize**2 <= i.poreSize**2:
                    i.parent = True
                    j.child = True

    # decide which pores are main pores and which are throats

    parents = [v for v in voidVoxels if not v.child]
    children = [v for v in voidVoxels if v.child]

    print("found parent and child voxels..")

    for num1, i in enumerate(parents):
        for num2, j in enumerate(parents):
            if i == j:
                pass
            dist = (i.pos[0] - j.pos[0])**2 + (i.pos[1] - j.pos[1])**2 +(i.pos[2] - j.pos[2])**2

            if dist < i.poreSize**2 + j.poreSize**2:
                i.throat = True
                j.throat = True

    mainPores = [v for v in parents if not v.throat]
    throatPores = [v for v in parents if v.throat]

    print("found main pores and throat pores..")

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
            plt.savefig(f"../runs/{conditions}/averagedfigs/poreDistribution")
            plt.close()

    return
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
    gridX = np.arange(boxMin[0], boxMax[0], voxelSize)
    gridY = np.arange(boxMin[1], boxMax[1], voxelSize)
    gridZ = np.arange(boxMin[2], boxMax[2], voxelSize)
    
    gx, gy, gz = np.meshgrid(gridX, gridY, gridZ, indexing = 'ij')
    coords = np.column_stack([gx.ravel(), gy.ravel(), gz.ravel()])


    allVoxels = [Voxel(pos = c) for c in coords]
    
    print("generated all voxels..")

    voxelCoords = np.array([v.pos for v in allVoxels])

    # find distance between pore centre and particles
    tree = KDTree(sphereCentres)
    filledVoxels = tree.query_ball_point(voxelCoords, r = parRadius, return_sorted = False)
    filledMask = np.array([len(id) > 0 for id in filledVoxels])

    print("generated mask..")

    for num, v in enumerate(allVoxels):
        if filledMask[num]:
            v.filled = True
        else:
            v.void = True
            dist, n = tree.query(v.pos, k = 1)
            v.poreSize = dist - parRadius

    voidVoxels = [v for v in allVoxels if v.void]

    print("found void voxels..")

    # filter through parent and child voxels
    for num1, i in enumerate(voidVoxels):
        for num2, j in enumerate(voidVoxels):
            if i == j:
                pass
            # is i a potential parent?
            elif i.poreSize > j.poreSize:
                dist = (i.pos[0] - j.pos[0])**2 + (i.pos[1] - j.pos[1])**2 +(i.pos[2] - j.pos[2])**2
                
                # does i fully encompass j?
                if dist + j.poreSize**2 <= i.poreSize**2:
                    i.parent = True
                    j.child = True

    # decide which pores are main pores and which are throats

    parents = [v for v in voidVoxels if not v.child]
    children = [v for v in voidVoxels if v.child]

    print("found parent and child voxels..")

    for num1, i in enumerate(parents):
        for num2, j in enumerate(parents):
            if i == j:
                pass
            dist = (i.pos[0] - j.pos[0])**2 + (i.pos[1] - j.pos[1])**2 +(i.pos[2] - j.pos[2])**2

            if dist < i.poreSize**2 + j.poreSize**2:
                i.throat = True
                j.throat = True

    mainPores = [v for v in parents if not v.throat]
    throatPores = [v for v in parents if v.throat]

    print("found main pores and throat pores..")

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
