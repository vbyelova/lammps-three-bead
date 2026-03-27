# maximal ball algorithm
# https://www.sciencedirect.com/science/article/pii/S0098300416305180?pes=vor&entityID=https%3A%2F%2Fpassport01.leeds.ac.uk%2Fidp%2Fshibboleth&utm_source=acs&getft_integrator=acs
# https://www.sciencedirect.com/science/article/pii/S037843710600464X?pes=vor&entityID=https%3A%2F%2Fpassport01.leeds.ac.uk%2Fidp%2Fshibboleth&utm_source=acs&getft_integrator=acs

import numpy as np
from scipy.spatial import KDTree
from scipy.spatial.distance import cdist

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
    allVoxels = np.column_stack([gx.ravel(), gy.ravel(), gz.ravel()])

    # find distance between pore centre and particles
    tree = KDTree(sphereCentres)
    dists, n = tree.query(allVoxels, k = 1)
    materialVoxels = allVoxels[dists < parRadius]
    voidVoxels = allVoxels[dists >= parRadius]

    # maximum pore size at each voxel
    maxBallRadius = dists[dists >= parRadius] - parRadius
    print(maxBallRadius)

    # filter through parent and child voxels
    parents = []
    children = []
    for i in maxBallRadius:
        for j in maxBallRadius:
            if i == j:
                pass
            elif i < j:
                parents.append(i)
                children.append(j)

    # decide which pores are main pores and which are throats

