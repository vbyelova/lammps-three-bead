# a python script to calculate the fractal dimension of a network consisting
# of three-bead molecules.
# Victoria Byelova

import matplotlib.pyplot as plt
import numpy as np
import pickle

from systemData import boxLength

def boxCounting(particles, boxLength):
    """A function to count the number of particles in given box intervals. First the number of voxels and
        their sizes is decided. Then bins are made to represent the voxels and digitize decides which voxel
        each particle belongs to. The unique voxels are counted and added to a list."""
    
    smallestBox = 0

    particleSize = 2**(1/6)
    halfLength = 0.5 * boxLength
    
    while 2**smallestBox < halfLength:
        smallestBox += 1

    numVoxels = [2**_ for _ in range(1, smallestBox)]
    voxelSizes = [boxLength / numDivs for numDivs in numVoxels] 
    totalUniqueVoxels = []
    print("initialised arrays")

    for vNum, v in enumerate(numVoxels):
        bins = np.linspace(- halfLength, halfLength, v + 1)
      
        uniqueVoxels = set()
        print("made bins for ", vNum)
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
        print(len(uniqueVoxels), v)
        totalUniqueVoxels.append([voxelSizes[vNum], len(uniqueVoxels)])

    return np.array(totalUniqueVoxels)

def calcFractalDimension(totalUniqueVoxels):
    logNumUniqueVoxels = np.array([np.log(i) for i in totalUniqueVoxels[:,1]])
    logInverseVoxelSizes = np.array([np.log(1 / i) for i in totalUniqueVoxels[:,0]])
    print(logNumUniqueVoxels)
    print(logInverseVoxelSizes)
    plt.plot(logInverseVoxelSizes, logNumUniqueVoxels, 'o')
    plt.xlabel("log inverse voxel sizes")
    plt.ylabel("log num unique voxels")

    trend = np.polyfit(logInverseVoxelSizes, logNumUniqueVoxels, 1)
    trendpoly = np.poly1d(trend)
    plt.plot(logInverseVoxelSizes, trendpoly(logInverseVoxelSizes))
    
    plt.show()
    print(trend[0])
    return

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



 #   for num, p in enumerate(particles):
        


    

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

with open("dillParticles.pkl", "rb") as f:
    particles = pickle.load(f)
totalUniqueVoxels = boxCounting(particles, boxLength)
calcFractalDimension(totalUniqueVoxels)

#sierpinskiTetrahedron()
