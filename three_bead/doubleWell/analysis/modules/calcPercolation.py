# python code to calculate percolation of three-bead network.

import pickle
import matplotlib.pyplot as plt
# import code by Dr. David Head
from .percPerBC import *
from .parseDump import *

def finalFrame(percFile, totalBonds, Npar):
    with open(percFile, "rb") as f:
        bonds = pickle.load(f)

    V = [n for n in range(0, Npar + 2)]


    totalFrames = len(totalBonds)
    E = bonds.get(totalFrames - 1)

    system = percPerBC(3)
    return print("Maximum spanning dimension of a component: {}".format(system.percolationDimension(V, E)))

def allFrames(percFile, totalBonds, Npar):
    dimsPercolated = []
    with open(percFile, "rb") as f:
        bonds = pickle.load(f)
    
    V = [n for n in range(0, Npar + 2)]

    system = percPerBC(3)

    totalFrames = len(totalBonds)
    for val in range(0, totalFrames + 1):
        E = bonds.get(val)
        dimsPercolated.append(system.percolationDimension(V, E))

    plt.plot(dimsPercolated, (range(0, totalFrames + 1)))
    plt.show()
    return

def testPercPerBC():
    V = [n for n in range(0, 13)]
    E = [
         [0, 1, [0, 0, 0]],
         [1, 2, [0, 0, 0]],
         [2, 3, [0, 0, 0]],
         [3, 4, [0, 0, 0]],
         [4, 0, [0, 1, 0]],
         [2, 5, [0, 0, 0]],
         [5, 6, [0, 0, 0]],
         [6, 7, [1, 0, 0]],
         [7, 8, [0, 0, 0]],
         [8, 2, [0, 0, 0]],
         [9, 2, [0, 0, 0]],
         [9, 10, [0, 0, 0]],
         [10, 11, [0, 0, 1]],
         [11, 12, [0, 0, 0]],
         [12, 2, [0, 0, 0]]
        ]

    system = percPerBC(3)
    print( "Maximum spanning dimension of a component: {}".format(system.percolationDimension(V,E) ) )
    print( system )

def testFrameByFrame():
    """tests multiple frames of percolation to see if dimensions ranging 1-3 can be achieved."""
    V = [n for n in range(0, 13)]
    nBonds = [12, 13, 14, 15]
    framesE = {0:[[0, 1, [0, 0, 0]],
                    [1, 2, [0, 0, 0]],
                    [2, 3, [0, 0, 0]],
                    [3, 4, [0, 0, 0]],
                    [2, 5, [0, 0, 0]],
                    [5, 6, [0, 0, 0]],
                    [7, 8, [0, 0, 0]],
                    [8, 2, [0, 0, 0]],
                    [9, 2, [0, 0, 0]],
                    [9, 10, [0, 0, 0]],
                    [11, 12, [0, 0, 0]],
                    [12, 2, [0, 0, 0]]], 
                1:[[0, 1, [0, 0, 0]],
                    [1, 2, [0, 0, 0]],
                    [2, 3, [0, 0, 0]],
                    [3, 4, [0, 0, 0]],
                    [4, 0, [0, 1, 0]],
                    [2, 5, [0, 0, 0]],
                    [5, 6, [0, 0, 0]],
                    [7, 8, [0, 0, 0]],
                    [8, 2, [0, 0, 0]],
                    [9, 2, [0, 0, 0]],
                    [9, 10, [0, 0, 0]],
                    [11, 12, [0, 0, 0]],
                    [12, 2, [0, 0, 0]]],
                2:[[0, 1, [0, 0, 0]],
                    [1, 2, [0, 0, 0]],
                    [2, 3, [0, 0, 0]],
                    [3, 4, [0, 0, 0]],
                    [4, 0, [0, 1, 0]],
                    [2, 5, [0, 0, 0]],
                    [5, 6, [0, 0, 0]],
                    [6, 7, [1, 0, 0]],
                    [7, 8, [0, 0, 0]],
                    [8, 2, [0, 0, 0]],
                    [9, 2, [0, 0, 0]],
                    [9, 10, [0, 0, 0]],
                    [11, 12, [0, 0, 0]],
                    [12, 2, [0, 0, 0]]],
                3:[[0, 1, [0, 0, 0]],
                    [1, 2, [0, 0, 0]],
                    [2, 3, [0, 0, 0]],
                    [3, 4, [0, 0, 0]],
                    [4, 0, [0, 1, 0]],
                    [2, 5, [0, 0, 0]],
                    [5, 6, [0, 0, 0]],
                    [6, 7, [1, 0, 0]],
                    [7, 8, [0, 0, 0]],
                    [8, 2, [0, 0, 0]],
                    [9, 2, [0, 0, 0]],
                    [9, 10, [0, 0, 0]],
                    [10, 11, [0, 0, 1]],
                    [11, 12, [0, 0, 0]],
                    [12, 2, [0, 0, 0]]]}
    
    with open("testfile.txt", "w") as f:
        for key, edgeList in framesE.items():
            f.write(f"{key}\n")
            f.write(f"{nBonds[key]}\n")
            f.write(f"")
            for E in edgeList:
                f.write(f"{E[0]}\n")
                f.write(f"{E[1]}\n")
                for val in E[2]:
                    f.write(f"{val}\n")
        f.close()            
    system = percPerBC(3)
    dimsPercolated = []
    for val in range(0, 4):
        E = framesE.get(val)
        dimsPercolated.append(system.percolationDimension(V, E))
    return print(dimsPercolated)

testFrameByFrame()
