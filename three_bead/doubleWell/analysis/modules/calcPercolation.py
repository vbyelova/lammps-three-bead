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

def getPercDims(barrier, refoldBarrier, runNum, Vf, numMol):
    
    conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{Vf}_mol{numMol}"
    filename = f"Run{runNum}_{conditions}"
    percDims = []
    with open(f"../runs/{conditions}/timesteps.pkl", "rb") as f:
        timesteps = pickle.load(f)
    with open(f"../runs/{conditions}/{filename}/analysis/percDimsPerFrame.txt", "r") as f:
        for line in f:
            percDims.append(int(line))
    
    return percDims

def plotPercolation(barrier, refoldBarrier, runNum, Vf, numMol, percDims):
    conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{Vf}_mol{numMol}"
    filename = f"Run{runNum}_{conditions}"
    with open(f"../runs/{conditions}/{filename}/timesteps.pkl", "rb") as f:
        timesteps = pickle.load(f)
    counter = 0
    while counter < len(percDims):
        print(timesteps[counter], percDims[counter])
        counter += 1
    print(timesteps, percDims)
    plt.plot(timesteps, percDims)
    plt.show()
    return

def plotAvPercolation(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength):
    allPercolation = defaultdict(list)
    avPercolation = {}
    avPercolationErr = {}

    fig, ax = plt.subplots()
    for barrier in unfoldBarriers:
        with open(f"../runs/{conditions}/timesteps.pkl", "rb") as f:
            timesteps = pickle.load(f)
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}"
        frames = len(timesteps)
        for runNum in range(numRuns):
            filename = f"Run{runNum}_{conditions}"
            with open(f"../runs/{conditions}/{filename}/analysis/percDims.pkl", "rb") as f:
                percDims = pickle.load(f)
            
    
            if 3 in percDims:
                while len(percDims) < frames:
                    percDims.append(3)
            allPercolation[barrier].append(percDims)

        percolationArray = np.array(allPercolation[barrier])
        avPercolation[barrier] = percolationArray.mean(axis = 0)
        avPercolationErr[barrier] = percolationArray.std(axis = 0)          

        ax.errorbar(timesteps, avPercolation[barrier], yerr = avPercolationErr[barrier],
                    label = f"barrier = {barrier}kT")

    ax.set_xlabel("simulation frame")
    ax.set_ylabel("dimensions of percolation")
    ax.legend()
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/avPercolation_vf{vf}.png")
    plt.show()
    return avPercolation, avPercolationErr
