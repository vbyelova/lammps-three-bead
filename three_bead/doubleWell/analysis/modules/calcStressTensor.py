import numpy as np
import re
import matplotlib.pyplot as plt

from collections import defaultdict
from .parseDump import *
from .threeBeadClasses import *

def plotForceMagnitude(barrier, refoldBarrier, runNum, vf, numMol, nBonds, bondInfo):
    """plots histograms of forces in the final frame of the simulation
        in each direction as well as a histogram of the force magnitudes."""
    conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}"
    filename = f"Run{runNum}_{conditions}"
    forceMagnitudes = []
    totalfx = []
    totalfy = []
    totalfz = []

    # properties[0] is bond index
    # properties[1] and [2] are atom indices
    # properties[3], [4] and [5] are fx, fy, fz
    # properties[6], [7] and [8] are dx, dy, dz
    # properties[9] is bond length

    finalFrame = len(nBonds) - 1
    for b in range(0, int(nBonds[finalFrame])):
        fx = bondInfo[finalFrame][b].properties[3]
        fy = bondInfo[finalFrame][b].properties[4]
        fz = bondInfo[finalFrame][b].properties[5]

        forceMag = np.sqrt(fx**2 + fy**2 + fz**2)
        forceMagnitudes.append(forceMag)
        totalfx.append(fx)
        totalfy.append(fy)
        totalfz.append(fz)

    fig, ax = plt.subplots(2, 2)
    fxhist = ax[0, 0].hist(totalfx, 20, color = "pink")
    fyhist = ax[1, 1].hist(totalfy, 20, color = "palevioletred")
    fzhist = ax[1, 0].hist(totalfz, 20, color = "mediumvioletred")
    maghist = ax[0, 1].hist(forceMagnitudes, 20, color = "rebeccapurple")
    ax[1, 1].set_xlabel("force")
    ax[0, 0].set_ylabel("number of bonds")
    ax[0, 1].set_xlabel("force magnitudes")
    ax[0, 1].semilogy()
    plt.tight_layout()
    plt.show()

def calcStressTensor(barrier, refoldBarrier, runNum, vf, numMol, nBonds, bondInfo, boxLength, timesteps):
    stressTensor = np.zeros((3, 3))
    avStressTensors = []
    bondCounter = 0
    frame = 0
    vol = boxLength ** 3
    while frame < len(timesteps):
        forces = np.array([bondInfo[frame][bondCounter].properties[3],
                            bondInfo[frame][bondCounter].properties[4],
                            bondInfo[frame][bondCounter].properties[5]])
        vectors = np.array([bondInfo[frame][bondCounter].properties[6],
                            bondInfo[frame][bondCounter].properties[7],
                            bondInfo[frame][bondCounter].properties[8]])
        stressTensor += np.outer(vectors, forces)
        bondCounter += 1
        if bondCounter == nBonds[frame]:
            stressTensor = - stressTensor / vol
            avStressTensors.append(stressTensor)
            bondCounter = 0
            frame += 1
    return avStressTensors

def calcPressure(barrier, refoldBarrier, runNum, vf, numMol, nBonds, avStressTensors, timesteps):
    pressure = []

    frames = [n for n in range(0, int(len(timesteps)))]
    for tensor in avStressTensors:
        pressure.append(np.trace(tensor) / 3)
    
    plt.plot(frames[:-2], pressure[:-2])
    plt.xlabel("simulation frame")
    plt.ylabel("pressure")
    plt.show()
    

def plotAvPressure(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, nBonds, avStressTensors, timesteps):
    pressure = defaultdict(list)
    pressureErr = defaultdict(list)
    avPressure = {}
    avPressureErr = []
    for barrier in unfoldBarriers:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}"

        for runNum in range(numRuns):
            filename = f"Run{runNum}_{conditions}"
            frames = [n for n in range(0, int(len(timesteps)))]
            for tensor in avStressTensors:
                press = np.trace(tensor) / 3
                pressure[barrier].append(press)
        
        avPressure[barrier] = sum(pressure[barrier]) / numRuns
        for val in avPressure[barrier]:
            pressureErr[barrier].append(val - avPressure[barrier] ** 2)


    plt.plot
