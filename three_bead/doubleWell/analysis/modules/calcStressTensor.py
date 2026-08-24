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

def calcStressTensor(barrier, refoldBarrier, runNum, vf, numMol, nBonds, bondInfo, boxLength, bondsPerAtom, prob, suffix):
    if bondsPerAtom == 2:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
    else:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
        if prob < 1:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}_prob{prob}{suffix}"
    filename = f"Run{runNum}_{conditions}"
    with open(f"../runs/{conditions}/timesteps.pkl", "rb") as f:
        stressTensor = np.zeros((3, 3))
    avStressTensors = []
    bondCounter = 0
    frame = 0
    vol = boxLength ** 3
    while frame < len(nBonds):
        forces = np.array([bondInfo[frame][bondCounter].properties[3],
                            bondInfo[frame][bondCounter].properties[4],
                            bondInfo[frame][bondCounter].properties[5]])
        vectors = np.array([bondInfo[frame][bondCounter].properties[6],
                            bondInfo[frame][bondCounter].properties[7],
                            bondInfo[frame][bondCounter].properties[8]])
        stressTensor += np.outer(vectors, forces)
        bondCounter += 1
        if bondCounter == nBonds[frame]:

            avStressTensors.append(- stressTensor / vol)
            stressTensor = np.zeros((3, 3))
            bondCounter = 0
            frame += 1
    return avStressTensors

def calcPressure(barrier, refoldBarrier, runNum, vf, numMol, nBonds, avStressTensors, bondsPerAtom, prob, suffix):

    if bondsPerAtom == 2:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
    else:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
        if prob < 1:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}_prob{prob}{suffix}"
    filename = f"Run{runNum}_{conditions}"
    with open(f"../runs/{conditions}/timesteps.pkl", "rb") as f:
        timesteps = pickle.load(f)
    pressure = []
    frames = [n for n in range(0, int(len(timesteps)))]
    for tensor in avStressTensors:
        pressure.append(np.trace(tensor) / 3)
    
    plt.plot(frames[:-2], pressure[:-2])
    plt.xlabel("simulation frame")
    plt.ylabel("pressure")
    plt.show()
    

def calcAvPressure(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, prob, suffixes):
    allPressure = defaultdict(list)
    avPressure = {}
    avPressureErr = {}

    for barrier, suffix in zip(unfoldBarriers, suffixes):
        if bondsPerAtom == 2:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
        else:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
            if prob < 1:
                conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}_prob{prob}{suffix}"
        with open(f"../runs/{conditions}/timesteps.pkl", "rb") as f:
            timesteps = pickle.load(f)
        frames = [n for n in range(0, int(len(timesteps)))]
        for runNum in range(numRuns):
            filename = f"Run{runNum}_{conditions}"
            
            with open(f"../runs/{conditions}/{filename}/analysis/nBonds.pkl", "rb") as f:
                nBonds = pickle.load(f)
            if len(nBonds) < 100: 
                print(f"skipping run{runNum} barrier{barrier}")
                continue
            with open(f"../runs/{conditions}/{filename}/analysis/bondInfo.pkl", "rb") as f:
                bondInfo = pickle.load(f)
            avStressTensors = calcStressTensor(barrier, refoldBarrier, runNum, vf, numMol,
                                               nBonds, bondInfo, boxLength, bondsPerAtom, prob, suffix)
            runPressure = [np.trace(tensor) / 3 for tensor in avStressTensors]
            allPressure[barrier].append(runPressure)
        
        pressureArray = np.array(allPressure[barrier])
        avPressure[barrier] = pressureArray.mean(axis = 0)
        avPressureErr[barrier] = pressureArray.std(axis = 0)

    return avPressure, avPressureErr

def plotAvPressure(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, prob, suffixes):
    
    colours = (["black", "dimgrey", "grey", "darkgrey", "silver"])
    percAt = []
    noperc = []
    perc = []
    with open(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/data/avPressure.pkl", "rb") as f:
        avPressure, avPressureErr = pickle.load(f)

    with open(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/data/avPercolation.pkl", "rb") as f:
        avPercolation, avPercolationErr = pickle.load(f)

    fig, ax = plt.subplots()

    for i, (barrier, suffix) in enumerate(zip(unfoldBarriers, suffixes)):
        if bondsPerAtom == 2:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
        else:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
            if prob < 1:
                conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}_prob{prob}{suffix}"
        with open(f"../runs/{conditions}/timesteps.pkl", "rb") as f:
            timesteps = pickle.load(f) 
        key = (barrier, suffix)
        if suffix != "":
            percAt.append(None)
            noperc.append((key))
        if suffix == "":
            perc.append(key)
            percAt.append(np.where(avPercolation[barrier] == 3)[0][0])

            ax.scatter(timesteps[percAt[i]], avPercolation[key][percAt[i]] / numMol * 100, color = "deeppink", marker = "x", zorder = 10)

    for i, (barrier, suffix) in enumerate(zip(unfoldBarriers, suffixes)):
        if bondsPerAtom == 2:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
        else:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
            if prob < 1:
                conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}_prob{prob}{suffix}"
        with open(f"../runs/{conditions}/timesteps.pkl", "rb") as f:
            timesteps = pickle.load(f)
        ax.errorbar(timesteps, avPressure[noperc[i]], yerr = avPressureErr[noperc[i]], color = colours[i],
                    label = f"barrier = {barrier}kT", linestyle = "dashdot")
        ax.errorbar(timesteps, avPressure[perc[i]], yerr = avPressureErr[perc[i]], color = colours[i],
                    label = f"barrier = {barrier}kT")

    ax.set_xlabel("timesteps")
    ax.set_ylabel("pressure")
    ax.legend()
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/pressure_vf{vf}.png")

    ax.set_xlabel("timesteps")
    ax.semilogx()
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/pressure_vf{vf}_semilog.png")
    plt.close()
    print("plotted average pressure..")
    return 

