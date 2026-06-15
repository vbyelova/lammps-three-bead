import matplotlib.pyplot as plt
import numpy as np
import pickle
from numba import jit
from collections import defaultdict

from .parseDump import *

@jit
def RDFnumba(positions, boxLength, dr):
    numBins = int(0.5 * boxLength / dr)
    rdfSet = np.zeros(numBins)
    for num1 in range(len(positions)):
        for num2 in range(len(positions)):
            # check if rdf is between central molecules
            if (num1 // 3 == (num1) // 3 and num1 // 3 ==(num1 - 1) // 3
                and num2 // 3 == (num2) // 3 and num2 // 3 ==(num2 - 1) // 3):
                continue
            if num1 == num2:
                continue
            dx = positions[num2, 0] - positions[num1, 0]
            dy = positions[num2, 1] - positions[num1, 1] 
            dz = positions[num2, 2] - positions[num1, 2]

            dx = dx - boxLength * np.round(dx / boxLength)
            dy = dy - boxLength * np.round(dy / boxLength)
            dz = dz - boxLength * np.round(dz / boxLength)
            sep = np.sqrt(dx**2 + dy**2 + dz**2)
            shellIndex = int(sep / dr)
            if shellIndex < numBins:
                rdfSet[shellIndex] += 1
                
    return rdfSet        

def posRDF(barrier, refoldBarrier, runNum, vf, numMol, boxLength, particles, bondsPerAtom, suffix):
    """calculates the positional radial distribution function. modified 
        so that the central particle of each molecule only is taken into account"""
    if bondsPerAtom == 2:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
    else:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
    filename = f"Run{runNum}_{conditions}"
    print(f"currently on run {runNum} for barrier {barrier}")
    with open(f"../runs/{conditions}/timesteps.pkl", "rb") as f:
        timesteps = pickle.load(f)
    parRadius = 2 ** (1/6)
    boxVol = boxLength ** 3
    dr = 0.005 * boxLength
    shellBounds = np.arange(0, 0.5 * boxLength, dr)
    frames = particles[0].properties.shape[0]
    frame = frames - 1
    rdfSet = np.zeros(len(shellBounds), dtype = np.float64)
    print("generated vals for final frame rdf..")

    positions = np.array([[p.properties[frame, 1],
                p.properties[frame, 2],
                p.properties[frame, 3]] for p in particles], dtype=np.float64)
    
    rdfSet = RDFnumba(positions, boxLength, dr)
    print(rdfSet)

    for i in range(len(shellBounds)):
        innerRadius = shellBounds[i]
        outerRadius = innerRadius + dr
        sliceVol = ((4/3) * np.pi * (outerRadius**3 - innerRadius**3))
        print(f"slice vol {sliceVol}")
        if sliceVol > 0:
            rdfSet[i] /= sliceVol
            rdfSet[i] /= (numMol)**2 / boxVol
    print("calculated rdf..")
    print(f"rdf set{i} {rdfSet[i]}")
    # fig, ax = plt.subplots()
    # ax.scatter(shellBounds, rdfSet)
    # ax.set_xlabel("r")
    # ax.set_ylabel("g(r)")
    # plt.savefig("./rdftest")
    # plt.close()
    # print("saved rdf plot..")
    return rdfSet, shellBounds

def calcPosRDF(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, suffixes):
    allRDFsets = defaultdict(list)
    avRDFs = {}
    avRDFsErr = {}
    fig, ax = plt.subplots()
    for barrier, suffix in zip(unfoldBarriers, suffixes):
        if bondsPerAtom == 2:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
        else:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
        with open(f"../runs/{conditions}/timesteps.pkl", "rb") as f:
            timesteps = pickle.load(f)
        for runNum in range(numRuns):
            filename = f"Run{runNum}_{conditions}"
            with open(f"../runs/{conditions}/{filename}/analysis/particleTraj.pkl", "rb") as f:
                particles = pickle.load(f)

            rdfSet, shellBounds = posRDF(barrier, refoldBarrier, runNum, vf, numMol, boxLength,
                            particles, bondsPerAtom, suffix)
            allRDFsets[barrier].append(rdfSet)

        rdfArray = np.array(allRDFsets[barrier])
        avRDFs[barrier] = rdfArray.mean(axis = 0)
        avRDFsErr[barrier] = rdfArray.std(axis = 0)
    return avRDFs, avRDFsErr, shellBounds

def plotAvRDF(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, suffixes):
    with open(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/data/avRDF.pkl", "rb") as f:
        avRDFs, avRDFsErr, shellBounds = pickle.load(f)
        
    fig, ax = plt.subplots()
    for barrier in unfoldBarriers:
        ax.errorbar(shellBounds, avRDFs[barrier], yerr = avRDFsErr[barrier],
                    label = f"barrier = {barrier}kT")
        
    ax.set_xlabel("r")
    ax.set_ylabel("g(r)")
    ax.legend()
    ax.set_title(f"RDF for vf = {vf}")
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/rdf{vf}.png")
    plt.show()
    print(f"plotted rdf for vf = {vf}..")
    return avRDFs, avRDFsErr, shellBounds
