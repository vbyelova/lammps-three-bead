import matplotlib.pyplot as plt
import numpy as np
import pickle
from numba import jit
from collections import defaultdict

from .parseDump import *

def parseLammpsPosRDF(barrier, refoldBarrier, runNum, Vf, numMol):

    conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{Vf}_mol{numMol}"
    filename = f"Run{runNum}_{conditions}"
    rdf = defaultdict(list)
    smolrdf = defaultdict(list)
    with open(f"../runs/{conditions}/{filename}/output/rdf.dat", "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rsplit()
            if len(parts) == 2:
                frame = int(parts[0])
            if len(parts) == 4:
                distance = float(parts[1])
                gofr = float(parts[2])
                smolrdf[frame].append([distance, gofr])
    with open(f"../runs/{conditions}/{filename}/output/bigRDF.dat", "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rsplit()
            if len(parts) == 2:
                frame = int(parts[0])
            if len(parts) == 4:
                distance = float(parts[1])
                gofr = float(parts[2])
                rdf[frame].append([distance, gofr])

    plt.plot(np.array(rdf[frame])[:, 0], np.array(rdf[frame])[:, 1], "-", color = "blue")
    plt.plot(np.array(smolrdf[frame])[:, 0], np.array(smolrdf[frame])[:, 1], "-", color = "red")
    print(np.array(rdf[frame]))
    plt.xlabel("distance")
    plt.ylabel("g(r)")
    plt.show()
    return

@jit
def RDFnumba(positions, boxLength, dr):
    numPar = len(positions)
    numBins = int(0.5 * boxLength / dr)
    rdfSet = np.zeros(numBins)
    for num1 in range(len(positions)):
        for num2 in range(len(positions)):
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

def posRDF(barrier, refoldBarrier, runNum, Vf, numMol, boxLength, particles):
    """calculates the positional radial distribution function"""

    conditions = f"unfold{barrier}_refold{2}_Vf{Vf}_mol{numMol}"
    filename = f"Run{runNum}_{conditions}"

    with open(f"../runs/{conditions}/timesteps.pkl", "rb") as f:
        timesteps = pickle.load(f)
    parRadius = 2 ** (1/6)
    boxVol = boxLength ** 3
    dr = 0.005 * boxLength
    shellBounds = np.arange(0, 0.5 * boxLength, dr)
    frame = len(timesteps) - 1
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
            rdfSet[i] /= (numMol * 3)**2 / boxVol
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

def avPosRDF(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, timesteps):
    allRDFsets = defaultdict(list)
    avRDFs = {}
    avRDFsErr = {}
    fig, ax = plt.subplots()
    for barrier in unfoldBarriers:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}"
        for runNum in range(numRuns):
            filename = f"Run{runNum}_{conditions}"
            with open(f"../runs/{conditions}/{filename}/analysis/particleTraj.pkl", "rb") as f:
                particles = pickle.load(f)

            rdfSet, shellBounds = posRDF(barrier, refoldBarrier, runNum, vf, numMol, boxLength,
                            particles, timesteps)
            allRDFsets[barrier].append(rdfSet)

        rdfArray = np.array(allRDFsets[barrier])
        avRDFs[barrier] = rdfArray.mean(axis = 0)
        avRDFsErr[barrier] = rdfArray.std(axis = 0)
        #print(f"shellbounds len {len(shellBounds)}, rdf stuff len {avRDFs[barrier]}")
        ax.errorbar(shellBounds, avRDFs[barrier], yerr = avRDFsErr[barrier],
                    label = f"barrier = {barrier}kT")
        
    ax.set_xlabel("r")
    ax.set_ylabel("g(r)")
    ax.legend()
    ax.set_title(f"RDF for vf = {vf}")
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/rdf{vf}.png")
    plt.show()
    print(f"plotted rdf for vf = {vf}..")
    return avRDFs, avRDFsErr, shellBounds
