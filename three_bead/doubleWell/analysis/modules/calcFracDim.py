# a python script to calculate the fractal dimension of a network consisting
# of three-bead molecules.
# Victoria Byelova

import matplotlib.pyplot as plt
import numpy as np
import pickle
import math
import pwlf


from collections import defaultdict
from .parseDump import *
from .calcPercolation import *

def boxCounting(filename, boxLength, timesteps, percDims, atPerc):
    """A function to count the number of particles in given boxLength intervals. First the number of voxels and
        their sizes is decided. Then bins are made to represent the voxels and digitize decides which voxel
        each particle belongs to. The unique voxels are counted and added to a list."""
    
    with open(filename, "rb") as f:
        particles = pickle.load(f)

    particleSize = 2**(1/6)
    halfLength = 0.5 * boxLength
    
    totalUniqueVoxels = []

    voxelSizes = np.logspace(np.log10(particleSize), np.log10(halfLength), num = 10)
    numVoxels = [int(boxLength/vSize) for vSize in voxelSizes]
    voxelSizes = [boxLength/ numDivs for numDivs in numVoxels]
    print("initialised arrays for box counting..")
    print(f"currently in {filename}")
    print(percDims)

    if atPerc == True:
        percAtFrame = percDims.index(3)
    elif atPerc == False:
        percAtFrame = particles[0].properties.shape[0] - 1
        print("false, num frames is ", percAtFrame)

    print(f"currently in {filename}")
    for vNum, v in enumerate(numVoxels):
        bins = np.linspace(- halfLength, halfLength, v + 1)
      
        uniqueVoxels = set()
        print("made bins for boxLength size", vNum + 1 )
        
        for pNum, p in enumerate(particles):
            x = p.properties[percAtFrame, 1]
            y = p.properties[percAtFrame, 2]
            z = p.properties[percAtFrame, 3]

            xVox = np.digitize(x, bins) - 1   #digitize returns 1-based, so -1
            yVox = np.digitize(y, bins) - 1
            zVox = np.digitize(z, bins) - 1

            uniqueVoxels.add((xVox, yVox, zVox))
        
        totalUniqueVoxels.append([voxelSizes[vNum], len(uniqueVoxels)])
    print("voxels size and number of particles per voxel: ", totalUniqueVoxels)
    return np.array(totalUniqueVoxels)
        

def calcFractalDimension(unfoldBarrier, refoldBarrier, runNum, vf, numMol, totalUniqueVoxels, boxLength, percDims, atPerc, bondsPerAtom, suffix):
    """calculates the fractal dimension for a single run"""
    if bondsPerAtom == 2:
        conditions = f"unfold{unfoldBarrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
    else:
        conditions = f"unfold{unfoldBarrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
    filename = f"Run{runNum}_{conditions}"
    with open(f"../runs/{conditions}/timesteps.pkl", "rb") as f:
        timesteps = pickle.load(f)
    if atPerc == True:
        atPercString = "AtPercolation"
    elif atPerc == False:
        atPercString = "SimEnd"
    logNumUniqueVoxels = np.array([np.log(i) for i in totalUniqueVoxels[:,1]])
    logInverseVoxelSizes = np.array([np.log(1 / i) for i in totalUniqueVoxels[:,0]])

    plt.plot(logInverseVoxelSizes, logNumUniqueVoxels, '*', color = "hotpink")
    plt.legend([f"vf = {vf}\nunfolding barrier  = {unfoldBarrier}kT\n n. molecules = {numMol}"])
    plt.title("Finding fractal dimension by boxLength counting method")
    plt.xlabel("log(1 / R)")
    plt.ylabel("log( N )")
    plt.title(f"box counting {atPercString}")


    # if boxLength > 30:
    model = pwlf.PiecewiseLinFit(logInverseVoxelSizes, logNumUniqueVoxels)
    breaks = model.fit_guess([-1.5])

    xHat = np.linspace(min(logInverseVoxelSizes), max(logInverseVoxelSizes), 100)
    yHat = model.predict(xHat)
    plt.plot(xHat, yHat, "-", color = "purple")
    plt.xlim(math.floor(min(logInverseVoxelSizes)), int(max(logInverseVoxelSizes)))
    plt.tight_layout()
    plt.savefig(f"../runs/{conditions}/{filename}/analysis/figs/boxCounting{atPercString}")
    plt.close()
    return 2, np.abs(model.slopes[1]), breaks[1]

    # if boxLength < 30:
    #     trend = np.polyfit(logInverseVoxelSizes, logNumUniqueVoxels, 1)
    #     trendpoly = np.poly1d(trend)
    #     plt.plot(logInverseVoxelSizes, trendpoly(logInverseVoxelSizes), color = "purple")
    #     plt.tight_layout()
    #     plt.savefig(f"../runs/{conditions}/{filename}/analysis/figs/boxCounting")
    #     plt.close()
    #     return 1, trend[0]

def findAllFractalDims(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, atPerc, bondsPerAtom, suffixes):
    """finds the fractal dimension for each run for each unfolding barrier and also calculates the average
        correlation length. Can find this for either at percolation point or in final frame of sim."""

    fractalDim = defaultdict(list)
    fractalDimError = defaultdict(list)

    fractalDimMean = {}
    finalFractalDims = []
    finalFractalDimsError = []

    corrLength = defaultdict(list)
    avCorrLength = {}
    avCorrLengthErr = {}

    for barrier, suffix in zip(unfoldBarriers, suffixes):
        if suffix != "":
            continue
        if bondsPerAtom == 2:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
        else:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
        with open(f"../runs/{conditions}/timesteps.pkl", "rb") as t:
            timesteps = pickle.load(t)
        validRuns = 0
        if atPerc == True:
            atPercString = "AtPerc"
        elif atPerc == False:
            atPercString = "SimEnd"
        with open(f"../runs/{conditions}/data/corrLengths{atPercString}.txt", "w") as f:
            f.write(f"# correlation length\t# standard deviation\n")

            for runNum in range(numRuns):
                filename = f"Run{runNum}_{conditions}"
                percDims = getPercDims(barrier, refoldBarrier, runNum, vf, numMol)
                if atPerc == True and 3 not in percDims:
                    print(f"{conditions} run {runNum} did not percolate, skipping" )
                    with open(f"../runs/{conditions}/analysisNotes.txt", "w") as s:
                        s.write(f"{conditions} run {runNum} did not percolate, skipped during fractal dim "
                                f"analysis.\n Maximum percolation dimension = {max(percDims)}\n")
                        s.close()
                    continue
                validRuns += 1
                totalUniqueVoxels = boxCounting(f"../runs/{conditions}/{filename}/analysis/particleTraj.pkl",
                                                boxLength, timesteps, percDims, atPerc)
                dimsAndBreak = calcFractalDimension(barrier, refoldBarrier, runNum, vf, numMol, 
                                                    totalUniqueVoxels, boxLength, percDims, atPerc, bondsPerAtom, suffix)
                if dimsAndBreak[0] == 2:
                    df = dimsAndBreak[1]
                    breakpoint = dimsAndBreak[2]
                    fractalDim[barrier].append(df)
                    print("breakpoint = ", breakpoint)
                    
                    corrLength[barrier].append(10 ** -breakpoint)

                if dimsAndBreak[0] == 1:
                    fractalDim[barrier].append(dimsAndBreak[1])

            fractalDimMean[barrier] = sum(fractalDim[barrier]) / validRuns
            for val in fractalDim[barrier]:
                fractalDimError[barrier].append((val - fractalDimMean[barrier])**2)
        
            finalFractalDimsError.append(np.sqrt(sum(fractalDimError[barrier]) / validRuns))  
            finalFractalDims.append(sum(fractalDim[barrier]) / validRuns)

            print(f"fractal dimensions for unfolding barrier = {barrier}: {fractalDim[barrier]}")
            
            corrLengthArray = np.array(corrLength[barrier])
            avCorrLength[barrier] = corrLengthArray.mean(axis = 0)
            avCorrLengthErr[barrier] = corrLengthArray.std(axis = 0)
                
            f.write(f"{avCorrLength[barrier]}\t{avCorrLengthErr[barrier]}\n")

    return finalFractalDims, finalFractalDimsError, avCorrLength, avCorrLengthErr

def plotAverageFractalDim(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom):

    """plots a scatter graph of what the average fractal dimension is with each unfolding barrier
        at percolation and at the end of the simulation."""
    with open(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/data/atPercFractalDims.pkl", "rb") as f:
        atPercFinalFractalDims, atPercFinalFractalDimsError, atPercAvCorrLength, atPercAvCorrLengthErr =  pickle.load(f)

    with open(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/data/simEndFractalDims.pkl", "rb") as f:
        simEndFinalFractalDims, simEndFinalFractalDimsError, simEndAvCorrLength, simEndAvCorrLengthErr =  pickle.load(f)  

    barWidth = 0.35
    x = np.arange(len(unfoldBarriers) - 1)
    fig, ax = plt.subplots()
    unfoldBarriers = sorted(set(unfoldBarriers))
    ax.errorbar(unfoldBarriers, atPercFinalFractalDims,
            yerr = atPercFinalFractalDimsError, label = "at gelation", marker = "x", color = "lightskyblue")
    ax.errorbar(unfoldBarriers, simEndFinalFractalDims,
            yerr = simEndFinalFractalDimsError, label = "at end of simulation", marker = "x", color = "steelblue")
    #ax.set_xticks(x)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
          fancybox=True, shadow=True, ncol=2, fontsize = 12)
    #ax.set_xticklabels(unfoldBarriers)
    ax.set_xlabel(f"$E_U$ (kT)", fontsize = 15)
    ax.set_ylabel(f"Fractal dimension $D_f$", fontsize = 15)
    plt.tight_layout()
    #ax.set_title("Fractal dimension through box counting method")
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/fractalDim_vf{vf}.png")
    plt.close()

    fig, ax = plt.subplots()
    ax.errorbar(unfoldBarriers, [atPercAvCorrLength[b] for b in unfoldBarriers],
            yerr = [atPercAvCorrLengthErr[b] for b in unfoldBarriers], color = "mediumseagreen", marker = "x", label = "at gelation")
    ax.errorbar(unfoldBarriers, [simEndAvCorrLength[b] for b in unfoldBarriers],
            yerr = [simEndAvCorrLengthErr[b] for b in unfoldBarriers], color = "darkgreen", marker = "x", label = "at end of simulation")
    
    #ax.set_xticks(x)
    #ax.set_xticklabels(unfoldBarriers)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
          fancybox=True, shadow=True, ncol=2, fontsize = 12)
    ax.set_xlabel("$E_U$ (kT)", fontsize = 15)
    ax.set_ylabel(r"correlation length $\xi$", fontsize = 15)
    plt.tight_layout()
    #ax.set_title("Cluster sizes at different simulation points")
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/corrLength_vf{vf}.png")
    plt.close()

    return (atPercFinalFractalDims, atPercFinalFractalDimsError, simEndFinalFractalDims, simEndFinalFractalDimsError,
            atPercAvCorrLength, atPercAvCorrLengthErr, simEndAvCorrLength, simEndAvCorrLengthErr)
