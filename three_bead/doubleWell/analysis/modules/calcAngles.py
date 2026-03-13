import numpy as np
import re
import matplotlib.pyplot as plt

from collections import defaultdict
from .parseDump import *

def calcAngles(barrier, refoldBarrier, runNum, vf, numMol):
    """saves angles of three bead molecules in each simulation frame."""

    conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}"
    filename = f"Run{runNum}_{conditions}"
    systemData = parseSystemData(f"../runs/{conditions}/{filename}/output/systemData.txt")
    boxLength, Npar, Nsteps, Nwrite, equilTime = systemData[0], systemData[1], systemData[2], systemData[3], systemData[4]
    
    angles = defaultdict(list)
    frame = 0
    counter = 0
    with open(f"../runs/{conditions}/{filename}/output/moleculeangles.dat", "r") as f:
        print("opened file")
        for line in f.readlines():
            line = line[:-1]
            #print(line)
            if re.search(r"\d+\s[\d.]+\s-?[\d.]+", line):
                #print(line, "matched")
                #angle = float(line.rsplit()[1]) * np.pi / 180
                angles[frame].append(float(line.rsplit()[1]))
                #print(line.rsplit()[1])
                counter += 1
                if counter == numMol:
                    #print("frame increased to ", frame)
                    frame += 1
                    counter = 0
    
    return angles

def unfoldedMolecules(barrier, refoldBarrier, runNum, vf, numMol, angles):
    """makes a list of molecules that are unfolded (have an angle of 120-180)"""
    # let's say a particle is unfolded if it's around 120-180 degrees
    # based on our bimodal distribution

    unfoldedMols = defaultdict(list)
    for frame in range(len(angles)):
        for molID, angle in enumerate(angles[frame]):
            if angle > 120:
                unfoldedMols[frame].append(molID)


    return unfoldedMols



def bimodalAngle(barrier, refoldBarrier, runNum, vf, numMol):
    """plots a histogram of the final angles of the molecules in the system."""

    conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}"
    filename = f"Run{runNum}_{conditions}"
    angles = calcAngles(barrier, refoldBarrier, runNum, vf, numMol)
    finalFrameAngles = angles[int(len(angles)-1)]
    hist, bins = np.histogram(finalFrameAngles, bins = 30)
    logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
    plt.hist(finalFrameAngles, bins = logbins, color = "pink")
    plt.legend([f"vf = {vf}\nunfolding barrier  = {barrier} kT\n num. mol. = {numMol}"])
    plt.xlabel("molecule angle \u03B8")
    plt.ylabel("log(number of molecules)")
    plt.title("Bimodal distribution of folded and unfolded molecules")
    plt.savefig(f"../runs/{conditions}/{filename}/analysis/figs/bimodalAngleDist")
    plt.close()
    return

def anglePopulation(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, timesteps):
    """plots a histogram of average distribution of folded vs unfolded molecules
        with time."""
    totalFolded = defaultdict(list)
    totalUnfolded = defaultdict(list)
    foldedMean = {}
    unfoldedMean = {}
    foldedError = []
    unfoldedError = []
    
    finalframe = len(timesteps) - 1
    
    for barrier in unfoldBarriers:
        for runNum in range(numRuns):
            angles = calcAngles(barrier, refoldBarrier, runNum, vf, numMol)
            unfoldedMol = unfoldedMolecules(barrier, refoldBarrier, runNum, vf, numMol, angles)
        
            numUnfolded = len(unfoldedMol[finalframe])
            numFolded = numMol - numUnfolded
            
            totalFolded[barrier].append(numFolded)
            totalUnfolded[barrier].append(numUnfolded)
            
        foldedArray = np.array(totalFolded[barrier])
        unfoldedArray = np.array(totalUnfolded[barrier])
            
        foldedMean[barrier] = foldedArray.mean()
        unfoldedMean[barrier] = unfoldedArray.mean()
        foldedError.append(foldedArray.std())
        unfoldedError.append(unfoldedArray.std())
        
    barWidth = 0.35
    x = np.arange(len(unfoldBarriers))
    
    fig, ax = plt.subplots()
    
    ax.bar(x - barWidth / 2, [foldedMean[b] for b in unfoldBarriers],
           width = barWidth, yerr = foldedError, color = "blue", label = "folded")
    ax.bar(x + barWidth / 2, [unfoldedMean[b] for b in unfoldBarriers],
           width = barWidth, yerr = unfoldedError, color = "red", label = "unfolded")
    
    ax.set_xticks(x)
    ax.set_xticklabels(unfoldBarriers)
    ax.set_xlabel("unfolding barrier (kT)")
    ax.set_ylabel("number of molecules")
    ax.legend()
    plt.savefig(f"./barrierAnglePop")
    plt.close()
    return

def avUnfoldPerFrame(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, timesteps):
    """plots a scatter graph of average unfolding events per timestep,
        compares against multiple unfolding barriers."""
    totalNewUnfoldedMol = defaultdict(list)
    avNewUnfoldedMol = {}
    avNewUnfoldedMolErr = {}
    fix, ax = plt.subplots()
    for barrier in unfoldBarriers:
        for runNum in range(numRuns):
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}"
            filename = f"Run{runNum}_{conditions}"
            nBonds = totalBonds(f"../runs/{conditions}/{filename}/output/nbonds.dat")
            
            newBonds = [0]
            for frame in range(1, len(nBonds)):
                newBonds.append(nBonds[frame] - nBonds[frame - 1])
            
            totalNewUnfoldedMol[barrier].append(newBonds)
            print(totalNewUnfoldedMol[barrier])
        
        newUnfoldedArray = np.array(totalNewUnfoldedMol[barrier])
        avNewUnfoldedMol[barrier] = newUnfoldedArray.mean(axis = 0)
        avNewUnfoldedMolErr[barrier] = newUnfoldedArray.std(axis = 0) / np.sqrt(numRuns)
        ax.errorbar(timesteps, avNewUnfoldedMol[barrier], yerr = avNewUnfoldedMolErr[barrier],
                label = f"unfolding barrier = {barrier}kT")
        
    ax.set_xlabel("simulation frame")
    ax.set_ylabel("number of new unfolding events")
    ax.legend()
    plt.savefig("./unfoldingEvents")
    plt.close()
    return