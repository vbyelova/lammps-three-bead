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
            if re.search(r"\d+\s[\d.]+\s-?[\d.]+", line):
                angles[frame].append(float(line.rsplit()[1]))
                counter += 1
                if counter == numMol:
                    frame += 1
                    counter = 0
    
    print("parsed angles..")
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

    print("generated list of unfolded molecules..")
    return unfoldedMols

def unfoldedOverTime(unfoldBarrier, refoldBarrier, numRuns, vf, numMol, boxLength, timesteps):
    """calculates number of molecules unfolded over course of simulation."""
    allUnfold = defaultdict(list)
    avUnfold = {}
    avUnfoldErr = {}

    fig, ax = plt.subplots()
    for barrier in unfoldBarrier:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}"
        for runNum in range(numRuns):
            filename = f"Run{runNum}_{conditions}"

            angles = calcAngles(barrier, refoldBarrier, runNum, vf, numMol)
            unfoldedMols = unfoldedMolecules(barrier, refoldBarrier, runNum, vf, numMol, angles)
            totalMols = [len(unfoldedMols[frame]) for frame in range(len(timesteps))]
            allUnfold[barrier].append(totalMols)

        avUnfoldArray = np.array(allUnfold[barrier])
        avUnfold[barrier] = avUnfoldArray.mean(axis = 0)
        avUnfoldErr[barrier] = avUnfoldArray.std(axis = 0)
    
        ax.errorbar(timesteps, avUnfold[barrier], yerr = avUnfoldErr[barrier],
                     label = f"barrier = {barrier}kT")
    ax.set_xlabel("simulation frame")
    ax.set_ylabel("number of unfolded molecules")
    ax.set_title("cumulative total of unfolded mol over simulation")
    ax.legend()
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/unfoldedMol_vf{vf}.png")
    ax.semilogx()
    ax.set_xlabel("simulation frame (semilog axis)")
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/unfoldedMol_vf{vf}_semilog.png")
    plt.close()

    print("plotted graphs of total unfolded molecules over time..")
    return avUnfold, avUnfoldErr

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
    print("plotted histogram of angle distribution in system..")
    return

def anglePopulation(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, timesteps):
    """plots a histogram of average distribution of folded vs unfolded molecules with
        energy barrier."""
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
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/barrierAnglePop_vf{vf}.png")
    plt.close()
    print("plotted folded vs unfolded population by simulation end..")
    return

def avUnfoldPerFrame(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength):
    """plots a scatter graph of average unfolding events per timestep,
        compares against multiple unfolding barriers."""
    totalNewUnfoldedMol = defaultdict(list)
    avNewUnfoldedMol = {}
    avNewUnfoldedMolErr = {}
    fix, ax = plt.subplots()
    for barrier in unfoldBarriers:
        with open(f"../runs/{conditions}/timesteps.pkl", "rb") as f:
            timesteps = pickle.load(f)
        for runNum in range(numRuns):
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}"
            filename = f"Run{runNum}_{conditions}"
            nBonds = totalBonds(f"../runs/{conditions}/{filename}/output/nbonds.dat")
            angles = calcAngles(barrier, refoldBarrier, runNum, vf, numMol)
            unfoldedMols = unfoldedMolecules(barrier, refoldBarrier, runNum, vf, numMol, angles)
            
            newUnfoldedCount = [0]
            print(nBonds)
            for frame in range(1, len(nBonds)):
                newUnfoldedCount.append(len(unfoldedMols[frame]) - len(unfoldedMols[frame - 1]))
            
            totalNewUnfoldedMol[barrier].append(newUnfoldedCount)
            print(totalNewUnfoldedMol[barrier])
        
        newUnfoldedArray = np.array(totalNewUnfoldedMol[barrier])
        avNewUnfoldedMol[barrier] = newUnfoldedArray.mean(axis = 0)
        avNewUnfoldedMolErr[barrier] = newUnfoldedArray.std(axis = 0)
        ax.errorbar(timesteps, avNewUnfoldedMol[barrier], yerr = avNewUnfoldedMolErr[barrier],
                label = f"unfolding barrier = {barrier}kT")

    ax.set_xlabel("simulation frame")
    ax.set_ylabel("number of new unfolding events")
    ax.legend()
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/unfoldingEvents_vf{vf}.png")

    print("plotted new unfolding events..")

    ax.semilogx()
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/unfoldingEvents_vf{vf}_semilog.png")
    plt.close()

    print("plotted unfolding events on semilog axis..")
    return avNewUnfoldedMol, avNewUnfoldedMolErr
