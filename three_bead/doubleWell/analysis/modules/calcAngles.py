import numpy as np
import re
import matplotlib.pyplot as plt

from collections import defaultdict
from .parseDump import *

def calcAngles(barrier, refoldBarrier, runNum, vf, numMol, bondsPerAtom, suffix):
    """saves angles of three bead molecules in each simulation frame."""

    if bondsPerAtom == 2:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
    else:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
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

def unfoldedMolecules(barrier, refoldBarrier, runNum, vf, numMol, angles, bondsPerAtom, suffix):
    """makes a list of molecules that are unfolded (have an angle of 120-180)"""
    # let's say a particle is unfolded if it's around 120-180 degrees
    # based on our bimodal distribution
    if bondsPerAtom == 2:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
    else:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"

    # with (f"../runs/{conditions}/angleInfo.in", "r") as f:
    #     lines = f.readlines()
    #     info = lines[2].rsplit()
    #     maximum = float(info[2]) * 180 / np.pi
    unfoldedMols = defaultdict(list)
    for frame in range(len(angles)):
        for molID, angle in enumerate(angles[frame]):
            if angle >= 120:
                unfoldedMols[frame].append(molID)
            # if angle > maximum:
            #     unfoldedMols[frame].append(molID)

    print("generated list of unfolded molecules..")
    return unfoldedMols

def unfoldedOverTime(unfoldBarrier, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, suffixes):
    """calculates number of molecules unfolded over course of simulation."""
    allUnfold = defaultdict(list)
    avUnfold = {}
    avUnfoldErr = {}

    with open(f"../runs/boxLength{boxLength}/vf{vf}/data/avPercolation.pkl", "rb") as f:
        avPercolation, avPercolationErr = pickle.load(f)

    for barrier, suffix in zip(unfoldBarrier, suffixes):
        key = (barrier, suffix)
        if bondsPerAtom == 2:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
        else:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
        with open(f"../runs/{conditions}/timesteps.pkl", "rb") as f:
            timesteps = pickle.load(f)
        for runNum in range(numRuns):
            filename = f"Run{runNum}_{conditions}"
            with open(f"../runs/{conditions}/{filename}/analysis/nBonds.pkl", "rb") as f:
                nBonds = pickle.load(f)
            angles = calcAngles(barrier, refoldBarrier, runNum, vf, numMol, suffix)
            unfoldedMols = unfoldedMolecules(barrier, refoldBarrier, runNum, vf, numMol, angles)
            totalMols = [len(unfoldedMols[frame]) for frame in range(len(timesteps))]
            allUnfold[key].append(totalMols)

        avUnfoldArray = np.array(allUnfold[key])
        avUnfold[key] = avUnfoldArray.mean(axis = 0)
        avUnfoldErr[key] = avUnfoldArray.std(axis = 0)
    
    return avUnfold, avUnfoldErr

def plotUnfoldedOverTime(unfoldBarrier, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, suffixes):
    percAt = []
    fig, ax = plt.subplots()
    with open(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}data/avPercolation.pkl", "rb") as f:
        avPercolation, avPercolationErr = pickle.load(f)
    with open(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}data/unfoldedOverTime.pkl", "rb") as f:
        avUnfold, avUnfoldErr = pickle.load(f)
    for i, (barrier, suffix) in enumerate(zip(unfoldBarrier, suffixes)):
        key = (barrier, suffix)
        if bondsPerAtom == 2:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
        else:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
        with open(f"../runs/{conditions}/timesteps.pkl", "rb") as f:
            timesteps = pickle.load(f)
        ax.errorbar(timesteps, avUnfold[key], yerr = avUnfoldErr[key],
                    label = f"barrier = {barrier}kT")
        if suffix != "":
            percAt.append(None)
        else:
            percAt.append(np.where(avPercolation[barrier] == 3)[0][0])
            colour = ax.get_lines()[i].get_color()
            #ax.axvline(timesteps[percAt[i]], color = colour, linestyle = "--", alpha = 0.5)
            ax.scatter(timesteps[percAt[i]], avUnfold[key][percAt[i]], color = "black", marker = "*")
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
    return 

def bimodalAngle(barrier, refoldBarrier, runNum, vf, numMol, bondsPerAtom, suffix):
    """plots a histogram of the final angles of the molecules in the system."""

    if bondsPerAtom == 2:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
    else:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
    filename = f"Run{runNum}_{conditions}"
    angles = calcAngles(barrier, refoldBarrier, runNum, vf, numMol, bondsPerAtom, suffix)
    finalFrameAngles = angles[int(len(angles)-1)]
    hist, bins = np.histogram(finalFrameAngles, bins = np.linspace(0,200, 31))
    #logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))
    plt.hist(finalFrameAngles, bins = bins, color = "pink")
    plt.legend([f"vf = {vf}\nunfolding barrier  = {barrier} kT\n num. mol. = {numMol}"])
    plt.xlabel("molecule angle \u03B8")
    plt.ylabel("number of molecules")
    plt.title("Bimodal distribution of folded and unfolded molecules")
    plt.savefig(f"../runs/{conditions}/{filename}/analysis/figs/bimodalAngleDist")
    plt.close()
    print("plotted histogram of angle distribution in system..")
    return hist, bins, finalFrameAngles

def plotAvBimodalAngle(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, suffixes):
    allHistograms = defaultdict(list)
    histogramsMean = {}
    histogramsErr = {}
    fig, axs = plt.subplots(2, 3, sharex = True, sharey = True, figsize = (15, 8))
    labels = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)]
    
    plt.ylabel("count of molecules")
    plt.xlabel("molecule angle")
    for i, (barrier, suffix) in enumerate(zip(unfoldBarriers, suffixes)):
        ax = axs[labels[i]]
        if bondsPerAtom == 2:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
        else:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
        for runNum in range(numRuns):
            filename = f"Run{runNum}_{conditions}"
            with open(f"../runs/{conditions}/{filename}/analysis/angleDist.pkl", "rb") as f:
                hist, bins, finalFrameAngles = pickle.load(f)
            allHistograms[barrier].append(hist)
        key = (barrier, suffix)
        allHistogramsArray = np.array(allHistograms[barrier])
        histogramsMean[key] = np.mean(allHistogramsArray, axis = 0)
        histogramsErr[key] = np.std(allHistogramsArray, axis = 0)
        
        if suffix != "":
            simtype = f"\nno intermol."
        else:
            simtype = ""

        bin_centers = (bins[:-1] + bins[1:]) / 2
        ax.bar(bin_centers, histogramsMean[key], width=np.diff(bins),
            yerr=histogramsErr[key], color="pink",
            label=f"E_b = {barrier}kT {simtype}")
        ax.legend()

  
    plt.tight_layout()
    # for i, ax_row in enumerate(axs):
    #     for j, ax in enumerate(ax_row):
    #         print(f"subplot ({i},{j}): patches={len(ax.patches)}, xlim={ax.get_xlim()}, ylim={ax.get_ylim()}")
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}avBimodalAngle.png")
    return

def anglePopulation(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, suffixes):
    """plots a histogram of average distribution of folded vs unfolded molecules with
        energy barrier."""
    totalFolded = defaultdict(list)
    totalUnfolded = defaultdict(list)
    foldedMean = {}
    unfoldedMean = {}
    foldedError = []
    unfoldedError = []    
    for barrier, suffix in zip(unfoldBarriers, suffixes):
        if bondsPerAtom == 2:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
        else:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"

        with open(f"../runs/{conditions}/timesteps.pkl", "rb") as f:
            timesteps = pickle.load(f)
        finalframe = len(timesteps) - 1

        for runNum in range(numRuns):
            filename = f"Run{runNum}_{conditions}"
            with open(f"../runs/{conditions}/{filename}/analysis/nBonds.pkl", "rb") as f:
                nBonds = pickle.load(f)
            if len(nBonds) < 100:
                print(f"skipping run {runNum} barrier {barrier}")
                continue
            angles = calcAngles(barrier, refoldBarrier, runNum, vf, numMol, bondsPerAtom, suffix)
            unfoldedMol = unfoldedMolecules(barrier, refoldBarrier, runNum, vf, numMol, angles, bondsPerAtom, suffix)
        
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
    

    ax.bar(x, [foldedMean[b] for b in unfoldBarriers],
          width = barWidth, yerr = foldedError, color = "blue", label = "folded")
    ax.bar(x, [unfoldedMean[b] for b in unfoldBarriers],
          width = barWidth, yerr = unfoldedError, color = "red", label = "unfolded")
    
    ax.set_xticks(x)
    ax.set_xticklabels(unfoldBarriers)
    ax.set_xlabel("unfolding barrier (kT)")
    ax.set_ylabel("number of molecules")
    ax.legend()
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}barrierAnglePop_vf{vf}.png")
    plt.show()
    plt.close()
    print("plotted folded vs unfolded population by simulation end..")
    return

def avUnfoldPerFrame(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, suffixes):
    """plots a scatter graph of average unfolding events per timestep,
        compares against multiple unfolding barriers."""
    totalNewUnfoldedMol = defaultdict(list)
    avNewUnfoldedMol = {}
    avNewUnfoldedMolErr = {}
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
            nBonds = totalBonds(f"../runs/{conditions}/{filename}/output/nbonds.dat")
            #print(f"len nBonds {len(nBonds)}, run{runNum} barrier{barrier}")
            if len(nBonds) < 100: 
                print(f"skipping run{runNum} barrier{barrier}")
                continue
            angles = calcAngles(barrier, refoldBarrier, runNum, vf, numMol, suffix)
            unfoldedMols = unfoldedMolecules(barrier, refoldBarrier, runNum, vf, numMol, angles)
            
            newUnfoldedCount = [0]
            #print(nBonds)
            for frame in range(1, len(nBonds)):
                newUnfoldedCount.append(len(unfoldedMols[frame]) - len(unfoldedMols[frame - 1]))
            
            totalNewUnfoldedMol[barrier].append(newUnfoldedCount)
            #print(totalNewUnfoldedMol[barrier])
        
        newUnfoldedArray = np.array(totalNewUnfoldedMol[barrier])
        avNewUnfoldedMol[barrier] = newUnfoldedArray.mean(axis = 0)
        avNewUnfoldedMolErr[barrier] = newUnfoldedArray.std(axis = 0)
        ax.errorbar(timesteps, avNewUnfoldedMol[barrier], yerr = avNewUnfoldedMolErr[barrier],
                label = f"unfolding barrier = {barrier}kT {suffix}")

    ax.set_xlabel("simulation frame")
    ax.set_ylabel("number of new unfolding events")
    ax.legend()
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/unfoldingEvents_vf{vf}.png")

    print("plotted new unfolding events..")

    ax.semilogx()
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/unfoldingEvents_vf{vf}_semilog.png")
    plt.close()

    print("plotted unfolding events on semilog axis..")
    return avNewUnfoldedMol, avNewUnfoldedMolErr
