import numpy as np
import re
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as mcolors
from itertools import cycle
from collections import defaultdict
from .parseDump import *
from scipy.optimize import curve_fit

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

    with open(f"../runs/{conditions}/angleInfo.in", "r") as f:
        lines = f.readlines()
        info = lines[2].rsplit()
        maximum = float(info[2]) * 180 / np.pi
    unfoldedMols = defaultdict(list)
    for frame in range(len(angles)):
        for molID, angle in enumerate(angles[frame]):
            # if angle >= 120:
            #    unfoldedMols[frame].append(molID)
             if angle > maximum:
                 unfoldedMols[frame].append(molID)

    print("generated list of unfolded molecules..")
    return unfoldedMols

def unfoldedOverTime(unfoldBarrier, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, suffixes):
    """calculates number of molecules unfolded over course of simulation."""
    allUnfold = defaultdict(list)
    avUnfold = {}
    avUnfoldErr = {}

    with open(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/data/avPercolation.pkl", "rb") as f:
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
            angles = calcAngles(barrier, refoldBarrier, runNum, vf, numMol, bondsPerAtom, suffix)
            unfoldedMols = unfoldedMolecules(barrier, refoldBarrier, runNum, vf, numMol, angles, bondsPerAtom, suffix)
            totalMols = [len(unfoldedMols[frame]) for frame in range(len(timesteps))]
            allUnfold[key].append(totalMols)

        avUnfoldArray = np.array(allUnfold[key])
        avUnfold[key] = avUnfoldArray.mean(axis = 0)
        avUnfoldErr[key] = avUnfoldArray.std(axis = 0)
    
    return avUnfold, avUnfoldErr

def plotUnfoldedOverTime(unfoldBarrier, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, suffixes):
    percAt = []
    perc = []
    noperc = []

    #colours = ["maroon", "firebrick", "indianred", "lightcoral", "lightpink"]
    colours = (["black", "dimgrey", "grey", "darkgrey", "silver"])

    fig, ax = plt.subplots()
    with open(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/data/avPercolation.pkl", "rb") as f:
        avPercolation, avPercolationErr = pickle.load(f)
    with open(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/data/unfoldedOverTime.pkl", "rb") as f:
        avUnfold, avUnfoldErr = pickle.load(f)

        
    for i, (barrier, suffix) in enumerate(zip(unfoldBarrier, suffixes)):
        if bondsPerAtom == 2:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
        else:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
        with open(f"../runs/{conditions}/timesteps.pkl", "rb") as f:
            timesteps = pickle.load(f) 
        key = (barrier, suffix)
        if suffix != "":
            percAt.append(None)
            noperc.append((key))
        if suffix == "":
            perc.append(key)
            percAt.append(np.where(avPercolation[barrier] == 3)[0][0])

            ax.scatter(timesteps[percAt[i]], avUnfold[key][percAt[i]] / numMol * 100, color = "deeppink", marker = "x", zorder = 10)
            print(f"barrier = {barrier}kT, percolates at {timesteps[percAt[i]]}")
    for i, barrier in enumerate(sorted(set(unfoldBarrier))):
       
        ax.errorbar(timesteps, avUnfold[perc[i]] / numMol * 100, yerr = avUnfoldErr[perc[i]] / numMol * 100,
                    label = f"E_U = {barrier}kT", color = colours[i])
        ax.errorbar(timesteps, avUnfold[noperc[i]] / numMol * 100, yerr = avUnfoldErr[noperc[i]] / numMol * 100,
                    color = colours[i], linestyle = "dashdot")

    ax.set_xlabel("simulation frame")
    ax.set_ylabel(f"number of unfolded molecules (%)")
    ax.legend()
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/unfoldedMol_vf{vf}.png")
    ax.semilogx()
    ax.set_xlabel("log(timesteps)")
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/unfoldedMol_vf{vf}_semilog.png")
    ax.semilogy()
    ax.set_ylabel("log(number of unfolded molecules)")
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/unfoldedMol_vf{vf}_log.png")

    plt.close()

    print("plotted graphs of total unfolded molecules over time..")
    return 

def collapseUnfoldedOverTime(unfoldBarrier, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, suffixes):
    percAt = []
    perc = []
    noperc = []
    percTimestep = []
    finalUnfolded = []

    colours = ["maroon", "firebrick", "indianred", "lightcoral", "lightpink"]

    fig, ax = plt.subplots()
    with open(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/data/avPercolation.pkl", "rb") as f:
        avPercolation, avPercolationErr = pickle.load(f)
    with open(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/data/unfoldedOverTime.pkl", "rb") as f:
        avUnfold, avUnfoldErr = pickle.load(f)


    # ax.set_prop_cycle(color = ["maroon", "firebrick", "indianred", "lightcoral", "lightpink"])

    for i, (barrier, suffix) in enumerate(zip(unfoldBarrier, suffixes)):
        if bondsPerAtom == 2:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
        else:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
        with open(f"../runs/{conditions}/timesteps.pkl", "rb") as f:
            timesteps = pickle.load(f)
        key = (barrier, suffix)
        if suffix != "":
            percAt.append(None)
            noperc.append((key))
        if suffix == "":
            perc.append(key)
            percAt.append(np.where(avPercolation[barrier] == 3)[0][0])
            percTimestep.append(timesteps[percAt[i]])
            ax.scatter(timesteps[percAt[i]]  - percTimestep[i], avUnfold[key][percAt[i]], color = "teal", marker = "*")
            print(timesteps[percAt[i]], percTimestep[i])
            
            finalUnfolded.append((avUnfold[key][-1]))

            print(f"barrier {barrier}, num unfolded at perc / num unfolded post-relax : {avUnfold[key][percAt[i]] / avUnfold[key][-1]}")
            print(f"barrier {barrier}, num unfolded mol at perc = {avUnfold[key][percAt[i]]}")
    print(finalUnfolded)
    avFinalUnfold = (np.array(finalUnfolded)).mean()
    print(avFinalUnfold)
    ax.scatter(unfoldBarrier, finalUnfolded)
    # for i, barrier in enumerate(sorted(set(unfoldBarrier))):        
    #     ax.errorbar([t - timesteps[percAt[i]] for t in timesteps], [val + avFinalUnfold * np.log(barrier) for val in avUnfold[perc[i]]], yerr = avUnfoldErr[perc[i]],
    #                 label = f"barrier = {barrier}kT", color = colours[i])
    #     # ax.errorbar(timesteps - percTimestep[i], avUnfold[noperc[i]], yerr = avUnfoldErr[noperc[i]],
        #             color = colours[i], linestyle = "dashdot")
    ax.set_xlabel("simulation frame")
    ax.set_ylabel("log(number of unfolded molecules)")
    ax.set_title("cumulative total of unfolded mol")
    ax.legend()
    #ax.semilogx()
    #ax.semilogy()
    plt.show()

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
    
    #colours = (["maroon", "firebrick", "indianred", "lightcoral", "lightpink"])
    colours = (["black", "dimgrey", "grey", "darkgrey", "silver"])
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
            yerr=histogramsErr[key], color = colours[i],
            label=f"E_b = {barrier}kT {simtype}")
        ax.legend()

  
    plt.tight_layout()

    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/avBimodalAngle.png")
    return

def plotBimodalAndDensityUnfoldedCorr(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, suffixes):
    allHistograms = defaultdict(list)
    histogramsMean = {}
    histogramsErr = {}
    fig, axs = plt.subplots(2, 5, sharex = "row", sharey = "row", figsize = (20, 8))
    labels1 = [(0, 0), (0, 1), (0, 2), (0, 3), (0, 4)]
    labels2 = [(1, 0), (1, 1), (1, 2), (1, 3), (1, 4)]


    allFilteredCounts = []
    allFilteredAngles = []

    colours = (["black", "dimgrey", "grey", "darkgrey", "silver"])
    for i, (barrier, suffix) in enumerate(zip(unfoldBarriers, suffixes)):
        ax = axs[labels1[barrier - 1]]
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
            yerr=histogramsErr[key], color = colours[barrier - 1],
            label=f"$E_U$ = {barrier}kT")
    
        
        ax.set_xlim([20, 185])
        ax.set_ylim([0, np.max(histogramsMean[key] + 500)])
        ax.axvline(60, linestyle = "--", color = "mediumblue")
        ax.axvline(180, linestyle = "--", color = "red")
        ax.set_title(f"$E_U$ = {barrier}kT", fontsize = 15)

        # check the difference between mode and mean of the raw distribution
        mode_angle = bin_centers[np.argmax(histogramsMean[key])]
        mean_angle = np.average(bin_centers, weights = histogramsMean[key])
        print(mode_angle, mean_angle)
    for barrier in unfoldBarriers:
        if bondsPerAtom == 2:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}"
        else:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}"
        ax = axs[labels2[barrier - 1]]
        with open(f"../runs/{conditions}/timesteps.pkl", "rb") as f:
            timesteps = pickle.load(f)
        for runNum in range(numRuns):
            filename = f"Run{runNum}_{conditions}"
            with open(f"../runs/{conditions}/{filename}/analysis/densityUnfoldedCorr.pkl", "rb") as f:
                allCounts, allMeanAngles, meanCorr, allCountsIncEmpty, allMeanAnglesIncEmpty = pickle.load(f)
            filtered = [(v, n) for v, n in zip(allCountsIncEmpty, allMeanAnglesIncEmpty)
                        if not (v == 0 and n == 0)]
            filteredCounts, filteredAngles = zip(*filtered)
            allFilteredCounts.extend(filteredCounts)
            allFilteredAngles.extend(filteredAngles)
            
        h = ax.hist2d(allFilteredAngles, allFilteredCounts, bins = (40, 40), range = [[20, 185], [0, 200]], norm = mcolors.PowerNorm(0.4), cmap = plt.cm.plasma)
        #ax.set_title(f"E_U = {barrier}kT")

    fig.text(0.08, 0.65, "count of molecules", va='center', rotation='vertical', fontsize=15)
    fig.text(0.08, 0.28, "count of molecules in voxel", va='center', rotation='vertical', fontsize=15)

    #fig.text(0.5, 0.52, "molecule angle", ha='center', fontsize=12)
    fig.text(0.5, 0.04, f"average molecule angle in voxel (\u03B8)", ha='center', fontsize=15)

    #plt.tight_layout()
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/bimodalAndDensityUnfoldedCorr.png")
    return

def anglePopulation(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, suffixes):
    """plots a histogram of average distribution of folded vs unfolded molecules with
        energy barrier."""
    totalFolded = defaultdict(list)
    totalUnfolded = defaultdict(list)
    foldedMean = {}
    unfoldedMean = {}
    foldedError = {}
    unfoldedError = {}
    plotLabels = []
    fig, ax = plt.subplots()
    barWidth = 0.4 
    perc = []
    noperc = []
    x = np.arange(len(unfoldBarriers))
    for barrier, suffix in zip(unfoldBarriers, suffixes):
        if suffix != "":
            plotLabels.append(f"{barrier} no intermol.")
            groupID = 0
            noperc.append((barrier, suffix))
        if suffix == "":
            plotLabels.append(f"{barrier}")
            groupID = 1
            perc.append((barrier, suffix))
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
            angles = calcAngles(barrier, refoldBarrier, runNum, vf, numMol, bondsPerAtom, suffix)
            unfoldedMol = unfoldedMolecules(barrier, refoldBarrier, runNum, vf, numMol, angles, bondsPerAtom, suffix)
        
            numUnfolded = len(unfoldedMol[finalframe])
            numFolded = numMol - numUnfolded
            
            totalFolded[(barrier, suffix)].append(numFolded)
            totalUnfolded[(barrier, suffix)].append(numUnfolded)
            
        foldedArray = np.array(totalFolded[(barrier, suffix)])
        unfoldedArray = np.array(totalUnfolded[(barrier, suffix)])
            
        foldedMean[(barrier, suffix)] = foldedArray.mean()
        unfoldedMean[(barrier, suffix)] = unfoldedArray.mean()
        foldedError[(barrier, suffix)] = (foldedArray.std())
        unfoldedError[(barrier, suffix)] =(unfoldedArray.std())

        print(f"unfolded percent is {unfoldedMean[barrier, suffix] / numMol * 100}, error is {unfoldedError[(barrier, suffix)] / numMol * 100}")

    x = np.arange(len(sorted(set(unfoldBarriers))))

    offset = barWidth * 0.5

    ax.bar(x - offset, [foldedMean[key] / numMol * 100 for key in perc],
          width = barWidth, yerr =  [foldedError[key] / numMol * 100 for key in perc], color = "darkblue", label = "folded")
    ax.bar(x - offset, [unfoldedMean[key] / numMol * 100 for key in perc], width = barWidth, yerr = [unfoldedError[key] / numMol * 100 for key in perc], color = "red", bottom = [foldedMean[key] / numMol * 100 for key in perc], label = "unfolded")
    ax.bar(x + offset, [foldedMean[key] / numMol * 100 for key in noperc],
          width = barWidth, yerr =  [foldedError[key] / numMol * 100 for key in noperc], color = "dimgrey", label = "folded (no intermol.)")
    ax.bar(x + offset, [unfoldedMean[key] / numMol * 100 for key in noperc], width = barWidth, yerr = [unfoldedError[key] / numMol * 100 for key in noperc], color = "lightgrey", bottom = [foldedMean[key] / numMol * 100 for key in noperc], label = "unfolded (no intermol.)")    

    # try replacing a bit with ""
    ax.set_xticks(x)
    ax.set_xticklabels(sorted(set(unfoldBarriers)))
    ax.set_xlabel(f"$E_U$ (kT)", fontsize = 12)
    ax.set_ylabel("molecule count (%)", fontsize = 12)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
        fancybox=True, shadow=True, ncol=2)
    plt.tight_layout()
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/barrierAnglePop_vf{vf}.png")
    # with open(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/data/barrierAnglePop.txt") as f:
    #     f.write(f"# number of folded or unfolded molecules at end of simulation")
    #     f.write(f"E_U\tInteracting?\tfolded\tfolded err\tunfolded\tunfolded err\n")
    #     for key in perc:
    #         f.write(f"{key[0]}\tyes\t{foldedMean[key]}")
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
            angles = calcAngles(barrier, refoldBarrier, runNum, vf, numMol, bondsPerAtom, suffix)
            unfoldedMols = unfoldedMolecules(barrier, refoldBarrier, runNum, vf, numMol, angles, bondsPerAtom, suffix)
            
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


def plotFolded(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, suffixes):

    def make_func(barrier, refoldBarrier):
        a = np.exp(- barrier)
        b = np.exp(- barrier) + np.exp(- refoldBarrier)

        def func(t, gamma):
            t = np.asarray(t)
            return (1 / b) * np.exp(- gamma * b * t)
            #return (b - a)/b * (1 - np.exp(- gamma * b * t))
        return func

    percAt = []
    perc = []
    noperc = []

    colours = (["black", "dimgrey", "grey", "darkgrey", "silver"])

    fig, ax = plt.subplots()
    with open(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/data/avPercolation.pkl", "rb") as f:
        avPercolation, avPercolationErr = pickle.load(f)
    with open(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/data/unfoldedOverTime.pkl", "rb") as f:
        avUnfold, avUnfoldErr = pickle.load(f)

        
    for i, (barrier, suffix) in enumerate(zip(unfoldBarriers, suffixes)):
        if bondsPerAtom == 2:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
        else:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
        with open(f"../runs/{conditions}/timesteps.pkl", "rb") as f:
            timesteps = pickle.load(f) 
        key = (barrier, suffix)
        if suffix != "":
            percAt.append(None)
            noperc.append((key))
        if suffix == "":
            perc.append(key)
            percAt.append(np.where(avPercolation[barrier] == 3)[0][0])

            # ax.scatter(timesteps[percAt[i]], avUnfold[key][percAt[i]] / numMol * 100, color = "deeppink", marker = "x", zorder = 10)
            # print(f"barrier = {barrier}kT, percolates at {timesteps[percAt[i]]}")

    
    for i, barrier in enumerate(sorted(set(unfoldBarriers))):
        func = make_func(barrier, refoldBarrier)
        ax.errorbar(timesteps, 1 - (avUnfold[perc[i]] / numMol ), yerr = avUnfoldErr[perc[i]] / numMol ,
                    label = f"E_U = {barrier}kT", color = colours[i])
        # ax.errorbar(timesteps, 100 - (avUnfold[noperc[i]] / numMol * 100), yerr = avUnfoldErr[noperc[i]] / numMol * 100,
        #             color = colours[i], linestyle = "dashdot")

        print(avUnfold[key][5:], numMol)
        popt, pcov  = curve_fit(func, timesteps, 1 - avUnfold[perc[i]] / numMol )
        ax.plot(timesteps, func(timesteps, *popt), color = colours[i], linestyle = "--", label = "fit") 
    ax.semilogy()
    ax.semilogx()
    plt.show()

    return
