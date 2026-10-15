# python code to parse data output from lammps simulations
# Victoria Byelova

import numpy as np
import math
import matplotlib.pyplot as plt
import pickle
import re
import configparser

from itertools import cycle
from collections import defaultdict
from .threeBeadClasses import *

def parseSystemData(name, simType):
    """parses the system information printed out previously to be used in analysis"""
    print("checking directory: ", name)
    config = configparser.ConfigParser()
    config.read(name)
    boxLength = int(config.get("systemData", "boxLength"))
    Npar = int(config.get("systemData", "Npar"))
    Nsteps = int(config.get("systemData", "Nsteps"))
    Nwrite = int(config.get("systemData", "Nwrite"))
    if simType == "step strain":
        equilTime = 0
    else:
        equilTime = int(config.get("systemData", "equilTime"))
    return boxLength, Npar, Nsteps, Nwrite, equilTime

def readData(name, numStep, everyN, numPar, equilTime = 0):
    """ A function to read in data from a lammps dump file and parse it into arrays."""
    # initialise counters
    frame = 0
    counter = 0
    cols = 6
    timesteps = []

    pattern = r"\d+\s+-?[\d.]+\s+-?[\d.]+\s+-?[\d.]+\s+-?[\d.]+\s+-?[\d.]+\s+-?[\d.]+" # digits sandwiched by spaces

    # read line by line to see if data or other stuff
    with open(name, 'r') as f:
        for line in f:
            if "ITEM: TIMESTEP" in line:
                timestep = (next(f)).split()[0]
                timesteps.append(int(timestep))

    particles = [Particle() for n in range(numPar)]
    for p in range(numPar):
        particles[p].properties = np.zeros((len(timesteps), cols))


    with open(name) as f:
        for line in f:
            if re.search(pattern, line):
                counter += 1
                parts = line.split()
                num = int(parts[0]) - 1 # change to 0 indexing for python

    # split up line and save into correct slot in array
                for x in range(0, cols):
                    particles[num].properties[frame, 0 + x] = parts[x]

    # loop once all particles have been accounted for
                if counter == numPar:
                    frame += 1
                    counter = 0
                  
    print("dump file processed..")
    return particles, timesteps

def parseBondInfo(barrier, refoldBarrier, runNum, vf, numMol, boxLength, bondsPerAtom, prob, suffix, simType, strainFraction):
    """gets stress tensor information for each bond as fx, fy, fz, dx, dy and dz."""
    if simType == "step strain":
        prefix = f"stepStrain{strainFraction}_"
    elif simType == "gelation":
        prefix = ""
    if bondsPerAtom == 2:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
    else:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
        if prob < 1:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}_prob{prob}{suffix}"
    filename = f"Run{runNum}_{conditions}"

    with open(f"../runs/{conditions}/{prefix}timesteps.pkl", "rb") as f:
        timesteps = pickle.load(f)
    with open(f"../runs/{conditions}/{filename}/{prefix}analysis/nBonds.pkl", "rb") as f:
        nBonds = pickle.load(f)

    bondInfo = {}
    dumpSteps = []
    with open(f"../runs/{conditions}/{filename}/{prefix}output/bondinfo.dat") as f:
        frame = -1
        for line in f:
            if line.startswith("ITEM: TIMESTEP"):
                dumpSteps.append(int(f.readline()))
            elif line.startswith("ITEM: NUMBER OF ENTRIES"):
                n = int(f.readline())
            elif line.startswith("ITEM: ENTRIES"):
                frame += 1
                rows = [f.readline().split() for _ in range(n)]
                bondInfo[frame] = np.array(rows, dtype=float)
                assert n == nBonds[frame], (frame, n, nBonds[frame])
    

    num = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
    pat = re.compile(rf"^\s*\d+(?:\s+{num}){{9}}\s*$")

    matched, missed = 0, []
    with open(f"../runs/{conditions}/{filename}/{prefix}output/bondinfo.dat") as f:
        for line in f:
            if line.startswith("ITEM") or len(line.split()) != 10:
                continue
            if pat.fullmatch(line.rstrip("\n")):
                matched += 1
            else:
                missed.append(line)
    print(matched, sum(nBonds), len(missed))
    print(missed[:5])
    return bondInfo


def totalBonds(filename):
    """counts total number of bonds in each frame of the simulation and saves to a list."""
    nBonds = []
    bondedAtoms = []
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            step, value = line.split()[:2]
            nBonds.append(int(float(value)))

    print("counted total bonds across the system in each frame..")
    return nBonds

def parseStress(barrier, refoldBarrier, runNum, vf, numMol, boxLength, bondsPerAtom, prob, suffix, simType, strainFraction):
    """extracts global stress calculated during lammps for use in calculating shear modulus"""

def boundaryCheck(val, boxLength):
    """checks if a coordinate is crossing a boundary condition and returns values for periodic boundaries."""
    if val < - 0.5 * boxLength:
        return -1
    if val > 0.5 * boxLength:
        return 1
    else:
        return 0

def wrapping(val, boxLength):
    """wraps a value around the periodic boundaries."""
    while val > 0.5 * boxLength or val < -0.5 * boxLength:
        if val > 0.5 * boxLength:
            val -= boxLength
        elif val < - 0.5 * boxLength:
            val += boxLength
    return val

def parseForPercolation(particles, filename, nBonds, boxLength, timesteps):
    """finds intermolecular bonds and tracks them as well as if they are across a boundary."""
    percolatedBonds = defaultdict(list)
    bondedAtoms = defaultdict(list)
    bondCounter = 0
    frame = 0
    with open(filename, "r") as f:
        for line in f:
            if re.search(r"\d+\s+\d+\s+\d+", line):
                atom1 = int(line.rsplit()[1]) - 1 # change to 0 indexing
                atom2 = int(line.rsplit()[2]) - 1
                dx = particles[atom1].properties[frame, 1] - particles[atom2].properties[frame, 1]
                dy = particles[atom1].properties[frame, 2] - particles[atom2].properties[frame, 2]
                dz = particles[atom1].properties[frame, 3] - particles[atom2].properties[frame, 3]
                overXbound = boundaryCheck(dx, boxLength)
                overYbound = boundaryCheck(dy, boxLength)
                overZbound = boundaryCheck(dz, boxLength)
                percolatedBonds[frame].append([atom1, atom2, [overXbound, overYbound, overZbound]])
                bondedAtoms[frame].append([atom1, atom2])
                bondCounter += 1
                if bondCounter == nBonds[frame]:
                    bondCounter = 0
                    frame += 1
                if frame >= len(nBonds):
                    break

    print("saved intermolecular bond information..")
    return percolatedBonds, bondedAtoms

def frameByFramePerc(barrier, refoldBarrier, runNum, vf, numMol, boxLength, bondsPerAtom, prob, suffix, simType, strainFraction):
    """finds bonded particles and saves them and over what boundary they are bonded. formatted
        for one value per line, for easier parsing in c++."""
    if simType == "step strain":
        prefix = f"stepStrain{strainFraction}_"
    elif simType == "gelation":
        prefix = ""
    bondCounter = 0
    frame = 0
    if bondsPerAtom == 2:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
    else:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
        if prob < 1:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}_prob{prob}{suffix}"
    filename = f"Run{runNum}_{conditions}"
    with open(f"../runs/{conditions}/{filename}/{prefix}analysis/nBonds.pkl", "rb") as bondsfile:
        nBonds = pickle.load(bondsfile)
    with open(f"../runs/{conditions}/{prefix}timesteps.pkl", "rb") as f:
        timesteps = pickle.load(f)
    with open(f"../runs/{conditions}/{filename}/{prefix}analysis/particleTraj.pkl", "rb") as parfile:
        particles = pickle.load(parfile)
    print(len(particles), " is len particles")
    print(numMol * 3, " is numPar")
    with open(f"../runs/{conditions}/{filename}/{prefix}output/bondinfo.dat", "r") as f:
        with open(f"../runs/{conditions}/{filename}/{prefix}analysis/percinfo.txt", "w") as w:
            w.write(f"{timesteps[frame]}\n")
            w.write(f"{nBonds[frame]}\n")
            for line in f:
                if re.search(r"\d+\s+\d+\s+\d+", line):
                    
                    info = line.rsplit()
                    atom1 = int(info[1]) - 1 # change to 0 indexing
                    atom2 = int(info[2]) - 1
                    #print(atom1, atom2)
                    dx = particles[atom1].properties[frame, 1] - particles[atom2].properties[frame, 1]
                    dy = particles[atom1].properties[frame, 2] - particles[atom2].properties[frame, 2]
                    dz = particles[atom1].properties[frame, 3] - particles[atom2].properties[frame, 3]
                    overXbound = boundaryCheck(dx, boxLength)
                    overYbound = boundaryCheck(dy, boxLength)
                    overZbound = boundaryCheck(dz, boxLength)                
                    w.write(f"{atom1}\n"
                            f"{atom2}\n"
                            f"{overXbound}\n"
                            f"{overYbound}\n"
                            f"{overZbound}\n")
                    bondCounter += 1
                    if bondCounter == nBonds[frame]:
                        print(frame, timesteps[frame], len(timesteps), len(nBonds))
                        bondCounter = 0
                        frame += 1
                        if frame >= len(nBonds):
                            break
                        w.write(f"{timesteps[frame]}\n")
                        #print(frame)
                        w.write(f"{nBonds[frame]}\n")

    print("found bonded particles and the boundaries they cross over..")
    return

def parseBondsForVis(filename, nBonds, timesteps):
    """generates a list of bonded atoms."""
    frame = 0
    totalFrames = len(timesteps)


    bondCounter = 0
    bondedAtoms = []

    with open(filename, "r") as f:
        for line in f:
            #print(f"frame: {frame} bondCounter: {bondCounter} lineID: {lineID}\n")
            if re.search(r"\d+\s+\d+\s+\d+", line):
                #print(f"frame {frame} lineID {lineID} bound count {bondCounter}")
                bondedAtoms.append(int(line.rsplit()[1]))
                bondedAtoms.append(int(line.rsplit()[2]))
                bondCounter += 1
                if bondCounter == nBonds[frame]:
                    bondCounter = 0
                    frame += 1

    print("got a list of bonded particles for processing..")
    return bondedAtoms



def coordination(barrier, refoldBarrier, runNum, vf, numMol, boxLength, bondsPerAtom, prob, suffix, simType, strainFraction):
    """a function that finds the coordination of the atoms belonging to a three-bead molecule.
        First finds whether the bond is intramolecular and intermolecular and only counts the
        intermolecular bonds"""
    if simType == "step strain":
        prefix = f"stepStrain{strainFraction}_"
    elif simType == "gelation":
        prefix = ""
    if bondsPerAtom == 2:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
    else:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
        if prob < 1:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}_prob{prob}{suffix}"
    filename = f"Run{runNum}_{conditions}"
    with open(f"../runs/{conditions}/{prefix}timesteps.pkl", "rb") as f:
        timesteps = pickle.load(f)
    with open(f"../runs/{conditions}/{filename}/analysis/nBonds.pkl", "rb") as f:
        nBonds = pickle.load(f)
    with open(f"../runs/{conditions}/{filename}/analysis/bondInfo.pkl", "rb") as f:
        bondInfo = pickle.load(f)
    frame = 0
    bondCounter = 0
    intermol = 0
    atomCoordination = np.zeros((len(timesteps), numMol * 3))
    moleculeCoordination = np.zeros((len(timesteps), numMol))
    trackedBonds = set()
    
    for frame in range(len(timesteps)):
        b = bondInfo[frame]
        a1 = b[:, 1].astype(int) - 1
        a2 = b[:, 2].astype(int) - 1
        inter = (a1 // 3) != (a2 // 3)
        pairs = np.unique(np.sort(np.stack([a1[inter], a2[inter]], axis=1), axis=1), axis=0)
        atomCoordination[frame] = np.bincount(pairs.ravel(), minlength=numMol * 3)
        moleculeCoordination[frame] = np.bincount(pairs.ravel() // 3, minlength=numMol)
        if len(b) - len(pairs) != 2 * numMol:
            print(f"frame {frame}: intra = {len(b) - len(pairs)}, expected {2*numMol}")


    centralAtomCoordination = atomCoordination[:, 1::3]
    print("counted intermolecular bonds per atom..")
    return moleculeCoordination, atomCoordination, centralAtomCoordination

def plotCentralCoordination(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, prob, suffixes, simType, strainFraction):
    """plots the average coordination number of the central bead of the molecule over time for every energetic barrier."""
    if simType == "step strain":
        prefix = f"stepStrain{strainFraction}_"
    elif simType == "gelation":
        prefix = ""
    fig, ax = plt.subplots()
    colours = (["black", "dimgrey", "grey", "darkgrey", "silver"])
    allCentralCoordination = defaultdict(list)
    percAt = []
    for i, (barrier, suffix) in enumerate(zip(unfoldBarriers, suffixes)):
        if suffix != "":
            continue
        if bondsPerAtom == 2:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
        else:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
            if prob < 1:
                conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}_prob{prob}{suffix}"
        with open(f"../runs/{conditions}/{prefix}timesteps.pkl", "rb") as f:
            timesteps = pickle.load(f)
        with open(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/{prefix}data/avPercolation.pkl", "rb") as f:
                avPercolation, avPercolationErr = pickle.load(f)

        for runNum in range(numRuns):
            filename = f"Run{runNum}_{conditions}"
            print(filename)
            with open(f"../runs/{conditions}/{filename}/analysis/percDims.pkl", "rb") as f:
                percDims = pickle.load(f)
            # with open(f"../runs/{conditions}/{filename}/analysis/nBonds.pkl", "rb") as f:
            #     nBonds = pickle.load(f)
            # with open(f"../runs/{conditions}/{filename}/analysis/bondInfo.pkl", "rb") as f:
            #     bondInfo = pickle.load(f)            
            # with open(f"../runs/{conditions}/{filename}/analysis/coordination.pkl", "rb") as f:
            #     coords = pickle.load(f)
            coords = coordination(barrier, refoldBarrier, runNum, vf, numMol, boxLength, bondsPerAtom, prob, suffix, simType, strainFraction)
            centralAtomCoordination = coords[2]
        
            allCentralCoordination[barrier].append(centralAtomCoordination.mean(axis = 1))
        coordArray = np.array(allCentralCoordination[barrier])
        avCentralCoord = coordArray.mean(axis = 0)
        avCentralCoordErr = coordArray.std(axis = 0)
        percAt.append(np.where(avPercolation[barrier] == 3)[0][0])
        ax.errorbar(timesteps, avCentralCoord, yerr = avCentralCoordErr, color = colours[i], label = f"$E_U$ = {barrier}kT")
        ax.scatter(timesteps[percAt[i]], avCentralCoord[percAt[i]], color = "deeppink", marker = "x", zorder = 10)
        print(f"barrier = {barrier}kT, coordination = {avCentralCoord[percAt[i]]} +- {avCentralCoordErr[percAt[i]]}, timestep of perc = {timesteps[percAt[i]]}")
    
    ax.set_xlabel("timesteps", fontsize = 12)
    ax.set_ylabel("central particle coordination", fontsize = 12)
    ax.legend()
    #plt.show()
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/centralCoord.jpg")
    ax.semilogx()
    ax.semilogy()
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/centralCoord_semilog.jpg")

    return
            # find percolation point
            # find coordination of central atoms at that point
            

def plotAvCoordination(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, prob, suffixes, simType, strainFraction):
    """a function to find the average molecule coordination per run"""
    if simType == "step strain":
        prefix = f"stepStrain{strainFraction}_"
    elif simType == "gelation":
        prefix = ""
    allCoordination = defaultdict(list)
    avCoordination = {}
    avCoordinationErr = {}
    percAt = []
    colours = (["black", "dimgrey", "grey", "darkgrey", "silver"])

    fig, ax = plt.subplots()
    for i, (barrier, suffix) in enumerate(zip(unfoldBarriers, suffixes)):
        if suffix != "":
            continue
        if bondsPerAtom == 2:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
        else:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
            if prob < 1:
                conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}_prob{prob}{suffix}"
        with open(f"../runs/{conditions}/{prefix}timesteps.pkl", "rb") as f:
            timesteps = pickle.load(f)
        with open(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/{prefix}data/avPercolation.pkl", "rb") as f:
            avPercolation, avPercolationErr = pickle.load(f)
        for runNum in range(numRuns):
            filename = f"Run{runNum}_{conditions}"
            print(filename)
            with open(f"../runs/{conditions}/{filename}/analysis/percDims.pkl", "rb") as f:
                percDims = pickle.load(f)
            # with open(f"../runs/{conditions}/{filename}/analysis/nBonds.pkl", "rb") as f:
            #     nBonds = pickle.load(f)
            # with open(f"../runs/{conditions}/{filename}/analysis/bondInfo.pkl", "rb") as f:
            #     bondInfo = pickle.load(f)            
            molCoordination = (coordination(barrier, refoldBarrier, runNum, vf, numMol, boxLength,
                                            bondsPerAtom, prob, suffix, simType, strainFraction))[0]
            avCoordinationPerFrame = molCoordination.mean(axis = 1)
            allCoordination[barrier].append(avCoordinationPerFrame)
        coordinationArray = np.array(allCoordination[barrier])
        avCoordination[barrier] = coordinationArray.mean(axis = 0)
        avCoordinationErr[barrier] = coordinationArray.std(axis = 0)
        print(f"minimum coordination is {min(avCoordination[barrier])}, maximum is {max(avCoordination[barrier])}, end is {avCoordination[barrier][-1]}")
        percAt.append(np.where(avPercolation[barrier] == 3)[0][0])
        ax.scatter(timesteps[percAt[i]], avCoordination[barrier][percAt[i]], color = "deeppink", marker = "x", zorder = 10)

        ax.errorbar(timesteps, avCoordination[barrier], yerr = avCoordinationErr[barrier],
                    color = colours[i], label = f"$E_U$ = {barrier}kT")
    ax.set_xlabel("simulation frame")
    ax.set_ylabel("molecule coordination")
    ax.legend()
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/avCoordination_vf{vf}.jpg")
    ax.semilogx()
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/avCoordination_vf{vf}_semilog.jpg")

    print("plotted average molecule coordination..")
    return avCoordination, avCoordinationErr

def plotTotalInterMol(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, prob, suffixes, simType, strainFraction):
    """plots scatter graph of average intermolecular bond count over simulation for each
        energy barrier."""
    if simType == "step strain":
        prefix = f"stepStrain{strainFraction}_"
    elif simType == "gelation":
        prefix = ""
    allBonds = defaultdict(list)
    avNBonds = {}
    avNBondsErr = {}
    percAt = []

    fig, ax = plt.subplots()

    colours = cycle(["maroon", "firebrick", "indianred", "lightcoral", "lightpink"])

    for i, (barrier, suffix) in enumerate(zip(unfoldBarriers, suffixes)):
        if bondsPerAtom == 2:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
        else:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
            if prob < 1:
                conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}_prob{prob}{suffix}"
        with open(f"../runs/{conditions}/{prefix}timesteps.pkl", "rb") as f:
            timesteps = pickle.load(f)
        with open(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/{prefix}data/avPercolation.pkl", "rb") as f:
            avPercolation, avPercolationErr = pickle.load(f)
        for runNum in range(numRuns):
            filename = f"Run{runNum}_{conditions}"
            with open(f"../runs/{conditions}/{filename}/{prefix}analysis/nBonds.pkl", "rb") as f:
                nBonds = pickle.load(f)
            
            allBonds[barrier].append(nBonds)
        
        allBondsArray = np.array(allBonds[barrier])
        avNBonds[barrier] = allBondsArray.mean(axis = 0)
        avNBondsErr[barrier] = allBondsArray.std(axis = 0)
        ax.errorbar(timesteps, avNBonds[barrier], yerr = avNBondsErr[barrier],
                    color = next(colours), label = f"E_U = {barrier}kT")
        percAt.append(np.where(avPercolation[barrier] == 3)[0][0])
        print(avNBonds[barrier][percAt[i]])
    
        ax.scatter(timesteps[percAt[i]], avNBonds[barrier][percAt[i]],
                   marker = "*", color = "black")
    plt.legend()
    plt.semilogx()
    plt.xlabel("timesteps")
    plt.ylabel("number of intermol. bonds")
    plt.show()
    return
