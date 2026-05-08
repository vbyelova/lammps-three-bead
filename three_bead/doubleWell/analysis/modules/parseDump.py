# python code to parse data output from lammps simulations
# Victoria Byelova

import numpy as np
import math
import matplotlib.pyplot as plt
import pickle
import re
import configparser

from collections import defaultdict
from .threeBeadClasses import *

def parseSystemData(name):
    """parses the system information printed out previously to be used in analysis"""
    print("checking directory: ", name)
    config = configparser.ConfigParser()
    config.read(name)
    boxLength = int(config.get("systemData", "boxLength"))
    Npar = int(config.get("systemData", "Npar"))
    Nsteps = int(config.get("systemData", "Nsteps"))
    Nwrite = int(config.get("systemData", "Nwrite"))
    equilTime = int(config.get("systemData", "equilTime"))
    return boxLength, Npar, Nsteps, Nwrite, equilTime

def readData(name, numStep, everyN, numPar, equilTime):
    """ A function to read in data from a lammps dump file and parse it into arrays."""
    # initialise counters
    frame = 0
    counter = 0
    cols = 6
    timesteps = []

    pattern = r"""\d+\s+\d+\s+\d+""" # digits sandwiched by spaces

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

def parseBondInfo(barrier, refoldBarrier, runNum, Vf, numMol, nBonds, boxLength):
    """gets stress tensor information for each bond as fx, fy, fz, dx, dy and dz."""
    conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{Vf}_mol{numMol}"
    filename = f"Run{runNum}_{conditions}"

    with open(f"../runs/{conditions}/timesteps.pkl", "rb") as f:
        timesteps = pickle.load(f)

    frame = 0
    bondCounter = 0
    bondInfo = defaultdict(list)

    for x, num in enumerate(nBonds):
        bondInfo[x] = [Bond() for n in range(nBonds[x])]
    
    print("saving bond force and direction.. ")
    with open(f"../runs/{conditions}/{filename}/output/bondinfo.dat", "r") as f:
        for line in f:
            if re.search(r"\d+\s+-?[\d.]+\s+-?[\d.]+\s+-?[\d.]+\s+-?[\d.]+\s+-?[\d.]+\s+-?[\d.]+", line):
                #print("matched")
                lineData = line.rsplit()
                num = int(lineData[0]) - 1 # change to 0 indexing
                for x in range(0, 9):
                    bondInfo[frame][num].properties[0 + x] = lineData[x]
                
                bondCounter += 1
                if bondCounter == nBonds[frame]:
                    bondCounter = 0
                    frame += 1
                    if frame == len(nBonds):
                        break
    
    print("len bond info ", len(bondInfo))
    print("len nbonds ",len(nBonds))
    print("num frames ", frame)

    return bondInfo



def checkClumping(unfoldBarriers, refoldBarrier, numRuns, Vf, numMol):
    """produces a histogram of the x coordinates of a simulation's final frame."""
    conditions = f"unfold{unfoldBarriers[0]}_refold{refoldBarrier}_Vf{Vf}_mol{numMol}"
    for n in range(0, numRuns):
        filename = f"Run{n}_{conditions}"
        with open(f"../runs/{conditions}/{filename}/analysis/dillParticles.pkl", "rb") as f:
            particles = pickle.load(f)
        xcoords = []
        for p in particles:
            xcoords.append(p.properties[-1, 1])

        hist, bins = np.histogram(xcoords, bins = 150)
        plt.hist(xcoords, bins, color = "palevioletred")
        plt.xlabel("x coordinate of particles")
        plt.ylabel("number of particles")
        plt.legend([f"unfoldBarrier = {unfoldBarriers[0]} kT\nnumMol = {numMol}\nnumber of bins = {len(bins)-1}"])
        plt.show()

def totalBonds(filename):
    """counts total number of bonds in each frame of the simulation and saves to a list."""
    nBonds = []
    bondedAtoms = []
    with open(filename, "r") as f:
        next(f)
        for line in f:
            line = int(line[:-1])
            nBonds.append(line)

    print("counted total bonds across the system in each frame..")
    return nBonds

def boundaryCheck(val, boxLength):
    """checks if a coordinate is crossing a boundary condition and returns values for periodic boundaries."""
    if val < - 0.5 * boxLength:
        return -1
    if val > 0.5 * boxLength:
        return 1
    else:
        return 0

def wrapping(val, boxLength):
    while val > 0.5 * boxLength or val < -0.5 * boxLength:
        if val > 0.5 * boxLength:
            val -= boxLength
        elif val < 0.5 * boxLength:
            val += boxLength
    return


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

def frameByFramePerc(barrier, refoldBarrier, runNum, vf, numMol, boxLength):
    """finds bonded particles and saves them and over what boundary they are bonded. formatted
        for one value per line, for easier parsing in c++."""
    bondCounter = 0
    frame = 0
    conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}"
    filename = f"Run{runNum}_{conditions}"
    with open(f"../runs/{conditions}/{filename}/analysis/nBonds.pkl", "rb") as bondsfile:
        nBonds = pickle.load(bondsfile)
    with open(f"../runs/{conditions}/timesteps.pkl", "rb") as f:
        timesteps = pickle.load(f)
    with open(f"../runs/{conditions}/{filename}/analysis/particleTraj.pkl", "rb") as parfile:
        particles = pickle.load(parfile)
    print(len(particles), " is len particles")
    print(numMol * 3, " is numPar")
    with open(f"../runs/{conditions}/{filename}/output/bondinfo.dat", "r") as f:
        with open(f"../runs/{conditions}/{filename}/analysis/percinfo.txt", "w") as w:
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


def checkBondLength(barrier, refoldBarrier, runNum, Vf, numMol, nBonds, boxLength):
    """checks if a bond length is unrealistically long, if so then this means that wrapping is not working"""
    bondCounter = 0
    frame = 0
    conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{Vf}_mol{numMol}"
    filename = f"threeBead_Run{runNum}_{conditions}"

    with open(f"../runs/{conditions}/{filename}/output/bondinfo.dat", "r") as f:
        for line in f:
            if re.search(r"^[1-9]*\s+\d+", line):
                bondID = line.rsplit()[0]
                bondLen = float(line.rsplit()[7])
                if bondLen > 0.5 * boxLength:
                    print(f"frame {frame}:bond length is too large at {bondLen} for bond {bondID}")
                bondCounter += 1

                if bondCounter == nBonds[frame]:
                    bondCounter = 0
                    frame += 1
                    print(f"moving to frame  {frame}")
        
    return print("bond length check completed.")

def coordination(barrier, refoldBarrier, runNum, Vf, numMol, nBonds, bondInfo):
    """a function that finds the coordination of the atoms belonging to a three-bead molecule.
        First finds whether the bond is intramolecular and intermolecular and only counts the
        intermolecular bonds"""
    conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{Vf}_mol{numMol}"
    filename = f"Run{runNum}_{conditions}"
    with open(f"../runs/{conditions}/timesteps.pkl", "rb") as f:
        timesteps = pickle.load(f)
    frame = 0
    bondCounter = 0
    atomCoordination = np.zeros((len(timesteps), numMol * 3))
    moleculeCoordination = np.zeros((len(timesteps), numMol))
    molCoordinationCount = defaultdict(list)
    trackedBonds = set()
    while frame < len(nBonds):
        bondedAtom1 = int(bondInfo[frame][bondCounter].properties[1]) - 1
        bondedAtom2 = int(bondInfo[frame][bondCounter].properties[2]) - 1

        mol1 = bondedAtom1 // 3
        mol2 = bondedAtom2 // 3
        if mol1 != mol2:
            # intermolecular bond case
            bondID = tuple(sorted([bondedAtom1, bondedAtom2]))
            
            if bondID not in trackedBonds:
                trackedBonds.add(bondID)
                atomCoordination[frame][bondedAtom1] += 1
                atomCoordination[frame][bondedAtom2] += 1
                moleculeCoordination[frame][mol1] += 1
                moleculeCoordination[frame][mol2] += 1
                
        
        bondCounter += 1
        if bondCounter == nBonds[frame]:
            trackedBonds.clear()
            bondCounter = 0 
            frame += 1
            print("frame ", frame)

    print("counted intermolecular bonds per atom..")
    return moleculeCoordination

def plotAvCoordination(unfoldBarriers, refoldBarrier, numRuns, numMol, Vf, boxLength):
    """a function to find the average molecule coordination per run"""
    allCoordination = defaultdict(list)
    avCoordination = {}
    avCoordinationErr = {}

    fig, ax = plt.subplots()
    for barrier in unfoldBarriers:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{Vf}_mol{numMol}"
        with open(f"../runs/{conditions}/timesteps.pkl", "rb") as f:
            timesteps = pickle.load(f)
        for runNum in range(numRuns):
            filename = f"Run{runNum}_{conditions}"
            with open(f"../runs/{conditions}/{filename}/analysis/percDims.pkl", "rb") as f:
                percDims = pickle.load(f)
            
                if 3 not in percDims:
                    print(f"skipping {filename}")
                    continue
            with open(f"../runs/{conditions}/{filename}/analysis/nBonds.pkl", "rb") as f:
                nBonds = pickle.load(f)
            with open(f"../runs/{conditions}/{filename}/analysis/bondInfo.pkl", "rb") as f:
                bondInfo = pickle.load(f)            
            molCoordination = coordination(barrier, refoldBarrier, runNum, Vf, numMol,
                                        nBonds, bondInfo)
            avCoordinationPerFrame = molCoordination.mean(axis = 1)
            allCoordination[barrier].append(avCoordinationPerFrame)
        coordinationArray = np.array(allCoordination[barrier])
        avCoordination[barrier] = coordinationArray.mean(axis = 0)
        avCoordinationErr[barrier] = coordinationArray.std(axis = 0)


        ax.errorbar(timesteps, avCoordination[barrier], yerr = avCoordinationErr[barrier],
                    label = f"barrier = {barrier}kT")
    ax.set_xlabel("simulation frame")
    ax.set_ylabel("molecule coordination")
    ax.legend()
    plt.savefig(f"../runs/boxLength{boxLength}/vf{Vf}/avCoordination_vf{Vf}.png")
    ax.semilogx()
    ax.set_xlabel("simulation frame (semilog)")
    plt.savefig(f"../runs/boxLength{boxLength}/vf{Vf}/avCoordination_vf{Vf}_semilog.png")

    plt.close()
    print("plotted average molecule coordination..")
    return avCoordination, avCoordinationErr

def plotTogether(unfoldBarriers, refoldBarrier, numRuns, numMol, vf, boxLength):
    with open(f"../runs/boxLength{boxLength}/vf{vf}/data/avNewUnfoldedMol.pkl", "rb") as f:
        avNewUnfoldedMol, avNewUnfoldedMolErr = pickle.load(f)
    with open(f"../runs/boxLength{boxLength}/vf{vf}/data/avPercolation.pkl", "rb") as f:
        avPercolation, avPercolationErr = pickle.load(f)
    # with open(f"../runs/boxLength{boxLength}/vf{vf}/data/avPressure.pkl", "rb") as f:
    #     avPressure, avPressureErr = pickle.load(f)
    with open(f"../runs/boxLength{boxLength}/vf{vf}/data/unfoldedOverTime.pkl", "rb") as f:
        avUnfoldOverTime, avUnfoldOverTimeErr = pickle.load(f)
    
    fig, (ax1, ax2, ax3) = plt.subplots(3, sharex = True, figsize = (15, 15))
    fig.subplots_adjust(wspace = 0, hspace = 0)
    for barrier in unfoldBarriers:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_numMol{numMol}"
        with open(f"../runs/{conditions}/timesteps.pkl", "rb") as f:
            timesteps = pickle.load(f)
        ax1.errorbar(timesteps, avUnfoldOverTime[barrier], yerr = avUnfoldOverTimeErr[barrier],
                     label = f"barrier = {barrier}kT")
        ax2.errorbar(timesteps, avNewUnfoldedMol[barrier], yerr = avNewUnfoldedMolErr[barrier],
                label = f"unfolding barrier = {barrier}kT")
        ax3.errorbar(timesteps, avPercolation[barrier], yerr = avPercolationErr[barrier],
                    label = f"barrier = {barrier}kT")
        # ax4.errorbar(timesteps, avPressure[barrier], yerr = avPressureErr[barrier],
        #             label = f"barrier = {barrier}kT")
        
    ax1.set_ylabel("number of unfolded molecules")
    ax2.set_ylabel("new unfolded molecules")
    ax3.set_ylabel("percolation dimension")
    ax1.legend()
    ax2.legend()
    ax3.legend()
#    ax4.set_ylabel("average pressure")
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/sharedaxis_vf{vf}.png")
    ax1.semilogx()
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/sharedaxis_vf{vf}_semilog.png")
    plt.close()

def plotTogetherShifted(unfoldBarriers, refoldBarrier, numRuns, numMol, vf, boxLength):
    with open(f"../runs/boxLength{boxLength}/vf{vf}/data/avNewUnfoldedMol.pkl", "rb") as f:
        avNewUnfoldedMol, avNewUnfoldedMolErr = pickle.load(f)
    # with open(f"../runs/boxLength{boxLength}/vf{vf}/data/avPercolation.pkl", "rb") as f:
    #     avPercolation, avPercolationErr = pickle.load(f)
    # with open(f"../runs/boxLength{boxLength}/vf{vf}/data/avPressure.pkl", "rb") as f:
    #     avPressure, avPressureErr = pickle.load(f)
    with open(f"../runs/boxLength{boxLength}/vf{vf}/data/unfoldedOverTime.pkl", "rb") as f:
        avUnfoldOverTime, avUnfoldOverTimeErr = pickle.load(f)
    fig, (ax1, ax2) = plt.subplots(2, sharex = True, figsize = (15, 15))
    fig.subplots_adjust(wspace = 0, hspace = 0)
    for barrier in unfoldBarriers:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}"
        with open(f"../runs/{conditions}/timesteps.pkl", "rb") as f:
            timesteps = pickle.load(f)
            timesteps = np.array(timesteps)
        percIndex = []
        for runNum in range(numRuns):
            filename = f"Run{runNum}_{conditions}"
            with open(f"../runs/{conditions}/{filename}/analysis/percDims.pkl", "rb") as f:
                percDims = pickle.load(f)
                if 3 in percDims:
                    percIndex.append(percDims.index(3))
                else:
                    print(f"no percolation in {conditions} run {runNum}")
                    continue
        avPercIndex = math.ceil(np.mean(percIndex))
        percTimestep = timesteps[avPercIndex]
        ax1.errorbar(timesteps - percTimestep, avUnfoldOverTime[barrier], yerr = avUnfoldOverTimeErr[barrier],
                     label = f"barrier = {barrier}kT")
        ax2.errorbar(timesteps - percTimestep, avNewUnfoldedMol[barrier], yerr = avNewUnfoldedMolErr[barrier],
                label = f"unfolding barrier = {barrier}kT")
        # ax3.errorbar(timesteps - percTimestep, avPercolation[barrier], yerr = avPercolationErr[barrier],
        #             label = f"barrier = {barrier}kT")
        # ax4.errorbar(timesteps, avPressure[barrier], yerr = avPressureErr[barrier],
        #             label = f"barrier = {barrier}kT")
        
    ax1.set_ylabel("number of unfolded molecules")
    ax2.set_ylabel("new unfolded molecules")
#    ax3.set_ylabel("percolation dimension")
    ax1.legend()
    ax2.legend()
#    ax3.legend()
#    ax4.set_ylabel("average pressure")
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/shifted_sharedaxis_vf{vf}.png")
    ax1.semilogx()
    ax1.set_xlabel("timesteps shifted")
    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/shifted_sharedaxis_vf{vf}_semilog.png")
    plt.close()
