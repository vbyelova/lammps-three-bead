# python code to parse data output from lammps simulations
# Victoria Byelova

import numpy as np
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
    Nsteps = eval(config.get("systemData", "Nsteps"))
    Nwrite = int(config.get("systemData", "Nwrite"))
    equilTime = int(config.get("systemData", "equilTime"))
    return boxLength, Npar, Nsteps, Nwrite, equilTime

def readData(name, numStep, everyN, numPar, equilTime):
    """ A function to read in data from a lammps dump file and parse it into arrays."""
    # initialise counters
    p = 0
    frame = 0
    counter = 0
    cols = 6
    timestep = int((numStep - equilTime)/everyN)

    pattern = r"""\d+\s\d+\s\d+\n""" # digits sandwiched by spaces

    # generate array of zeros to store data in each particle
    particles = [Particle() for n in range(numPar)]
    while p < numPar:
        particles[p].properties = np.zeros((timestep + 1, cols))
        p += 1
    
    # read line by line to see if data or other stuff
    with open(name) as f:
        for line in f:
            if re.search(pattern, line):
                line = line[:-1]
                counter += 1
                num = int(line.rsplit()[0]) - 1

    # split up line and save into correct slot in array
                for x in range(0, cols):
                    particles[num].properties[frame, 0 + x] = line.rsplit()[x]

    # loop once all particles have been accounted for
                if counter == numPar:
                    frame += 1
                    counter = 0
                  
    print("dump file processed..")
    print([particles[0].properties[0]])
    return particles

def parseBondInfo(barrier, refoldBarrier, runNum, Vf, numMol, nBonds):
    """gets stress tensor information for each bond as fx, fy, fz, dx, dy and dz."""
    conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{Vf}_mol{numMol}"
    filename = f"Run{runNum}_{conditions}"
    systemData = parseSystemData(f"../runs/{conditions}/{filename}/output/systemData.txt")
    boxLength, Npar, Nsteps, Nwrite, equilTime = systemData[0], systemData[1], systemData[2], systemData[3], systemData[4]

    vol = boxLength ** 3

    frame = 0
    bondCounter = 0
    bondInfo = defaultdict(list)

    print("saving bond force and direction.. ")
    with open(f"../runs/{conditions}/{filename}/output/bondinfo.dat", "r") as f:
        for line in f:
            if bondCounter == 0:
                bondInfo[frame] = [Bond() for n in range(nBonds[frame])]

            if re.search(r"\d+\s+-?[\d.]+\s+-?[\d.]+\s+-?[\d.]+\s+-?[\d.]+\s+-?[\d.]+\s+-?[\d.]+", line):
                #print("matched")
                lineData = line.rsplit()
                num = int(lineData[0]) - 1
                for x in range(0, 9):
                    bondInfo[frame][num].properties[0 + x] = lineData[x]
                
                bondCounter += 1
                if bondCounter == nBonds[frame]:
                    bondCounter = 0
                    frame += 1
    
    print("len bond info ", len(bondInfo))
    print("len nbonds ",len(nBonds))
    print("num frames ", frame)
    print(bondInfo[1][0].properties)

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
    """checks if a coordinate is crossing a boundary condition and wraps around the box if it is."""
    if val < - 0.5 * boxLength or val > 0.5 * boxLength:
        return 1
    else:
        return 0

def parseForPercolation(particles, filename, nBonds, boxLength):
    """finds intermolecular bonds and tracks them as well as if they are across a boundary."""
    percolatedBonds = defaultdict(list)
    bondedAtoms = defaultdict(list)
    totalFrames = len(nBonds) 
    bondCounter = 0
    frame = 0
    with open(filename, "r") as f:
        for line in f:
            if re.search(r"\d+\s+\d+\s+\d+", line):
                atom1 = int(line.rsplit()[1]) - 1
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

    print("saved intermolecular bond information..")
    return percolatedBonds, bondedAtoms

def FrameByFramePerc(particles, filename, writeto, nBonds, boxLength):
    """finds bonded particles and saves them and over what boundary they are bonded. formatted
        for use in a c++ file."""
    totalFrames = len(nBonds)
    bondCounter = 0
    frame = 0
    with open(filename, "r") as f:
        with open(writeto, "w") as w:
            w.write(f"frame 0\n")
            for line in f:
                if re.search(r"\d+\s+\d+\s+\d+", line):
                    info = line.rsplit()
                    atom1 = int(info[0]) - 1
                    atom2 = int(info[1]) - 1
                    print(atom1, atom2)
                    dx = particles[atom1].properties[frame, 1] - particles[atom2].properties[frame, 1]
                    dy = particles[atom1].properties[frame, 2] - particles[atom2].properties[frame, 2]
                    dz = particles[atom1].properties[frame, 3] - particles[atom2].properties[frame, 3]
                    overXbound = boundaryCheck(dx, boxLength)
                    overYbound = boundaryCheck(dy, boxLength)
                    overZbound = boundaryCheck(dz, boxLength)                
                    w.write(f"{atom1} {atom2} {overXbound} {overYbound} {overZbound}\n")
                    bondCounter += 1
                    if bondCounter == nBonds[frame]:
                        bondCounter = 0
                        frame += 1
                        w.write(f"frame {frame}\n")
        
    return

def parseBondsForVis(filename, nBonds):
    """generates a list of bonded atoms."""
    frame = 0
    totalFrames = len(nBonds)


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

def writePSFfiles(totalBonds, bondedAtoms, Npar, dillParticles):
    """"generates psf files to draw shorter bonds in vmd. i wish this worked."""
    with open(dillParticles, "rb") as f:
        particles = pickle.load(f)
    frame = 0
    lineCount = 0
    bCount = 0
    nums = [_ for _ in range(0, Npar + 1)]
    tCount = 1
    seg = "SYS"
    resid = 1
    resname = "MOL"
    name = "A"
    typeB = "B"
    charge = 0
    mass = 1

    while frame < len(totalBonds):
        with open(f"../output/psfFiles/datafile%05d.psf" % frame, "w") as f:
            f.write(f"  1 !NTITLE\n")
            f.write(f"PSF\n\n")
            f.write(f"0 !NATOM\n")
            # f.write(f"{Npar} !NATOM\n")
            # for num, p in enumerate(particles):
            #     f.write(f"{int(particles[num].properties[frame, 0]):>8d} {seg:<4s} "
            #             f"{resid:>4d} {resname:<4s} {name:4s} {typeB:<4s} "
            #             f"{charge:>10.6f} {mass:13.4f}       0\n")
            
            f.write(f"\n{totalBonds[frame]} !NBONDS: bonds")
            while bCount < totalBonds[frame]:

                if bCount % 4 == 0:
                    f.write(f"\n")
                f.write(f"{bondedAtoms[lineCount]:>8d} {bondedAtoms[lineCount + 1]:>8d} ")
                lineCount += 2
                bCount += 1

            f.write(f"\n\n{Npar} !NTHETA: angles\n")
            while tCount < Npar:
                f.write(f"{nums[tCount]} {nums[tCount + 1]} {nums[tCount + 2]}\n")
                tCount += 3
            tCount = 1
            f.write(f"\n       0 !NPHI: dihedrals\n\n")
            f.write(f"    0 !NPHI: impropers\n\n")
            f.write(f"    0 !NCRTERM: cross-terms\n\n")
        bCount = 0
        frame += 1
    
    print("psf files made for visualising bonds..")
    return

            

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

def coordination(numMol, nBonds, bondInfo):
    """a function that finds the coordination of the atoms belonging to a three-bead molecule.
        First finds whether the bond is intramolecular and intermolecular and only counts the
        intermolecular bonds"""
    frame = 0
    bondCounter = 0
    atomCoordination = np.zeros((int(len(nBonds)) - 1, numMol * 3))
    moleculeCoordination = np.zeros((int(len(nBonds)) - 1, numMol))
    while frame < len(nBonds) - 1:
        bondedAtom1 = int(bondInfo[frame][bondCounter].properties[1]) - 1
        bondedAtom2 = int(bondInfo[frame][bondCounter].properties[2]) - 1

        mol1 = bondedAtom1 // 3
        mol2 = bondedAtom2 // 3
        if mol1 == mol2:
            # intramolecular bond case
            pass
        else:
            # intermolecular bond case
            atomCoordination[frame][bondedAtom1] += 1
            atomCoordination[frame][bondedAtom2] += 1
            moleculeCoordination[frame][mol1] += 1
            moleculeCoordination[frame][mol2] += 1
        bondCounter += 1
        if bondCounter == nBonds[frame]:
            bondCounter = 0
            frame += 1

    print("counted intermolecular bonds per atom..")
    avAtomCoord = np.mean(atomCoordination, axis = 1)
    xvals = np.arange(len(nBonds) - 1)
    avMolCoord = np.mean(moleculeCoordination, axis = 1)
    plt.scatter(xvals, avMolCoord)
    plt.semilogx()
    plt.xlabel("time frame (semilog axis)")
    plt.ylabel("average molecule coordination")
    plt.show()
    return
