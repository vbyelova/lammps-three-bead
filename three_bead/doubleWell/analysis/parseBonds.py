import pickle
import re
from collections import defaultdict

from parseDump import * 

def totalBonds(filename):
    totalBonds = []
    bondedAtoms = []
    with open(filename, "r") as f:
        next(f)
        for line in f.readlines():
            line = int(line[:-1])
            totalBonds.append(line)

    print("counted total bonds across the system in each frame..")
    return totalBonds

def boundaryCheck(val, boxLength):
    if val < - 0.5 * boxLength or val > 0.5 * boxLength:
        return 1
    else:
        return 0

def parseForPercolation(picklefile, filename, totalBonds, boxLength):
    with open(picklefile, "rb") as f:
        particles = pickle.load(f)
    percolatedBonds = defaultdict(list)
    totalFrames = len(totalBonds)
    bondCounter = 0
    frame = 0
    lineID = 0
    with open(filename, "r") as f:
        lines = f.readlines()
        while frame < totalFrames and lineID < len(lines):
            while bondCounter < totalBonds[frame]:
                if re.search(r"\d+\s+\d+\s+\d+", lines[lineID]):
                    atom1 = int(lines[lineID].rsplit()[1]) - 1
                    atom2 = int(lines[lineID].rsplit()[2]) - 1
                    dx = particles[atom1].properties[frame, 1] - particles[atom2].properties[frame, 1]
                    dy = particles[atom1].properties[frame, 2] - particles[atom2].properties[frame, 2]
                    dz = particles[atom1].properties[frame, 3] - particles[atom2].properties[frame, 3]
                    overXbound = boundaryCheck(dx, boxLength)
                    overYbound = boundaryCheck(dy, boxLength)
                    overZbound = boundaryCheck(dz, boxLength)
                    percolatedBonds[frame].append([atom1, atom2, [overXbound, overYbound, overZbound]])
                    lineID += 1
                    bondCounter += 1
                else:
                    lineID += 1
            bondCounter = 0
            frame += 1
            
    print("saved intermolecular bond information..")
    return percolatedBonds


def parseBondsForVis(filename, totalBonds):
    frame = 0
    totalFrames = len(totalBonds)

    lineID = 0
    bondCounter = 0
    bondedAtoms = []

    with open(filename, "r") as f:
        lines = f.readlines()
        while frame < totalFrames and lineID < len(lines):
            while bondCounter < totalBonds[frame]:
                #print(f"frame: {frame} bondCounter: {bondCounter} lineID: {lineID}\n")
                if re.search(r"\d+\s+\d+\s+\d+", lines[lineID]):
                    #print(f"frame {frame} lineID {lineID} bound count {bondCounter}")
                    bondedAtoms.append(int(lines[lineID].rsplit()[1]))
                    bondedAtoms.append(int(lines[lineID].rsplit()[2]))
                    lineID += 1
                    bondCounter += 1
                else:
                    lineID += 1
            bondCounter = 0
            frame += 1

    print("got a list of bonded particles for processing..")
    return bondedAtoms

def writePSFfiles(totalBonds, bondedAtoms, Npar, dillParticles):
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

def bondStressTensor(filename, totalBonds):
    frame = 0
    totalFrames = len(totalBonds)
    
    lineID = 0
    bondCounter = 0

    with open(filename, "r") as f:
        lines = f.readlines()
        while frame < totalFrames and lineID < len(lines):
            while bondCounter < totalBonds[frame]:
                if re.search(r"\d+\s+\d+\s+\d+\d+\s+\d+\s+\d+", lines[lineID]):
                    lineID += 1
                    bondCounter += 1
                else:
                    lineID += 1
            bondCounter = 0
            frame += 1
            

def checkBondLength(unfoldBarriers, refoldBarrier, numRuns, Vf, numMol):
    bondCounter = 0
    frame = 0
    conditions = f"unfold{unfoldBarriers[0]}_refold{refoldBarrier}_Vf{Vf}_mol{numMol}"
    filename = f"threeBead_Run0_{conditions}"
    nBonds = totalBonds(f"../runs/{conditions}/{filename}/output/nbonds.dat")

    systemData = parseSystemData(f"../runs/{conditions}/{filename}/output/systemData.txt")
    boxLength = systemData[0]
    with open(f"../runs/{conditions}/{filename}/output/bonddistance.dat", "r") as f:
        for line in f.readlines():
            if re.search(r"^[1-9]*\s+\d+", line):
                bondID = line.rsplit()[0]
                bondLen = float(line.rsplit()[1])
                if bondLen > 0.5 * boxLength:
                    print(f"frame {frame}:bond length is too large at {bondLen} for bond {bondID}")
                bondCounter += 1

                if bondCounter == nBonds[frame]:
                    bondCounter = 0
                    frame += 1
                    print(f"moving to frame  {frame}")
        
    return print("bond length check completed.")
                        