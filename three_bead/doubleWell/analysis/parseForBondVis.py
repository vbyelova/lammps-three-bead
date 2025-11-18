import pickle
import re
from collections import defaultdict
from systemData import Npar, boxLength

def totalBonds():
    totalBonds = []
    bondedAtoms = []
    with open("../output/nbonds.dat", "r") as f:
        next(f)
        for line in f.readlines():
            line = int(line[:-1])
            totalBonds.append(line)

    return totalBonds

def boundaryCheck(val):
    if val < - 0.5 * boxLength or val > 0.5 * boxLength:
        return 1
    else:
        return 0

def parseForPercolation(particles, totalBonds):
    percolatedBonds = defaultdict(list)
    totalFrames = len(totalBonds)
    bondCounter = 0
    frame = 0
    lineID = 0
    with open("../output/bondedatoms.dat", "r") as f:
        lines = f.readlines()
        while frame < totalFrames and lineID < len(lines):
            while bondCounter < totalBonds[frame]:
                if re.search(r"\d+\s+\d+\s+\d+", lines[lineID]):
                    atom1 = int(lines[lineID].rsplit()[1]) - 1
                    atom2 = int(lines[lineID].rsplit()[2]) - 1
                    dx = particles[atom1].properties[frame, 1] - particles[atom2].properties[frame, 1]
                    dy = particles[atom1].properties[frame, 2] - particles[atom2].properties[frame, 2]
                    dz = particles[atom1].properties[frame, 3] - particles[atom2].properties[frame, 3]
                    overXbound = boundaryCheck(dx)
                    overYbound = boundaryCheck(dy)
                    overZbound = boundaryCheck(dz)
                    percolatedBonds[frame].append([atom1, atom2, [overXbound, overYbound, overZbound]])
                    lineID += 1
                    bondCounter += 1
                else:
                    lineID += 1
            bondCounter = 0
            frame += 1
            
    return percolatedBonds


def parseBondsForVis(totalBonds):
    frame = 0
    totalFrames = len(totalBonds)

    lineID = 0
    bondCounter = 0
    bondedAtoms = []

    with open("../output/bondedatoms.dat", "r") as f:
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


    while frame < totalBonds[frame]:
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
            
with open("dillParticles.pkl", "rb") as f:
    particles = pickle.load(f)

totalBonds = totalBonds()
percolatedBonds = parseForPercolation(particles, totalBonds)
bondedAtoms = parseBondsForVis(totalBonds)
writePSFfiles(totalBonds, bondedAtoms, Npar, "dillParticles.pkl")

with open("dillPercolation.pkl", "wb") as f:
    pickle.dump(percolatedBonds, f)
