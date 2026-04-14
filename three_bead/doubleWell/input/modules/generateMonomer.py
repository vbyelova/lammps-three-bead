# a python code to generate coordinates for monomers

import numpy as np
from numpy import random

bondLen = 1.112462048
mod = 0.9    # used to generate particles within a certain distance of boundary to avoid breaking

class Particle():
    """ An object that has x,y,z coordinates and an atom ID/type to be assigned."""
    def __init__(self):
        self.x = 0
        self.y = 0
        self.z = 0
        self.atomType = 0
        self.atomID = 0

def wrapping(dr, boxLength):
    """Check if separation of two particles is within a boxlength and wraps the coordinates if 
        such is the case."""
    if dr >= 0.5 * boxLength:
        dr -= boxLength
    elif dr <= - 0.5 * boxLength:
        dr += boxLength
    return dr

def equilateral(bondLen, p1, p2, p3):
    """generate particles and position them in a triangle."""

    p1.x = p2.x - bondLen * 0.5
    p1.y = p2.y - np.sqrt(3)/2 * bondLen
    p1.z = p2.z
    p3.x = p2.x + 0.5 * bondLen
    p3.y = p2.y - np.sqrt(3)/2 * bondLen
    p3.z = p2.z


    # type 1 = sticker / reactive
    # type 2 = hinge / non-reactive
    p1.atomType = 1
    p2.atomType = 2
    p3.atomType = 1

    return

def dimer(bondLen, p1, p2):
    p1.x = p2.x + bondLen
    p1.y = p2.y + bondLen
    p1.z = p2.z + bondLen

    return

def writeFile(filename, particles, numMol, boxLength):

    """ a function that allocates IDs to atoms and bonds then writes them to an output file."""

    atomID = 1
    bondCounter = 1
    molID = 1
    atomTypes = []
    molIDs = []

    #angles = list(range(1, numMol * 3 + 1))

    numPar = numMol * 3

    for num, p in enumerate(particles):
        p.atomID = atomID
        atomID += 1

    while bondCounter < numPar + 1:
        for n in range(0, 1):
            molIDs.append(molID)
        bondCounter += 1
        atomTypes.append(1)      
        # create a list of which atoms are bonded to which
        molID += 1

    
    with open(filename, "w") as f:

        f.write("LAMMPS config file for N monomers\n")

        f.write(f"\n{numPar} atoms\n")
        f.write(f"0 bonds\n")
        f.write(f"0 angles\n")
        f.write("1 atom types\n")
        f.write("2 bond types\n")

        f.write("\n#Sim box size\n")
        f.write(f"-{0.5 * boxLength} {0.5 * boxLength} xlo xhi\n")
        f.write(f"-{0.5 * boxLength} {0.5 * boxLength} ylo yhi\n")
        f.write(f"-{0.5 * boxLength} {0.5 * boxLength} zlo zhi\n")

        f.write("\nMasses\n")
        f.write("#AtomID mass\n")
        f.write("1  1\n")

        f.write("\nAtoms\n")
        f.write("#atomID moleculeID atomType X Y Z\n")
        for num, p in enumerate(particles):
            f.write(f"{p.atomID} {molIDs[num]} {atomTypes[num]} {p.x} {p.y} {p.z} \n")


def generateMonomer(filename, numMol, boxLength):
    particles = [Particle() for _ in range(numMol * 3)]
    for p in particles:
        p.x = random.uniform(- mod * 0.5 * boxLength, mod * 0.5 * boxLength)
        p.y = random.uniform(- mod * 0.5 * boxLength, mod * 0.5 * boxLength)
        p.z = random.uniform(- mod * 0.5 * boxLength, mod * 0.5 * boxLength)
    print("initialised particles...")

    for num,p in enumerate(particles):
        wrapping(p.x, boxLength)
        wrapping(p.y, boxLength)
        wrapping(p.z, boxLength)
    print("wrapping around periodic box edges...")
    print("writing to file...")
    writeFile(filename, particles, numMol, boxLength)
    print("generated molecule input file!")
