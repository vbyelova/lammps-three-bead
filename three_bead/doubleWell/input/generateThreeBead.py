# a python code to generate coordinates for N three bead molecules.
# modified so that only sticker-hinge bonds exist for double-well potential.

import numpy as np
from numpy import random

# modify these two for different volume fraction
numMol = 10
boxLength = 20

bondLen = 1.112462048
mod = 0.8      # used to generate particles within a certain distance of boundary to avoid breaking

class Particle():
    """ An object that has x,y,z coordinates and an atom ID/type to be assigned."""
    def __init__(self):
        self.x = random.uniform(- mod * 0.5 * boxLength, mod * 0.5 * boxLength)
        self.y = random.uniform(- mod * 0.5 * boxLength, mod * 0.5 * boxLength)
        self.z = random.uniform(- mod * 0.5 * boxLength, mod *  0.5 * boxLength)
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


def writeFile(fileName, particles, numMol):

    """ a function that allocates IDs to atoms and bonds then writes them to an output file."""

    atomID = 1
    bondID = 1
    bondCounter = 1
    molID = 1
    atomTypes = []
    molIDs = []
    bondIDs = []
    bondTypes = []
    bondedAtom1 = []
    bondedAtom2 = []
    angles = list(range(1, numMol * 3 + 1))

    numPar = numMol * 3

    for num, p in enumerate(particles):
        p.atomID = atomID
        bondIDs.append(bondID)
        atomID += 1
        bondID += 1
    
    while bondCounter < numPar:
        for n in range(0, 3):
            molIDs.append(molID)

        bondTypes.append(1)
        bondTypes.append(1)

        atomTypes.append(1)
        atomTypes.append(1)
        atomTypes.append(1)      
        # create a list of which atoms are bonded to which
        bondedAtom1.append(bondCounter)
        bondedAtom2.append(bondCounter + 1)
        bondedAtom1.append(bondCounter + 1)
        bondedAtom2.append(bondCounter + 2)
        molID += 1
        bondCounter += 3

    
    with open(fileName, "w") as f:

        f.write("LAMMPS config file for N 3 bead models\n")

        f.write(f"\n{numPar} atoms\n")
        f.write(f"{int(numPar * 2/3)} bonds\n")
        f.write(f"{numMol} angles\n")
        f.write("1 atom types\n")
        f.write("2 bond types\n")
        f.write("1 angle types\n")

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

        f.write("\nBonds\n")
        f.write("#bondID bondType particlesBonding\n")
        for num in range(0, int(len(particles)*(2/3))):
            f.write(f"{bondIDs[num]} {bondTypes[num]} {bondedAtom1[num]} {bondedAtom2[num]} \n")

        f.write("\nAngles\n")
        f.write("#angleID angleType particle1 particle2 particle3\n")
        num = 0
        counter = 1
        while num < numMol * 3:
            f.write(f"{counter} 1 {angles[num]} {angles[num + 1]} {angles[num + 2]}\n")
            counter += 1
            num += 3

particles = [Particle() for _ in range(numMol * 3)]

counter = 0
while counter < numMol:
    equilateral(bondLen, particles[counter], particles[counter + 1], particles[counter + 2])
    counter += 3

for num,p in enumerate(particles):
    wrapping(p.x, boxLength)
    wrapping(p.y, boxLength)
    wrapping(p.z, boxLength)

writeFile("network.in", particles, numMol)
