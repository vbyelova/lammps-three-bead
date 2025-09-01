# a python code to generate coordinates for N three bead molecules.
import numpy as np
from numpy import random


numMol = 1
bondLen = 2
boxLength = 5
mod = 0.9        # used to generate particles within a certain distance of boundary to avoid breaking

class Particle():
    """ An object that has x,y,z coordinates and an atom ID/type to be assigned."""
    def __init__(self):
        self.x = random.uniform(- mod * boxLength, mod * boxLength)
        self.y = random.uniform(- mod * boxLength, mod * boxLength)
        self.z = random.uniform(- mod * boxLength, mod * boxLength)
        self.atomType = 0
        self.atomID = 0

def equilateral(bondLen, p1, p2, p3):
    """generate particles and position them in a triangle."""

    p1.x = p2.x - bondLen # * 0.5
    p1.y = p2.y # + np.sqrt(3)/2 * bondLen
    p1.z = p2.z
    p3.x = p2.x + 0.5 * bondLen
    p3.y = p2.y + np.sqrt(3)/2 * bondLen
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
    numPar = numMol * 3

    for num, p in enumerate(particles):
        p.atomID = atomID
        bondIDs.append(bondID)
        atomID += 1
        bondID += 1
    
    while bondCounter < numPar:
        for n in range(0, 3):
            molIDs.append(molID)
        # configure whether the bond can break or not
        # bond type 1 = sticker-hinge
        # bond type 2 = sticker=sticker, used for unfolding
        bondTypes.append(1)
        bondTypes.append(1)
        bondTypes.append(2)
        # atom type 1 = sticker
        # atom type 2 = hinge
        atomTypes.append(1)
        atomTypes.append(2)
        atomTypes.append(1)      
        # create a list of which atoms are bonded to which
        bondedAtom1.append(bondCounter)
        bondedAtom2.append(bondCounter + 1)
        bondedAtom1.append(bondCounter + 1)
        bondedAtom2.append(bondCounter + 2)
        bondedAtom1.append(bondCounter)
        bondedAtom2.append(bondCounter + 2)
        molID += 1
        bondCounter += 3

    
    with open(fileName, "w") as f:

        f.write("LAMMPS config file for N 3 bead models\n")

        f.write(f"\n{numPar} atoms\n")
        f.write(f"{numPar} bonds\n")
        f.write("2 atom types\n")
#        f.write("3 atom types\n")
        f.write("5 bond types\n")

        f.write("\n#Sim box size\n")
        f.write(f"-{boxLength} {boxLength} xlo xhi\n")
        f.write(f"-{boxLength} {boxLength} ylo yhi\n")
        f.write(f"-{boxLength} {boxLength} zlo zhi\n")

        f.write("\nMasses\n")
        f.write("#AtomID mass\n")
        f.write("1  1\n")
        f.write("2  1\n")
#        f.write("3  1\n")

        f.write("\nAtoms\n")
        f.write("#atomID moleculeID atomType X Y Z\n")
        for num, p in enumerate(particles):
            f.write(f"{p.atomID} {molIDs[num]} {atomTypes[num]} {p.x} {p.y} {p.z} \n")

        f.write("\nBonds\n")
        f.write("#bondID bondType particlesBonding\n")
        for num, p in enumerate(particles):
            f.write(f"{bondIDs[num]} {bondTypes[num]} {bondedAtom1[num]} {bondedAtom2[num]} \n")

particles = [Particle() for _ in range(numMol * 3)]

counter = 0
while counter < numMol:
    equilateral(bondLen, particles[counter], particles[counter + 1], particles[counter + 2])
    counter += 3

writeFile("network.in", particles, numMol)
