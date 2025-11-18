import numpy as np
import pickle
import re

from collections import defaultdict
from systemData import boxLength
from threeBeadClasses import Particle

def readTotalBonds(name):
    """A function to read in data about how many total bonds there are across
        the system.
    """

    totalBonds = []

    with open(name) as text:
        next(text)
        for line in text.readlines():
            line = line[:-1]
            sumBonds = int(line.rsplit()[1]) + int(line.rsplit()[2]) + int(line.rsplit()[3])
            # get total number of bonds in system per timestep
            totalBonds.append(sumBonds)
            #print(line)
    return totalBonds

def readBondedAtomData(name, particles, totalBonds):
    """ A function to read in data from a lammps dump file and parse it into arrays.
        This function processes data that identifies which two atoms are in a bond.
        Data is read in as: bond index atom1 atom2
        The function then saves intermolecular bonds in a dictionary, with each key
        corresponding to a timestep.
     """
    # initialise counter
    row = 0
    counter = 0
    
    bonds = defaultdict(list)
    bondForces = defaultdict(list)
    #molecules = [Molecule() for _ in range(int(Npar / 3) + 1)]
    pattern = r"""\d+\s\d+\s\d+""" # digits sandwiched by spaces
    # generate array of zeros to store data in each particle

    # read line by line to see if data or other stuff
    with open(name) as text:
        for line in text.readlines():
            if re.search(pattern, line):
                #print("pattern matched")
                line = line[:-1]
                # subtract 1 for indexing purposes 
                atom1 = int(line.rsplit()[1]) - 1
                atom2 = int(line.rsplit()[2]) - 1
                mol1 = int(particles[atom1].properties[row, 5] - 1)
                mol2 = int(particles[atom2].properties[row, 5] - 1)
                forces = [float(line.rsplit()[3]), float(line.rsplit()[4]), float(line.rsplit()[5])]
                bondForces[row].append(forces)
                if mol1 != mol2:
                    dx = particles[atom1].properties[row, 1] - particles[atom2].properties[row, 1]
                    dy = particles[atom1].properties[row, 2] - particles[atom2].properties[row, 2]
                    dz = particles[atom1].properties[row, 3] - particles[atom2].properties[row, 3]
                    if dx < - 0.5 * boxLength or dx > 0.5 * boxLength:
                        overXboundary = 1
                    else:
                        overXboundary = 0
                    if dy < - 0.5 * boxLength or dy > 0.5 * boxLength:
                        overYboundary = 1
                    else:
                        overYboundary = 0
                    if dz < - 0.5 * boxLength or dz > 0.5 * boxLength:
                        overZboundary = 1
                    else:
                        overZboundary = 0
                    
                    overBoundaries = [overXboundary, overYboundary, overZboundary]

                    bonds[row].append([atom1, atom2, mol1, mol2, overBoundaries])
                elif mol1 == mol2:
                    pass
                counter += 1
                if counter == totalBonds[row]:
                    counter = 0
                    row += 1
                    if row >= len(totalBonds):
                        break
                
    return bonds, bondForces


totalBonds = readTotalBonds("../output/nbonds.dat")
with open("dillParticles.pkl", "rb") as f:
    particles = pickle.load(f)
bondData = readBondedAtomData("../output/bondedatoms.dat", particles, totalBonds)
bonds = bondData[0]
bondForces = bondData[1]
#print(bondForces)
with open("intermolBondsPerTimestep.pkl", "wb") as f:
    pickle.dump(bonds, f)
