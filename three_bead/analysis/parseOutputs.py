# a python code to parse data output from lammps simulations


import numpy as np
import pickle
import re

from systemData import Nsteps, Nwrite, Npar


class Particle():
    """ A  particle with id, type, molecule id and number of intramolecular
        sticky bonds."""
    def __init__(self):
        self.properties = []

def readData(name, particles, numStep, everyN, numPar):
    """ A function to read in data from a lammps dump file and parse it into arrays."""
    # initialise counters
    p = 0
    row = 0
    counter = 0
    cols = 7
    timestep = int(numStep/everyN)

    # generate array of zeros to store data in each particle
    while p < numPar:
        particles[p].properties = np.zeros((timestep + 1, cols))
        p += 1
    
    # read line by line to see if data or other stuff
    with open(name) as text:
        for line in text.readlines():
            if re.search(pattern, line):
                line = line[:-1]
                counter += 1
                num = int(line.rsplit()[0]) - 1

    # split up line and save into correct slot in array
                for x in range(0, cols):
                    particles[num].properties[row, 0 + x] = line.rsplit()[x]

    # loop once all particles have been accounted for
                if counter == numPar:
                    row += 1
                    counter = 0
                  
    return particles

def dillParticles(particles):
    with open("dillParticles.pkl", "wb") as f:
        pickle.dump(particles,f)

pattern = r"""\d+\s\d+\s\d+\s\d+\n""" # digits sandwiched by spaces

particles = [Particle() for n in range(Npar)]
readData("../output/dump.lammpstrj", particles, Nsteps, Nwrite, Npar)
dillParticles(particles)



# only use those where atom 1 and atom 3 both have 0 ssintra bonds
# iterate every 3?
 
# find separation of these molecules

# load in av. free energy of molecule
# can i do it per each molecule? 
