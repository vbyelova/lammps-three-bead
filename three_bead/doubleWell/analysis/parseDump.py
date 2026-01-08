# a python code to parse data output from lammps simulations
# Victoria Byelova

import numpy as np
import matplotlib.pyplot as plt
import pickle
import re
import configparser

from threeBeadClasses import Particle

def parseSystemData(name):
    print("checking directory: ", name)
    config = configparser.ConfigParser()
    config.read(name)
    boxLength = int(config.get("systemData", "boxLength"))
    Npar = int(config.get("systemData", "Npar"))
    Nsteps = eval(config.get("systemData", "Nsteps"))
    Nwrite = int(config.get("systemData", "Nwrite"))
    equilTime = int(config.get("systemData", "equilTime"))
    return boxLength, Npar, Nsteps, Nwrite, equilTime

def readData(name, particles, numStep, everyN, numPar, equilTime):
    """ A function to read in data from a lammps dump file and parse it into arrays."""
    # initialise counters
    p = 0
    frame = 0
    counter = 0
    cols = 6
    timestep = int((numStep - equilTime)/everyN)

    pattern = r"""\d+\s\d+\s\d+\n""" # digits sandwiched by spaces

    # generate array of zeros to store data in each particle
    while p < numPar:
        particles[p].properties = np.zeros((timestep + 1, cols))
        p += 1
    
    # read line by line to see if data or other stuff
    with open(name) as f:
        for line in f.readlines():
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
    return particles

def dillParticles(name, particles):
    with open(name, "wb") as f:
        pickle.dump(particles, f)
    print(f"saved dump info as {name}.")

def checkClumping(unfoldBarriers, refoldBarrier, numRuns, Vf, numMol):
    conditions = f"langevin_10_unfold{unfoldBarriers[0]}_refold{refoldBarrier}_Vf{Vf}_mol{numMol}"
    for n in range(0, numRuns):
        filename = f"Run{n}_{conditions}"
        with open(f"../runs/{conditions}/{filename}/analysis/dillParticles.pkl", "rb") as f:
            particles = pickle.load(f)
        xcoords = []
        for p in particles:
            xcoords.append(p.properties[-1, 1])

        hist, bins = np.histogram(xcoords, bins = 80)
        plt.hist(xcoords, bins, color = "palevioletred")
        plt.title(f"nve + overdamped langevin (10)")
        plt.xlabel("x coordinate of particles")
        plt.ylabel("number of particles")
        plt.legend([f"unfoldBarrier = {unfoldBarriers[0]} kT\nnumMol = {numMol}\nnumber of bins = {len(bins)-1}"])
        plt.show()
