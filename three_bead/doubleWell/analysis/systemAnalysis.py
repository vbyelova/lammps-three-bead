import numpy as np
import pickle
import re
import os
from collections import defaultdict

from threeBeadClasses import Particle
from parseDump import *
from parseBonds import *

# let's get the system data first
unfoldBarriers = 3
refoldBarrier = 2
numRuns = 3
Vf = 0.07
numMol = 3938

conditions = f"langevin_10_unfold{unfoldBarriers}_refold{refoldBarrier}_Vf{Vf}_mol{numMol}"
for n in range(0, numRuns):
    filename = f"Run{n}_{conditions}"
    systemData = parseSystemData(f"../runs/{conditions}/{filename}/output/systemData.txt")

    boxLength, Npar, Nsteps, Nwrite, equilTime = systemData[0], systemData[1], systemData[2], systemData[3], systemData[4]

    # let's parse the dump file first and save it for later

    particles = [Particle() for n in range(Npar)]
    checking = readData(f"../runs/{conditions}/{filename}/output/dump.lammpstrj", particles, Nsteps, Nwrite, Npar, equilTime)
    dillParticles(f"../runs/{conditions}/{filename}/analysis/dillParticles.pkl", particles)

# count up total bonds in the system, both inter and intramolecular

#nBonds = totalBonds(f"../runs/{conditions}/{filename}/output/nbonds.dat")

# process the bonds to later check if they are percolated or not

#percolatedBonds = parseForPercolation(f"../runs/{conditions}/{filename}/analysis/dillParticles.pkl", f"../runs/{conditions}/{filename}/output/bondedatoms.dat", nBonds, boxLength)

#with open(f"../runs/{conditions}/{filename}/analysis/dillPercolation.pkl", "wb") as f:
#    pickle.dump(percolatedBonds, f)

# get all the bonded atoms and then write psf files for visualisation

#bondedAtoms = parseBondsForVis(f"../runs/{conditions}/{filename}/output/bondedatoms.dat", nBonds)
#writePSFfiles(nBonds, bondedAtoms, Npar, "dillParticles.pkl")

print("done :D")
