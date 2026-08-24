import numpy as np
import pickle
import re
import os
import subprocess
import matplotlib.pyplot as plt
from collections import defaultdict

from modules.parseDump import *
from modules.calcAngles import *
from modules.calcFracDim import * 
from modules.calcPercolation import *
from modules.calcStressTensor import *
from modules.threeBeadClasses import *
from modules.parseRDF import *
from modules.calcPorosity import *

realRadius = 30e-10
kT = 298 * 1.38e-23 # J / K
viscosity = 0.89e-3 # Pa * s, https://wiki.anton-paar.com/uk-en/water/ 

diffReal = kT / (6 * np.pi * viscosity * realRadius)

parRadius = 2 ** (1/6)
simRadius = parRadius + (0.5 * parRadius) / (np.sqrt(3) / 2)
scaleFactor = realRadius / simRadius

def MSD(particles, frames, timesteps):
    msd = []

    for frame in range(0, frames):
        info_t = particles[0].properties[frame]
        xt, yt, zt = info_t[1], info_t[2], info_t[3]
        info_0 = particles[0].properties[0]
        x0, y0, z0 = info_0[1], info_0[2], info_0[3]
        dx = (xt - x0)
        dy = (yt - y0)
        dz = (zt - z0)
        newMsd = dx**2 + dy**2 + dz**2
        msd.append(newMsd)

    return msd

barrier = 10
refoldBarrier = 2
vf = 1e-05
numMol = 1
numRuns = 20

dt = 0.001
allMsds = []
for runNum in range(0, numRuns):

    conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}"

    filename = f"Run{runNum}_{conditions}"
    systemData = parseSystemData(f"../runs/{conditions}/{filename}/output/systemData.txt")

    boxLength, nPar, nSteps, nWrite, equilTime = systemData[0], systemData[1], systemData[2], systemData[3], systemData[4]

    data = readData(f"../runs/{conditions}/{filename}/output/dump.lammpstrj", nSteps, nWrite, nPar, equilTime)
    particles, timesteps = data[0], data[1]

    frames = len(timesteps)
    newMsd = MSD(particles, frames, timesteps)
    allMsds.append(newMsd)

allMsds = np.array(allMsds)
avMsd = allMsds.mean(axis = 0)
avMsdErr = allMsds.std(axis = 0) / np.sqrt(numRuns)

time = np.array(timesteps) * dt

plt.scatter(time, avMsd)
plt.show()

start = int(0.2 * len(time))
end = int(0.9 * len(time))
slope, intercept = np.polyfit(time[start:end], avMsd[start:end], 1)
diffSim = slope / 6

timescale = diffSim * scaleFactor **2 / diffReal
realTime = time[-1] * timescale
print(realTime)
