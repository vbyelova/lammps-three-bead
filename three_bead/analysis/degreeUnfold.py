import matplotlib.pyplot as plt
import numpy as np
import pickle

from operator import add
from parseOutputs import Particle
from systemData import Nsteps, Nwrite

class Molecule():
    def __init__(self):
        self.molNum = 0
        self.separation = []

def checkUnfolded(particles):
    """Checks if a particle underwent unfolding during the simulation and adds to a list.
        Syntax for columns is id x y z type mol intramol_bonds."""
    unfoldedPar = []
    for num, p in enumerate(particles):
        if particles[num].properties[0][6] != particles[num].properties[-1][6]:
            unfoldedPar.append(particles[num])
            print(particles[num].properties)
 
    return unfoldedPar

def getSep(particles):
    """Gets separation between sticker particles in three bead molecules at each time step."""
    p = 0
    m = 0
    n = 0
    row = 0
    loopNum = 0

    # there are 2 particles per molecule
    molecules = [Molecule for _ in range(0, int(len(particles) * 0.5))]

    # every 2 rows in the data is 1 timestep
    if len(particles) == 0:
        raise ValueError("No particles in the system have undergone unfolding")
    timesteps = int(particles[0].properties.shape[0])


    #create array for storing separation for each molecule
    while n < len(molecules):
        molecules[n].separation = np.zeros((timesteps, 1))
        n += 1

    # get separation of sticker particles over time. double-check they are same molecule
    while m < len(molecules):
        print("LOOP NUMBER ", loopNum)
        loopNum += 1
        if row == timesteps:
            break
        if particles[p].properties[row][5] != particles[p + 1].properties[row][5]:
            raise ValueError("Particles in the list have not been ordered according to their molecule")
        elif particles[p].properties[row][5] == particles[p + 1].properties[row][5]:
            print("setting up")
            molecules[m].molNum = particles[p].properties[row][5]
            print("mol num", molecules[m].molNum)
            dx = particles[p].properties[row][1] - particles[p + 1].properties[row][1]
            dy = particles[p].properties[row][2] - particles[p + 1].properties[row][2]
            dz = particles[p].properties[row][3] - particles[p + 1].properties[row][3]
            molecules[m].separation[row] = np.sqrt(dx**2 + dy**2 + dz**2)
            print("found sep", molecules[m].separation[row])
            p += 2
            m += 1
            row += 1
            if p == len(particles):
                p = 0
                m = 0

    return molecules

def plotAvUnfold(molecules, Nsteps, Nwrite):
    timesteps = int(Nsteps/Nwrite)
    simTime = []
    for t in range(0, timesteps + 1):
        simTime.append(1000 * t)
    data = np.zeros((timesteps + 1, int(len(molecules))))
    for num, m in enumerate(molecules):
    #    print("mol sep", m.separation)
        data[:, num] = m.separation.flatten()
    sumSep = np.sum(data, axis = 1)
    avSep = [sep / int(len(molecules)) for sep in sumSep]

    plt.plot(simTime, avSep)
    return plt.show()

with open("dillParticles.pkl", "rb") as f:
    particles = pickle.load(f)
unfoldedPar = checkUnfolded(particles)
molecules = getSep(unfoldedPar)
data = plotAvUnfold(molecules, Nsteps, Nwrite)

