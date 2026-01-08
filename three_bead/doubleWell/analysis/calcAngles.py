import numpy as np
import re
import matplotlib.pyplot as plt

from collections import defaultdict
from parseBonds import totalBonds
from parseDump import *


def calcAngles(barrier, refoldBarrier, Vf, numMol):
    conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{Vf}_mol{numMol}"
    filename = f"threeBead_Run0_{conditions}"
    #totalBonds = totalBonds(f"../runs/{conditions}/{filename}/output/nbonds.dat")
    systemData = parseSystemData(f"../runs/{conditions}/{filename}/output/systemData.txt")
    boxLength, Npar, Nsteps, Nwrite, equilTime = systemData[0], systemData[1], systemData[2], systemData[3], systemData[4]
    
    angles = defaultdict(list)
    frame = 0
    counter = 0
    with open(f"../runs/{conditions}/{filename}/output/moleculeangles.dat", "r") as f:
        print("opened file")
        for line in f.readlines():
            line = line[:-1]
            #print(line)
            if re.search(r"\d+\s[\d.]+\s-?[\d.]+", line):
                #print(line, "matched")
                #angle = float(line.rsplit()[1]) * np.pi / 180
                angles[frame].append(float(line.rsplit()[1]))
                #print(line.rsplit()[1])
                counter += 1
                if counter == numMol:
                    #print("frame increased to ", frame)
                    frame += 1
                    counter = 0
    

    #with open(picklename, "wb") as f:
    #    pickle.dump(f)

    # print(angles[frame - 1])
    # plt.hist(angles[frame - 1], bins = 18)
    # plt.ylabel("number of molecules")
    # plt.xlabel("angle")
    # plt.title("Unfolding barrier = 5")
    # plt.show()

    return angles[frame - 1]


