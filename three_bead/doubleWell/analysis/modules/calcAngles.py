import numpy as np
import re
import matplotlib.pyplot as plt

from collections import defaultdict
from .parseDump import *

def calcAngles(barrier, refoldBarrier, Vf, numMol):
    """saves angles of three bead molecules in each simulation frame and returns
        the final frame."""

    conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{Vf}_mol{numMol}"
    filename = f"Run0_{conditions}"
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
    
    return angles

def unfoldedMolecules(barrier, refoldBarrier, runNum, Vf, numMol, angles):
    """makes a list of molecules that are unfolded (have an angle of 120-180)"""
    # let's say a particle is unfolded if it's around 120-180 degrees
    # based on our bimodal distribution
    counter = 0
    unfoldedMols = defaultdict(list)
    while counter < len(angles):
        for n in angles[counter]:
            if n > 120:
                unfoldedMols[counter].append(n)
            counter += 1




def bimodalAngle(barrier, refoldBarrier, runNum, Vf, numMol):
    """plots a histogram of the final angles of the molecules in the system."""

    conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{Vf}_mol{numMol}"
    filename = f"Run{runNum}_{conditions}"
    angles = calcAngles(3, refoldBarrier, Vf, numMol)
    finalFrameAngles = angles[int(len(angles)-1)]
    hist, bins = np.histogram(finalFrameAngles, bins = 30)
    logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
    plt.hist(finalFrameAngles, bins = logbins, color = "pink")
    plt.legend([f"Vf = {Vf}\nunfolding barrier  = {barrier} kT\n num. mol. = {numMol}"])
    plt.xlabel("molecule angle \u03B8")
    plt.ylabel("log(number of molecules)")
    plt.title("Bimodal distribution of folded and unfolded molecules")
    plt.savefig(f"../runs/{conditions}/{filename}/analysis/figs/bimodalAngleDist")

    return

def unfoldingRDF(barrier, refoldBarrier, runNum, Vf, numMol, boxLength):
    """find where unfolding is happening through a discrete radial distribution function."""

    # for each molecule:
    #   find how many unfolded mol are dr away and divide that number by volume of shell
    #   divide that by the number of mol in sim / volume of box
    #   repeat for multiple dr
    # plot

def unfoldingWithPerc(barrier, refoldBarrier, runNum, Vf, numMol, nBonds):
    """plots number of unfolded particles against number of bonds in system"""
