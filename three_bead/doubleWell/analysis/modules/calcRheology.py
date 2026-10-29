import numpy as np
import re
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as mcolors
import matplotlib.patches as patches
from itertools import cycle
from collections import defaultdict
from .parseDump import *
from scipy.optimize import curve_fit
import pickle
import pandas as pd

def calcAvShearModulus(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, prob, suffixes, simType, strainFraction):
    """gets out the average stress and strain, uses to calculate shear modulus at end of sim"""
    allShearModuli = defaultdict(list)
    allShearModuliErr = defaultdict(list)
    avShearModulus = {}
    avShearModulusErr = {}
    if simType == "step strain":
            prefix = f"stepStrain{strainFraction}_"
    elif simType == "gelation":
        prefix = ""
    for barrier, suffix in zip(unfoldBarriers, suffixes):
        key = (barrier, suffix)
        if bondsPerAtom == 2:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
        else:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
            if prob < 1:
                conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}_prob{prob}{suffix}"
        with open(f"../runs/{conditions}/{prefix}timesteps.pkl", "rb") as f:
            timesteps = pickle.load(f)

        with open(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/{prefix}data/avStressAndStrain.pkl", "rb") as f:
            avStresses, avStressesErr, avStrains, avStrainsErr = pickle.load(f)
        plt.plot(timesteps, avStresses[barrier])
        plt.show()
        for w in range(int(len(timesteps)/2), len(timesteps)):
            shearModulus = avStresses[barrier][w] / avStrains[barrier][w]
            frameErr = np.abs(shearModulus) * np.sqrt((avStressesErr[barrier][w] / avStresses[barrier][w])**2 + (avStrainsErr[barrier][w] / avStrains[barrier][w])**2)
            allShearModuli[barrier].append(shearModulus)
            allShearModuliErr[barrier].append(frameErr)
            print(w, shearModulus)

        shearModuliArray = np.array(allShearModuli[barrier])
        shearModuliErrArray = np.array(allShearModuliErr[barrier])
        avShearModulus[barrier] = np.mean(shearModuliArray)
        avShearModulusErr[barrier] = np.mean(shearModuliErrArray)
        print(avStresses[barrier][-1], avStrains[barrier][-1])
        print(f"shear modulus in final frame for strain = {strainFraction} and E_U = {barrier}kT: {avShearModulus[barrier]} +- {avShearModulusErr[barrier]}")
    return avShearModulus, avShearModulusErr


def plotShearModuli(unfoldBarriers, refoldBarrier, numRuns, vf, numMol, boxLength, bondsPerAtom, prob, suffixes, simType, strainFraction):
    if simType == "step strain":
        prefix = f"stepStrain{strainFraction}_"
    elif simType == "gelation":
        prefix = ""
    fig, ax = plt.subplots()
    for barrier, suffix in zip(unfoldBarriers, suffixes):
        key = (barrier, suffix)
        if bondsPerAtom == 2:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
        else:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
            if prob < 1:
                conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}_prob{prob}{suffix}"
        with open(f"../runs/{conditions}/{prefix}timesteps.pkl", "rb") as f:
            timesteps = pickle.load(f)
        with open(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/{prefix}data/shearModuli.pkl", "rb") as f:
            avShearModuli, avShearModuliErr = pickle.load(f)

        ax.errorbar(barrier, avShearModuli[barrier], yerr = avShearModuliErr[barrier], label = f"$E_U$ = {barrier}kT")

    plt.savefig(f"../runs/boxLength{boxLength}/vf{vf}/bondsPerAtom{bondsPerAtom}/prob{prob}/shearModuli.jpg")
