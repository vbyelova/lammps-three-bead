import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

def parseLammpsPosRDF(barrier, refoldBarrier, runNum, Vf, numMol):

    conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{Vf}_mol{numMol}"
    filename = f"Run{runNum}_{conditions}"
    rdf = defaultdict(list)
    smolrdf = defaultdict(list)
    with open(f"../runs/{conditions}/{filename}/output/rdf.dat", "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rsplit()
            if len(parts) == 2:
                frame = int(parts[0])
            if len(parts) == 4:
                distance = float(parts[1])
                gofr = float(parts[2])
                smolrdf[frame].append([distance, gofr])
    with open(f"../runs/{conditions}/{filename}/output/bigRDF.dat", "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rsplit()
            if len(parts) == 2:
                frame = int(parts[0])
            if len(parts) == 4:
                distance = float(parts[1])
                gofr = float(parts[2])
                rdf[frame].append([distance, gofr])

    plt.plot(np.array(rdf[frame])[:, 0], np.array(rdf[frame])[:, 1], "-", color = "blue")
    plt.plot(np.array(smolrdf[frame])[:, 0], np.array(smolrdf[frame])[:, 1], "-", color = "red")
    print(np.array(rdf[frame]))
    plt.xlabel("distance")
    plt.ylabel("g(r)")
    plt.show()
    return

def posRDF(barrier, refoldBarrier, runNum, Vf, numMol, particles):
    """calculates the positional radial distribution function"""
    return


def calcAvPosRDF(unfoldBarriers, refoldBarrier):
    return
