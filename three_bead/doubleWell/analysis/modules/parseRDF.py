import matplotlib.pyplot as plt
import numpy as np
import pickle
import re
from collections import defaultdict

def parsePosRDF(barrier, refoldBarrier, runNum, Vf, numMol):
    conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{Vf}_mol{numMol}"
    filename = f"Run{runNum}_{conditions}"
    rdf = defaultdict(list)
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
                rdf[frame].append([distance, gofr])

    plt.scatter(np.array(rdf[frame])[:, 0], np.array(rdf[frame])[:, 1])
    plt.show()
    return
