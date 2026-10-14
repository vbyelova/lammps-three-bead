import multiprocessing as mp

from modules.calcPercolation import *
import subprocess
import os

def process_single_run(args):
    """function used to process as many percolation calculations as possible. Calls on a C++ function."""
    barrier, refoldBarrier, runNum, vf, numMol, boxLength, bondsPerAtom, prob, suffix, simType, strainFraction  = args
    if simType == "step strain":
        prefix = f"stepStrain{strainFraction}_"
    elif simType == "gelation":
        prefix = ""
    if bondsPerAtom == 2:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
    else:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
        if prob < 1:
            conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}_prob{prob}{suffix}"
    
    filename = f"Run{runNum}_{conditions}"
    output_path = f"../runs/{conditions}/{filename}/{prefix}analysis/percDimsPerFrame.txt"
    
    if not os.path.exists(output_path):
        if suffix == "NOPERCOLATION":
            return
        
        frameByFramePerc(barrier, refoldBarrier, runNum, vf, numMol, boxLength, bondsPerAtom, prob, suffix, simType, strainFraction)
        
        systemDataFile = f"../runs/{conditions}/{filename}/{prefix}output/systemData.txt"
        percInfoFile   = f"../runs/{conditions}/{filename}/{prefix}analysis/percinfo.txt"
        subprocess.run(["./addingData/addingData", systemDataFile, percInfoFile, output_path])
