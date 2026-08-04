import multiprocessing as mp

from modules.calcPercolation import *
import subprocess
import os

def process_single_run(args):
    barrier, suffix, refoldBarrier, vf, numMol, boxLength, bondsPerAtom, runNum = args
    
    if bondsPerAtom == 2:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}{suffix}"
    else:
        conditions = f"unfold{barrier}_refold{refoldBarrier}_Vf{vf}_mol{numMol}_bondsPerAtom{bondsPerAtom}{suffix}"
    
    filename = f"Run{runNum}_{conditions}"
    output_path = f"../runs/{conditions}/{filename}/analysis/percDimsPerFrame.txt"
    
    if not os.path.exists(output_path):
        if suffix == "NOPERCOLATION":
            return
        
        frameByFramePerc(barrier, refoldBarrier, runNum, vf, numMol, boxLength, bondsPerAtom)
        
        systemDataFile = f"../runs/{conditions}/{filename}/output/systemData.txt"
        percInfoFile   = f"../runs/{conditions}/{filename}/analysis/percinfo.txt"
        subprocess.run(["./addingData/addingData", systemDataFile, percInfoFile, output_path])
