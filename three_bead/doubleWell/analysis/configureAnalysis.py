import numpy as np
import pickle
import re
import os
import subprocess
from collections import defaultdict

from modules.parseDump import *
from modules.calcAngles import *
from modules.calcFracDim import * 
from modules.calcPercolation import *
from modules.calcStressTensor import *
from modules.calcCorrelation import *
from modules.threeBeadClasses import *
from modules.parseRDF import *
from modules.calcPorosity import *


#generate optimal values for quartic function

# first input value is unfolding barrier height
# second input value is minimum energy for unfolded state

barrierToUnfold = [1, 2, 3, 4, 5]
suffixes = ["", "", "", "", ""]
barrierToRefold = 2

numRuns = 10
vf = 0.07

boxLength = 100
prob = 1 # probability of unfolding
bondsPerAtom = 3

numMol = 31512
#generate molecule input file


for barrier, suffix in zip(barrierToUnfold, suffixes):
    conditions = f"unfold{barrier}_refold{barrierToRefold}_Vf{vf}_mol{numMol}{suffix}"

    for n in range (0, numRuns):
        filename = f"Run{n}_{conditions}"

        if suffix != "":
            continue
        else:
            if not os.path.exists(f"../runs/{conditions}/{filename}/analysis/parsePerc.py"):
                with open(f"../runs/{conditions}/{filename}/analysis/parsePerc.py", "w") as f:
                    f.write(f"if not os.path.exists(f'../runs/{conditions}/filename/analysis/percDimsPerFrame.txt'):\n")
                    f.write("f")



    with open(f"../runs/{conditions}/runScripts.py", "w") as f:
        f.write(f"# run this script from {conditions} to run the simulations.\n")
        f.write(f"import subprocess\n"
                f"import os\n\n")
        for filename in files:
            f.write(f"if not os.path.exists(f'../runs/{conditions}/{filename}/input/log.lammps'):\n")
            f.write(f"\tprint('running simulation {filename}..')\n")
            f.write(f"\tsubprocess.run(['lmp','-in','in.lammps'], cwd = '{filename}/input')\n\n")

