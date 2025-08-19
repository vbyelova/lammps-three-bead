# a code to see if the bell model really will work for my code.

import numpy as np
import pandas as pd
import MDAnalysis as md
import matplotlib.pyplot as plt

# load in molecule data 
timestep = 20000
increase = 1000 
skiplines = 9
columns = ["id", "type", "mol", "c_ssintra"]
# id type mol c_ssintra
df = pd.read_csv("../output/ssintra.dump", skipinitialspace = True, sep = "\s+", usecols = columns, skiprows = skiplines)
print(df.id)

for n in df.values:
    

# only use those where atom 1 and atom 3 both have 0 ssintra bonds
# iterate every 3?
# have additional column for time step 
# find separation of these molecules

# load in av. free energy of molecule
# can i do it per each molecule? 
