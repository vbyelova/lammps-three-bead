# python code setting up MDAnalysis for post-mortem analysis for LAMMPS

import MDAnalysis as mda
import numpy as np

u = mda.Universe("../output/dump.lammpstrj", convert_units =  False, format = "LAMMPSDUMP", timeunit = None,
                  lengthunit = None, additional_columns = ['type', 'mol', 'c_ssintra'])
print("WARNING: Time and length units are set to None so that MDAnalysis doesn't get angry. LAMMPS LJ sims are unitless.")

for ts in u.trajectory:
    types = ts.data['type']
    n_mol = ts.data['mol']
    n_intra = ts.data['c_ssintra']

print(types, n_mol, n_intra)
