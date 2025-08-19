# python code setting up MDAnalysis for post-mortem analysis for LAMMPS

import MDAnalysis as mda
from MDAnalysis.tests.datafiles import PSF, DCD, GRO, XTC
import numpy as np

#u = mda.Universe(LAMMPSDATA, atom_style = "id type mol x y z")

print(mda.Universe(PSF, DCD))
print("Using MDANalysis version", mda.__version__)
