# python code to calculate percolation of three-bead network.


import pickle

# import code by Dr. David Head
from percPerBC import percPerBC
from systemData import Npar


with open("intermolBondsPerTimestep.pkl", "rb") as f:
    bonds = pickle.load(f)

finalFrame = len(bonds)
# vertices are the particles. this needs to be 1D.
V = [n for n in range(1, Npar + 1)]

# edges are the bonds. these are 3D where the first two values are
# the vertices and the third value is [x, y, z] where 0 = not bonded
# over boundary and 1 = bonded over boundary.

E = [n[-3:] for n in bonds.get(finalFrame)]

system = percPerBC(3)
print("Maximum spanning dimension of a component: {}".format(system.percolationDimension(V, E)))
#print(system)
