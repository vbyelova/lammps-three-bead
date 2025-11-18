# python code to calculate percolation of three-bead network.


import pickle

# import code by Dr. David Head
from percPerBC import percPerBC
from systemData import Npar
from parseForBondVis import totalBonds

with open("dillPercolation.pkl", "rb") as f:
    bonds = pickle.load(f)

# vertices are the particles. this needs to be 1D.
V = [n for n in range(0, Npar + 2)]

# edges are the bonds. these are 3D where the first two values are
# the vertices and the third value is [x, y, z] where 0 = not bonded
# over boundary and 1 = bonded over boundary.
totalFrames = len(totalBonds)
E = bonds.get(totalFrames - 1)

def testPercPerBC():
    V = [n for n in range(0, 13)]
    E = [
         [0, 1, [0, 0, 0]],
         [1, 2, [0, 0, 0]],
         [2, 3, [0, 0, 0]],
         [3, 4, [0, 0, 0]],
         [4, 0, [0, 1, 0]],
         [2, 5, [0, 0, 0]],
         [5, 6, [0, 0, 0]],
         [6, 7, [1, 0, 0]],
         [7, 8, [0, 0, 0]],
         [8, 2, [0, 0, 0]],
         [9, 2, [0, 0, 0]],
         [9, 10, [0, 0, 0]],
         [10, 11, [0, 0, 1]],
         [11, 12, [0, 0, 0]],
         [12, 2, [0, 0, 0]]
        ]

    system = percPerBC(3)
    print( "Maximum spanning dimension of a component: {}".format(system.percolationDimension(V,E) ) )
    print( system )

system = percPerBC(3)
print("Maximum spanning dimension of a component: {}".format(system.percolationDimension(V, E)))
#print(system)

#testPercPerBC()
