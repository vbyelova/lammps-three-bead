import matplotlib.pyplot as plt
import numpy as np
from numpy import random
from scipy.optimize import minimize, root

from doubleWellCoeffs import *
from generateThreeBead import *

#generate optimal values for quartic function

# first input value is unfolding barrier height
# second input value is minimum energy for unfolded state

barrierToUnfold = 3
barrierToRefold = 2

quarticVals = calcDoubleWellCoeffs(barrierToUnfold,barrierToRefold)

angle0, k2, k3, k4 = quarticVals[0], quarticVals[1], quarticVals[2], quarticVals[3]

writeAngleFile(angle0, k2, k3, k4)

#generate molecule input file

#first value is number of molecules
#second value is length of box

numMol = 360
boxLength = 20
generateThreeBead(numMol, boxLength)
