# input desired parameters to see number of particles required for a certain volume fraction.
# comment out the functions that you will not use

import numpy as np



def whatNumPar():
    """spits out number of particles required for certain volume fraction"""

def whatSigma():
    """gives value for sigma required given a certain volume fraction and number of particles"""

def whatVolFrac():
    """gives the volume fraction as a result of all parameters specified"""

option = input("Select an option:\n \
                (1) Find number of particles for given volume fraction\n \
                (2) Find value of sigma for given vol. frac. and particle num.\n \
                (3) Find volume fraction for given sigma and particle num.\n ")

if option == 1:
    whatNumPar()
elif option == 2:
    whatSigma()
elif option == 3:
    whatVolFrac()
