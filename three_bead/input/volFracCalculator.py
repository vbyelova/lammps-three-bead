# input desired parameters to see number of particles required for a certain volume fraction.
# comment out the functions that you will not use

import numpy as np



def whatNumMol():
    sigma = float(input("Enter a value for sigma: "))
    volFrac = float(input("Enter the desired volume fraction: "))
    boxLen = float(input("Enter the box length: "))

    numPar = (boxLen**3 * volFrac) / ((sigma * 2 **(1 / 6))**3 * (np.pi / 6))
    numMol = numPar / 3

    return print("Number of MOLECULES in system should be ", numMol)

def whatVolFrac():
    """gives the volume fraction as a result of all parameters specified"""
    sigma = float(input("Enter a value for sigma: "))
    boxLen = float(input("Enter the box length: "))
    numMol = float(input("Enter the number of molecules: "))

    volFrac = ((sigma * 2 **(1 / 6))**3 * (np.pi / 6) * numMol * 3) / boxLen**3

    return print("Volume fraction of the system will be ", volFrac)

option = input("Select an option:\n \
                (1) Find number of particles for given volume fraction\n \
                (2) Find volume fraction for given sigma and particle num.\n ")

if option == "1":
    whatNumMol()

elif option == "2":
    whatVolFrac()
