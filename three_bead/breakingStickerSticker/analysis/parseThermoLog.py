# code to parse log file from lammps sim to get the thermodynamic info
# Victoria Byelova

import pickle
import re

def readLogFile(name):
    pattern = something
    with open(name) as text:
        for line in text.readlines():
            if re.search(pattern, line):
                line = line[:-1]
