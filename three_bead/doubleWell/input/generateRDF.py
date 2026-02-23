# code to generate input files for lammps. this one is for calculating a bigger rdf

import numpy as np
from numpy import random
from random import randint

def generateRDF(conditions, filename, boxLength, numMol, prob, bondsPerAtom):
    Nsteps = 700000 
    equilTime = 120000
    r_cutoff = 1.112462048
    seed = randint(100000, 999999)
    comm = 0.55 * boxLength
    numPar = numMol * 3
    tstep = 0.001
    halfLength = int(0.5 * boxLength)

    # actually write the lammps input file now
    with open(f"../runs/{conditions}/{filename}/input/rdf.lammps", "w") as f:
        # set some simulation variables
        f.write(f"# simulation variables\n\n")
        f.write(f"variable Nsteps         equal {Nsteps}\n")
        f.write(f"variable Nwrite         equal 1000\n")
        f.write(f"variable equilTime      equal {equilTime}\n")
        f.write(f"variable seed           equal {seed}\n")
        f.write(f"variable tstep          equal {tstep}\n")
        f.write(f"variable boxLength      equal lx\n")
        f.write(f"variable r_cutoff       equal {r_cutoff}\n")
        f.write(f"variable prob           equal {prob}\n\n")

        # now some configuration stuff
        f.write(f"# configure simulation\n\n")
        f.write(f"units lj\n")
        f.write(f"boundary p p p\n")
        f.write(f"atom_style angle\n\n")
        f.write(f"read_data {filename}.in"
                " extra/bond/per/atom 2 extra/special/per/atom 200\n")
        f.write(f"log log_rerun.lammps\n\n")

        # equilibration stuff
        f.write(f"# equilibration parameters\n\n")
        f.write(f"bond_style harmonic\n")
        f.write(f"bond_coeff 1 3 {r_cutoff} #intramolecular\n")
        f.write(f"bond_coeff 2 3 {r_cutoff} #intermolecular\n\n")
        f.write(f"fix nvelim all nve/limit 1.0\n")
        f.write(f"fix thermoequil all langevin 1 1 10 {seed}\n\n")
        f.write(f"angle_style harmonic\n")
        f.write(f"angle_coeff 1 1000 60\n")
        f.write(f"comm_modify cutoff {comm}\n\n")

        # now defining main simulation parameters
        f.write(f"# main simulation parameters\n\n")
        f.write(f"pair_style lj/cut {0.5 * boxLength}\n")
        f.write(f"pair_coeff * * 1.0 1.0\n")


        f.write(f"include ../../angleInfo.in\n")

        # let's run the main sim
        f.write(f"# main simulation\n\n")
        f.write(f"timestep {tstep}\n")
        f.write(f"neighbor 0.5 bin\n")
        f.write(f"neigh_modify page 200000 one 20000\n")
        f.write(f"read_dump ../output/dump.lammpstrj 0 x y z box yes\n")
        f.write(f"reset_timestep 0\n")
        f.write(f"compute myrdf all rdf 100 1 1\n")
        f.write(f"fix a all ave/time 500 1 500 c_myrdf[*] file "
        f"../output/bigRDF.dat "
        f"mode vector\n")
        f.write(f"rerun ../output/dump.lammpstrj first 0 "
        f"every 500 dump x y z box yes format native\n")
        f.write(f"unfix a\n")


        
    print(f"file written to ../runs/{conditions}/{filename}")
