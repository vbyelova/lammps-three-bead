# code to generate input files for lammps.

import numpy as np
from numpy import random
from random import randint

def generateLammpsInput(conditions, filename, boxLength, numMol, prob, bondsPerAtom):
    Nsteps = 700000 
    equilTime = 120000
    r_cutoff = 1.112462048
    seed = randint(100000, 999999)
    comm = 10
    numPar = numMol * 3
    tstep = 0.001
    halfLength = int(0.5 * boxLength)

    print("writing input file..")
    # write system parameters to file
    with open(f"../runs/{conditions}/{filename}/output/systemData.txt", "w") as f:
        f.write(f"[systemData]\n")
        f.write(f"boxLength = {boxLength}\n")
        f.write(f"Npar = {numPar}\n")
        f.write(f"Nsteps = {Nsteps}\n")
        f.write("Nwrite = 1000\n")
        f.write(f"equilTime = {equilTime}\n")
        f.write(f"prob = {prob}\n")

    # actually write the lammps input file now
    with open(f"../runs/{conditions}/{filename}/input/in.lammps", "w") as f:
        # set some simulation variables
        f.write(f"# simulation variables\n\n")
        f.write(f"variable Nsteps         equal {Nsteps}\n")
        f.write(f"variable Nwrite         equal 1000\n")
        f.write(f"variable equilTime      equal {equilTime}\n")
        f.write(f"variable seed           equal {seed}\n")
        f.write(f"variable tstep          equal {tstep}\n")
        f.write(f"variable boxLength      equal lx\n")
        f.write(f"variable r_cutoff       equal {r_cutoff}\n")
        f.write(f"variable s              equal logfreq3(10,100,{Nsteps})\n")
        f.write(f"variable prob           equal {prob}\n\n")

        # now some configuration stuff
        f.write(f"# configure simulation\n\n")
        f.write(f"units lj\n")
        f.write(f"boundary p p p\n")
        f.write(f"atom_style angle\n\n")
        f.write(f"read_data {filename}.in"
                f" extra/bond/per/atom {bondsPerAtom} extra/special/per/atom {bondsPerAtom}00\n")
        f.write(f"log input_log.lammps\n")
        f.write("min_style fire\n\n")

        # print thermodynamic info
        f.write(f"thermo 1000\n")

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
        f.write(f"pair_style soft {r_cutoff}\n")
        f.write(f"pair_coeff * * 10\n\n")

        f.write(f"minimize 1.0e-4 1.0e-6 1000 10000\n")
        f.write(f"velocity all create 1 {seed} mom yes rot yes dist gaussian\n\n")

        f.write(f"dump 1 all custom 1000 ../output/dump.lammpstrj"
                f" id x y z type mol\n")
        f.write(f"dump_modify 1 sort id\n\n")

        f.write(f"timestep 0.0001\n")
        f.write(f"run {equilTime}\n")
        f.write(f"unfix nvelim\n")
        f.write(f"unfix thermoequil\n\n")

        f.write(f"fix pressequil all npt temp 1 1 10 iso 1 1 10\n")
        f.write(f"run {equilTime}\n")
        f.write(f"unfix pressequil\n\n")

        # now defining main simulation parameters
        f.write(f"# main simulation parameters\n\n")
        f.write(f"pair_style hybrid/overlay zero 5 lj/cut {r_cutoff}\n")
        f.write(f"pair_coeff * * lj/cut 1.0 1.0\n")
        f.write(f"pair_coeff * * zero 5\n")
        f.write(f"pair_modify shift yes\n\n")

        f.write(f"include ../../angleInfo.in\n")
        f.write(f"compute angles all angle/local theta\n")
        f.write(f"compute angle_pe all angle/local eng\n\n")

        f.write(f"fix bd_nve all nve\n")
        f.write(f"fix bd_langevin all langevin 1 1 10 {seed}\n")
        f.write(f"fix percolate all bond/create 1 1 1 {r_cutoff} 2 iparam {bondsPerAtom}"
                f" 1 jparam {bondsPerAtom} 1 molecule inter prob {prob} {seed}\n\n")
        
        # a bunch of computes
        f.write(f"# compute commands for analysis\n\n")
        f.write(f"compute intra all nbond/atom bond/type 1\n")
        f.write(f"compute inter all nbond/atom bond/type 2\n")
        f.write(f"compute totalbonds all count/type bond\n")
        f.write(f"variable totalintra equal c_totalbonds[1]\n")
        f.write(f"variable totalinter equal c_totalbonds[2]\n")
        f.write(f"variable allbonds equal v_totalintra+v_totalinter\n\n")

        f.write(f"compute bondedatomid1 all property/local batom1\n")
        f.write(f"compute bondedatomid2 all property/local batom2\n")
        f.write(f"compute btype all property/local btype\n")
        f.write(f"compute bondforce all bond/local fx fy fz\n")
        f.write(f"compute bonddir all bond/local dx dy dz\n")
        f.write(f"compute bondlen all bond/local dist\n")
        f.write(f"variable vol equal lx*ly*lz\n\n")

        # more dumps

        f.write(f"# dump commands\n\n")
        f.write(f"undump 1\n")
        f.write(f"reset_timestep 0\n")
        f.write(f"dump 1 all custom 1000 ../output/dump.lammpstrj"
                f" id x y z type mol\n")
        f.write(f"dump_modify 1 sort id\n")
        f.write(f"dump_modify 1 every v_s\n\n")

        f.write(f"compute myrdf all rdf {halfLength} 1 1\n")
        f.write(f"fix rdfstuff all ave/time 1000 1 1000 c_myrdf[*] "
                f"file ../output/rdf.dat mode vector\n\n")

        f.write(f'fix nbondsfile all print v_s "${{allbonds}} " file'
                f' ../output/nbonds.dat screen no\n')
        f.write(f"dump bondinfo all local 1000 ../output/bondinfo.dat"
                f" index c_bondedatomid1 c_bondedatomid2"
                f" c_bondforce[*] c_bonddir[*] c_bondlen\n")

        f.write(f"dump angles all local 1000 ../output/moleculeangles.dat"
                f" index c_angles c_angle_pe\n\n")
        
        f.write(f"dump_modify bondinfo every v_s\n\n")
        f.write(f"dump_modify angles every v_s\n\n")


        # let's run the main sim
        f.write(f"# main simulation\n\n")
        f.write(f"timestep {tstep}\n")
        f.write(f"neighbor 0.5 bin\n")
        f.write(f"neigh_modify every 1 delay 2 check yes\n")
        f.write(f"run {Nsteps}\n")

    print(f"file written to ../runs/{conditions}/{filename}")
