### Three-bead particle simulation in LAMMPS

Changes to code can be done in text editor or any coding GUI.
Code is run on the command line.

***Start with input directory***
**generateThreeBead.py**
1. Change box size and number of particles to desired volume fraction
2. Run program
 `$ python generateThreeBead.py`


**in.lammps**
1. Change simulation run time
2. Change sigma to desired volume fraction
3. Run in.lammps
`$ lmp -in in.lammps` 
NOTE: this line may be different depending on your LAMMPS build.

***Move to analysis directory***
1. Parse outputs so that further analysis can be done in Python
`$ python parseOutputs.py`
