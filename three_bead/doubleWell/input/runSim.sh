#!/bin/bash
#SBATCH --job-name=threebeadnetwork   # Job name
#SBATCH --time=10:00:00         # Request runtime (hh:mm:ss)
#SBATCH --mem=1G                # Request memory
#SBATCH --ntasks=1              # Number of tasks
#SBATCH --cpus-per-task=1       # Number of cores per task

# Load any necessary modules

# give omp resource allocation

# Execute your application
module load miniforge openmpi gcc

cd ../runs/unfold1_refold2_Vf0.07_mol3939
python runScripts.py
cd ../unfold1_refold2_Vf0.1_mol5627
python runScripts.py
cd ../unfold1_refold2_Vf0.15_mol8441
python runScripts.py
