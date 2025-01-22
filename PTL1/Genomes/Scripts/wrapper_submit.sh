## Last updated Jan 2025 by Auden Cote-L'Heureux

## This script is intended to be used to process genomic CDS with EukPhylo part 1 on an HPC that uses the Slurm workload manager.
## The first part of the script are Slurm-specific parameters that should be adjusted by users to fit their resource allocation
## needs and restrictions, followed by some example commands taken from the GitHub Wiki, more detail for which can be found
## here: https://github.com/Katzlab/EukPhylo/wiki/EukPhylo-Part-1:-GF-assignment


#!/bin/bash

## Slurm specific code

#SBATCH --job-name=EukPhylo
#SBATCH --output=EukPhylo.%j.out # Stdout (%j expands to jobId)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=64 # #change to number of srun when running multiple instances
#SBATCH --mem=160G

module purge       #Cleans up any loaded modules
module load slurm
module load tqdm
module load Biopython/1.75-foss-2019b-Python-3.7.4
module load BLAST+/2.9.0-gompi-2019b
module load DIAMOND/0.9.30-GCC-8.3.0

path='/Your/Home/Folder'

## Example run command

# Start at script 1 and go through script 5 (the final script) using the Universal genetic code
srun -D ${path}Scripts python3 ${path}Scripts/wrapper.py -1 1 -2 5 --cds ${path}Input -o ${path}Output --genetic_code Universal --databases ${path}Databases 

