#!/bin/bash
#
#SBATCH --job-name=PTL1_genome
#SBATCH --output=PTL1.%j.out # Stdout (%j expands to jobId)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=64 ##change to number of srun when running multiple instances
#SBATCH --mem=160G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOUREMAIL@smith.edu

module purge       #Cleans up any loaded modules
module use /gridapps/modules/all    #make sure module locations is loaded

module load slurm
module load Biopython/1.75-foss-2019b-Python-3.7.4
module load BLAST+
module load DIAMOND/0.9.30-GCC-8.3.0

export PATH=$PATH:/Users/katzlab/scratch/katzlab/grid_phylotol_setup/programs/standard-RAxML-master
export PATH=$PATH:/Users/katzlab/scratch/katzlab/grid_vsearch_setup/vsearch-2.15.1-linux-x86_64/bin

python wrapper.py -1 1 -2 5 --cds ../TestData --genetic_code Universal --databases ../Databases