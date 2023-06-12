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
module load tqdm
module load Biopython/1.75-foss-2019b-Python-3.7.4
module load BLAST+/2.9.0-gompi-2019b
module load DIAMOND/0.9.30-GCC-8.3.0

path='/beegfs/fast/katzlab/PTL1/Genomes/'

srun -D ${path}Scripts python3 ${path}Scripts/wrapper.py -1 1 -2 5 --cds ${path}PTL1GenomesBatches/PTL1GenomesBatch2 -o ${path}Output/PTL1Genomes_OutputBatch2 --genetic_code Universal --databases ${path}Databases &
wait
