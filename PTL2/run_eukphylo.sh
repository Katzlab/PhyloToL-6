## Last updated Jan 2025 by Auden Cote-L'Heureux

## This shell script is used for running EukPhylo part 2, and includes a general setup for use on an HPC that uses 
## the Slurm workload manager. It also includes several example run commands, which correspond to examples explained in more detail in the 
## EukPhylo Wiki (https://github.com/Katzlab/EukPhylo/wiki/EukPhylo-Part-2:-MSAs,-trees,-and-contamination-loop).
## These run commands can also be copied and run in the terminal / command line separately, without a shell script.


#!/bin/bash

## SLURM-SPECIFIC SETUP BELOW

#SBATCH --job-name=EukPhylo # Job name
#SBATCH --output=Run_EukPhylo.%j.out # Stdout (%j expands to jobId)
#SBATCH --nodes=1
#SBATCH --ntasks=10 ## On the Smith College HPC (Grid), we have to change this to be double the number of task/batches you want to launch
#SBATCH -c 24 # Number of Cores per Task
#SBATCH --mem=125G # Requested Memory
#SBATCH -q long # Partition. Only use this on certain HPCs (e.g., Unity at UMass).
#SBATCH -t 336:00:00 # Job time limit

module purge       #Cleans up any loaded modules
module use /gridapps/modules/all    #make sure module locations is loaded
module load slurm
module load ETE
module load Biopython/1.79-foss-2021b
module load DIAMOND/2.0.13-GCC-11.2.0
module load MAFFT
module load RAxML
module load IQ-TREE/2.1.2-gompi-2021b
module load tqdm/4.64.1-GCCcore-12.2.0
module load Python/3.9.6-GCCcore-11.2.0
export PATH=$PATH:/Path/To/Executable/Files

parent='/Your/Home/Folder/' # The folder where you are running EukPhylo (this should contain the Scripts and input data folders

## RUN COMMANDS BELOW

# A simple run of part 2, starting from ReadyToGo files and running through tree building
srun --exact -n 1 -D ${parent} python3 ${parent}Scripts/eukphylo.py --start raw --end trees --gf_list  ${parent}listofOGs.txt --taxon_list ${parent}taxon_list.txt --data ${parent}Input_folder --output ${parent}Output_folder

# Another example run starting from ReadyToGo files and running through tree building, with the commonly used similarity filter cutoff, blacklist, and  "sim_taxa_list" arguments (see Wiki)
srun --exact -n 1 -D ${parent} python3 ${parent}Scripts/eukphylo.py --start raw --end trees --gf_list  ${parent}listofOGs.txt --taxon_list ${parent}taxon_list.txt --data ${parent}Input_folder --output ${parent}Output_folder --similarity_filter --blacklist ${parent}Blacklist.txt --sim_cutoff 0.99 --sim_taxa sim_taxa_list.txt

# An example of running just the concatenation step of part 2, starting from trees
srun --exact -n 1 -D ${parent} python eukphylo.py --start trees --concatenate --concat_target_taxa Sr_rh --data Output

# See the Wiki (https://github.com/Katzlab/EukPhylo/wiki/EukPhylo-Part-2:-MSAs,-trees,-and-contamination-loop) for more details!




