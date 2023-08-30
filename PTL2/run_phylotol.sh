#!/bin/bash
#SBATCH --job-name=YOUR_PROJECT_NAME ##change this to a shortened name of your project
#SBATCH --output=Run_phylotol.%j.out # Stdout (%j expands to jobId)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=4 ##change to number of srun when running multiple instances
#SBATCH --mem=150G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=email@smith.edu ##add your email address

module purge       #Cleans up any loaded modules

module use /gridapps/modules/all    #make sure module locations is loaded

module load slurm
module load ETE
module load Biopython
module load DIAMOND/0.9.30-GCC-8.3.0
module load MAFFT
module load BioPerl
module load RAxML
module load IQ-TREE/2.1.2-foss-2020a
export PATH=$PATH:/beegfs/fast/katzlab/grid_phylotol_setup/programs/standard-RAxML-master

parent='/beegfs/fast/katzlab/<fill in with your path>/' #do not edit this line before katzlab/, add your path starting with the name of your folder

#if you are running batches, you need an srun line for each batch!
srun -D ${parent} python3 ${parent}Scripts/phylotol.py --start raw --end trees --gf_list listofOGs.txt --taxon_list taxon_list.txt --data OutgroupR2Gs --output ${parent}Output_folder > Output_folder.out &
#srun -D ${parent} python3 ${parent}Scripts/phylotol.py --start raw --end trees --gf_list mei_batch2.txt --taxon_list taxon_list.txt --data OutgroupR2Gs --output ${parent}mei_output2 > mei_output2.out &
#srun -D ${parent} python3 ${parent}Scripts/phylotol.py --start raw --end trees --gf_list mei_batch3.txt --taxon_list taxon_list.txt --data OutgroupR2Gs --output ${parent}mei_output3 > mei_output3.out &
#srun -D ${parent} python3 ${parent}Scripts/phylotol.py --start raw --end trees --gf_list mei_batch4.txt --taxon_list taxon_list.txt --data OutgroupR2Gs --output ${parent}mei_output4 > mei_output4.out &

wait
