#!/bin/bash
#SBATCH --job-name=meta033 ##change this to a shortened name of your project
#SBATCH --output=Run_phylotol.%j.out # Stdout (%j expands to jobId)
#SBATCH --nodes=1
#SBATCH --ntasks=10 ##change this to be double the number of task/batches you want to launch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=youremail@smith.edu ##add your email address

module purge       #Cleans up any loaded modules

module use /gridapps/modules/all    #make sure module locations is loaded

module load slurm
module load ETE
module load Biopython/1.79-foss-2021b
module load DIAMOND/2.0.13-GCC-11.2.0
module load MAFFT
module load BioPerl
module load RAxML
module load IQ-TREE/2.1.2-gompi-2021b
module load tqdm/4.64.1-GCCcore-12.2.0
module load Python/3.9.6-GCCcore-11.2.0
export PATH=$PATH:/beegfs/fast/katzlab/grid_phylotol_setup/programs/standard-RAxML-master

parent='/beegfs/fast/katzlab/Adri/p2PTL/033_meta/B1_meta_033/' #add your path starting with the name of your folder, should begin with /beegfs/fast/katzlab/

#if you are running batches, you need an srun line for each batch!
srun --exact -n 1 -D ${parent} python3 ${parent}Scripts/phylotol.py --similarity_filter --sim_cutoff 0.95 --sim_taxa sim_taxa.txt --blacklist GuidanceRemovedSeqs_allConservedRuns_ML_nov_dec_2023.txt --start raw --end trees --gf_list B1_listofOGs.txt --taxon_list taxon_list.txt --data OutgroupR2Gs --output ${parent}Output_folder_B1 > Output_folder_B1.out &
srun --exact -n 1 -D ${parent} python3 ${parent}Scripts/phylotol.py --similarity_filter --sim_cutoff 0.95 --sim_taxa sim_taxa.txt --blacklist GuidanceRemovedSeqs_allConservedRuns_ML_nov_dec_2023.txt --start raw --end trees --gf_list B2_listofOGs.txt --taxon_list taxon_list.txt --data OutgroupR2Gs --output ${parent}Output_folder_B2 > Output_folder_B2.out &
srun --exact -n 1 -D ${parent} python3 ${parent}Scripts/phylotol.py --similarity_filter --sim_cutoff 0.95 --sim_taxa sim_taxa.txt --blacklist GuidanceRemovedSeqs_allConservedRuns_ML_nov_dec_2023.txt --start raw --end trees --gf_list B3_listofOGs.txt --taxon_list taxon_list.txt --data OutgroupR2Gs --output ${parent}Output_folder_B3 > Output_folder_B3.out &
srun --exact -n 1 -D ${parent} python3 ${parent}Scripts/phylotol.py --similarity_filter --sim_cutoff 0.95 --sim_taxa sim_taxa.txt --blacklist GuidanceRemovedSeqs_allConservedRuns_ML_nov_dec_2023.txt --start raw --end trees --gf_list B4_listofOGs.txt --taxon_list taxon_list.txt --data OutgroupR2Gs --output ${parent}Output_folder_B4 > Output_folder_B4.out &
srun --exact -n 1 -D ${parent} python3 ${parent}Scripts/phylotol.py --similarity_filter --sim_cutoff 0.95 --sim_taxa sim_taxa.txt --blacklist GuidanceRemovedSeqs_allConservedRuns_ML_nov_dec_2023.txt --start raw --end trees --gf_list B5_listofOGs.txt --taxon_list taxon_list.txt --data OutgroupR2Gs --output ${parent}Output_folder_B5 > Output_folder_B5.out &
wait
