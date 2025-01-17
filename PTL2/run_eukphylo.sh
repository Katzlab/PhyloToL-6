###GA 11/11/24
###Updated run script to include grid and unity commands. 
###The first block of code is specific to the grid. The second block is specific to unity. Pick one and delete the other.

#!/bin/bash
#SBATCH --job-name=GA1 # Job name
#SBATCH --output=Run_phylotol.%j.out # Stdout (%j expands to jobId)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=youremail@smith.edu ##add your email address



###Grid start
#SBATCH --nodes=1
#SBATCH --ntasks=10 ##change this to be double the number of task/batches you want to launch

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
export PATH=$PATH:/beegfs/fast/katzlab/grid_phylotol_setup/programs/standard-RAxML-master
parent='/beegfs/fast/katzlab/' #add your path starting with the name of your folder, should begin with /beegfs/fast/katzlab/
#if you are running batches, you need an srun line for each batch!
srun --exact -n 1 -D ${parent} python3 ${parent}Scripts/eukphylo.py --similarity_filter --sim_cutoff 0.95 --sim_taxa sim_taxa.txt --blacklist GuidanceRemovedSeqs_allConservedRuns_ML_nov_dec_2023.txt --start raw --end trees --gf_list B1_listofOGs.txt --taxon_list taxon_list.txt --data OutgroupR2Gs --output ${parent}Output_folder_B1 > Output_folder_B1.out &
wait
###Grid end


###Unity start
#SBATCH -c 24 # Number of Cores per Task
#SBATCH --mem=125G # Requested Memory
#SBATCH -q long # Partition
#SBATCH -t 336:00:00 # Job time limit
module purge       #Cleans up any loaded modules
module load miniconda/22.11.1-1
module load mafft/7.481
module load conda/latest
conda activate /work/pi_lkatz_smith_edu/Conda_PTL6p2/envs/PTL/
parent='/work/pi_lkatz_smith_edu/' #add your path startin>
#if you are running batches, you need an srun line for each batch!
srun -D ${parent} python3 ${parent}Scripts/eukphylo.py --similarity_filter --sim_cutoff 0.99 -->
wait
###Unity end
