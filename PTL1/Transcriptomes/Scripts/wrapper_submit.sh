#!/bin/bash
#
#SBATCH --job-name=PTL1_GBF
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
module load tqdm/4.62.3-GCCcore-11.2.0
module load Biopython/1.79-foss-2021b
module load BLAST+/2.12.0-gompi-2021b
module load DIAMOND/2.0.13-GCC-11.2.0
module load VSEARCH/2.22.1-GCC-11.3.0

parent='/beegfs/fast/katzlab/becky/PTL1/Transcriptomes/Forams/'

srun -D ${parent}Scripts python3 ${parent}Scripts/wrapper.py -1 1 -2 7 -x --assembled_transcripts ${parent}AssembledTranscripts -o ${parent}  -n ${parent}Conspecific.txt --genetic_code Universal &
#srun -D ${parent}Scripts python3 ${parent}Scripts/wrapper.py -1 1 -2 7 -x --assembled_transcripts ${parent}Assembled_Transcripts -o ${parent}  -n ${parent}Conspecific.txt --genetic_code ${parent}Gcodes.txt > log.out &
#srun -D ${parent}HQ/Scripts python3 ${parent}HQ/Scripts/wrapper.py -1 2 -2 7 -x --assembled_transcripts ${parent}Plate7/Assembled_Transcripts -o ${parent}Plate7  -n ${parent}Plate7/Conspecific.txt --genetic_code ${parent}Plate7/Gcodes.txt &
#srun -D ${parent}HQ/Scripts python3 ${parent}HQ/Scripts/wrapper.py -1 1 -2 7 -x --assembled_transcripts ${parent}Plate11/Assembled_Transcripts -o ${parent}Plate11  -n ${parent}Plate11/Conspecific.txt --genetic_code ${parent}Plate11/Gcodes.txt &
#srun -D ${parent}HQ/Scripts python3 ${parent}HQ/Scripts/wrapper.py -1 2 -2 7 -x --assembled_transcripts ${parent}Plate18/Assembled_Transcripts -o ${parent}Plate18  -n ${parent}Plate18/Conspecific.txt --genetic_code ${parent}Plate18/Gcodes.txt &
wait
