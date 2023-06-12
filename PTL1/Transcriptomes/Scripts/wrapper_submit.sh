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
module load tqdm
module load Biopython/1.75-foss-2019b-Python-3.7.4
module load BLAST+/2.9.0-gompi-2019b
module load DIAMOND/0.9.30-GCC-8.3.0
module load VSEARCH/2.21.1-GCC-10.3.0

parent='/beegfs/fast/katzlab/becky/PTL1/Transcriptomes/Forams/'

srun -D ${parent}Scripts python3 ${parent}Scripts/wrapper.py -1 1 -2 7 -x --assembled_transcripts ${parent}AssembledTranscripts -o ${parent}  -n ${parent}Conspecific.txt --genetic_code Universal &
#srun -D ${parent}HQ/Scripts python3 ${parent}HQ/Scripts/wrapper.py -1 2 -2 7 -x --assembled_transcripts ${parent}Plate4/Assembled_Transcripts -o ${parent}Plate4  -n ${parent}Plate4/Conspecific.txt --genetic_code ${parent}Plate4/Gcodes.txt &
#srun -D ${parent}HQ/Scripts python3 ${parent}HQ/Scripts/wrapper.py -1 2 -2 7 -x --assembled_transcripts ${parent}Plate7/Assembled_Transcripts -o ${parent}Plate7  -n ${parent}Plate7/Conspecific.txt --genetic_code ${parent}Plate7/Gcodes.txt &
#srun -D ${parent}HQ/Scripts python3 ${parent}HQ/Scripts/wrapper.py -1 1 -2 7 -x --assembled_transcripts ${parent}Plate11/Assembled_Transcripts -o ${parent}Plate11  -n ${parent}Plate11/Conspecific.txt --genetic_code ${parent}Plate11/Gcodes.txt &
#srun -D ${parent}HQ/Scripts python3 ${parent}HQ/Scripts/wrapper.py -1 2 -2 7 -x --assembled_transcripts ${parent}Plate18/Assembled_Transcripts -o ${parent}Plate18  -n ${parent}Plate18/Conspecific.txt --genetic_code ${parent}Plate18/Gcodes.txt &
wait
