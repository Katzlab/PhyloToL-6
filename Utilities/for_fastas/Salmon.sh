#!/bin/bash

## Last updated on Jan 9th 2024 by Auden Cote-L'Heureux

#Intent: Calculate TPM for assembled transcripts
#Dependencies: None
#Inputs: Must be in a folder along with a folder called 'Transcriptomes', containing assembled transcripts as output by rnaSpades (transcripts.fasta), 
## and a folder called 'RawReads' containing the fwd and rev reads prior to assembly, with the same file prefixes as the corresponding assembled transcript files
#Outputs: A folder, containing a 'quant' file which has TPM data.

## If running on an HPC, include parameters here! For example, on a Slurm system you might use


#SBATCH --job-name=tpm
#SBATCH --output=Salmon.%j.out # Stdout (%j expands to jobId)
#SBATCH --nodes=1
#SBATCH --ntasks=60
#SBATCH --mem=60G

mkdir Indices

## First, build transcript indices

cd Transcriptomes

IFS='/'
for TRANS in *; do
        #read -a trapsplit <<<"$TRANS"
        #traf=${trapsplit[1]}
        tax=${TRANS:0:10}

	./../salmon-1.9.0_linux_x86_64/bin/salmon index -t $TRANS -i ../Indices/$tax
done

## Now calculate TPM

cd Indices

IFS='/'
for TRANS in *; do
        read -a trapsplit <<<"$TRANS"
        tax=${TRANS:0:10}
        fpe='NA'; rpe='NA'; fpesub="FPE"; rpesub="RPE"
        for TRIM in ../RawReads/*; do
                read -a tripsplit <<<"$TRIM"
                trif=${tripsplit[2]}
                if [ "${trif:0:10}" == "$tax" ]; then
                        if [[ "$trif" == *"$fpesub"* ]]; then
                                fpe=$trif
                        fi
                        if [[ "$trif" == *"$rpesub"* ]]; then
                                rpe=$trif
                        fi
                fi
        done
        
        if [ "$rpe" != 'NA' ]; then
                ./../salmon-1.9.0_linux_x86_64/bin/salmon quant -i $TRANS -l A -1 ../RawReads/$fpe -2 ../RawReads/$rpe --validateMappings -o ../quants/$tax
        fi

        if [ "$rpe" == 'NA' ]; then
                ./../salmon-1.9.0_linux_x86_64/bin/salmon quant -i $TRANS -l A -r ../RawReads/$fpe --validateMappings -o ../quants/$tax
        fi

done