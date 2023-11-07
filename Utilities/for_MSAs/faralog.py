#Author, date: Godwin Ani, last updated Nov 7th 2023
#Motivation: Produce gappyness stats for paralogs in an alignment.
#Intent: To use the gappyness stats in filtering paralogs.
#Inputs: A folder of alignment files
#Outputs: A spreadsheet of the gappyness stats
#Example: python faralog.py --alignment /Path/to/alignments  --code ten digit code

#Dependencies
import os, sys, re
from Bio import SeqIO
from tqdm import tqdm
import itertools
import pandas as pd
import argparse



parser = argparse.ArgumentParser(
                    prog = 'F/paralog assessment script',
                    description = "Updated November 7, 2023 by Godwin Ani.")

parser.add_argument('-a', '--alignment', help = 'The alignment files folder')
parser.add_argument('-c', '--code', help = 'The ten digit code for f/paralog assessment')
args = parser.parse_args()


def faralog_gaps():
    df = pd.DataFrame()
    #looping through the fasta files in the folder.   
    for file in tqdm(os.listdir(args.alignment)):
        if file.endswith('.fasta'):
            #creating empty lists to store the id, sequence, and splitted sequences of the fasta files.
            name =[]
            seq = []
            seq_len = []
            outer = []
            total = []
            split = []
            #reading the fasta files with Biopython (looping each sequence in a file and populating the empty name and seq lists).
            for x in SeqIO.parse(args.alignment + '/' + file, "fasta"):
                if x.id.startswith(args.code):
                    name.append(x.id)
                    seq.append(x.seq)
                    seq_len.append(len(x.seq))
                    total.append(x.seq.count('-'))
            #looping through the seq list, converting the Biopython seq object to string, and splitting it with all possible bases to get a list of gaps.
            for x in seq:
                split.append(re.split(r'[A-Z]', str(x)))
            #looping through the gap lists and summing the first splitted and last splitted gaps in the sequence to get the terminal gaps length. 
            for x in split:
                outer.append(len(x[0] + x[-1]))
            internal = [total - outer for total, outer in zip(total, outer)]
            data = {'Sequence': name, 'Terminal_gaps': outer, 'Internal_gaps':internal,'Aln_length':seq_len} 
            df1 = pd.DataFrame(data)
            df1['Total'] = df1['Terminal_gaps'] + df1['Internal_gaps']
            df1['Total'] = df1['Total']/df1['Aln_length']
            df1['OG'] = df1['Sequence'].str[-10:]
            df1['Avg_terminal_gaps_for_OG'] = df1['Terminal_gaps'].mean()
            df1['Avg_internal_gaps_for_OG'] = df1['Internal_gaps'].mean()
            df1['#term/avgTerm'] = df1['Terminal_gaps']/df1['Avg_terminal_gaps_for_OG']
            df1['#int/avgint'] = df1['Internal_gaps']/df1['Avg_internal_gaps_for_OG']
            best = df1[df1.Total == df1.Total.min()]
            best_terminal_gaps = best['Terminal_gaps'].iloc[0]  
            if best_terminal_gaps == 0:
                best_terminal_gaps += 1
            best_internal_gaps = best['Internal_gaps'].iloc[0]
            if best_internal_gaps == 0:
                best_internal_gaps += 1
            df1['Seq_term/best_seq_term'] = df1['Terminal_gaps']/best_terminal_gaps
            df1['Seq_int/best_seq_int'] = df1['Internal_gaps']/best_internal_gaps
            df1['Seq_term/avgTerm'] = df1['Terminal_gaps']/df1['Avg_terminal_gaps_for_OG']
            df1['Seq_int/avgInt'] =  df1['Internal_gaps']/df1['Avg_internal_gaps_for_OG']
            cond = best['Total'].iloc[0]
            df1['Best_seq'] = ['Yes' if x == cond else 'No' for x in df1['Total']]
            df = pd.concat([df,df1])
    df = df[['OG', 'Sequence', 'Best_seq', 'Terminal_gaps', 'Avg_terminal_gaps_for_OG','Seq_term/best_seq_term','Seq_term/avgTerm','Internal_gaps','Avg_internal_gaps_for_OG','Seq_int/best_seq_int','Seq_int/avgInt','Aln_length', 'Total']]
    df = df.round(3)
    df.to_csv(args.code + '.csv', index = False)
    
    
    
                
faralog_gaps()                


