'''
#Author, date: Godwin Ani, 24th- Nov - 2023.
#Dependencies: Python3, Biopython
#Inputs: A folder of containing ReadyToGo files, PerSequenceStatSummaries, rules file and the script.
#The csv file should contain 10 digit codes and the limits(column headers are name, lower, and upper).
#Outputs: a folder of curated ready to go files.
#Example: python GC_identifier.py -i (input folder of r2gs), -r (rules csv file), -s (folder of per seq stats) 
'''


import os
import argparse
from Bio import SeqIO
os.makedirs('output', exist_ok= True)

def process_OG6_A_G(input_dir, rules_csv_path, stats_csv_dir):
    # Read rules from CSV into a dictionary
    rules = {}
    with open(rules_csv_path, encoding='utf-8') as rules_file:
        next(rules_file)  # Skip the header
        for line in rules_file:
            name, lower, upper = line.strip().split(',')
            rules[name] = {'lower': float(lower), 'upper': float(upper)}

    # Initialize dictionaries for storing sequences
    tot_neutral = {}
    tot_A = {}
    tot_G = {}

    # Process CSV files in 'PerSequenceStatSummaries' directory
    for file in os.listdir(stats_csv_dir):
        if file.endswith('.csv'):
            rule = rules.get(file[:10], {'lower': 0.0, 'upper': 1.0})
            lower, upper = rule['lower'], rule['upper']

            neutral_sequences, A_sequences, G_sequences = [], [], []
            with open(os.path.join(stats_csv_dir, file), encoding='utf-8') as csv_file:
                header = next(csv_file).strip().split(',')
                sequence_index, gc3_degen_index = header.index('Sequence'), header.index('GC3-Degen')
                for line in csv_file:
                    parts = line.strip().split(',')
                    sequence, gc3_degen = parts[sequence_index], float(parts[gc3_degen_index])
                    if lower <= gc3_degen <= upper:
                        neutral_sequences.append(sequence)
                    elif gc3_degen < lower:
                        A_sequences.append(sequence)
                    elif gc3_degen > upper:
                        G_sequences.append(sequence)
            tot_neutral[file[:10]] = neutral_sequences
            tot_A[file[:10]] = A_sequences
            tot_G[file[:10]] = G_sequences

    # Process sequence files in 'Input' directory
    for file in os.listdir(input_dir):
        sequences = []
        with open(os.path.join(input_dir, file), encoding='latin-1') as fasta_file:
            for record in SeqIO.parse(fasta_file, 'fasta'):
                file_key = file[:10]
                if record.id in tot_neutral[file_key]:
                    record.description = record.id
                    sequences.append(record)
                elif record.id in tot_A[file_key]:
                    # Create a new record for 'OGA' version
                    record_oga = record[:]
                    record_oga.id = f"{record.id.split('_OG6')[0]}_OGA{record.id.split('_OG6')[1]}"
                    record_oga.description = record_oga.id
                    sequences.append(record_oga)
                elif record.id in tot_G[file_key]:
                    # Create a new record for 'OGG' version
                    record_ogg = record[:]
                    record_ogg.id = f"{record.id.split('_OG6')[0]}_OGG{record.id.split('_OG6')[1]}"
                    record_ogg.description = record_ogg.id
                    sequences.append(record_ogg)
        SeqIO.write(sequences, 'output/' + file, 'fasta')

def main():
    parser = argparse.ArgumentParser(description='Process files based on rules.')
    parser.add_argument('-i', '--input_dir', required=True, help='Path to the input directory containing sequence files.')
    parser.add_argument('-r','--rules_csv_path', required=True, help='Path to the rules CSV file.')
    parser.add_argument('-s','--stats_csv_dir', required=True, help='Path to the directory containing statistics CSV files.')

    args = parser.parse_args()

    process_OG6_A_G(args.input_dir, args.rules_csv_path, args.stats_csv_dir)

if __name__ == "__main__":
    main()
