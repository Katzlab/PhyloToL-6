'''
Author: Auden Cote-L'Heureux, Laura Katz and ChatGPT
Last updated: June 9th, 2024
Motivation: Count the number of occurrences of each taxa in each OG in a post guidance file
Dependencies: Bio python, os, sys
Inputs: Directory of postguidance files 
Optional: use the --minor flag and include a file named focal_minors.txt in same folder (do not put file name in command line). This file should be csv of targets (Am_tu, Sr_rh, Sr_ci)
Optional: OR use the --major flag and include a file named focal_majors.txt in same folder (do not put file name in command line). This file should be csv of targets (Am, Sr)
Outputs: CSV file tallying all the counts of taxa in each OG file plus minor and major clade tallies
Command line: python3 CountTaxonOccurence_faster_minor.py  --input <dir of postguidance files> --minor
'''


import os
import sys
from Bio import SeqIO
import argparse


def get_args():
    parser = argparse.ArgumentParser(
        prog = 'Taxon occurrence counting script',
        description = "Updated June 9, 2024"
    )

    parser.add_argument('-i', '--input', type = str, required = True, help = 'Path to the folder containing the aligned/unaligned fasta files')
    parser.add_argument('--minor', action='store_true', help = 'Flag to use focal minor clades from focal_minors.txt')
    parser.add_argument('--major', action='store_true', help = 'Flag to use focal major clades from focal_majors.txt')
    args = parser.parse_args()
        
    if(args.input.endswith('/')):
        args.input = args.input[:-1]
        
    if(not os.path.isdir(args.input)):
        print('\nThe input folder (--input) could not be found. Make sure you have given the correct path.\n')
        exit()
                
    return args.input, args.minor, args.major


def count_tips(in_dir, use_focal_minors, use_focal_majors):
    
    focal_minors = []
    if use_focal_minors:
        try:
            with open('focal_minors.txt', 'r') as f:
                focal_minors = f.read().strip().split(',')
                focal_minors = [minor.strip() for minor in focal_minors]
        except FileNotFoundError:
            print('A file called focal_minors.txt must be included. This file should have a csv of target minor clades such as "Am_tu, Sr_ci, Sr_rh"')
            exit()
    
    focal_majors = []
    if use_focal_majors:
        try:
            with open('focal_majors.txt', 'r') as f:
                focal_majors = f.read().strip().split(',')
                focal_majors = [major.strip() for major in focal_majors]
        except FileNotFoundError:
            print('A file called focal_majors.txt must be included. This file should have a list of target major clades such as "Am, Sr"')
            exit()
    
    count_data = {}
    major_clades = set()
    minor_clades = set()
    
    for file in os.listdir(in_dir):
        if file.split('.')[-1] in ('fasta', 'fas', 'faa', 'fna'):
            fname = os.path.join(in_dir, file)
            
            count_data[file] = {}
            
            # Use SeqIO.index for faster access
            seq_index = SeqIO.index(fname, "fasta")
            tips = [seq_index[record_id].id[:10] for record_id in seq_index]
            
            for tip in tips:
                major_clade = tip[:2]
                minor_clade = tip[:5]
                
                if use_focal_minors and minor_clade not in focal_minors:
                    continue
                
                if use_focal_majors and major_clade not in focal_majors:
                    continue
                
                major_clades.add(major_clade)
                minor_clades.add(minor_clade)
                
                if tip not in count_data[file]:
                    count_data[file][tip] = 0
                count_data[file][tip] += 1
    
    if use_focal_minors:
        # Filter major and minor clades based on focal minors
        major_clades = sorted({minor[:2] for minor in focal_minors})
        minor_clades = sorted(focal_minors)
    elif use_focal_majors:
        # Filter major and minor clades based on focal majors
        major_clades = sorted(focal_majors)
        minor_clades = sorted({clade for clade in minor_clades if clade[:2] in focal_majors})
    else:
        major_clades = sorted(major_clades)
        minor_clades = sorted(minor_clades)
    
    taxa = sorted({tax for file_data in count_data.values() for tax in file_data})
    
    # Generate output file name based on input folder name
    folder_name = os.path.basename(in_dir)
    output_file_name = folder_name + "_TaxonOccurrence.csv"
    
    with open(output_file_name, 'w') as o:
        o.write(',' + ','.join(major_clades + minor_clades + taxa) + '\n')
        for file in count_data:
            o.write(file)
            
            major_clade_counts = {clade: 0 for clade in major_clades}
            minor_clade_counts = {clade: 0 for clade in minor_clades}
            
            for tax, count in count_data[file].items():
                major_clade_counts[tax[:2]] += count
                minor_clade_counts[tax[:5]] += count
            
            for clade in major_clades:
                o.write(',' + str(major_clade_counts[clade]))
                
            for clade in minor_clades:
                o.write(',' + str(minor_clade_counts[clade]))
            
            for tax in taxa:
                if tax in count_data[file]:
                    o.write(',' + str(count_data[file][tax]))
                else:
                    o.write(',0')
            o.write('\n')


def main():
    in_dir, use_focal_minors, use_focal_majors = get_args()
    count_tips(in_dir, use_focal_minors, use_focal_majors)
    

main()
