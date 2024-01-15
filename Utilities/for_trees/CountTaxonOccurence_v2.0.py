'''
Author: Auden Cote-L'Heureux
Last updated: 10/24/23
Motivation: Count the number of occurences of each taxa in each OG in a post guidance file
Dependencies: Bio python, os, sys
Inputs: Directory of postguidance files
Outputs: CSV file tallying all the counts of taxa in each OG file
Command line: python count_tips.py --input <dir of postguidance files>
'''

import os
import sys
from Bio import SeqIO
import argparse


def get_args():

	parser = argparse.ArgumentParser(
		prog = 'Taxon occurence counting script',
		description = "Updated Oct 24th, 2023 by Auden Cote-L'Heureux."
	)

	parser.add_argument('-i', '--input', type = str, required = True, help = 'Path to the folder containing the aligned/unaligned fasta files')
	args = parser.parse_args()
		
	if(args.input.endswith('/')):
		args.input = args.input[:-1]
		
	if(not os.path.isdir(args.input)):
		print('\nThe input folder (--input) could not be found. Make sure you have given the correct path.\n')
		exit()
				
	return args.input


def count_tips(in_dir):

	count_data = { }
	for file in os.listdir(in_dir):
		if file.split('.')[-1] in ('fasta', 'fas', 'faa', 'fna'):
			fname = in_dir + '/' + file
			
			count_data.update({ file : { } })
			tips = [record.id[:10] for record in SeqIO.parse(in_dir+'/'+file, 'fasta')]
					
			for tip in tips:
				tip = tip.strip()
				if(tip[:10] not in count_data[file]):
					count_data[file].update({ tip[:10] : 0 })
				count_data[file][tip[:10]] += 1
				
	taxa = sorted(list(dict.fromkeys([tax for file in count_data for tax in count_data[file]])))
				
	with open('TaxonOccurrence.csv', 'w') as o:
		o.write(',' + ','.join(taxa) + '\n')
		for file in count_data:
			o.write(file)
			for tax in taxa:
				if(tax in count_data[file]):
					o.write(',' + str(count_data[file][tax]))
				else:
					o.write(',0')
			o.write('\n')	
	
	
def main():

	in_dir = get_args()

	count_tips(in_dir)
	

main()
