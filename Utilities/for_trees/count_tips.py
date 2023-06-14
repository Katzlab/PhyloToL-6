'''
Author: Auden Cote-L'Heureux
Last updated: 06/14/23 by Elinor
Motivation: Count the number of occurences of each taxa in each OG in a post guidance file
Dependencies: Bio python, os, sys
Inputs: Directory of postguidance files
Outputs: CSV file tallying all the counts of taxa in each OG file
Command line: python count_tips.py --input <dir of postguidance files>
'''

import os
import sys
from Bio import SeqIO


def get_args():

	in_dir = ''
	
	if('--input' in sys.argv or '-i' in sys.argv):
		try:
			if('--input' in sys.argv):
				in_dir = sys.argv[sys.argv.index('--input') + 1]
			else:
				in_dir = sys.argv[sys.argv.index('-i') + 1]
		except IndexError:
			print('\nError: Something went wrong went parsing the arguments... maybe you forgot to input an input directory of trees?\n')
			print('\nPlease input a folder of .tre files:\n\n\tpython count_tips.py --input <path/to/folder>\n')
			exit()
	else:
		print('\nPlease input a folder of .tre files:\n\n\tpython count_tips.py --input <path/to/folder>\n')
		exit()
		
		
	if(in_dir.endswith('/')):
		in_dir = in_dir[:-1]
		
	if(not os.path.isdir(in_dir)):
		print('\nPlease input a folder of .tre files:\n\n\tpython count_tips.py --input <path/to/folder>\n')
		exit()
				
	return in_dir


def count_tips(in_dir):

	count_data = { }
	for file in os.listdir(in_dir):
		if(file.endswith('.fas')):
			fname = in_dir + '/' + file
			
			count_data.update({ file : { } })
			tips = [record.id[:10] for record in SeqIO.parse(in_dir+'/'+file, 'fasta')]
					
			for tip in tips:
				tip = tip.strip()
				if(tip[:10] not in count_data[file]):
					count_data[file].update({ tip[:10] : 0 })
				count_data[file][tip[:10]] += 1
				
	taxa = sorted(list(dict.fromkeys([tax for file in count_data for tax in count_data[file]])))
				
	with open('tip_count_data.csv', 'w') as o:
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