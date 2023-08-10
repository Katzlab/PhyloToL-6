#Author, date: Auden Cote-L'Heureux, Aug 10 2023
#Motivation: Understand variation at nucleotide level between sequences
#Intent: Create a nucleotide alignment from an amino acid alignment
#Dependencies: Python3, Biopython
#Inputs: An aligned amino acid fasta file or folder of aligned amino acid fasta files, and a nucleotide fasta file or a folder of nucleotide fasta files
#Outputs: An aligned nucleotide file or folder of aligned nucleotide files
#Example: python3 BacktranslateAlignment.py -a AminoAcidAlignment.fasta -n Nucleotides.fasta


import os, sys
from Bio import SeqIO
import argparse

universal_6fold = {
        	'GCT': ['A', 'four', 0], 'GCC': ['A', 'four', 0], 'GCA': ['A', 'four', 0],
        	'GCG': ['A', 'four', 0], 'CGT': ['R', 'six', 0], 'CGC': ['R', 'six', 0],
        	'CGG': ['R', 'six', 0], 'CGA': ['R', 'six', 0], 'AGA': ['R', 'six', 0],
        	'AGG': ['R', 'six', 0], 'AAT': ['N', 'two', 0], 'AAC': ['N', 'two', 0],
        	'GAT': ['D', 'two', 0], 'GAC': ['D', 'two', 0], 'TGT': ['C', 'two', 0],
        	'TGC': ['C', 'two', 0], 'CAA': ['Q', 'two', 0], 'CAG': ['Q', 'two', 0],
        	'GAA': ['E', 'two', 0], 'GAG': ['E', 'two', 0], 'GGT': ['G', 'four', 0],
        	'GGC': ['G', 'four', 0], 'GGA': ['G', 'four', 0], 'GGG': ['G', 'four', 0],
        	'CAT': ['H', 'two', 0], 'CAC': ['H', 'two', 0], 'ATT': ['I', 'three', 0],
        	'ATC': ['I', 'three', 0], 'ATA': ['I', 'three', 0], 'ATG': ['M', 'one', 0],
        	'TTA': ['L', 'six', 0], 'TTG': ['L', 'six', 0], 'CTT': ['L', 'six', 0],
        	'CTC': ['L', 'six', 0], 'CTA': ['L', 'six', 0], 'CTG': ['L', 'six', 0],
        	'AAA': ['K', 'two', 0], 'AAG': ['K', 'two', 0], 'TTT': ['F', 'two', 0],
        	'TTC': ['F', 'two', 0], 'CCT': ['P', 'four', 0], 'CCC': ['P', 'four', 0],
        	'CCA': ['P', 'four', 0], 'CCG': ['P', 'four', 0], 'TCT': ['S', 'six', 0],
        	'TCC': ['S', 'six', 0], 'TCA': ['S', 'six', 0], 'TCG': ['S', 'six', 0],
        	'AGT': ['S', 'six', 0], 'AGC': ['S', 'six', 0], 'ACT': ['T', 'four', 0],
        	'ACC': ['T', 'four', 0], 'ACA': ['T', 'four', 0], 'ACG': ['T', 'four', 0],
        	'TGG': ['W', 'one', 0], 'TAT': ['Y', 'two', 0], 'TAC': ['Y', 'two', 0],
        	'GTT': ['V', 'four', 0], 'GTC': ['V', 'four', 0], 'GTA': ['V', 'four', 0],
        	'GTG': ['V', 'four', 0], 'TAA': ['*', 'none', 0], 'TGA': ['*', 'none', 0],
        	'TAG': ['*', 'none', 0]}

aas = list(dict.fromkeys([universal_6fold[codon][0] for codon in universal_6fold]))
codons_per_aa = { aa : [codon for codon in universal_6fold if universal_6fold[codon][0] == aa] for aa in aas }


def get_args():

	parser = argparse.ArgumentParser(
		prog = 'Backtranslating script, Version 1.0',
		description = "Updated Aug 10th, 2023 by Auden Cote-L'Heureux"
	)

	parser.add_argument('-a', '--amino', type = str, required = True, help = 'Path to a fasta file or folder containing multiple fasta files with aligned amino acid sequences')
	parser.add_argument('-n', '--nucl', type = str, required = True, help = 'Path to a fasta file or folder containing multiple fasta files with nucleotide sequences. Every amino acid sequence considered must have a nucleotide sequence in this (these) file(s) of matching length with the same sequence identifier.')

	return parser.parse_args()


def backtranslate(aa, nucls, output):

	o = open(output, 'w')

	for rec in aa:
		if(rec.id in nucls):
			if len(nucls[rec.id]) == len(str(rec.seq).replace('-', '')) * 3:
				running_seq = ''; c = 0; fail = False
				for i, char in enumerate(str(rec.seq)):
					if(char == '-'):
						running_seq += '---'
					else:
						codon = nucls[rec.id][c:c+3]
						if char == 'X' or codon in codons_per_aa[char]:
							running_seq += codon
							c += 3
						else:
							fail = True

				if fail:
					print('\nWARNING: The nucleotide sequence ' + rec.id + ' does not match the corresponding amino acid sequence. This sequence will be missing from the alignment.\n')
				else:
					o.write('>' + rec.id + '\n')
					o.write(running_seq + '\n\n')
			else:
				print('\nWARNING: The nucleotide sequence ' + rec.id + ' is not 3x the length of the corresponding amino acid sequence. This sequence will be missing from the alignment.\n')
		else:
			print('\nWARNING: There is no nucleotide sequence for the amino acid sequence ' + rec.id + '. This sequence will be missing from the alignment.\n')

	o.close()


if __name__ == '__main__':
	
	args = get_args()

	if os.path.isfile(args.nucl):
		if args.nucl.split('.')[-1] in ('fasta', 'fas', 'fna'):
			try:
				nucls = { rec.id : str(rec.seq).replace('-', '') for rec in SeqIO.parse(args.nucl, 'fasta') }
			except:
				print('\nERROR: It appears that a single file of nucleotide sequences was input but is improperly formatted. Make sure this file has the extension fasta, fas, or fna and contains unaligned nucleotide sequences.\n')
				exit()
		else:
			print('\nERROR: It appears that a single file of nucleotide sequences was input but is improperly formatted. Make sure this file has the extension fasta, fas, or fna and contains unaligned nucleotide sequences.\n')
			exit()
	elif os.path.isdir(args.nucl):
		try:
			nucls = { rec.id : str(rec.seq).replace('-', '') for file in os.listdir(args.nucl) for rec in SeqIO.parse(args.nucl + '/' + file, 'fasta') if file.split('.')[-1] in ('fasta', 'fas', 'fna') }
		except:
			print('\nERROR: It appears that a folder of nucleotide files was input but one or more files is improperly formatted. Make sure the files have the extension fasta, fas, or fna and contain unaligned nucleotide sequences.\n')
			exit()
	else:
		print('\nERROR: Invalid path to a nucleotide file or folder of nucleotide files.\n')
		exit()

	if len(nucls) > 0:
		if os.path.isfile(args.amino):
			if args.amino.split('.')[-1] in ('fasta', 'fas', 'faa'):
				try:
					aa = [rec for rec in SeqIO.parse(args.amino, 'fasta')]
				except:
					print('\nERROR: It appears that a single amino acid file was input but is improperly formatted. Make sure the file has the extension fasta, fas, or faa and contains aligned amino acid sequences.\n')
					exit()

				backtranslate(aa, nucls, '.'.join(args.amino.split('.')[:-1]) + '_NTDAligned.fasta')
			else:
				print('\nERROR: It appears that a single amino acid file was input but is improperly formatted. Make sure the file has the extension fasta, fas, or faa and contains aligned amino acid sequences.\n')
				exit()
		elif os.path.isdir(args.amino):

			if not os.path.isdir('NTDAligned'):
				os.mkdir('NTDAligned')

			for file in os.listdir(args.amino):
				if file.split('.')[-1] in ('fasta', 'fas', 'faa'):
					try:
						aa = [rec for rec in SeqIO.parse(args.amino + '/' + file, 'fasta')]
					except:
						print('\nERROR: the amino acid file ' + file + ' could not be read. Make sure the file has the extension fasta, fas, or faa and contains aligned amino acid sequences.\n')
						exit()

					backtranslate(aa, nucls, 'NTDAligned/' + '.'.join(file.split('.')[:-1]) + '_NTDAligned.fasta')
		else:
			print('\nERROR: Invalid path to a amino acid file or folder of amino acid files.\n')
			exit()
	else:
		print('\nERROR: No nucleotide sequences were read from the input file(s). Make sure these files are properly formatted and not empty.\n')
		exit()


















