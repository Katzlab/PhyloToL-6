'''
Written March 2023 by Elinor (esterner27@gmail.com) to plot length, coverage and GC of assembled transcripts

This script will iterate through all assembled files (named as 10 digit code _assembledTranscripts) with and gather GC, length and coverage. With that data, it plots R scripts

Input: 
	Folder called Renamed_assembled_files of previously renamed files (if this is the case, put -r or --renamed in the command line)
	tsv file of LKH number and new names formatted like this: LKHxxx\tLKHxxx-10_digit_code\tdescriptor_of_taxon called new_names.tsv

To run:
	python assess_transcriptomes.py -input <pathway to directory of spades output>

Output: csv file of length, GC, coverage of each transcript, and multiple R plots, faceted by taxon and a csv file of data. It plots GC by length, and distributions of coverage, length and GC

'''
import os
import sys
from pathlib import Path
import sys
from Bio.SeqUtils import GC
from Bio import SeqIO

def script_help():

	print('\nThis script grabs and plots GC, length and coverage of transcriptomes. \n\nInput:\ntsv file of tab separated LKH number, ten digit code and taxon info (taxonomy, lifestage, etc).\nAND\nfolder of the folders output by spades named with LKH number \n(LKH999 or WTALKH999)\nOR\ndirectory of renamed assemblies in this format: ten_digit_code_assembledTranscripts.fasta\n\nOutput is multiple R plots, faceted by taxon and a csv file of data. \n\nIt plots GC by length, and distributions of coverage, length and GC.\n\n To run: \n\n\tpython assess_transcriptomes.py <pathway to directory of spades output>\n\n-r or --renamed if your assemblies are already renamed to this format: ten_digit_code_assembledTranscript.fasta/nand this command if they are not yet named: --raw <path to directory of spades output folders>\n\n-h or --help for this message\n\n')

def get_args():
	#this parses user arguments. Checks if the files are renamed already or not (--renamed or --raw), and gets the directory of those files.
	
	if('--help' in sys.argv or '-h' in sys.argv):#check for help function in command line
		script_help()
		exit()


	if ('--input'in sys.argv or '-i' in sys.argv):#check for renamed parameter
		renamed = True
		try:
			if('--input' in sys.argv):
				input_dir = sys.argv[sys.argv.index('--input') + 1]
			else:
				input_dir = sys.argv[sys.argv.index('-i') + 1]
		except IndexError:
			print('\nSomething went wrong went parsing the arguments. Did you input a directory of assemblies?\n')

	
		make_dirs(input_dir)

	
def make_dirs(input_dir):

	Path(f'plots').mkdir(parents=True, exist_ok=True)#makes output folder for r plots
	assess_transcriptomes(input_dir)#skip renaming funtion if theyre already renamed
	
	

def assess_transcriptomes(input_dir):
	
	#iterate through renamed fasta files
	for file in os.listdir(input_dir):
		if file.endswith('fasta'):
			print(f'gathering data from {file}..\n')
			data = []#initiate list of data per LKH
			records = list(SeqIO.parse(f'{input_dir}/{file}', 'fasta'))#parse the fasta files
			
			taxon_info = get_taxon_info(file)#send to get_taxon_info function
			#parse output from that function
			ten_digit_code = taxon_info[0]
			taxon = taxon_info[1]

		#	extract all data for each transcript
			for record in records:
				gc = GC(record.seq)
				length = len(record.seq)
				iden = record.id
				cov = record.id.split('_')[5]
		
				transcript_data = f'{iden}, {length}, {gc}, {cov}, {ten_digit_code}, {taxon}'#list of data per transcriptome
				data.append(transcript_data)#data for each LKH
		
			all_data.update({file.split('_assembledTranscripts')[0]: data})#make dict of ten digit code: all information

	write_to_file()

def get_taxon_info(file):#parse the info in the new_names.tsv file to get info on the taxon
	with open('new_names.tsv', 'r') as o:
		cell_info = o.readlines()
		for line in cell_info:
			if '_'.join(file.split('_')[0:3]) in line:
				ten_digit_code = line.split('\t')[0]
				taxon = line.split('\t')[1]
	
	try:
		return ten_digit_code, taxon
	except UnboundLocalError:
		print(f'no taxon information in new_names.tsv for {file}')

	print('done getting taxon info from new_names.tsv')


def write_to_file():

	with open('assembly_assessment.csv', 'w') as o:#create output csv
		print('writing data to assembly_assessment.csv\n')
		o.write('seqID, length, GC, cov, ten_digit_code, taxon_info\n')#write header
		for lkh, data in all_data.items():
			for transcript in data:
				o.write(f'{transcript.strip()}\n')#write each line

	plot_assessment()

def plot_assessment():
	print('Plotting your data\n')
	os.system('Rscript plot_assemblies.R')
	for file in os.listdir():
		if file.endswith('.png') or file.endswith('.pdf'):
			os.system(f'mv {file} plots/{file}')

	print('\n\n\n\nDone! All your plots are saved to the folder called plots\n\n')

if __name__ == '__main__':
	all_data = {}
	
	try:
		f = open("new_names.tsv")

	except FileNotFoundError:
		print('\n\n Did you include a tsv file of your cells? It should be called new_names.tsv and formatted like this: 10_digit_code\tdescriptor_of_taxon\n\n')	
		exit()

	else:
		get_args()
