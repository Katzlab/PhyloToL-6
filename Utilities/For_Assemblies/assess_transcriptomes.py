'''
Written March 2023 by Elinor (esterner27@gmail.com) to plot length, coverage and GC of assembled transcripts

This script will rename the spades output to new names in the txt file, then iterate through them all and gather GC, length and coverage. With that data, it plots R scripts

Input: 
	Directory of directories output by rnaSpades OR folder called Renamed_assembled_files of previously renamed files (if this is the case, put -r or --renamed in the command line)
	txt file of LKH number and new names formatted like this: LKHxxx\tLKHxxx-10_digit_code-descriptor_of_taxon

To run if your files are already renamed:
	python assess_transcriptomes.py <pathway to directory of spades output>
	python assess_transcriptomes.py --renamed

Output: csv file of length, GC, coverage of each transcript, and multiple R plots, faceted by taxon and a csv file of data. It plots GC by length, and distributions of coverage, length and GC

'''
import os
import sys
from pathlib import Path
import sys
from Bio.SeqUtils import GC
from Bio import SeqIO

def script_help():

	print('\nThis script grabs and plots GC, length and coverage of transcriptomes. \n\nInput:\ntxt file of tab separated LKH number, ten digit code and taxon info (taxonomy, lifestage, etc).\nAND\nfolder of the folders output by spades named with LKH number \n(LKH999 or WTALKH999)\nOR\ndirectory of renamed assemblies in this format: ten_digit_code_assembledTranscripts.fasta\n\nOutput is multiple R plots, faceted by taxon and a csv file of data. \n\nIt plots GC by length, and distributions of coverage, length and GC.\n\n To run: \n\n\tpython assess_transcriptomes.py <pathway to directory of spades output>\n\nOptional parameters:\n\n-r or --renamed if your assemblies are already renamed to this format: ten_digit_code_assembledTranscript.fasta\n\n-h or --help for this message\n\n')

def get_args():
	#this parses user arguments. Checks if the files are renamed already or not (--renamed or --raw), and gets the directory of those files.
	
	renamed = False
	if('--help' in sys.argv or '-h' in sys.argv):#check for help function in command line
		script_help()
		exit()


	if ('--renamed'in sys.argv or '-r' in sys.argv):#check for renamed parameter
		renamed = True
		try:
			if('--renamed' in sys.argv):
				input_dir = sys.argv[sys.argv.index('--renamed') + 1]
			else:
				input_dir = sys.argv[sys.argv.index('-r') + 1]
		except IndexError:
			print('\nSomething went wrong went parsing the arguments. Did you input a directory of assemblies?\n')

	if ('--raw' in sys.argv or '-w' in sys.argv):#check for renamed parameter 
		renamed = False
		try:
			if('--raw' in sys.argv):
				input_dir = sys.argv[sys.argv.index('--raw') + 1]
			else:
				input_dir = sys.argv[sys.argv.index('-r') + 1]
		except IndexError:
			print('\nSomething went wrong went parsing the arguments. Did you input a directory of assemblies?\n')
	make_dirs(renamed, input_dir)

	
def make_dirs(renamed, input_dir):

	Path(f'plots').mkdir(parents=True, exist_ok=True)#makes output folder for r plots
	
	if renamed == True:
		assess_transcriptomes(input_dir)#skip renaming funtion if theyre already renamed
	elif renamed == False:
		Path(f'Renamed_assembled_files').mkdir(parents=True, exist_ok=True)#makes output folder for renamed fasta files
		#print('renaming')
		#print(input_dir)
		rename_assembled_transcriptomes(input_dir)#send to renaming function if they're not renamed yet.
	else:
		print('\nplease specify if your files are already renamed (--renamed) or if they are not yet renamed (--raw)\n')


def rename_assembled_transcriptomes(input_dir):#this function grabs and renames assembled transcripts files to the format required for this script and for phylotolv6 part1 

	names = { line.split('\t')[0].strip() : line.split('\t')[1].strip() for line in open('new_names.txt')}#make dictionary of original name: new name to replace
	for dir in os.listdir(input_dir):
		for old_name, new_name in names.items():
			if(old_name in dir):
				if(os.path.isfile(f'{input_dir}/{dir}/transcripts.fasta')):#go into spades output directory
					print(f'{old_name} is being renamed to {new_name}\n')
					os.system(f'cp {input_dir}/{dir}/transcripts.fasta Renamed_assembled_files/{new_name}_assembledTranscripts.fasta')#rename and copy to renamed_assembled files
	input_dir = 'Renamed_assembled_files'

	assess_transcriptomes(input_dir)

def assess_transcriptomes(input_dir):
	
	#iterate through renamed fasta files
	for file in os.listdir(input_dir):
		if file.endswith('fasta'):
			print(f'gathering data from {file}..\n')
			data = []#initiate list of data per LKH
			records = list(SeqIO.parse(f'{input_dir}/{file}', 'fasta'))#parse the fasta files
			
			taxon_info = get_taxon_info(file)#send to get_taxon_info function
			#parse output from that function
			lkh = taxon_info[0]
			ten_digit_code = taxon_info[1]
			taxon = taxon_info[2]

		#	extract all data for each transcript
			for record in records:
				gc = GC(record.seq)
				length = len(record.seq)
				iden = record.id
				cov = record.id.split('_')[5]
		
				transcript_data = f'{iden}, {lkh}, {length}, {gc}, {cov}, {ten_digit_code}, {taxon}'#list of data per transcriptome
				data.append(transcript_data)#data for each LKH
		
			all_data.update({file.split('_assembledTranscripts')[0]: data})#make dict of ten digit code: all information

	write_to_file()

def get_taxon_info(file):#parse the info in the new_names.txt file to get info on the taxon
	with open('new_names.txt', 'r') as o:
		cell_info = o.readlines()
		for line in cell_info:
			if '_'.join(file.split('_')[0:3]) in line:
				lkh = line.split('\t')[0]
				ten_digit_code = line.split('\t')[1]
				taxon = line.split('\t')[2]
	
	try:
		return lkh, ten_digit_code, taxon
	except UnboundLocalError:
		print(f'no taxon information in new_names.txt for {file}')

	print('done with getting taxon info')


def write_to_file():

	with open('assembly_assessment.csv', 'w') as o:#create output csv
		print('writing data to assembly_assessment.csv\n')
		o.write('seqID, lkh, length, GC, cov, ten_digit_code, taxon_info\n')#write header
		for lkh, data in all_data.items():
			for transcript in data:
				o.write(f'{transcript}')#write each line

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
	get_args()
