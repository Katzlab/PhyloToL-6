import os, sys, re
from Bio import SeqIO


def check_cds(params):

	if os.path.isdir(params.cds):
		for file in os.listdir(params.cds):
			if file[10:] != '_GenBankCDS.fasta' and 'DS_Store' not in file:
				print('\nERROR: The file ' + file + ' in the give folder of assembled transcripts is incorrectly formatted. The files must start with a ten digit taxon identifier and then be named like Op_me_Hsap_GenBankCDS.fasta\n')
				exit()
	else:
		print('\nERROR: CDS folder could not be found. Please ensure the given path is correct.\n')
		exit()


def check_databases(params):

	if os.path.isdir(params.databases):
		if os.path.isdir(params.databases + '/db_OG'):
			fasta = [file for file in os.listdir(params.databases + '/db_OG') if file.endswith('.fasta')]
			dmnd = [file for file in os.listdir(params.databases + '/db_OG') if file.endswith('.dmnd')]

			if len(fasta) == 0:
				print('\nERROR: No Hook fasta file found in the Databases/db_OG folder\n')
				exit()
			elif len(fasta) > 1:
				print('\nERROR: More than one Hook fasta file found in the Databases/db_OG folder. Please delete all except for the correct file.\n')
				exit()
			else:
				for rec in SeqIO.parse(params.databases + '/db_OG/' + fasta[0], 'fasta'):
					try:
						og_number = re.split('OG.{1}_', rec.id)[-1][:6]
						og_prefix = rec.id.split(og_number)[-2][-4:]
						og = og_prefix + og_number
						
						if rec.id[-10:] != og:
							print('\nError: The sequence name ' + rec.id + ' in the given Hook database fasta file is incorrectly formatted. Each sequence ID should start with a ten-digit taxon identifier and end with a ten-digit gene family identifier (which must start with OGX_, with "X" being any digit. E.g. Op_me_Hsap_0_OG6_110767)\n')
							exit()
					except IndexError:
						print('\nError: The sequence name ' + rec.id + ' in the given Hook database fasta file is incorrectly formatted. Each sequence ID should start with a ten-digit taxon identifier and end with a ten-digit gene family identifier (which must start with OGX_, with "X" being any digit. E.g. Op_me_Hsap_0_OG6_110767)\n')
						exit()
			if len(dmnd) == 0:
				print('\nERROR: No Hook Diamond database (.dmnd) file found in the Databases/db_OG folder.\n')
				exit()
			elif len(dmnd) > 1:
				print('\nERROR: No Hook Diamond database (.dmnd) file found in the Databases/db_OG folder. Please delete all except for the correct file.\n')
				exit()
		else:
			print('\nERROR: The db_OG folder could not be found in the databases folder.\n')
			exit()
	else:
		print('\nERROR: Databases folder could not be found. Please ensure the given path is correct.\n')
		exit()


def run(params):

	print('\nChecking the input files and scripts setup...\n')

	check_cds(params)

	check_databases(params)

	print('\nAll checks passed!\n')






