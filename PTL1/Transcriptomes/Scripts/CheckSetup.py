import os, sys, re
from Bio import SeqIO


def check_transcripts(params):

	taxa = []

	if os.path.isdir(params.assembled_transcripts):
		for file in os.listdir(params.assembled_transcripts):
			if file[10:] == '_assembledTranscripts.fasta':
				taxa.append(file[:10])

				for rec in SeqIO.parse(params.assembled_transcripts + '/' + file, 'fasta'):
					if 'node_' not in rec.id.lower() or '_length_' not in rec.id.lower() or '_cov_' not in rec.id.lower():
						print('\nERROR: The sequence record ' + rec.id + ' from taxon ' + file[:10] + ' is incorrectly formatted. All sequence names must have be formatted in the style of rnaSpades output (with a NODE ID, length, and coverage)\n')
						exit()
			elif 'DS_Store' not in file:
				print('\nERROR: The file ' + file + ' in the give folder of assembled transcripts is incorrectly formatted. The files must start with a ten digit taxon identifier and then be named like Am_tu_Hp01_assembledTranscripts.fasta\n')
				exit()
	else:
		print('\nERROR: Assembled transcripts folder could not be found. Please ensure the given path is correct.\n')
		exit()


	return taxa


def check_databases(params):

	if os.path.isdir(params.databases):
		if os.path.isdir(params.databases + '/db_BvsE'):
			if not os.path.isfile(params.databases + '/db_BvsE/eukout.dmnd'):
				print('\nERROR: eukout.dmnd could not be found in the Databses/db_BvsE folder')
				exit()
			if not os.path.isfile(params.databases + '/db_BvsE/micout.dmnd'):
				print('\nERROR: micout.dmnd could not be found in the Databses/db_BvsE folder')
				exit()
			if not os.path.isfile(params.databases + '/db_BvsE/SSULSUdb.nhr'):
				print('\nERROR: SSULSUdb.nhr could not be found in the Databses/db_BvsE folder')
				exit()
			if not os.path.isfile(params.databases + '/db_BvsE/SSULSUdb.nin'):
				print('\nERROR: SSULSUdb.nin could not be found in the Databses/db_BvsE folder')
				exit()
			if not os.path.isfile(params.databases + '/db_BvsE/SSULSUdb.nsq'):
				print('\nERROR: SSULSUdb.nsq could not be found in the Databses/db_BvsE folder')
				exit()
		else:
			print('\nERROR: The db_BvsE folder could not be found in the databases folder.\n')
			exit()
		
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
		
		if os.path.isdir(params.databases + '/db_StopFreq'):
			if not os.path.isfile(params.databases + '/db_StopFreq/RepEukProts.dmnd'):
				print('\nERROR: The RepEukProts.dmnd file could not be found in the Databases/db_StopFreq folder.\n')
		else:
			print('\nERROR: The db_StopFreq folder could not be found in the databases folder.\n')
			exit()
	else:
		print('\nERROR: Databases folder could not be found. Please ensure the given path is correct.\n')
		exit()


def check_conspecifics(params, transcript_taxa):

	if params.xplate_contam:
		if params.conspecific_names != None:
			if os.path.isfile(params.conspecific_names):
				lines = [line for line in open(params.conspecific_names)]

				seen_taxa = []
				for line in lines:
					if len(line.split('\t')) > 0 and not len(line.split('\t')) > 2:
						taxon = line.split('\t')[0]
						if taxon not in seen_taxa:
							seen_taxa.append(taxon)
						else:
							print('\nERROR: The taxon ' + taxon + ' is listed in the conspecific names file twice. Please make sure this file is correct.\n')
							exit()
						
						if taxon not in transcript_taxa:
							print('\nERROR: The taxon ' + taxon + ' is listed in the conspecific names file but there is no assembled transcripts file that corresponds to this taxon.\n')
							exit()
					else:
						print('\nERROR: Unclear how to parse the line ' + line + ' in the given conspecific names file. This file should have two tab-separated columns.\n')
						exit()

				missing = [tax for tax in transcript_taxa if tax not in seen_taxa]
				if len(missing) > 0:
					print('\nERROR: The following taxa have assembled transcripts but are not listed in the conspecific names file:' + '\n\t'.join(missing) + '\n')
					exit()
			else:
				print('\nERROR: The given file of conspecific names (' + params.conspecific_names + ') could not be found.\n')
				exit()
		else:
			print('\nERROR: If running cross-plate contamination, a file with conspecific names is required.\n')
			exit()


def check_gcodes(params, transcript_taxa):

	valid_codes = ['bleph','blepharisma','chilo','chilodonella','condy', 'condylostoma','none','eup','euplotes','peritrich','vorticella','ciliate','universal','taa','tag','tga','mesodinium']

	if params.genetic_code != None:
		if os.path.isfile(params.genetic_code):
			lines = [line for line in open(params.genetic_code)]

			seen_taxa = []
			for line in lines:
				if len(line.split('\t')) > 0 and not len(line.split('\t')) > 2:
					taxon = line.split('\t')[0]

					if taxon not in seen_taxa:
						seen_taxa.append(taxon)
					else:
						print('\nERROR: The taxon ' + taxon + ' is listed in the genetic codes file twice. Please make sure this file is correct.\n')
						exit()
					
					if taxon not in transcript_taxa:
						print('\nERROR: The taxon ' + taxon + ' is listed in the genetic codes file but there is no assembled transcripts file that corresponds to this taxon.\n')
						exit()

					if line.split('\t')[1].strip().lower() not in valid_codes:
						print('\nERROR: The code ' + line.split('\t')[1].strip() + ' is not a valid genetic code. Make sure you input only accepted genetic codes and that the genetic codes file is properly formatted.\n')
						exit()
				else:
					print('\nERROR: Unclear how to parse the line ' + line + ' in the given genetic codes file. This file should have two tab-separated columns.\n')
					exit()

			missing = [tax for tax in transcript_taxa if tax not in seen_taxa]
			if len(missing) > 0:
				print('\nERROR: The following taxa have assembled transcripts but are not listed in the genetic codes names file:' + '\n\t'.join(missing) + '\n')
				exit()
		else:
			if params.genetic_code.lower() not in valid_codes:
				print('\nERROR: The code ' + params.genetic_code + ' is not a valid genetic code. Make sure you input only accepted genetic codes and that the genetic codes file is properly formatted.\n')
				exit()

	else:
		print('\nERROR: An input genetic code is required.\n')
		exit()


def run(params):

	print('\nChecking the input files and scripts setup...\n')

	transcript_taxa = check_transcripts(params)

	check_databases(params)

	check_conspecifics(params, transcript_taxa)

	check_gcodes(params, transcript_taxa)

	print('\nAll checks passed!\n')






