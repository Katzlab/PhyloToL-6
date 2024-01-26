# Last updated Sept 2023
# Author: Auden Cote-L'Heureux

# This script is a WRAPPER for the PhyloToL Part 1 GENOMES pipeline. Users should
# use this script to run the pipeline, rather than running any of the sub-scripts (number 1a through 5b)
# independently. To run an individual step in the pipeline, use --script X where X is the number (1 through 5).
# To run multiple sets (usually all of them), use --first script 1 --last_script 5, or whichever first
# and last scripts are desired. Run "python wrapper.py --help" for details on how to run this script. Before
# running this script ensure that the databases are correctly located and named, and that input CDS are named 
# in the format Op_me_Hsap_GenBankCDS.fasta, where Op_me_Hsap can be replaced with any 10-digit sample 
# identifier.


import os, sys, re
import argparse
import CheckSetup


def get_args():

	parser = argparse.ArgumentParser(
                prog = 'PhyloToL v6.0 Part 1 for GenBank Genomes',
                description = "Updated January 19th, 2023 by Auden Cote-L'Heureux. Link to GitHub: https://github.com/AudenCote/PhyloToL_v6.0"
                )

	parser.add_argument('-s', '--script', default = -1, type = int, choices = { 1, 2, 3, 4, 5 }, help = 'Script to run if you are only running one script')
	parser.add_argument('-1', '--first_script', default = -1, type = int, choices = { 1, 2, 3, 4 }, help = 'First script to run')
	parser.add_argument('-2', '--last_script', default = -1, type = int, choices = { 2, 3, 4, 5 }, help = 'First script to run')
	parser.add_argument('-c', '--cds', type = str, help = 'Path to a folder of nucleotide CDS. Each file name should start with a unique 10 digit code, and end in "_GenBankCDS.fasta", E.g. Op_me_hsap_GenBankCDS.fasta')
	parser.add_argument('-o', '--output', default = '../', type = str, help = 'An "Output" folder will be created at this directory to contain all output files. By default this folder will be created at the parent directory of the Scripts folder')
	parser.add_argument('-g', '--genetic_code', type = str, help = 'If all of your taxa use the same genetic code, you may enter it here (to be used in script 4). Otherwise, stop after script 3 and fill in "gcode_output.tsv" before running script 4')
	parser.add_argument('-d', '--databases', type = str, default = '../Databases', help = 'Path to databases folder (which should contain db_OG)')

	return parser.parse_args()


def script_one(args, ten_digit_codes):

	CheckSetup.run(args)

	for file in os.listdir(args.cds):
		if file[10:] == '_GenBankCDS.fasta' and file[:10] in ten_digit_codes:
			os.system('python 1_RenameCDS.py -in ' + args.cds + '/' + file + ' -s GenBank -o ' + args.output + '/Output')


def script_two(args):
	
	valid_codes = ['bleph','blepharisma','chilo','chilodonella','condy', 'condylostoma','none','eup','euplotes','peritrich','vorticella','ciliate','universal','taa','tag','tga','mesodinium']

	for folder in os.listdir(args.output + '/Output'):
		if os.path.isfile(args.output + '/Output/' + folder + '/' + folder + '_GenBankCDS.Prepped.fasta'):
			os.system('python 2_GCodeEval.py --input_file ' + args.output + '/Output/' + folder + '/' + folder + '_GenBankCDS.Prepped.fasta')

	gcode_info = []
	for folder in os.listdir(args.output + '/Output'):
		if os.path.isfile(args.output + '/Output/' + folder + '/' + folder + '_GenBankCDS.Prepped.GeneticCode.txt'):
				with open(args.output + '/Output/' + folder + '/' + folder + '_GenBankCDS.Prepped.GeneticCode.txt') as f:
					gcode_temp = [folder]
					for line in f:
						line_sep = line.strip().split('\t')
						if line_sep[0] == 'TGA':
							gcode_temp.append(line_sep[1])
						elif line_sep[0] == 'TAG':
							gcode_temp.append(line_sep[1])
						elif line_sep[0] == 'TAA':
							gcode_temp.append(line_sep[1])

					gcode_info.append(gcode_temp)

	stop = False
	gcode_file = { }
	if args.genetic_code.endswith('.txt') or args.genetic_code.endswith('.tsv'):
		if os.path.isfile(args.genetic_code):
			for line in open(args.genetic_code):
				try:
					if line.split('\t')[1].strip().lower() in valid_codes:
						gcode_file.update({ line.split('\t')[0] : line.split('\t')[1].strip() })
					else:
						print('Genetic code ERROR -- ' + line.split('\t')[1].strip() + ' is not a valid genetic code. Please fill out the "gcode_output.tsv" file and continue with script 3.')
				except IndexError:
					print('\nGenetic code ERROR -- it looks like you tried to enter a .txt/.tsv file, but it is improperly formatted. Stopping after script 2; you may fill out the file gcode_output.tsv and continue with script 3.\n')
					stop = True
		else:
			print('\nGenetic code ERROR -- it looks like you tried to enter a .txt/.tsv file, but it could not be found. Stopping after script 2; you may fill out the file gcode_output.tsv and continue with script 3.\n')
			stop = True

	with open(args.output + '/Output/gcode_output.tsv', 'w') as g:
		g.writelines('10 Digit Code\tIn-frame TAG Density\tIn-frame TGA Density\tIn-frame TAA Density\tGenetic Code\n')
		for row in gcode_info:
			if args.genetic_code == None:
				g.writelines(row[0] + '\t' + row[1] + '\t' + row[2] + '\t' + row[3] + '\t\n')
			elif args.genetic_code.lower() in valid_codes:
				g.writelines(row[0] + '\t' + row[1] + '\t' + row[2] + '\t' + row[3] + '\t' + args.genetic_code + '\n')
			elif args.genetic_code.endswith('.txt') or args.genetic_code.endswith('.tsv'):
				try:
					g.writelines(row[0] + '\t' + row[1] + '\t' + row[2] + '\t' + row[3] + '\t' + gcode_file[row[0]] + '\n')
				except KeyError:
					g.writelines(row[0] + '\t' + row[1] + '\t' + row[2] + '\t' + row[3] + '\t\n')
					print('\nGenetic code ERROR -- it looks like you tried to enter a .txt/.tsv file, but a taxon is missing. Stopping after script 2; you may fill out the file gcode_output.tsv and continue with script 3.\n')
					stop = True
			else:	
				stop = True

	if stop or args.genetic_code == None:
		print('\nStopping after script 2 because genetic code information is incomplete. Please fill out the file "gcode_output.tsv" and continue with script 3.\n')
		exit()


def script_three(args):

	valid_codes = ['bleph','blepharisma','chilo','chilodonella','condy', 'condylostoma','none','eup','euplotes','peritrich','vorticella','ciliate','universal','taa','tag','tga','mesodinium']
	
	lines = [line.strip().split('\t') for line in open(args.output + '/Output/gcode_output.tsv', 'r')]
	with open(args.output + '/Output/gcode_output.tsv', 'r') as g:
		for folder in os.listdir(args.output + '/Output'):
			if os.path.isfile(args.output + '/Output/' + folder + '/' + folder + '_GenBankCDS.Prepped.fasta'):
				for line in lines:
					if line[0] == folder and line[-1].lower() in valid_codes:
						os.system('python 3_GCodeTranslate.py --input_file ' + args.output + '/Output/' + folder + '/' + folder + '_GenBankCDS.Prepped.fasta --genetic_code ' + line[-1])
					elif line[-1].lower() not in valid_codes and line[-1] != 'Genetic Code':
						print('\n' + line[-1] + ' is not a valid genetic code. Skipping taxon ' + folder + '.\n')


def script_four(args):

	valid_codes = ['universal', 'blepharisma', 'chilodonella', 'condylostoma', 'euplotes', 'peritrich', 'vorticella', 'mesodinium', 'tag', 'tga', 'taa', 'none']
	
	gcode_by_folder = { line.strip().split('\t')[0] : line.strip().split('\t')[-1] for line in open(args.output + '/Output/gcode_output.tsv', 'r') }
	for folder in os.listdir(args.output + '/Output'):
		if os.path.isdir(args.output + '/Output/' + folder):
			gcode_formatted = gcode_by_folder[folder][0].upper() + gcode_by_folder[folder].lower()[1:]
			if os.path.isfile(args.output + '/Output/' + folder + '/' + folder + '_GenBankCDS.' + gcode_formatted + '.AA.fasta'):
				os.system('python 4_CountOGsDiamond.py -in ' + args.output + '/Output/' + folder + '/' + folder + '_GenBankCDS.' + gcode_formatted + '.AA.fasta -t 30 --databases ' + args.databases + ' --evalue 1e-15')
		


def script_five(args):

	gcode_by_folder = { line.strip().split('\t')[0] : line.strip().split('\t')[-1] for line in open(args.output + '/Output/gcode_output.tsv', 'r') }
	for folder in os.listdir(args.output + '/Output'):
		if os.path.isdir(args.output + '/Output/' + folder):
			gcode_formatted = gcode_by_folder[folder][0].upper() + gcode_by_folder[folder].lower()[1:]
			if os.path.isfile(args.output + '/Output/' + folder + '/' + folder + '_GenBankCDS.Renamed.' + gcode_formatted + '.AA.fasta'):
				step5_cmd = 'python 5a_FinalizeName.py -in ' + args.output + '/Output/' + folder + '/DiamondOG/' + folder + '_GenBankCDS.Renamed.' + gcode_formatted + '.AA.fasta -n ' + folder
				os.system(step5_cmd)

	os.mkdir(args.output + '/Output/Intermediate')

	for file in os.listdir(args.output + '/Output'):
		if file != 'ReadyToGo' and file != 'Intermediate':
			os.system('mv ' + args.output + '/Output/' + file + ' ' + args.output + '/Output/Intermediate')

	os.system('python 5b_SummaryStats.py -i ' + args.output + '/Output -d ' + args.databases)


if __name__ == "__main__":

	args = get_args()

	if (args.first_script == 1 or args.script == 1) and not os.path.isdir(args.cds):
		print('\nIf starting at the first script, a valid path to a folder of nucleotide CDS files (which must end in .fasta) should be input using the --cds argument')
		exit()

	ten_digit_codes = []
	if args.first_script == 1 or args.script == 1:
		for file in os.listdir(args.cds):
			if file[10:] == '_GenBankCDS.fasta' and '.DS_Store' not in file:
				ten_digit_codes.append(file[:10])
	else:
		if not os.path.isdir(args.output + '/Output'):
			print('\nA folder called "Output" is not found at the given output path. Enter the correct path for --output or start from script 1.\n')
			exit()

	if(len(ten_digit_codes) > len(list(dict.fromkeys(ten_digit_codes)))):
		print('\nDuplicate 10-digit codes are not allowed. Aborting.\n')
		exit()

	for code in ten_digit_codes:
		for c, char in enumerate(code):
			if (c != 2 and c != 5 and char not in 'qwertyuiopasdfghjklzxcvbnmQWERTYUIOPASDFGHJKLZXCVBNM1234567890') or ((c == 2 or c == 5) and char != '_'):
				print('\n' + code + ' is an invalid 10-digit code sample identifier. It must of the format Op_me_hsap (Homo sapiens for example). Please ask for help if this does not make sense.\n')
				exit()

	if os.path.isdir(args.output + '/Output') and (args.first_script == 1 or args.script == 1):
		print('\nAn "Output" folder already exists at the given path. Please delete or rename this folder and try again.\n')
		exit()
	elif os.path.isdir(args.output + '/Output/Intermediate'):
		print('\nIt looks like this run is already complete. Try deleting/renaming the Output folder and try again.\n')
		exit()
	elif not os.path.isdir(args.output + '/Output'):
		os.mkdir(args.output + '/Output')
	
	scripts = [0, script_one, script_two, script_three, script_four, script_five]

	if args.script == -1:
		if args.first_script < args.last_script:
			for i in range(1 + args.last_script - args.first_script):
				print('\nRunning script ' + str(i + args.first_script) + '...\n')
				if i + args.first_script == 1:
					scripts[i + args.first_script](args, ten_digit_codes)
				else:
					scripts[i + args.first_script](args)
		else:
			print('\nInvalid script combination: the first script must be less than the last script. If you want to use only once script, use the --script argument.\n')
			exit()
	else:
		if args.script == 1:
			scripts[args.script](args, ten_digit_codes)
		else:
			scripts[args.script](args)

















