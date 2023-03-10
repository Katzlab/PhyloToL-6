import os, sys, re
import argparse


def get_args():

	parser = argparse.ArgumentParser(
                prog = 'PhyloToL v6.0 Part 1 for GenBank Genomes',
                description = "Updated January 19th, 2023 by Auden Cote-L'Heureux. Link to GitHub: https://github.com/AudenCote/PhyloToL_v6.0"
                )

	parser.add_argument('-s', '--script', default = -1, type = int, choices = { 1, 2, 3, 4, 5, 6 }, help = 'Script to run if you are only running one script')
	parser.add_argument('-1', '--first_script', default = -1, type = int, choices = { 1, 2, 3, 4, 5 }, help = 'First script to run')
	parser.add_argument('-2', '--last_script', default = -1, type = int, choices = { 2, 3, 4, 5, 6 }, help = 'First script to run')
	parser.add_argument('-c', '--cds', type = str, help = 'Path to a folder of nucleotide CDS. Each file name should start with a unique 10 digit code, and end in "_GenBankCDS.fasta", E.g. Op_me_hsap_GenBankCDS.fasta')
	parser.add_argument('-o', '--output', default = '../', type = str, help = 'An "Output" folder will be created at this directory to contain all output files. By default this folder will be created at the parent directory of the Scripts folder')
	parser.add_argument('-x', '--xplate_contam', action = 'store_true', help = 'Run cross-plate contamination removal (includes all files)')
	parser.add_argument('-g', '--genetic_code', type = str, help = 'If all of your taxa use the same genetic code, you may enter it here (to be used in script 4). Otherwise, stop after script 3 and fill in "gcode_output.tsv" before running script 4')
	parser.add_argument('-l', '--minlen', type = int, default = 200, help = 'Minimum CDS length')
	parser.add_argument('-d', '--databases', type = str, default = '../Databases', help = 'Path to databases folder (which should contain db_OG)')

	return parser.parse_args()


def script_one(args, ten_digit_codes):

	for file in os.listdir(args.cds):
		if file[10:] == '_GenBankCDS.fasta' and file[:10] in ten_digit_codes:
			os.system('python 1_RenameCDS.py -in ' + args.cds + '/' + file + ' -s GenBank -o ' + args.output + '/Output')


def script_two(args):

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
								
	with open(args.output + '/Output/gcode_output.tsv', 'w') as g:
		g.writelines('10 Digit Code\tIn-frame TGA Density\tIn-frame TAG Density\tIn-frame TAA Density\tGenetic Code\n') 
		for row in gcode_info:
			g.writelines(row[0] + '\t' + row[1] + '\t' + row[2] + '\t' + row[3] + '\n')


def script_three(args):

	valid_codes = ['universal', 'blepharisma', 'chilodonella', 'condylostoma', 'euplotes', 'peritrich', 'vorticella', 'mesodinium', 'tag', 'tga', 'taa', 'none']

	if args.genetic_code != None and args.genetic_code.lower() in valid_codes:
		for folder in os.listdir(args.output + '/Output'):
			if os.path.isfile(args.output + '/Output/' + folder + '/' + folder + '_GenBankCDS.Prepped.fasta'):
				os.system('python 3_GCodeTranslate.py -in ' + args.output + '/Output/' + folder + '/' + folder + '_GenBankCDS.Prepped.fasta -g ' + args.genetic_code.lower())
	else:
		lines = [line.strip().split('\t') for line in open(args.output + '/Output/gcode_output.tsv', 'r')]
		with open(args.output + '/Output/gcode_output.tsv', 'r') as g:
			for folder in os.listdir(args.output + '/Output'):
				if os.path.isfile(args.output + '/Output/' + folder + '/' + folder + '_GenBankCDS.Prepped.fasta'):
					for line in lines:
						if line[0] == folder and line[-1].lower() in valid_codes:
							os.system('python 3_GCodeTranslate.py -in ' + args.output + '/Output/' + folder + '/' + folder + '_GenBankCDS.Prepped.fasta -g ' + line[-1])
						elif line[-1].lower() not in valid_codes and line[-1] != 'Genetic Code':
							print('\n' + line[-1] + ' is not a valid genetic code. Skipping taxon ' + folder + '.\n')


def script_four(args):

	for folder in os.listdir(args.output + '/Output'):
		if os.path.isfile(args.output + '/Output/' + folder + '/' + folder + '_GenBankCDS.Universal.AA.fasta'):
			os.system('python 4_CountOGsDiamond.py -in ' + args.output + '/Output/' + folder + '/' + folder + '_GenBankCDS.Universal.AA.fasta -t 30 --databases ' + args.databases + ' --evalue 1e-15')
	


def script_five(args):

	for folder in os.listdir(args.output + '/Output'):
		if os.path.isfile(args.output + '/Output/' + folder + '/' + folder + '_GenBankCDS.Renamed.Universal.AA.fasta'):
			step5_cmd = 'python 5_FinalizeName.py -in ' + args.output + '/Output/' + folder + '/DiamondOG/' + folder + '_GenBankCDS.Renamed.Universal.AA.fasta -n ' + folder
			os.system(step5_cmd)

	os.mkdir(args.output + '/Output/Intermediate')

	for file in os.listdir(args.output + '/Output'):
		if file != 'ReadyToGo' and file != 'Intermediate':
			os.system('mv ' + args.output + '/Output/' + file + ' ' + args.output + '/Output/Intermediate')



if __name__ == "__main__":

	args = get_args()

	if (args.first_script == 1 or args.script == 1) and not os.path.isdir(args.cds):
		print('\nIf starting at the first script, a valid path to a folder of nucleotide CDS files (which must end in .fasta) should be input using the --cds argument')
		quit()

	ten_digit_codes = []
	if args.first_script == 1 or args.script == 1:
		for file in os.listdir(args.cds):
			if file[10:] == '_GenBankCDS.fasta':
				ten_digit_codes.append(file[:10])
	else:
		if not os.path.isdir(args.output + '/Output'):
			print('\nA folder called "Output" is not found at the given output path. Enter the correct path for --output or start from script 1.\n')

	if(len(ten_digit_codes) > len(list(dict.fromkeys(ten_digit_codes)))):
		print('\nDuplicate 10-digit codes are not allowed. Aborting.\n')
		quit()

	for code in ten_digit_codes:
		for c, char in enumerate(code):
			if (c != 2 and c != 5 and char not in 'qwertyuiopasdfghjklzxcvbnmQWERTYUIOPASDFGHJKLZXCVBNM1234567890') or ((c == 2 or c == 5) and char != '_'):
				print('\n' + code + ' is an invalid 10-digit code sample identifier. It must of the format Op_me_hsap (Homo sapiens for example). Please ask for help if this does not make sense.\n')
				quit()

	if os.path.isdir(args.output + '/Output') and (args.first_script == 1 or args.script == 1):
		print('\nAn "Output" folder already exists at the given path. Please delete or rename this folder and try again.\n')
		quit()
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
			quit()
	else:
		if args.script == 1:
			scripts[args.script](args, ten_digit_codes)
		else:
			scripts[args.script](args)

















