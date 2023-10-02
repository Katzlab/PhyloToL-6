#Dependencies
import os, sys, re
import shutil
import argparse


def get_args():

	parser = argparse.ArgumentParser(
                prog = 'PhyloToL v6.0 Part 1 for Transcriptomes',
                description = "Updated September 29th, 2023 by Auden Cote-L'Heureux. Link to GitHub: https://github.com/AudenCote/PhyloToL_v6.0"
                )

	parser.add_argument('-s', '--script', default = -1, type = int, choices = { 1, 2, 3, 4, 5, 6, 7 }, help = 'Script to run if you are only running one script')
	parser.add_argument('-n', '--conspecific_names', type = str, help = 'A .txt or .tsv file with two tab-separated columns; the first should have 10 digit codes, the second species or other identifying names. This is used to determine which sequences to remove (only between "species") in cross-plate contamination assessment')
	parser.add_argument('-1', '--first_script', default = -1, type = int, choices = { 1, 2, 3, 4, 5, 6 }, help = 'First script to run')
	parser.add_argument('-2', '--last_script', default = -1, type = int, choices = { 2, 3, 4, 5, 6, 7 }, help = 'First script to run')
	parser.add_argument('-a', '--assembled_transcripts', type = str, help = 'Path to a folder of assembled transcripts, assembled by rnaSPAdes. Each assembled transcript file name should start with a unique 10 digit code, and end in "_assembledTranscripts.fasta", E.g. Op_me_hsap_assembledTranscripts.fasta')
	parser.add_argument('-o', '--output', default = '../', type = str, help = 'An "Output" folder will be created at this directory to contain all output files. By default this folder will be created at the parent directory of the Scripts folder')
	parser.add_argument('-x', '--xplate_contam', action = 'store_true', help = 'Run cross-plate contamination removal (includes all files)')
	parser.add_argument('-g', '--genetic_code', type = str, help = 'If all of your taxa use the same genetic code, you may enter it here (to be used in script 5). Alternatively, if you need to use a variety of genetic codes but know which codes to use, you may fill give here the path to a .txt or .tsv with two tab-separated columns, the first with the ten-digit codes and the second column with the corresponding genetics codes. Otherwise, stop at script 4 and fill in "gcode_output.tsv" before running script 5')
	parser.add_argument('-min', '--minlen', type = int, default = 200, help = 'Minimum transcript length')
	parser.add_argument('-max', '--maxlen', type = int, default = 12000, help = 'Maximum transcript length')
	parser.add_argument('-d', '--databases', type = str, default = '../Databases', help = 'Path to databases folder')

	return parser.parse_args()

		
#running the first script on all the bare files
def script_one(args, ten_digit_codes):
	for file in os.listdir(args.assembled_transcripts):
		if file[10:] == '_assembledTranscripts.fasta' and file[:10] in ten_digit_codes:
			os.system('python 1a_TranscriptLengthFilter.py --input_file ' + args.assembled_transcripts + '/' + file + ' --output_file ' + args.output + '/Output/' + file[:10] + ' --minLen ' + str(args.minlen) + ' --maxLen ' + str(args.maxlen) + ' --spades') #SPADES ARGUMENT??

	if args.xplate_contam:
		if not os.path.isfile(args.conspecific_names):
			print('\nERROR: If you are running cross-plate contamination, a file designating species assignments is required for the --conspecific_names argument\n')
			exit()
		else:
			os.system('python 1b_CrossPlateContamination.py ' + args.output + '/Output/XlaneBleeding ' + str(args.minlen) + ' ' + args.conspecific_names)


def script_two(args):

	for folder in os.listdir(args.output + '/Output/'):
		if os.path.isfile(args.output + '/Output/' + folder + '/SizeFiltered/' + folder + '.' + str(args.minlen) + 'bp.fasta'):
			os.system('python 2a_Identify_rRNA.py --input_file ' + args.output + '/Output/' + folder + '/SizeFiltered/' + folder + '.' + str(args.minlen) + 'bp.fasta --databases ' + args.databases)

			fasta_withBact = args.output + '/Output/' + folder + '/' + folder + '_NorRNAseqs.fasta'
			os.system('python 2b_Identify_Proks.py --input_file ' + fasta_withBact + ' --databases ' + args.databases)

#NEED TO SORT OUT FILE NAMES ETC. BELOW HERE

#running the third script
def script_three(args):

	for folder in os.listdir(args.output + '/Output'):
		if os.path.isfile(args.output + '/Output/' + folder + '/' + folder + '_WTA_EPU.fasta'):
			os.system('python 3_AssignOGs.py --input_file ' + args.output + '/Output/' + folder + '/' + folder + '_WTA_EPU.fasta --evalue 1e-15 --databases ' + args.databases)


#running the fourth script
def script_four(args):			
	for folder in os.listdir(args.output + '/Output'):
		if os.path.isfile(args.output + '/Output/' + folder + '/' + folder + '_WTA_EPU.Renamed.fasta'):
				os.system('python 4_InFrameStopCodonEstimator.py --input_file ' + args.output + '/Output/' + folder + '/' + folder + '_WTA_EPU.Renamed.fasta --databases ' + args.databases)
				
	#putting all of the gcode summaries produced by fourth script into a spreadsheet
	gcode_info = []
	for folder in os.listdir(args.output + '/Output'):
		if os.path.isfile(args.output + '/Output/' + folder + '/' + folder + '_WTA_EPU.Renamed_StopCodonStats.tsv'):
				with open(args.output + '/Output/' + folder + '/' + folder + '_WTA_EPU.Renamed_StopCodonStats.tsv') as f:
					for line in f:
						line_sep = line.split('\t')
						if line_sep[0] == 'Summary':
							gcode_info.append([folder, line_sep[6], line_sep[7], line_sep[8][:-1]])

	valid_codes = ['bleph','blepharisma','chilo','chilodonella','condy', 'condylostoma','none','eup','euplotes','peritrich','vorticella','ciliate','universal','taa','tag','tga','mesodinium']
	
	stop = False
	gcode_file = { }
	if args.genetic_code.endswith('.txt') or args.genetic_code.endswith('.tsv'):
		if os.path.isfile(args.genetic_code):
			for line in open(args.genetic_code):
				try:
					if line.split('\t')[1].strip().lower() in valid_codes:
						gcode_file.update({ line.split('\t')[0] : line.split('\t')[1].strip() })
					else:
						print('Genetic code ERROR -- ' + line.split('\t')[1].strip() + ' is not a valid genetic code. Please fill out the "gcode_output.tsv" file and continue with script 5.')
				except IndexError:
					print('\nGenetic code ERROR -- it looks like you tried to enter a .txt/.tsv file, but it is improperly formatted. Stopping after script 4; you may fill out the file gcode_output.tsv and continue with script 5.\n')
					stop = True
		else:
			print('\nGenetic code ERROR -- it looks like you tried to enter a .txt/.tsv file, but it could not be found. Stopping after script 4; you may fill out the file gcode_output.tsv and continue with script 5.\n')
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
					g.writelines(row[0] + '\t' + row[1] + '\t' + row[2] + '\t' + row[3]+ '\t' + gcode_file[row[0]] + '\n')
				except KeyError:
					g.writelines(row[0] + '\t' + row[1] + '\t' + row[2] + '\t' + row[3]+ '\t' + 'Universal' + '\n')
					print('\nDefaulting to Universal genetic code for taxon ' + row[0] + '\n')
					#print('\nGenetic code ERROR -- it looks like you tried to enter a .txt/.tsv file, but a taxon is missing. Stopping after script 4; you may fill out the file gcode_output.tsv and continue with script 5.\n')
					#stop = True
			else:	
				stop = True

	if stop or args.genetic_code == None:
		print('\nStopping after script 4 because genetic code information is incomplete. Please fill out the file "gcode_output.tsv" and continue with script 5.\n')
		exit()


def script_five(args):
	
	valid_codes = ['bleph','blepharisma','chilo','chilodonella','condy', 'condylostoma','none','eup','euplotes','peritrich','vorticella','ciliate','universal','taa','tag','tga','mesodinium']
	
	lines = [line.strip().split('\t') for line in open(args.output + '/Output/gcode_output.tsv', 'r')]
	with open(args.output + '/Output/gcode_output.tsv', 'r') as g:
		for folder in os.listdir(args.output + '/Output'):
			if os.path.isfile(args.output + '/Output/' + folder + '/' + folder + '_WTA_EPU.Renamed.fasta'):
				for line in lines:
					if line[0] == folder and line[-1].lower() in valid_codes:
						os.system('python 5_GCodeTranslate.py --input_file ' + args.output + '/Output/' + folder + '/' + folder + '_WTA_EPU.Renamed.fasta --genetic_code ' + line[-1])
					elif line[-1].lower() not in valid_codes and line[-1] != 'Genetic Code':
						print('\n' + line[-1] + ' is not a valid genetic code. Skipping taxon ' + folder + '.\n')


def script_six(args):

	prefixes = []
	for file in os.listdir(args.output + '/Output'):
		if file.endswith('_AA.ORF.fasta'): 
			prefixes.append(file[:10])
	
	unique_prefixes = list(dict.fromkeys(prefixes))

	hook_fasta = ''
	for file in os.listdir(args.databases + '/db_OG'):
		if file.split('.')[-1] in ('fasta', 'fas', 'fa', 'faa'):
			hook_fasta = args.databases + '/db_OG/' + file

	if hook_fasta == '':
		print('\nNo .fasta file could be found containing Hook sequences. This should be supplied along with the .dmnd-formatted database file in the Databases/db_OG folder. Quitting before script 6.\n')

	for prefix in unique_prefixes:
		os.system('python 6_FilterPartials.py --file_prefix ' + args.output + '/Output/' + prefix + ' --hook_fasta ' + hook_fasta)
		

def script_seven(args):

	for file in os.listdir(args.output + '/Output/ToRename'):
		if '.AA.ORF.fasta' in file:
			os.system('python 7a_FinalizeName.py --input_file ' + args.output + '/Output/ToRename/' + file + ' --name ' + file[:10])

	os.mkdir(args.output + '/Output/Intermediate')

	for file in os.listdir(args.output + '/Output'):
		if file != 'ReadyToGo' and file != 'Intermediate':
			os.system('mv ' + args.output + '/Output/' + file + ' ' + args.output + '/Output/Intermediate')

	os.system('python 7b_SummaryStats.py -i ' + args.output + '/Output -d ' + args.databases)


if __name__ == "__main__":

	args = get_args()

	if (args.first_script == 1 or args.script == 1) and (args.assembled_transcripts == None or not os.path.isdir(args.assembled_transcripts)):
		print('\nERROR: If starting at the first script, a valid path to a folder of assembled transcript files (which must end in .fasta, .fa, or .fna) should be input using the --assembled_transcripts argument')
		quit()

	if args.genetic_code == None and args.script == -1:
		if args.first_script < 5 and args.last_script >= 5:
			print('\nERROR: You cannot run script 5 without giving a genetic code! If all of the taxa in the run use the same genetic code, then use the --genetic_code argument (e.g. -g Universal). Otherwise, stop after script 4, fill out the spreadsheet called "gcode_translate.tsv," and then run scripts 5-7. If this does not make sense, please ask for help.')
			quit()

	ten_digit_codes = []
	if args.first_script == 1 or args.script == 1:
		for file in os.listdir(args.assembled_transcripts):
			if file[10:] == '_assembledTranscripts.fasta':
				ten_digit_codes.append(file[:10])
	else:
		if not os.path.isdir(args.output + '/Output'):
			print('\nERROR: A folder called "Output" is not found at the given output path. Enter the correct path for --output or start from script 1.\n')
			quit()

	if(len(ten_digit_codes) > len(list(dict.fromkeys(ten_digit_codes)))):
		print('\nERROR: Duplicate 10-digit codes are not allowed.\n')
		quit()

	for code in ten_digit_codes:
		for c, char in enumerate(code):
			if (c != 2 and c != 5 and char not in 'qwertyuiopasdfghjklzxcvbnmQWERTYUIOPASDFGHJKLZXCVBNM1234567890') or ((c == 2 or c == 5) and char != '_'):
				print('\nERROR: ' + code + ' is an invalid 10-digit code sample identifier. It must of the format Op_me_hsap (Homo sapiens for example). Please ask for help if this does not make sense.\n')
				quit()

	if os.path.isdir(args.output + '/Output') and (args.first_script == 1 or args.script == 1):
		print('\nERROR: An "Output" folder already exists at the given path. Please delete or rename this folder and try again.\n')
		quit()
	elif os.path.isdir(args.output + '/Output/Intermediate'):
		print('\nIt looks like this run is already complete. Try deleting/renaming the Output folder and try again.\n')
		quit()
	elif not os.path.isdir(args.output + '/Output'):
		os.mkdir(args.output + '/Output')

	scripts = [0, script_one, script_two, script_three, script_four, script_five, script_six, script_seven]

	if args.script == -1:
		if args.first_script < args.last_script:
			for i in range(1 + args.last_script - args.first_script):
				print('\nRunning script ' + str(i + args.first_script) + '...\n')
				if i + args.first_script == 1:
					if len(ten_digit_codes) == 0:
						print('\nNo properly-named assembled transcripts files found.\n')
						quit()
					else:
						scripts[i + args.first_script](args, ten_digit_codes)
				else:
					scripts[i + args.first_script](args)
		else:
			print('\nERROR: Invalid script combination: the first script must be less than the last script. If you want to use only once script, use the --script argument.\n')
			quit()
	else:
		if args.script == 1:
			if len(ten_digit_codes) == 0:
				print('\nNo properly-named assembled transcripts files found.\n')
				quit()
			else:
				scripts[args.script](args, ten_digit_codes)
		else:
			scripts[args.script](args)

		
		
						
				



















