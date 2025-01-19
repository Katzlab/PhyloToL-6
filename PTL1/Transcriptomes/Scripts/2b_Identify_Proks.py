# Last updated Sept. 2023
# Authors: Xyrus Maurer-Alcala and Auden Cote-L'Heureux

# This script is intended to identify likely prokarotic (contaminant) sequences. It does 
# this by similarity-searching against a reference database of eukaryote and prokaryote
# sequences, and it labels the output sequences with an "E" (likely eukaryotic), "P" (likely
# prokaryotic) or "U" (Unknown) in the sequence ID. This is done by comparing e-values: if
# a sequence hits a eukaryotic sequence with an e-value >100 times that of its best hit
# to a prokaryotic sequence, it is labeled with an "E"; if it's best hit to a prokaryotic
# sequence has an e-value >1000 times that of its best hit to a eukaryotic sequence, it is
# labeled with a "P". Anything else gets a "U". This script should be run as part of the 
# EukPhylo Part 1 pipeline using the script wrapper.py.

# Prior to running this script, ensure that you have run scripts 1a (and optionally
# script 1b) and 2a, and that your prokaryote and reference databases (or the default 
# ones provided on the GitHub) is in the proper database folder 
# (Databases/BvsE/eukout.dmnd and micout.dmnd).


import argparse, os, sys
from argparse import RawTextHelpFormatter,SUPPRESS
from distutils import spawn
from Bio import SeqIO


#------------------------------ Colors For Print Statements ------------------------------#
class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   ORANGE = '\033[38;5;214m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'


#------------------------------- Main Functions of Script --------------------------------#

###########################################################################################
###---------------------------- UPDATE DIAMOND PATH BELOW! -----------------------------###
###########################################################################################
	## IF Diamond is IN YOUR PATH then no updating is needed...

def check_diamond_path():

	diamond_path = ''

	if diamond_path == '':
		diamond_path = spawn.find_executable("diamond")
		#diamond_path = '/path/to/diamond'
	else:
		pass

	if diamond_path == None:
		print (color.BOLD + '\n\nPlease open this script and check that you have included'\
		+' the PATH to the'+color.BLUE+' "Diamond" '+color.END+color.BOLD+'executable.\n\n'+color.END)
		print (color.BOLD+color.BLUE+'LOOK FOR:\n\n'+color.RED\
		+'#------------------------------ UPDATE DIAMOND PATH BELOW! -------------------------------#'\
		+color.BLUE+'\n\nThis is somewhere around lines 50 - 80...\n\n'+color.END)

		sys.exit()
	else:
		pass
	
	return diamond_path
	
	
###########################################################################################
###--------------------- Parses and Checks Command-Line Arguments ----------------------###
###########################################################################################
	
def check_args():

	parser = argparse.ArgumentParser(description=
	color.BOLD + '\nThis script will categorize Contigs as'+color.ORANGE+' STRONGLY '+color.END\
	+color.BOLD+color.RED+'Eukaryotic \nOR Prokaryotic'+color.END+color.BOLD+' using a set of Proteins'\
	' from diverse\n'+color.ORANGE+'Eukaryotes, Bacteria and Archaea'+color.END\
	+color.BOLD+'.'+color.END+usage_msg(), usage=SUPPRESS,formatter_class=RawTextHelpFormatter)
	
	required_arg_group = parser.add_argument_group(color.ORANGE+color.BOLD+'Required Options'+color.END)
 
	required_arg_group.add_argument('--input_file','-in', action='store',
	help=color.BOLD+color.GREEN+'Fasta file of Nucleotide sequences (with rRNAs removed)'+color.END)
	required_arg_group.add_argument('--databases','-d', action='store',
	help=color.BOLD+color.GREEN+"Path to databases"+color.END)

	optional_arg_group = parser.add_argument_group(color.ORANGE+color.BOLD+'Options'+color.END)
	optional_arg_group.add_argument('-author', action='store_true',
	help=color.BOLD+color.GREEN+' Print author contact information\n'+color.END)

	if len(sys.argv[1:]) == 0:
		print (parser.description)
		print ('\n')
		sys.exit()

	args = parser.parse_args()
	
	quit_eval = return_more_info(args)
	if quit_eval > 0:
		sys.exit()
		
	return args


###########################################################################################
###------------------------------- Script Usage Message --------------------------------###
###########################################################################################

def usage_msg():
	return (color.BOLD+color.RED+'\n\nExample usage:'+color.CYAN+' python 2b_remove_Bact.py --input_file'\
	' ../Op_me_Xxma/Op_me_Xxma_NorRNAseqs.fasta'+color.END)


##########################################################################################
###-------- Storage for LARGE (Annoying) Print Statements for Flagged Options ---------###
##########################################################################################

def return_more_info(args):

	valid_arg = 0

	author = (color.BOLD+color.ORANGE+'\n\n\tQuestions/Comments? Email Xyrus (author) at'\
	' maurerax@gmail.com\n\n'+color.END)

	if args.author == True:
		print (author)
		valid_arg += 1

	print(args.input_file)

	if args.input_file != None:
		if os.path.isfile(args.input_file) != False:
			if args.input_file.split('/')[-1] not in os.listdir('/'.join(args.input_file.split('/')[:-1])):
				print (color.BOLD+color.RED+'\nError:'+color.END+color.BOLD+' The provided Fasta file '\
				'('+color.DARKCYAN+args.input_file.split('/')[-1]+color.END+color.BOLD+')\ndoes not'\
				' exist or is incorrectly formatted.\n\nDouble-check then try again!\n\n'+color.END) 
				valid_arg += 1
			elif args.input_file.endswith('NorRNAseqs.fasta') != True:
				print (color.BOLD+'\n\nInvalid Fasta File! Only Fasta Files that were processed'\
				' with '+color.GREEN+'2a_remove_rRNA.py '+color.END+color.BOLD+'are valid\n\n'\
				'However, to bypass that issue, Fasta Files MUST end with '+color.CYAN+\
				'"NorRNAseqs.fas"\n\n'+color.END)
				valid_arg += 1
		else:
			print (color.BOLD+color.RED+'\nError:'+color.END+color.BOLD+' The provided Fasta file '\
			'('+color.DARKCYAN+args.input_file.split('/')[-1]+color.END+color.BOLD+')\ndoes not'\
			' exist or is incorrectly formatted.\n\nDouble-check then try again!\n\n'+color.END) 
			valid_arg += 1

	if os.path.isdir(args.databases + '/db_BvsE') != True:
		print (color.BOLD+color.RED+'\nError:'+color.END+color.BOLD+' Cannot find the '\
		+color.ORANGE+'db_BvsE Folder!\n\n'+color.END+color.BOLD+'Ensure that this folder '\
		'can be found in the main '+color.ORANGE+'Databases Folder'+color.END+color.BOLD\
		+'\n\nThen try once again.')
		valid_arg += 1
	elif os.path.isfile(args.databases + '/db_BvsE/eukout.dmnd') != True:
		print (color.BOLD+color.RED+'\nError:'+color.END+color.BOLD+' Cannot find the '\
		'Diamond formatted '+color.ORANGE+'Eukaryotic Protein database!\n\n'+color.END+color.BOLD+\
		'Ensure that it can be found in the '+color.ORANGE+'db_BvsE folder'+color.END+\
		color.BOLD+',\nwhich can be found in the main '+color.ORANGE+'Databases Folder'+\
		color.END+color.BOLD+'\n\nThen try once again.'+color.END)
	elif os.path.isfile(args.databases + '/db_BvsE/micout.dmnd') != True:
		print (color.BOLD+color.RED+'\nError:'+color.END+color.BOLD+' Cannot find the '\
		'Diamond formatted '+color.ORANGE+'Bacterial/Archaeal Protein database!\n\n'+color.END+color.BOLD+\
		'Ensure that it can be found in the '+color.ORANGE+'db_BvsE folder'+color.END+\
		color.BOLD+',\nwhich can be found in the main '+color.ORANGE+'Databases Folder'+\
		color.END+color.BOLD+'\n\nThen try once again.'+color.END)

		valid_arg += 1

	return valid_arg
	
	
###########################################################################################
###--------------------------- Does the Inital Folder Prep -----------------------------###
###########################################################################################

def prep_folders(args):

	BvE_folder = '/'.join(args.input_file.split('/')[:-1]) + '/BvE/'
	
	if os.path.isdir(BvE_folder) != True:
		os.system('mkdir '+BvE_folder)

		
###########################################################################################
###---------------- Runs Diamond on Bact and Euk small RefSeq Databases ----------------###
###########################################################################################

def ublast_BvE(args, diamond_path):

	BvE_folder = '/'.join(args.input_file.split('/')[:-1]) + '/BvE/'
	mic_output = args.input_file.split('/')[-1]+'micresults.'
	euk_output = args.input_file.split('/')[-1]+'eukresults.'

	print(color.BOLD+'\n\n"BLAST"-ing against PROK database using DIAMOND: ' + color.DARKCYAN + 'micout.dmnd' + color.END + '\n\n')

	Prok_diamond_cmd = diamond_path + ' blastx -q ' + args.input_file + ' --max-target-seqs 1 -d ' + args.databases + '/db_BvsE/micout.dmnd --evalue 1e-5 --threads 60 --outfmt 6 -o ' + BvE_folder + 'allmicresults.tsv'

	os.system(Prok_diamond_cmd)

	print(color.BOLD+'\n\n"BLAST"-ing against EUK database using DIAMOND: ' + color.DARKCYAN + 'eukout.dmnd' + color.END + '\n\n')

	Euk_diamond_cmd = diamond_path + ' blastx -q ' + args.input_file + ' --max-target-seqs 1 -d ' + args.databases + '/db_BvsE/eukout.dmnd --evalue 1e-5 --threads 60 --outfmt 6 -o ' + BvE_folder + 'alleukresults.tsv'

	os.system(Euk_diamond_cmd)



###########################################################################################
###---------------- Compares Bacterial and Euk Hits for Classification -----------------###
###########################################################################################

def compare_hits(args):

	BvE_folder = '/'.join(args.input_file.split('/')[:-1]) + '/BvE/'
	
	EukDict = {}
	ProkDict = {}
	CompDict = {}	
	
	inFasta = [seq_rec for seq_rec in SeqIO.parse(args.input_file,'fasta')]
	
	for seq_rec in inFasta:
			EukDict[seq_rec.description] = ''
			ProkDict[seq_rec.description] = ''
			CompDict[seq_rec.description] = []

	inEukHits = [i for i in open(BvE_folder + 'alleukresults.tsv').readlines()]
	inEukHits.sort(key=lambda x: (float(x.split('\t')[-2]), -int(x.split('\t')[3])))
	
	inProkHits = [i for i in open(BvE_folder + 'allmicresults.tsv').readlines()]
	inProkHits.sort(key=lambda x: (float(x.split('\t')[-2]), -int(x.split('\t')[3])))
	
	for i in inEukHits:
		if EukDict[i.split('\t')[0]] == '':
			EukDict[i.split('\t')[0]] = float(i.split('\t')[-2])
	
	for i in inProkHits:
		if ProkDict[i.split('\t')[0]] == '':
			ProkDict[i.split('\t')[0]] = float(i.split('\t')[-2])
	
	for k in CompDict.keys():
		if EukDict[k] != '':
			CompDict[k].append(EukDict[k])
		else:
			CompDict[k].append('no hit')
		if ProkDict[k] != '':
			CompDict[k].append(ProkDict[k])
		else:
			CompDict[k].append('no hit')

	for k, v in CompDict.items():

	### Contigs lacking STRONG Eukaryotic OR Prokaryotic Hits
		if v[0] == 'no hit' and v[1] == 'no hit':
			CompDict[k].append('UNDETERMINED')

	### Contigs lacking STRONG Eukaryotic with a Prokaryotic Hit	
		elif v[0] != 'no hit' and v[1] == 'no hit':
			CompDict[k].append('EUKARYOTIC')

	### Contigs with a Eukaryotic but without a Prokaryotic Hit		
		elif v[0] == 'no hit' and v[1] != 'no hit':
			CompDict[k].append('PROKARYOTIC')

	### Uses Basic math to determine if contigs with are MORE Eukaryotic than Prokaryotic				
		else:
			try:
				prok_euk_ratio = float(v[1])/float(v[0])
				euk_prok_ratio = float(v[0])/float(v[1])
	
				if prok_euk_ratio >= 100:
					CompDict[k].append('EUKARYOTIC')
	
				elif  euk_prok_ratio >= 1000:
					CompDict[k].append('PROKARYOTIC')
	
				else:
					CompDict[k].append('UNDETERMINED')
	
			except:
				CompDict[k].append('divide by zero')
		
	with open(BvE_folder + 'comparisons.txt','w+') as w:
		for k, v in CompDict.items():
			w.write(k+':'+':'.join([str(i) for i in v])+'\n')

	BvE_folder = '/'.join(args.input_file.split('/')[:-1]) + '/BvE/'
	BvE_output_base = BvE_folder+args.input_file.split('/')[-1].split('.fas')[0]
	
### Gathers the sequences and categorizes them	
	Euk_Fasta = sorted((seq_rec for seq_rec in inFasta if CompDict[seq_rec.description][-1] == 'EUKARYOTIC'), key=lambda x: -int(len(x.seq)))
	Prok_Fasta = sorted((seq_rec for seq_rec in inFasta if CompDict[seq_rec.description][-1] == 'PROKARYOTIC'), key=lambda x: -int(len(x.seq)))
	Und_Fasta = sorted((seq_rec for seq_rec in inFasta if CompDict[seq_rec.description][-1] == 'UNDETERMINED'), key=lambda x: -int(len(x.seq)))
	Zero_Fasta = sorted((seq_rec for seq_rec in inFasta if CompDict[seq_rec.description][-1] == 'divide by zero'), key=lambda x: -int(len(x.seq)))


### Writes out all of the categorized sequences
	with open(args.input_file.split('NorRNA')[0] + 'WTA_EPU.fasta', 'w') as epu:
		with open(BvE_output_base+'.Not_Bact.fasta','w+') as nb:
			for euk_seq in Euk_Fasta:
				nb.write('>' + euk_seq.description + '\n' + str(euk_seq.seq) + '\n')
				epu.write('>' + euk_seq.description + '_E' + '\n' + str(euk_seq.seq) + '\n')


		with open(BvE_output_base+'.Bact_Hit.fasta','w+') as pr:
			for prok_seq in Prok_Fasta:
				pr.write('>' + prok_seq.description + '\n' + str(prok_seq.seq) + '\n')		
				epu.write('>' + prok_seq.description + '_P' + '\n' + str(prok_seq.seq) + '\n')

		with open(BvE_output_base+'.Undetermined.fasta','w+') as und:
			for und_seq in Und_Fasta:
				und.write('>' + und_seq.description + '\n' + str(und_seq.seq) + '\n')
				epu.write('>' + und_seq.description + '_U' + '\n' + str(und_seq.seq) + '\n')

		if len(Zero_Fasta) != 0:
			with open(BvE_output_base+'.DivideByZero.fasta','w+') as w:
				for zero_seq in Zero_Fasta:
					w.write('>' + zero_seq.description + '\n' + str(zero_seq.seq) + '\n')
					epu.write('>' + zero_seq.description + '_U' + '\n' + str(zero_seq.seq) + '\n')
		else:
			pass
		
	return str(len(Euk_Fasta)), str(len(Prok_Fasta)), str(len(Und_Fasta))	

###########################################################################################
###-------------------------------- Next Script Message --------------------------------###
###########################################################################################

def next_script(args):

	print (color.BOLD+'\nLook for '+color.DARKCYAN+args.input_file.split('/')[-1]\
	.split('NorRNA')[0]+'WTA_EPU.fasta'+color.END+color.BOLD+' in the '\
	+args.input_file.split('/')[1]+' Folder\n\n' + color.END)
	print (color.BOLD + 'Next Script is: ' + color.GREEN + '3_CountOGsDiamond.py\n\n'+ color.END)


##########################################################################################
###--------------------- Cleans up the Folder and Moves Final Files -------------------###
##########################################################################################

def clean_up(args):

	home_folder = '/'.join(args.input_file.split('/')[:-1])
	
	os.system('cp '+home_folder+'/*WTA_EPU.fasta '+home_folder+'/BvE/')
	os.system('mv '+home_folder+'/*NorRNA*fasta '+home_folder+'/rRNA_Removal/')


##########################################################################################
###--------------- Checks Command Line Arguments and Calls on Functions ---------------###
##########################################################################################

def main():

	usearch_path = check_diamond_path()

	args = check_args()

	prep_folders(args)
	
	ublast_BvE(args, usearch_path)

	Euk_Contigs, Prok_Contigs, Und_Contigs = compare_hits(args)
		
	clean_up(args)
	
	next_script(args)

main()
