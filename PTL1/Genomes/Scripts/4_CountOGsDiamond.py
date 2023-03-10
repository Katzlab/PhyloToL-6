#!/usr/bin/env python3.5

##__Updated__: 19_09_2017
##__Author__: Xyrus Maurer-Alcala; maurerax@gmail.com
##__Usage__: python 3g_GCodeTranslate.py --help

##############################################################################
##                                                                          ##
## This scrip will categorize TRANSLATED CDSs into Homologous Gene Families ##
##                                                                          ##
##     Questions about Gene Family Binning/Source? SEE NOTES at Bottom!     ##
##                                                                          ##
##      E-mail Xyrus (author) for help if needed: maurerax@gmail.com        ##
##                                                                          ##
##############################################################################

import argparse, os, re, sys
from argparse import RawTextHelpFormatter, SUPPRESS
from distutils import spawn
from Bio import SeqIO


#----------------------------- Colors For Print Statements ------------------------------#
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
   

#------------------------------ UPDATE DIAMOND PATH BELOW! -------------------------------#
def check_diamond_path():
	### IF Diamond is IN YOUR PATH then no updating is needed...
	diamond_path = ''

	if diamond_path == '':
		diamond_path = spawn.find_executable("diamond")
		#diamond_path = /path/to/diamond
	else:
		pass

	if diamond_path == None:
		print (color.BOLD + '\n\nPlease open this script and check that you have included'\
		+ ' the PATH to the' + color.BLUE + ' "diamond" '+ color.END + color.BOLD\
		+ 'executable.\n\n' + color.END)
		print (color.BOLD + color.BLUE + 'LOOK FOR:\n\n' + color.RED\
		+'#------------------------------ UPDATE DIAMOND PATH BELOW! -------------------------------#'\
		+ color.BLUE + '\n\nThis is somewhere around lines 55 - 80...\n\n' + color.END)

		sys.exit()
	else:
		pass

	return diamond_path

#------------------------------- Main Functions of Script --------------------------------#

###########################################################################################
###--------------------- Parses and Checks Command-Line Arguments ----------------------###
###########################################################################################

def check_args():

	parser = argparse.ArgumentParser(description=
	color.BOLD + '\n\nThis script will categorize Contigs into'+color.ORANGE+' "Homologous" '\
	+color.END+color.BOLD+'Gene Families (OGs)\nbased on '+color.RED+'OrthoMCL'+color.END\
	+color.BOLD+"'s Gene Family Grouping\n\n\nNotes on this script and "+color.GREEN+\
	'OrthoMCL Families'+color.END+color.BOLD+' can be found\nat the bottom of '+color.GREEN\
	+'THIS script (4_CountOGsDiamond.py)\n'+color.END+usage_msg(), usage=SUPPRESS,
	formatter_class=RawTextHelpFormatter)
	
	required_arg_group = parser.add_argument_group(color.ORANGE+color.BOLD+'Required Options'+color.END)
 
	required_arg_group.add_argument('--input_file','-in', action='store',
	help=color.BOLD+color.GREEN+'Fasta file of Nucleotide sequences enriched \nwith'\
	' Eukaryotic protein coding transcripts'+color.END)
	required_arg_group.add_argument('--databases','-d', action='store',
	help=color.BOLD+color.GREEN+'Path to folder containing db_OG'+color.END)
	required_arg_group.add_argument('--evalue','-e', action='store',
	help=color.BOLD+color.GREEN+'Maximum OG assignment e-value'+color.END)

	optional_arg_group = parser.add_argument_group(color.ORANGE+color.BOLD+'Options'+color.END)

	optional_arg_group.add_argument('--threads','-t', default='2',
	help=color.BOLD+color.GREEN+' Number of threads to use for BLAST\n (default = 2)\n'+color.END)

	optional_arg_group.add_argument('-author', action='store_true',
	help=color.BOLD+color.GREEN+' Prints author contact information\n'+color.END)

	if len(sys.argv[1:]) == 0:
		print (parser.description)
		print ('\n')
		sys.exit()

	args = parser.parse_args()
	
	quit_eval = return_more_info(args)
	if quit_eval > 0:
		sys.exit()

	args.diamond = check_diamond_path()
	
	args.home_folder = '/'.join(args.input_file.split('/')[:-1]) + '/'
	
	args.tsv_out = args.home_folder + args.input_file.split('/')[-1].replace('CDS','CDS.Renamed').replace('.AA.fasta','_allOGCleanresults.tsv')
	
	args.aa_out = args.home_folder + args.input_file.split('/')[-1].replace('CDS','CDS.Renamed')
	args.ntd_out = args.home_folder + args.input_file.split('/')[-1].replace('CDS','CDS.Renamed').replace('AA','NTD')

	return args

		
###########################################################################################
###------------------------------- Script Usage Message --------------------------------###
###########################################################################################

def usage_msg():
	return (color.BOLD+color.RED+'\n\nExample usage:'+color.CYAN+' python 4_CountOGsDiamond.py'\
	' --input_file ../Stentor_coeruleus.WGS.CDS.Prep/Stentor_coeruleus.WGS.CDS.Universal.AA.fasta'+color.END)


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

	if args.input_file != None:
		if os.path.isfile(args.input_file) != False:
			if args.input_file.split('/')[-1] not in os.listdir('/'.join(args.input_file.split('/')[:-1])):
				print (color.BOLD+color.RED+'\nError:'+color.END+color.BOLD+' The provided Fasta file '\
				'('+color.DARKCYAN+args.input_file.split('/')[-1]+color.END+color.BOLD+')\ndoes not'\
				' exist or is incorrectly formatted.\n\nDouble-check then try again!\n\n'+color.END) 
				valid_arg += 1
			elif args.input_file.endswith('AA.fasta') != True:
				print (color.BOLD+'\n\nInvalid Fasta File! Only Fasta Files that were processed'\
				' with '+color.GREEN+'3g_GCodeTranslate.py '+color.END+color.BOLD+'are valid\n\n'\
				'However, to bypass that issue, Fasta Files MUST end with '+color.CYAN+\
				'"AA.fasta"\n\n'+color.END)
				valid_arg += 1
		else:
			print (color.BOLD+color.RED+'\nError:'+color.END+color.BOLD+' The provided Fasta file '\
			'('+color.DARKCYAN+args.input_file.split('/')[-1]+color.END+color.BOLD+')\ndoes not'\
			' exist or is incorrectly formatted.\n\nDouble-check then try again!\n\n'+color.END) 
			valid_arg += 1

	if os.path.isdir(args.databases + '/db_OG') != True:
		print (color.BOLD+color.RED+'\nError:'+color.END+color.BOLD+' Cannot find the '\
		+color.ORANGE+'db_OG Folder!\n\n'+color.END+color.BOLD+'Ensure that this folder '\
		'can be found in the main '+color.ORANGE+'Databases Folder'+color.END+color.BOLD\
		+'\n\nThen try once again\n\n.'+color.END)
		valid_arg += 1

	elif os.path.isfile(args.databases + '/db_OG/OGSout.dmnd') != True:
		print (color.BOLD+color.RED+'\nError:'+color.END+color.BOLD+' Cannot find the '\
		'Diamond formatted '+color.ORANGE+'Gene Family databases!\n\n'+color.END+color.BOLD+\
		'Ensure that they can be found in the '+color.ORANGE+'db_OG folder'+color.END+\
		color.BOLD+',\nwhich can be found in the main '+color.ORANGE+'Databases Folder'+\
		color.END+color.BOLD+'\n\nThen try once again.\n\n'+color.END)
		valid_arg += 1

	return valid_arg


###########################################################################################
###--------------------------- Does the Inital Folder Prep -----------------------------###
###########################################################################################

def prep_folders(args):

	OG_folder = '/'.join(args.input_file.split('/')[:-1])+'/DiamondOG/'
	
	if os.path.isdir(OG_folder) != True:
		os.system('mkdir '+OG_folder)		

		
###########################################################################################
###--------------------- Runs Diamond on Split OrthoMCL Databases ----------------------###
###########################################################################################

def OG_ublast(args):

	OG_diamond_cmd = args.diamond + ' blastp -q ' + args.input_file + ' -d ' + args.databases + '/db_OG/OGSout.dmnd --evalue ' + args.evalue + ' --subject-cover 0.5 --threads ' + args.threads + ' --outfmt 6 -o ' + args.input_file.split('.fas')[0] + '_allOGresults'	
	os.system(OG_diamond_cmd)


###########################################################################################
###--------------- Keeps the Single BEST Hit (HSP-score) Per Transcript ----------------###
###########################################################################################

def keep_best(args):
	print (color.BOLD+color.PURPLE+'\n\nProcessing OG-database results to keep only the BEST'\
	'\nmatch for each transcript\n\n'+color.END)
	
	inTSV = [i for i in open(args.input_file.split('.fas')[0]+'_allOGresults').read().split('\n') if i != '']
	
	inTSV.sort(key = lambda x: -float(x.split('\t')[-1]))
	
	keep = []
	for i in inTSV:
		if any(i.split('\t')[0] in j for j in keep) != True:
			keep.append(i)

	updated_lines = list(set([line.split('\t')[0]+'_'+'_'.join(line.split('\t')[1].split('_')[-2:])+\
	'\t'+'\t'.join(line.split('\t')[1:])+'\n' for line in keep]))
		
	with open(args.tsv_out, 'w+') as w:
		for i in updated_lines:
			w.write(i+'\n')
	

###########################################################################################
###-------- Copies and Updates Names of Transcripts With OG Hits to New Fasta ----------###
###########################################################################################

def update_fasta(args):

	print (color.BOLD+color.PURPLE+'Updating Sequence Names with their BEST OG hits\n\n'+color.END)

	keep = [i for i in open(args.tsv_out).read().split('\n') if i != '']

	keep_dict = {line.split('\t')[0].split('_OG5')[0]:line.split('\t')[0].split('_OG5')[0]+\
	'_OG5_'+line.split('\t')[1].split('_')[-1] for line in keep if 'OG5' in line.split('\t')[1]}
	
	protFasta = [seq_rec for seq_rec in SeqIO.parse(args.input_file,'fasta')]
	
	ntdFasta = [seq_rec for seq_rec in SeqIO.parse(args.input_file.replace('.AA.','.NTD.'),'fasta')]

	updated_prot_name = ['>'+keep_dict[i.description]+'\n'+str(i.seq).rstrip('*')+'\n' for i in protFasta if i.description in keep_dict.keys()]
	updated_ntd_name = ['>'+keep_dict[i.description]+'\n'+str(i.seq).rstrip('*')+'\n' for i in ntdFasta if i.description in keep_dict.keys()]

	with open(args.aa_out,'w+') as w:
		for i in updated_prot_name:
			w.write(i)

	with open(args.ntd_out,'w+') as x:
		for i in updated_ntd_name:
			x.write(i)			


##########################################################################################
###--------------------- Cleans up the Folder and Moves Final Files -------------------###
##########################################################################################

def clean_up(args):

	os.system('mv '+args.input_file.replace('.fasta','_allOGresults')+' '+args.home_folder+\
	'/DiamondOG')

	os.system('cp '+args.aa_out+' '+args.home_folder+'/DiamondOG/')
	os.system('cp '+args.ntd_out+' '+args.home_folder+'/DiamondOG/')
	os.system('cp '+args.tsv_out+' '+args.home_folder+'/DiamondOG/')


##########################################################################################
###----------------------------- Calls on Above Functions -----------------------------###
##########################################################################################	
				
def main():

	args = check_args()
	
	prep_folders(args)

	OG_ublast(args)
	
	keep_best(args)

	update_fasta(args)

	clean_up(args)

	print (color.BOLD+'Next Script is: '+color.GREEN+'5g_FinalizeName.py\n\n'+color.END)

main()

#----------------------------------------- NOTES -----------------------------------------#
#
# This script uses a "BLAST"-based approach to identify ANCIENT homologous gene families.
#
# Gene family designations were taken from OrthoMCL.org and serve as the database for 
# this script's gene family assignments. These gene family assignments are NON-EXHAUSTIVE
# and most Lineage-Specific families will be missed!
#
# If you have any questions contact Xyrus (author): maurerax@gmail.com