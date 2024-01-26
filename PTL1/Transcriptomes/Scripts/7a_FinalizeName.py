# Last updated Sept 29th, 2023
# Authors: Xyrus Maurer-Alcala and Auden Cote-L'Heureux

# This script does not process sequence data in any way. It only renames the outputs of 
# script 6 to the 10-digit taxon code which prefixes the file names, and then moves output
# 'ReadyToGo' files into a separate folder. It is intended to be run as part of the PhyloToL
# 6 Part 1 pipeline using the script wrapper.py.

import argparse, os, sys
from argparse import RawTextHelpFormatter,SUPPRESS

#----------------------- Solely to Make Print Statements Colorful -----------------------#

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
###--------------------- Parses and Checks Command-Line Arguments ----------------------###
###########################################################################################

def check_args():

	parser = argparse.ArgumentParser(description=
	color.BOLD + '\n\nThis script is intended to '+color.RED+'Rename '+color.END\
	+color.BOLD+'the core set of '+color.PURPLE+'ORFS\n'+color.END+color.BOLD+'with a valid '\
	+color.RED+'10-character code'+color.END+color.BOLD+' for use in the KatzLab\nPhylogenomic Pipeline'\
	+usage_msg(), usage=SUPPRESS, formatter_class=RawTextHelpFormatter)
	
	required_arg_group = parser.add_argument_group(color.ORANGE+color.BOLD+'Required Options'+color.END)
 
	required_arg_group.add_argument('--input_file','-in', action='store',
	help=color.BOLD+color.GREEN+' One of the Fasta files that is to be renamed\n'+color.END)
	required_arg_group.add_argument('--name','-n', action='store',
	help=color.BOLD+color.GREEN+' A valid 10-Character code for updating the data\n'+color.END)


	optional_arg_group = parser.add_argument_group(color.ORANGE+color.BOLD+'Options'+color.END)

	optional_arg_group.add_argument('-author', action='store_true',
	help=color.BOLD+color.GREEN+' Prints author contact information\n'+color.END)

	if len(sys.argv[1:]) == 0:
		print (parser.description)
		print ('\n')
		sys.exit()

	args = parser.parse_args()
	
	quit_eval = return_more_info(args)
	if quit_eval > 0:
		print ('\n')
		sys.exit()

	args.all_output_folder = '/'.join(args.input_file.split('/')[:-2])
	
	args.file_prefix = args.input_file.split('/')[-1].split('_Filtered.Final')[0]
	if 'fasta' in args.file_prefix:
		args.file_prefix = args.name

	args.r2g_aa = args.all_output_folder + '/ReadyToGo/ReadyToGo_AA/'
	args.r2g_ntd = args.all_output_folder + '/ReadyToGo/ReadyToGo_NTD/'
	args.r2g_tsv = args.all_output_folder + '/ReadyToGo/ReadyToGo_TSV/'
	

	return args


###########################################################################################
###------------------------------- Script Usage Message --------------------------------###
###########################################################################################

def usage_msg():
	return (color.BOLD+color.RED+'\n\nExample usage:'+color.CYAN+' python 7_FinalizeName.py'\
	' --input_file ../ToRename/Op_me_Xxma_Filtered.Final.AA.ORF.fasta --name Op_me_Xxma'+color.END)


##########################################################################################
###-------- Storage for LARGE (Annoying) Print Statements for Flagged Options ---------###
##########################################################################################

def return_more_info(args):

	valid_args = 0

	author = (color.BOLD+color.ORANGE+'\n\n\tQuestions/Comments? Email Xyrus (author) at'\
	' maurerax@gmail.com\n\n'+color.END)

	if args.author == True:
		print (author)
		valid_args += 1

	print(args.input_file)

	if args.input_file.endswith('AA.ORF.fasta'):
		args.input_NTD = args.input_file.replace('AA.ORF.fasta','NTD.ORF.fasta')
		args.input_AA = args.input_file
# 		args.input_TSV = ('/').join(args.input_file.split('/')[:-1])+'/SpreadSheets/'+args.input_file.split('/')[-1].replace('AA.ORF.fasta','allOGCleanresults.tsv')
		args.input_TSV = args.input_file.replace('AA.ORF.fasta','allOGCleanresults.tsv')

	elif args.input_file.endswith('NTD.ORF.fasta'):
		args.input_NTD = args.input_file
		args.input_AA = args.input_file.replace('NTD.ORF.fasta','AA.ORF.fasta')
# 		args.input_TSV = ('/').join(args.input_file.split('/')[:-1])+'/SpreadSheets/'+args.input_file.split('/')[-1].replace('NTD.ORF.fasta','allOGCleanresults.tsv')
		args.input_TSV = args.input_file.replace('AA.ORF.fasta','allOGCleanresults.tsv')
		print(args.input_TSV)

	if os.path.isfile(args.input_NTD) != True:
		print (color.BOLD+color.RED+'\nError:'+color.END+color.BOLD+' The provided Nucleotide '\
		'Fasta file ('+color.DARKCYAN+args.input_NTD.split('/')[-1]+color.END+color.BOLD+')\ndoes not'\
		' exist or is incorrectly formatted.\n\nDouble-check then try again!\n\n'+color.END) 
		valid_args += 1

	if os.path.isfile(args.input_AA) != True:
		print (color.BOLD+color.RED+'\nError:'+color.END+color.BOLD+' The provided Protein '\
		'Fasta file ('+color.DARKCYAN+args.input_AA.split('/')[-1]+color.END+color.BOLD+')\ndoes not'\
		' exist or is incorrectly formatted.\n\nDouble-check then try again!\n\n'+color.END) 
		valid_args += 1

	if os.path.isfile(args.input_TSV) != True:
		print (color.BOLD+color.RED+'\nError:'+color.END+color.BOLD+' The provided TSV '\
		' file ('+color.DARKCYAN+args.input_TSV.split('/')[-1]+color.END+color.BOLD+')\ndoes not'\
		' exist or is incorrectly formatted.\n\nDouble-check then try again!\n\n'+color.END) 
		valid_args += 1

	return valid_args

###########################################################################################
###-------------------- Double Checks Format for 10-Character Code ---------------------###
###########################################################################################

def check_code(args):
	
	check_name = args.name.split('_')
	
	if len(args.name) != 10:
		print (color.BOLD+'\n\nNew Species Prefix is not 10 characters long\n\n')
		print ('Three examples below:\n'+color.CYAN+'\n\tSr_ci_Cunc\n\n\tOp_me_Hsap\n\n\t'\
		'Am_ar_Ehis\n\n'+color.END) 
		sys.exit()

	elif args.name.count('_') != 2:
		print (color.BOLD+'\n\nCheck the format of your Species Prefix!\n\n')
		print ('Three examples below:\n'+color.CYAN+'\n\tSr_ci_Cunc\n\n\tOp_me_Hsap\n\n\t'\
		'Am_ar_Ehis\n\n'+color.END) 

		sys.exit()

	if len(check_name[0]) == 2 and len(check_name[1]) == 2 and len(check_name[2]) == 4:
		print (color.BOLD+"\n\nRenaming "+color.ORANGE+args.input_file.split('/')[-1]\
		.split('_Filtered')[0]+color.END+color.BOLD+"'s files with the following 10-character\n"\
		"code: "+color.CYAN+args.name+color.END+'\n')
	else:
		print (color.BOLD+'\n\nCheck the format of your Species Prefix!\n\n')
		print ('Three examples below:\n'+color.CYAN+'\n\tSr_ci_Cunc\n\n\tOp_me_Hsap\n\n\t'\
		'Am_ar_Ehis\n\n'+color.END) 
		sys.exit()

			
##########################################################################################
###------------------------- Creates Folders For Storing Data -------------------------###
##########################################################################################

def prep_folders(args):
	

	if os.path.isdir(args.all_output_folder + '/ReadyToGo/') != True:
		os.system('mkdir ' + args.all_output_folder + '/ReadyToGo')


	if os.path.isdir(args.r2g_ntd) != True:
		os.system('mkdir ' + args.r2g_ntd)
	if os.path.isdir(args.r2g_aa) != True:
		os.system('mkdir ' + args.r2g_aa)
	if os.path.isdir(args.r2g_tsv) != True:
		os.system('mkdir ' + args.r2g_tsv)

	if os.path.isdir(args.all_output_folder + '/' + args.file_prefix + '/Renamed') != True:
		os.system('mkdir ' + args.all_output_folder + '/' + args.file_prefix + '/Renamed')

###########################################################################################
###----------- Renames the NTD and AA CDSs with the Given 10-Character Code ------------###
###########################################################################################

def rename_paralogs(args):

	home_folder = args.all_output_folder + '/' + args.file_prefix + '/Renamed/'

	print (color.BOLD+'\nRenaming Translated (Protein) '+color.PURPLE+'ORFs\n'+color.END)
	renamed_Final_Prots = open(args.input_AA).read().replace('>','>'+args.name+'_XX_')
	
	print (color.BOLD+'\nRenaming Nucleotide '+color.PURPLE+'ORFs\n'+color.END)
	renamed_Final_Nucs = open(args.input_NTD).read().replace('>','>'+args.name+'_XX_')

	
	print (color.BOLD+'\nUpdating CDS Names in the Spreadsheet'+color.END)
	if '\n\n' in open(args.input_TSV).read():
		renamed_Final_tsv = args.name+'_XX_'+open(args.input_TSV).read().rstrip('\n')\
		.replace('\n\n','\n'+args.name+'_XX_')
	else:
		renamed_Final_tsv = args.name+'_XX_'+open(args.input_TSV).read().rstrip('\n')\
		.replace('\n','\n'+args.name+'_XX_')
		
	with open(home_folder+args.name+'_XX_'+args.input_AA.split('/')[-1],'w+') as w:
		w.write(renamed_Final_Prots)

	with open(home_folder+args.name+'_XX_'+args.input_NTD.split('/')[-1],'w+') as x:
		x.write(renamed_Final_Nucs)

	
	with open(home_folder+args.name+'_XX_'+args.input_TSV.split('/')[-1],'w+') as y:
		y.write(renamed_Final_tsv)


##########################################################################################
###-------------------- Cleans up the Folder and Moves Final Files --------------------###
##########################################################################################
def clean_up(args):

	home_folder = args.all_output_folder + '/' + args.file_prefix + '/Renamed/'

	os.system('cp '+home_folder+'*tsv '+args.r2g_tsv)

	os.system('cp '+home_folder+'*_XX_*AA.ORF.fasta '+args.r2g_aa)
	os.system('cp '+home_folder+'*_XX_*NTD.ORF.fasta '+args.r2g_ntd)

	os.system('cp '+home_folder+'*_XX_*tsv ' + args.all_output_folder + '/' + args.file_prefix)
	os.system('cp '+home_folder+'*_XX_*AA.ORF.fasta ' + args.all_output_folder + '/' + args.file_prefix)
	os.system('cp '+home_folder+'*_XX_*NTD.ORF.fasta ' + args.all_output_folder + '/' + args.file_prefix)
	
	os.system('rm ' + args.all_output_folder + '/ToRename/*'+args.file_prefix+'*')
	
	if os.path.isdir(args.all_output_folder + '/Finished/') != True:
		os.system('mkdir ' + args.all_output_folder + '/Finished')
	
	os.system('mv ' + args.all_output_folder + '/' + args.file_prefix + ' ' + args.all_output_folder + '/Finished')

###########################################################################################
###-------------------------------- Next Script Message --------------------------------###
###########################################################################################

def next_script(args):

	print (color.BOLD+'\nThere is no next script!\n\n'+color.END)

##########################################################################################
###--------------- Checks Command Line Arguments and Calls on Functions ---------------###
##########################################################################################
			
def main():

	args = check_args()
	
	check_code(args)
	
	prep_folders(args)
	
	rename_paralogs(args)
		
	clean_up(args)
	
	next_script(args)
	
main()
