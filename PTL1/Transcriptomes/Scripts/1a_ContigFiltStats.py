#!/usr/bin/env python3.6

##__Updated__: 01_04_2023
##__Author__: Xyrus Maurer-Alcala; maurerax@gmail.com
##__Usage__: python 1_ContigFiltStats.py
##__Options__: python 1_ContigFiltStats.py --help

##########################################################################################
## This script is intended to remove small transcripts or small contigs below a given   ##
## minimum size from a transcriptome assembly.                                          ##
##                                                                                      ##
## Prior to running this script, ensure the following:                                  ##
## 1. You have assembled your transcriptome and COPIED the 'assembly' file              ##
##    (contigs.fasta, or scaffolds.fasta) to the PostAssembly Folder                    ##
##                                                                                      ##
##                              COMMAND Example Below                                   ##
##                                                                                      ##
##          E-mail Xyrus (author) for help if needed: maurerax@gmail.com                ##
##                                                                                      ##
##                          Next Script(s) to Run:                                      ##
##  AutoBactVsEuk.py (removes SSU then Bact) or 2a_removeSSU.py then 2b_removeBact.py   ##
##                                                                                      ##
##########################################################################################


import argparse, os, sys
from argparse import RawTextHelpFormatter,SUPPRESS
from Bio import SeqIO
from Bio.SeqUtils import GC


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


#------------------------------- Main Functions of Script --------------------------------#

###########################################################################################
###--------------------- Parses and Checks Command-Line Arguments ----------------------###
###########################################################################################

def check_args():

	parser = argparse.ArgumentParser(description=
	color.BOLD+'\nThis script will remove Contigs (and provide a summary of statistics)'\
	+'\nfrom your Assembly that are shorter than a given length.'+color.ORANGE+\
	'\n\nA good minimum length to start with is 200bp.'+color.END+color.BOLD+\
	'\n\nThe minimum length value should be adjusted for your data sets.\n'+color.END+usage_msg(), 
	usage=SUPPRESS,formatter_class=RawTextHelpFormatter)
	
	required_arg_group = parser.add_argument_group(color.ORANGE+color.BOLD+'Required Options'+color.END)

	required_arg_group.add_argument('--input_file','-in', action='store',
	help=color.BOLD+color.GREEN+" Fasta file of Protein/Nucleotide sequences\n"+color.END)

	required_arg_group.add_argument('--output_file','-out',
	help=color.BOLD+color.GREEN+" Desired Output Name\n\n"+color.END)

	required_arg_group.add_argument('--minLen','-min', default=200, type=int,
	help=color.BOLD+color.GREEN+" Minimum number of base pairs for contigs\n (default = 200)"+color.END)
	required_arg_group.add_argument('--maxLen','-max', default=15000, type=int,
	help=color.BOLD+color.GREEN+" Minimum number of base pairs for contigs\n (default = 15000)"+color.END)

	optional_arg_group = parser.add_argument_group(color.ORANGE+color.BOLD+'Options'+color.END)

	optional_arg_group.add_argument('--spades','-spades', action='store_true',
	help=color.BOLD+color.GREEN+'rnaSPAdes transcriptome assembly\n'+color.END)

	optional_arg_group.add_argument('--genbank','-gb', action='store_true',
	help=color.BOLD+color.GREEN+'Assembly from Genbank\n (Will include Accession Number in'\
	' contig name)\n'+color.END)

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

	args = parser.parse_args()
	
	return args


###########################################################################################
###------------------------------- Script Usage Message --------------------------------###
###########################################################################################

def usage_msg():
	return color.BOLD+color.RED+'\n\nExample usage:'+color.CYAN+' python 1_ContigFiltStats.py'\
	' --input_file ../Op_me_Xxma_rnaSPAdes_scaffolds_15_05.fasta --output_file '\
	'Op_me_Xxma --minLen 200 --spades'+color.END


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
		else:
			print (color.BOLD+color.RED+'\nError:'+color.END+color.BOLD+' The provided Fasta file '\
			'('+color.DARKCYAN+args.input_file.split('/')[-1]+color.END+color.BOLD+')\ndoes not'\
			' exist or is incorrectly formatted.\n\nDouble-check then try again!\n\n'+color.END) 
			valid_arg += 1
	
	if args.output_file == None:
		valid_arg += 1

	return valid_arg

###########################################################################################
###--------------------------- Does the Inital Folder Prep -----------------------------###
###########################################################################################

def prep_folders(args):
	Home_folder_name = args.output_file
	
	if os.path.isdir(args.output_file) != True:
		os.system('mkdir ' + args.output_file)

	if os.path.isdir(args.output_file + '/OriginalFasta/') != True:
		os.system('mkdir ' + args.output_file +'/OriginalFasta/')

	if os.path.isdir(args.output_file + '/SizeFiltered/') != True:
		os.system('mkdir ' + args.output_file +'/SizeFiltered/')

	if os.path.isdir('/'.join(args.output_file.split('/')[:-1]) + '/XlaneBleeding/') != True:
		os.system('mkdir ' + '/'.join(args.output_file.split('/')[:-1]) + '/XlaneBleeding/')


###########################################################################################
###---------- Renames the Contigs, Writes them out, and Calculates Basic Info ----------###
###########################################################################################

def rename_Transcriptome(args):
	
	home_folder = args.output_file + '/SizeFiltered/'

	print (color.BOLD+'\n\nPrepping '+color.GREEN+args.input_file.split('/')[-1]+color.END)

	inFasta = [i for i in SeqIO.parse(args.input_file,'fasta') if len(i.seq) >= args.minLen and len(i.seq) <= args.maxLen]
	inFasta.sort(key=lambda seq_rec: -len(seq_rec.seq))

	renamed_seqs = []
	seq_code_dict = {}

	count = 1

	seq_name_start = 'Contig'

	if args.genbank == True:
		for seq_rec in inFasta:
			seq_code_dict.setdefault(seq_rec.id,[]).append(seq_rec.id.split('_')[-1].split('.')[0]+'_Contig_'+str(count)+'_Len'+str(len(seq_rec.seq)))
			seq_code_dict.setdefault(seq_rec.id,[]).append(str(seq_rec.seq).upper())
			renamed_seqs.append('>'+seq_rec.id.split('_')[-1].split('.')[0]+'_Contig_'+str(count)+'_Len'+str(len(seq_rec.seq))+'\n'+str(seq_rec.seq).upper())
			count += 1
	elif args.spades == True:
		for seq_rec in inFasta:
			seq_code_dict.setdefault(seq_rec.description,[]).append(seq_name_start+'_'+str(count)+'_Len'+str(len(seq_rec.seq))+'_Cov'+str(int(round(float(seq_rec.description.split('_')[-3])))))
			seq_code_dict.setdefault(seq_rec.description,[]).append(seq_rec.description.split('_')[5])
			seq_code_dict.setdefault(seq_rec.description,[]).append(str(seq_rec.seq).upper())
			renamed_seqs.append('>'+seq_name_start+'_'+str(count)+'_Len'+str(len(seq_rec.seq))+'_Cov'+str(int(round(float(seq_rec.description.split('_')[-3]))))+'\n'+str(seq_rec.seq).upper())
			count += 1
	else:
		for seq_rec in inFasta:
			seq_code_dict.setdefault(seq_rec.description,[]).append(seq_name_start+'_'+str(count)+'_Len'+str(len(seq_rec.seq)))
			seq_code_dict.setdefault(seq_rec.description,[]).append(str(seq_rec.seq).upper())
			renamed_seqs.append('>'+seq_name_start+'_'+str(count)+'_Len'+str(len(seq_rec.seq))+'\n'+str(seq_rec.seq).upper())
			count += 1
	
	print (color.BOLD+'\n\nThere are '+color.RED+str(len(renamed_seqs))+' contigs > '+str(args.minLen)\
	+color.END+color.BOLD+' in '+color.DARKCYAN+args.input_file.split('/')[-1]+color.END)
	
	with open(home_folder + args.output_file.split('/')[-1] + '.' + str(args.minLen)+'bp.fasta','w+') as w:
		for seq in renamed_seqs:
			w.write(seq+'\n')
	
	if args.spades != True:
		with open(home_folder + args.output_file.split('/')[-1] + '.' + str(args.minLen) + 'bp.SeqCodes.tsv','w+') as w:
			w.write('Original Name\tNew Name\tSeq Length\t Seq GC\n')
			for k, v in seq_code_dict.items():
				w.write(k+'\t'+v[0]+'\t'+str(len(v[1]))+'\t'+str(GC(v[1]))+'\n')
	else:
		with open(home_folder + args.output_file.split('/')[-1] + '.' + str(args.minLen) + 'bp.SeqCodes.tsv','w+') as w:
			w.write('Original Name\tNew Name\tSeq Length\tSeq GC\tSeq Coverage\n')
			for k, v in seq_code_dict.items():
				w.write(k+'\t'+v[0]+'\t'+str(len(v[2]))+'\t'+str(GC(v[2]))+'\t'+str(v[1])+'\n')

	
###########################################################################################
###-------------------------- Cleans Up the PostAssembly Folder ------------------------###
###########################################################################################

def clean_up(args):
		
	os.system('cp ' + args.input_file + ' ' + args.output_file + '/OriginalFasta/' + args.input_file.split('/')[-1].replace('.fasta', '.Original.fasta'))
	
	os.system('cp ' + args.output_file + '/SizeFiltered/' + args.output_file.split('/')[-1] + '.' + str(args.minLen)+'bp.fasta ' + '/'.join(args.output_file.split('/')[:-1]) + '/XlaneBleeding/')


###########################################################################################
###-------------------------------- Next Script Message --------------------------------###
###########################################################################################

def next_script(args):

	print (color.BOLD+'\n\nLook for '+color.DARKCYAN+args.output_file+'.'+str(args.minLen)+\
	'bp.fasta'+color.END+color.BOLD+'\n\n')

	print ('Next Script is: '+color.GREEN+'2_Auto_rRNA_BvE.py'+color.END\
	+color.BOLD+'\n(Alternatively'+color.GREEN+' 2a_remove_rRNA.py followed by 2b_remove_Bact.py'\
	+color.END+color.BOLD+')\n\n'+ color.END)


##########################################################################################
###--------------- Checks Command Line Arguments and Calls on Functions ---------------###
##########################################################################################

def main():

	args = check_args()

	prep_folders(args)
	
	temp = rename_Transcriptome(args)
	
	clean_up(args)
	
	next_script(args)
	
main()
	
