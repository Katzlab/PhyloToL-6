# Last updated Sept 2017
# Author: Xyrus Maurer-Alcala

# This script takes in CDS as formatted by GenBank for genome assemblies,
# removes input CDS that are shorter than 30bp, renames CDSs for ease of processing
# in later steps of the pipeline, and creates the general output folder structure.
# Input CDS should be named as Op_me_Hsap_GenBankCDS.fasta, with Op_me_Hsap replaced with
# a unique 10-digit sample identifier for each input file. This script is intended to be
# run as part of the PhyloToL 6 Part 1 pipeline using the script wrapper.py.

from Bio import SeqIO
from Bio.SeqUtils import GC
import argparse, os, sys, time
from argparse import RawTextHelpFormatter,SUPPRESS


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
	color.BOLD + '\n\nThis script is intended to extract '+color.RED+'Annotated '+\
	color.PURPLE+'ORFS\n'+color.END+color.BOLD+'from a provided Genbank formatted file.'\
	+usage_msg(), usage=SUPPRESS, formatter_class=RawTextHelpFormatter)
	
	required_arg_group = parser.add_argument_group(color.ORANGE+color.BOLD+'Required Options'+color.END)
 
	required_arg_group.add_argument('--input_file','-in', action='store',
	help=color.BOLD+color.GREEN+' Fasta file with CDSs\n'+color.END)

	required_arg_group.add_argument('--output_dir','-o', action='store',
	help=color.BOLD+color.GREEN+' Output directory\n'+color.END)

	optional_arg_group = parser.add_argument_group(color.ORANGE+color.BOLD+'Options'+color.END)

	optional_arg_group.add_argument('--source','-s', action='store', default='GenBank',
	help=color.BOLD+color.GREEN+' Data Source of CDSs (default = "GenBank")\n'+color.END)

	optional_arg_group.add_argument('--list_source','-lsrc', action='store_true',
	help=color.BOLD+color.GREEN+' Lists supported data sources\n'+color.END)

	optional_arg_group.add_argument('-author', action='store_true',
	help=color.BOLD+color.GREEN+' Prints author contact information\n'+color.END)


	if len(sys.argv[1:]) == 0:
		print (parser.description)
		print ('\n')
		sys.exit()

	args = parser.parse_args()
	
	more_info = return_more_info(args)
	if more_info != None:
		print (parser.description)
		print (more_info)
		sys.exit()
	
	args.folder = args.output_dir + '/' + args.input_file.split('/')[-1][:10]
		
	return args
	
	
###########################################################################################
###------------------------------- Script Usage Message --------------------------------###
###########################################################################################

def usage_msg():
	return (color.BOLD+color.RED+'\n\nExample usage:'+color.CYAN+' python 1g_RenameCDS.py'\
	' --input_file ../Stentor_coeruleus.WGS.CDS.Prep/Stentor_coeruleus.WGS.CDS.fasta --source'\
	' GenBank'+color.END)


##########################################################################################
###-------- Storage for LARGE (Annoying) Print Statements for Flagged Options ---------###
##########################################################################################

def return_more_info(args):

	acceptable_sources = ['in-house', 'in-lab', 'GenBank', 'gb', 'NCBI']

	author = (color.BOLD+color.ORANGE+'\n\n\tQuestions/Comments? Email Xyrus (author) at'\
	' maurerax@gmail.com\n\n'+color.END)

	if args.author == True:
		return author

	if args.list_source == True:
		print (color.BOLD+color.RED+'\nThese are the currently supported data sources.\n'+color.END)
		print (color.BOLD+color.ORANGE+'\n'.join(acceptable_sources)+'\n\n'+color.END)
		sys.exit()

	if args.source.lower() not in [i.lower() for i in acceptable_sources]:
		print (color.BOLD+color.RED+'\nUnsupported source was provided.\n\nEnsure that '\
		'you are providing a valid data source (see below).\n'+color.END)
		print (color.BOLD+color.ORANGE+'\n'.join(acceptable_sources)+'\n'+color.END)
		sys.exit()
	
	if args.input_file != None:
		if args.input_file.split('/')[-1] not in os.listdir('/'.join(args.input_file.split('/')[:-1])):
			print (color.BOLD+color.RED+'\nError:'+color.END+color.BOLD+' The provided Fasta file '\
			'('+color.DARKCYAN+args.input_file.split('/')[-1]+color.END+color.BOLD+')\ndoes not'\
			' exist or is incorrectly formatted.\n\nDouble-check then try again!\n\n'+color.END) 
			sys.exit()
	

###########################################################################################
###--------------------------- Does the Inital Folder Prep -----------------------------###
###########################################################################################

def prep_folders(args):

	if os.path.isdir(args.folder) != True:
		os.system('mkdir '+args.folder)
		os.system('cp '+args.input_file+' '+args.folder)
		args.input_file = args.folder+'/'+args.input_file.split('/')[-1]
		
	if os.path.isdir(args.folder+'/Original') != True:
		os.system('mkdir '+args.folder+'/Original')

	os.system('cp '+args.input_file+' '+args.folder+'/Original/')

###########################################################################################
###------------- Renames Protein-Coding CDS Sequences to Standard Format ---------------###
###########################################################################################

def renamed_GenomeCDS(args):
		
	print (color.BOLD+'\n\nPrepping to rename '+color.GREEN+args.input_file.split('/')[-1]+\
	color.END+color.BOLD+"'s CDS sequences"+color.END)
	inFasta = sorted((i for i in SeqIO.parse(args.input_file,'fasta')),key=lambda seq_rec: -len(seq_rec.seq))

	renamed_seqs = []
	seq_code_dict = {}

	count = 1
	for seq_rec in inFasta:
		seq_code_dict.setdefault(seq_rec.description,[]).append('Contig_'+str(count)+'_Len'+str(len(seq_rec.seq)))
		seq_code_dict[seq_rec.description].append(str(seq_rec.seq).upper())
		renamed_seqs.append('>Contig_' + str(count) + '_Len' + str(len(seq_rec.seq)) + '\n' + str(seq_rec.seq).upper())
		count += 1

	## keeps only CDSs that are greater than 30 bp (10 AA --> This is a cut-off in the 
	## phylogenomic pipeline too!)
	renamed_seqs = [i for i in renamed_seqs if len(i.split('\n')[-1]) > 30]
	
	print (color.BOLD+'\n\nFor '+color.DARKCYAN+args.input_file.split('/')[-1]+color.END+\
	color.BOLD+', '+color.RED+str(len(renamed_seqs))+' CDS sequences\n'+color.END+color.BOLD+
	'were renamed while preserving the '+color.ORANGE+args.source+color.END+color.BOLD+' formatting'\
	+color.END+'\n')
	
	with open(args.input_file.replace('.fasta','.Prepped.fasta'),'w+') as w:
		w.write('\n'.join(renamed_seqs))

	with open(args.input_file.split('/')[-1].replace('.fasta','.SeqCodes.tsv'),'w+') as w:
		w.write('Original Name\tNew Name\tSeq Length\t Seq GC\n')
		for k, v in seq_code_dict.items():
			w.write(k+'\t'+v[0]+'\t'+str(len(v[1]))+'\t'+str(GC(v[1]))+'\n')


###########################################################################################
###--------------------- Cleans up the Folder and Moves Final Files --------------------###
###########################################################################################

def clean_up(args):

#	os.system('rm '+args.input_file)
	os.system('mv ' + args.input_file.split('/')[-1].replace('.fasta','.SeqCodes.tsv') + ' ' + args.folder + '/Original/')
	os.system('mv ' + args.input_file + ' ' + args.folder + '/Original/')


###########################################################################################
###-------------------------------- Next Script Message --------------------------------###
###########################################################################################

def next_script(args):

	print (color.BOLD+'\nLook for '+color.DARKCYAN+args.input_file.split('/')[-1].replace('.fasta','.Renamed.fasta')\
	+'.fasta'+color.END+color.BOLD+'\nin the '+color.ORANGE+args.folder.split('/')[-1]+\
	' Folder\n\n'+color.END+color.BOLD)

	print ('Next Script(s) are:\n\n'+color.PURPLE+'2g_GCodeEval.py'+color.END+color.BOLD\
	+' (if Genetic Code is '+color.RED+'Unknown'+color.END+color.BOLD+')\n\nOtherwise:\n\n'+\
	color.PURPLE+'3g_GCodeTranslate.py\n\n'+color.END)
	
	
##########################################################################################
###----------------------------- Calls on Above Functions -----------------------------###
##########################################################################################

def main():

	args = check_args()

	prep_folders(args)
		
	renamed_GenomeCDS(args)
	
	clean_up(args)
	
	next_script(args)
	
main()
