#!/usr/bin/env python3.5

##__Updated__: 2023-09-27 by Auden Cote-L'Heureux
##__Author__: Xyrus Maurer-Alcala; maurerax@gmail.com; xyrus.maurer-alcala@izb.unibe.ch
##__Usage__: python 6_FilterPartials.py --help


##################################################################################################
## This script is intended to remove incomplete transcripts that have a more complete mate		##
##																								##
## Prior to running this script, ensure the following:											##
##																								##
## 1. You have assembled your transcriptome and COPIED the 'assembly' file 						##
##    (contigs.fasta, or scaffolds.fasta) to the PostAssembly Folder							##
## 2. Removed small sequences (usually sequences < 200bp)			##
## 3. Removed SSU/LSU sequences from your Fasta File											##
## 4. Classified your sequences as Strongly Prokaryotic/Eukaryotic or Undetermined				##
## 5. Classified sequences into OGs 								##
## 6. You either know (or have inferred) the genetic code of the organism						##
## 7. You have translated the sequences and checked for the data in the RemovePartials folder	##
##																								##
## 					E-mail Xyrus (author) for help if needed: maurerax@gmail.com				##
##																								##
##										Next Script(s) to Run: 									##
##						 	  			  7_FinalRename.py										##
##																								##
##################################################################################################

from Bio import SeqIO
from Bio.Seq import Seq
from statistics import mean

from distutils import spawn
import argparse, os, sys, time, re
from argparse import RawTextHelpFormatter,SUPPRESS

from tqdm import tqdm


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
		+' the PATH to the'+color.BLUE+' "diamond" '+color.END+color.BOLD+'executable.\n\n'+color.END)
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
	color.BOLD + '\n\nThis script is intended to '+color.RED+'Identify and Collapse '+color.END\
	+color.BOLD+'partial '+color.PURPLE+'ORFS\n'+color.END+color.BOLD+'present within a '\
	+color.RED+'Given'+color.END+color.BOLD+' transcriptome (or replicate) transcriptome(s)'\
	+usage_msg(), usage=SUPPRESS, formatter_class=RawTextHelpFormatter)

	required_arg_group = parser.add_argument_group(color.ORANGE+color.BOLD+'Required Options'+color.END)

	required_arg_group.add_argument('--file_prefix','-fp', action='store',
	help=color.BOLD+color.GREEN+' File prefix that is unique (or common)\n to the files '\
	'to be processed\n'+color.END)


	optional_arg_group = parser.add_argument_group(color.ORANGE+color.BOLD+'Options'+color.END)

	optional_arg_group.add_argument('--identity','-id', type=float, action='store', default=0.98,
	help=color.BOLD+color.GREEN+' Identity threshold for identifying \n "partials" to larger'\
	' contigs\n (default = 0.98)\n'+color.END)
	optional_arg_group.add_argument('-author', action='store_true',
	help=color.BOLD+color.GREEN+' Prints author contact information\n'+color.END)
	optional_arg_group.add_argument('--hook_fasta','-f', help='Path to the fasta file of the Hook DB in the Databases/db_OG folder')

	if len(sys.argv[1:]) == 0:
		print (parser.description)
		print ('\n')
		sys.exit()

	args = parser.parse_args()

	args.id_print = str(int(float(args.identity)*100))

	args.all_output_folder = '/'.join(args.file_prefix.split('/')[:-1]) + '/'
	args.file_prefix = args.file_prefix.split('/')[-1]

	args.file_listNTD = [args.all_output_folder + i for i in os.listdir(args.all_output_folder) if args.file_prefix in i and i.endswith('NTD.ORF.fasta')]

	args.file_listAA = [args.all_output_folder + i for i in os.listdir(args.all_output_folder) if args.file_prefix in i and i.endswith('AA.ORF.fasta')]

	args.file_listTSV = [args.all_output_folder + i for i in os.listdir(args.all_output_folder) if args.file_prefix in i and i.endswith('results.tsv')]

	quit_eval = return_more_info(args)
	if quit_eval > 0:
		print ('\n')
		sys.exit()

	return args


###########################################################################################
###------------------------------- Script Usage Message --------------------------------###
###########################################################################################

def usage_msg():
	return (color.BOLD+color.RED+'\n\nExample usage:'+color.CYAN+' python 6_RemovePartials.py'\
	' --file_prefix Op_me_Xxma'+color.END)


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

	if args.file_listNTD == []:
		print (color.BOLD+'\n\nNo '+color.ORANGE+'Nucleotide Fasta Files'+color.END+color.BOLD+\
		' found!\n\nCheck that your'+color.GREEN+' File Prefix'+color.END+color.BOLD+\
		'is present in\nthe files of interest')
		valid_arg += 1

	if args.file_listAA == []:
		print (color.BOLD+'\n\nNo '+color.ORANGE+'Protein Fasta Files'+color.END+color.BOLD+\
		' found!\n\nCheck that your'+color.GREEN+' File Prefix'+color.END+color.BOLD+\
		'is present in\nthe files of interest')
		valid_arg += 1

	if args.file_listTSV == []:
		print (color.BOLD+'\n\nNo '+color.ORANGE+'OG-Assignment Spreadsheets'+color.END+color.BOLD+\
		' found!\n\nCheck that your'+color.GREEN+' File Prefix'+color.END+color.BOLD+\
		'is present in\nthe files of interest')
		valid_arg += 1

	if len(args.file_listNTD) == len(args.file_listAA) == len(args.file_listTSV):
		pass
	else:
		print (color.BOLD+color.RED+'\n\nError:'+color.END+color.BOLD+' Unequal numbers of'\
		' input files found.\n\nDouble-check that there are:'+color.CYAN+'SINGLE'+color.END\
		+color.BOLD+' Nucleotide and Protein fasta files and OG-assignment Spreadsheet for'\
		' each transcriptome\n\nThen try once again.'+color.END)
		valid_arg += 1

	return valid_arg


##########################################################################################
###------------------------- Creates Folders For Storing Data -------------------------###
##########################################################################################

def prep_folders(args):

	if os.path.isdir(args.all_output_folder + 'ToRename') != True:
		os.system('mkdir ' + args.all_output_folder + 'ToRename')

	if os.path.isdir(args.all_output_folder + args.file_prefix) != True:
		os.system('mkdir ' + args.all_output_folder + args.file_prefix)

	if os.path.isdir(args.all_output_folder + args.file_prefix + '/Original') != True:
		os.system('mkdir ' + args.all_output_folder + args.file_prefix + '/Original')
		os.system('mkdir ' + args.all_output_folder + args.file_prefix + '/Original/SpreadSheets')
		os.system('mkdir ' + args.all_output_folder + args.file_prefix + '/Original/Concatenated/')
		os.system('mkdir ' + args.all_output_folder + args.file_prefix + '/Original/Concatenated/SpreadSheets')

	if os.path.isdir(args.all_output_folder + args.file_prefix + '/Processed') != True:
		os.system('mkdir ' + args.all_output_folder + args.file_prefix + '/Processed')
		os.system('mkdir ' + args.all_output_folder + args.file_prefix + '/Processed/SpreadSheets')


##########################################################################################
###-------------------- Merges Fasta Files When Replicates Present --------------------###
##########################################################################################

def merge_fasta_replicates(args, type):

	cat_folder = args.all_output_folder + args.file_prefix + '/Original/Concatenated/'

	count = 0
	fasta_to_merge = []

	if type == 'NTD':
		fasta_list = args.file_listNTD
	else:
		fasta_list = args.file_listAA

	for file in fasta_list:
		fasta_to_merge += ['>'+str(count)+'_'+i for i in open(file).read().split('>') if i != '']
		count += 1

	with open(cat_folder+args.file_prefix+'.'+type+'.Concatenated.fasta','w+') as w:
		w.write(''.join(fasta_to_merge))

	time.sleep(.75)


##########################################################################################
###--------------------- Merges TSV Files When Replicates Present ---------------------###
##########################################################################################

def merge_tsv_replicates(args):

	cat_folder = args.all_output_folder + args.file_prefix + '/Original/Concatenated/SpreadSheets/'

	count = 0
	tsv_to_merge = []

	for file in args.file_listTSV:
		tsv_to_merge += [str(count)+'_'+i for i in open(file).read().split('\n') if i != '']
		count += 1

	with open(cat_folder+args.file_prefix+'_Concatenated.allOGCleanresults.tsv','w+') as w:
		w.write('\n'.join(tsv_to_merge))

	time.sleep(.75)


##########################################################################################
###------------------ Calls on the other Merge Functions by Data Type -----------------###
##########################################################################################

def merge_relevant_data(args):

	print (color.BOLD+'\n\nMerging Transcriptome data together.'+color.END)

	merge_fasta_replicates(args, 'NTD')
	merge_fasta_replicates(args, 'AA')
	merge_tsv_replicates(args)


##########################################################################################
###-------------------- Removes Nearly Identical ORFs from Data Set -------------------###
##########################################################################################

def filter_NTD_data(args, OGLenDB):

	cat_folder = args.all_output_folder + args.file_prefix + '/Original/Concatenated/'
	proc_folder = args.all_output_folder + args.file_prefix + '/Processed/'

	print (color.BOLD+'\n\nRemoving Partial '+color.PURPLE+'ORFs'+color.END+color.BOLD+\
	' with >'+args.id_print+'% Nucleotide Identity over >33% of\ntheir length when '\
	'compared to more complete '+color.PURPLE+'ORFs '+color.END+color.BOLD+'from: '\
	+color.CYAN+args.file_prefix+'\n\n'+color.END)

	#Read in the NTD and AA sequences output by script 5 and concatenated across files for the taxon above
	starting_NTD_seqs = { rec.id : str(rec.seq) for rec in SeqIO.parse(cat_folder+args.file_prefix+'.NTD.Concatenated.fasta', 'fasta') }
	starting_AA_seqs = { rec.id : str(rec.seq) for rec in SeqIO.parse(cat_folder+args.file_prefix+'.AA.Concatenated.fasta', 'fasta') }

	#Creating a record of short sequences left over from translation to be removed
	short_from_translation = []
	if os.path.isfile('ShortTranscripts_FromTranslation.txt'):
		for line in open('ShortTranscripts_FromTranslation.txt'):
			short_from_translation.append(line.strip())

	#Remove sequences <33% or >150% the average length of their OG in the Hook database
	good_NTD_seqs = []; good_AA_seqs = []
	for rec in starting_NTD_seqs:
		og_number = re.split('OG.{1}_', rec)[-1][:6]
		og_prefix = rec.split(og_number)[0][-4:]
		og = og_prefix + og_number

		if len(starting_AA_seqs[rec]) <= 1.5*OGLenDB[og] and len(starting_AA_seqs[rec]) >= 0.33*OGLenDB[og] and '-'.join(rec.split('_')[1:]) not in short_from_translation:
			good_NTD_seqs.append((rec, starting_NTD_seqs[rec]))
			good_AA_seqs.append((rec, starting_AA_seqs[rec]))

	#Write out all sequences after removal of ORFs based on comparison to mean Hook OG length
	with open(proc_folder + args.file_prefix + '.NTD.HookLenFiltered_NotPartialFiltered.fasta', 'w') as o:
		for seq in good_NTD_seqs:
			o.write('>' + seq[0] + '\n' + seq[1] + '\n\n')

	with open(proc_folder + args.file_prefix + '.AA.HookLenFiltered_NotPartialFiltered.fasta', 'w') as o:
		for seq in good_AA_seqs:
			o.write('>' + seq[0] + '\n' + seq[1] + '\n\n')

	#BLAST-ing all nucleotide sequences that survived the Hook-relative length filter against each other
	db_cmd = 'makeblastdb -in ' + proc_folder + args.file_prefix + '.NTD.HookLenFiltered_NotPartialFiltered.fasta -dbtype nucl -parse_seqids -out ' + proc_folder + args.file_prefix + '.NTD.HookLenFiltered_NotPartialFiltered'
	blastn_cmd = 'blastn -query ' + proc_folder + args.file_prefix + '.NTD.HookLenFiltered_NotPartialFiltered.fasta -db ' + proc_folder + args.file_prefix + '.NTD.HookLenFiltered_NotPartialFiltered -perc_identity 98 -outfmt "6 std qcovs" -out ' + proc_folder + '/SpreadSheets/All_NTD_SelfBLAST_Results.tsv'

	os.system(db_cmd)
	os.system(blastn_cmd)

	#Creating a record of query, subject pairs for each OG
	pairs_per_og = { }
	for line in open(proc_folder + '/SpreadSheets/All_NTD_SelfBLAST_Results.tsv'):

		#IMPORTANT: This line is where the query coverage threshold is determined, and it might be helpful to adjust this when trying to optimally filter chimeras. By default it is set to 20%
		if line.split('\t')[0] != line.split('\t')[1] and line.split('\t')[0][-10:] == line.split('\t')[1][-10:] and int(line.split('\t')[-1].strip()) > 20:
			if line.split('\t')[0][-10:] not in pairs_per_og:
				pairs_per_og.update({ line.split('\t')[0][-10:] : [] })

			pairs_per_og[line.split('\t')[0][-10:]].append([line.split('\t')[0], line.split('\t')[1]])

	#For each OG, iterating through master sequences by decreasing score (cov*len), removing sequences that hit each master
	partials_to_remove = []; pairs_with_query_removed = []
	for og in pairs_per_og:

		#Sorting the sequences by score (cov*len)
		seqs = sorted(list(dict.fromkeys([seq for pair in pairs_per_og[og] for seq in pair])), key = lambda x : -(int(x.split('Len')[-1].split('_')[0]) * int(x.split('Cov')[-1].split('_')[0])))

		#Iterating through masters and determining which sequences to remove
		for master in seqs:
			if master not in partials_to_remove:
				for pair in pairs_per_og[og]:
					if pair[1] == master:
						partials_to_remove.append(pair[0])
						pairs_with_query_removed.append(pair)

	#Writing out a record of the BLAST hits of all relevant pairs (i.e. when removed sequences hit a master sequence)
	with open(args.all_output_folder + args.file_prefix + '/'+args.file_prefix+'_SeqPairsAbove98.txt','w') as w:
		for line in open(proc_folder + '/SpreadSheets/All_NTD_SelfBLAST_Results.tsv'):
			if [line.split('\t')[0], line.split('\t')[1]] in pairs_with_query_removed:
				w.write(line)

	####################################################################
	## Finalized Outputs are Summarized and Written Out to New Fastas ##
	####################################################################

	# print (color.BOLD+'There were '+color.CYAN+str(len(inFasta_NTD_rawLen))+color.END+color.BOLD\
	# +color.PURPLE+' ORFs '+color.END+color.BOLD+'originally, with '+color.ORANGE+\
	# str(nuc_tsv_100)+color.END+color.BOLD+' Partial '+color.PURPLE+'ORFs'+color.END+\
	# color.BOLD+' that\nwere '+color.RED+'100% Identical'+color.END+color.BOLD+' to larger'\
	# +color.PURPLE+' ORFs.\n\n'+color.END)

	# print(color.BOLD+'Of the '+color.CYAN+str(len(inFasta_NTD_rawLen))+color.END+color.BOLD\
	# +' original'+color.PURPLE+' ORFs'+color.END+color.BOLD+', '+color.ORANGE+str(len(set(seqs_to_toss)))+\
	# color.END+color.BOLD+' are '+color.PURPLE+'Partial ORFs '+color.END+color.BOLD+'(e.g. '+\
	# color.RED+'> '+args.id_print+'%'+color.END+color.BOLD+'\nNUCLEOTIDE identity) to larger'\
	# +color.PURPLE+' ORFs'+color.END+color.BOLD+' with '+color.ORANGE+str(too_short+too_long)\
	# +color.END+color.BOLD+' additional'+color.PURPLE+' ORFs\n'+color.END+color.BOLD+'that were either '+\
	# color.RED+'TOO LONG or SHORT.\n\n'+color.END)

	# print (color.BOLD+'Overall, there are '+color.GREEN+str(len(good_NTD_seqs))+' Unique ORFs'\
	# +color.END+color.BOLD+' for '+color.CYAN+args.file_prefix+'\n'+color.END)

	with open(proc_folder+args.file_prefix+'_Filtered.Final.NTD.ORF.fasta','w+') as w:
		for i in good_NTD_seqs:
			if i[0] not in partials_to_remove:
				w.write('>' + i[0] + '\n' + str(i[1]) + '\n')

	with open(proc_folder+args.file_prefix+'_Filtered.Final.AA.ORF.fasta','w+') as x:
		for i in good_AA_seqs:
			if i[0] not in partials_to_remove:
				x.write('>' + i[0] + '\n' + str(i[1]) + '\n')

	good_seq_names = [i[0] for i in good_NTD_seqs]

	with open(proc_folder + '/SpreadSheets/' + args.file_prefix + '_Filtered.Final.allOGCleanresults.tsv', 'w') as t:
		for line in open(cat_folder + '/SpreadSheets/' + args.file_prefix + '_Concatenated.allOGCleanresults.tsv'):
			if line.split('\t')[0] in good_seq_names and line.split('\t')[0] not in partials_to_remove:
				t.write(line)


##########################################################################################
###--------------------- Cleans up the Folder and Moves Final Files -------------------###
##########################################################################################

def clean_up(args):

	for i in args.file_listNTD:
		os.system('mv ' + i + ' ' + args.all_output_folder + args.file_prefix + '/Original/')
		os.system('mv ' + i.replace('NTD.ORF.fasta','AA.ORF.fasta') + ' ' + args.all_output_folder + args.file_prefix + '/Original/')
		os.system('mv ' + i.split('named')[0]+'named*allOGCleanresults.tsv ' + args.all_output_folder + args.file_prefix + '/Original/SpreadSheets/')

	os.system('cp ' + args.all_output_folder + args.file_prefix + '/Processed/*ORF.fasta ' + args.all_output_folder + '/ToRename/')
	os.system('cp ' + args.all_output_folder + args.file_prefix + '/Processed/SpreadSheets/*allOGCleanresults.tsv ' + args.all_output_folder + '/ToRename/')


###########################################################################################
###-------------------------------- Next Script Message --------------------------------###
###########################################################################################

def next_script():

	print(color.BOLD+'\nNext Script is: '+color.GREEN+'7_FinalizeName.py\n\n'+color.END)


##########################################################################################
###------------------- Checks Command Line Arguments and Calls Steps ------------------###
##########################################################################################

def main():

	diamond_path = check_diamond_path()

	args = check_args()

	prep_folders(args)

	merge_relevant_data(args)

	OGLenDB = {}
	for rec in SeqIO.parse(args.hook_fasta, 'fasta'):
		if rec.id[-10:] not in OGLenDB:
			OGLenDB.update({ rec.id[-10:] : [] })

		OGLenDB[rec.id[-10:]].append(len(str(rec.seq)))

	for og in OGLenDB:
		OGLenDB[og] = mean(OGLenDB[og])

	filter_NTD_data(args, OGLenDB)

	clean_up(args)

	next_script()

main()
