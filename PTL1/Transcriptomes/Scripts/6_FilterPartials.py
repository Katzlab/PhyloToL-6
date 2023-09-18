#!/usr/bin/env python3.5

##__Updated__: 2023-09-18
##__Author__: Xyrus Maurer-Alcala; maurerax@gmail.com; xyrus.maurer-alcala@izb.unibe.ch
##__Usage__: python 6_FilterPartials.py --help


##################################################################################################
## This script is intended to remove incomplete transcripts that have a more complete mate		##
##																								##
## Prior to running this script, ensure the following:											##
##																								##
## 1. You have assembled your transcriptome and COPIED the 'assembly' file 						##
##    (contigs.fasta, or scaffolds.fasta) to the PostAssembly Folder							##
## 2. Removed small sequences (usually sequences < 300bp) with ContigFilterPlusStats.py			##
## 3. Removed SSU/LSU sequences from your Fasta File											##
## 4. Classified your sequences as Strongly Prokaryotic/Eukaryotic or Undetermined				##
## 5. Classified the Non-Strongly Prokaryotic sequences into OGs 								##
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
###------------------- Uses Diamond to perform Self-vs-Self "BLAST" -------------------###
##########################################################################################

def self_blast(args, diamond_path):

	cat_folder = args.all_output_folder + args.file_prefix + '/Original/Concatenated/'

	diamond_makedb = diamond_path + ' makedb --in ' + cat_folder + args.file_prefix + '.AA.Concatenated.fasta -d ' + cat_folder + args.file_prefix + '.AA.Concatenated'

	diamond_self = diamond_path + ' blastp -q ' + cat_folder + args.file_prefix + '.AA.Concatenated.fasta -d ' + cat_folder + args.file_prefix + '.AA.Concatenated --strand plus --no-self-hits --id '+str(args.identity)+\
	' --query-cover 0.7 --evalue 1e-15 --threads 60 --outfmt 6 -o ' + cat_folder + 'SpreadSheets/' + args.file_prefix + '.Concatenated.Self.'+str(args.id_print)+'ID.tsv'

	print (color.BOLD+'\n\nBinning ALL '+color.ORANGE+'Nucleotide ORFs'+color.END+color.BOLD\
	+' for '+color.GREEN+args.file_prefix+color.END+color.BOLD+' at '+args.id_print\
	+'% identity.\n\n'+color.END)

	os.system(diamond_makedb)
	os.system(diamond_self)

	return cat_folder+'SpreadSheets/'+args.file_prefix+'.Concatenated.Self.'+str(args.id_print)+'ID.tsv'


##########################################################################################
###------------------- Uses USearch to perform Self-vs-Self "BLAST" -------------------###
##########################################################################################

def check_Self_vs_Self(tsv_file):

	evaluation = ''

	tsv_in = [i for i in open(tsv_file).read().split('\n') if i != '']

	if len(tsv_in) == 0:
		evaluation = 'empty'
		with open(tsv_file,'w+') as w:
			w.write('No Self-vs-Self hits were found')
	else:
		evaluation = 'continue'

	return evaluation



##########################################################################################
###-------------------- Removes Nearly Identical ORFs from Data Set -------------------###
##########################################################################################

def filter_NTD_data(args, OGLenDB):

	cat_folder = args.all_output_folder + args.file_prefix + '/Original/Concatenated/'
	proc_folder = args.all_output_folder + args.file_prefix + '/Processed/'

	##########################################
	## Set-up Useful Lists and Dictionaries ##
	##########################################

	nuc_Above98_hit = {}
	seqs_to_toss = []
	prepped_NTD = []
	prepped_AA = []

	nuc_tsv_100 = 0

	orf_lens = { rec.id : len(str(rec.seq)) for file in args.file_listNTD for rec in SeqIO.parse(file, 'fasta') }

	replicates = ''

	if len(args.file_listNTD) > 1:
		replicates = 'yes'
	else:
		replicates = 'nope'

	print (color.BOLD+'\n\nRemoving Partial '+color.PURPLE+'ORFs'+color.END+color.BOLD+\
	' with >'+args.id_print+'% Nucleotide Identity over >70% of\ntheir length when '\
	'compared to more complete '+color.PURPLE+'ORFs '+color.END+color.BOLD+'from: '\
	+color.CYAN+args.file_prefix+'\n\n'+color.END)

	#####################################################################
	## Self-v-self BLAST Output Parsing - first checks for Seq-length! ##
	#####################################################################

	nuc_tsv_raw = [i.rstrip('\n') for i in open(cat_folder+'SpreadSheets/'+args.file_prefix\
	+'.Concatenated.Self.'+str(args.id_print)+'ID.tsv').readlines() if i != '\n']

	too_long = 0
	for line in nuc_tsv_raw:
		og_number = re.split('OG.{1}_', line)[-1][:6]
		og_prefix = line.split(og_number)[0][-4:]
		og = og_prefix + og_number

		if og in OGLenDB.keys():
			if orf_lens['_'.join(line.split('\t')[1].split('_')[1:])] > 4.5*OGLenDB[og] or orf_lens['_'.join(line.split('\t')[1].split('_')[1:])] < 1.5*OGLenDB[og]:
				seqs_to_toss.append(line.split('\t')[1])
				too_long += 1

	nuc_tsv = [line for line in nuc_tsv_raw if line.split('\t')[1] not in seqs_to_toss]

	if len(nuc_tsv) > 0:
		if 'Cov' in nuc_tsv[0].split('\t')[0].split('_')[-3]:
			nuc_tsv.sort(key=lambda x: (-int(x.split('\t')[1].split('Len')[-1].split('_')[0]),-int(x.split('\t')[1].split('Cov')[-1].split('_')[0])))
		else:
			nuc_tsv.sort(key=lambda x: -int(x.split('\t')[1].split('Len')[-1].split('_')[0]))

	for line in nuc_tsv:
		if line.split('\t')[1] not in seqs_to_toss:
			nuc_Above98_hit.setdefault(line.split('\t')[1],[]).append(line.split('\t')[0])
			seqs_to_toss.append(line.split('\t')[0])
			if line.split('\t')[2] == '100.0':
				nuc_tsv_100 += 1

	seqs_to_toss = list(set(seqs_to_toss))
	inFasta_NTD_rawLen = [i for i in SeqIO.parse(cat_folder+args.file_prefix+'.NTD.Concatenated.fasta', 'fasta') if i.description]
	inFasta_NTD = [i for i in inFasta_NTD_rawLen if i.description not in seqs_to_toss]
	inFasta_AA = [i for i in SeqIO.parse(cat_folder+args.file_prefix+'.AA.Concatenated.fasta','fasta') if i.description not in seqs_to_toss]

	if replicates != '':
		for i in inFasta_NTD:
			if i.description not in nuc_Above98_hit.keys():
				prepped_NTD.append('>'+'_'.join(i.description.split('_')[1:])+'_Trans1\n'+str(i.seq))
			else:
				Rep_Num = str(len(set([i.description.split('_')[0]]+[j.split('_')[0] for j in nuc_Above98_hit[i.description]])))
				prepped_NTD.append('>'+'_'.join(i.description.split('_')[1:])+'_Trans'+Rep_Num+'\n'+str(i.seq))
		for i in inFasta_AA:
			if i.description not in nuc_Above98_hit.keys():
				prepped_AA.append('>'+'_'.join(i.description.split('_')[1:])+'_Trans1\n'+str(i.seq).replace('*','X'))
			else:
				Rep_Num = str(len(set([i.description.split('_')[0]]+[j.split('_')[0] for j in nuc_Above98_hit[i.description]])))
				prepped_AA.append('>'+'_'.join(i.description.split('_')[1:])+'_Trans'+Rep_Num+'\n'+str(i.seq).replace('*','X'))
	else:
		for i in inFasta_NTD:
			if i.description not in nuc_Above98_hit.keys():
				prepped_NTD.append('>'+i.description+'\n'+str(i.seq))
			else:
				prepped_NTD.append('>'+i.description+'\n'+str(i.seq))
		for i in inFasta_AA:
			if i.description not in nuc_Above98_hit.keys():
				prepped_AA.append('>'+i.description+'\n'+str(i.seq).replace('*','X'))
			else:
				prepped_AA.append('>'+i.description+'\n'+str(i.seq).replace('*','X'))

	with open(args.all_output_folder + args.file_prefix + '/'+args.file_prefix+'_SeqPairsAbove98.txt','w+') as w:
		for k, v in nuc_Above98_hit.items():
			w.write(k+'\t'+'\t'.join(v)+'\n')

	############################################################################################
	## Check for abnormally short and long sequences for the taxon for every Gene Family (OG) ##
	##																													   ##	
	##	IMPORTANT: As of September 2023, ACL and LAK think that the following filter, which	   ##
	## filters ORFs by length relative to the length of other ORFs from the taxon, may be     ##
	## worth revising/removing in the next version of PhyloToL Part 1.							   ##
	##																													   ##
	############################################################################################

	print (color.BOLD+'Removing Abnormally Short (30% length) OR Long (300% length)'\
	+color.PURPLE+' ORFs'+color.END+color.BOLD+'\ncompared to typical '+color.ORANGE+'Gene '\
	'Family '+color.END+color.BOLD+'member length for: '+color.CYAN+args.file_prefix+'\n\n'+color.END)

	self_OGLenDB={} ##
	seqs_to_toss = [] ##
	too_long = too_short = 0 ##

	for i in prepped_NTD:
		og_number = re.split('OG.{1}_', i.split('\n')[0])[-1][:6]
		og_prefix = i.split('\n')[0].split(og_number)[0][-4:]
		og = og_prefix + og_number

		self_OGLenDB.setdefault(og,[]).append(len(i.split('\n')[-1]))

	good_NTD_names = []
	for i in prepped_NTD:
		og_number = re.split('OG.{1}_', i.split('\n')[0])[-1][:6]
		og_prefix = i.split('\n')[0].split(og_number)[0][-4:]
		og = og_prefix + og_number

		if (0.3*sum(self_OGLenDB[og])/float(len(self_OGLenDB[og]))) <= len(i.split('\n')[-1]) <= (3*sum(self_OGLenDB[og])/float(len(self_OGLenDB[og]))):
			good_NTD_names.append(i.split('\n')[0])

	good_NTD_seqs = [i for i in prepped_NTD if i.split('\n')[0] in good_NTD_names]
	good_AA_seqs = [i for i in prepped_AA if i.split('\n')[0] in good_NTD_names]

	too_short = len(prepped_NTD) - len(good_NTD_names)

	####################################################################
	## Finalized Outputs are Summarized and Written Out to New Fastas ##
	####################################################################

	print (color.BOLD+'There were '+color.CYAN+str(len(inFasta_NTD_rawLen))+color.END+color.BOLD\
	+color.PURPLE+' ORFs '+color.END+color.BOLD+'originally, with '+color.ORANGE+\
	str(nuc_tsv_100)+color.END+color.BOLD+' Partial '+color.PURPLE+'ORFs'+color.END+\
	color.BOLD+' that\nwere '+color.RED+'100% Identical'+color.END+color.BOLD+' to larger'\
	+color.PURPLE+' ORFs.\n\n'+color.END)

	print(color.BOLD+'Of the '+color.CYAN+str(len(inFasta_NTD_rawLen))+color.END+color.BOLD\
	+' original'+color.PURPLE+' ORFs'+color.END+color.BOLD+', '+color.ORANGE+str(len(set(seqs_to_toss)))+\
	color.END+color.BOLD+' are '+color.PURPLE+'Partial ORFs '+color.END+color.BOLD+'(e.g. '+\
	color.RED+'> '+args.id_print+'%'+color.END+color.BOLD+'\nNUCLEOTIDE identity) to larger'\
	+color.PURPLE+' ORFs'+color.END+color.BOLD+' with '+color.ORANGE+str(too_short+too_long)\
	+color.END+color.BOLD+' additional'+color.PURPLE+' ORFs\n'+color.END+color.BOLD+'that were either '+\
	color.RED+'TOO LONG or SHORT.\n\n'+color.END)

	print (color.BOLD+'Overall, there are '+color.GREEN+str(len(good_NTD_seqs))+' Unique ORFs'\
	+color.END+color.BOLD+' for '+color.CYAN+args.file_prefix+'\n'+color.END)

	with open(proc_folder+args.file_prefix+'_Filtered.Final.NTD.ORF.fasta','w+') as w:
		for i in good_NTD_seqs:
			w.write(i+'\n')
	with open(proc_folder+args.file_prefix+'_Filtered.Final.AA.ORF.fasta','w+') as x:
		for i in good_AA_seqs:
			x.write(i+'\n')

	return good_NTD_names


##########################################################################################
###------------------- Updates SpreadSheet with Update Sequence Names -----------------###
##########################################################################################

def update_tsv(args, NTD_list_names):

	cat_folder = args.all_output_folder + args.file_prefix + '/Original/Concatenated/SpreadSheets/'
	proc_folder = args.all_output_folder + args.file_prefix + '/Processed/'

	inTSV = {'_'.join(i.split('\t')[0].split('_')[1:]):'\t'.join(i.split('\t')[1:]) for i in open(cat_folder+\
	args.file_prefix+'_Concatenated.allOGCleanresults.tsv').readlines() if i != '\n'}

	Updated_inTSV = [i.strip('>')+'\t'+inTSV[i.split('_Trans')[0].strip('>')] for i in NTD_list_names]

	with open(proc_folder+'/SpreadSheets/'+args.file_prefix+'_Filtered.Final.allOGCleanresults.tsv','w+') as w:
		for line in Updated_inTSV:
			w.write(line+'\n')


############################################################################################
##																													   ##	
##	IMPORTANT: As of September 2023, ACL and LAK think that the following filter, which	   ##
## filters ORFs by length that do NOT hit other ORFs from the taxon by self-BLAST, may be ##
## worth revising/removing in the next version of PhyloToL Part 1.							   ##
##																													   ##
############################################################################################


def no_partials_present(args, OGLenDB):

	print (color.BOLD+color.RED+'\n\nWarning:'+color.END+color.BOLD+' No partial sequences'\
	' were found with > '+str(args.id_print)+'% nucleotide identity.\n\nThe data will still be '\
	'checked for ORFs that are unexpectedly '+color.ORANGE+'Short'+color.END+color.BOLD+' or'\
	+color.ORANGE+' Long.\n\n'+color.END)

	cat_folder = args.all_output_folder + args.file_prefix + '/Original/Concatenated/'
	proc_folder = args.all_output_folder + args.file_prefix + '/Processed/'

	NTD_file = cat_folder+args.file_prefix+'.NTD.Concatenated.fasta'
	AA_file = cat_folder+args.file_prefix+'.AA.Concatenated.fasta'
	TSV_file = cat_folder+'/SpreadSheets/'+args.file_prefix+'_Concatenated.allOGCleanresults.tsv'

	self_OGLenDB = {}
	seqs_to_toss = []
	too_long, too_short = 0, 0

    ## Small changes in this section for Auden (ought to work now)
    ## Lists -> Dictionaries and some data curation steps

	inFasta = {i.description:str(i.seq) for i in SeqIO.parse(NTD_file,'fasta')}

	for k,v in inFasta.items():
		og_number = re.split('OG.{1}_', k)[-1][:6]
		og_prefix = k.split(og_number)[0][-4:]
		og = og_prefix + og_number

		if len(v) >= 4.5*OGLenDB[og]:
			seqs_to_toss.append(k)
			too_long+= 1

	prepped_NTD = [i for i in inFasta if i not in seqs_to_toss]

	print (color.BOLD+'Removing Abnormally Short (30% length) OR Long (300% length)'\
	+color.PURPLE+' ORFs'+color.END+color.BOLD+'\ncompared to typical '+color.ORANGE+'Gene '\
	'Family '+color.END+color.BOLD+'member length for: '+color.CYAN+args.file_prefix+'\n\n'+color.END)

    ## toss those sequences from the sequence dictonary (less headache)
	for crap_seq in seqs_to_toss:
		del inFasta[crap_seq]

	for k, v in inFasta.items():
		og_number = re.split('OG.{1}_', k)[-1][:6]
		og_prefix = k.split(og_number)[0][-4:]
		og = og_prefix + og_number

		self_OGLenDB.setdefault(og,[]).append(len(v))

	self_OGLenDB_Final = {k:sum(v)/len(v) for k, v in self_OGLenDB.items()}

	good_NTD_data = { }
	for k, v in inFasta.items():
		og_number = re.split('OG.{1}_', k)[-1][:6]
		og_prefix = k.split(og_number)[0][-4:]
		og = og_prefix + og_number

		if 0.3*self_OGLenDB_Final[og] <= len(v) <= 3*self_OGLenDB_Final[og]:
			good_NTD_data.update({ k : v })

	good_AA_data = {i.description:str(i.seq) for i in SeqIO.parse(AA_file,'fasta') if i.description in good_NTD_data.keys()}

	good_TSV_data = [i for i in open(cat_folder+'/SpreadSheets/'+args.file_prefix+'_Concatenated.allOGCleanresults.tsv')\
		.read().split('\n') if i != '' and i.split('\t')[0] in good_NTD_data.keys()]

	renamed_TSV_data = [i.split('\t')[0]+'\t'+'\t'.join(i.split('\t')[1:]) for i in good_TSV_data]

	with open(proc_folder+args.file_prefix+'_Filtered.Final.NTD.ORF.fasta','w+') as w:
		for k,v in good_NTD_data.items():
			w.write('>'+k+'\n'+v+'\n')

	with open(proc_folder+args.file_prefix+'_Filtered.Final.AA.ORF.fasta','w+') as x:
		for k, v in good_AA_data.items():
			x.write('>'+k+'\n'+v+'\n')

	with open(proc_folder+'/SpreadSheets/'+args.file_prefix+'_Filtered.Final.allOGCleanresults.tsv','w+') as y:
		y.write('\n'.join(renamed_TSV_data))


############################################################################################
##																													   ##	
##	IMPORTANT: The following function was added in September 2023 by ACL and LAK to        ##
##	replace the pre-Guidance filter that removes sequences <50% and >150% the average      ##
## length of the OG in the Hook. It may be worth revisiting in the next version of        ##
## PhyloToL part 1.							   																##
##																													   ##
############################################################################################


def filter_all_by_hook(args, OGLenDB):

	proc_folder = args.all_output_folder + args.file_prefix + '/Processed/'

	remaining_NTD_seqs = { rec.id : str(rec.seq) for rec in SeqIO.parse(proc_folder+args.file_prefix+'_Filtered.Final.NTD.ORF.fasta', 'fasta') }
	remaining_AA_seqs = { rec.id : str(rec.seq) for rec in SeqIO.parse(proc_folder+args.file_prefix+'_Filtered.Final.AA.ORF.fasta', 'fasta') }

	short_from_translation = []
	if os.path.isfile('UnexpexctedShortStuffBlameXyrus.txt'):
		for line in open('UnexpexctedShortStuffBlameXyrus.txt'):
			short_from_translation.append(line.strip())

	good_NTD_seqs = []; good_AA_seqs = []
	for rec in remaining_NTD_seqs:
		og_number = re.split('OG.{1}_', rec)[-1][:6]
		og_prefix = rec.split(og_number)[0][-4:]
		og = og_prefix + og_number

		if len(remaining_AA_seqs[rec]) <= 1.5*OGLenDB[og] and len(remaining_AA_seqs[rec]) >= 0.33*OGLenDB[og] and rec.split('_Trans')[0] not in short_from_translation:
			good_NTD_seqs.append((rec, remaining_NTD_seqs[rec]))
			good_AA_seqs.append((rec, remaining_AA_seqs[rec]))

	with open(proc_folder+args.file_prefix+'_Filtered.Final.NTD.ORF.fasta','w') as w:
		for rec in good_NTD_seqs:
			w.write('>' + rec[0] + '\n' + rec[1] + '\n\n')

	with open(proc_folder+args.file_prefix+'_Filtered.Final.AA.ORF.fasta','w') as x:
		for rec in good_AA_seqs:
			x.write('>' + rec[0] + '\n' + rec[1] + '\n\n')


##########################################################################################
###--------------------- Cleans up the Folder and Moves Final Files -------------------###
##########################################################################################

def clean_up(args):

	for i in args.file_listNTD:
		os.system('mv ' + i + ' ' + args.all_output_folder + args.file_prefix + '/Original/')
		os.system('mv ' + i.replace('NTD.ORF.fasta','AA.ORF.fasta') + ' ' + args.all_output_folder + args.file_prefix + '/Original/')
		os.system('mv ' + i.split('named')[0]+'named*allOGCleanresults.tsv ' + args.all_output_folder + args.file_prefix + '/Original/SpreadSheets/')


###########################################################################################
###-------------------------------- Next Script Message --------------------------------###
###########################################################################################

def next_script():

	print(color.BOLD+'\nNext Script is: '+color.GREEN+'6b_update_cov_post_removepartials.py\n\n'+color.END)


##########################################################################################
###------------------- Checks Command Line Arguments and Calls Steps ------------------###
##########################################################################################

def main():

	diamond_path = check_diamond_path()

	args = check_args()

	prep_folders(args)

	merge_relevant_data(args)

	self_BLAST_out = self_blast(args, diamond_path)

	evaluation = check_Self_vs_Self(self_BLAST_out)

	OGLenDB = {}
	for rec in SeqIO.parse(args.hook_fasta, 'fasta'):
		if rec.id[-10:] not in OGLenDB:
			OGLenDB.update({ rec.id[-10:] : [] })

		OGLenDB[rec.id[-10:]].append(len(str(rec.seq)))

	for og in OGLenDB:
		OGLenDB[og] = mean(OGLenDB[og])

	if evaluation != 'empty':
		NTD_names = filter_NTD_data(args, OGLenDB)
		update_tsv(args, NTD_names)
	else:
		no_partials_present(args, OGLenDB)

	filter_all_by_hook(args, OGLenDB)

	clean_up(args)

	next_script()

main()
