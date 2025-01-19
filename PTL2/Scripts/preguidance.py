# Last updated Nov 2023
# Authors: Auden Cote-L'Heureux, Mario Ceron-Romero, Godwin Ani

# This script is only run when --start = unaligned. This typically means that a user
# is inputting ReadyToGo files as output by EukPhylo part 1. The script contains two optional
# filters. One filter aims to remove sequences outside silent-site GC content ranges set by 
# the user, and relies on the output of the utility script ‘GC_Identifier_v1.0.py.’ See the manual
# for details on using this filter. Sequence filtration by composition is set using the --og_prefix
# parameter; If no OG prefix is input, all sequences will be taken. The second optional filter
# is intended to remove highly similar sequences that might represent ingroup paralogs or 
# redundant transcripts. For each taxon, and for OG in that taxon, all sequences with that OG are 
# similarity-searched against the longest sequence with the OG, and sequences that hit this longest 
# sequence with a percent identity greater --sim_cutoff (usually very high, e.g. 99%) are removed. 
# If users are concerned about only a particular taxon or set of taxa having highly redundant sequences 
# (e.g. human genome), or want to remove highly similar sequences from all but a focal group of taxa, 
# they can input a list of taxa on which the similarity filter is to be exclusively applied (--sim_taxa).

# One other optional filter is a simple blacklist of sequences (--blacklist). Any sequences with IDs that are in this 
# (text) file will be removed from the output Pre-Guidance files.

# The script also reads in the --gf_list (to select a certain list of OGs from input ReadyToGo files) and
# --taxon_list (a list of all taxon names corresponding to the first ten digits of the ReadyToGo files to use)
# and reorganizes the sequence data by OG rather than by taxa. It does this whether the above filters are used
# or not.

#Dependencies
import os, sys, re
from Bio import SeqIO

#This function is called ONLY in eukphylo.py.
def run(params):

	#Reading in the list of gene families to use (--gf_list)
	try:
		ogs = list(dict.fromkeys([line.strip() for line in open(params.gf_list)]))
	except (FileNotFoundError, TypeError) as e:
		print('\nERROR: Unable to read GF list file. Please make sure that the path is correct and that the file is formatted correctly.\n\n' + str(e) + '\n')
		exit()

	#Reading in the list of taxa to use (--taxon_list)
	try:
		taxa = list(dict.fromkeys([line.strip() for line in open(params.taxon_list)]))
	except (FileNotFoundError, TypeError) as e:
		print('\nERROR: Unable to read taxon list file. Please make sure that the path is correct and that the file is formatted correctly.\n\n' + str(e) + '\n')
		exit()

	#Reading in the list of taxa to which to apply the similarity filter
	if params.sim_taxa != None:
		try:
			sim_taxa = list(dict.fromkeys([line.strip() for line in open(params.sim_taxa)]))
		except (FileNotFoundError, TypeError) as e:
			print('\nERROR: Unable to read similarity taxa list file. Please make sure that the path is correct and that the file is formatted correctly.\n\n' + str(e) + '\n')
			exit()
	else:
		sim_taxa = 'all'

	#Reading in any black-listed sequences.
	if params.blacklist != None:
		try:
			blacklist_seqs = list(dict.fromkeys([line.strip() for line in open(params.blacklist)]))
		except (FileNotFoundError, TypeError) as e:
			print('\nERROR: Unable to read blacklist file. Please make sure that the path is correct and that the file is formatted correctly.\n\n' + str(e) + '\n')
			exit()
	else:
		blacklist_seqs = []

	#Looking for input data
	if not os.path.isdir(params.data):
		print('\nInput amino-acid data files not found. Please make sure that the given path (--data) is correct.\n')

	aa_files = [f for f in os.listdir(params.data) if f[:10] in taxa if f.endswith('.faa') or f.endswith('.fa') or f.endswith('.fasta')]

	missing_taxa = [tax for tax in taxa if tax not in [f[:10] for f in aa_files]]
	if(len(missing_taxa) > 0):
		print('\nWARNING: The following taxa in the taxon list are missing amino-acid files in ' + params.data + ':\n' + '\n'.join(['\t' + t for t in missing_taxa]) + '\n')

	os.mkdir(params.output + '/Output/Intermediate/SF_Diamond')
	
	removed_file = open(params.output + '/Output/Pre-Guidance/SimFilter_removed.txt', 'w')

	#Applying similarity filter to each OG and taxon.
	for og in ogs:
		print('\nProcessing ' + og + '\n')
		with open(params.output + '/Output/Pre-Guidance/' + og + '_preguidance.fasta', 'w') as preguidance_file:
			for taxon_file in aa_files:
				recs = []
				#Sorting the records by length
				for rec in sorted([rec for rec in SeqIO.parse(params.data + '/' + taxon_file, 'fasta') if rec.id[-10:] == og and rec.id not in blacklist_seqs and rec.id[-10:].startswith(params.og_identifier)], key=lambda x: -len(x.seq)):
					if(rec.id == rec.description):
						recs.append(rec)
					else:
						print('\n\tThe sequence ID ' + rec.description + ' is invalid. Please make sure that sequence IDs contain no spaces, tabs, etc. This sequence is being excluded.\n')

				#Getting the list of taxa to apply similarity filter to
				if sim_taxa == 'all':
					use_taxon = True
				else:
					if taxon_file[:10] in sim_taxa:
						use_taxon = True
					else:
						use_taxon = False

				masters = []; removed = 0; flag = 0; cycle = 0
				if params.similarity_filter and use_taxon:
					if len(recs) > 1:
						while flag == 0:
							#Creating output files to use in similarity searching
							master_file_name = params.output + '/Output/Intermediate/SF_Diamond/' + og + '_' + taxon_file[:10] + '_master_' + str(cycle)
							query_file_name = params.output + '/Output/Intermediate/SF_Diamond/' + og + '_' + taxon_file[:10] + '_queries_' + str(cycle) + '.fasta'
							diamond_out_name = params.output + '/Output/Intermediate/SF_Diamond/' + og + '_' + taxon_file[:10] + '_diamond_results_' + str(cycle) + '.tsv'

							#Writing out the master (longest) sequence in the OG/taxon
							open(master_file_name + '.faa', 'w').write('>' + recs[0].id + '\n' + str(recs[0].seq) + '\n\n')
							masters.append(recs[0])

							#Writing out all other (query) sequences
							with open(query_file_name, 'w') as queries:
								for rec in recs[1:]:
									queries.write('>' + rec.id + '\n' + str(rec.seq) + '\n\n')

							#Similarity searching all query sequences against the master sequence
							os.system('diamond makedb --in ' + master_file_name + '.faa -d ' + master_file_name)
							os.system('diamond blastp -d ' + master_file_name + '.dmnd -q ' + query_file_name + ' --outfmt 6 -o ' + diamond_out_name)

							#Reading the result
							diamond_out = open(diamond_out_name).readlines()
							recs_to_remove = []
							for line in diamond_out:
								line = line.strip().split('\t')

								#If a sequence hits above the --sim_cutoff identity filter, remove it
								if float(line[2])/100 >= params.sim_cutoff:
									recs_to_remove.append(line[0]); removed =+ 1

							if len([rec for rec in recs[1:] if rec.id not in recs_to_remove]) < 2:
								recs = [rec for rec in recs[1:] if rec.id not in recs_to_remove]
								flag = 1
							else:
								recs = [rec for rec in recs[1:] if rec.id not in recs_to_remove]
								cycle += 1
							
							for item in recs_to_remove:
								removed_file.write(f"{item}\n")

					print('\n\t' + str(removed) + ' sequence(s) removed by the similarity filter (' + str(cycle + 1) + ' iterations) from ' + taxon_file[:10] + '\n')
				
				#Write out the final Pre-Guidance file.
				for rec in recs + masters:
					preguidance_file.write('>' + rec.id + '\n' + str(rec.seq) + '\n\n')
    
	removed_file.close()

	#Remove intermediate files if not specified that they should be kept (--keep_temp)
	if(not params.keep_temp):
		os.system('rm -r ' + params.output + '/Output/Intermediate/SF_Diamond')

