import os, sys, re
from logger import Logger
from Bio import SeqIO



def run(params):

	try:
		ogs = list(dict.fromkeys([line.strip() for line in open(params.gf_list)]))
	except (FileNotFoundError, TypeError) as e:
		Logger.Error('Unable to read GF list file. Please make sure that the path is correct and that the file is formatted correctly.\n\n' + str(e))

	try:
		taxa = list(dict.fromkeys([line.strip() for line in open(params.taxon_list)]))
	except (FileNotFoundError, TypeError) as e:
		Logger.Error('Unable to read taxon list file. Please make sure that the path is correct and that the file is formatted correctly.\n\n' + str(e))

	if not os.path.isdir(params.data):
		Logger.Error(Logger.Error('Input amino-acid data files not found. Please make sure that the given path (--data) is correct.'))

	aa_files = [f for f in os.listdir(params.data) if f[:10] in taxa if f.endswith('.faa') or f.endswith('.fa') or f.endswith('.fasta')]

	missing_taxa = [tax for tax in taxa if tax not in [f[:10] for f in aa_files]]
	if(len(missing_taxa) > 0):
		Logger.Warning('The following taxa in the taxon list are missing amino-acid files in ' + params.data + ':\n' + '\n'.join(['\t' + t for t in missing_taxa]))

	os.mkdir(params.output + '/Output/Temp/OF-SF_Diamond')

	for og in ogs:
		Logger.Message('Processing ' + og)
		with open(params.output + '/Output/Pre-Guidance/' + og + '_preguidance.faa', 'w') as preguidance_file:
			for taxon_file in aa_files:
				recs = []
				for rec in sorted([rec for rec in SeqIO.parse(params.data + '/' + taxon_file, 'fasta') if rec.id[-10:] == og], key=lambda x: -len(x.seq)):
					if(rec.id == rec.description):
						recs.append(rec)
					else:
						Logger.Warning('\tThe sequence ID ' + rec.description + ' is invalid. Please make sure that sequence IDs contain no spaces, tabs, etc. This sequence is being excluded.')

				masters = []; removed = 0; flag = 0; cycle = 0
				if len(recs) > 1:
					while flag == 0:
						master_file_name = params.output + '/Output/Temp/OF-SF_Diamond/' + og + '_' + taxon_file[:10] + '_master_' + str(cycle)
						query_file_name = params.output + '/Output/Temp/OF-SF_Diamond/' + og + '_' + taxon_file[:10] + '_queries_' + str(cycle) + '.faa'
						diamond_out_name = params.output + '/Output/Temp/OF-SF_Diamond/' + og + '_' + taxon_file[:10] + '_diamond_results_' + str(cycle) + '.tsv'

						open(master_file_name + '.faa', 'w').write('>' + recs[0].id + '\n' + str(recs[0].seq) + '\n\n')
						masters.append(recs[0])
						
						with open(query_file_name, 'w') as queries:
							for rec in recs[1:]:
								queries.write('>' + rec.id + '\n' + str(rec.seq) + '\n\n')

						os.system('diamond makedb --in ' + master_file_name + '.faa -d ' + master_file_name)
						os.system('diamond blastp -d ' + master_file_name + '.dmnd -q ' + query_file_name + ' --outfmt 6 -o ' + diamond_out_name)

						diamond_out = open(diamond_out_name).readlines()
						recs_to_remove = []
						for line in diamond_out:
							line = line.strip().split('\t')
							alignment_length = int(line[3]); gaps = int(line[5]); seq = str(line[0]); identity = float(line[2])
							
							if ((alignment_length - gaps) < params.overlap_cutoff * len(recs[0].seq) and cycle == 0) or identity > params.sim_cutoff:
								recs_to_remove.append(seq); removed =+ 1

						if len([rec for rec in recs[1:] if rec.id not in recs_to_remove]) < 2:
							recs = [rec for rec in recs[1:] if rec.id not in recs_to_remove]
							flag = 1
						else:
							recs = [rec for rec in recs[1:] if rec.id not in recs_to_remove]
							cycle += 1

				Logger.Message('\t' + str(removed) + ' sequence(s) removed by the overlap/similarity filters (' + str(cycle + 1) + ' iterations) from ' + taxon_file[:10])
				for rec in recs + masters:
					preguidance_file.write('>' + rec.id + '\n' + str(rec.seq) + '\n\n')

	if(not params.keep_temp):
		os.system('rm -r ' + params.output + '/Output/Temp/OF-SF_Diamond')









































