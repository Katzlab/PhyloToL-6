# Last updated Sept 2023
# Authors: Auden Cote-L'Heureux and Mario Ceron-Romero

# This script runs Guidance in an iterative fashion for more both MSA construction 
# and more rigorous homology assessment than what is offered in PhyloToL 6 part 1.
# Guidance runs until the input number of iterations (--guidance_iters, default = 5) 
# has been reached, or until there are no sequences below the sequence score cutoff.
# All sequences below the score cutoff (--seq_cutoff, default = 0.3) are removed at
# each iteration. By default, PhyloToL does not remove residues that fall below the 
# given residue cutoff (--res_cutoff) and columns that fall below the given column 
# cutoff (--col_cutoff, defaults are 0), though this can be turned on by adjusting 
# these parameters. Outputs at this point are found in the “Guidance_NotGapTrimmed” 
# output folder. We then run MSAs through TrimAl to remove all sites in the alignment 
# that are at least 95% gaps (or --gap_trim_cutoff) generating files in the “Guidance” 
# output folder.

# This step is either intended to be run starting with --start = unaligned (but not raw)
# inputs, meaning one amino acid alignment per OG. It can also be run directly after the
# preguidance step. The run() function is called in two places: in phylotol.py generally,
# and in contamination.py if the contamination loop is using Guidance as the re-alignment
# method.

#Dependencies
import os, sys, re
from Bio import SeqIO

#Called in phylotol.py and contamination.py
def run(params):

	if params.start == 'raw' or params.start == 'unaligned':
		#Checking that pre-Guidance has been run or that unaligned files per OG are provided.
		if params.start == 'raw':
			preguidance_path = params.output + '/Output/Pre-Guidance'
		else:
			preguidance_path = params.data

		if not os.path.isdir(preguidance_path):
			print('\nERROR: The path ' + preguidance_path + ' could not be found when trying to locate pre-Guidance (unaligned) files. Make sure that the --start and --data parameters are correct and/or that the pre-Guidance step ran successfully.\n')
			exit()
		
		if len([f for f in os.listdir(preguidance_path) if f.endswith('.fa') or f.endswith('.faa') or f.endswith('.fasta')]) == 0:
			print('\nERROR: No pre-Guidance (unaligned) files could be found at the path ' + preguidance_path + '. Make sure that the --start and --data parameters are correct, that the pre-Guidance step ran successfully, and that the unaligned files are formatted correctly (they must have the file extension .faa, .fa, or .fasta).\n')
			exit()

		#Creating intermedate folders that will later be deleted unless running with --keep_temp
		os.mkdir(params.output + '/Output/Intermediate/Guidance')
		os.mkdir(params.output + '/Output/Intermediate/Guidance/Input')
		os.mkdir(params.output + '/Output/Intermediate/Guidance/Output')

		guidance_input = params.output + '/Output/Intermediate/Guidance/Input/'
		os.system('cp -r ' + preguidance_path + '/* ' + guidance_input)

		guidance_removed_file = open(params.output + '/Output/GuidanceRemovedSeqs.txt', 'w')
		guidance_removed_file.write('Sequence\tScore\n')

		#For each unaligned AA fasta file
		for file in [f for f in os.listdir(guidance_input) if f.endswith('.fa') or f.endswith('.faa') or f.endswith('.fasta')]:
			tax_guidance_outdir = params.output + '/Output/Intermediate/Guidance/Output/' + file.split('.')[0].split('_preguidance')[0]
			os.mkdir(tax_guidance_outdir)

			fail = False
			#For each iteration
			for i in range(params.guidance_iters):
				n_recs = len([r for r in SeqIO.parse(guidance_input + '/' + file, 'fasta')])

				#Guidance can't handle inputs with fewer than 4 sequences
				if n_recs < 4:
					print('\nWARNING: Gene famiily ' + file.split('.')[0].split('_preguidance')[0] + ' contains fewer than 4 sequences after ' + str(i) + ' Guidance iterations, therefore no alignment will be produced for this gene family.\n')
					os.system('rm -rf ' + tax_guidance_outdir)
					if i == 0:
						fail = True
					break

				#Determining MAFFT algorithm based on the number of input sequences
				if n_recs < 200:
					mafft_alg = 'genafpair'
				else:
					mafft_alg = 'auto'

				#Running Guidance (one per OG per iteration)
				os.system('Scripts/guidance.v2.02/www/Guidance/guidance.pl --seqFile ' + guidance_input + '/' + file + ' --msaProgram MAFFT --seqType aa --outDir ' + tax_guidance_outdir + ' --seqCutoff ' + str(params.seq_cutoff) + ' --colCutoff ' + str(params.col_cutoff) + " --outOrder as_input --bootstraps 10 --MSA_Param '\\--" + mafft_alg + " --maxiterate 1000 --thread " + str(params.guidance_threads) + " --bl 62 --anysymbol' > " + params.output + '/Output/Intermediate/Guidance/Output/' + file[:10] + '/log.txt')

				#Checking for a sequence score file; if not available, Guidance failed.
				if os.path.isfile(tax_guidance_outdir + '/MSA.MAFFT.Guidance2_res_pair_seq.scr_with_Names'):
					#All sequences below score cutoff
					seqs_below = len([line for line in open(tax_guidance_outdir + '/MSA.MAFFT.Guidance2_res_pair_seq.scr_with_Names').readlines()[1:-1] if float(line.split()[-1]) < params.seq_cutoff])
					#If fewer than four were above the cutoff, this OG is done iterating.
					if n_recs - seqs_below < 4:
						print('\nWARNING: Gene famiily ' + file.split('.')[0].split('_preguidance')[0] + ' contains fewer than 4 sequences after ' + str(i + 1) + ' Guidance iterations, therefore no alignment will be produced for this gene family.\n')
						os.system('rm -rf ' + tax_guidance_outdir)
						break
					#If all sequences were above the cutoff, this OG is done iterating.
					if seqs_below == 0 or i == params.guidance_iters - 1:
						print('\nGuidance complete after ' + str(i + 1) + ' iterations for gene family ' + file.split('.')[0].split('_preguidance')[0] + '\n')
						break
					#Recording list of sequences removed by Guidance.
					for line in [line for line in open(tax_guidance_outdir + '/MSA.MAFFT.Guidance2_res_pair_seq.scr_with_Names').readlines()[1:-1] if float(line.split()[-1]) < params.seq_cutoff]:
						guidance_removed_file.write(line)
					#Copying over the old file with the new results
					os.system('cp ' + tax_guidance_outdir + '/Seqs.Orig.fas.FIXED.Without_low_SP_Seq.With_Names ' + guidance_input + '/' + file)
					#Cleaning up the intermediate files for the iteration.
					os.system('rm -r ' + tax_guidance_outdir + '/*')
				else:
					fail = True
					break

			#After all iterations, THEN apply residue and column cutoffs
			if not fail:
				#Getting a list of sequences to keep
				seqs2keep = [rec.description for rec in SeqIO.parse(tax_guidance_outdir + '/Seqs.Orig.fas.FIXED.Without_low_SP_Seq.With_Names', 'fasta')]
				orig_seqs = [rec.description for rec in SeqIO.parse(tax_guidance_outdir + '/MSA.MAFFT.aln.With_Names', 'fasta')]
				running_aln = { rec.description : str(rec.seq) for rec in SeqIO.parse(tax_guidance_outdir + '/MSA.MAFFT.aln.With_Names', 'fasta') if rec.description in seqs2keep }

				#Residues that fall below the confidence cutoff (--res_cutoff) are replaced with 'X'
				for site in [(int(line.split()[1]), int(line.split()[0]) - 1) for line in open(tax_guidance_outdir + '/MSA.MAFFT.Guidance2_res_pair_seq.scr').readlines()[1:-1] if float(line.split(' ')[-1].strip()) < params.res_cutoff]:
					if(orig_seqs[site[0]] in seqs2keep):
						running_aln[orig_seqs[site[0]]][site[1]] = 'X'

				#Removing columns below the --col_cutoff
				cols2remove = [int(line.split()[0]) - 1 for line in open(tax_guidance_outdir + '/MSA.MAFFT.Guidance2_res_pair_col.scr').readlines()[1:-1] if float(line.split(' ')[-1].strip()) < params.col_cutoff]
				for seq in running_aln:
					running_aln[seq] = ''.join([running_aln[seq][i] for i in range(len(running_aln[seq])) if i not in cols2remove])

				with open(tax_guidance_outdir + '/postGuidance_preTrimAl_unaligned.fasta', 'w') as o:
					for seq in running_aln:
						o.write('>' + seq + '\n' + str(running_aln[seq]).replace('-', '') + '\n\n')

				#Aligning one last time after removing the final set of sequences and applying the res and col cutoffs
				print('mafft ' + tax_guidance_outdir + '/postGuidance_preTrimAl_unaligned.fasta > ' + tax_guidance_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '_postGuidance_preTrimAl_aligned.fasta')
				os.system('mafft ' + tax_guidance_outdir + '/postGuidance_preTrimAl_unaligned.fasta > ' + tax_guidance_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '.postGuidance_preTrimAl_aligned.fasta')

				#Gap trimming
				os.system('Scripts/trimal-trimAl/source/trimal -in ' + tax_guidance_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '.postGuidance_preTrimAl_aligned.fasta -out ' + tax_guidance_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '.95gapTrimmed.fasta -gapthreshold 0.05 -fasta')

				#Copying over final aligments (pre and post gap trimming) into output folder.
				os.system('cp ' + tax_guidance_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '.95gapTrimmed.fasta ' + params.output + '/Output/Guidance/' + file.split('.')[0].split('_preguidance')[0] + '.95gapTrimmed.fasta')
				os.system('cp ' + tax_guidance_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '.postGuidance_preTrimAl_aligned.fasta ' + params.output + '/Output/NotGapTrimmed/' + file.split('.')[0].split('_preguidance')[0] + '.postGuidance_preTrimAl_aligned.fasta')
				
				#Removing intermediate files if not --keep_temp
				if not params.keep_temp:
					for gdir_file in os.listdir(tax_guidance_outdir):
						if gdir_file not in ('MSA.MAFFT.Guidance2_res_pair_seq.scr_with_Names', 'MSA.MAFFT.aln.With_Names', 'MSA.MAFFT.Guidance2_res_pair_col.scr', 'log', 'postGuidance_preTrimAl_unaligned.fasta'):
							os.system('rm -r ' + tax_guidance_outdir + '/' + gdir_file)
						else:
							if gdir_file == 'MSA.MAFFT.aln.With_Names':
								os.system('mv ' + tax_guidance_outdir + '/' + gdir_file + ' ' + tax_guidance_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '_' + gdir_file + '.aln')
							else:
								os.system('mv ' + tax_guidance_outdir + '/' + gdir_file + ' ' + tax_guidance_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '_' + gdir_file)

		guidance_removed_file.close()

























