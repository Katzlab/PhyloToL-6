# Last updated March 2024 by Godwin Ani
# Authors: Auden Cote-L'Heureux and Mario Ceron-Romero



#Dependencies
import os, sys, re
import argparse
from Bio import SeqIO

#Reading arguments
parser = argparse.ArgumentParser(
                    prog = 'Guidance Wrapper v2.1',
                    description = "Updated July 21, 2023 by Auden Cote-L'Heureux"
                    )
parser.add_argument('--input', '-i', required = True, type = str, help = 'Path to folder of unaligned amino acid fasta files to align. File extensions must be fasta, fa, fas, or faa. Try using the absolute rather than relative path if working on the Grid and having trouble')
parser.add_argument('--output', '-o', default = '.', type = str, help = 'Path to folder where files will be created')
parser.add_argument('--codon', action = 'store_true', help = 'Run on nucleotide files by translating codons')
parser.add_argument('--iterations', '-n', default = 5, type = int, help = 'Number of Guidance iterations (default = 5)')
parser.add_argument('--guidance_path', '-p', default = 'guidance.v2.02', type = str, help = 'Path to the guidance_v2.02 folder')
parser.add_argument('--seq_cutoff', '-s', default = 0.3, type = float, help = 'Taxa are removed if their score is below this cutoff')
parser.add_argument('--col_cutoff', '-c', default = 0.0, type = float, help = 'Columns are removed if their score is below this cutoff')
parser.add_argument('--res_cutoff', '-r', default = 0.0, type = float, help = 'During guidance, residues are removed if their score is below this cutoff')
parser.add_argument('--force', '-f', action = 'store_true', help = 'Delete existing output folder at given output path')
parser.add_argument('--keep_temp', '-k', action = 'store_true', help = 'Keep all Guidance intermediate files (beware, some can be very large)')
parser.add_argument('--guidance_threads', '-t', default = 20, type = int, help = 'Number of threads to allocate to Guidance')
parser.add_argument('--keep_iter', '-z', action = 'store_true', help = 'Keep all Guidance intermediate files (beware, some can be very large)')

args = parser.parse_args()

#Creating intermedate folders that will later be deleted unless running with --keep_temp

os.makedirs(args.output + '/Output/Intermediate/Guidance')
os.mkdir(args.output + '/Output/Intermediate/Guidance/Input')
os.mkdir(args.output + '/Output/Intermediate/Guidance/Output')
os.mkdir(args.output + '/Output/Guidance')
os.mkdir(args.output + '/Output/NotGapTrimmed')
guidance_input = args.output + 'Output/Intermediate/Guidance/Input/'
os.system('cp -r ' + args.input + '/* ' + guidance_input)

guidance_removed_file = open(args.output + '/Output/GuidanceRemovedSeqs.txt', 'w')
guidance_removed_file.write('Sequence\tScore\n')

#For each unaligned AA fasta file
for file in [f for f in os.listdir(guidance_input) if f.endswith('.fa') or f.endswith('.faa') or f.endswith('.fasta')]:
	tax_guidance_outdir = args.output + '/Output/Intermediate/Guidance/Output/' + file.split('.')[0].split('_preguidance')[0]
	os.mkdir(tax_guidance_outdir)
	fail = False
    #For each iteration
	for i in range(args.iterations):
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
		#Determining if its nucleotide or AA	
		if args.codon:
			seqtype = 'codon'
		else:
			seqtype = 'aa'

		#Run Guidance - new
		os.system(args.guidance_path + '/www/Guidance/guidance.pl --seqFile ' + guidance_input + '/' + file + ' --msaProgram MAFFT --seqType ' + seqtype + ' --outDir ' + tax_guidance_outdir + ' --seqCutoff ' + str(args.seq_cutoff) + ' --colCutoff ' + str(args.col_cutoff) + " --outOrder as_input --bootstraps 10 --MSA_Param '\\--" + mafft_alg + " --maxiterate 1000 --thread " + str(args.guidance_threads) + " --bl 62 --anysymbol' > " + args.output + '/Output/Intermediate/Guidance/Output/' + file[:10] + '/log.txt')
                
              
		#Checking for a sequence score file; if not available, Guidance failed.
		if os.path.isfile(tax_guidance_outdir + '/MSA.MAFFT.Guidance2_res_pair_seq.scr_with_Names'):
		#All sequences below score cutoff
			seqs_below = len([line for line in open(tax_guidance_outdir + '/MSA.MAFFT.Guidance2_res_pair_seq.scr_with_Names').readlines()[1:-1] if float(line.split()[-1]) < args.seq_cutoff])
			#If fewer than four were above the cutoff, this OG is done iterating.
			if n_recs - seqs_below < 4:
				print('\nWARNING: Gene famiily ' + file.split('.')[0].split('_preguidance')[0] + ' contains fewer than 4 sequences after ' + str(i + 1) + ' Guidance iterations, therefore no alignment will be produced for this gene family.\n')
				os.system('rm -rf ' + tax_guidance_outdir)
				break
			#If all sequences were above the cutoff, this OG is done iterating.
			if seqs_below == 0 or i == args.iterations - 1:
				print('\nGuidance complete after ' + str(i + 1) + ' iterations for gene family ' + file.split('.')[0].split('_preguidance')[0] + '\n')
				break
			#Recording list of sequences removed by Guidance.
			for line in [line for line in open(tax_guidance_outdir + '/MSA.MAFFT.Guidance2_res_pair_seq.scr_with_Names').readlines()[1:-1] if float(line.split()[-1]) < args.seq_cutoff]:
				guidance_removed_file.write(line)
			#Copying over the old file with the new results
			os.system('cp ' + tax_guidance_outdir + '/Seqs.Orig.fas.FIXED.Without_low_SP_Seq.With_Names ' + guidance_input + '/' + file)
			#Handling intermediate files for each iteration.	
			if args.keep_iter:
				if i +1 < args.iterations:
					os.makedirs(args.output + '/Output/iterations/', exist_ok = True)
					os.makedirs(args.output + '/Output/iterations/' + str(i+1)+'/', exist_ok = True)
					os.makedirs(args.output + '/Output/iterations/' + str(i+1) + '/' + file.split('.')[0].split('_preguidance')[0], exist_ok = True)
					iteration_folder = args.output + '/Output/iterations/'+ str(i +1) + '/' + file.split('.')[0].split('_preguidance')[0]
					os.system('cp -r ' + tax_guidance_outdir + '/* ' + iteration_folder)
					os.system('rm -r ' + tax_guidance_outdir + '/*')
					if not args.keep_temp:
						for gdir_file in os.listdir(iteration_folder):
							if gdir_file not in ('MSA.MAFFT.Guidance2_res_pair_seq.scr_with_Names', 'MSA.MAFFT.aln.With_Names', 'MSA.MAFFT.Guidance2_res_pair_col.scr', 'log', 'postGuidance_preTrimAl_unaligned.fasta'):
								os.system('rm -r ' + iteration_folder + '/' + gdir_file)
							else:
								if gdir_file == 'MSA.MAFFT.aln.With_Names':
									os.system('mv ' + iteration_folder + '/' + gdir_file + ' ' + iteration_folder + '/' + file.split('.')[0].split('_preguidance')[0] + '_' + gdir_file + '.aln')
								else:
									os.system('mv ' + iteration_folder + '/' + gdir_file + ' ' + iteration_folder + '/' + file.split('.')[0].split('_preguidance')[0] + '_' + gdir_file)
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
		for site in [(int(line.split()[1]), int(line.split()[0]) - 1) for line in open(tax_guidance_outdir + '/MSA.MAFFT.Guidance2_res_pair_seq.scr').readlines()[1:-1] if float(line.split(' ')[-1].strip()) < args.res_cutoff]:
			if(orig_seqs[site[0]] in seqs2keep):
				running_aln[orig_seqs[site[0]]][site[1]] = 'X'

		#Removing columns below the --col_cutoff
		cols2remove = [int(line.split()[0]) - 1 for line in open(tax_guidance_outdir + '/MSA.MAFFT.Guidance2_res_pair_col.scr').readlines()[1:-1] if float(line.split(' ')[-1].strip()) < args.col_cutoff]
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
		os.system('cp ' + tax_guidance_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '.95gapTrimmed.fasta ' + args.output + '/Output/Guidance/' + file.split('.')[0].split('_preguidance')[0] + '.95gapTrimmed.fasta')
		os.system('cp ' + tax_guidance_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '.postGuidance_preTrimAl_aligned.fasta ' + args.output + '/Output/NotGapTrimmed/' + file.split('.')[0].split('_preguidance')[0] + '.postGuidance_preTrimAl_aligned.fasta')
				
		#Removing intermediate files if not --keep_temp
		if not args.keep_temp:
			for gdir_file in os.listdir(tax_guidance_outdir):
				if gdir_file not in ('MSA.MAFFT.Guidance2_res_pair_seq.scr_with_Names', 'MSA.MAFFT.aln.With_Names', 'MSA.MAFFT.Guidance2_res_pair_col.scr', 'log', 'postGuidance_preTrimAl_unaligned.fasta'):
					os.system('rm -r ' + tax_guidance_outdir + '/' + gdir_file)
				else:
					if gdir_file == 'MSA.MAFFT.aln.With_Names':
						os.system('mv ' + tax_guidance_outdir + '/' + gdir_file + ' ' + tax_guidance_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '_' + gdir_file + '.aln')
					else:
						os.system('mv ' + tax_guidance_outdir + '/' + gdir_file + ' ' + tax_guidance_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '_' + gdir_file)

guidance_removed_file.close()
