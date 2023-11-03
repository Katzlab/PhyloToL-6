import os, sys, re
from Bio import SeqIO
from logger import Logger

def run(params):

	if params.start == 'raw' or params.start == 'unaligned':
		if params.start == 'raw':
			preguidance_path = params.output + '/Output/Pre-Guidance'
		else:
			preguidance_path = params.data

		if not os.path.isdir(preguidance_path):
			Logger.Error('The path ' + preguidance_path + ' could not be found when trying to locate pre-Guidance (unaligned) files. Make sure that the --start and --data parameters are correct and/or that the pre-Guidance step ran successfully.')

		if len([f for f in os.listdir(preguidance_path) if f.endswith('.fa') or f.endswith('.faa') or f.endswith('.fasta')]) == 0:
			Logger.Error('No pre-Guidance (unaligned) files could be found at the path ' + preguidance_path + '. Make sure that the --start and --data parameters are correct, that the pre-Guidance step ran successfully, and that the unaligned files are formatted correctly (they must have the file extension .faa, .fa, or .fasta).')

		os.mkdir(params.output + '/Output/Temp/Guidance')
		os.mkdir(params.output + '/Output/Temp/Guidance/Input')
		os.mkdir(params.output + '/Output/Temp/Guidance/Output')

		guidance_input = params.output + '/Output/Temp/Guidance/Input/'
		os.system('cp -r ' + preguidance_path + '/* ' + guidance_input)

		guidance_removed_file = open(params.output + '/Output/GuidanceRemovedSeqs.txt', 'w')
		guidance_removed_file.write('Sequence\tScore\n')

		for file in [f for f in os.listdir(guidance_input) if f.endswith('.fa') or f.endswith('.faa') or f.endswith('.fasta')]:
			tax_guidance_outdir = params.output + '/Output/Temp/Guidance/Output/' + file.split('.')[0].split('_preguidance')[0]
			os.mkdir(tax_guidance_outdir)

			fail = False
			for i in range(params.guidance_iters):
				n_recs = len([r for r in SeqIO.parse(guidance_input + '/' + file, 'fasta')])

				if n_recs < 4:
					Logger.Warning('Gene famiily ' + file.split('.')[0].split('_preguidance')[0] + ' contains fewer than 4 sequences after ' + str(i) + ' Guidance iterations, therefore no alignment will be produced for this gene family.')
					os.system('rm -rf ' + tax_guidance_outdir)
					if i == 0:
						fail = True
					break

				if n_recs < 200:
					mafft_alg = 'genafpair'
				else:
					mafft_alg = 'auto'

				os.system('Scripts/guidance.v2.02/www/Guidance/guidance.pl --seqFile ' + guidance_input + '/' + file + ' --msaProgram MAFFT --seqType aa --outDir ' + tax_guidance_outdir + ' --seqCutoff ' + str(params.seq_cutoff) + ' --colCutoff ' + str(params.col_cutoff) + " --outOrder as_input --bootstraps 10 --MSA_Param '\\--" + mafft_alg + " --maxiterate 1000 --thread " + str(params.guidance_threads) + " --bl 62 --anysymbol' > " + params.output + '/Output/Temp/Guidance/Output/' + file[:10] + '/log.txt')

				if os.path.isfile(tax_guidance_outdir + '/MSA.MAFFT.Guidance2_res_pair_seq.scr_with_Names'):
					seqs_below = len([line for line in open(tax_guidance_outdir + '/MSA.MAFFT.Guidance2_res_pair_seq.scr_with_Names').readlines()[1:-1] if float(line.split()[-1]) < params.seq_cutoff])

					if n_recs - seqs_below < 4:
						Logger.Warning('Gene famiily ' + file.split('.')[0].split('_preguidance')[0] + ' contains fewer than 4 sequences after ' + str(i + 1) + ' Guidance iterations, therefore no alignment will be produced for this gene family.')
						os.system('rm -rf ' + tax_guidance_outdir)
						break

					if seqs_below == 0 or i == params.guidance_iters - 1:
						Logger.Message('Guidance complete after ' + str(i + 1) + ' iterations for gene family ' + file.split('.')[0].split('_preguidance')[0])
						break

					for line in seqs_below:
						guidance_removed_file.write(line)

					os.system('cp ' + tax_guidance_outdir + '/Seqs.Orig.fas.FIXED.Without_low_SP_Seq.With_Names ' + guidance_input + '/' + file)

					os.system('rm -r ' + tax_guidance_outdir + '/*')
				else:
					fail = True
					break

			if not fail:
				seqs2keep = [rec.description for rec in SeqIO.parse(tax_guidance_outdir + '/Seqs.Orig.fas.FIXED.Without_low_SP_Seq.With_Names', 'fasta')]
				orig_seqs = [rec.description for rec in SeqIO.parse(tax_guidance_outdir + '/MSA.MAFFT.aln.With_Names', 'fasta')]
				running_aln = { rec.description : str(rec.seq) for rec in SeqIO.parse(tax_guidance_outdir + '/MSA.MAFFT.aln.With_Names', 'fasta') if rec.description in seqs2keep }

				for site in [(int(line.split()[1]), int(line.split()[0]) - 1) for line in open(tax_guidance_outdir + '/MSA.MAFFT.Guidance2_res_pair_seq.scr').readlines()[1:-1] if float(line.split(' ')[-1].strip()) < params.res_cutoff]:
					if(orig_seqs[site[0]] in seqs2keep):
						running_aln[orig_seqs[site[0]]][site[1]] = 'X'

				cols2remove = [int(line.split()[0]) - 1 for line in open(tax_guidance_outdir + '/MSA.MAFFT.Guidance2_res_pair_col.scr').readlines()[1:-1] if float(line.split(' ')[-1].strip()) < params.col_cutoff]
				for seq in running_aln:
					running_aln[seq] = ''.join([running_aln[seq][i] for i in range(len(running_aln[seq])) if i not in cols2remove])

				with open(tax_guidance_outdir + '/postGuidance_preTrimAl_unaligned.fasta', 'w') as o:
					for seq in running_aln:
						o.write('>' + seq + '\n' + str(running_aln[seq]).replace('-', '') + '\n\n')

				print('mafft ' + tax_guidance_outdir + '/postGuidance_preTrimAl_unaligned.fasta > ' + tax_guidance_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '_postGuidance_preTrimAl_aligned.fasta')
				os.system('mafft ' + tax_guidance_outdir + '/postGuidance_preTrimAl_unaligned.fasta > ' + tax_guidance_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '.postGuidance_preTrimAl_aligned.fasta')

				os.system('Scripts/trimal-trimAl/source/trimal -in ' + tax_guidance_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '.postGuidance_preTrimAl_aligned.fasta -out ' + tax_guidance_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '.95gapTrimmed.fasta -gapthreshold 0.05 -fasta')

				os.system('cp ' + tax_guidance_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '.95gapTrimmed.fasta ' + params.output + '/Output/Guidance/' + file.split('.')[0].split('_preguidance')[0] + '.95gapTrimmed.fasta')
				os.system('cp ' + tax_guidance_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '.postGuidance_preTrimAl_aligned.fasta ' + params.output + '/Output/NotGapTrimmed/' + file.split('.')[0].split('_preguidance')[0] + '.postGuidance_preTrimAl_aligned.fasta')
				
				
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

























