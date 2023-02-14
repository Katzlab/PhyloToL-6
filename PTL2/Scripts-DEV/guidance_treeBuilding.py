import os, sys, re
from Bio import SeqIO


def run(params):

	if params.start == 'raw':
		preguidance_path = params.output + '/Output/Pre-Guidance'
	elif params.start == 'unaligned':
		preguidance_path = params.data

	if not os.path.isdir(preguidance_path):
		Logger.Error('The path ' + preguidance_path + ' could not be found when trying to locate pre-Guidance (unaligned) files. Make sure that the --start and --data parameters are correct and that the pre-Guidance step ran successfully.')

	if len([f for f in os.listdir(preguidance_path) if f.endswith('.fa') or f.endswith('.faa') or f.endswith('.fasta')]) == 0:
		Logger.Error('No pre-Guidance (unaligned) files could be found at the path ' + preguidance_path + '. Make sure that the --start and --data parameters are correct, that the pre-Guidance step ran successfully, and that the unaligned files are formatted correctly (they must have the file extension .faa, .fa, or .fasta).')

	os.mkdir(params.output + '/Output/Temp/Guidance')
	os.mkdir(params.output + '/Output/Temp/Guidance/Input')
	os.mkdir(params.output + '/Output/Temp/Guidance/Output')

	guidance_input = params.output + '/Output/Temp/Guidance/Input/'
	os.system('cp -r ' + preguidance_path + '/* ' + guidance_input)

	for file in [f for f in os.listdir(guidance_input) if f.endswith('.fa') or f.endswith('.faa') or f.endswith('.fasta')]:
		tax_guidance_outdir = params.output + '/Output/Temp/Guidance/Output/' + file.split('.')[0].split('_preguidance')[0]
		os.mkdir(tax_guidance_outdir)

		for i in range(params.guidance_iters):
			n_recs = len([r for r in SeqIO.parse(guidance_input + '/' + file, 'fasta')])

			if n_recs < 4:
				Logger.Warning('Gene famiily ' + file.split('.')[0].split('_preguidance')[0] + ' contains fewer than 4 sequences after ' + str(i) + ' Guidance iterations, therefore no alignment will be produced for this gene family.')
				os.system('rm -rf ' + tax_guidance_outdir)
				break

			if n_recs < 200:
				mafft_alg = 'genafpair'
			else:
				mafft_alg = 'auto'

			os.system('guidance.v2.02/www/Guidance/guidance.pl --seqFile ' + guidance_input + '/' + file + ' --msaProgram MAFFT --seqType aa --outDir ' + tax_guidance_outdir + ' --seqCutoff ' + str(params.seq_cutoff) + ' --colCutoff ' + str(params.col_cutoff) + " --outOrder as_input --bootstraps 10 --MSA_Param '\\--'" + mafft_alg + " --maxiterate 1000 --thread " + str(params.guidance_threads) + " --bl 62 --anysymbol' > " + params.output + '/Output/Temp/Guidance/Output/' + file[:10] + '/log.txt')

			exit()
			seqs_below = len([line for line in open(tax_guidance_outdir + '/MSA.MAFFT.Guidance2_res_pair_seq.scr_with_Names').readlines()[1:] if float(line.split()[-1]) < params.seqcutoff])

			if n_recs - seqs_below < 4:
				Logger.Warning('Gene famiily ' + file.split('.')[0].split('_preguidance')[0] + ' contains fewer than 4 sequences after ' + str(i + 1) + ' Guidance iterations, therefore no alignment will be produced for this gene family.')
				os.system('rm -rf ' + tax_guidance_outdir)
				break

			if len(seqs_below) == 0 or i == params.guidance_iters - 1:
				Logger.Message('Guidance complete after ' + str(i + 1) + ' iterations for gene family ' + file.split('.')[0].split('_preguidance')[0])
				break

			os.system('cp ' + tax_guidance_outdir + '/Seqs.Orig.fas.FIXED.Without_low_SP_Seq.With_Names ' + guidance_input + '/' + file)

		seqs2keep = [rec.description for rec in SeqIO.parse('Seqs.Orig.fas.FIXED.Without_low_SP_Seq.With_Names')]
		orig_seqs = [rec.description for rec in SeqIO.parse('MSA.MAFFT.aln.With_Names', 'fasta')]
		running_aln = { rec.description : str(rec.seq) for rec in SeqIO.parse('MSA.MAFFT.aln.With_Names', 'fasta') if rec.description in seqs2keep }

		for site in [(int(line.split()[1]), int(line.split()[0]) - 1) for line in open(file.split('.')[0].split('_preguidance')[0] + '/MSA.MAFFT.Guidance2_res_pair_seq.scr').readlines()[1:] if float(line.split([2])) < params.res_cutoff]:
			if(orig_seqs[site[0]] in seqs2keep):
				running_aln[orig_seqs[site[0]]][site[1]] = 'X'

		cols2remove = [int(line.split()[0] - 1) for line in open(file.split('.')[0].split('_preguidance')[0] + '/MSA.MAFFT.Guidance2_res_pair_col.scr').readlines()[1:] if float(line.split([1])) < params.col_cutoff]
		for seq in running_aln:
			running_aln[seq] = ''.join([running_aln[seq][i] for i in range(len(running_aln[seq])) if i not in cols2remove])

		with open(tax_guidance_outdir + '/postGuidance_preTrimAl.fasta', 'w') as o:
			for seq in running_aln:
				o.write('>' + seq + '\n' + str(running_aln[seq]) + '\n\n')

		os.system('trimal-trimAl/source/trimal -in ' + tax_guidance_outdir + '/postGuidance_preTrimAl.fasta -out ' + tax_guidance_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '.95gapTrimmed.fas -gapthreshold 0.05 -fasta')
		os.system('trimal-trimAl/source/trimal -in ' + tax_guidance_outdir + '/postGuidance_preTrimAl.fasta -out ' + tax_guidance_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '.95gapTrimmed.phy -gapthreshold 0.05 -phylip')

		os.system('cp ' + tax_guidance_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '.95gapTrimmed.fas ' + params.output + '/Output/Guidance/' + file.split('.')[0].split('_preguidance')[0] + '.95gapTrimmed.fas')

		if params.end == 'trees':
			if params.tree_method == 'iqtree':
				if not os.path.isdir(params.output + '/Output/Temp/IQTree'):
					os.mkdir(params.output + '/Output/Temp/IQTree')
				tax_iqtree_outdir = params.output + '/Output/Temp/IQTree/' + file.split('.')[0].split('_preguidance')[0]
				os.mkdir(tax_iqtree_outdir)

				#RUN IQ-TREE

			elif params.tree_method == 'raxml':
				if not os.path.isdir(params.output + '/Output/Temp/RAxML'):
					os.mkdir(params.output + '/Output/Temp/RAxML')
				tax_raxml_outdir = params.output + '/Output/Temp/RAxML/' + file.split('.')[0].split('_preguidance')[0]
				os.mkdir(tax_raxml_outdir)

				#RUN RAXML



































