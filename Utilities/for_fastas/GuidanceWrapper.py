#Author: Auden Cote-L'Heureux
#Last update: 03/01/2023; added replacement of non-AA characters '*' and 'U' by 'X' after discussion with LAK and ML on 02/28/2023


import os, sys
from Bio import SeqIO

batchdir = sys.argv[1]

seq_cutoff = 0.3
col_cutoff = 0.4

for file in os.listdir(batchdir):
	if(not os.path.isdir('/beegfs/fast/katzlab/acl/Hook-DEV/Guidance_G1NS/GuidanceOutput/' + file.split('.')[0])):
		os.mkdir('/beegfs/fast/katzlab/acl/Hook-DEV/Guidance_G1NS/GuidanceOutput/' + file.split('.')[0])

		recs = [r for r in SeqIO.parse(batchdir + '/' + file, 'fasta')]
		with open(batchdir + '/' + file, 'w') as o:
			for rec in recs:
				o.write('>' + rec.description + '\n' + str(rec.seq).replace('*', 'X').replace('U', 'X') + '\n\n')

		n_recs = len(recs)

		for i in range(5):

			if n_recs < 4:
				print('Fewer than 4 sequences in gene family ' + file.split('.')[0])
				#os.system('rm -rf ' + tax_guidance_outdir)
				break

			if n_recs < 200:
				mafft_alg = 'genafpair'
			else:
				mafft_alg = 'auto'

			os.system('/beegfs/fast/katzlab/acl/Hook-DEV/Guidance_G1NS/guidance.v2.02/www/Guidance/guidance.pl --seqFile ' + batchdir + '/' + file + ' --msaProgram MAFFT --seqType aa --outDir /beegfs/fast/katzlab/acl/Hook-DEV/Guidance_G1NS/GuidanceOutput/' + file.split('.')[0] + ' --seqCutoff ' + str(seq_cutoff) + ' --colCutoff ' + str(col_cutoff) + " --outOrder as_input --bootstraps 10 --MSA_Param '\\--" + mafft_alg + " --maxiterate 1000 --thread 20 --bl 62 --anysymbol'")

			seqs_below = len([line for line in open('GuidanceOutput/' + file.split('.')[0] + '/MSA.MAFFT.Guidance2_res_pair_seq.scr_with_Names').readlines()[1:-1] if float(line.split()[-1]) < seq_cutoff])

			if n_recs - seqs_below < 4:
				print('Fewer than 4 sequences left after ' + str(i + 1) + ' iterations in gene family ' + file.split('.')[0])
				#os.system('rm -rf GuidanceOutput/' + file.split('.')[0])
				break

			if seqs_below == 0 or i == 4:
				print('Guidance complete after ' + str(i + 1) + ' iterations for gene family ' + file.split('.')[0])
				break

			os.system('cp GuidanceOutput/' + file.split('.')[0] + '/Seqs.Orig.fas.FIXED.Without_low_SP_Seq.With_Names ' + batchdir + '/' + file)
			os.system('rm -r GuidanceOutput/' + file.split('.')[0])