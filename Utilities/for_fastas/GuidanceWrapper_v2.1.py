#Author, date: Auden Cote-L'Heureux, July 17 2023
#Motivation: Run Guidance on multiple files
#Intent: To be used in place of PhyloToL part 2 if there is a good reason not to use the guidance-only mode of PhyloToL
#Dependencies: Python3, BioPython
#Inputs: A folder of unaligned fasta files
#Outputs: Folder for each input file with selected Guidance output files, renamed with the original file name at the front
#Example: python3 GuidanceWrapper_v2.1.py -i FastaFiles


#Dependencies
import os, sys
import argparse
from Bio import SeqIO

#Reading arguments
parser = argparse.ArgumentParser(
                    prog = 'Guidance Wrapper v2.1',
                    description = "Updated July 21, 2023 by Auden Cote-L'Heureux"
                    )
parser.add_argument('--input', '-i', required = True, type = str, help = 'Path to folder of unaligned amino acid fasta files to align. File extensions must be fasta, fa, fas, or faa. Try using the absolute rather than relative path if working on the Grid and having trouble')
parser.add_argument('--output', '-o', default = '.', type = str, help = 'Path to folder where a folder named "GuidanceOutput" will be created')
parser.add_argument('--codon', action = 'store_true', help = 'Run on nucleotide files by translating codons')
parser.add_argument('--iterations', '-n', default = 5, type = int, help = 'Number of Guidance iterations (default = 5)')
parser.add_argument('--guidance_path', '-p', default = 'guidance.v2.02', type = str, help = 'Path to the guidance_v2.02 folder')
parser.add_argument('--seq_cutoff', '-s', default = 0.3, type = float, help = 'Taxa are removed if their score is below this cutoff')
parser.add_argument('--col_cutoff', '-c', default = 0.0, type = float, help = 'Columns are removed if their score is below this cutoff')
parser.add_argument('--res_cutoff', '-r', default = 0.0, type = float, help = 'During guidance, residues are removed if their score is below this cutoff')
parser.add_argument('--force', '-f', action = 'store_true', help = 'Delete existing output folder at given output path')
parser.add_argument('--keep_temp', '-k', action = 'store_true', help = 'Keep all Guidance intermediate files (beware, some can be very large)')
parser.add_argument('--guidance_threads', '-t', default = 20, type = int, help = 'Number of threads to allocate to Guidance')

args = parser.parse_args()

#Creating an output folder
if not os.path.isdir(args.output + '/GuidanceOutput'):
	os.mkdir(args.output + '/GuidanceOutput')
elif not args.force:
	print('\nERROR: An output folder already exists at the given path. Delete it or run in --force mode\n')
	quit()

#For each file in the input folder
for file in os.listdir(args.input):
	#If it has an appropriate extension for an amino acid fasta file
	if file.split('.')[-1] in ('fasta', 'fa', 'fas', 'faa'):
		#Create the output folder
		tax_guidance_outdir = args.output + '/GuidanceOutput/' + file.split('.')[0]
		if(not os.path.isdir(tax_guidance_outdir)):
			os.mkdir(tax_guidance_outdir)

			#Read in the sequence data and replace any * and U characters with 'X'
			recs = [r for r in SeqIO.parse(args.input + '/' + file, 'fasta')]
			with open(args.input + '/' + file, 'w') as o:
				for rec in recs:
					o.write('>' + rec.description + '\n' + str(rec.seq).replace('*', 'X').replace('U', 'X') + '\n\n')

			fail = False
			#For each Guidance iteration (default = 5)
			for i in range(args.iterations):
				#Get the number sequences in the unaligned file
				n_recs = len([r for r in SeqIO.parse(args.input + '/' + file, 'fasta')])

				#Cancel the Guidance run if there are fewer than 4 sequences
				if n_recs < 4:
					print('\nGene famiily ' + file.split('.')[0].split('_preguidance')[0] + ' contains fewer than 4 sequences after ' + str(i) + ' Guidance iterations, therefore no alignment will be produced for this gene family\n')
					os.system('rm -rf ' + tax_guidance_outdir)
					if i == 0:
						fail = True
					break

				#MAFFT parameters depend on sequence count
				if n_recs < 200:
					mafft_alg = 'genafpair'
				else:
					mafft_alg = 'auto'

				if args.codon:
					seqtype = 'codon'
				else:
					seqtype = 'aa'

				#Run Guidance
				os.system(args.guidance_path + '/www/Guidance/guidance.pl --seqFile ' + args.input + '/' + file + ' --msaProgram MAFFT --seqType ' + seqtype + ' --outDir ' + tax_guidance_outdir + ' --seqCutoff ' + str(args.seq_cutoff) + ' --colCutoff ' + str(args.col_cutoff) + " --outOrder as_input --bootstraps 10 --MSA_Param '\\--" + mafft_alg + " --maxiterate 1000 --thread " + str(args.guidance_threads) + " --bl 62 --anysymbol' > " + tax_guidance_outdir + '/log.txt')

				#If it ran successfully
				if os.path.isfile(tax_guidance_outdir + '/MSA.MAFFT.Guidance2_res_pair_seq.scr_with_Names'):
					sep = ' '
					if '\t' in open(tax_guidance_outdir + '/MSA.MAFFT.Guidance2_res_pair_seq.scr_with_Names').readlines()[1]:
						sep = '\t'
					
					#Create a record of sequences below the sequence score cutoff
					seqs_below = len([line for line in open(tax_guidance_outdir + '/MSA.MAFFT.Guidance2_res_pair_seq.scr_with_Names').readlines()[1:-1] if float(line.split(sep)[-1]) < args.seq_cutoff])

					#If the number of remaining sequences is less than 4, then stop iterating
					if n_recs - seqs_below < 4:
						print('\nGene famiily ' + file.split('.')[0].split('_preguidance')[0] + ' contains fewer than 4 sequences after ' + str(i + 1) + ' Guidance iterations, therefore no alignment will be produced for this gene family.\n')
						os.system('rm -rf ' + tax_guidance_outdir)
						break

					#If there are no sequences below the cutoff or the number of Guidance iterations has been reached, then stop iterating
					if seqs_below == 0 or i == args.iterations - 1:
						print('\nGuidance complete after ' + str(i + 1) + ' iterations for gene family ' + file.split('.')[0].split('_preguidance')[0] + '\n')
						break

					#If another iteration is needed, then replace the original input file with the sequence-filtered run to input to Guidance for the next iteration
					os.system('cp ' + tax_guidance_outdir + '/Seqs.Orig.fas.FIXED.Without_low_SP_Seq.With_Names ' + args.input + '/' + file)

					#If another iteration is needed, remove the Guidance output from the last iteration
					os.system('rm -r ' + tax_guidance_outdir + '/*')
				else:
					fail = True
					break

			#If Guidance ran successfully
			if not fail:
				#Create a record of sequences above the sequence score cutoff
				seqs2keep = [rec.description for rec in SeqIO.parse(tax_guidance_outdir + '/Seqs.Orig.fas.FIXED.Without_low_SP_Seq.With_Names', 'fasta')]
				orig_seqs = [rec.description for rec in SeqIO.parse(tax_guidance_outdir + '/MSA.MAFFT.aln.With_Names', 'fasta')]

				#Put this sequence information in a dictionary and write out the surviving sequences to an unaligned .fasta file
				running_aln = { rec.description : str(rec.seq) for rec in SeqIO.parse(tax_guidance_outdir + '/MSA.MAFFT.aln.With_Names', 'fasta') if rec.description in seqs2keep }
				with open(tax_guidance_outdir + '/postGuidance_preMAFFT_unaligned.fasta', 'w') as o:
					for seq in running_aln:
						o.write('>' + seq + '\n' + str(running_aln[seq]).replace('-', '') + '\n\n')

				#Align the file with MAFFT
				os.system('mafft ' + tax_guidance_outdir + '/postGuidance_preMAFFT_unaligned.fasta > ' + tax_guidance_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '.postGuidance_MAFFT_realigned.fasta')

				#Read in the MAFFT alignment
				running_aln = { rec.description : str(rec.seq) for rec in SeqIO.parse(tax_guidance_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '.postGuidance_MAFFT_realigned.fasta', 'fasta') }

				sep = ' '
				if '\t' in open(tax_guidance_outdir + '/MSA.MAFFT.Guidance2_res_pair_seq.scr').readlines()[1]:
					sep = '\t'
				
				#Apply residue cutoff per site per sequence
				for site in [(int(line.split(sep)[1]), int(line.split(sep)[0]) - 1) for line in open(tax_guidance_outdir + '/MSA.MAFFT.Guidance2_res_pair_seq.scr').readlines()[1:-1] if float(line.split(sep)[-1].strip()) < args.res_cutoff]:
					if(orig_seqs[site[0]] in seqs2keep):
						running_aln[orig_seqs[site[0]]][site[1]] = 'X'

				sep = ' '
				if '\t' in open(tax_guidance_outdir + '/MSA.MAFFT.Guidance2_res_pair_col.scr').readlines()[1]:
					sep = '\t'

				#Apply column cutoff per column
				cols2remove = [int(line.split(sep)[0]) - 1 for line in open(tax_guidance_outdir + '/MSA.MAFFT.Guidance2_res_pair_col.scr').readlines()[1:-1] if float(line.split(sep)[-1].strip()) < args.col_cutoff]
				for seq in running_aln:
					running_aln[seq] = ''.join([running_aln[seq][i] for i in range(len(running_aln[seq])) if i not in cols2remove])

				#Write out the final filtered alignment
				with open(tax_guidance_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '.postGuidance_MAFFT_realigned.fasta', 'w') as o:
					for seq in running_aln:
						o.write('>' + seq + '\n' + str(running_aln[seq]) + '\n\n')				
				
				#If the user hasn't specified to keep all intermediate files, then delete the extraneous ones and rename the rest per the input file name
				if not args.keep_temp:
					for gdir_file in os.listdir(tax_guidance_outdir):
						if gdir_file not in ('MSA.MAFFT.Guidance2_res_pair_seq.scr_with_Names', 'MSA.MAFFT.aln.With_Names', 'MSA.MAFFT.Guidance2_res_pair_col.scr', 'log', file.split('.')[0].split('_preguidance')[0] + '.postGuidance_MAFFT_realigned.fasta'):
							os.system('rm -r ' + tax_guidance_outdir + '/' + gdir_file)
						else:
							if gdir_file == 'MSA.MAFFT.aln.With_Names':
								os.system('mv ' + tax_guidance_outdir + '/' + gdir_file + ' ' + tax_guidance_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '_' + gdir_file + '.aln')
							else:
								os.system('mv ' + tax_guidance_outdir + '/' + gdir_file + ' ' + tax_guidance_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '_' + gdir_file)









