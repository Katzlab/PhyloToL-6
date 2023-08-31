import os, sys, re
from Bio import SeqIO
from color import color

def run(params):

	if params.start == 'aligned':
		guidance_path = params.data
	else:
		guidance_path = params.output + '/Output/Guidance'

	if not os.path.isdir(guidance_path):
		print('\nThe path ' + guidance_path + ' could not be found when trying to locate Guidance (aligned) files. Make sure that the --start and --data parameters are correct and/or that the Guidance step ran successfully.\n')

	if len([f for f in os.listdir(guidance_path) if f.endswith('.fa') or f.endswith('.faa') or f.endswith('.fasta') or f.endswith('.fas') or f.endswith('.aln')]) == 0:
		print('\nNo Guidance (unaligned) files could be found at the path ' + guidance_path + '. Make sure that the --start and --data parameters are correct, that the Guidance step ran successfully, and that the aligned files are formatted correctly (they must have the file extension .faa, .fa, .aln, .fas, or .fasta).\n')

	for file in [f for f in os.listdir(guidance_path)  if f.endswith('.fa') or f.endswith('.faa') or f.endswith('.fasta') or f.endswith('.fas') or f.endswith('.aln')]:
		
		if params.tree_method == 'iqtree':
			if not os.path.isdir(params.output + '/Output/Intermediate/IQTree'):
				os.mkdir(params.output + '/Output/Intermediate/IQTree')

			tax_iqtree_outdir = params.output + '/Output/Intermediate/IQTree/' + '.'.join(file[:-1])
			os.mkdir(tax_iqtree_outdir)

			os.system('iqtree2 -s ' + guidance_path + '/' + file + ' -m LG+G --prefix ' + tax_iqtree_outdir + '/' + file + '_IQTree')
			
			if os.path.isfile(tax_iqtree_outdir + '/' + '.'.join(file[:-1]) + '_IQTree.treefile'):
				os.system('cp ' + tax_iqtree_outdir + '/' + '.'.join(file[:-1]) + '_IQTree.treefile ' + params.output + '/Output/Trees/' + '.'.join(file[:-1]) + '_IQTree.tree')
				#color(params.output + '/Output/Trees/' + file.split('.')[0].split('_preguidance')[0] + '_IQTree.tree')
			else:
				print('\nNo tree file created by IQ-Tree for OG ' + file[:10] + '\n')

		elif params.tree_method == 'raxml':
			if not os.path.isdir(params.output + '/Output/Intermediate/RAxML'):
				os.mkdir(params.output + '/Output/Intermediate/RAxML')
				
			tax_raxml_outdir = params.output + '/Output/Intermediate/RAxML/' + '.'.join(file[:-1])
			os.mkdir(tax_raxml_outdir)

			os.system('./Scripts/trimal-trimAl/source/trimal -in ' + guidance_path + '/' + file + ' -phylip -out ' + tax_raxml_outdir + '/aligned.phy')

			print('raxmlHPC-PTHREADS-AVX2 -s ' + tax_raxml_outdir + '/aligned.phy -m PROTGAMMALG -f d -p 12345 -# 1 -w ' + tax_raxml_outdir + ' -n ' + '.'.join(file[:-1]) + '_RAxML -T ' + str(params.guidance_threads))
			os.system('raxmlHPC-PTHREADS-AVX2 -s ' + tax_raxml_outdir + '/aligned.phy -m PROTGAMMALG -f d -p 12345 -# 1 -w ' + tax_raxml_outdir + ' -n ' + '.'.join(file[:-1]) + '_RAxML -T ' + str(params.guidance_threads))

			if os.path.isfile(tax_raxml_outdir + '/RAxML_bestTree.' + '.'.join(file[:-1]) + '_RAxML'):
				os.system('cp ' + tax_raxml_outdir + '/RAxML_bestTree.' + '.'.join(file[:-1]) + '_RAxML ' + params.output + '/Output/Trees/' + '.'.join(file[:-1]) + '_RAxML.tree')
				#color(params.output + '/Output/Trees/' + file.split('.')[0].split('_preguidance')[0] + '_RAxML.tree')
			else:
				print('\nNo tree file created by RAxML for OG ' + file[:10] + '\n')




