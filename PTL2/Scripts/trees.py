import os, sys, re
from Bio import SeqIO
from logger import Logger
from color import color

def run(params):

	if params.start == 'aligned':
		guidance_path = params.data
	else:
		guidance_path = params.output + '/Output/Guidance'

	if not os.path.isdir(guidance_path):
		Logger.Error('The path ' + guidance_path + ' could not be found when trying to locate Guidance (aligned) files. Make sure that the --start and --data parameters are correct and/or that the Guidance step ran successfully.')

	if len([f for f in os.listdir(guidance_path) if f.endswith('.fa') or f.endswith('.faa') or f.endswith('.fasta') or f.endswith('.fas') or f.endswith('.aln')]) == 0:
		Logger.Error('No Guidance (unaligned) files could be found at the path ' + guidance_path + '. Make sure that the --start and --data parameters are correct, that the Guidance step ran successfully, and that the aligned files are formatted correctly (they must have the file extension .faa, .fa, .aln, .fas, or .fasta).')

	for file in [f for f in os.listdir(guidance_path)  if f.endswith('.fa') or f.endswith('.faa') or f.endswith('.fasta') or f.endswith('.fas') or f.endswith('.aln')]:
		
		if params.tree_method == 'iqtree':
			if not os.path.isdir(params.output + '/Output/Intermediate/IQTree'):
				os.mkdir(params.output + '/Output/Intermediate/IQTree')

			tax_iqtree_outdir = params.output + '/Output/Intermediate/IQTree/' + file.split('.')[0].split('_preguidance')[0]
			os.mkdir(tax_iqtree_outdir)

			os.system('iqtree2 -s ' + guidance_path + '/' + file + ' -m LG+G -T AUTO --prefix ' + tax_iqtree_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '_IQTree')
			
			if os.path.isfile(tax_iqtree_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '_IQTree.treefile'):
				os.system('cp ' + tax_iqtree_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '_IQTree.treefile ' + params.output + '/Output/Trees/' + file.split('.')[0].split('_preguidance')[0] + '_IQTree.tree')
				#color(params.output + '/Output/Trees/' + file.split('.')[0].split('_preguidance')[0] + '_IQTree.tree')
			else:
				Logger.Warning('No tree file created by IQ-Tree for OG ' + file[:10])

		elif params.tree_method == 'raxml':
			if not os.path.isdir(params.output + '/Output/Intermediate/RAxML'):
				os.mkdir(params.output + '/Output/Intermediate/RAxML')
				
			tax_raxml_outdir = params.output + '/Output/Intermediate/RAxML/' + file.split('.')[0].split('_preguidance')[0]
			os.mkdir(tax_raxml_outdir)

			os.system('./Scripts/trimal-trimAl/source/trimal -in ' + guidance_path + '/' + file + ' -phylip -out ' + tax_raxml_outdir + '/aligned.phy')

			print('raxmlHPC -s ' + tax_raxml_outdir + '/aligned.phy -m PROTGAMMALG -f d -p 12345 -# 10 -n ' + file.split('.')[0].split('_preguidance')[0] + '_RAxML -T ' + str(params.guidance_threads))
			os.system('raxmlHPC -s ' + tax_raxml_outdir + '/aligned.phy -m PROTGAMMALG -f d -p 12345 -# 10 -n ' + file.split('.')[0].split('_preguidance')[0] + '_RAxML -T ' + str(params.guidance_threads))

			if os.path.isfile(tax_raxml_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '_IQTree.treefile'):
				os.system('cp ' + tax_raxml_outdir + '/RAxML_bestTree.' + file.split('.')[0].split('_preguidance')[0] + '_RAxML ' + params.output + '/Output/Trees/' + file.split('.')[0].split('_preguidance')[0] + '_RAxML.tree')
				#color(params.output + '/Output/Trees/' + file.split('.')[0].split('_preguidance')[0] + '_RAxML.tree')
			else:
				Logger.Warning('No tree file created by RAxML for OG ' + file[:10])



































