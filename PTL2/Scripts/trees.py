# Last updated Sept 2023
# Authors: Auden Cote-L'Heureux and Mario Ceron-Romero

# This is a relatively simple script that only runs trees, using either IQ-Tree
# or RAxML. The run() function is called in two places: both in phylotol.py, and
# in contamination.py, where it is used to re-build trees. When starting at this
# step, users must input one aligned amino acid fasta file per OG. Otherwise, if 
# starting at the pre-Guidance or Guidance steps, this step will be run if --end = trees.

#Dependencies
import os, sys, re
from Bio import SeqIO
from color import color

#Called in phylotol.py and contamination.py
def run(params):

	#Checking whether aligned files were input, or it should just start with the Guidance outputs from the previous step.
	if params.start == 'aligned':
		guidance_path = params.data
	else:
		guidance_path = params.output + '/Output/Guidance'

	#Looking for input files
	if not os.path.isdir(guidance_path):
		print('\nERROR: The path ' + guidance_path + ' could not be found when trying to locate Guidance (aligned) files. Make sure that the --start and --data parameters are correct and/or that the Guidance step ran successfully.\n')
		exit()
		
	if len([f for f in os.listdir(guidance_path) if f.endswith('.fa') or f.endswith('.faa') or f.endswith('.fasta') or f.endswith('.fas') or f.endswith('.aln')]) == 0:
		print('\nERROR: No Guidance (unaligned) files could be found at the path ' + guidance_path + '. Make sure that the --start and --data parameters are correct, that the Guidance step ran successfully, and that the aligned files are formatted correctly (they must have the file extension .faa, .fa, .aln, .fas, or .fasta).\n')
		exit()

	#For each input alignment
	for file in [f for f in os.listdir(guidance_path)  if f.endswith('.fa') or f.endswith('.faa') or f.endswith('.fasta') or f.endswith('.fas') or f.endswith('.aln')]:

		#Run IQ-Tree
		if params.tree_method == 'iqtree' or params.tree_method == 'iqtree_fast':
			#Make intermediate folders
			if not os.path.isdir(params.output + '/Output/Intermediate/IQTree'):
				os.mkdir(params.output + '/Output/Intermediate/IQTree')
				
			tax_iqtree_outdir = params.output + '/Output/Intermediate/IQTree/' + file.split('.')[0].split('_preguidance')[0]
			os.mkdir(tax_iqtree_outdir)

			#Run IQ-Tree
			if params.tree_method == 'iqtree':
				os.system('iqtree2 -s ' + guidance_path + '/' + file + ' -m LG+G -T 10 --prefix ' + tax_iqtree_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '.IQTree')
			elif params.tree_method == 'iqtree_fast':
				os.system('iqtree2 -s ' + guidance_path + '/' + file + ' -m LG+G -T 10 --fast --prefix ' + tax_iqtree_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '.IQTree')
				
			#Copy over the final output
			if os.path.isfile(tax_iqtree_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '.IQTree.treefile'):
				os.system('cp ' + tax_iqtree_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '.IQTree.treefile ' + params.output + '/Output/Trees/' + file.split('.')[0].split('_preguidance')[0] + '.IQTree.tree')
				#color(params.output + '/Output/Trees/' + file.split('.')[0].split('_preguidance')[0] + '.IQTree.tree')
			else:
				print('\nWARNING: No tree file created by IQ-Tree for OG ' + file[:10] + '\n')

		#If not IQ-Tree, then run RAxML
		elif params.tree_method == 'raxml':
			#Make intermediate folders
			if not os.path.isdir(params.output + '/Output/Intermediate/RAxML'):
				os.mkdir(params.output + '/Output/Intermediate/RAxML')
				
			tax_raxml_outdir = params.output + '/Output/Intermediate/RAxML/' + file.split('.')[0].split('_preguidance')[0]
			os.mkdir(tax_raxml_outdir)

			#Reformat the alignment as phylip
			os.system('./Scripts/trimal-trimAl/source/trimal -in ' + guidance_path + '/' + file + ' -phylip -out ' + tax_raxml_outdir + '/aligned.phy')

			#Run RAxML
			os.system('raxmlHPC -s ' + tax_raxml_outdir + '/aligned.phy -m PROTGAMMALG -f d -p 12345 -# 10 -n ' + file.split('.')[0].split('_preguidance')[0] + '_RAxML -T ' + str(params.guidance_threads))
			
			#Copy over final output
			if os.path.isfile(tax_raxml_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '.IQTree.treefile'):
				os.system('cp ' + tax_raxml_outdir + '/RAxML_bestTree.' + file.split('.')[0].split('_preguidance')[0] + '_RAxML ' + params.output + '/Output/Trees/' + file.split('.')[0].split('_preguidance')[0] + '_RAxML.tree')
				#color(params.output + '/Output/Trees/' + file.split('.')[0].split('_preguidance')[0] + '_RAxML.tree')
			else:
				print('\nWARNING: No tree file created by RAxML for OG ' + file[:10] + '\n')



































