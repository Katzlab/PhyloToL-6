#Author, date: Auden Cote-L'Heureux, last updated Aug 1st 2023
#Motivation: Select robust sequences from trees
#Intent: Select clades of interest from large trees using taxonomic specifications
#Dependencies: Python3, ete3, Biopython
#Inputs: A folder of trees and corresponding unaligned .fasta files
#Outputs: A folder of grabbed clades and filtered unaligned fasta files
#Example: python CladeGrabbing_v2.1.py --input /Path/to/trees --target Sr_rh --min_presence 20

#Dependencies
import os, re, sys
from Bio import SeqIO
import ete3
import argparse


def get_args():

	parser = argparse.ArgumentParser(
		prog = 'Clade grabber, Version 2.1',
		description = "Updated Aug 1st, 2023 by Auden Cote-L'Heureux"
	)

	parser.add_argument('-i', '--input', type = str, required = True, help = 'Path to a folder containing input trees (which must have the file extension .tre, .tree, .treefile, or .nex)')
	parser.add_argument('-t', '--target', type = str, required = True, help = 'A comma-separated list of any number of digits/characters to describe focal taxa (e.g. Sr_ci_S OR Am_t), or a file with the extension .txt containing a list of complete or partial taxon codes. All sequences containing the complete/partial code will be identified as belonging to target taxa.')
	parser.add_argument('-m', '--min_presence', type = int, required = True, help = 'Minimum number of target taxa present in clade for it to be selected')
	parser.add_argument('-a', '--at_least', type = str, default = '', help = 'A comma-separated list of any number of digits/characters (e.g. Sr_ci_S OR Am_t), or a file with the extension .txt containing a list of complete or partial taxon codes, to describe taxa that MUST be present in a clade for it to be selected (e.g. you may want at least one whole genome).')
	parser.add_argument('-na', '--n_at_least', type = int, default = 0, help = 'The number of species belonging to taxa in the --at_least list that must be present in the clade. Default is 1.')
	parser.add_argument('-o', '--outgroup', type = str, default = '', help = 'A comma-separated list of any number of digits/characters (e.g. Sr_ci_S OR Am_t), or a file with the extension .txt containing a list of complete or partial taxon codes, to describe taxa that will be included as outgroups in the output unaligned fasta files (which will contain only sequences from a single selected clade, and all outgroup sequences in the tree captured by this argument).')
	parser.add_argument('-c', '--contaminants', type = float, default = 2, help = 'The number of non-ingroup contaminants allowed in a clade, or if less than 1 the proportion of sequences in a clade that can be non-ingroup (i.e. presumed contaminants). Default is to allow 2 contaminants.')
	
	return parser.parse_args()

	
def get_newick(fname):
	
	newick = ''
	for line in open(fname):
		line = line.split(' ')[-1]
		if(line.startswith('(') or line.startswith('tree1=')):
			newick = line.split('tree1=')[-1].replace("'", '').replace('\\', '')

	return newick


#This function reroots the tree on the largest Ba/Za clade. If there is no prokaryote clade,
#it roots on the largest Op clade, then Pl, then Am, then Ex, then Sr.
def reroot(tree):

	#This nested function returns the largest clade of a given taxonomic group
	def get_best_clade(taxon):

		best_size = 0; best_clade = []; seen_leaves = []
		#Traverse all nodes
		for node in tree.traverse('levelorder'):
			#If the node is big enough and not subsumed by a node we've already accepted
			if len(node) >= 3 and len(list(set(seen_leaves) & set([leaf.name for leaf in node]))) == 0:
				leaves = [leaf.name for leaf in node]
				
				#Create a record of leaves that belong to the taxonomic group
				target_leaves = set()
				for leaf in leaves[::-1]:
					if leaf[:2] in taxon:
						target_leaves.add(leaf[:10])
						leaves.remove(leaf)

				#If this clade is better than any clade we've seen before, grab it
				if len(target_leaves) > best_size and len(leaves) <= 2:
					best_clade = node
					best_size = len(target_leaves)
					seen_leaves.extend([leaf.name for leaf in node])

		return best_clade

	#Get the biggest clade for each taxonomic group (stops once it finds one)
	for taxon in [('Ba', 'Za'), ('Op'), ('Pl'), ('Am'), ('Ex'), ('Sr')]:
		clade = get_best_clade(taxon)
		if len([leaf for leaf in clade if leaf.name[:2] in taxon]) > 3:
			tree.set_outgroup( clade)

			break

	return tree
	
	
def get_subtrees(args, file):

	newick = get_newick(args.input + '/' + file)	

	tree = ete3.Tree(newick)

	majs = list(dict.fromkeys([leaf.name[:2] for leaf in tree]))

	#Only try to reroot trees with more than 2 major clades. This was added to fix the ETE3 "Cannot set myself as outgroup" error
	if len(majs) > 2:
		tree = reroot(tree)								

	#Getting a clean list of all target taxa
	if '.' in args.target:
		try:
			target_codes = [l.strip() for l in args.target.open().readlines() if l.strip() != '']
		except AttributeError:
			print('\n\nError: invalid "target" argument. This must be a comma-separated list of any number of digits/characters to describe focal taxa (e.g. Sr_ci_S OR Am_t), or a file with the extension .txt containing a list of complete or partial taxon codes. All sequences containing the complete/partial code will be identified as belonging to target taxa.\n\n')
	else:
		target_codes = [code.strip() for code in args.target.split(',') if code.strip() != '']

	#Getting a clean list of all "at least" taxa
	if '.' in args.at_least:
		try:
			at_least_codes = [l.strip() for l in args.at_least.open().readlines() if l.strip() != '']
		except AttributeError:
			print('\n\nError: invalid "at_least" argument. This must be a comma-separated list of any number of digits/characters (e.g. Sr_ci_S OR Am_t), or a file with the extension .txt containing a list of complete or partial taxon codes, to describe taxa that MUST be present in a clade for it to be selected (e.g. you may want at least one whole genome).\n\n')
	else:
		at_least_codes = [code.strip() for code in args.at_least.split(',') if code.strip() != '']

	target_codes = list(dict.fromkeys(target_codes + at_least_codes))

	#Creating a record of selected subtrees, and all of the leaves in those subtrees
	selected_nodes = []; seen_leaves = []

	#Iterating through all nodes in tree, starting at "root" then working towards leaves
	for node in tree.traverse('levelorder'):
		#If a node is large enough and is not contained in an already selected clade
		if len(node) >= args.min_presence and len(list(set(seen_leaves) & set([leaf.name for leaf in node]))) == 0:
			leaves = [leaf.name for leaf in node]

			#Accounting for cases where e.g. one child is a contaminant, and the other child is a good clade with 1 fewer than the max number of contaminants
			children_keep = 0
			for child in node.children:
				for code in target_codes:
					for leaf in child:
						if leaf.name.startswith(code):
							children_keep += 1
							break

			if children_keep == len(node.children):
				#Creating a record of all leaves belonging to the target/"at least" group of taxa, and any other leaves are contaminants
				target_leaves = set(); at_least_leaves = set()
				for code in target_codes:
					for leaf in leaves[::-1]:
						if leaf.startswith(code):
							target_leaves.add(leaf[:10])

							if code in at_least_codes:
								at_least_leaves.add(leaf[:10])

							leaves.remove(leaf)

				#Grab a clade as a subtree if 1) it has enough target taxa; 2) it has enough "at least" taxa; 3) it does not have too many contaminants
				if len(target_leaves) >= args.min_presence and len(at_least_leaves) >= args.n_at_least and ((args.contaminants < 1 and len(leaves) < args.contaminants * len(target_leaves)) or len(leaves) < args.contaminants):
					selected_nodes.append(node)
					seen_leaves.extend([leaf.name for leaf in node])

	#Write the subtrees to output .tre files
	for i, node in enumerate(selected_nodes[::-1]):
		with open('Subtrees/' + '.'.join(file.split('.')[:-1]) + '_' + str(i) + '.tre', 'w') as o:
			o.write(node.write())
	

def make_new_unaligned(args):

	if not os.path.isdir('Subtrees_unaligned'):
		os.mkdir('Subtrees_unaligned')

	#Getting a clean list of outgroup taxa
	if '.' in args.outgroup:
		try:
			outgroup_codes = [l.strip() for l in args.outgroup.open().readlines() if l.strip() != '']
		except AttributeError:
			print('\n\nError: invalid "target" argument. This must be a comma-separated list of any number of digits/characters (e.g. Sr_ci_S OR Am_t), or a file with the extension .txt containing a list of complete or partial taxon codes, to describe taxa that will be included as outgroups in the output unaligned fasta files (which will contain only sequences from a single selected clade, and all outgroup sequences in the tree captured by this argument).\n\n')
	else:
		outgroup_codes = [code.strip() for code in args.outgroup.split(',') if code.strip() != '']

	for tree_file in os.listdir('Subtrees'):
		if tree_file.endswith('.tre'):
			og = tree_file[:10]

			tree = ete3.Tree('Subtrees/' + tree_file)
			
			#Get the fasta (aligned or unaligned, but if aligned not columns removed) for each subtree
			for f, file in enumerate(os.listdir(args.input)):
				if file.startswith(og) and file.split('.')[-1] in ('fa', 'faa', 'fna', 'fasta', 'fas'):
					with open('Subtrees_unaligned/' + tree_file.replace('.tre', '.fasta'), 'w') as o:
						#For all records in the fasta file
						for rec in SeqIO.parse(args.input + '/' + file, 'fasta'):
							keep = False
							#Keep it if it's in the outgroup taxa list
							for code in outgroup_codes:
								if rec.id.startswith(code):
									keep = True
									break

							#Or keep it if it's in the subtree
							if rec.description in [leaf.name for leaf in tree]:
								keep = True

							#Write the sequence to the output file if kept
							if keep:
								o.write('>' + rec.description + '\n' + str(rec.seq).replace('-', '') + '\n\n')
						break

def main():
	
	args = get_args()
	
	if(not os.path.isdir('Subtrees')):
		os.mkdir('Subtrees')
	
	f = 0
	for file in os.listdir(args.input):
		if file.split('.')[-1] in ('tre', 'tree', 'treefile', 'nex'):
			print(str(f + 1) + '. ' + file)
			f += 1
				
			get_subtrees(args, file)
						
	make_new_unaligned(args)
													
	
main()

	