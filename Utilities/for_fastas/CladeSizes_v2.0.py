#Author, date: Auden Cote-L'Heureux, last updated Dec 1st 2023
#Motivation: Understand the topology of trees
#Intent: Describe clade sizes for different taxonomic groups
#Dependencies: Python3, ete3
#Inputs: A folder of trees
#Outputs: a spreadsheet describing clade sizes
#Example: python CladeSizes_v2.0.py -i /path/to/trees


import os, sys, re
import ete3
import argparse
from tqdm import tqdm


#Extract Newick string
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
	
	
def get_clades(file):

	newick = get_newick(file)	

	tree = ete3.Tree(newick)

	#Root tree if possible
	try:
		tree = reroot(tree)
	except:
		print('\nUnable to re-root the tree ' + file + ' (maybe it had only 1 major clade, or an inconvenient polytomy). Skipping this step and continuing to try to grab robust clades from the tree.\n')					

	majs = list(dict.fromkeys([leaf.name[:2] for leaf in tree]))
	mins = list(dict.fromkeys([leaf.name[:5] for leaf in tree]))

	#For each major and minor clade, find all monophyletic clades of that taxon
	clades_per_tax = { }
	for clade_code in majs + mins:
		clades_per_tax.update({ clade_code : [] })

		seen_leaves = []
		for node in tree.traverse('levelorder'):
			if len(list(set(seen_leaves) & set([leaf.name for leaf in node]))) == 0:
				leaves = [leaf.name for leaf in node]

				taxa = list(dict.fromkeys([leaf[:10] for leaf in leaves]))

				target_leaves = set()
				for leaf in leaves[::-1]:
					if leaf.startswith(clade_code):
						target_leaves.add(leaf[:10])

				if len(target_leaves) == len(taxa):
					clades_per_tax[clade_code].append(len(list(dict.fromkeys([l[:10] for l in target_leaves]))))
					seen_leaves.extend([leaf.name for leaf in node])

	clades_per_tax = { clade : sorted(clades_per_tax[clade], key = lambda x : -x ) for clade in clades_per_tax }

	return clades_per_tax


#Write the output spreadsheet
def write_output(clades_per_tax_per_file):

	clades = sorted(list(dict.fromkeys([clade for file in clades_per_tax_per_file for clade in clades_per_tax_per_file[file] if len(clade) == 2]))) + sorted(list(dict.fromkeys([clade for file in clades_per_tax_per_file for clade in clades_per_tax_per_file[file] if len(clade) == 5])))
	with open('CladeSizesPerTaxon.csv', 'w') as o:
		o.write('OG,' + ','.join(clades) + '\n')
		for tree_file in clades_per_tax_per_file:
			o.write(tree_file)
			for clade in clades:
				if clade in clades_per_tax_per_file[tree_file]:
					if len(clades_per_tax_per_file[tree_file][clade]) > 1:
						o.write(',' + ' '.join(['(' + str(size) + ')' for size in clades_per_tax_per_file[tree_file][clade]]))
					elif len(clades_per_tax_per_file[tree_file][clade]) == 1:
						o.write(',' + str(clades_per_tax_per_file[tree_file][clade][0]))
					else:
						o.write(',')
				else:
					o.write(',')

			o.write('\n')


if __name__ == '__main__':

	parser = argparse.ArgumentParser(
                    prog = 'Clade sizes utility version 2.0',
                    description = "Updated Dec 1, 2023 by Auden Cote-L'Heureux"
                    )

	parser.add_argument('-i', '--input', required = True, help = 'Path to folder of tree files')

	args = parser.parse_args()
	
	clades_per_tax_per_file = { }
	for tree_file in tqdm(os.listdir(args.input)):
		if tree_file.split('.')[-1] in ('tre', 'tree', 'treefile', 'nex'):
			clades_per_tax = get_clades(args.input + '/' + tree_file)
			clades_per_tax_per_file.update({ tree_file.split('.')[0] : clades_per_tax })

	write_output(clades_per_tax_per_file)





































