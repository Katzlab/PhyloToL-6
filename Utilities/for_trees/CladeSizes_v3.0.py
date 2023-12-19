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
	
	
def get_clades(file, args):

	newick = get_newick(file)	

	tree = ete3.Tree(newick)

	#Root tree if possible
	try:
		tree = reroot(tree)
	except:
		print('\nUnable to re-root the tree ' + file + ' (maybe it had only 1 major clade, or an inconvenient polytomy). Skipping this step and continuing to try to grab robust clades from the tree.\n')					

	if args.taxon_group != None:
		majs = []; mins = []
		majs.append([line.strip() for line in open(args.taxon_group)])
	else:
		majs = list(dict.fromkeys([leaf.name[:2] for leaf in tree]))
		mins = list(dict.fromkeys([leaf.name[:5] for leaf in tree]))


	#For each major and minor clade, find all monophyletic clades of that taxon
	clades_per_tax = { }; majs_per_clade = { }; mins_per_clade = { }; seen_leaves_majs = []
	for clade_code in majs + mins:

		if type(clade_code) is str:
			clades_per_tax.update({ clade_code : [] })
		elif type(clade_code) is list:
			clades_per_tax.update({ args.taxon_group : [] })
			majs_per_clade.update({ args.taxon_group : [] })
			mins_per_clade.update({ args.taxon_group : [] })

		seen_leaves = []
		for node in tree.traverse('levelorder'):
			if len(list(set(seen_leaves) & set([leaf.name for leaf in node]))) == 0 and len([leaf for leaf in node]) >= args.size:
				if (len(clade_code) == 5 and args.single_mins_only and len(list(set(seen_leaves_majs) & set([leaf.name for leaf in node]))) == 0) or (len(clade_code) == 2 or not args.single_mins_only):
					
					leaves = [leaf.name for leaf in node]

					taxa = list(dict.fromkeys([leaf[:10] for leaf in leaves]))

					target_leaves = set()
					for leaf in leaves[::-1]:
						if type(clade_code) is str:
							if leaf.startswith(clade_code):
								target_leaves.add(leaf[:10])
						elif type(clade_code) is list:
							for code in clade_code:
								if leaf.startswith(code):
									target_leaves.add(leaf[:10])
									break

					minors = len(list(dict.fromkeys([l[:5] for l in target_leaves])))

					if len(target_leaves) == len(taxa) and len(target_leaves) >= args.size:
						if type(clade_code) is str and ((len(clade_code) == 2 and args.majs_req_mult_mins and minors > 1) or (len(clade_code) == 2 and not args.majs_req_mult_mins) or len(clade_code) == 5):
							clades_per_tax[clade_code].append(len(list(dict.fromkeys([l[:10] for l in target_leaves]))))
							seen_leaves.extend([leaf.name for leaf in node])

							if len(clade_code) == 2:
								seen_leaves_majs.extend(seen_leaves)

						elif type(clade_code) is list:
							clades_per_tax[args.taxon_group].append(len(list(dict.fromkeys([l[:10] for l in target_leaves]))))
							majs_per_clade[args.taxon_group].append(len(list(dict.fromkeys([l[:2] for l in target_leaves]))))
							mins_per_clade[args.taxon_group].append(minors)

							seen_leaves.extend([leaf.name for leaf in node])

					

	clades_per_tax = { clade : sorted(clades_per_tax[clade], key = lambda x : -x ) for clade in clades_per_tax }


	return clades_per_tax, majs_per_clade, mins_per_clade



#Write the output spreadsheet
def write_output(clades_per_tax_per_file, args, majs_per_clade = None, mins_per_clade = None):

	clades = sorted(list(dict.fromkeys([clade for file in clades_per_tax_per_file for clade in clades_per_tax_per_file[file] if len(clade) == 2]))) + sorted(list(dict.fromkeys([clade for file in clades_per_tax_per_file for clade in clades_per_tax_per_file[file] if len(clade) == 5])))
	
	if args.taxon_group != None:
		clades = clades + [args.taxon_group]

	sumfilt = { }
	if args.summary_filter != None:
		if os.path.isfile(args.summary_filter):
			for line in open(args.summary_filter):
				line = line.strip().split('\t')
				if len(line) == 2:
					try:
						sumfilt.update({ line[0] : int(line[1]) })
					except ValueError:
						print('\nError: the line ' + '\t'.join(line) + ' in the input summary filter file is not properly formatted.\n')
						exit()
		else:
			print('\nError: it looks like you tried to input a summary filter file, but this file could not be found.\n')
			exit()

	sumclades = sorted(list([key for key in sumfilt.keys() if len(key) == 2])) + sorted(list([key for key in sumfilt.keys() if len(key) == 5]))

	with open('CladeSizesPerTaxon.csv', 'w') as o:
		o.write('OG,')

		if len(sumclades) > 0:
			o.write(','.join(sumclades))
			sum_iter_group = sumclades
		else:
			o.write(','.join(clades))
			sum_iter_group = clades

		if args.taxon_group == None:
			o.write(',' + ','.join(clades) + '\n')
		else:
			o.write(',' + ','.join(clades) + ',Majs,Mins\n')

		for tree_file in clades_per_tax_per_file:
			o.write(tree_file)
			
			for clade in sum_iter_group:
				if clade in clades_per_tax_per_file[tree_file]:
					minsize = 0
					if len(sumclades) > 0:
						minsize = sumfilt[clade]

					nclades = len([size for size in clades_per_tax_per_file[tree_file][clade] if size >= minsize])

					o.write(',' + str(nclades))
				else:
					o.write(',0')

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

			if args.taxon_group != None:
				o.write(',' + ' '.join(['(' + str(size) + ')' for size in majs_per_clade[tree_file][args.taxon_group]]))
				o.write(',' + ' '.join(['(' + str(size) + ')' for size in mins_per_clade[tree_file][args.taxon_group]]))

			o.write('\n')


if __name__ == '__main__':

	parser = argparse.ArgumentParser(
                    prog = 'Clade sizes utility version 2.0',
                    description = "Updated Dec 1, 2023 by Auden Cote-L'Heureux"
                    )

	parser.add_argument('-i', '--input', required = True, help = 'Path to folder of tree files')
	parser.add_argument('-t', '--taxon_group', help = 'Optional. Path to a file with a list of taxa (e.g. phototrophs) that will be considered a single group in clade grabbing.')
	parser.add_argument('-s', '--size', default = 1, type = int, help = 'Only count clades larger than a certain size')
	parser.add_argument('--majs_req_mult_mins', action = 'store_true', help = 'Flag: if set, major clades will only be counted that contan multiple minor clades.')
	parser.add_argument('--single_mins_only', action = 'store_true', help = 'Flag: if set, clades will only be counted under the minor clade category if they contain a single minor clade.')
	parser.add_argument('--summary_filter', type = str, help = 'Path to a file with two tab-separated columns to filter the "summary" section of the output (which reports the number of clades above a certain given size for a given set of taxonomic groups) the first column should list all taxonomic groups to summarize, and the second column the minimum size of the clades in order to be included in the summary for each group.')


	args = parser.parse_args()
	
	clades_per_tax_per_file = { }; majs_per_clade_per_file = { }; mins_per_clade_per_file = { }
	for tree_file in tqdm(os.listdir(args.input)):
		if tree_file.split('.')[-1] in ('tre', 'tree', 'treefile', 'nex'):
			clades_per_tax, majs_per_clade, mins_per_clade = get_clades(args.input + '/' + tree_file, args)
			clades_per_tax_per_file.update({ tree_file.split('.')[0] : clades_per_tax })
			majs_per_clade_per_file.update({ tree_file.split('.')[0] : majs_per_clade })
			mins_per_clade_per_file.update({ tree_file.split('.')[0] : mins_per_clade })

	write_output(clades_per_tax_per_file, args, majs_per_clade = majs_per_clade_per_file, mins_per_clade = mins_per_clade_per_file)





































