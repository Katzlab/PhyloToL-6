#Intent: Summarize the taxonomic distribution of sister sequences for each taxon in a set of trees
#Requirements: A folder of trees
#Example script call: python ContaminationBySisters_v2.1.py --input /Path/to/trees
#For more info: python ContaminationBySisters_v2.0.py --help

#Dependencies
import os, sys, re
import ete3
import argparse
from statistics import mean


def get_args():

	parser = argparse.ArgumentParser(
		prog = 'Contamination by Sisters, Version 2.1',
		description = "Updated July 5th, 2023 by Auden Cote-L'Heureux."
	)

	parser.add_argument('-i', '--input', type = str, required = True, help = 'Path to the folder containing the input trees')
	parser.add_argument('-l', '--level', choices = ['major', 'minor', 'species'], default = 'species', help = 'The level at which to report sister relationships. Options are major, minor, and species.')
	parser.add_argument('-b', '--branch_length_filter', type = float, help = 'Branch length filter to use: this number will be multiplied by the mean branch length in the tree, and only branches smaller than that will be considered (e.g. enter 1.5 to only report branches smaller than 150 percent the average branch length). If you skip this argument, branches will not be filtered for length')
	parser.add_argument('-1', '--single_sister_only', action = 'store_true', help = 'Only include sequences in the final summary file with one sister taxon')
	parser.add_argument('-q', '--query_clades', nargs = '+', help = 'A list of 2, 4, 5, 7, 8, or 10 digit codes specifying which taxa for which you would like sisters reported (e.g. -s Am,Ba,Pl_gr will report sisters to taxa that are Amoebozoa, Bacteria, or green algae), separated by a comma. Alternatively, input a file containing a list of 10 digit codes of taxa for sisters to represent if there are a lot')
	parser.add_argument('-s', '--sister_clades', nargs = '+', help = 'A list of 2, 4, 5, 7, 8, or 10 digit codes specifying which sister taxa to report (e.g. -s Am,Ba,Pl_gr will report sisters that are Amoebozoa, Bacteria, or green algae), separated by a comma. Alternatively, input a file containing a list of 10 digit codes of taxa for sisters to represent if there are a lot')

	args = parser.parse_args()

	def reset_arg(arg):
		if arg != None:
			temp = []
			for val in arg:
				if os.path.isfile(val):
					temp.extend([line.strip() for line in open(val) if line.strip() != ''])
				else:
					temp.append(val.strip())
			
			return temp
		else:
			return None

	args.query_clades = reset_arg(args.query_clades)
	args.sister_clades = reset_arg(args.sister_clades)

	return args
	
	
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
	

def write_all_data(args):

	report = open('PerSequenceData' + args.run_info_string + '.csv', 'w')
	report.write('Tree,Taxon,Sequence,Result,Sisters,BranchLength,BranchLength.Mean\n')
				
	#Iterating over all the input trees
	for file in os.listdir(args.input):	
		if file.split('.')[-1] in ('tre', 'tree', 'treefile', 'nex'):
			
			#Read the tree using ete3 and reroot it using the above function
			newick = get_newick(args.input + '/' + file)
			tree = ete3.Tree(newick)
			tree = reroot(tree)

			#Get the average branch length at terminal nodes
			mean_bl = mean([leaf.dist for leaf in tree])

			#For each sequence
			for leaf in tree:

				#Test whether it is in once of the input query clades
				consider = False
				if args.query_clades == None:
					consider = True
				else:
					for clade in args.query_clades:
						if leaf.startswith(clade):
							consider = True
							break

				if consider:

					#This loop will keep moving towards the root of the tree until it finds a node that
					#has leaves from a cell other than the one for which we are looking for sisters
					parent_node = leaf; sister_taxa = {leaf.name[:10]}
					while len(sister_taxa) == 1:
						parent_node = parent_node.up
						for l2 in parent_node:
							sister_taxa.add(l2.name[:10])

					#Create a record of the sister sequences
					sisters = [sister for sister in parent_node if sister.name[:10] != leaf.name[:10]]
					#And the sisters' minor clades
					sister_minors = list(dict.fromkeys([sister.name[:5] for sister in sisters]))

					#Classify the taxonomic distribution of sisters
					if len(sister_minors) == 1:
						if sister_minors[0] == leaf.name[:5]:
							result = 'same_minor'
						else:
							result = sister_minors[0]
					else:
						result = 'non-monophyletic'

					#Write to the output file
					report.write (file + ',' + leaf.name[:10] + ',' + leaf.name + ',' +  result + ',' + ' '.join(list(dict.fromkeys([sister.name for sister in sisters]))) + ',' + str(leaf.dist) + ',' + str(mean_bl) + '\n')

	report.close()

	
#This summarizes the huge table initially generated, which isn't very human friendly
def summarize(args):
	
	#Reading the output file from above and creating a summary file
	report = open('PerSequenceData' + args.run_info_string + '.csv', 'r').readlines()[1:]
	output = open('PerTaxonSummary' + args.run_info_string + '.csv', 'w')
	
	#A list of all the trees and taxa
	trees = sorted(list(dict.fromkeys([line.split(",")[0] for line in report])))
	queries = sorted(list(dict.fromkeys([line.split(",")[1] for line in report])))
	
	sisters_by_taxon = { }; stats_by_taxon = { }
	#For each line in the huge table
	for line in report:
		line = line.split(',')
		#Apply the branch-length filter if activated
		if (args.branch_length_filter != None and float(line[5]) < float(line[6]) * args.branch_length_filter) or (args.branch_length_filter == None):
			#Get a record of sister taxa for this sequence
			sister_taxa = list(dict.fromkeys([sis[:10] for sis in line[4].split(' ') if sis.strip() != '']))
			
			#If applicable, only consider the sequence if it has a single sister
			if not args.single_sister_only or len(sister_taxa) == 1:
				if(line[1] not in sisters_by_taxon):
					sisters_by_taxon.update({ line[1] : { } })
					stats_by_taxon.update({ line[1] : { 'same' : 0, 'total' : 0, 'nm' : 0} })

				#Counting up the number of "results" per taxon
				if line[3] == 'same_minor':
					stats_by_taxon[line[1]]['same'] += 1
				elif line[3] == 'non-monophyletic':
					stats_by_taxon[line[1]]['nm'] += 1

				stats_by_taxon[line[1]]['total'] += 1

				#For each sister species, add to the count
				for sis in sister_taxa:

					if args.level == 'major':
						if sis[:2] not in sisters_by_taxon[line[1]]:
							sisters_by_taxon[line[1]].update({ sis[:2] : 0 })

						sisters_by_taxon[line[1]][sis[:2]] += 1
					elif args.level == 'minor':
						if sis[:2] not in sisters_by_taxon[line[1]]:
							sisters_by_taxon[line[1]].update({ sis[:2] : 0 })

						if sis[:5] not in sisters_by_taxon[line[1]]:
							sisters_by_taxon[line[1]].update({ sis[:5] : 0 })

						sisters_by_taxon[line[1]][sis[:5]] += 1
						sisters_by_taxon[line[1]][sis[:2]] += 1
					elif args.level == 'species':
						if sis[:2] not in sisters_by_taxon[line[1]]:
							sisters_by_taxon[line[1]].update({ sis[:2] : 0 })

						if sis not in sisters_by_taxon[line[1]]:
							sisters_by_taxon[line[1]].update({ sis : 0 })

						sisters_by_taxon[line[1]][sis] += 1
						sisters_by_taxon[line[1]][sis[:2]] += 1

	#Organize the list of sister taxonomic groups that were counted
	sisters_to_write = sorted(list(dict.fromkeys([sister for query in sisters_by_taxon for sister in sisters_by_taxon[query]])))

	#Write the output headers
	output.write('Taxon,SameMinor,NM,Total,PropSM,' + ','.join([sis for sis in sisters_to_write if len(sis) == 2]) + ',' + ','.join([sis for sis in sisters_to_write if len(sis) > 2]) + '\n')

	#Write out the counts for sister relationships
	for taxon in sisters_by_taxon:
		output.write(taxon + ',' + str(stats_by_taxon[taxon]['same']) + ',' + str(stats_by_taxon[taxon]['nm']) + ',' + str(stats_by_taxon[taxon]['total']) + ',' + str(round(stats_by_taxon[taxon]['same']/stats_by_taxon[taxon]['total'], 2)) + ',')

		for sister in [sis for sis in sisters_to_write if len(sis) == 2]:
			if sister in sisters_by_taxon[taxon]:
				output.write(str(sisters_by_taxon[taxon][sister]) + ',')
			else:
				output.write('0,')
		for sister in [sis for sis in sisters_to_write if len(sis) > 2]:
			if sister in sisters_by_taxon[taxon]:
				output.write(str(sisters_by_taxon[taxon][sister]) + ',')
			else:
				output.write('0,')
		
		output.write('\n')
	

#A wrapper to call all above functions
def main():
	
	#Getting the parameters
	args = get_args()

	run_info_string = '_'
	if args.single_sister_only:
		run_info_string = run_info_string + 'SingleSisterOnly_'
	else:
		run_info_string = run_info_string + 'SisterClades_'

	if args.branch_length_filter == None:
		run_info_string = run_info_string + 'NoBL_'
	else:
		run_info_string = run_info_string + 'BLFilter' + str(args.branch_length_filter) + 'x_'

	args.run_info_string = run_info_string + args.level
	
	#Writing the big not-human-friendly spreadsheet
	write_all_data(args)
	
	#Writing the summary spreadsheet
	summarize(args)
	
	
#Calling the main wrapper function
main()








			