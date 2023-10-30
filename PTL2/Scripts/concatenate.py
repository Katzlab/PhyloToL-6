import os, sys
from Bio import SeqIO
import ete3
import argparse
from tqdm import tqdm


#Small utility function to extract newick strings from nexus file
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


#Function to select sequences to use per tree
def remove_paralogs(params):

	seqs_per_og = { }
	for file in tqdm(os.listdir(params.output + '/Output/Guidance')):
		if file.split('.')[-1] in ('fasta', 'fas', 'faa'):

			prefix = '.'.join(file.split('.')[:-1])
			tre_f = [t for t in os.listdir(params.output + '/Output/Trees') if t.startswith(prefix)]
			if len(tre_f) == 0:
				tre_f = [t for t in os.listdir(params.output + '/Output/Trees') if t.startswith(prefix.split('.')[0])]

				if len(tre_f) == 0:

					tre_f = [t for t in os.listdir(params.output + '/Output/Trees') if t.startswith(file[:10])]
					
					if len(tre_f) == 0:
						print('\nNo tree file found for alignment ' + file + '. Skipping this gene family.\n')
						continue
					elif len(tre_f) > 1:
						print('\nMore than one tree file found matching the alignment file ' + file + '. Please make your file names unique: there should be one alignment file for every tree file, with a matching unique prefix (everything before the first "."). Skipping this gene family.\n')
						continue
				elif len(tre_f) > 1:
					print('\nMore than one tree file found matching the alignment file ' + file + '. Please make your file names unique: there should be one sequence file for every tree file, with a matching unique prefix (everything before the first "."). Skipping this gene family.\n')
					continue
			elif len(tre_f) > 1:
				print('\nMore than one tree file found matching the alignment file ' + file + '. Please make your file names unique: there should be one sequence file for every tree file, with a matching unique prefix (everything before the first "."). Skipping this gene family.\n')
				continue

			seqs_per_og.update({ file : [] })
			og_recs = { rec.id : rec for rec in SeqIO.parse(params.output + '/Output/Guidance/' + file, 'fasta')}

			newick = get_newick(params.output + '/Output/Trees/' + tre_f[0])
			tree = ete3.Tree(newick)

			try:
				tree = reroot(tree)
			except:
				print('\nUnable to re-root the tree ' + file + ' (maybe it had only 1 major clade, or an inconvenient polytomy). Skipping this step and continuing to try to grab robust clades from the tree.\n')					

			#Getting a clean list of all target taxa

			if type(params.concat_target_taxa) is list:
				target_codes = [code.strip() for code in params.concat_target_taxa if code.strip() != '']
			elif params.concat_target_taxa != None:
				if os.path.isfile(params.concat_target_taxa):
					try:
						target_codes = [l.strip() for l in open(params.concat_target_taxa).readlines() if l.strip() != '']
					except AttributeError:
						print('\n\nError: invalid "concat_target_taxa" argument. This must be a comma-separated list of any number of digits/characters to describe focal taxa (e.g. Sr_ci_S OR Am_tu), or a file with the extension .txt containing a list of complete or partial taxon codes. All sequences containing the complete/partial code will be identified as belonging to target taxa.\n\n')
			else:
				print('\nERROR: missing --concat_target_taxa argument. When concatenating, you need to give the taxonomic group (sequence prefix), groups, or a file containing a list of groups (multiple prefixes) for which to select sequences to construct a concatenated alignment\n')
				exit()

			monophyletic_clades = { }

			#Create list of relevant major/minor clades for clade grabbing
			for taxon in target_codes:
				if len(taxon) < 5 and taxon[:2] not in monophyletic_clades:
					monophyletic_clades.update({ taxon : [] })
				elif len(taxon) >= 5 and taxon[:5] not in monophyletic_clades:
					monophyletic_clades.update({ taxon[:5] : [] })

			#Grab best clades from all target groups
			seen_leaves = []
			for clade in monophyletic_clades:
				for node in tree.traverse('levelorder'):
					#If the node is big enough and not subsumed by a node we've already accepted
					if len(list(set(seen_leaves) & set([leaf.name for leaf in node]))) == 0:
						leaves = [leaf.name for leaf in node]

						#Create a record of leaves that belong to the taxonomic group
						target_leaves = set()
						for leaf in leaves[::-1]:
							if leaf[:2] == clade or leaf[:5] == clade:
								target_leaves.add(leaf[:10])
								leaves.remove(leaf)

						#If the clade is monophyletic
						if len(leaves) == 0:
							monophyletic_clades[clade].append(node)
							seen_leaves.extend([leaf.name for leaf in node])

			#Get all target taxa in the alignment
			taxa = []
			for seq in tree:
				for code in target_codes:
					if code in seq.name:
						taxa.append(seq.name[:10])
						break

			taxa = list(dict.fromkeys(taxa))

			#For each taxon, get its best sequence
			for tax in taxa:

				#Get all sequences belonging to the taxon
				taxseqs = [seq.name for seq in tree if seq.name[:10] == tax]

				score = False

				#If there's more than one sequence
				if len(taxseqs) > 1:

					clades = { }

					#Get the size of the clade in which each sequence falls (at minor clade level if available, otherwise major clade)
					if tax[:5] in monophyletic_clades:
						clades = { seq : len([leaf for clade in monophyletic_clades[tax[:5]] for leaf in clade if seq in [l.name for l in clade]]) for seq in taxseqs }
					elif tax[:2] in monophyletic_clades:
						clades = { seq : len([leaf for clade in monophyletic_clades[tax[:2]] for leaf in clade if seq in [l.name for l in clade]]) for seq in taxseqs }

					#If there's more than one sequence that falls in a robust clade
					if len(clades) > 0:

						#Filter the list of sequences to those that fall in clades
						taxseqs = [seq for seq in taxseqs if seq in clades]

						#Get the largest clade in which a sequence from the taxon falls
						best_size = max(list(clades.values()))

						#Get a list of sequences in a clade of that size
						best_seqs = [seq for seq in taxseqs]# if clades[seq] == best_size]

						#If there is only one sequence in the best-sized clade, take it and finish
						if len(best_seqs) == 1:
							seqs_per_og[file].append(og_recs[best_seqs[0]])
						#Otherwise, need to take the sequence with the best score that falls into a clade of that size
						else:
							taxseqs = best_seqs
							score = True
					#Otherwise, of all sequences that don't fall in any clade, take the one with the best score
					else:
						score = True
				#If there's only one sequence for the taxon, no problem
				elif len(taxseqs) == 1:
					seqs_per_og[file].append(og_recs[taxseqs[0]])

				#If scoring is necessary, do it on the filter set of sequences for the taxon and keep the best
				if score:
					use_cov = True
					for seq in taxseqs:
						if 'Cov' not in seq[10:]:
							use_cov = False
							break

					if use_cov:
						taxseqs = sorted(taxseqs, key = lambda x : -len(og_recs[x].seq.replace('-', '')) * float(x.split('Cov')[-1].split('_')[0]))
					else:
						taxseqs = sorted(taxseqs, key = lambda x : -len(og_recs[x].seq.replace('-', '')))

					seqs_per_og[file].append(og_recs[taxseqs[0]])

	return seqs_per_og


#Function to concatenate all the selected sequences into one alignment file
def concat(seqs_per_og, params):

	taxa = list(dict.fromkeys([rec.id[:10] for og in seqs_per_og for rec in seqs_per_og[og]]))

	seqs_per_og = { og : { rec.id[:10] : str(rec.seq) for rec in seqs_per_og[og] } for og in seqs_per_og }

	concat_seqs_per_tax = { tax : '' for tax in taxa }
	for taxon in taxa:
		for og in seqs_per_og:
			if taxon in seqs_per_og[og]:
				concat_seqs_per_tax[taxon] += seqs_per_og[og][taxon]
			else:
				print(list(seqs_per_og[og].values()))
				print(og)
				concat_seqs_per_tax[taxon] += ''.join(['-' for i in range(len(list(seqs_per_og[og].values())[0]))])

	with open(params.output + '/Output/ConcatenatedAlignment.fasta', 'w') as o:
		for tax in concat_seqs_per_tax:
			o.write('>' + tax + '\n' + concat_seqs_per_tax[tax] + '\n\n')
			

#wrapper
def run(params):
	
	if not os.path.isdir(params.output + '/Output/Guidance') or not os.path.isdir(params.output + '/Output/Trees'):
		print('\nERROR in concatenation: cannot find alignments and/or trees (looking in ' + params.output + '/Output/Guidance and ' + params.output + '/Output/Trees)')
		exit()
	else:
		seqs_per_og = remove_paralogs(params)

		concat(seqs_per_og, params)



		

























