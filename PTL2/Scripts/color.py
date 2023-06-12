import os, sys
import ete3

	
def get_newick(fname):
	
	newick = ''
	for line in open(fname):
		line = line.split(' ')[-1]
		if(line.startswith('(') or line.startswith('tree1=')):
			newick = line.split('tree1=')[-1].replace("'", '').replace('\\', '')

	return newick


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


def write_lines(o, newick, taxa_and_colors, tree_font_size):
	ntax = str(len(taxa_and_colors))
			
	o.write('#NEXUS\n')	
	o.write('begin taxa;\n')
	o.write('\tdimensions ntax=' + ntax + ';\n')
	o.write('\ttaxlabels\n')
			
	for taxon in taxa_and_colors:
		o.write('\t' + taxon + '\n')
			
	o.write(';\nend;\n\n')
			
	o.write('begin trees;\n')
	o.write('\ttree tree_1 = [&R]\n')
	o.write(newick)
	o.write('end;\n\n')
			
	with open('figtree_format.txt', 'r') as ff:
		for line in ff:
			if('.fontSize' in line):
				o.write(line.replace('8', tree_font_size))
			else:
				o.write(line)


def write_nexus(newick, leaf_colors, params):
		
	with open(out_path, 'w') as o:
		write_lines(o, newick, taxa_and_colors, tree_font_size)
	

def color(file, params):

	colors = { 'Ba' : '[&!color=#000000]', 'Za' : '[&!color=#808080]', 'Sr' : '[&!color=#7b2516]', 'Op' : '[&!color=#12aaff]', 'Pl' : '[&!color=#006300]', 'Ex' : '[&!color=#ffa100]', 'EE' : '[&!color=#ff6288]', 'Am' : '[&!color=#aa00ff]' }

	newick = get_newick(file)	
	tree = ete3.Tree(newick)
	tree = reroot(tree)
	tree.ladderize()

	leaf_colors = [leaf + colors[leaf[:2]] for leaf in tree]

	with open(params.output + '/ColoredTrees/' + file.split('.tree')[0] + '_Colored.tree', 'w') as o:
		write_lines(o, newick, leaf_colors, params.tree_font_size)
				
