#Author, date: Auden Cote-L'Heureux & Godwin Ani, last updated Nov 27th 2023
#Motivation: Visualize placement of taxa by taxonomic group in trees
#Intent: Color tip labels in trees by taxonomic group
#Dependencies: Python3, ete3
#Inputs: A folder of trees
#Outputs: a folder of colored trees
#Example: python ColorByClade_v2.1.py -i /path/to/trees


import os, sys
import ete3
import argparse


#Needed for communicating with Figtree program
figtree_format = '''begin figtree;
	set appearance.backgroundColorAttribute="Default";
	set appearance.backgroundColour=#ffffff;
	set appearance.branchColorAttribute="User selection";
	set appearance.branchColorGradient=false;
	set appearance.branchLineWidth=1.0;
	set appearance.branchMinLineWidth=0.0;
	set appearance.branchWidthAttribute="Fixed";
	set appearance.foregroundColour=#000000;
	set appearance.hilightingGradient=false;
	set appearance.selectionColour=#2d3680;
	set branchLabels.colorAttribute="User selection";
	set branchLabels.fontName="sansserif";
	set branchLabels.fontSize=8;
	set branchLabels.fontStyle=0;
	set branchLabels.isShown=false;
	set branchLabels.significantDigits=4;
	set layout.expansion=0;
	set layout.layoutType="RECTILINEAR";
	set layout.zoom=0;
	set legend.attribute=null;
	set legend.fontSize=10.0;
	set legend.isShown=false;
	set legend.significantDigits=4;
	set nodeBars.barWidth=4.0;
	set nodeBars.displayAttribute=null;
	set nodeBars.isShown=false;
	set nodeLabels.colorAttribute="User selection";
	set nodeLabels.displayAttribute="Node ages";
	set nodeLabels.fontName="sansserif";
	set nodeLabels.fontSize=8;
	set nodeLabels.fontStyle=0;
	set nodeLabels.isShown=false;
	set nodeLabels.significantDigits=4;
	set nodeShape.colourAttribute="User selection";
	set nodeShape.isShown=false;
	set nodeShape.minSize=10.0;
	set nodeShape.scaleType=Width;
	set nodeShape.shapeType=Circle;
	set nodeShape.size=4.0;
	set nodeShape.sizeAttribute="Fixed";
	set polarLayout.alignTipLabels=false;
	set polarLayout.angularRange=0;
	set polarLayout.rootAngle=0;
	set polarLayout.rootLength=100;
	set polarLayout.showRoot=true;
	set radialLayout.spread=0.0;
	set rectilinearLayout.alignTipLabels=false;
	set rectilinearLayout.curvature=0;
	set rectilinearLayout.rootLength=100;
	set scale.offsetAge=0.0;
	set scale.rootAge=1.0;
	set scale.scaleFactor=1.0;
	set scale.scaleRoot=false;
	set scaleAxis.automaticScale=true;
	set scaleAxis.fontSize=8.0;
	set scaleAxis.isShown=false;
	set scaleAxis.lineWidth=1.0;
	set scaleAxis.majorTicks=1.0;
	set scaleAxis.origin=0.0;
	set scaleAxis.reverseAxis=false;
	set scaleAxis.showGrid=true;
	set scaleBar.automaticScale=true;
	set scaleBar.fontSize=10.0;
	set scaleBar.isShown=true;
	set scaleBar.lineWidth=1.0;
	set scaleBar.scaleRange=0.0;
	set tipLabels.colorAttribute="User selection";
	set tipLabels.displayAttribute="Names";
	set tipLabels.fontName="sansserif";
	set tipLabels.fontSize=12;
	set tipLabels.fontStyle=0;
	set tipLabels.isShown=true;
	set tipLabels.significantDigits=4;
	set trees.order=false;
	set trees.orderType="increasing";
	set trees.rooting=false;
	set trees.rootingType="User Selection";
	set trees.transform=false;
	set trees.transformType="cladogram";
end;'''


def get_args():

	parser = argparse.ArgumentParser(
		prog = 'Tree-coloring script, Version 2.1',
		description = "Updated Nov 27th, 2023"
	)

	parser.add_argument('-i', '--input', type = str, required = True, help = 'Path to a folder containing input trees (which must have the file extension .tre, .tree, .treefile, or .nex)')
	parser.add_argument('-k', '--keyfile', type = str, help = 'Path to a text file with two tab-separated columns; the first a set of keys and the second a color for each key in hex-code format. Any sequence starting with a particular key will be assigned the color corresponding to that key in this file.')
	
	return parser.parse_args()


#Function to extract newick string from either newick or nexus file
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


def fix_node_labels(newick):

	out_newick = ''
	for chunk in newick.split(':'):
		if ';' in chunk:
			out_newick += chunk
		elif ')' in chunk:
			out_newick += chunk.split(')')[0] + ')[&NULL_LABEL=' + chunk.split(')')[-1] + ']:'
		else:
			out_newick += chunk + ':'

	return out_newick


def write_lines(o, newick, taxa_and_colors, tree_font_size):
	ntax = str(len(taxa_and_colors))

	newick = fix_node_labels(newick)
	
	#writes the header to the tree file		
	o.write('#NEXUS\n')	
	o.write('begin taxa;\n')
	o.write('\tdimensions ntax=' + ntax + ';\n')
	o.write('\ttaxlabels\n')

	#write out all taxa 	
	for taxon in taxa_and_colors:
		o.write('\t' + taxon + '\n')
			
	o.write(';\nend;\n\n')
			
	o.write('begin trees;\n')
	o.write('\ttree tree_1 = [&R]\n')
	o.write(newick)
	o.write('end;\n\n')
			
	
	for line in figtree_format:
		if('.fontSize' in line):
			o.write(line.replace('8', tree_font_size))
		else:
			o.write(line)


def write_nexus(newick, leaf_colors, params):
		
	with open(out_path, 'w') as o:
		write_lines(o, newick, taxa_and_colors, tree_font_size)
	

def tree_formatting_wrapper(file):
	try:
		newick = get_newick(file)
		tree = ete3.Tree(newick)
		
		majs = list(dict.fromkeys([leaf.name[:2] for leaf in tree]))
		#Only try to reroot trees with more than 2 major clades. This was added to fix the ETE3 "Cannot set myself as outgroup" error
		
		if len(majs) > 2: 
			tree = reroot(tree)	
		tree.ladderize(direction = 1)
		tree.write(outfile = 'ColoredTrees/' + file.split('/')[-1].split('.tree')[0] + '_Colored.tree')
	except Exception as e:
		print(f" {file.split('/')[-1]} has {e} error ")


	
def color(file, args):

	if args.keyfile != None:
		if os.path.isfile(args.keyfile):
			try:
				colors = { line.split('\t')[0] : line.split('\t')[1].strip() for line in open(args.keyfile) if len(line.split('\t')) == 2 }
			except:
				print('\nERROR: your keyfile is incorrectly formatted\n')
				exit()
		else:
			print('\nERROR: your input keyfile could not be found\n')
	else:
		colors = { 'Ba' : '[&!color=#000000]', 'Za' : '[&!color=#808080]', 'Sr' : '[&!color=#B4A26D]', 'Op' : '[&!color=#1260CC]', 'Pl' : '[&!color=#026736]', 'Ex' : '[&!color=#E63B60]', 'EE' : '[&!color=#C76A6A]', 'Am' : '[&!color=#29C5F6]', 'EE_cr' : '[&!color=#08B461]', 'EE_ha' : '[&!color=#03EA74]', 'Sr_ci' : '[&!color=#A97533]', 'Sr_ap' : '[&!color=#D4BA99]', 'Sr_rh' : '[&!color=#8A3324]', 'Sr_st' : '[&!color=#E97451]', 'Sr_di' : '[&!color=#492815]' }
		
	newick = get_newick(file)

	tree = ete3.Tree(newick)

	leaf_colors = []
	for leaf in tree:
		keys = sorted([key for key in colors if leaf.name.startswith(key)], key = lambda x : -len(x))

		# the line below allows you to have keys anywhere within name and not just start of name.. to use, you have to # the line above
		#keys = sorted([key for key in colors if key in leaf.name], key=lambda x: -len(x))
		
		if len(keys) > 0:
			if '[&!color=' in colors[keys[0]]:
				leaf_colors.append(leaf.name + colors[keys[0]])
			else:
				leaf_colors.append(leaf.name + '[&!color=' + colors[keys[0]] + ']')
		else:
			leaf_colors.append(leaf.name)
	
	with open('ColoredTrees/' + file.split('/')[-1].split('.tree')[0] + '.tree', 'w') as o:
		write_lines(o, newick, leaf_colors, str(12))#change tree font size here (right now it is 12)

if __name__ == '__main__':

	args = get_args()

	if not os.path.isdir('ColoredTrees'):
		os.mkdir('ColoredTrees')

	for tree in os.listdir(args.input):
		if tree.split('.')[-1] in ('tree', 'tre', 'treefile', 'nex'):
			tree_formatting_wrapper(args.input + '/' + tree)
        
	for tree in os.listdir('ColoredTrees'):
		if tree.split('.')[-1] in ('tree', 'tre', 'treefile', 'nex'):
			color('ColoredTrees/' + tree, args)
        














