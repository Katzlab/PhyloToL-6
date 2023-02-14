import os, sys, re
import ete3
from logger import Logger


def get_best_clade(params, tree):

	if params.target_taxa_file != None:
		try:
			target_taxa_list = [line.strip() for line in open(PathtoFiles + '/' + target_taxa_list)]
		except (FileNotFoundError, TypeError):
			Logger.Error('The --target_taxa_file could not be found or was incorrectly formatted.')

	if params.at_least_file != None:
		try:
			at_least_list = [line.strip() for line in open(PathtoFiles + '/' + at_least_file)]
		except (FileNotFoundError, TypeError):
			Logger.Error('The --at_least_file could not be found or was incorrectly formatted.')
	else:
		at_least_list = []


	### FROM HERE BELOW IN THIS FUNCTION NEEDS REPLACING WITH ETE3
		
	forbidden_nodes = [node for node in nodes_to_exclude]
	for node in nodes_to_exclude:
		for num in tree.getNodeNumsAbove(node):
			forbidden_nodes.append(tree.node(num))
	
	best_node = None
	best_size = 0
	for node in tree.iterNodesNoRoot():
		if(node not in forbidden_nodes):
			leaves = tree.getAllLeafNames(node)
							
			num = 0.0; dem = 0.0;
		
			non_minor = []
			for leaf in leaves:
				if(leaf[:2] != target_minor and leaf[:4] != target_minor):
					num += 1.0;
					non_minor.append(leaf[:10])
					
			if(target_taxa_list == 'na' or target_taxa_list == '' or target_taxa_list == 'NA'):
				n_targets = len(list(dict.fromkeys([tip[:10] for tip in leaves if(tip[:2] == target_clade or tip[:3] == target_clade or tip[:4] == target_clade or tip[:5] == target_clade or tip[:7] == target_clade or tip[:8] == target_clade)])))
			else:
				n_targets = len(list(dict.fromkeys([tip[:10] for tip in leaves if((tip[:2] == target_clade or tip[:3] == target_clade or tip[:4] == target_clade or tip[:5] == target_clade or tip[:7] == target_clade or tip[:8] == target_clade) and (tip[:10] in target_taxa_list or tip[:8] in target_taxa_list))])))
							
			at_least_taxa = len(list(dict.fromkeys([leaf[:10] for leaf in tree.getAllLeafNames(node) if leaf[:10] in at_least_list])))

			if(num <= cont_num_contams and n_targets > best_size and n_targets >= min_presence and at_least_taxa >= num_at_least):
				best_node = node
				best_size = n_targets
			
	return best_node