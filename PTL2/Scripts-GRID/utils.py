import os, sys, re
import argparse
from logger import Logger
import shutil


def get_params():

	parser = argparse.ArgumentParser(
                    prog = 'PhyloToL v6.0',
                    description = "Updated December 9, 2022 by Auden Cote-L'Heureux. Link to GitHub: https://github.com/AudenCote/PhyloToL_v6.0"
                    )

	common = parser.add_argument_group('Commonly adjusted parameters')
	common.add_argument('--start', default = 'raw', choices = {'raw', 'unaligned', 'aligned', 'trees'}, help = 'Stage at which to start running PhyloToL.')
	common.add_argument('--end', default = 'trees', choices = {'unaligned', 'aligned', 'trees'}, help = 'Stage until which to run PhyloToL. Options are "unaligned" (which will run up to but not including guidance), "aligned" (which will run up to but not including RAxML), and "trees" which will run through RAxML')
	common.add_argument('--gf_list', default = None, help = 'Path to the file with the GFs of interest. Only required if starting from the raw dataset.')
	common.add_argument('--taxon_list', default = None, help = 'Path to the file with the taxa (10-digit codes) to include in the output.')
	common.add_argument('--data', help = 'Path to the input dataset. The format of this varies depending on your --start parameter')
	common.add_argument('--output', default = '../', help = 'Directory where the output folder should be created. If not given, the folder will be created in the parent directory of the folder containing the scripts.')
	common.add_argument('--force', action = 'store_true', help = 'Overwrite all existing files in the "Output" folder.')
	common.add_argument('--tree_method', default = 'iqtree', choices = {'iqtree, raxml, all'}, help = 'Program to use for tree-building')

	core = parser.add_argument_group('Core parameters (rarely altered from the defaults)')
	core.add_argument('--blast_cutoff', default = 1e-20, type = float, help = 'Blast e-value cutoff')
	core.add_argument('--len_cutoff', default = 10, type = int, help = 'Amino acid length cutoff for removal of very short sequences after column removal in Guidance.')
	core.add_argument('--sim_cutoff', default = 0.98, type = float, help = 'Sequences from the same taxa that are assigned to the same OG are removed if they are more similar than this cutoff')
	core.add_argument('--overlap_cutoff', default = 0.35, type = float, help = 'A sequence is removed if its alignment length to the longest sequence in its OG & taxon is less than this proportion of the length of the longest sequence')
	core.add_argument('--guidance_iters', default = 5, type = int, help = 'Number of Guidance iterations for sequence removal')
	core.add_argument('--seq_cutoff', default = 0.3, type = float, help = 'During guidance, taxa are removed if their score is below this cutoff')
	core.add_argument('--col_cutoff', default = 0.4, type = float, help = 'During guidance, columns are removed if their score is below this cutoff')
	core.add_argument('--res_cutoff', default = 0.0, type = float, help = 'During guidance, residues are removed if their score is below this cutoff')
	core.add_argument('--guidance_threads', default = 20, type = int, help = 'Number of threads to allocate to Guidance')

	CL = parser.add_argument_group('Contamination loop parameters')
	CL.add_argument('--contamination_loop', default = 'clade', choices = {'seq', 'clade'}, help = 'Remove sequences by looking at the sisters of each sequence in a rules file or by picking the best clades')
	CL.add_argument('--nloops', default = 5, type = int, help = 'The maximum number of contamination-removal loops')
	
	CL.add_argument('--sister_rules', default = None, help = 'Path to a file of rules, nly used in "seq" mode. Sequences in the rules file with specified contaminants will be removed if sister only to those contaminants')
	CL.add_argument('--remove_nm_sister', action = 'store_true', help = 'Only used in "seq" mode. When set to true, sequences sister to a non-monophyletic group of taxa that all fall into categories designated as contaminants in the rules file will be removed')

	CL.add_argument('--target_clade', nargs = '+', default = None, help = 'Only used in "clade" mode. Selected clades can have no more than num_contams (below) sequences that are not of this clade (can be 2, 4, 5, 7, 8, or 10 digits)')
	CL.add_argument('--num_contams', type = int, default = 2, help = 'Only used in "clade" mode. Selected clades can have no more than this number of sequences that are not of the target clade')
	CL.add_argument('--min_target_presence', type = int, default = 8, help = 'Only used in "clade" mode. The minimum number of species belonging to a target clade allowed in a selected clade')
	CL.add_argument('--target_taxa_file', default = None, help = 'Only used in "clade" mode. Path to a text file containing 2, 4, 5, 7, 8, or 10 digit codes that are your target taxa; you must have at least min_target_presence (above) of ten-digit codes that match these criteria in your selected clade')
	CL.add_argument('--at_least_file', default = None, help = 'Only used in "clade" mode. A file containing 2, 4, 5, 7, 8, or 10 digit codes; any selected clade must have at least at_least_sisters_num of taxa that match these criteria; this is used to require the presence of certain sister lineages')
	CL.add_argument('--at_least_sisters_num', type = int, default = 1, help = 'Only used in "clade" mode. See above.')
	CL.add_argument('--coverage_diff_OM', type = int, default = None, help = 'The orders of magnitude by which the coverages of two sequences in the same clade can differ for the seq with lower coverage to be removed (default is to NOT filter for differential coverage)')
	CL.add_argument('--coverage_diff_abs', type = int, default = None, help = 'The absolute number by which the coverages of two sequences in the same clade can differ for the seq with lower coverage to be removed (default is to NOT filter for differential coverage)')

	sisters = parser.add_argument_group('Parameters for the sister report')
	CL.add_argument('--query_clades', nargs = '+', default = None, help = 'A list of 2, 4, 5, 7, 8, or 10 digit codes specifying which taxa for which to count sisters, separated by a comma. Alternatively, input a file containing a list of 10 digit codes of taxa for which to return sisters if there are a lot')
	CL.add_argument('--sister_clades', nargs = '+', help = 'A list of 2, 4, 5, 7, 8, or 10 digit codes specifying which taxa which sisters to represent in the spreadsheet, separated by a comma. Alternatively, input a file containing a list of 10 digit codes of taxa for sisters to represent if there are a lot')
	CL.add_argument('--break_up', nargs = '+', default = None, help = 'A list of major clades for which to break up the sister report into minor clades')
	CL.add_argument('--branch_length_filter', default = 'avg', choices = {'avg', int, float}, help = 'Filter tips to represent by branch length')
	CL.add_argument('--single_sister_only', action = 'store_true', help = 'Whether or not you only want to report sister relationships when there is only a single taxon sister to a sequence')


	other = parser.add_argument_group('Other arguments')
	other.add_argument('--concatenate', action = 'store_true', help = 'Remove paralogs and generate an alignment for concatenation')
	other.add_argument('--tree_font_size', default = 12, help = "Change this if you're not quite happy with the font size in the output trees. If you want smaller font in your trees, you can lower this value; and if you want larger font in your trees, you can raise this value. Some common values are 8, 10, and 12. Size 16 font is pretty big, and size 4 font is probably too small for most purposes. Iconoclasts use size 9, 11, or 13 font.")
	other.add_argument('--keep_temp', action = 'store_true', help = "Change this if you're not quite happy with the font size in the output trees. If you want smaller font in your trees, you can lower this value; and if you want larger font in your trees, you can raise this value. Some common values are 8, 10, and 12. Size 16 font is pretty big, and size 4 font is probably too small for most purposes. Iconoclasts use size 9, 11, or 13 font.")

	return parser.parse_args()


def clean_up(params):

	if not os.path.isdir(params.output + '/Output'):
		os.mkdir(params.output + '/Output')
	elif os.path.isdir(params.output + '/Output') and params.force == False:
		Logger.Error('An "Output" folder already exists at the given path. Please delete or rename this folder and try again.')
	elif params.force and len([d for d in os.listdir(params.output + '/Output') if d != '.DS_Store']) > 0:
		Logger.Warning('An "Output" folder already exists at the given path, but all contents were deleted in --force mode.')
		os.system('rm -r ' + params.output + '/Output/*')

	os.mkdir(params.output + '/Output/Temp')

	def copy_input():
		if os.path.isdir(params.data):
			input_files = [f for f in os.listdir(params.data) if f.endswith('.faa') or f.endswith('.fasta') or f.endswith('.fa')]
			if len(input_files) > 0:
				for f in input_files:
					if(params.start == 'unaligned'):
						shutil.copyfile(params.data + '/' + f, params.output + '/Output/Pre-Guidance/' + f)
					elif(params.start == 'aligned'):
						shutil.copyfile(params.data + '/' + f, params.output + '/Output/Guidance/' + f)
					elif(params.start == 'trees'):
						shutil.copyfile(params.data + '/' + f, params.output + '/Output/Trees/' + f)
			else:
				Logger.Warning('The given path to a folder of ' + params.start.strip('s') + ' files was located, but no ' + params.start.strip('s') + ' files were found. Make sure the file extensions are .fasta, .fa, or .faa.')
		else:
			Logger.Error('Input ' + params.start.strip('s') + ' data files not found. Please make sure that the given path (--data) is correct or set --start to "raw".')
	
	os.mkdir(params.output + '/Output/Pre-Guidance')
	if params.start == 'unaligned':
		copy_input()

	if params.start in ('unaligned', 'aligned') or params.end in ('aligned', 'trees', None):
		os.mkdir(params.output + '/Output/Guidance')
		if(params.start == 'aligned'):
			copy_input()

	if params.end == 'trees' or params.contamination_loop != None:
		os.mkdir(params.output + '/Output/Trees')
		if(params.start == 'trees'):
			copy_input()

	



































