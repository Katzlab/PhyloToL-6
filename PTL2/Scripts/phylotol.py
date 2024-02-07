# Last updated: Jan 2024
# Author: Auden Cote-L'Heureux

# This script is what users should call when running any or all components of
# PhyloToL 6 part 2. It briefly determines which parts of the pipeline should be
# run (pre-Guidance, Guidance, tree building, contamination loop, and/or
# concatenation) based on the --start and --end parameters, and then runs all 
# of these components. Each component is actually run by the run() function in 
# that component's respective script, named accordingly, which are loaded into 
# this script as Python libraries at the top.

#!/usr/bin/python3
import os, sys, re
import contamination
import utils
import preguidance
import guidance
import trees
import concatenate


if __name__ == '__main__':

	#Reading in all input parameters
	params = utils.get_params()

	#First, checking for and cleaning up existing outputs from previous runs,
	#then setting up an empty output folder structure for the new run
	if not (params.concatenate and params.start == 'trees'):
		print('\nCleaning up existing files and organizing output folder\n')
		utils.clean_up(params)

	#Running pre-Guidance (preguidance.py)
	if params.start == 'raw':
		print('\nRunning preguidance\n')
		preguidance.run(params)

	#Running Guidance (guidance.py)
	if params.start in ('unaligned', 'raw') and params.end in ('aligned', 'trees'):
		print('\nRunning guidance\n')
		guidance.run(params)

	#Building trees (trees.py)
	if params.start != 'trees' and params.end == 'trees':
		print('\nBuilding trees\n')
		trees.run(params)

	#Running the contamination loop (contamination.py)
	if params.contamination_loop != None:
		print('\nRunning contamination loop\n')
		contamination.run(params)

	#Running concatenation (concatenate.py)
	if params.concatenate:
		print('\nChoosing orthologs and concatenating alignments...\n')
		concatenate.run(params)
	
