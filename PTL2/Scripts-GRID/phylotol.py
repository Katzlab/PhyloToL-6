#!/usr/bin/python3
import os, sys, re

import utils
import preguidance
import guidance_treeBuilding
from logger import Logger

if __name__ == '__main__':

	params = utils.get_params()
	
	Logger.Message('Cleaning up existing files and organizing output folder', Logger.BOLD)
	utils.clean_up(params)

	if params.start == 'raw':
		Logger.Message('Running preguidance', Logger.BOLD)
		preguidance.run(params)
	
	if params.start in ('unaligned', 'raw') and params.end in ('aligned', 'trees'):
		Logger.Message('Running guidance', Logger.BOLD)
		guidance_treeBuilding.run(params)

	#ADD HERE: If only running tree-building...

	if params.contamination_loop != None:
		Logger.Message('Running contamination loop', Logger.BOLD)

	if not params.keep_temp:
		os.system('rm -r ' + params.output + '/Output/Temp')
	
