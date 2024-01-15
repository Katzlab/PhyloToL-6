'''
#Author, date: Uploaded by Adri Grow, 2023
#Intent: To get the NCBI taxonomic classification of organisms. 
#Dependencies: Python3, esearch
#Inputs: Spreadsheet with ten digit codes in the first column and the genus and species names in the second column.
#Outputs: Spreadsheet with taxonomy.
#Example: python GetTaxonomy_v1.0.py --input_file <path to .csv file>
'''

import os
import sys
from subprocess import check_output
from time import sleep


output_handle = 'output_taxonomies.csv'


def bad_script_call():
	
	print('\nPlease input a spreadsheet with ten digit codes in the first column and the genus and species names in the second column. Preferably, the genus and species name will be separated by a space and there will be no extraneous characters in the second column\n\n\tpython get_taxonomy.py --input_file <path to .csv file>\n')
	exit()


def get_args():
	
	input_path = ''
	
	try:
		if(sys.argv[1] == '--input_file'):
			input_path= sys.argv[2]
		else:
			bad_script_call()
	except IndexError:
		bad_script_call()
		
	if(input_path != ''):
		return input_path
	else:
		bad_script_call()
		

#Deterministic function to return the genus + species name given a line in the input spreadsheet
def get_name(line):
	name = line.split(',')[1]
	if(len(name.split(' ')) > 1):
		if(name.split(' ')[1].strip() == 'sp' or name.split(' ')[1].strip() == 'sp.'):
			name = name.split(' ')[0]	
		else:
			name = name.split(' ')[0] + ' ' + name.split(' ')[1]
	elif(len(name.split('_')) > 1):
		if(name.split('_')[1].strip() == 'sp' or name.split('_')[1].strip() == 'sp.'):
			name = name.split('_')[0]
		else:
			name = name.split('_')[0] + ' ' + name.split('_')[1]
	elif(len(name.split('-')) > 1):
		if(name.split('-')[1].strip() == 'sp' or name.split('-')[1].strip() == 'sp.'):
			name = name.split('-')[0]
		else:
			name = name.split('-')[0] + ' ' + name.split('-')[1]
			
	return name.split('\n')[0].split('\r')[0].split('(')[0]


#Returning all of the genus + species names in the input spreadsheet
def parse_taxaselection(input_path):
	taxa = {}
	for lines in open(input_path, 'r'):
		for i, line in enumerate(lines.split('\r')):
			if(len(line.split(',')) > 1 and not ('Clade' in line and 'Taxa' in line)):	
				taxa.update({ get_name(line) : line.split(',')[0].strip() })
				
	return taxa
	

#This function queries Entrez Search with the genus and species name and returns the taxonomy for each name if available
def get_taxonomy(taxa):
	
	#A rough list of most common clades for reference when multiple taxonomies are returned. Feel free to add to this, just keep the format of the names.
	clades = { 'op_' : 'opisthokont', 'op_fu_' : 'fung', 'op_me_' : 'metazoa', 'pl_' : 'archaeplastid', 'pl_gr_' : 'green alga', 'pl_rh_' : 'red alga', 'pl_gl_' : 'glaucophyt', 'sr_ap_' : 'apicomplexa', 'sr_ci_' : 'ciliat', 'sr_rh_' : 'rhizari', 'sr_di_' : 'dinoflagell', 'sr_st_' : 'stramenopil', 'ba_' : 'bacteria', 'ba_cy_' : 'cyanobacteria', 'ba_ad_' : 'acidobacteria', 'ba_pa_' : 'alphaproteobacteria', 'ba_pb_' : 'betaproteobacteria', 'ba_pg_' : 'gammaproteobacteria', 'ba_pd' : 'deltaproteobacteria', 'za_' : 'archea', 'ex_' : 'excavata', 'am_' : 'amoeb', 'am_tu_' : 'tubilinea', 'am_di_' : 'discosea', 'op_fb_' : 'fung', 'op_mb_' : 'metazoa', 'pl_rb_' : 'red alga', 'sr_ab_' : 'apicomplexa', 'sr_cb_' : 'ciliat', 'sr_rb_' : 'rhizari', 'sr_db_' : 'dinoflagell', 'sr_sb_' : 'stramenopil', 'am_tb' : 'tubilinea', 'am_db' : 'discosea' }

	taxonomies = {}
	
	for taxon in taxa:
		code = taxa[taxon].lower()
	
		print('\nFetching taxonomy for taxon ' + taxon + '\n')
		try:
			output = str(check_output('esearch -db taxonomy -query "' + taxon + '" | efetch -format xml', shell = True, executable='/bin/bash'))
			#sleep(2)
		
			taxonomy_strings = []
			
			for line in output.split('<'):
				if('lineage' in line.lower() and 'lineageex' not in line.lower()):
					if('cellular organisms' in line.split('>')[1].split('\n')[0]):
						taxonomy_strings.append(line.split('>')[1].split('\n')[0])
						
		except:
			continue
						
		
		
		hits = 0
		#If multiple taxonomies were returned
		if(len(taxonomy_strings) > 1):
			relevant_taxonomies = []
			#Check whether the relevant minor clade has a keyword
			if(code[:6] in clades):
				keyword = clades[code[:6]]
				for taxonomy in taxonomy_strings:
					#Get the taxonomy that contains the keyword
					if(keyword in taxonomy.lower()):
						hits += 1
						relevant_taxonomies.append(taxonomy)
			#If the minor clade does not have a keyword, check whether the major clade has a keyword
			elif(code[:3] in clades):
				keyword = clades[code[:3]]
				for taxonomy in taxonomy_strings:
					#Get the taxonomy that contains the keyword
					if(keyword in taxonomy.lower()):
						hits += 1
						relevant_taxonomies.append(taxonomy)
			else:
				hits = 2
				
			if(hits > 1):
				final_string = 'MULTIPLE TAXONOMIES: ' + ' | '.join(relevant_taxonomies)
			else:
				final_string = relevant_taxonomies[0]
		elif(len(taxonomy_strings) == 0):
			final_string = 'NO TAXONOMY RETURNED'
		else:
			final_string = taxonomy_strings[0]
								
		print(final_string)
						
		taxonomies.update({ taxon : final_string })
    
	return taxonomies
	

#Writing the output spreadsheet with taxonomies using the input spreadsheet as a template
def write_spreadsheet(taxonomies, input_path):
	
	with open(output_handle, 'w') as o:
		for lines in open(input_path, 'r'):
			for i, line in enumerate(lines.split('\r')):
				if(len(line.split(',')) > 1 and not ('Clade' in line and 'Taxa' in line)):	
					if(get_name(line) in taxonomies):
						line_split = line.strip().split(',')
						line_split.append(taxonomies[get_name(line)])
						
						o.write(','.join(line_split) + '\n')
					else:
						o.write(line + '\n')
			

#This is a one-time-use function to add taxonomy to a specific spreadsheet (R2G and BB counts)
def add_to_seq_count_spreadsheet():

	used = []
	
	with open('seq_counts_taxaselection.csv', 'w') as o:
		for lines in open(output_handle, 'r'):
			for count_lines in open('seq_counts_selection.csv', 'r'):
				for count_line in count_lines.split('\r'):
					for i, line in enumerate(lines.split('\r')):
						if(count_line.split(',')[0] in line and count_line.split(',')[0] not in used):
							used.append(count_line.split(',')[0])
							o.write(count_line + ',' + line.split(',')[1] + ',' + line.split(',')[6] + '\n')
    		
	
def main():

	input_path = get_args()

	taxa = parse_taxaselection(input_path)
			
	taxonomies = get_taxonomy(taxa)
	
	write_spreadsheet(taxonomies, input_path)
	
	#You will probably never need to use this function
	#add_to_seq_count_spreadsheet()

	
main()
