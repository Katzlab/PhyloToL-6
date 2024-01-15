'''
#Author, date: Elinor Sterner Jan-26-2023, updated Feb-12-2023.
#Intent: To get the unique taxa from a taxonomic classification. 
#Dependencies: Python3
#Inputs: text file of taxonomies. make sure each taxonomic level is separated with '; ' (semicolon space).
#Outputs: Spreadsheet with unique taxa. If you run the script multiple times, DELETE THE PREVIOUS OUTPUT. 
#Example: python GetUniqueTaxa_v1.0.py 
'''

import os
from pathlib import Path

Path(f'unique_taxon_lists').mkdir(parents=True, exist_ok=True)#makes output folder

with open ('all_taxa.txt') as t:
	names = t.readlines()
	to_uniquify = []
	skipped = []

	for name in names: #iterate through each line of txt file


		#Adding short names to keep
		if len(name.split('; ')) <= 2:
			if ';' in name:#this removes things that are only one word (ex 'no taxID')
				to_uniquify.append(name)#add short taxonomies to list to uniquify because we don't want to remove any taxonomic information


		#for names over 2 parts
		else:
			if ' ' in name.split('; ')[-1]:#check if it has a species name
				if len(name.split('; ')) > 3:#remove sp and genus from taxonomies with over 3 parts
					short_tax = name.split('; ')[:-2]
					to_uniquify.append('; '.join(short_tax))
				
				#if short, only remove -1
				if len(name.split('; ')) < 4:#remove sp only from taxa with 3 or fewer parts
					short_tax = name.split('; ')[:-1]
					to_uniquify.append('; '.join(short_tax))

			#for names without species names
			else: #if there is no species name, remove genus only 
				if len(name.split('; ')) > 3:#remove sp and genus from taxonomies with over 3 parts
					short_tax = name.split('; ')[:-1]
					to_uniquify.append('; '.join(short_tax))
				else:
					to_uniquify.append(name)

unique_taxonomies = {}#initialize dictionary
count = 0
for i in to_uniquify:#iterate through names that have species removed
	item = i.split('; ')#divide up each word in the line of the taxonomy
	for name in item:#iterate through words in each
		count+=1#count number of words total
		unique_taxonomies.update({name.strip(): item[0]})#write to dictionary with key as word (so they automatically get uniquified) and value = MC



print(f'You have {len(unique_taxonomies)} unique taxonomies\n\n')


#write to csvs by major clade
for taxa, mc in unique_taxonomies.items():

	file = open(f'unique_taxon_lists/{mc.strip()}.csv', 'a+')#open/create a file with the name of the major clade of each taxonomy
	file.write(f'{taxa}\n')#write each taxonomy to the file named by its major clade
	file.close#close the file i dont know what happens if you dont do this but its probably bad

with open('unique_taxon_lists/all_unique_terms.csv', 'w') as o:
	for taxa, mc in unique_taxonomies.items():
		o.write(f'{taxa}\n')

	

