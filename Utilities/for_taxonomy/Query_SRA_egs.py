'''
#Author, date: Elinor Sterner, Feb 2023
#Intent: To grab recent assemblies (since 2020) and GCA codes. 
#Dependencies: Python3, Biopython
#Inputs: Folder named 'unique_taxon_lists' with files of keywords by major clade (separated by new lines).
#Outputs: File of species, IDs, and GCA or SRR codes AND a file with uniquified codes.
#Example: python Query_SRA_egs.py -t (transcriptome, SRA db) or -g (genome, assembly db)
'''

from Bio import Entrez
from Bio import SeqIO
import os
import sys


def get_args():
	
	Entrez.email = "@smith.edu"#CHANGE UR EMAIL
	Entrez.tool = "Biopython_NCBI_Entrez_downloads.ipynb"

	if len(sys.argv) < 2:
		print(f'enter -t or -g in command line to choose genomes (-g) or transcriptomes (-t)')
	if '-t' in sys.argv:
		data_type = False
	elif '-g' in sys.argv:
		data_type = True

	with open('RecentIDs.csv', 'w') as o:#starts output file and writes header
		o.write('major clade, keyword, species, ID, experiment, sequencing technology, GCA/SRR,\n')

	get_keywords(data_type)


def get_keywords(data_type):

	for file in os.listdir('unique_taxon_lists'):
		if file.endswith('_unique.csv'):#put name of file to look at here. or only .csv to look at all of them
			with open(f'unique_taxon_lists/{file}', 'r') as lines:#read each file
				mc = file.split("_unique.csv")[0]
				print(f'Searching taxonomic names in {mc}\n\n')
				
				for line in lines.readlines():#iterate file
					keyword = line.strip()#keyword for genbank search is each word in the files
					
					if data_type == False:
						fetch_SRA(mc, keyword)
					if data_type == True:
						fetch_CDS(mc, keyword)

	write_unique_codes()


def fetch_CDS(mc, keyword):#searches your keywords in the assembly database
	print(f'\nGrabbing recent CDSs')
	all_stuff = []#initiate list, will put genbank codes into this

	#get IDs of assemblies for keyword since 2020. Returns multiple IDs
	handle = Entrez.esearch(db="assembly", term=keyword + "[Organism:exp]" + "2020 [SeqReleaseDate]:3000", retmax=100)
	id_record = Entrez.read(handle)
	print(f'There are {len(id_record["IdList"])} assemblies labeled as {keyword} in genbank since 2020\nFetching IDs and GCAs\n')

	#Iterate through list of IDs given above, seach for their associated GCAs. Only one GCA for each ID, and each corresponds to 1 individual sequenced
	for tax_id in id_record['IdList']:
		handle = Entrez.esummary(db="assembly", id=tax_id, retmode="text")
		gca_records = Entrez.read(handle, validate=False)
		handle.close()
		
		#parse the output (its really awful. pythonic turduken: dict(list(str(dict))) type deal)
		for record in gca_records['DocumentSummarySet']['DocumentSummary']:
			sp = record['Organism']
			gca=record['AssemblyAccession']
			stuff = f'{mc}, {keyword}, {sp},{tax_id}, , ,{gca}'
			all_stuff.append(stuff)
	
	write_to_csv(all_stuff)#send this new info to be added to output sheet


def fetch_SRA(mc, keyword):#searches your keywords in the SRA db
	
	all_stuff = []
	print(f'\nGrabbing recent SRRs')
	# get IDs from taxonomies
	handle = Entrez.esearch(db="sra", term=keyword + "[Organism:exp]"+ " 2020:2023[PDAT]", retmax=100)
	id_record = Entrez.read(handle, validate = False)
	print(f'There are {len(id_record["IdList"])} SRAs labeled as {keyword} in genbank since 2020\nFetching SRAs\n')


	#get SRRs for taxonomy
	for rec in id_record['IdList']:#iterates through all of the IDs for the taxonomy
		tax_id = rec
		handle = Entrez.esummary(db="sra", id=tax_id)
		srr_records = Entrez.read(handle)#parse genbank info

		#parse out all information needed from genbank info
		sp = srr_records[0]['ExpXml'].split('ScientificName="')[1].split('"')[0]#extract species from genbank info
		srr = srr_records[0]['Runs'].split('"')[1]#extract srr from genbank info
		seq_type = srr_records[0]['ExpXml'].split('<LIBRARY_STRATEGY>')[1].split('</LIBRARY_STRATEGY>')[0]#parse to "library_strategy" parameter to check if its amplicon
		machine = srr_records[0]['ExpXml'].split('<Platform instrument_model="')[1].split('">')[0]#get the type of sequencing machine used
		if 'AMPLICON' not in seq_type:
			stuff = f'{mc}, {keyword}, {sp}, {tax_id}, {seq_type}, {machine}, {srr},'#write to comma separated string
			all_stuff.append(stuff)


	write_to_csv(all_stuff)

def write_to_csv(data):#writes the output from fetch_SRA or fetch_CDS to a csv
	with open('RecentIDs.csv', 'a+') as o:
		for i in data:
			o.write(f'{i}\n')

def write_unique_codes():#uniquify the list of IDs that the scipt grabbed. since we are searching all taxanomic levels, we query many repeats so this removes them.

	#Writing unique files
	with open('RecentIDs.csv', 'r') as o:# read file of all data
		taxa = o.readlines()
		print(f'\nThere are {len(taxa)} codes before uniquifying\n\n')
		unique_lines = {line.split(', ')[-1] : line.split(', ')[0:-1] for line in taxa}#makes dictionary of SRR/GCA:other info to uniquify the codes
		print(f'\nYou have {len(unique_lines)} unique codes... writing them to unique_taxa.csv')

	with open ('unique_taxa.csv', 'w') as o:#start csv of unique codes
		for gca, other in unique_lines.items():#parse uniquified dictionary
			o.write(f'{(", ").join(other)}, {gca}')#write out (use join to convert the list containing other info to a string)



if __name__ == '__main__':
	get_args()
