#Author, date: ACL June 8 2023
#Motivation: Get record of OG presence across taxa from ReadyToGo files
#Intent: Create a spreadsheet summarizing OG presence
#Inputs: A folder of ReadyToGo files
#Outputs: Spreadsheet
#Example: Python SharedOGs_v1.0.py ReadyToGo_AA


import os, sys
from Bio import SeqIO
from tqdm import tqdm


input_dir = sys.argv[1]

print('\nCreating a record of taxa per OG...')

taxa_by_og = { }
for file in tqdm(os.listdir(input_dir)):
	if file.split('.')[-1] in ('fasta', 'faa', 'fna', 'fa'):
		tax = file[:10]
		for rec in SeqIO.parse(input_dir + '/' + file, 'fasta'):
			if rec.id[-10:] not in taxa_by_og:
				taxa_by_og.update({ rec.id[-10:] : [] })

			taxa_by_og[rec.id[-10:]].append(tax)


print('\nWriting output file...')

all_taxa = sorted(list(dict.fromkeys([tax for og in taxa_by_og for tax in taxa_by_og[og]])))
all_maj = sorted(list(dict.fromkeys([tax[:2] for og in taxa_by_og for tax in taxa_by_og[og]])))
with open('OGSharedness.csv', 'w') as o:
	o.write('OG,Sequences,Species,Paralogness,MinorClades,MajorClades,' + ','.join(all_maj) + ',' + ','.join(all_taxa) + '\n')
	for og in tqdm(taxa_by_og):

		og_majs = list(dict.fromkeys([tax[:2] for tax in taxa_by_og[og]]))
		og_taxa = list(dict.fromkeys(taxa_by_og[og]))
		
		o.write(og + ',' + str(len(taxa_by_og[og])) + ',' + str(len(list(dict.fromkeys(taxa_by_og[og])))) + ',' + str(len(taxa_by_og[og])/len(list(dict.fromkeys(taxa_by_og[og])))) + ',' + str(len(list(dict.fromkeys([tax[:5] for tax in taxa_by_og[og]])))) + ',' + str(len(list(dict.fromkeys([tax[:2] for tax in taxa_by_og[og]])))))
		for maj in all_maj:
			if maj in og_majs:
				o.write(',' + str(len([tax for tax in og_taxa if tax[:2] == maj])))
			else:
				o.write(',0')

		for tax in all_taxa:
			if tax in taxa_by_og[og]:
				o.write(',' + str(taxa_by_og[og].count(tax)))
			else:
				o.write(',0')
		o.write('\n')




