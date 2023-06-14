![Lab image](https://github.com/Katzlab/PhyloToL-6/blob/621517812a5ed4256f19ba26498a665bdd5060af/Katzlab.png)
# Utilities
> This folder contains many useful tools for analyzing sequence data.

## For Taxonomies dir:
### Query_SRA_egs.py:

**Purpose** Gives a spreadsheet of all SRA codes or GCA codes from Genbank associated with taxonomic terms in an input csv file

**Input** Input is folder 'unique_taxon_lists' with files of keywords by major clade (separated by new lines) by get_unique_taxa.py or manually 

**Output** all SRAs or GCA since 2020 (can be adjusted by modifying script). For SRAs, the script also gives sequecing technology used (pacbio, miseq, etc) and experiment type. It excludes all SRAs that include the word 'amplicon'.

**Usage** -t (transcriptome, searches SRA db) or -g (genome, searches assembly db) in the command line to specify data type. 
>  Example command line: `python Query_SRA_egs.py -t OR -g`

### get_unique_taxa.py: 
Written by Elinor 1/26, updated 2/12

**Purpose** make lists of unique taxonomy from phylotol master taxonomy column (genbank taxonomy for each taxa in the pipeline). These lists are the intended input for Query_SRA_egs.py. This cuts off the genus (and species if there is one), uniquifies the list and writes them out to files by the first word of the taxonomy

**Input** text file of taxonomies called `all_taxa.txt`. make sure each taxonomic level is separated with `; ` (semicolon space) or the script will not parse the names right

**Output** txt file of _all_ unique names found, and a directory of txt files of unique names sorted by major clade (the first word in the line of input taxonomy

**Usage**
>`python get_unique_taxa.py`

WARNING: if you run the script multiple times, DELETE THE PREVIOUS OUTPUT. this is because it appends lines to the 
end of files so you will have many duplicates

### get_taxonomy.py:

**Purpose** 

Queries Entrez Search with the genus and species name associated with 10 digit codes and returns the taxonomy for each name if available.

**Input**

Spreadsheet with ten digit codes in the first column and the genus and species names in the second column (csv).

**Output**

CSV file called `output_taxonomies.csv` with 10 digit codes and genbank taxonomy.

**Usage**
Input a spreadsheet with ten digit codes in the first column and the genus and species names in the second column. Preferably, the genus and species name will be separated by a space and there will be no extraneous characters in the second column.

>`python get_taxonomy.py --input_file <path to .csv file>`


## For Assemblies dir:
### assess_transcriptomes.py:
Written March 2023 by Elinor (esterner27@gmail.com) to plot length, coverage and GC of assembled transcripts

**Purpose** This script will iterate through all assembled files (named as 10 digit code plus assembledTranscripts) with and gather GC, length and coverage. With that data, it plots R scripts


**Input**
	folder of renamed assemblies (if this is the case, put `-i` or --input in the command line)
	tsv file of ten digit code and taxon info (lifestage, species, etc) formatted like this: 10_digit_code\tdescriptor_of_taxon
	R script plot_assemblies.R, which is called from within this python script

**Usage**

To run:
>`python assess_transcriptomes.py --input  <pathway to directory of renamed assemblies>`

**Output** csv file of length, GC, coverage of each transcript, and multiple R plots, faceted by taxon and a csv file of data. It plots GC by length, and distributions of coverage, length and GC content across the whole transcript

## For Trees dir:
### count_tips.py

Written Auden Cote-L'Heureux, updated by Elinor (esterner27@gmail.com) to count occurrences of taxa in each OG in a batch of post guidance files.

**Purpose** This script will iterate through a directory of postguidance files and tally the number of times each taxa is found in an alignment.

**Input**
	folder of post guidance files (if this is the case, put `-i` or --input in the command line)

**Usage**

To run:
>`python count_tips.py --input  <pathway to directory of post guidance files>`

**Output** csv file of each taxa by each OG file, and tallies of the number of times each taxa is found in each OG



### Katz lab
>[About Katz Lab](https://www.science.smith.edu/katz-lab/)  &nbsp; \| &nbsp;
📧[Mail](lkatz@smith.edu) &nbsp; \| &nbsp; 📞 Call : (413) 585-3825 &nbsp;   \|   
:office: Address: Burton Hall 201, 46 College Lane,
Smith College, Northampton Massachusetts.
