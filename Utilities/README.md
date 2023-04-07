![Lab image](https://github.com/Katzlab/PhyloToL-6/blob/621517812a5ed4256f19ba26498a665bdd5060af/Katzlab.png)
# Utilities
> This folder contains many useful tools for analyzing sequence data.

## For taxonomy dir:
### Query_SRA_egs.py:

**Purpose** Gives a spreadsheet of all SRA codes or GCA codes from Genbank associated with taxonomic terms in an input csv file

**Input** Input is folder 'unique_taxon_lists' with files of keywords by major clade (separated by new lines) by get_unique_taxa.py or manually 

**Output** all SRAs or GCA since 2020 (can be adjusted by modifying script). For SRAs, the script also gives sequecing technology used (pacbio, miseq, etc) and experiment type. It excludes all SRAs that include the word 'amplicon'.

**Usage** -t (transcriptome, searches SRA db) or -g (genome, searches assembly db) in the command line to specify data type. 
>  Example command line: `python Query_SRA_egs.py -t OR -g`

### get_unique_taxa.py: 
Written by Elinor 1/26, updated 2/12

**Purpose** make lists of unique taxonomy from phylotol master taxonomy column (genbank taxonomy for each taxa in the pipeline). These lists are the intended input for Query_SRA_egs.py. This cuts off the genus (and species if there is one), uniquifies the list and writes them out to files by the first word of the taxonomy

**Input** text file of taxonomies. make sure each taxonomic level is separated with '; ' (semicolon space) or the script will not parse the names right


WARNING: if you run the script multiple times, DELETE THE PREVIOUS OUTPUT. this is because it appends lines to the 
end of files so you will have many duplicates

>  Example command line: `python get_unique_taxa.py`

### Katz lab
>[About Katz Lab](https://www.science.smith.edu/katz-lab/)  &nbsp; \| &nbsp;
📧[Mail](lkatz@smith.edu) &nbsp; \| &nbsp; 📞 Call : (413) 585-3825 &nbsp;   \|   
:office: Address: Burton Hall 201, 46 College Lane,
Smith College, Northampton Massachusetts.
