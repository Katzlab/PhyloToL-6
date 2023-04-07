![Lab image](https://github.com/Katzlab/PhyloToL-6/blob/621517812a5ed4256f19ba26498a665bdd5060af/Katzlab.png)
# Utilities
> This folder contains many useful tools for analyzing sequence data.

## For taxonomy dir:
### Query_SRA_egs.py: 
input list of taxonomic names and output all SRAs or GCA since 2020 (can be adjusted by modifying script). For SRAs, the script also gives sequsncing technology used (pacbio, miseq, etc) and experiment type. It excludes all SRAs that include the word 'amplicon'.  Input is folder 'unique_taxon_lists' with files of keywords by major clade (separated by new lines). Put -t (transcriptome, SRA db) or -g (genome, assembly db) in the command line to specify data type. 
>  Example command line: `python Query_SRA_egs.py -t OR -g`


### Katz lab
>[About Katz Lab](https://www.science.smith.edu/katz-lab/)  &nbsp; \| &nbsp;
📧[Mail](lkatz@smith.edu) &nbsp; \| &nbsp; 📞 Call : (413) 585-3825 &nbsp;   \|   
:office: Address: Burton Hall 201, 46 College Lane,
Smith College, Northampton Massachusetts.
