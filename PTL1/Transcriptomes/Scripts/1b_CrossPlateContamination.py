#!/usr/bin/python3

__author__ = 'Jean-David Grattepanche'
__version__ = 'ACL fixed sequence naming issue Feb 23, 2022'
__email__ = 'jeandavid.grattepanche@gmail.com'



import sys
import os
import re
import time
import string
import os.path
from Bio import SeqIO
from sys import argv
listtaxa=[]
toosim = 0.99
seqcoverage = 0.7

def merge_files(folder, minlen, conspecific_names):
	mergefile = open('/'.join(folder.split('/')[:-1]) + '/forclustering.fasta','w+')
	print("MERGE following files")
	for taxafile in os.listdir(folder):
		if taxafile[0] != ".":
			listtaxa.append(taxafile.split('.' + str(minlen) + 'bp')[0])

			for line2 in SeqIO.parse(folder+'/'+taxafile, 'fasta'):
				if int(len(str(line2.seq))) >= int(minlen):
					mergefile.write('>'+taxafile.split('.' + str(minlen) + 'bp')[0] + '_' + line2.description + '\n' + str(line2.seq) + '\n')
				else:
					print(line2, " is too short")
	mergefile.close()

	sort_cluster(folder, listtaxa, minlen, conspecific_names)


def sort_cluster(folder, listtaxa, minlen, conspecific_names):
	if not os.path.exists('/'.join(folder.split('/')[:-1]) + '/clusteringresults_vsearch/'):
		os.makedirs('/'.join(folder.split('/')[:-1]) + '/clusteringresults_vsearch/')

	fastalist = []; fastadict= {}
	conspecific_names_dict = { line.split('\t')[0] : line.split('\t')[1].strip() for line in open(conspecific_names) }

	print('CREATE a dictionnary of sequences')
	for record in SeqIO.parse(open('/'.join(folder.split('/')[:-1]) + '/forclustering.fasta','r'),'fasta'):
		if record.id[:10] not in conspecific_names_dict:
			print('\nError in cross-plate contamination assessment: the ten-digit code ' + record.id[:10] + ' is not found in the conspecific names file. Please check that this file is correct and try again.\n')
			quit()

		IDL  = record.description, int(record.description.split('_Cov')[1].replace('\n',''))
		fastalist.append(IDL)
		fastadict[record.description] = record.seq

	print("CLUSTER sequences that overlap at least 70%")
	os.system('vsearch --cluster_fast ' + '/'.join(folder.split('/')[:-1]) + '/forclustering.fasta --strand both --query_cov '+str(seqcoverage)+' --id '+str(toosim) +' --uc ' + '/'.join(folder.split('/')[:-1]) + '/clusteringresults_vsearch/results_forclustering.uc --threads 60' )
	
	#input2 = open('/'.join(folder.split('/')[:-1]) + '/clusteringresults_vsearch/results_forclustering.uc','r')
	#input2 = open('/Output_PostClusterBackup/clusteringresults_vsearch/results_forclustering.uc','r')
	cluster_output = '/'.join(folder.split('/')[:-1]) + '/clusteringresults_vsearch/results_forclustering.uc'
	out2 = open('/'.join(folder.split('/')[:-1]) + '/fastatokeep.fas','w+')
	out3 = open('/'.join(folder.split('/')[:-1]) + '/fastatoremoved.fas','w+')
	out4 = open('/'.join(folder.split('/')[:-1]) + '/fastatoremoved.uc','w+')
	print("CREATE a dictionary with clustering results")
	clustdict= {}; clustlist = []; allseq = []; clustline = {}; list= []; i=0; j=0
	for row2 in open(cluster_output, 'r'):
		if row2.split('\t')[0] == 'C' and int(row2.split('\t')[2]) < 2: # keep all unique sequences
			out2.write('>'+row2.split('\t')[8] + '\n' + str(fastadict[row2.split('\t')[8]])+ '\n')
		if row2.split('\t')[0] == 'C' and int(row2.split('\t')[2]) > 1: # create another dictionary
#			print("create dico: ", row2.split('\t')[8])
			clustdict.setdefault(row2.split('\t')[8], [row2.split('\t')[8]])
			clustlist.append(row2.split('\t')[8])

	for row3 in open(cluster_output, 'r'):
		if row3.split('\t')[0] == 'H':
#			print("add dico: ", row3.split('\t')[9], row3.split('\t')[8])
			clustdict[row3.split('\t')[9].replace('\n','')].append(row3.split('\t')[8].replace('\n',''))
			clustline[row3.split('\t')[8].replace('\n','')] = row3.replace('\n','')
			clustline[row3.split('\t')[9].replace('\n','')] = row3.replace('\n','')

	print("PARSE the clusters: keep seed sequences (highest coverage) for each cluster")
	for clust in clustlist:
		list = sorted(clustdict[clust], reverse = True, key=lambda x: int(x.split('_Cov')[1]))
		master = list[0]
		Covmaster = int(list[0].split('_Cov')[1])
		master8dig = ('_').join(list[0].split('_')[0:3])[:-2]
		for seq in list:
			clustered =  seq.replace('\n','')
			Covclustered = int(clustered.split('_Cov')[1])
			clustered8dig = ('_').join(clustered.split('_')[0:3])[:-2]
# 			print(master8dig, Covmaster, '//', clustered8dig, Covclustered)
			if float(Covmaster/Covclustered) < 10:
				out2.write('>'+clustered + '\n' + str(fastadict[clustered])+ '\n')
				i +=1
			elif conspecific_names_dict[master[:10]] == conspecific_names_dict[clustered[:10]]:
				out2.write('>'+clustered + '\n' + str(fastadict[clustered])+ '\n')
				i +=1
			elif Covclustered >= 50:				
				out2.write('>'+clustered + '\n' + str(fastadict[clustered])+ '\n')
				i +=1
			else:
				j +=1
				out4 = open('/'.join(folder.split('/')[:-1]) + '/fastatoremoved.uc','a')
				out3.write('>'+clustered + '\n' + str(fastadict[clustered])+ '\n')
				print(clustline[clustered],'\t' , master )
				out4.write(clustline[clustered]+ '\t' + master + '\n')
				out4.close()


	print('there are ', str(i),' sequences kept and ',str(j),' sequences removed')
	
	out2.close()
	out3.close()

	splittaxa(folder, listtaxa, minlen)
	
def splittaxa(folder, listtaxa, minlen):
	for taxa in listtaxa:
		tax_sf_path = '/'.join(folder.split('/')[:-1]) + '/' + taxa + '/SizeFiltered/'
		os.system('mv ' + tax_sf_path + taxa + '.' + str(minlen) + 'bp.fasta' + ' ' + tax_sf_path + taxa + '.' + str(minlen) + 'bp.preXPlate.fasta')

		with open(tax_sf_path + taxa + '.' + str(minlen) + 'bp.fasta','w') as o:
			for kept in SeqIO.parse('/'.join(folder.split('/')[:-1]) + '/fastatokeep.fas','fasta'):
				if taxa in kept.description:
					o.write('>' + kept.description.replace(taxa + '_', '') + '\n' + str(kept.seq) + '\n')

	os.system('mv ' + '/'.join(folder.split('/')[:-1]) + '/fastatokeep.fas ' + '/'.join(folder.split('/')[:-1]) + '/clusteringresults_vsearch/')
	os.system('mv ' + '/'.join(folder.split('/')[:-1]) + '/fastatoremoved.fas ' + '/'.join(folder.split('/')[:-1]) + '/clusteringresults_vsearch/')
	os.system('mv ' + '/'.join(folder.split('/')[:-1]) + '/fastatoremoved.uc ' + '/'.join(folder.split('/')[:-1]) + '/clusteringresults_vsearch/')
	os.system('mv ' + '/'.join(folder.split('/')[:-1]) + '/forclustering.fasta ' + '/'.join(folder.split('/')[:-1]) + '/clusteringresults_vsearch/')
	
def main():

	script, folder, minlen, conspecific_names = argv
	merge_files(folder, minlen, conspecific_names)

main()
















