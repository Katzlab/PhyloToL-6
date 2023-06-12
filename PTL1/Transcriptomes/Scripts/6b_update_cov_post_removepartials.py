#!/usr/bin/python
from __future__ import print_function

__author__ = "Jean-David Grattepanche"
__version__ = "2, August 28, 2017"
__email__ = "jeandavid.grattepanche@gmail.com"


import sys
import os
import re
import time
import string
import os.path
from Bio import SeqIO
from sys import argv

def Addcoverage(code):
	seqfolder = code
	all_output_folder = '/'.join(code.split('/')[:-1])
	code = code.split('/')[-1]
	
	covupd = {}
	for seqcoll in open(seqfolder + '/' + code + '_SeqPairsAbove98.txt','r'):
		CL = 0
		for transc in seqcoll.split('\t'):
			if CL == 0:
				reftrans = ('_').join(transc.split('_')[1:])
			coverage = int(transc.split('Cov')[1].split('_')[0])
			Length = int(transc.split('Len')[1].split('_')[0])
			CL += coverage * Length
		covupd[reftrans] = CL 

	if os.path.isdir(seqfolder + '/Updated_Coverage/') != True:
		os.system('mkdir ' + seqfolder + '/Updated_Coverage/')
	if os.path.isdir(seqfolder + '/Updated_Coverage/SpreadSheets/') != True:
		os.system('mkdir ' + seqfolder + '/Updated_Coverage/SpreadSheets/')

	for spreadsh in os.listdir(seqfolder + '/Processed/SpreadSheets/'):
		if spreadsh.endswith('.tsv'):
			outtsvtokeep = open(seqfolder + '/Updated_Coverage/SpreadSheets/' + spreadsh.split('Final')[0] + 'UC.Final' + spreadsh.split('Final')[1],'w+')
			for row in open(seqfolder + '/Processed/SpreadSheets/'+ spreadsh, 'r'):
				if row.split('_Trans')[0] in covupd:
					og_number = re.split('OG.{1}_', row)[-1][:6]
					og_prefix = row.split(og_number)[0][-4:]
					og = og_prefix + og_number

					newcov2 = round(covupd[row.split('_Trans')[0]] / int(row.split('_Len')[1].split('_')[0]))
					outtsvtokeep.write(row.split('Cov')[0]+'Cov'+str(newcov2)+'_' + og_prefix +row.split(og_prefix)[1].split('_Trans')[0] +'\t' +('\t').join(row.split('\t')[1:]))
				else: 
					if 'Trans' in row:
						outtsvtokeep.write(row.split('_Trans')[0]+ '\t' +('\t').join(row.split('\t')[1:]))
					else:
						outtsvtokeep.write(row)
			outtsvtokeep.close()

	for seqfile in os.listdir(seqfolder + '/Processed'):
		if seqfile.endswith('.fasta'):
			outseqtokeep = open(seqfolder + '/Updated_Coverage/' + seqfile.split('Final')[0] + 'UC.Final' + seqfile.split('Final')[1],'w+')
			for Seq in SeqIO.parse(seqfolder + '/Processed/' + seqfile ,'fasta'):
				if Seq.description.split('_Trans')[0] not in covupd:
					outseqtokeep.write('>'+Seq.description.split('_Trans')[0]+ '\n'+str(Seq.seq) +'\n')
				else:
					og_number = re.split('OG.{1}_', Seq.description)[-1][:6]
					og_prefix = Seq.description.split(og_number)[0][-4:]
					og = og_prefix + og_number

					newcov = round(covupd[Seq.description.split('_Trans')[0]] / int(Seq.description.split('_Len')[1].split('_')[0]))
					outseqtokeep.write('>'+Seq.description.split('Cov')[0]+'Cov'+str(newcov)+'_' + Seq.description.split(og)[0][-2:] + og + '\n'+str(Seq.seq) +'\n')
			outseqtokeep.close()

	if os.path.isdir(all_output_folder + '/ToRename') != True:
		os.system('mkdir ' + all_output_folder + '/ToRename')

	os.system('cp ' + seqfolder + '/Updated_Coverage/*fasta ' + all_output_folder + '/ToRename/')
	os.system('cp ' + seqfolder + '/Updated_Coverage/SpreadSheets/*tsv ' + all_output_folder + '/ToRename/')
			
			
def main():
	script, code = argv
	Addcoverage(code)
main()






