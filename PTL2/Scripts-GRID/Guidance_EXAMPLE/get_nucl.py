import os, sys
from Bio import SeqIO


nucls = { rec.id : str(rec.seq) for rec in SeqIO.parse('AllCuratedNTDSeqs_020722_previous.fasta', 'fasta') }

with open('AllCuratedNTDSeqs_022222_Seq0.7.fasta', 'w') as o:
	for file in os.listdir('SeqFiltered_ToAlign/Preguidance_0.7_Aligned'):
		if('.fas' in file):
			for rec in SeqIO.parse('SeqFiltered_ToAlign/Preguidance_0.7_Aligned/' + file, 'fasta'):
				o.write('>' + rec.id + '\n' + nucls[rec.id] + '\n\n')
