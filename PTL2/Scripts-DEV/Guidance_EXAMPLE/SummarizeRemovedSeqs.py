import os, sys
from Bio import SeqIO
from tqdm import tqdm


prev_recs = [rec.id for rec in SeqIO.parse('AllCuratedNTDSeqs_020722_previous.fasta', 'fasta')]

recs5 = { rec.id : file[:12] for file in os.listdir('Seq0.5') for rec in SeqIO.parse('Seq0.5/' + file, 'fasta') }
recs7 = { rec.id : file[:12] for file in os.listdir('Seq0.5') for rec in SeqIO.parse('Seq0.7/' + file, 'fasta') }

#print(len(prev_recs))
print(recs5)
#print(len(recs7))

n_file = open('NumRemovedSeqs.csv', 'w')
seqs5_file = open('AllRemovedSeqs_0.5.csv', 'w')
seqs7_file = open('AllRemovedSeqs_0.7.csv', 'w')

n_by_clade_5 = { }; n_by_clade_7 = { }
for rec in prev_recs:
	if(recs5[rec] not in n_by_clade_5):
		n_by_clade_5.update({ recs5[rec] : 0 })
		n_by_clade_7.update({ recs5[rec] : 0 })
	if(rec not in recs5):
		seqs5_file.write(rec + '\n')
		n_by_clade_5[recs5[rec]] += 1
	if(rec not in recs7):
		seqs7_file.write(rec + '\n')
		n_by_clade_7[recs5[rec]] += 1


seqs5_file.close()


