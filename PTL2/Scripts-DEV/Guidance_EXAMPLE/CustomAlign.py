import os, sys
from Bio import SeqIO


score_by_seq = { line.split('\t')[0] : float(line.split('\t')[1]) for line in open('ScoreBySeq.tsv') }

remove3 = [seq for seq in score_by_seq if score_by_seq[seq] < 0.3]
remove5 = [seq for seq in score_by_seq if score_by_seq[seq] < 0.5]
remove7 = [seq for seq in score_by_seq if score_by_seq[seq] < 0.7]

for clade in os.listdir('GuidanceOutput'):
	if('OG5_' in clade):
		with open('SeqFiltered_ToAlign/Preguidance_0.7/' + clade[:12] + '.fasta', 'w') as o:
			for rec in SeqIO.parse('GuidanceOutput/' + clade + '/Seqs.Orig.fas', 'fasta'):
				if(rec.id not in remove7):
					o.write('>' + rec.id + '\n' + str(rec.seq) + '\n\n')