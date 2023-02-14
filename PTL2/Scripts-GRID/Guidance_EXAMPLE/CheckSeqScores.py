import os,sys



with open('ScoreBySeq.tsv', 'w') as o:
	for clade in os.listdir('GuidanceOutput'):
		if(os.path.isfile('GuidanceOutput/' + clade + '/MSA.MAFFT.Guidance2_res_pair_seq.scr_with_Names')):
			for line in open('GuidanceOutput/' + clade + '/MSA.MAFFT.Guidance2_res_pair_seq.scr_with_Names'):
				if('SEQUENCE_NAME' not in line):
					o.write(line)




















