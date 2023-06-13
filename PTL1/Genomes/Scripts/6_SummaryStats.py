import os, sys
import argparse
from Bio import SeqIO
import CUB
from statistics import mean
from math import ceil, floor
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np


def get_args():

	parser = argparse.ArgumentParser(
		prog = 'PTL6p1 Script 8: Stat Summary',
		description = "Updated March 31th, 2023 by Auden Cote-L'Heureux"
	)

	parser.add_argument('-i', '--input', type = str, required = True, help = 'Input path to the "Output" folder produced by PhyloToL Part 1. This folder should contain both the "ReadyToGO" and "Intermediate" folders.')
	parser.add_argument('-d', '--databases', type = str, default = '../Databases', help = 'Path to databases folder')
	parser.add_argument('-r', '--r2g_jf', action = 'store_true', help = 'Create ReadyToGo files filtered to only include sequences between the 25th and 75th percentile of silent-site GC content. Please be aware that these are not necessarily the correct or non-contaminant sequences; examine the GC3xENc plots carefully before using these data.')

	#Curate genetic code

	return parser.parse_args()


def hook_lens(args):

	print('\nGetting average OG lengths in the Hook DB...')

	len_by_og = { }
	for file in os.listdir(args.databases + '/db_OG'):
		if file.endswith('.fasta') and os.path.isfile(args.databases + '/db_OG/' + file.replace('.fasta', '.dmnd')):
			for rec in tqdm(SeqIO.parse(args.databases + '/db_OG/' + file, 'fasta')):
				if rec.id[-10:] not in len_by_og:
					len_by_og.update({ rec.id[-10:] : [] })

				len_by_og[rec.id[-10:]].append(len(str(rec.seq)))

	for og in len_by_og:
		len_by_og[og] = mean(len_by_og[og])

	return len_by_og


def aa_comp_lengths(args, gcodes):

	print('\nGetting amino acid composition data from ReadyToGo files...')

	r2g_lengths = { }; aa_comp = { }; recid_by_contig_n = { }
	for file in tqdm([f for f in os.listdir(args.input + '/ReadyToGo/ReadyToGo_AA')]):
		if file.endswith('.fasta') and file[:10] in gcodes:
			for rec in SeqIO.parse(args.input + '/ReadyToGo/ReadyToGo_AA/' + file, 'fasta'):
				r2g_lengths.update({ rec.id : len(str(rec.seq)) * 3 })

				fymink = 0; garp = 0; other = 0; total = 0; x = 0
				for char in str(rec.seq):
					if char in 'FYMINKfymink':
						fymink += 1
					elif char in 'GARPgarp':
						garp += 1
					elif char in 'Xx':
						x += 1
					else:
						other += 1

					total += 1

				aa_comp.update({ rec.id : { 'FYMINK' : fymink/total, 'GARP' : garp/total, 'Other' : other/total, 'X' : x/total } })

				recid_by_contig_n.update({ rec.id.split('Contig_')[-1].split('_')[0] : rec.id })

	print('\nGetting transcript sequence data from original assembled transcript files...')

	transcripts = { }; transcript_id_corr = { }
	for tax in tqdm([f for f in os.listdir(args.input + '/Intermediate/')]):
		if os.path.isdir(args.input + '/Intermediate/' + tax + '/Original'):
			for file in os.listdir(args.input + '/Intermediate/' + tax + '/Original'):
				if file.endswith('_GenBankCDS.fasta'):
					for rec in SeqIO.parse(args.input + '/Intermediate/' + tax + '/Original/' + file, 'fasta'):
						transcripts.update({ rec.id : (file[:10], str(rec.seq)) })
						if rec.id.split('NODE_')[-1].split('_')[0] in recid_by_contig_n:
							transcript_id_corr.update({ recid_by_contig_n[rec.id.split('NODE_')[-1].split('_')[0]] : rec.id})

	return aa_comp, transcripts, r2g_lengths, transcript_id_corr


def get_nuc_comp(args, gcodes):

	print('\nGetting nucleotide composition data from ReadyToGo files...')

	nuc_comp = { }
	for file in tqdm([f for f in os.listdir(args.input + '/ReadyToGo/ReadyToGo_NTD')]):
		if file.endswith('.fasta') and file[:10] in gcodes:
			cub_out = CUB.CalcRefFasta(args.input + '/ReadyToGo/ReadyToGo_NTD/' + file, gcodes[file[:10]].lower())[0]
			for k in cub_out:
				nuc_comp.update({ k : cub_out[k] })

	return nuc_comp


def per_seq(args, nuc_comp, aa_comp, all_transcripts, r2g_lengths, transcript_id_corr, og_mean_lens):

	if not os.path.isdir(args.input + '/PerSequenceStatSummaries'):
		os.mkdir(args.input + '/PerSequenceStatSummaries')

	taxa = list(dict.fromkeys([seq[:10] for seq in nuc_comp]))

	for taxon in taxa:
		with open(args.input + '/PerSequenceStatSummaries/' + taxon + '.csv', 'w') as o:
			o.write('Sequence,Taxon,OG,OrigName,OrigLength,R2GLength,AvgLengthOGinHook,AmbiguousCodons,GC-Overall,GC1,GC2,GC3,GC3-Degen,ExpWrightENc,ObsWrightENc_6Fold,ObsWrightENc_No6Fold,ObsWeightedENc_6Fold,ObsWeightedENc_No6Fold,FYMINK,GARP,OtherAA,N.Xs\n')
			for rec in nuc_comp:
				if rec[:10] == taxon:
					o.write(rec + ',' + rec[:10] + ',' + rec[-10:])

					try:
						o.write(',' + transcript_id_corr[rec] + ',' + str(len(all_transcripts[transcript_id_corr[rec]][1])))
					except KeyError:
						o.write(',NA,NA')

					o.write(',' + str(r2g_lengths[rec]) + ',' + str(round(og_mean_lens[rec[-10:]], 2)))

					v = nuc_comp[rec]
					gcs = [str(round(v.gcOverall, 2)), str(round(v.gc1, 2)), str(round(v.gc2, 2)), str(round(v.gc3, 2)), str(round(v.gc4F, 2))]
					ENc = [str(round(v.expENc, 2)), str(round(v.obsENc_6F, 2)), str(round(v.obsENc_No6F, 2)), str(round(v.SunENc_6F, 2)),str(round(v.SunENc_No6F, 2))]
					o.write(',' + ','.join([str(v.amb_cdn)] + gcs + ENc))

					o.write(',' + str(round(aa_comp[rec]['FYMINK'], 2)) + ',' + str(round(aa_comp[rec]['GARP'], 2)) + ',' + str(round(aa_comp[rec]['Other'], 2)) + ',' + str(round(aa_comp[rec]['X'], 2)) + '\n')


def per_tax(args, nuc_comp, aa_comp, all_transcripts, r2g_lengths, gcodes, og_mean_lens):

	taxa = list(dict.fromkeys([seq[:10] for seq in nuc_comp]))

	with open(args.input + '/PerTaxonSummary.csv', 'w') as o:
		o.write('Taxon,OrigSeqs,Orig_MedianGC,Orig_GCWidth_5-95Perc,Orig_MedianLen,Orig_IQRLen,R2GSeqs,R2GOGs,R2GMedian_GC3,R2G_5Perc_GC3,R2G_95Perc_GC3,R2G_GC3Width_5-95Perc,R2G_MedianENc,R2G_IQRENc,R2G_MedianLen,R2G_IQRLen,R2G_Prop.G1.5_OGAvg,R2G_Prop.L0.5_OGAvg,R2G_MeanXs,GeneticCode\n')

		for taxon in taxa:
			try:
				o.write(taxon)

				transcripts = [all_transcripts[seq][1].upper() for seq in all_transcripts if all_transcripts[seq][0] == taxon]
				o.write(',' + str(len(transcripts)))

				transcript_gcs = []
				for transcript in transcripts:
					transcript_gcs.append((transcript.count('G') + transcript.count('C'))/len(transcript))

				transcript_gcs = sorted(transcript_gcs)
				o.write(',' + str(round(transcript_gcs[floor(len(transcripts)*0.5)], 2)))
				o.write(',' + str(round(transcript_gcs[floor(len(transcripts)*0.95)] - transcript_gcs[floor(len(transcripts)*0.05)], 2)))

				transcript_lens = sorted([len(transcript) for transcript in transcripts])
				o.write(',' + str(round(transcript_lens[floor(len(transcripts)*0.5)], 2)))
				o.write(',' + str(round(transcript_lens[floor(len(transcripts)*0.75)] - transcript_lens[floor(len(transcripts)*0.25)], 2)))

				r2g_ntds = [nuc_comp[seq] for seq in nuc_comp if seq[:10] == taxon]
				o.write(',' + str(len(r2g_ntds)))
				r2g_ogs = list(dict.fromkeys([seq[-10:] for seq in nuc_comp if seq[:10] == taxon]))
				o.write(',' + str(len(r2g_ogs)))

				r2g_gc3s = sorted([seq.gc4F for seq in r2g_ntds])
				o.write(',' + str(round(r2g_gc3s[floor(len(r2g_ntds)*0.5)], 2)))
				o.write(',' + str(round(r2g_gc3s[floor(len(r2g_gc3s)*0.05)], 2)))
				o.write(',' + str(round(r2g_gc3s[floor(len(r2g_gc3s)*0.95)], 2)))
				o.write(',' + str(round(r2g_gc3s[floor(len(r2g_gc3s)*0.95)] - r2g_gc3s[floor(len(r2g_gc3s)*0.05)], 2)))

				r2g_encs = sorted([seq.obsENc_6F for seq in r2g_ntds])
				o.write(',' + str(round(r2g_encs[floor(len(r2g_encs)*0.5)], 2)))
				o.write(',' + str(round(r2g_encs[floor(len(r2g_encs)*0.75)] - r2g_encs[floor(len(r2g_encs)*0.25)], 2)))

				tax_r2g_lens = sorted([r2g_lengths[seq] for seq in r2g_lengths if seq[:10] == taxon])
				o.write(',' + str(round(tax_r2g_lens[floor(len(tax_r2g_lens)*0.5)], 2)))
				o.write(',' + str(round(tax_r2g_lens[floor(len(tax_r2g_lens)*0.75)] - tax_r2g_lens[floor(len(tax_r2g_lens)*0.25)], 2)))

				prop_len_g = len([seq for seq in r2g_lengths if seq[:10] == taxon and r2g_lengths[seq] > 4.5 * og_mean_lens[seq[-10:]]])/len(tax_r2g_lens)
				prop_len_l = len([seq for seq in r2g_lengths if seq[:10] == taxon and r2g_lengths[seq] < 1.5 * og_mean_lens[seq[-10:]]])/len(tax_r2g_lens)

				o.write(',' + str(round(prop_len_g, 2)) + ',' + str(round(prop_len_l, 2)))

				o.write(',' + str(mean([aa_comp[seq]['X'] for seq in aa_comp if seq[:10] == taxon])))

				o.write(',' + gcodes[taxon] + '\n')
			except:
				pass


def r2g_jf(args, nuc_comp, gcodes):

	#Q: should there be an maximum IQR cutoff at which we do NOT produce a file here?

	if not os.path.isdir(args.input + '/ReadyToGo/ReadyToGo_NTD_JF'):
		os.mkdir(args.input + '/ReadyToGo/ReadyToGo_NTD_JF')

	if not os.path.isdir(args.input + '/ReadyToGo/ReadyToGo_AA_JF'):
		os.mkdir(args.input + '/ReadyToGo/ReadyToGo_AA_JF')

	for file in os.listdir(args.input + '/ReadyToGo/ReadyToGo_NTD'):
		if file.endswith('.fasta') and file[:10] in gcodes:
			taxon = file[:10]

			r2g_ntds = [nuc_comp[seq] for seq in nuc_comp if seq[:10] == taxon]
			r2g_gc3s = sorted([seq.gc4F for seq in r2g_ntds])

			with open(args.input + '/ReadyToGo/ReadyToGo_NTD_JF/' + file.replace('.fasta', '.JF.fasta'), 'w') as o:
				for rec in SeqIO.parse(args.input + '/ReadyToGo/ReadyToGo_NTD/' + file, 'fasta'):
					if nuc_comp[rec.id].gc4F > r2g_gc3s[floor(len(r2g_gc3s)*0.25)] and nuc_comp[rec.id].gc4F < r2g_gc3s[floor(len(r2g_gc3s)*0.75)]:
						o.write('>' + rec.id + '\n' + str(rec.seq) + '\n\n')

			with open(args.input + '/ReadyToGo/ReadyToGo_AA_JF/' + file.replace('.fasta', '.JF.fasta').replace('NTD', 'AA'), 'w') as o:
				for rec in SeqIO.parse(args.input + '/ReadyToGo/ReadyToGo_AA/' + file.replace('NTD', 'AA'), 'fasta'):
					if nuc_comp[rec.id].gc4F > r2g_gc3s[floor(len(r2g_gc3s)*0.25)] and nuc_comp[rec.id].gc4F < r2g_gc3s[floor(len(r2g_gc3s)*0.75)]:
						o.write('>' + rec.id + '\n' + str(rec.seq) + '\n\n')


def plot_jf(args, nuc_comp):

	if not os.path.isdir(args.input + '/GC3xENc_Plots'):
		os.mkdir(args.input + '/GC3xENc_Plots')

	taxa = list(dict.fromkeys([rec[:10] for rec in nuc_comp]))

	gc3_null = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100]
	enc_null = [31, 31.5958, 32.2032, 32.8221, 33.4525, 34.0942, 34.7471, 35.411, 36.0856, 36.7707, 37.4659, 38.1707, 38.8847, 39.6074, 40.3381, 41.0762, 41.8208, 42.5712, 43.3264, 44.0854, 44.8471, 45.6102, 46.3735, 47.1355, 47.8949, 48.65, 49.3991, 50.1406, 50.8725, 51.593, 52.3, 52.9916, 53.6656, 54.32, 54.9525, 55.561, 56.1434, 56.6975, 57.2211, 57.7124, 58.1692, 58.5898, 58.9723, 59.3151, 59.6167, 59.8757, 60.0912, 60.2619, 60.3873, 60.4668, 60.5, 60.4668, 60.3873, 60.2619, 60.0912, 59.8757, 59.6167, 59.3151, 58.9723, 58.5898, 58.1692, 57.7124, 57.2211, 56.6975, 56.1434, 55.561, 54.9525, 54.32, 53.6656, 52.9916, 52.3, 51.593, 50.8725, 50.1406, 49.3991, 48.65, 47.8949, 47.1355, 46.3735, 45.6102, 44.8471, 44.0854, 43.3264, 42.5712, 41.8208, 41.0762, 40.3381, 39.6074, 38.8847, 38.1707, 37.4659, 36.7707, 36.0856, 35.411, 34.7471, 34.0942, 33.4525, 32.8221, 32.2032, 31.5958, 31]

	for taxon in taxa:
		comp_data = [(nuc_comp[rec].gc4F, nuc_comp[rec].obsENc_6F) for rec in nuc_comp if rec[:10] == taxon]

		plt.figure()
		plt.plot(np.array(gc3_null), np.array(enc_null), color = 'black', linewidth=2)
		plt.scatter(np.array([val[0] for val in comp_data]), np.array([val[1] for val in comp_data]), s = 1)
		plt.xlabel("GC content (3rd pos, 4-fold sites)")
		plt.ylabel("Observed Wright ENc (6 Fold)")
		plt.savefig(args.input + '/GC3xENc_Plots/' + taxon + '.png')
	
if __name__ == "__main__":
	args = get_args()

	valid_codes = ['universal', 'blepharisma', 'chilodonella', 'condylostoma', 'euplotes', 'peritrich', 'vorticella', 'mesodinium', 'tag', 'tga', 'taa', 'none']

	gcodes = { }
	if os.path.isfile(args.input + '/Intermediate/gcode_output.tsv'):
		for line in open(args.input + '/Intermediate/gcode_output.tsv'):
			if len(line.split('\t')) == 5 and line.split('\t')[4].strip().lower() in valid_codes:
				gcodes.update({ line.split('\t')[0] : line.split('\t')[4].strip() })
			elif line.split('\t')[4].strip().lower() != '':
				print('\nInvalid genetic code assignment for taxon ' + line.split('\t')[0] + '. Skipping this taxon in script 6 (summary statistics)\n')
	else:
		print('\nGenetic code assignment file (Output/Intermediate/gcode_output.tsv) not found. Quitting script 6 (summary statistics).\n')
		exit()

	aa_comp, transcripts, r2g_lengths, transcript_id_corr = aa_comp_lengths(args, gcodes)
	nuc_comp = get_nuc_comp(args, gcodes)
	og_mean_lens = hook_lens(args)

	per_tax(args, nuc_comp, aa_comp, transcripts, r2g_lengths, gcodes, og_mean_lens)
	per_seq(args, nuc_comp, aa_comp, transcripts, r2g_lengths, transcript_id_corr, og_mean_lens)

	if args.r2g_jf:
		r2g_jf(args, nuc_comp, gcodes)

	plot_jf(args, nuc_comp)





















