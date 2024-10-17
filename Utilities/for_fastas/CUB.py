#Author, date: Xyrus (last modified by him Sept 17 2020), most recently updated by Auden on October 17 2024
#Motivation: Generate lots of codon usage statistics to aid in identifying useful characteristics for de novo ORF calling
#Intent: Summarize nucleotide composition statistics for a fasta file or folder of fasta files
#Dependencies: Python3, numpy, BioPython
#Inputs: Fasta file or folder of fasta files
#Outputs: A fasta file filtered for properly formatted sequences and several spreadsheets summarizing GC, ENc, RSCU, etc.
#Example: python3 CUB.py -i seqs.fasta 
#Note: Use "python3 CUB.py -i seqs.fasta --require_start --require_stop" when using the on R2G files.

#Dependencies
import os
import re
import sys
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC
import argparse

class CalcCUB:
    """
    Returns the Effective Number of Codons used (observed and expected)
    following the equations originally from Wright 1990.
    """
    def expWrightENc(gc3):
        # Calculates the expected ENc from a sequence's GC3 under Wright 1990
        if gc3 > 1:
            # If GC3 looks as though it is > 1 (e.g. 100%), converts to a float ≤ 1.
            # Calculations expect a value between 0 and 1
            gc3 = gc3/100
        exp_enc = 2+gc3+(29/((gc3**2)+(1-gc3)**2))
        return round(exp_enc, 4)

    def nullENcGC3():
        # Calculates the expected ENc from the null distribution of GC3
        # values (0, 100% GC)
        null = [CalcCUB.expWrightENc(n) for n in np.arange(0,.51,0.01)]
        null += null[:-1][::-1]
        return [str(i)+'\t'+str(j) for i, j in zip([n for n in range(0, 101)],null)]


    def calcWrightENc(cdnTable):
        # Follows Wright's (1990) calculations for determining ENc scores.

        def faCalcWright(aa_counts):
            # Returns the codon homozygosity (fa) for a given "type" of AA (e.g.
            # 2-fold degeneracy).
            counts = [i[2] for i in aa_counts]
            # n_aa --> number of this particular AA
            n_aa = sum(counts)
            # fa --> codon homozygosity
            try:
                fa = (((n_aa*sum([(i/float(n_aa))**2 for i in counts]))-1)/(n_aa-1))
            except:
                fa = 0
            return fa

        def ENcWright_by_Degen(fa_data):
            # Same as used in Wright 1990, averages the homozygosity across all codons
            # of a given class (e.g. 2-fold degeneracy)

            # Codons without any degeneracy (e.g. ATG == M) have 100% homozygosity
            # and provide a "base" for the ENc score
            enc = 2
            for k, v in fa_data.items():
                non_zero_vals, non_zero_sum = len([i for i in v if i != 0]), sum([i for i in v if i != 0])
                try:
                    f_aa = non_zero_sum/non_zero_vals
                except:
                    f_aa = 1
                enc += k/f_aa
            return enc

        # Determines the number of degenerate groups to use (i.e. whether 6-Fold
        # degeneracy is present).
        degen_cdns = {}
        for k, v in cdnTable.items():
            if v[1] not in degen_cdns.keys():
                degen_cdns[v[1]] = [v[0]]
            else:
                if v[0] not in degen_cdns[v[1]]:
                    degen_cdns[v[1]] += [v[0]]

        # Calculates codon homozygosity (fa) for each amino acid. Groups the
        # resulting values based on the amino acids degeneracy (e.g. 'two-fold').
        fa_cdns = {len(v):[] for k, v in degen_cdns.items() if 'one' not in k}

        for k, v in degen_cdns.items():
            # Skip codons lacking degeneracy
            if 'one' in k:
                continue

            for aa in v:
                aa_counts = [cdnTable[k] for k in cdnTable.keys() if cdnTable[k][0] == aa]
                fa_cdns[len(v)] += [faCalcWright(aa_counts)]
        enc_val = min(61, round(ENcWright_by_Degen(fa_cdns),4))
        return enc_val

    def SunEq5(cdnTable):
        def calcFcf(aa_counts):
            counts = [i[2] for i in aa_counts]
            pseudocounts = [i+1 for i in counts]
            na = sum(pseudocounts)
            fcf = sum([(i/float(na))**2 for i in pseudocounts]), sum(pseudocounts)
            return fcf

        ENcWeightedPsuedo = 0
        degen_cdns = {}

        for k, v in cdnTable.items():
            if v[1] == 'none':
                continue
            if v[1] not in degen_cdns.keys():
                degen_cdns[v[1]] = [v[0]]
            else:
                if v[0] not in degen_cdns[v[1]]:
                    degen_cdns[v[1]] += [v[0]]
        for k, v in degen_cdns.items():
            fcf_nc = []
            for aa in v:
                aa_counts = [cdnTable[k] for k in cdnTable.keys() if cdnTable[k][0] == aa]
                fcf_nc.append(calcFcf(aa_counts))
            weightedENc = (len(fcf_nc) /
                (sum([i[0]*i[1] for i in fcf_nc]) /
                sum([i[1] for i in fcf_nc])))
            ENcWeightedPsuedo += weightedENc
        return round(ENcWeightedPsuedo,4)

    def calcRCSU(cdnTbl):
        rscu = {k:[v[0]] for k, v in cdnTbl.items() if v[0].isalpha()}
        for k, v in rscu.items():
            try:
                aa_info = [(key, val[-1]) for key, val in cdnTbl.items() if val[0] == v[0]]
                aa_cnts = [x[1] for x in aa_info]
                cdn_rscu = (cdnTbl[k][-1]*len(aa_cnts))/sum(aa_cnts)
                rscu[k] += [str(round(cdn_rscu,4))]
            except:
                rscu[k] += ['NA']
        return rscu


class GenUtil(object):
    """
    "Overflow" of functions for now. Just a precaution to make the code a
    little cleaner/easier to manage.

    This class inclues means to normalize/check the user-provided genetic code,
    which if not valid will default to the "universal" genetic code.

    Similarly, This class will return the appropriate
    codon count table and provides a function to update its values.
    """
    def convertGenCode(gCode):
        # Will interpret the user provided genetic code (gcode) and checks that
        # it is currently available for use with the NCBI/biopython
        # supported translation tables. Default is universal.
        # Dictionary of the possible/functional genetic codes that are supported.
        # --- Chilodonella and condylostoma are to come!
        transTable = {'universal':1, 'blepharisma':4,
            'ciliate':6, 'euplotes':10, 'mesodinium':29, 'myrionecta':29, 'peritrich':30,
            '1':1, '4':4, '6':6, '10':10, '29':29, '30':30, 'chilo':'chilo'}

        if str(gCode).lower() not in transTable:
            print("\nWarning: Provided genetic code is not supported (yet).\n")
            print("Currently running using the UNIVERSAL genetic code.\n\n")
            print("Alternative genetic codes are as follows (Note: numbers "\
                "correspond to NCBI genetic code tables):\n")
            print('\n'.join(list(transTable.keys()))+'\n')
            return 'Universal',1
        else:
            return gCode,transTable[str(gCode).lower()]

    def getCDNtable(gCode):
        # Returns the appropriate codon table to be used for the ENc calculations.
        # Universal codon table, with 6-fold degenerate codons split
        # into four-fold and two-fold groups.
        universal_no6fold = {
        	'GCT': ['A', 'four', 0], 'GCC': ['A', 'four', 0], 'GCA': ['A', 'four', 0],
        	'GCG': ['A', 'four', 0], 'CGT': ['R', 'four', 0], 'CGC': ['R', 'four', 0],
        	'CGG': ['R', 'four', 0], 'CGA': ['R', 'four', 0], 'AGA': ['R_', 'two', 0],
        	'AGG': ['R_', 'two', 0], 'AAT': ['N', 'two', 0], 'AAC': ['N', 'two', 0],
        	'GAT': ['D', 'two', 0], 'GAC': ['D', 'two', 0], 'TGT': ['C', 'two', 0],
        	'TGC': ['C', 'two', 0], 'CAA': ['Q', 'two', 0], 'CAG': ['Q', 'two', 0],
        	'GAA': ['E', 'two', 0], 'GAG': ['E', 'two', 0], 'GGT': ['G', 'four', 0],
        	'GGC': ['G', 'four', 0], 'GGA': ['G', 'four', 0], 'GGG': ['G', 'four', 0],
        	'CAT': ['H', 'two', 0], 'CAC': ['H', 'two', 0], 'ATT': ['I', 'three', 0],
        	'ATC': ['I', 'three', 0], 'ATA': ['I', 'three', 0], 'ATG': ['M', 'one', 0],
        	'TTA': ['L_', 'two', 0], 'TTG': ['L_', 'two', 0], 'CTT': ['L', 'four', 0],
        	'CTC': ['L', 'four', 0], 'CTA': ['L', 'four', 0], 'CTG': ['L', 'four', 0],
        	'AAA': ['K', 'two', 0], 'AAG': ['K', 'two', 0], 'TTT': ['F', 'two', 0],
        	'TTC': ['F', 'two', 0], 'CCT': ['P', 'four', 0], 'CCC': ['P', 'four', 0],
        	'CCA': ['P', 'four', 0], 'CCG': ['P', 'four', 0], 'TCT': ['S', 'four', 0],
        	'TCC': ['S', 'four', 0], 'TCA': ['S', 'four', 0], 'TCG': ['S', 'four', 0],
        	'AGT': ['S_', 'two', 0], 'AGC': ['S_', 'two', 0], 'ACT': ['T', 'four', 0],
        	'ACC': ['T', 'four', 0], 'ACA': ['T', 'four', 0], 'ACG': ['T', 'four', 0],
        	'TGG': ['W', 'one', 0], 'TAT': ['Y', 'two', 0], 'TAC': ['Y', 'two', 0],
        	'GTT': ['V', 'four', 0], 'GTC': ['V', 'four', 0], 'GTA': ['V', 'four', 0],
        	'GTG': ['V', 'four', 0], 'TAA': ['*', 'none', 0], 'TGA': ['*', 'none', 0],
        	'TAG': ['*', 'none', 0], 'XXX': ['_missing', 'none', 0]}

        # Universal codon table, with 6-fold degenerate codons kept
        # whole, no splitting! Traditional Universal codon table.
        universal_6fold = {
        	'GCT': ['A', 'four', 0], 'GCC': ['A', 'four', 0], 'GCA': ['A', 'four', 0],
        	'GCG': ['A', 'four', 0], 'CGT': ['R', 'six', 0], 'CGC': ['R', 'six', 0],
        	'CGG': ['R', 'six', 0], 'CGA': ['R', 'six', 0], 'AGA': ['R', 'six', 0],
        	'AGG': ['R', 'six', 0], 'AAT': ['N', 'two', 0], 'AAC': ['N', 'two', 0],
        	'GAT': ['D', 'two', 0], 'GAC': ['D', 'two', 0], 'TGT': ['C', 'two', 0],
        	'TGC': ['C', 'two', 0], 'CAA': ['Q', 'two', 0], 'CAG': ['Q', 'two', 0],
        	'GAA': ['E', 'two', 0], 'GAG': ['E', 'two', 0], 'GGT': ['G', 'four', 0],
        	'GGC': ['G', 'four', 0], 'GGA': ['G', 'four', 0], 'GGG': ['G', 'four', 0],
        	'CAT': ['H', 'two', 0], 'CAC': ['H', 'two', 0], 'ATT': ['I', 'three', 0],
        	'ATC': ['I', 'three', 0], 'ATA': ['I', 'three', 0], 'ATG': ['M', 'one', 0],
        	'TTA': ['L', 'six', 0], 'TTG': ['L', 'six', 0], 'CTT': ['L', 'six', 0],
        	'CTC': ['L', 'six', 0], 'CTA': ['L', 'six', 0], 'CTG': ['L', 'six', 0],
        	'AAA': ['K', 'two', 0], 'AAG': ['K', 'two', 0], 'TTT': ['F', 'two', 0],
        	'TTC': ['F', 'two', 0], 'CCT': ['P', 'four', 0], 'CCC': ['P', 'four', 0],
        	'CCA': ['P', 'four', 0], 'CCG': ['P', 'four', 0], 'TCT': ['S', 'six', 0],
        	'TCC': ['S', 'six', 0], 'TCA': ['S', 'six', 0], 'TCG': ['S', 'six', 0],
        	'AGT': ['S', 'six', 0], 'AGC': ['S', 'six', 0], 'ACT': ['T', 'four', 0],
        	'ACC': ['T', 'four', 0], 'ACA': ['T', 'four', 0], 'ACG': ['T', 'four', 0],
        	'TGG': ['W', 'one', 0], 'TAT': ['Y', 'two', 0], 'TAC': ['Y', 'two', 0],
        	'GTT': ['V', 'four', 0], 'GTC': ['V', 'four', 0], 'GTA': ['V', 'four', 0],
        	'GTG': ['V', 'four', 0], 'TAA': ['*', 'none', 0], 'TGA': ['*', 'none', 0],
        	'TAG': ['*', 'none', 0], 'XXX': ['_missing', 'none', 0]}

        # Blepharisma (table 4) genetic code codon table, with 6-fold degenerate
        # codons kept whole, no splitting!
        blepharisma_6fold = {**universal_6fold,
            'TGA': ['W', 'two', 0], 'TGG': ['W', 'two', 0],
            'TAA': ['*', 'two', 0], 'TAG': ['*', 'two', 0]}

        # Blepharisma (table 4) genetic code codon table, with 6-fold degenerate
        # codons split into four-fold and two-fold groups.
        blepharisma_no6fold = {**universal_no6fold,
            'TGA': ['W', 'two', 0], 'TGG': ['W', 'two', 0],
            'TAA': ['*', 'two', 0], 'TAG': ['*', 'two', 0]}

        # Chilodonella genetic code codon table, with 6-fold degenerate
        # codons kept whole, no splitting!
        chilo_6fold = {**universal_6fold,
            'CAA': ['Q', 'four', 0], 'CAG': ['Q', 'four', 0],
            'TAA': ['*', 'one', 0], 'TAG': ['Q', 'four', 0],
            'TGA': ['Q', 'four', 0]}

        # Chilodonella genetic code codon table, with 6-fold degenerate
        # codons split into four-fold and two-fold groups.
        # Note that this also splits four-fold degenerate codons that OUGHT to
        # be in "different" functional categories (e.g. CAG =/= TAG)
        chilo_no6fold = {**universal_no6fold,
            'TAA': ['*', 'one', 0], 'TAG': ['Q_', 'one', 0],
            'TGA': ['Q_', 'one', 0]}

        # Ciliate (table 6) genetic code codon table, with 6-fold degenerate
        # codons kept whole, no splitting! Traditional ciliate codon table.
        ciliate_6fold = {**universal_6fold,
            'CAA': ['Q', 'four', 0], 'CAG': ['Q', 'four', 0],
            'TAA': ['Q', 'four', 0], 'TAG': ['Q', 'four', 0],
            'TGA': ['*', 'one', 0]}

        # Ciliate (table 6) genetic code codon table, with 6-fold degenerate
        # codons split into four-fold and two-fold groups.
        # Note that this also splits four-fold degenerate codons that OUGHT to
        # be in "different" functional categories (e.g. CAA =/= TAA)
        ciliate_no6fold = {**universal_no6fold,
            'TAA': ['Q_', 'two', 0], 'TAG': ['Q_', 'two', 0],
            'TGA': ['*', 'one', 0]}

        # Euplotes codon table, with 6-fold degenerate codons kept
        # whole, no splitting! Traditional Universal codon table.
        euplotes_6fold = {**universal_6fold,
            'TGA': ['C', 'three', 0], 'TGT': ['C', 'three', 0],
            'TGC': ['C', 'three', 0], 'TAA': ['*', 'two', 0],
            'TAG': ['*', 'two',0]}

        # Euplotes genetic code codon table, with 6-fold degenerate codons
        # split into four-fold and two-fold groups.
        euplotes_no6fold = {**universal_no6fold,
            'TGA': ['C', 'three', 0], 'TGT': ['C', 'three', 0],
            'TGC': ['C', 'three', 0], 'TAA': ['*', 'two', 0],
            'TAG': ['*', 'two',0]}

        # Mesodinium/Myrionecta (table 29) genetic code codon table, with 6-fold
        # degenerate codons kept whole, no splitting! Traditional ciliate codon table.
        mesodinium_6fold = {**universal_6fold,
            'TAA': ['Y', 'four', 0], 'TAT': ['Y', 'four', 0],
            'TAG': ['Y', 'four', 0], 'TAC': ['Y', 'four', 0],
            'TGA': ['*', 'one', 0]}

        # Mesodinium/Myrionecta (table 29) genetic code codon table, with 6-fold
        # degenerate codons split into four-fold and two-fold groups.
        mesodinium_no6fold = {**universal_no6fold,
            'TAA': ['Y', 'four', 0], 'TAT': ['Y', 'four', 0],
            'TAG': ['Y', 'four', 0], 'TAC': ['Y', 'four', 0],
            'TGA': ['*', 'one', 0]}

        # Peritrich (table 30) genetic code codon table, with 6-fold degenerate
        # codons kept whole, no splitting! Traditional ciliate codon table.
        peritrich_6fold = {**universal_6fold,
            'GAA': ['E', 'four', 0], 'GAG': ['E', 'four', 0],
            'TAA': ['E', 'four', 0], 'TAG': ['E', 'four', 0],
            'TGA': ['*', 'one', 0]}

        # Peritrich (table 30) genetic code codon table, with 6-fold degenerate
        # codons split into four-fold and two-fold groups.
        # Note that this also splits four-fold degenerate codons that OUGHT to
        # be in "different" functional categories (e.g. CAA =/= TAA)
        peritrich_no6fold = {**universal_no6fold,
            'TAA': ['E_', 'two', 0], 'TAG': ['E_', 'two', 0],
            'TGA': ['*', 'one', 0]}

        cdnTableDict = {1:[universal_no6fold,universal_6fold],
            4:[blepharisma_no6fold, blepharisma_6fold],
            6:[ciliate_no6fold,ciliate_6fold],
            10:[euplotes_no6fold,euplotes_6fold],
            29:[mesodinium_no6fold,mesodinium_6fold],
            30:[peritrich_no6fold,peritrich_6fold],
            'chilodonella':[chilo_no6fold,chilo_6fold],
            'chilo':[chilo_no6fold,chilo_6fold]}
        return cdnTableDict[gCode]

    def mapCdns(seq, cdnTable):
        # Updates the codon counts for a given sequence to the respective codon
        # count table (e.g. with or without 6-fold degeneracy).
        codons = [seq[n:n+3] for n in range(0, len(seq)-len(seq)%3, 3)]
        amb_cdn = 0
        for c in codons:
            try:
                cdnTable[c][-1] += 1
            except:
                amb_cdn += 1
        if cdnTable['TCC'][1] == 'six':
            return cdnTable, amb_cdn
        else:
            return cdnTable

class GCeval():
    """
    Returns %GC values from DNA sequences of various types.
    """
    def gcTotal(seq):
        # This function returns global GC content
        return round(GC(seq), 4)

    def gc1(seq):
        # This function return the GC content of the first position of a codon
        return round(GC(''.join([seq[n] for n in range(0, len(seq), 3)])), 4)

    def gc2(seq):
        # This function return the GC content of the second position of a codon
        return round(GC(''.join([seq[n] for n in
            range(1, len(seq)-len(seq[1:]) % 3, 3)])), 4)

    def gc3(seq):
        # This function return the GC content of the third position of a codon
        return round(GC(''.join([seq[n] for n in
            range(2, len(seq)-len(seq[2:]) % 3, 3)])), 4)

    def gc3s(cdnTbl):
        # This function return the GC content of the third position of a codon excluding Tryp and Met
        syn = round(GC(''.join([k[-1]*v[-1] for k, v in cdnTbl.items() if v[0] != 'W' and v[0] != 'M'])), 4)
        
        return syn

    def gc3_4F(cdnTbl):
    #     # This function return the GC content of the third position of four-fold
    #     # degenerate codons

        FrFold = round(GC(''.join([k[-1]*v[-1] for k, v in cdnTbl.items() if v[1] == 'four'])), 4)
        
        return FrFold

class SeqInfo(object):
    """
    Provides a means to harbor the data for each individual contig/gene in a
    given fasta file.
    This includes GC content (various types), Effective Number of codons
    (ENc; again various calculations), Relative Synonymous Codon Usage (RSCU).
    """
    def __init__(self,seq,gcode='universal'):
        self.ntd = str(seq)
        self.gcode, self.transTable = GenUtil.convertGenCode(gcode)
        # Dictionary of the GC-related functions/calculations
        self.gcFuncs = {'gcOverall':GCeval.gcTotal,'gc1':GCeval.gc1,'gc2':GCeval.gc2,'gc3':GCeval.gc3}

    def countCodons(self):
        # Stores the different codon tables and updates their codon counts
        cdnTbls = GenUtil.getCDNtable(self.transTable)
        self.cdnCounts_6F,self.amb_cdn = GenUtil.mapCdns(self.ntd, cdnTbls[1])
        self.cdnCounts_No6F = GenUtil.mapCdns(self.ntd, cdnTbls[0])

    def ENcStats(self):
        # Stores the various Effective Number of Codons calculations in the class
        self.expENc = CalcCUB.expWrightENc(self.gc3s)
        self.obsENc_6F = CalcCUB.calcWrightENc(self.cdnCounts_6F)
        self.obsENc_No6F = CalcCUB.calcWrightENc(self.cdnCounts_No6F)
        self.SunENc_6F = CalcCUB.SunEq5(self.cdnCounts_6F)
        self.SunENc_No6F = CalcCUB.SunEq5(self.cdnCounts_No6F)

    def GCstats(self):
        # Stores the various GC-stats in the class
        for k, v in self.gcFuncs.items():
            setattr(self,k,v(self.ntd))
        self.gc4F = GCeval.gc3_4F(self.cdnCounts_No6F)
        self.gc3s = GCeval.gc3s(self.cdnCounts_No6F)


    def RSCUstats(self):
        self.rscu_No6Fold = CalcCUB.RSCU(self.cdnCounts_No6F)
        self.rscu_6Fold = CalcCUB.RSCU(self.cdnCounts_6F)


def prepFolders(outName):
    if os.path.isdir(outName) == False:
        os.mkdir(outName)
    if os.path.isdir(outName+'/SpreadSheets') == False:
        os.mkdir(outName+'/SpreadSheets')


def CalcRefFasta(fasta, gCode):
    seqDB = {i.description:SeqInfo(i.seq, gCode) for i in SeqIO.parse(fasta,'fasta')}
    GenCDNtable = {}
    for k, v in seqDB.items():
        v.countCodons()
        v.GCstats()
        v.ENcStats()
        for k, v in v.cdnCounts_6F.items():
            if k.isalpha() and k not in GenCDNtable .keys():
                GenCDNtable[k] = [v[0],v[-1]]
            else:
                GenCDNtable[k][-1] += v[-1]
    RSCU = CalcCUB.calcRCSU(GenCDNtable)
    return seqDB, RSCU


def WriteWrightOut(seqData, outName, comp):
    if comp == False:
        with open(outName+'/SpreadSheets/'+outName.split('/')[-1]+'.ENc.Raw.tsv','w+') as w:
            w.write('SequenceID\tAmbiguousCodons\tGC-Overall\tGC1\tGC2\tGC3\tGC3S\t'
                'GC3-Degen\tExpWrightENc\tObsWrightENc_6Fold\tObsWrightENc_No6Fold\t'
                'ObsWeightedENc_6Fold\tObsWeightedENc_No6Fold\n')
            for k, v in seqData.items():
                name = [k]
                gcs = [str(v.gcOverall),str(v.gc1),str(v.gc2),str(v.gc3),str(v.gc3s),str(v.gc4F)]
                ENc = [str(v.expENc),str(v.obsENc_6F),str(v.obsENc_No6F),
                str(v.SunENc_6F),str(v.SunENc_No6F)]
                w.write('\t'.join(name+[str(v.amb_cdn)]+gcs+ENc)+'\n')
    else:
        with open(outName+'/SpreadSheets/'+outName.split('/')[-1]+'.CompTrans.ENc.Raw.tsv','w+') as w:
            w.write('SequenceID\tAmbiguousCodons\tGC-Overall\tGC1\tGC2\tGC3\tGC3S\t'
                'GC3-Degen\tExpWrightENc\tObsWrightENc_6Fold\tObsWrightENc_No6Fold\t'
                'ObsWeightedENc_6Fold\tObsWeightedENc_No6Fold\n')
            for k, v in seqData.items():
                name = [k]
                gcs = [str(v.gcOverall),str(v.gc1),str(v.gc2),str(v.gc3),str(v.gc3s),str(v.gc4F)]
                ENc = [str(v.expENc),str(v.obsENc_6F),str(v.obsENc_No6F),
                str(v.SunENc_6F),str(v.SunENc_No6F)]
                w.write('\t'.join(name+[str(v.amb_cdn)]+gcs+ENc)+'\n')


def getCompFasta(fasta, gCode, require_start, require_stop):
    stopCDNs = {'1':['TAA','TAG','TGA'], '4':['TAA','TAG'], '6':['TGA'], '10':['TAA','TAG'],
    '29':['TGA'], '30':['TGA'], 'universal':['TAA','TAG','TGA'], 'blepharisma':['TAA','TAG'],
    'ciliate':['TGA'],'euplotes':['TAA','TAG'], 'mesodinium':['TGA'], 'peritrich':['TGA'],
    'chilo':['TAA']}
    if gCode.lower() not in stopCDNs.keys():
        print('\nWARNING: Genetic code not recognized. Defaulting to universal.')
        stops = stopCDNs['1']
    else:
        stops = stopCDNs[gCode]

    with open(fasta.replace('.fasta','.Comp.fasta'),'w+') as w:
        for i in SeqIO.parse(fasta,'fasta'):
            if len(i.seq) % 3 == 0:
                if ((require_start and str(i.seq).upper().startswith('ATG')) or not require_start) and ((require_stop and str(i.seq).upper()[-3:] in stops) or not require_stop):
                    w.write('>'+i.description+'\n'+str(i.seq)+'\n')

    return fasta.replace('.fasta','.Comp.fasta')


def WriteNullENcOut(outName):
    with open(outName+'/SpreadSheets/' + outName.split('/')[-1] + '.ENc.Null.tsv','w+') as w:
        w.write('GC3S\tENc\n')
        w.write('\n'.join(CalcCUB.nullENcGC3()))


def WriteRSCUtbl(RSCUtbl, outName):
    with open(outName+'/SpreadSheets/' + outName.split('/')[-1] + '.RSCU.tsv','w+') as w:
        w.write('Codon\tAmino Acid\tRSCU\n')
        for k,v in RSCUtbl.items():
            w.write(k+'\t'+'\t'.join(v)+'\n')


def get_args():

    parser = argparse.ArgumentParser(
                prog = 'Codon Usage Bias Statistic Calculator version 2.1',
                description = "Written by Xyrus Maurer-Alcala, updated July 19, 2023 by Auden Cote-L'Heureux"
                )

    parser.add_argument('--input', '-i', required = True, default = None, type = str, help = 'Path to a fasta file of nucleotide sequences OR to a folder containing multiple fasta files')
    parser.add_argument('--genetic_code', '-g', default = 'universal', choices = { 'universal', 'blepharisma', 'ciliate', 'euplotes', 'mesodinium', 'peritrich', 'chilo'}, help = 'Genetic code to use (be careful if you are working with ciliates). Available codes are:  "universal" (TAA, TAG, TGA), "blepharisma" (TAA, TAG), "ciliate" (TGA), "euplotes" (TAA, TAG), "mesodinium" (TGA), "peritrich" (TGA) and "chilo" (TAA). Ask someone if you are not sure what this means.')
    parser.add_argument('--require_start', action = 'store_true', help = 'Filter for sequences that begin with a start codon')
    parser.add_argument('--require_stop', action = 'store_true', help = 'Filter for sequences that end with a stop codon')

    return parser.parse_args()


if __name__ == "__main__":

    args = get_args()

    #If only one file was input
    if os.path.isfile(args.input) and args.input.split('.')[-1] in ('fasta', 'fas', 'fa', 'fna'):
        outName = args.input.split('.')[0]

        compFasta = getCompFasta(args.input, args.genetic_code, args.require_start, args.require_stop)
        prepFolders(outName)
        fastaDataRaw, RSCUtbl = CalcRefFasta(args.input, args.genetic_code)
        fastaDataComp, RSCUtbl = CalcRefFasta(compFasta, args.genetic_code)
        WriteWrightOut(fastaDataRaw, outName, comp=False)
        WriteWrightOut(fastaDataComp, outName, comp=True)
        WriteNullENcOut(outName)
        WriteRSCUtbl(RSCUtbl, outName)
        os.system('mv ' + compFasta + ' ' + outName + '/')
    #If a folder of files was input
    elif os.path.isdir(args.input):
        if not os.path.isdir('TemporaryCUBOutput'):
            os.mkdir('TemporaryCUBOutput')
        for file in os.listdir(args.input):
            if file != '.DS_Store':
                if os.path.isfile(args.input + '/' + file) and file.split('.')[-1] in ('fasta', 'fas', 'fa', 'fna'):
                    outName = 'TemporaryCUBOutput/' + file.split('.')[0]

                    compFasta = getCompFasta(args.input + '/' + file, args.genetic_code, args.require_start, args.require_stop)
                    prepFolders(outName)
                    fastaDataRaw, RSCUtbl = CalcRefFasta(args.input + '/' + file, args.genetic_code)
                    fastaDataComp, RSCUtbl = CalcRefFasta(compFasta, args.genetic_code)
                    WriteWrightOut(fastaDataRaw, outName, comp=False)
                    WriteWrightOut(fastaDataComp, outName, comp=True)
                    WriteNullENcOut(outName)
                    WriteRSCUtbl(RSCUtbl, outName)
                    os.system('mv ' + compFasta + ' ' + outName + '/')
                else:
                    print('\nWARNING: One of the files in your input folder (' + file + ') is not formatted correctly. These should be nucleotide fasta files with a file extension fasta, fas, fa, or fna. Skipping this file.\n')

        folders = ['TemporaryCUBOutput/' + folder for folder in os.listdir('TemporaryCUBOutput') if os.path.isfile('TemporaryCUBOutput/' + folder + '/SpreadSheets/' + folder + '.CompTrans.ENc.Raw.tsv')]

        #Combining information from all input files into one set of outputs
        if len(folders) > 0:
            if not os.path.isdir('CUBOutput'):
                os.mkdir('CUBOutput')

            if not os.path.isdir('CUBOutput/SpreadSheets'):
                os.mkdir('CUBOutput/SpreadSheets')

            os.system('cat TemporaryCUBOutput/*/*.Comp.fasta > CUBOutput/AllFiles.Comp.fasta')
            os.system('cp ' + folders[0] + '/SpreadSheets/*.ENc.Null.tsv CUBOutput/SpreadSheets/ENc.Null.tsv')

            with open('CUBOutput/SpreadSheets/RSCU.tsv', 'w') as o:
                o.write('File\tCodon\tAmino Acid\tRSCU\n')
                for folder in folders:
                    for line in open(folder + '/SpreadSheets/' + folder.split('/')[-1] + '.RSCU.tsv'):
                        if 'Amino Acid' not in line:
                            o.write(folder.split('/')[-1] + '\t' + line)

            with open('CUBOutput/SpreadSheets/ENc.Raw.tsv', 'w') as o:
                o.write('File\tSequenceID\tAmbiguousCodons\tGC-Overall\tGC1\tGC2\tGC3\tGC3S\tGC3-Degen\tExpWrightENc\tObsWrightENc_6Fold\tObsWrightENc_No6Fold\tObsWeightedENc_6Fold\tObsWeightedENc_No6Fold\n')
                for folder in folders:
                    for line in open(folder + '/SpreadSheets/' + folder.split('/')[-1] + '.ENc.Raw.tsv'):
                        if 'SequenceID' not in line:
                            o.write(folder.split('/')[-1] + '\t' + line)

            with open('CUBOutput/SpreadSheets/CompTrans.ENc.Raw.tsv', 'w') as o:
                o.write('File\tSequenceID\tAmbiguousCodons\tGC-Overall\tGC1\tGC2\tGC3\tGC3S\tGC3-Degen\tExpWrightENc\tObsWrightENc_6Fold\tObsWrightENc_No6Fold\tObsWeightedENc_6Fold\tObsWeightedENc_No6Fold\n')
                for folder in folders:
                    for line in open(folder + '/SpreadSheets/' + folder.split('/')[-1] + '.CompTrans.ENc.Raw.tsv'):
                        if 'SequenceID' not in line:
                            o.write(folder.split('/')[-1] + '\t' + line)

            os.system('rm -r TemporaryCUBOutput')
        else:
            print('\nERROR: No composition information was created -- something probably went wrong with the formatting and/or filtering of the input fasta files. Make sure that this contains properly formatted nucleotide sequences (divisible by 3, etc).')
    else:
        print('\nERROR: Invalid --input. This should be a nucleotide fasta file with a file extension fasta, fas, fa, or fna or a folder of such files.\n')











