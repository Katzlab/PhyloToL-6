#Author, date: Written by Claude, modified by Ozan, September 30 2026
#Intent: Calculate pairwise uncorrected p-distances between all sequences in an aligned fasta file or folder of aligned fasta files
#Dependencies: Python3, numpy, BioPython
#Inputs: Aligned fasta file (amino acid or nucleotide) or folder of aligned fasta files
#Outputs: A spreadsheet (tsv) of pairwise p-distances with three columns: first sequence, second sequence, pairwise distance
#Example: python3 PairwiseDistance.py -i alignment.fasta (or alignment_DIR)
#Note: p-distance = number of differing sites/number of compared sites (no substitution model or multiple-hit correction)
#Note: Sequence type (AA vs NTD) is detected automatically; use --seq_type to override if detection is wrong 
#Note: Gaps and ambiguous characters are handled by pairwise deletion by default; use --complete_deletion to drop any column with a gap/ambiguity in ANY sequence
#Note: Distances are reported as NA when two sequences share no comparable (non-gap, non-ambiguous) sites
#Note: Files with duplicate sequence headers are skipped, since the output would be ambiguous
#Note: Use "python3 PairwiseDistance.py -i GeneFolder --concatenate" to also get distances for the concatenated alignment of all genes. Header lines should be identical between gene files to be concatenated (PairwiseOutput/Concatenated.PairwiseDist.tsv). Concatenated distances are summed differences/summed compared sites across genes (NOT an average of per-gene distances).

#Dependencies
import os
import sys
import shutil
import argparse
import tempfile
from collections import Counter
import numpy as np
from Bio import SeqIO


class SeqUtil(object):
    """
    Functions to read in the alignment, check that headers are unique and that
    sequences are aligned, and determine whether the alignment is composed of
    amino acid (AA) or nucleotide (NTD) sequences.
    """
    # Characters counted as nucleotides when auto-detecting the sequence type
    NTD_CHARS = set('ACGTUN')
    # Characters ignored when auto-detecting the sequence type (gaps/missing/stops)
    IGNORE_CHARS = set('-.?*')
    # Characters treated as valid (comparable) sites. Everything else (gaps,
    # IUPAC ambiguity codes, X, etc.) is treated as missing data.
    VALID_NTD = 'ACGT'
    VALID_AA = 'ACDEFGHIKLMNPQRSTVWY'

    @staticmethod
    def readAlignment(fasta):
        # Returns a list of header and sequence, keeping the full header
        # (tabs in headers are replaced with spaces to keep the tsv columns intact)
        return [(i.description.replace('\t', ' '), str(i.seq).upper()) for i in SeqIO.parse(fasta, 'fasta')]

    @staticmethod
    def checkDuplicates(records):
        # Returns a list of headers that appear more than once
        counts = Counter([r[0] for r in records])
        return [k for k, v in counts.items() if v > 1]

    @staticmethod
    def checkAligned(records):
        # All sequences in an alignment should be the same length
        return len(set([len(r[1]) for r in records])) == 1

    @staticmethod
    def detectSeqType(records, threshold=0.9):
        # Calls the alignment NTD if >= 90% of the non-gap characters are
        # A, C, G, T, U or N; otherwise AA
        chars = Counter(''.join([r[1] for r in records]))
        total = sum([v for k, v in chars.items() if k not in SeqUtil.IGNORE_CHARS])
        if total == 0:
            return None
        ntd = sum([v for k, v in chars.items() if k in SeqUtil.NTD_CHARS])
        return 'ntd' if ntd/float(total) >= threshold else 'aa'


class DistCalc(object):
    """
    Returns pairwise uncorrected p-distances between all sequences of an
    alignment: the proportion of compared sites at which two sequences differ.
    """
    @staticmethod
    def encodeAlignment(records, seqType):
        # Converts the sequences into a numpy array (one row per sequence) and
        # builds a matching mask of valid (comparable) sites
        alphabet = SeqUtil.VALID_NTD if seqType == 'ntd' else SeqUtil.VALID_AA
        names = [r[0] for r in records]
        seqs = [r[1].replace('U', 'T') if seqType == 'ntd' else r[1] for r in records]
        arr = np.vstack([np.frombuffer(s.encode('ascii', 'replace'), dtype=np.uint8) for s in seqs])
        valid = np.isin(arr, np.frombuffer(alphabet.encode(), dtype=np.uint8))
        return names, arr, valid

    @staticmethod
    def completeDeletion(valid):
        # Masks every column that has a gap/ambiguity in any sequence
        return valid & valid.all(axis=0)

    @staticmethod
    def calcPairwise(names, arr, valid):
        # Compares each sequence against all following sequences at once.
        # Returns (name1, name2, p-distance, differences, compared sites)
        results = []

        for i in range(len(names)-1):
            # Sites where both sequences have a valid residue (pairwise deletion)
            shared = valid[i] & valid[i+1:]
            nSites = shared.sum(axis=1)
            # Shared sites where the two residues differ (Hamming distance)
            nDiffs = ((arr[i] != arr[i+1:]) & shared).sum(axis=1)

            # p-distance; pairs with no shared sites give 0/0 = nan (reported as NA)
            with np.errstate(divide='ignore', invalid='ignore'):
                dist = nDiffs/nSites

            for j, d in enumerate(dist):
                results.append((names[i], names[i+1+j], d, int(nDiffs[j]), int(nSites[j])))

        return results

    @staticmethod
    def addToConcat(concat, distData, geneName):
        # Adds one gene's differences and compared sites to the running totals
        # for each pair. Pairs are stored in sorted order so that (A, B) and
        # (B, A) from different genes are summed together.
        for seq1, seq2, d, nDiffs, nSites in distData:
            if nSites == 0:
                continue
            key = tuple(sorted((seq1, seq2)))
            if key not in concat:
                concat[key] = [0, 0, []]
            concat[key][0] += nDiffs
            concat[key][1] += nSites
            concat[key][2].append(geneName)
        return concat

# Function for creating the output folder
def prepFolders(outName):
    if os.path.isdir(outName) == False:
        os.mkdir(outName)

# Function to format distances for output; pairs with no shared comparable sites become NA
def formatDist(d):
    # Pairs with no shared sites (nan) become NA
    if not np.isfinite(d):
        return 'NA'
    # Fixed-point formatting avoids scientific notation (e.g. 1e-06) for small distances
    return f'{float(d):.6f}'

# Function to calculate pairwise distances
def CalcDistFasta(fasta, args):
    # Runs all checks and calculations for a single fasta file. Returns the
    # list of pairwise distances, or None if the file should be skipped.
    records = SeqUtil.readAlignment(fasta)

    if len(records) < 2:
        print('\nWARNING: ' + fasta + ' contains fewer than two sequences. Skipping this file.\n')
        return None

    dups = SeqUtil.checkDuplicates(records)
    if len(dups) > 0:
        print('\nERROR: ' + fasta + ' contains duplicate sequence headers. Please make the following headers unique. Skipping this file.\n')
        print('\n'.join(dups) + '\n')
        return None

    if not SeqUtil.checkAligned(records):
        print('\nERROR: Sequences in ' + fasta + ' are not all the same length. Is this an alignment? Skipping this file.\n')
        return None

    if args.seq_type == 'auto':
        seqType = SeqUtil.detectSeqType(records)
        if seqType is None:
            print('\nERROR: ' + fasta + ' contains no sequence characters. Skipping this file.\n')
            return None
        print('\nDetected sequence type for ' + fasta + ': ' + seqType.upper())
    else:
        seqType = args.seq_type

    names, arr, valid = DistCalc.encodeAlignment(records, seqType)
    if args.complete_deletion:
        valid = DistCalc.completeDeletion(valid)
        if valid.sum() == 0:
            print('\nWARNING: No columns remain in ' + fasta + ' after complete deletion; all distances will be NA.\n')

    return DistCalc.calcPairwise(names, arr, valid), seqType


def WriteDistOut(distData, outName):
    with open(outName+'/'+outName.split('/')[-1]+'.PairwiseDist.tsv','w+') as w:
        w.write('Sequence1\tSequence2\tDistance\n')
        for seq1, seq2, d, nDiffs, nSites in distData:
            w.write(seq1+'\t'+seq2+'\t'+formatDist(d)+'\n')


def WriteConcatOut(concat, outFile):
    # Distance for the concatenated alignment = summed differences / summed
    # compared sites across all genes in which both sequences are present
    with open(outFile,'w+') as w:
        w.write('Sequence1\tSequence2\tDistance\tDifferences\tComparedSites\tNumGenes\tGenes\n')
        for (seq1, seq2), (nDiffs, nSites, genes) in sorted(concat.items()):
            w.write('\t'.join([seq1, seq2, formatDist(nDiffs/float(nSites)),
                str(nDiffs), str(nSites), str(len(genes)), ','.join(genes)])+'\n')

# Function to parse arguments
def get_args():

    parser = argparse.ArgumentParser(
                prog = 'Pairwise Distance Calculator',
                description = "Calculates pairwise uncorrected p-distances for AA or NTD alignments."
                )

    parser.add_argument('--input', '-i', required = True, type = str, help = 'Path to an aligned fasta file OR to a folder containing multiple aligned fasta files')
    parser.add_argument('--seq_type', '-t', default = 'auto', choices = {'auto', 'aa', 'ntd'}, help = 'Sequence type. Default "auto" detects AA vs NTD from the sequence composition.')
    parser.add_argument('--concatenate', '-c', action = 'store_true', help = 'Folder input only: also calculate distances for the concatenation of all gene alignments, by summing differences and compared sites per pair across genes')
    parser.add_argument('--complete_deletion', action = 'store_true', help = 'Remove every alignment column containing a gap/ambiguous character in ANY sequence (default is pairwise deletion)')

    return parser.parse_args()


def main():

    args = get_args()
    extensions = ('fasta', 'fas', 'fa', 'fna', 'faa', 'afa')

    #If only one file was input
    if os.path.isfile(args.input) and args.input.split('.')[-1] in extensions:
        if args.concatenate:
            print('\nWARNING: --concatenate needs a folder of gene alignments; ignoring it for a single file.\n')

        outName = os.path.splitext(args.input)[0]

        out = CalcDistFasta(args.input, args)
        if out is None:
            sys.exit(1)
        prepFolders(outName)
        WriteDistOut(out[0], outName)
    #If a folder of files was input
    elif os.path.isdir(args.input):
        # A fresh temporary folder for every run, so leftovers from a previous
        # (e.g. crashed) run are never merged into the output
        tmpDir = tempfile.mkdtemp(prefix = 'PairwiseTmp_')
        # Running totals for the concatenated alignment (--concatenate)
        concat, seqTypes = {}, set()
        for file in os.listdir(args.input):
            if file != '.DS_Store':
                if os.path.isfile(args.input + '/' + file) and file.split('.')[-1] in extensions:
                    outName = tmpDir + '/' + os.path.splitext(file)[0]

                    out = CalcDistFasta(args.input + '/' + file, args)
                    if out is None:
                        continue
                    distData, seqType = out
                    prepFolders(outName)
                    WriteDistOut(distData, outName)
                    if args.concatenate:
                        seqTypes.add(seqType)
                        concat = DistCalc.addToConcat(concat, distData, os.path.splitext(file)[0])
                else:
                    print('\nWARNING: One of the files in your input folder (' + file + ') is not formatted correctly. These should be aligned fasta files with a file extension ' + ', '.join(extensions) + '. Skipping this file.\n')

        folders = [tmpDir + '/' + folder for folder in os.listdir(tmpDir) if os.path.isfile(tmpDir + '/' + folder + '/' + folder + '.PairwiseDist.tsv')]

        #Combining information from all input files into one set of outputs
        if len(folders) > 0:
            if not os.path.isdir('PairwiseOutput'):
                os.mkdir('PairwiseOutput')

            with open('PairwiseOutput/PairwiseDist.tsv', 'w') as o:
                o.write('File\tSequence1\tSequence2\tDistance\n')
                for folder in folders:
                    with open(folder + '/' + folder.split('/')[-1] + '.PairwiseDist.tsv') as f:
                        for n, line in enumerate(f):
                            if n > 0:
                                o.write(folder.split('/')[-1] + '\t' + line)

            for folder in folders:
                shutil.copy(folder + '/' + folder.split('/')[-1] + '.PairwiseDist.tsv', 'PairwiseOutput/')

            #Distances for the concatenated alignment of all genes
            if args.concatenate:
                if len(seqTypes) > 1:
                    print('\nWARNING: Your input folder contains both AA and NTD alignments. The concatenated distances mix both types.\n')
                if len(concat) == 0:
                    print('\nWARNING: No sequence pairs share comparable sites across genes; no concatenated distances were written. Check that headers match between gene files.\n')
                else:
                    WriteConcatOut(concat, 'PairwiseOutput/Concatenated.PairwiseDist.tsv')

            shutil.rmtree(tmpDir)
        else:
            shutil.rmtree(tmpDir)
            print('\nERROR: No pairwise distances were calculated -- something probably went wrong with the formatting of the input fasta files. Make sure these are aligned sequences with unique headers.')
            sys.exit(1)
    else:
        print('\nERROR: Invalid --input. This should be an aligned fasta file with a file extension ' + ', '.join(extensions) + ' or a folder of such files.\n')
        sys.exit(1)


if __name__ == "__main__":
    main()
