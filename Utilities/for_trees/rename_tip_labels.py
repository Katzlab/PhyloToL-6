
''' Author, Date : Godwin Ani, 15th - September - 2023.
Motivation : To make phylogenetic trees more presentable.
Intent : Rename the tip labels of phylogenetic trees.
Dependencies : Python3, Pandas
Inputs : A folder containing trees and a csv file(with headers).
The first column of the csv is the 10 digit code and other columns are the information to be added to the tip labels.
Outputs : A folder of trees with renamed tips.
python3 rename_tip_labels.py -i input to_folder_of_trees -s to_spreadsheet
'''

import os, re, sys, argparse
import pandas as pd
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input')
parser.add_argument('-s', '--spreadsheet')
args = parser.parse_args()
os.makedirs(args.input + '/renamed', exist_ok = True)
df = pd.read_csv(args.spreadsheet, index_col = 0)
df = df.astype(str)
df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)
df['Merged'] = df.apply(lambda row: '_'.join(row), axis=1)
for file in os.listdir(args.input):
    if file.endswith('.tree') or file.endswith('.tre'):
        with open(args.input + '/' + file, 'r') as tree:
            tree = tree.read()
            tree = tree.replace('Len_', 'L')
            tree = tree.replace('Cov_', 'Cv')
            tree = tree.replace('Contig_', 'Ct')
            tree = tree.replace('Len', 'L')
            tree = tree.replace('Cov', 'Cv')
            tree = tree.replace('Contig', 'Ct')
            search_strings = df.index.tolist()
            replacement_strings = df['Merged'].tolist()
            for search, replace in zip(search_strings, replacement_strings):
                tree = tree.replace(search,replace)
                with open(args.input + '/renamed/' + file, 'w') as o:
                    o.write(tree)


