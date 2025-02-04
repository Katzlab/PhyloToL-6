
''' Author, Date : Godwin Ani, 10 - July - 2024.
Motivation : To make phylogenetic trees more presentable.
Intent : Shorten the tip labels of phylogenetic trees.
Dependencies : Python3, ete3
Inputs : A folder containing trees
Outputs : A folder of trees with shortened tips.
python3 RenameTips_v1.0.py -i input to_folder_of_trees 
'''


import os, re, sys, argparse, string
import ete3



parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input')
args = parser.parse_args()
os.makedirs(args.input + '/renamed', exist_ok = True)

def get_newick(fname):
    newick = ''
    for line in open(fname):
        line = line.split(' ')[-1]
        if(line.startswith('(') or line.startswith('tree1=')):
            newick = line.split('tree1=')[-1].replace("'", '').replace('\\', '')
    return newick


def tree_formatting_wrapper(file):
    newick = get_newick(file)
    tree = ete3.Tree(newick)    
    any_letter = tuple(string.ascii_letters)
    for leaf in tree:
        if leaf.name.startswith(any_letter):
            leaf.name = str(leaf.name).split('_Len')[0]
            leaf.name = str(leaf.name).replace('Contig_', 'Ct')
            leaf.name = str(leaf.name).replace('_XX_0', '')
    tree.write(format=1, outfile=args.input + '/renamed/' +file.split('/')[-1] + '.tree')


for tree in os.listdir(args.input):
    if tree.split('.')[-1] in ('tree', 'tre', 'treefile', 'nex'):
        tree_formatting_wrapper(args.input + '/' + tree)





