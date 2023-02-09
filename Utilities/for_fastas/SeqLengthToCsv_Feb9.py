'''
Professor L. Katz and Godwin Ani
9th-Feb-2023
Seq_length_to_csv is a program that exports the length of DNA sequences excluding gaps and missing data to a csv file.
'''

import os, sys, re
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import pandas as pd
from tqdm import tqdm

'''Ensure that DNA fasta files are stored in a root directory folder named "Seq_length". This folder will automatically be created if none exists.
	The result of the process is saved into a sub-folder within the "Seq_length" folder named "csv" which is created automatically if none exists. '''
	
	
def Seq_length_to_csv():
    if not os.path.isdir('Seq_length'):
        os.mkdir('Seq_length')
    
    if not os.path.isdir('Seq_length/csv/'):
        os.mkdir('Seq_length/csv/')
        
    list = []
    for file in os.listdir('Seq_length'):
        if file.endswith('.fasta'):
            list.append(file)
            
    for x in tqdm(list):
        name =[]
        seq_length = []
        file_name = x.split('fasta')[0] + '.csv'
        for x in SeqIO.parse('Seq_length/' + x, "fasta"):
            name.append(x.id)
            seq = x.seq
            seq_length.append(seq.count("A") + seq.count("a") + seq.count("T") + seq.count("t") + seq.count("G") + seq.count("g") + seq.count("C") + seq.count("c"))
            a = np.array([seq_length])
            mean_length = round(np.mean(a), 2)
            data = {'Name of sequence' : name, 'Length of sequence' : seq_length, 'Average sequence length' : mean_length}
            df = pd.DataFrame(data)
            df.index += 1
            df.to_csv('Seq_length/csv/' + file_name)


Seq_length_to_csv()


