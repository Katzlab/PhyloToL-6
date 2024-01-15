#Professor L. Katz and Godwin Ani
#9th- Feb- 2023
#A clustering program that accepts and validates users input

import os, sys, re
from tqdm import tqdm


def cluster():
    ''' This function takes an amino acid or DNA sequence stored in a fasta format and clusters it.
It uses two nested functions (cluster_AA and cluster_DNA) to perform this operation.
Ensure that fasta files to be clustered are stored in a root directory folder named "to_cluster". This folder will automatically be created if none exists.
The result of the clustering process is saved into a root directory folder named "Clustered" which is created automatically if none exists. '''
  
    if not os.path.isdir('to_cluster'):
        os.mkdir('to_cluster')
        
    if not os.path.isdir('Clustered'):
        os.mkdir('Clustered')
    # Nested function for amino acids clustering.    
    def cluster_AA():
    #input validation for the sequence identity threshold.
        while True:
            try:
                val1 =input( 'Amino Acids Sequence Identity Threshold (e.g. 0.99, 0.95)? : ')
                integer, fractional = val1.split('.')
                val1 = float(val1)
                if int(integer)== 0 and len(fractional) == 2:
                    break
            except ValueError:
                pass
            print('ERROR! Use format 0.## for Amino acids sequence identity threshold.')
 
      #Input validation for the overlap value.
        while True:
            try:
                val2 =input( 'Amino Acids Alignment Overlap Value (e.g. 0.67, 0.75)? : ')
                integer, fractional = val2.split('.')
                val2 = float(val2)
                if int(integer)== 0 and len(fractional) == 2:
                    break
            except ValueError:
                pass
            print('ERROR! Use format 0.## for Amino acids sequence alignment overlap value')
            
  #Selects amino acids fasta files in "to_cluster" folder and clusters them with CD-HIT.
        for file in tqdm(os.listdir('to_cluster')):
            if file.endswith('.fasta'):
                os.system('cd-hit -i to_cluster/' + file + ' -o Clustered/' + file + ' -c %s -d 0 -aS %s' %( val1, val2))
     #Renaming the result of the clustering.
        for file in os.listdir('Clustered'):
            if file.endswith('.clstr'):
                os.rename('Clustered/' + file, 'Clustered/' + file.split('FILE')[0] + 'Clustered.txt')


    #Nested function for DNA clustering.
    def cluster_DNA():
        #Input validation for the sequence identity threshold.
        while True:
            try:
                val1 =input( 'DNA Sequence Identity Threshold (e.g. 0.99, 0.95)? : ')
                integer, fractional = val1.split('.')
                val1 = float(val1)
                if int(integer)== 0 and len(fractional) == 2:
                    break
            except ValueError:
                pass
            print('ERROR! Use format 0.## for DNA sequence identity threshold.')
      #Input validation for the overlap value. 
        while True:
            try:
                val2 =input( 'DNA Sequence Alignment Overlap Value (e.g. 0.67, 0.75)? : ')
                integer, fractional = val2.split('.')
                val2 = float(val2)
                if int(integer)== 0 and len(fractional) == 2:
                    break
            except ValueError:
                pass
            print('ERROR! Use format 0.## for DNA sequence alignment overlap value.')
            
            #Selects DNA fasta files in "to_cluster" folder and clusters them with CD-HIT.
        for file in tqdm(os.listdir('to_cluster')):
            if file.endswith('.fasta'):
                os.system('cd-hit-est -i to_cluster/' + file + ' -o Clustered/' + file + ' -c %s -d 0 -aS %s' %( val1, val2))

        #Renaming the result of the clustering.
        for file in os.listdir('Clustered'):
            if file.endswith('.clstr'):
                os.rename('Clustered/' + file, 'Clustered/' + file.split('FILE')[0] + 'Clustered.txt')

    # Prompts for user input and function call.
    choice_1 = input('Are you clustering Amino Acids or DNA? (AA or DNA) : ')
    if choice_1 in ['AA', 'Aa', 'aa']:
        cluster_AA()
    elif choice_1 in ['DNA', 'Dna', 'dna']:
        cluster_DNA()
    else:
        print('Sorry. This program can only cluster Amino Acids and DNA')

cluster()

