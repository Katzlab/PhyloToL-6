'''
#Author, date: Uploaded by Adri Grow, 2023
#Intent: map a group of trimmed reads to a reference. 
#Dependencies: Python3, hisat2, samtools, sambamba
#Inputs: Folder named 'TrimmedReads' containing all the trimmed reads.
#Outputs: Folders with the names of the LKHs containing the sam/bam files.
#Example: python readmapping.py
'''

import os
from Bio import SeqIO

#this first command builds your reference with Hisat.
#If you've already done this, DON'T run this command! Instead, comment it out (use a # in front of it).
#It will output several files. Don't worry about them, Hisat will know what to do.
os.system("hisat2-build NEW_Foram_only_053023.fasta Foram_Index") #change to your reference.fasta and rename the index

folder = os.listdir("TrimmedReads") #Insert the name of the folder which has your trimmed reads inside the quotes
folder.sort() #This sorts the folder so that all the LKHs are in order.

for x in folder:
	if "LKH" in x and "FPE" in x: #assigning a variable to forward reads. Make sure you have both forward and reverse reads for each cell!
		FPE = x
	if "LKH" in x and "RPE" in x: #assigning a variable to reverse reads.
		RPE = x

		if(FPE[:7] == RPE[:7]): 
			#The next few lines are several Hisat commands that will create new files.
			#EDIT the name of the index and the name of the trimmed reads folder in the first command below
			os.system("hisat2 -x Foram_Index -1 TrimmedReads/" +FPE+ " -2 TrimmedReads/" +RPE+ " -S sample.sam") 
			os.system("samtools view -bS sample.sam > sample.bam")
			os.system("samtools fixmate -O bam sample.bam  fixmate_sample.bam")
			os.system("samtools sort -O bam -o sorted_sample.bam fixmate_sample.bam")
			os.system("sambamba markdup -r sorted_sample.bam sorted_sample.dedup.bam")
			os.system("samtools view -h -b -q 40 sorted_sample.dedup.bam > sorted_sample.q40.bam")
			os.system("samtools view -h -b -q 20 sorted_sample.dedup.bam > sorted_sample.q20.bam")
			os.system("samtools view -h -F 4 -b sorted_sample.dedup.bam > defaultparameters_sample.bam")

			if not os.path.isdir(x[:7]):
				os.mkdir(x[0:7]) #making folders with the names of the LKHs
				
			for file in os.listdir('.'): #These lines move the sam/bam files that Hisat creates into the new LKH folders.
				if(file.endswith('.sam') or file.endswith('.bam')):
					os.rename(file,x[:7] + '/' + file)
				
print("~~~~~~~~~~~:>~") #When the snake appears, your script has run!
		
	
	
