'''
#Author, date: Auden Cote-L'Heureux, last updated in September 2023 by Adri
#Intent: To make a simple jellyfish plot from a dataset like the one produced by the CUB_v2.1.py script, also in the Github Utilities folder. 
#Input: CompTrans.ENc.Raw.tsv and ENc.Null.tsv
#Output: Jellyfish plots
'''

#load necessary packages
library(tidyverse)

#Change to the path of the directory you're working from
#Use "getwd()" in console below to get path
setwd("/Users/agrow/Desktop/NewCUB_meta/CUBOutput/SpreadSheets")

#You may need to change name of data frame below
#if you used the CUB_v2.1.py script from Github, it should be in the 
#CUBOutput folder and then inside the SpreadSheets folder
#you are looking for the CompTrans.ENc.Raw.tsv
gc3 <- data.frame(read_tsv('CompTrans.ENc.Raw.tsv'))%>%
  mutate(taxon = paste(substr(SequenceID, 1, 4), substr(SequenceID,6,10), sep = '')) #this line reads in your 10-digit codes to a column in the data frame called taxon

#This .tsv is generated by the CUB script and will be in the same folder as the .tsv above
#This generates the null expectation line
enc_null <- data.frame(read_tsv('ENc.Null.tsv'))

#change data in first line here to what you want plotted
#you need as.numeric to ensure R is reading the variable correctly
gc3_plot <- ggplot(gc3, aes(as.numeric(GC3.Degen), as.numeric(ObsWrightENc_No6Fold)))+
  geom_point(size = 0.1)+
  geom_line(data = enc_null, aes(GC3, ENc))+
  theme_classic()+
  labs(x = 'GC3 Degen', y = 'ObsWrightENc_No6fold')+
  theme(legend.position = 'none')+
  ggtitle("Metatranscriptomics R2G NTD files")+
  theme(axis.text.x=element_text(colour="black"))+
  theme(axis.text.y=element_text(colour="black"))+
  facet_wrap(vars(taxon), scales = 'free_x')
gc3_plot