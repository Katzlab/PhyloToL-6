library(tidyverse)

#Change this path
setwd('/Your/Working/Directory')

#You will need to change name of data frame below
#if you used the CUB.py script, it should be in the spreadsheets folder
#in the output and end in CompTrans.ENc.Raw.tsv
gc3 <- data.frame(read_tsv('CompTrans.ENc.Raw.tsv'))

gc3$GC3.Degen <- as.numeric(gc3$GC3.Degen)
gc3$ObsWrightENc_6Fold <- as.numeric(gc3$ObsWrightENc_6Fold)

#The data for the null expectation curve will be in the same folder as above
enc_null <- data.frame(read_tsv('ENc.Null.tsv'))

gc3_plot <- ggplot(data = gc3, mapping = aes(as.numeric(GC3.Degen), as.numeric(ObsWrightENc_6Fold))) +
   geom_point() +
   geom_line(data = enc_null, aes(GC3, ENc)) +
   theme_classic() +
   labs(x = '%GC at 3rd-pos 4-fold sites', y = 'Observed Wright ENc (6Fold)') +
   theme(
      legend.position = 'none'
   )

gc3_plot
