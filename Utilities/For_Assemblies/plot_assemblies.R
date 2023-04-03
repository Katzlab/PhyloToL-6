library(ggplot2)
library(tidyverse)
library(scales)

#read in data output by python script
transcriptome_data <- read.csv('assembly_assessment.csv')

#make a column combining ten digit code and taxon information
transcriptome_data_expert_level <- transcriptome_data|>
  mutate(taxon_id = paste(ten_digit_code, taxon_info, sep = "-"))

#plot length vs coverage, faceted and labeled by taxon
len_cov <- ggplot(transcriptome_data, aes(x = GC, y = length, color = taxon_info))+
  geom_point(size = .5)+
  ylim(200, 15000)+
  facet_wrap(~ten_digit_code + taxon_info)+
  ggtitle('GC% and length of Allogromia assembled transcripts')+
  theme(strip.text = element_text(size = 7))

ggsave('len_cov.png', device = 'png', width = 8.5, height = 6)#save plot


#plot distribution of GC
dist_gc <- ggplot(transcriptome_data_expert_level, aes(x = GC, fill = taxon_info))+
  geom_histogram()+
  facet_wrap(~factor(taxon_id, levels = unique(taxon_id)))+#  facet_grid(~factor(my_variable, levels=c('val1', 'val2', 'val3', ...)))
  ggtitle('Distribution of GC')+
  theme(strip.text = element_text(size = 5))
dist_gc
ggsave('gc_dist.png', device = 'png', width = 8, height = 6)


#plot distribution of length
dist_len <- ggplot(transcriptome_data_expert_level, aes(x = reorder(taxon_id,length), y = length, color = taxon_info, fill= taxon_info))+
  geom_violin()+
  geom_boxplot(color = 'Black', fill = 'NA', outlier.shape = NA, width = 0.5)+
  ggtitle('Length Distribution')+
  scale_y_continuous(labels = comma, trans = 'log10')+
  ylab('Length (log10)')+
  xlab('Ten digit code, taxon info')
dist_len + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave('length_violin.png', device = 'png',width = 7.5, height = 6)


#plot distribution of coverage
dist_cov <- ggplot(transcriptome_data_expert_level, aes(x = reorder(taxon_id,cov), y = cov, color = taxon_info, fill= taxon_info))+
  geom_violin()+
  geom_boxplot(color = 'Black', fill = 'NA', outlier.shape = NA, width = 0.5)+
  ggtitle('Distribution of coverage')+
  scale_y_continuous(labels = comma, trans = 'log10')+
  ylab('Coverage (log10)')+
  xlab('Ten digit code, taxon info')
dist_cov+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave('cov_violin.png', device = 'png', width = 7.5, height = 6)
