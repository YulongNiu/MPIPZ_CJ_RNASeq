
## originally by Yulong Niu
## yulong.niu@hotmail.com

##############################tile plot###############################
library('tidyverse')
library('foreach')
library('RColorBrewer')

setwd('/extDisk1/RESEARCH/MPIPZ_CJ_RNASeq/geneset')

intereGO <- read_csv('Interesting_GO_terms.csv')

GOfiles <- dir(pattern = 'mergeDay8.*BP.csv')
intereGOMerge <- foreach(i = seq_along(GOfiles), .combine = bind_rows) %do% {

  eachGO <- read_csv(GOfiles[i]) %>%
    select(-X1, -under_represented_pvalue, -ontology) %>%
    mutate(pvalue = -log2(over_represented_pvalue)) %>%
    mutate(cl = str_extract(GOfiles[i], '(?<=cluster).*?(?=_BP)') %>%
             as.numeric) %>%
    left_join(intereGO, .)
}

intereGOMerge %>%
  mutate(pvalue = pvalue %>% if_else(. > -log2(0.05), ., 0)) %>%
  mutate(cl = cl %>%
           paste0('cluster_', .) %>%
           factor(levels = paste0('cluster_', sort(unique(cl))))) %>%
  ggplot(aes(x = Annotation, y = cl, fill = pvalue)) +
  geom_raster() +
  scale_fill_gradientn(colours = colorRampPalette(brewer.pal(n = 8, name = 'RdYlBu') %>% rev %>% {c('#FFFFFF', .[6:8])})(20), name = '-log2(p-value)') +
  xlab('Gene ontology - Biological process') +
  ylab('cluster') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave('Interesting_GO_terms.jpg', width = 10)
ggsave('Interesting_GO_terms.pdf', width = 10)
######################################################################
