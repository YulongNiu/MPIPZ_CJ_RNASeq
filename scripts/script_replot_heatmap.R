########################replot heatmap with DEGs and cluster############
library('tidyverse')
library('ComplexHeatmap')

setwd('/extDisk1/RESEARCH/MPIPZ_CJ_RNASeq/results/')

wholeDEG <- read_csv('eachGroup_mergeDay8.csv')
kmeansRes <- read_csv('../results_sig/kmeans_10_mergeDay8.csv') %>%
  select(ID, cl)

fcsig <- wholeDEG %>%
  select(ends_with('FoldChange')) %>%
  transmute_all(list(~ case_when(. > log2(1.5) ~ 1,
                                 . < -log2(1.5) ~ 1,
                                 TRUE ~ 0)))
padjsig <- wholeDEG %>%
  select(ends_with('padj')) %>%
  abs %>%
  `<`(0.05) %>%
  as_tibble %>%
  transmute_all(list(~ if_else(is.na(.), FALSE, .)))

heatsig <- (padjsig * fcsig) %>%
  as_tibble %>%
  rowSums %>%
  {. >= 1} %>%
  which %>%
  slice(wholeDEG, .) %>%
  inner_join(kmeansRes)

write_csv(heatsig, 'eachGroup_mergeDay8_deg_sig.csv')
#######################################################################