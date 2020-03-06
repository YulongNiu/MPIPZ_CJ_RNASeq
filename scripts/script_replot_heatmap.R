########################replot heatmap with DEGs and cluster############
library('tidyverse')
library('ComplexHeatmap')
library('limma')
library('DESeq2')
library('RColorBrewer')

setwd('/extDisk1/RESEARCH/MPIPZ_CJ_RNASeq/results/')

load('eachGroup_mergeDay8.RData')

##~~~~~~~~~~~~~~~~~~~~~~useful funcs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
meanFe <- function(v) {

  require('magrittr')

  res <- v %>%
    split(rep(1 : 8, each = 6)) %>%
    sapply(mean, na.rm = TRUE)

  return(res)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~select DEGs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
wholeDEG <- read_csv('eachGroup_mergeDay8.csv')
kmeansRes <- read_csv('kmeans_10_mergeDay8.csv') %>%
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
  dplyr::slice(wholeDEG, .) %>%
  inner_join(kmeansRes)

heatsig %>%
  mutate_at(c('Gene', 'Description'), .funs = list(~if_else(is.na(.), '', .))) %>%
  write_csv('eachGroup_mergeDay8_deg_sig.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~heatmap~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## rlog transformed
rawC <- rldData %>%
  as.data.frame %>%
  .[, c(13:24, 1:12, 37:48, 25:36)] %>%
  rownames_to_column('ID') %>%
  as_tibble %>%
  inner_join(heatsig %>% select(ID, cl))

scaleC <- rawC %>%
  select(contains('Rep')) %>%
  t %>%
  scale %>%
  t %>%
  as_tibble %>%
  bind_cols(rawC %>% select(ID, cl))

cairo_pdf('kmeans_10_mergeDay8_heatmap_sig2.pdf')
iron <- HeatmapAnnotation(Iron = rep(rep(c('EDTA', 'FeCl3'), 2),
                                     each = 12),
                          col = list(Iron = c('EDTA' = 'white', 'FeCl3' = 'grey')),
                          gp = gpar(col = 'black'))
syncom <- HeatmapAnnotation(SynCom = rep(rep(c('HK', 'Live'), 4), each = 6),
                            col = list(SynCom = c('HK' = 'white', 'Live' = 'grey')),
                            gp = gpar(col = 'black'))
Heatmap(matrix = scaleC %>% select(contains('Rep')),
        name = 'Scaled Counts',
        row_order = order(scaleC$cl) %>% rev,
        row_split = scaleC$cl,
        row_gap = unit(2, "mm"),
        column_order = 1 : 48,
        column_split = rep(c('Col-0', 'f6\'h1'), each = 24),
        show_column_names = FALSE,
        col = colorRampPalette(rev(brewer.pal(n = 10, name = 'Spectral'))[c(-3, -4, -7, -8)])(10),
        top_annotation = c(iron, syncom))
dev.off()
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~box plot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sampleN <- c('Col0_FeEDTA_HK', 'Col0_FeEDTA_Live',
             'Col0_FeCl3_HK', 'Col0_FeCl3_Live',
             'f6h1_FeEDTA_HK', 'f6h1_FeEDTA_Live',
             'f6h1_FeCl3_HK', 'f6h1_FeCl3_Live')

boxplotData <- rldData %>%
  .[, c(13:24, 1:12, 37:48, 25:36)] %>%
  t %>%
  scale %>%
  t %>%
  apply(1, meanFe) %>%
  t %>%
  set_colnames(sampleN) %>%
  as.data.frame %>%
  rownames_to_column('ID') %>%
  as_tibble %>%
  inner_join(kmeansRes)

for (i in 1:10) {
  boxplotData %>%
    filter(cl == i) %>%
    select(-ID, -cl) %>%
    gather(key = 'Conditions', value = 'ScaleCounts') %>%
    mutate(Conditions = Conditions %>% factor(levels = sampleN)) %>%
    ggplot(aes(x = Conditions, y = ScaleCounts)) +
    geom_boxplot() +
    ylab('Scaled counts') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90),
          plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
          legend.text.align = 0,
          axis.text = element_text(size = 13),
          axis.title = element_text(size = 14),
          legend.text=element_text(size= 13),
          legend.title = element_text(size = 14))

  ggsave(paste0('kmeans_10_mergeDay8_boxplot_cluster', i, '.pdf'))
  ggsave(paste0('kmeans_10_mergeDay8_boxplot_cluster', i, '.jpeg'))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#######################################################################
