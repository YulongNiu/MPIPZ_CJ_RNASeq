
## originally by Yulong Niu
## yulong.niu@hotmail.com

######################hierarchical clustering####################
setwd('/extDisk1/RESEARCH/MPIPZ_CJ_RNASeq/results_sig/')

library('readr')
library('magrittr')
library('tibble')
library('gplots')
library('dendextend')
library('dynamicTreeCut')
library('ggplot2')
library('tidyr')
library('DESeq2')
library('dplyr')
library('RColorBrewer')
library('gridExtra')
library('cluster')
library('scales')

load('eachGroup_mergeDay8.RData')
deganno <- read_csv('eachGroup_mergeDay8.csv',
                   col_types = cols(Chromosome = col_character()))

##~~~~~~~~~~~~~~~~~~~~~~useful funcs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
meanFe <- function(v) {

  require('magrittr')

  res <- v %>%
    split(rep(1 : 8, each = 6)) %>%
    sapply(mean, na.rm = TRUE)

  return(res)
}

##p value calculation from WGCNA
corPvalueStudent <- function(cor, nSamples) {

  ## ref: https://courses.lumenlearning.com/introstats1/chapter/testing-the-significance-of-the-correlation-coefficient/
  T <- sqrt(nSamples - 2) * cor / sqrt(1 - cor^2)

  p <- apply(T, 1:2, function(x) {
    if (x < 0) {
      eachp <- 1 -  pt(x, nSamples - 2, lower.tail = FALSE)
    } else {
      eachp <- pt(x, nSamples - 2, lower.tail = FALSE)
    }
  })

  return(p)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~prepare counts~~~~~~~~~~~~~~~~~~~~~~~~~~
rawCount <- rlog(degresRmZero) %>% assay

## mean value of normalized count
sampleN <- colnames(degresRmZero) %>% substring(first = 1, last = nchar(.) - 5) %>% unique
meanCount <- rawCount %>%
  apply(1, meanFe) %>%
  t
colnames(meanCount) <- sampleN

scaleCount <- meanCount %>%
  t %>%
  scale %>%
  t
scaleCount %<>% .[complete.cases(.), ]
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~K-means cluster~~~~~~~~~~~~~~~~~~~~~~
z_var <- apply(meanCount, 1, var)
z_mean <- apply(meanCount, 1, mean)
plot(z_mean, z_var, pch = '.')
abline(h = 1, col='red')
abline(v = 1, col='red')
text(x = 13,
     y = 23,
     labels = 'variance > 1 &\n mean > 1',
     col = 'red')

## filter
## meanCount %<>% .[which(z_var > 0 & z_mean > 0), ]

## choose groups
## 1. sum of squared error
wss <- (nrow(scaleCount) - 1) * sum(apply(scaleCount, 2, var))

for (i in 2:20) {
  wss[i] <- sum(kmeans(scaleCount,
                       centers=i,
                       algorithm = 'MacQueen')$withinss)
}

ggplot(tibble(k = 1:20, wss = wss), aes(k, wss)) +
  geom_point(colour = '#D55E00', size = 3) +
  geom_line(linetype = 'dashed') +
  xlab('Number of clusters') +
  ylab('Sum of squared error')
ggsave('kmeans_sse_mergeDay8.pdf')
ggsave('kmeans_sse_mergeDay8.jpg')

## 2. Akaike information criterion
kmeansAIC = function(fit){
  m = ncol(fit$centers)
  n = length(fit$cluster)
  k = nrow(fit$centers)
  D = fit$tot.withinss
  return(D + 2*m*k)
}

aic <- numeric(20)
for (i in 1:20) {
  fit <- kmeans(x = scaleCount, centers = i, algorithm = 'MacQueen')
  aic[i] <- kmeansAIC(fit)
}

ggplot(tibble(k = 1:20, aic = aic), aes(k, wss)) +
  geom_point(colour = '#009E73', size = 3) +
  geom_line(linetype = 'dashed') +
  xlab('Number of clusters') +
  ylab('Akaike information criterion')
ggsave('kmeans_AIC_mergeDay8.pdf')
ggsave('kmeans_AIC_mergeDay8.jpg')

## cluster
kClust10 <- kmeans(scaleCount, centers = 10, algorithm= 'MacQueen', nstart = 1000, iter.max = 20)

## selected genes
read_csv('../results/selected_genes.csv') %>%
  mutate(cl = ID %>% kClust10$cluster[.]) %>%
  filter(!is.na(cl)) %>%
  inner_join(deganno, c('ID' = 'ID')) %>%
  select(ID, Gene.x, Gene_Symbol, cl, Col0_FeCl3_Live_vs_Col0_FeCl3_HK_log2FoldChange) %>%
  ## filter(Col0_FeCl3_Live_vs_Col0_FeCl3_HK_log2FoldChange >= log2(1.5) | Col0_FeCl3_Live_vs_Col0_FeCl3_HK_log2FoldChange <= -log2(1.5)) %>%
  write_csv('selected_genes_cluster.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


##~~~~~~~~~~~~~~~~~~~~~~~~~~plot patterns~~~~~~~~~~~~~~~~~~~~~~~~
cl <- kClust10$cluster
prefix <- 'kmeans_10_mergeDay8'

clusterGene <- scaleCount %>%
  as.data.frame %>%
  rownames_to_column(var = 'ID') %>%
  as_tibble %>%
  {
    cl <- as.data.frame(cl) %>%
      rownames_to_column(var = 'ID')
    inner_join(., cl)
  }

## plot core cluster
clusterCore <- clusterGene %>%
  group_by(cl) %>%
  summarise_at(-1, mean, na.rm = TRUE) %>% ## mean of each cluster
  mutate(cl = paste0('cluster_', cl) %>%
           factor(levels = paste0('cluster_', cl))) %>%
  gather(Sample, NorExpress, -1) %>%
  mutate(Sample = Sample %>% factor(levels = sampleN, ordered = TRUE))

## ## choose colors
## show_col(hue_pal()(4))
ggplot(clusterCore, aes(Sample, NorExpress, col = cl, group = cl)) +
  geom_point() +
  geom_line() +
  facet_wrap(. ~ cl, ncol = 2) +
  ylab('Scaled counts') +
  scale_color_manual(values = hue_pal()(10),
    breaks = kClust10$cluster %>%
                       table %>%
                       names %>%
                       paste0('cluster_', .),
    labels = kClust10$cluster %>%
                       table %>%
                       {paste0('cluster_', names(.), ' ', .)},
    guide = guide_legend(title = 'kmeans (k = 10)')) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste0(prefix, '.pdf'))
ggsave(paste0(prefix, '.jpg'))

## plot all genes
clusterGenePlot <- clusterGene %>%
  gather(Sample, NorExpress, -ID, -cl) %>%
  mutate(cl = paste0('cluster_', cl) %>%
           factor(levels = paste0('cluster_', sort(unique(cl))))) %>%
  mutate(Sample = Sample %>% factor(levels = sampleN, ordered = TRUE))

clusterCorePlot <- clusterCore %>% dplyr::mutate(ID = 1 : nrow(clusterCore))
ggplot(clusterGenePlot, aes(Sample, NorExpress, group = ID)) +
  geom_line(color = 'grey30', alpha = 0.01) +
  facet_wrap(. ~ cl, ncol = 2) +
  geom_point(data = clusterCorePlot, aes(Sample, NorExpress, col = cl, group = ID)) +
  geom_line(data = clusterCorePlot, aes(Sample, NorExpress, group = cl, col = cl)) +
  ylab('Scaled counts') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  guides(colour = guide_legend(title = 'kmeans (k=10)'))
ggsave(paste0(prefix, '_genes.pdf'), width = 10, dpi = 320)
ggsave(paste0(prefix, '_genes.jpg'), width = 10, dpi = 320)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~cluster cor phenotype~~~~~~~~~~~~~~~~~
traits <- data.frame(FeEDTA = c(0, 0, 1, 1, 0, 0, 1, 1),
                     Live = c(0, 1, 0, 1, 0, 1, 0, 1),
                     Col0 = c(1, 1, 1, 1, 0, 0, 0, 0),
                     IronStarvation = c(1, 0, 0, 0, 1, 1, 0, 0))

cores <- clusterGene %>%
  group_by(cl) %>%
  summarise_at(-1, mean, na.rm = TRUE) %>%
  mutate(cl = cl %>% paste0('cluster_', .)) %>%
  column_to_rownames(var = 'cl') %>%
  t

moduleTraitCor <- cor(cores, traits, use = 'p')
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(traits))

traitPPlot <- moduleTraitPvalue %>%
  as.data.frame %>%
  rownames_to_column('cluster') %>%
  gather(trait, pvalue, -1) %>%
  as_tibble

traitCorPlot <- moduleTraitCor %>%
  as.data.frame %>%
  rownames_to_column('cluster') %>%
  gather(trait, correlation, -1) %>%
  as_tibble %>%
  mutate(x = rep(0 : (ncol(traits) - 1), each = ncol(cores))) %>%
  mutate(y = rep((ncol(cores) - 1) : 0, ncol(traits))) %>%
  inner_join(traitPPlot) %>%
  mutate(addtext = paste0(round(correlation, digit = 2),
                          '\n',
                          '(',
                          round(pvalue, digit = 2),
                          ')'))

ggplot(traitCorPlot, aes(x = x, y = y, fill = correlation)) +
  geom_tile() +
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name = 'RdYlBu')))(100),
                       breaks = seq(-1, 1, 0.5),
                       labels = format(seq(-1, 1, 0.5)),
                       limits = c(-1, 1)) +
  geom_text(aes(label = addtext)) +
  scale_x_continuous(breaks = 0 : (ncol(traits) - 1), labels = colnames(traits)) +
  scale_y_continuous(breaks = 0 : (ncol(cores) - 1), labels = paste0('cluster_', (ncol(cores)):1)) +
  xlab('Trait') +
  ylab('Cluster')
ggsave(paste0(prefix, '_trait.jpg'))
ggsave(paste0(prefix, '_trait.pdf'))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~heat map~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
scaleC <- rawCount %>%
  t %>%
  scale %>%
  t %>%
  as.data.frame %>%
  rownames_to_column('ID') %>%
  as_tibble %>%
  rename_at(-1, .funs = list(~paste0('Scale_', .)))

rawC <- rawCount %>%
  as.data.frame %>%
  rownames_to_column('ID') %>%
  as_tibble %>%
  rename_at(-1, .funs = list(~paste0('Raw_', .)))

degresC <- deganno %>%
  select(ID, Col0_FeCl3_HK_vs_Col0_FeEDTA_HK_pvalue : f6h1_FeCl3_Live_vs_f6h1_FeCl3_HK_log2FoldChange)

heatPlot <- rawC %>%
  inner_join(scaleC) %>%
  inner_join(degresC) %>%
  {
    cl <- as.data.frame(cl) %>%
      rownames_to_column(var = 'ID')
    inner_join(., cl)
  } %T>%
  {(sum(names(cl) == .$ID) == nrow(.)) %>% print} %>% ## check cl names and degresC row names
  slice(cl %>% order)

inner_join(deganno, heatPlot) %>%
  mutate_at(c('Gene', 'Description'), .funs = list(~if_else(is.na(.), '', .))) %>%
  write_csv(paste0(prefix, '.csv'))

heatRawPlot <- heatPlot %>%
  select(ID, starts_with('Raw')) %>%
  gather(sample, raw, -1) %>%
  mutate(x = rep(0 : 47, each = nrow(heatPlot))) %>%
  mutate(y = rep(0 : (nrow(heatPlot) - 1), 48))

heatScalePlot <- heatPlot %>%
  select(ID, starts_with('Scale')) %>%
  gather(sample, scale, -1) %>%
  mutate(x = rep(0 : 47, each = nrow(heatPlot))) %>%
  mutate(y = rep(0 : (nrow(heatPlot) - 1), 48))

heatlog2FCPlot <- heatPlot %>%
  select(ID, ends_with('FoldChange')) %>%
  gather(sample, log2FC, Flg22_vs_Mock_log2FoldChange : Flg22_SynCom35_vs_Mock_log2FoldChange) %>%
  mutate(x = rep(0 : 2, each = nrow(heatPlot))) %>%
  mutate(y = rep(0 : (nrow(heatPlot) - 1), 3))

## sig |FC| > 1 and padj < 0.05
fcsig <- heatPlot %>%
  select(ends_with('FoldChange')) %>%
  transmute_all(list(~ case_when(. > 1 ~ 1,
                                . < -1 ~ -1,
                                TRUE ~ 0)))

padjsig <- heatPlot %>%
  select(ends_with('padj')) %>%
  abs %>%
  `<`(0.05) %>%
  as_tibble %>%
  transmute_all(list(~ if_else(is.na(.), FALSE, .)))

heatsigPlot <- (padjsig * fcsig) %>%
  as_tibble %>%
  gather(sample, sig, 1:3) %>%
  mutate(x = rep(0 : 2, each = nrow(heatPlot))) %>%
  mutate(y = rep(0 : (nrow(heatPlot) - 1), 3))

heatGroupPlot <- heatPlot %>%
  select(ID, cluster = cl) %>%
  mutate(x = 0) %>%
  mutate(y = 0 : (nrow(heatPlot) - 1))

theme_flg22 <- function(...) {
  theme_bw() %+replace%
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks.length = unit(0, 'mm'),
          axis.line = element_blank(),
          panel.spacing = unit(0, 'mm'),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), 'line'),
          legend.spacing = unit(0, 'mm'),
          ...)
}

ggplot(heatRawPlot, aes(x = x, y = y, fill = log2(raw))) +
  geom_tile() +
  scale_fill_gradientn(colours = colorRampPalette(brewer.pal(n = 7, name = 'GnBu'))(100), name = 'log2(count)') +
  scale_x_continuous(breaks = 0 : 47,
                     labels = rep(sampleN, each = 6) %>%
                       paste(rep(1 : 6, 8), sep = '_')) +
  theme_flg22(legend.position = 'left',
              axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste0(prefix, '_heatmap_raw.jpg'))
ggsave(paste0(prefix, '_heatmap_raw.pdf'))

ggplot(heatScalePlot, aes(x = x, y = y, fill = scale)) +
  geom_tile() +
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 10, name = 'Spectral')))(100), name = 'scale(count)') +
  scale_x_continuous(breaks = 0 : 47,
                     labels = rep(sampleN, each = 6) %>%
                       paste(rep(1 : 6, 8), sep = '_')) +
  theme_flg22(legend.position = 'left',
              axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste0(prefix, '_heatmap_scale.jpg'))
ggsave(paste0(prefix, '_heatmap_scale.pdf'))

ggplot(heatlog2FCPlot, aes(x = x, y = y, fill = log2FC)) +
  geom_tile() +
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 8, name = 'PiYG')))(100), name = 'log2(FoldChange)') +
  scale_x_continuous(breaks = 0 : 2,
                     labels = paste(c('flg22', 'flg22_SynCom33', 'flg22_SynCom35'), 'vs. Mock')) +
  theme_flg22(legend.position = 'left',
              axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste0(prefix, '_heatmap_logFC.jpg'))
ggsave(paste0(prefix, '_heatmap_logFC.pdf'))

ggplot(heatsigPlot, aes(x = x, y = y)) +
  geom_tile(aes(fill = factor(sig))) +
  scale_fill_manual(name = 'Significant', labels = c('Down', 'No', 'Up'), values = c('green', 'grey90', 'red')) +
  scale_x_continuous(breaks = 0 : 2,
                     labels = paste(c('flg22', 'flg22_SynCom33', 'flg22_SynCom35'), 'vs. Mock')) +
  theme_flg22(legend.position = 'left',
              axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste0(prefix, '_heatmap_sig.jpg'))
ggsave(paste0(prefix, '_heatmap_sig.pdf'))

ggplot(heatGroupPlot, aes(x = x, y = y)) +
  geom_tile(aes(fill = factor(cluster))) +
  scale_x_continuous(breaks = 0,
                     labels = 'group') +
  theme_flg22(legend.position = 'left',
              axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste0(prefix, '_heatmap_group.jpg'))
ggsave(paste0(prefix, '_heatmap_group.pdf'))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~merge all plots~~~~~~~~~~~~~~~~~~~~
groupe <- ggplot(heatGroupPlot, aes(x = x, y = y)) +
  geom_tile(aes(fill = factor(cluster))) +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(expand = c(0, 0), breaks = NULL) +
  scale_x_continuous(expand = c(0, 0), breaks = NULL) +
  theme_flg22(title = element_blank(),
              legend.position = 'none')

groupne <- heatGroupPlot %>%
  group_by(cluster) %>%
  summarise(y = median(y)) %>%
  mutate(x = 0, cluster = cluster %>% paste0('cluster', .)) %>%
  ggplot(aes(x = x, y = y, label = cluster)) +
  geom_text() +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(expand = c(0, 0),  limits = c(0, nrow(heatGroupPlot)), breaks = NULL) +
  scale_x_continuous(expand = c(0, 0), breaks = NULL) +
  theme_flg22(title = element_blank(),
              legend.position = 'none')

rawe <- ggplot(heatRawPlot, aes(x = x, y = y, fill = log2(raw))) +
  geom_tile() +
  scale_fill_gradientn(colours = colorRampPalette(brewer.pal(n = 7, name = 'GnBu'))(100), name = 'log2(count)') +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(expand = c(0, 0), breaks = NULL) +
  scale_x_continuous(expand = c(0, 0), breaks = NULL) +
  theme_flg22(title = element_blank(),
              legend.position = 'none')

scalee <- ggplot(heatScalePlot, aes(x = x, y = y, fill = scale)) +
  geom_tile() +
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 10, name = 'Spectral')))(10), name = 'scale(count)') +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(expand = c(0, 0), breaks = NULL) +
  scale_x_continuous(expand = c(0, 0), breaks = NULL) +
  theme_flg22(title = element_blank(),
              legend.position = 'none')

fce <- ggplot(heatlog2FCPlot, aes(x = x, y = y, fill = log2FC)) +
  geom_tile() +
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 8, name = 'PiYG')))(100), name = 'log2(FoldChange)') +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(expand = c(0, 0), breaks = NULL) +
  scale_x_continuous(expand = c(0, 0), breaks = NULL) +
  theme_flg22(title = element_blank(),
              legend.position = 'none')

sige <- ggplot(heatsigPlot, aes(x = x, y = y)) +
  geom_tile(aes(fill = factor(sig))) +
  scale_fill_manual(name = 'Significant', labels = c('Down', 'No', 'Up'), values = c('green', 'grey90', 'red')) +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(expand = c(0, 0), breaks = NULL) +
  scale_x_continuous(expand = c(0, 0), breaks = NULL) +
  theme_flg22(title = element_blank(),
              legend.position = 'none')

sigte <- heatGroupPlot %>%
  select(cluster, y) %>%
  inner_join(heatsigPlot) %>%
  select(sample, sig, x, y, cluster) %>%
  group_by(sample, cluster) %>%
  count(sig) %>%
  spread(sig, n) %>%
  {
    loc <- heatGroupPlot %>%
      group_by(cluster) %>%
      summarise(y = median(y))
    inner_join(., loc)
  } %>%
  rename('down' = `-1`, 'no' = `0`, 'up' = `1`) %>%
  ungroup %>%
  mutate_at(c('down', 'no', 'up'), .funs = list(~if_else(is.na(.), 0L, .))) %>%
  mutate(x = rep(c(0.2, 0.4, 0), each = max(cl))) %>%
  mutate(signum = paste0(down, '/', no, '/', up)) %>%
  select(signum, x, y) %>%
  ggplot(aes(x = x, y = y, label = signum)) +
  geom_text() +
  labs(x = NULL, y = NULL) +
  scale_y_continuous(expand = c(0, 0),  limits = c(0, nrow(heatGroupPlot)), breaks = NULL) +
  scale_x_continuous(expand = c(0, 0), limits = c(-0.1, 0.5), breaks = NULL) +
  theme_flg22(title = element_blank(),
              legend.position = 'none')


blanke <- ggplot(tibble(x = 0, y = 0 : (nrow(heatPlot) - 1)),
                 aes(x = x, y = y)) +
  geom_tile(colour = 'white') +
  scale_y_continuous(expand = c(0, 0), breaks = NULL) +
  scale_x_continuous(expand = c(0, 0), breaks = NULL) +
  theme_flg22(title = element_blank(),
              legend.position = 'none')

g <- grid.arrange(groupne,
                  groupe,
                  blanke,
                  rawe,
                  blanke,
                  scalee,
                  blanke,
                  fce,
                  blanke,
                  sige,
                  nrow = 1,
                  ncol = 10,
                  widths = c(3.5, 1, 0.5, 13, 0.5, 13, 0.5, 3, 0.5, 3) %>% {. / sum(.)})
ggsave(file = paste0(prefix, '_heatmap_merge_1stadd.pdf'), plot = g)
ggsave(file = paste0(prefix, '_heatmap_merge_1stadd.jpg'), plot = g)

g <- grid.arrange(groupne,
                  groupe,
                  blanke,
                  scalee,
                  nrow = 1,
                  ncol = 4,
                  widths = c(3.5, 1, 0.5, 13) %>% {. / sum(.)})
ggsave(file = paste0(prefix, '_heatmap_all.pdf'), plot = g)
ggsave(file = paste0(prefix, '_heatmap_all.jpg'), plot = g)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## write the cluster file
inner_join(deganno, heatPlot) %>%
  mutate_at(c('Gene', 'Description'), .funs = list(~if_else(is.na(.), '', .))) %>%
  write_csv(paste0(prefix, '.csv'))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################################################################

##########################plot genes############################
setwd('/extDisk1/RESEARCH/MPIPZ_CJ_RNASeq/results/')

library('ggplot2')
library('readr')
library('dplyr')
library('magrittr')
library('tidyr')

anno <- read_csv('eachGroup_vs_Col0_FeCl3_HK_Day8_raw_k.csv',
                 col_types = cols(Chromosome = col_character())) %>%
  mutate(Gene = Gene %>% {if_else(is.na(.), '', .)}) %>%
  mutate(Description = Description %>% {if_else(is.na(.), '', .)})

cgenes <- read_csv('selected_genes.csv',
                   col_types = cols(Chromosome = col_character()))

genePlot <- anno %>%
  inner_join(cgenes) %>%
  mutate(Gene = paste(Gene_Symbol, ID, sep = '_')) %>%
  select(Gene, Col0_FeCl3_HK_Day8_Rep1 : f6h1_FeEDTA_Live_Day15_Rep3) %>%
  gather(ID, NormCount, -1) %>%
  mutate(ID = ID %>%
           substring(1, nchar(.) - 5)) %>%
  mutate(Time = ID %>%
           strsplit(split = '_', fixed = TRUE) %>%
           sapply('[[', 4) %>%
           factor(levels = c('Day8', 'Day15'))) %>%
  mutate(Condi = ID %>%
           strsplit(split = '_', fixed = TRUE) %>%
           sapply(function(x) {paste(x[1:3], collapse = '_')})) %>%
  mutate(ID = factor(ID)) %>%
  mutate(Gene = factor(Gene))

genes <- genePlot$Gene %>% unique %>% as.character

for(i in genes) {

  genePlot %>%
    filter(Gene %in% i) %>%
    ggplot(aes(x = Condi, y = NormCount)) +
    geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.7) +
    stat_summary(fun.y = mean, geom = 'point', color='red', size = 3) +
    facet_wrap(. ~ Time, ncol = 2) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

  ggsave(file = paste0('selectgenes/', i, '.pdf'))
  ggsave(file = paste0('selectgenes/', i, '.jpg'))
}
##################################################################

########################separate DEGs############################
library('tibble')
library('readr')
library('dplyr')

setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/')
savepath <- '/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/kmeanssig'

## cgenes <- c('AT1G14550.1', 'AT2G30750.1', 'AT2G19190.1')
## kres %>%
##   slice(which(kres$ID %in% cgenes))

kres <- read_csv('kmeans_10.csv',
                 col_types = cols(Chromosome = col_character()))

## base columns
g <- c('Flg22', 'Flg22_SynCom33', 'Flg22_SynCom35')

for (i in g) {

  eachcols <- c(paste(i, 1:3, sep = '_'),
                paste0(i, c('_vs_Mock_pvalue', '_vs_Mock_padj', '_vs_Mock_log2FoldChange')))

  eachres <- kres %>%
    select(ID : Mock_3, eachcols, cl) %>%
    filter(eachcols[6] %>% get %>% abs %>% {. > 1}) %>%
    arrange(eachcols[6] %>% get %>% desc)

  cls <- eachres$cl %>% unique

  for (j in cls) {
    fname <- i %>% paste0('_vs_Mock_cluster', j, '.csv') %>% file.path(savepath, .)
    eachres %>%
      filter(cl == j) %>%
      mutate(Gene = Gene %>% coalesce('')) %>%
      mutate(Description = Description %>% coalesce('')) %>%
      write_csv(fname)
  }
}
#################################################################
