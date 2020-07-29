
## originally by Yulong Niu
## yulong.niu@hotmail.com

###########################DEGs##################################
##~~~~~~~~~~~~~~~~~~~~~~~~~~~useful funcs~~~~~~~~~~~~~~~~~~~~~~~~
checkZeros <- function(v, threshold) {
  res <- ifelse(sum(v == 0) > threshold, FALSE, TRUE)
  return(res)
}

checkFe <- function(v, threshold) {

  require('magrittr')

  res <- v %>%
    split(rep(1 : 16, each = 3)) %>%
    sapply(checkZeros, threshold) %>%
    all

  return(res)
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library('tximport')
library('rhdf5')
library('magrittr')
library('DESeq2')
library('tidyverse')
library('foreach')

anno <- read_csv('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/Ensembl_ath_Anno.csv',
                 col_types = cols(Chromosome = col_character())) %>%
  mutate(Gene = Gene %>% {if_else(is.na(.), '', .)}) %>%
  mutate(Description = Description %>% {if_else(is.na(.), '', .)})


##~~~~~~~~~~~~~~~~~~~~load k alignments~~~~~~~~~~~~~~~~~~~~~~~~~~
wd <- '/extDisk1/RESEARCH/MPIPZ_CJ_RNASeq/bamfiles'
setwd(wd)

labelanno <- read_csv('../results/LibraryIDs.csv') %>%
  mutate_all(list(~ str_replace(., 'Day14', 'Day15'))) %>%
  mutate(Days = rep(c(8, 15, 8, 15), each = 24)) %>%
  filter(Timepoint %>% str_detect('Day8')) %>%
  arrange(Anno) %>%
  arrange(Days)

slabel <- labelanno$Sample %>%
  paste0('_ath_kallisto')

files <- file.path(wd, slabel, 'abundance.h5')
names(files) <- labelanno$Anno

kres <- tximport(files, type = 'kallisto', txOut = TRUE)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~normalization~~~~~~~~~~~~~~~~~
setwd('/extDisk1/RESEARCH/MPIPZ_CJ_RNASeq/results/')

## sampleTable
condi <- labelanno$Anno %>% substring(first = 1, last = nchar(.) - 5) %>% unique

sampleTable <- factor(rep(condi, each = 6), levels = condi) %>%
  as.character %>%
  strsplit(split = '_', fixed = TRUE) %>%
  do.call(rbind, .) %>%
  set_colnames(c('Genotype', 'Iron', 'SynCom', 'Time')) %>%
  as_tibble %>%
  mutate(Conditions = paste(Genotype, Iron, SynCom, sep = '_') %>%
           factor,
         Treatment = paste(Iron, SynCom, sep = '_') %>%
           factor,
         Genotype = Genotype %<>% factor,
         SynCom = SynCom %<>% factor,
         Iron = Iron %<>% factor,
         Batch = c(rep(rep(1:2, each = 3), 8))) %>%
  as.data.frame %>%
  set_rownames(colnames(kres$counts))

degres <- DESeqDataSetFromTximport(kres, sampleTable, ~ Conditions)

degres <- DESeq(degres)
## resultsNames(degres)

## count transformation
rld <- rlog(degres)
ntd <- normTransform(degres)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~hidden batch effect~~~~~~~~~~~~~~~~~~~~~
library('sva')
library('ggplot2')

dat <- rld %>%
  assay %>%
  {.[rowMeans(.) > 1, ]}
mod <- model.matrix(~ Conditions, colData(degres))
mod0 <- model.matrix(~ 1, colData(degres))

## manual detect surrogate variance
svnum <- 4
svseq <- svaseq(dat, mod, mod0, n.sv = svnum)

## auto detect sv
svobj <- sva(dat, mod, mod0)
svnum <- svobj$sv %>% ncol

svobj$sv %>%
  set_colnames(paste0('sv', seq_len(svnum))) %>%
  as_tibble %>%
  gather(key = 'sv', value = 'value') %>%
  mutate(condition = colData(degres) %>%
           .$Conditions %>%
           rep(svnum) %>%
           as.character,
         sample = rep(colnames(degres), svnum)) %>%
  mutate(group = paste(sv, condition, sep = '_')) %>%
  ggplot(aes(sample, value, colour = sv, group = group)) +
  geom_point() +
  geom_line() +
  theme(axis.text.x = element_text(angle = 90))
ggsave('auto_sv_day8.jpg')
ggsave('auto_sv_day8.pdf')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~pair-wise comparison~~~~~~~~~~~~~~~~~~~~~~~~~
degres$sv1 <- svobj$sv[, 1]
degres$sv2 <- svobj$sv[, 2]
degres$sv3 <- svobj$sv[, 3]
degres$sv4 <- svobj$sv[, 4]
design(degres) <- ~sv1 + sv2 + sv3 + sv4 + Conditions

degres <- DESeq(degres)

cond <- list(c('Col0_FeCl3_HK', 'Col0_FeEDTA_HK'),
             c('Col0_FeCl3_Live', 'Col0_FeEDTA_Live'),
             c('f6h1_FeCl3_HK', 'f6h1_FeEDTA_HK'),
             c('f6h1_FeCl3_Live', 'f6h1_FeEDTA_Live'),
             c('Col0_FeEDTA_Live', 'Col0_FeEDTA_HK'),
             c('Col0_FeCl3_Live', 'Col0_FeCl3_HK'),
             c('f6h1_FeEDTA_Live', 'f6h1_FeEDTA_HK'),
             c('f6h1_FeCl3_Live', 'f6h1_FeCl3_HK'))

resRaw <- lapply(cond,
                 function(x) {
                   degres %>%
                     results(contrast = c('Conditions', x)) %T>%
                     summary %>%
                     as_tibble %>%
                     select(pvalue, padj, log2FoldChange) %>%
                     rename_all(.funs =  list(~paste0(paste(x, collapse = '_vs_'), '_', .)))
                 }) %>%
  bind_cols

res <- cbind.data.frame(as.matrix(mcols(degres)[, 1:10]), assay(rld), stringsAsFactors = FALSE) %>%
  rownames_to_column(., var = 'ID') %>%
  as_tibble %>%
  bind_cols(resRaw) %>%
  inner_join(anno, by = 'ID') %>%
  select(ID, Gene : Description, Col0_FeCl3_HK_Day8_Rep1 : f6h1_FeCl3_Live_vs_f6h1_FeCl3_HK_log2FoldChange) %>%
  arrange(Col0_FeCl3_HK_vs_Col0_FeEDTA_HK_padj)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PCA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library('directlabels')
library('ggplot2')
library('RColorBrewer')
library('limma')
library('sva')

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Day8~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cols <- c('#1b9e77', '#d95f02')

group <- sampleTable$Conditions
design <- model.matrix(~ group)
rldData <- dat %>%
  removeBatchEffect(covariates = svobj$sv,
                    design = design)

pca <- prcomp(t(rldData))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
percentVar <- round(100 * percentVar)
pca1 <- pca$x[,1]
pca2 <- pca$x[,2]

pcaData <- colData(degres) %>%
  as_tibble %>%
  mutate(Colours = factor(Genotype, labels = cols)) %>%
  mutate(Exp = paste(Iron, SynCom, sep = '_') %>% factor) %>%
  select(Conditions, Genotype, Colours, Exp, Batch) %>%
  mutate(PC1 = pca1, PC2 = pca2)

ggplot(pcaData, aes(x = PC1, y = PC2, colour = Colours)) +
  geom_point(aes(shape = Exp), size = 4) +
  scale_colour_manual(values = levels(pcaData$Colours),
                      name = 'Genotype',
                      labels = expression('Col-0',  italic("f6\'h1"))) +
  scale_shape_manual(values = c(2, 17, 1, 16),
                     name = 'Experimental\nConditions',
                     labels = expression(FeCl[3]+HK, FeCl[3]+Live, FeEDTA+HK, FeEDTA+Live)) +
  stat_ellipse(aes(x = PC1, y = PC2, group = Conditions), type = 't', linetype = 2) +
  coord_fixed(1) +
  xlab(paste0('PC1: ',percentVar[1],'% variance')) +
  ylab(paste0('PC2: ',percentVar[2],'% variance')) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
        legend.text.align = 0,
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 14),
        legend.text=element_text(size= 13),
        legend.title = element_text(size = 14))

ggsave('PCA_mergeDay8_sva.pdf', width = 12)
ggsave('PCA_mergeDay8_sva.jpg', width = 12)

write_csv(res, 'eachGroup_mergeDay8.csv')
save(degres, rldData, file = 'eachGroup_mergeDay8.RData')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################
