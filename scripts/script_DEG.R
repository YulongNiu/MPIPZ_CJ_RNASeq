###########################DEGs##################################
##~~~~~~~~~~~~~~~~~~~~~~~~~~~useful funcs~~~~~~~~~~~~~~~~~~~~~~~~
checkZeros <- function(v, threshold) {
  res <- ifelse(sum(v == 0) > threshold, FALSE, TRUE)
  return(res)
}

checkFe <- function(v, threshold) {

  require('magrittr')

  res <- v %>%
    split(rep(1 : 8, each = 3)) %>%
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
  mutate(Days = rep(c(8, 15, 8, 14), each = 24)) %>%
  filter(Sample %>% str_detect('3989')) %>%
  arrange(Anno) %>%
  arrange(Days) %>%
  filter(Timepoint == 'Day15')

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
sampleTable <- data.frame(condition = factor(rep(condi, each = 3), levels = condi))
rownames(sampleTable) <- colnames(kres$counts)

degres <- DESeqDataSetFromTximport(kres, sampleTable, ~ condition)

## DEGs
degres %<>%
  estimateSizeFactors %>%
  counts(normalized = TRUE) %>%
  apply(1, checkFe, 1) %>%
  degres[., ]

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
mod <- model.matrix(~ condition, colData(degres))
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
           .$condition %>%
           rep(svnum) %>%
           as.character,
         sample = rep(colnames(degres), svnum)) %>%
  mutate(group = paste(sv, condition, sep = '_')) %>%
  ggplot(aes(sample, value, colour = sv, group = group)) +
  geom_point() +
  geom_line() +
  theme(axis.text.x = element_text(angle = 90))
ggsave('auto_sv_day15.jpg')
ggsave('auto_sv_day15.pdf')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~pair-wise comparison~~~~~~~~~~~~~~~~~
degres$sv1 <- svobj$sv[, 1]
degres$sv2 <- svobj$sv[, 2]
degres$sv3 <- svobj$sv[, 3]
degres$sv4 <- svobj$sv[, 4]
design(degres) <- ~sv1 + sv2 + sv3 + sv4 + condition

cond <- list(c('Col0_FeCl3_HK_Day15', 'Col0_FeEDTA_HK_Day15'),
             c('Col0_FeCl3_Live_Day15', 'Col0_FeEDTA_Live_Day15'),
             c('f6h1_FeCl3_HK_Day15', 'f6h1_FeEDTA_HK_Day15'),
             c('f6h1_FeCl3_Live_Day15', 'f6h1_FeEDTA_Live_Day15'),
             c('Col0_FeEDTA_Live_Day15', 'Col0_FeEDTA_HK_Day15'),
             c('Col0_FeCl3_Live_Day15', 'Col0_FeCl3_HK_Day15'),
             c('f6h1_FeEDTA_Live_Day15', 'f6h1_FeEDTA_HK_Day15'),
             c('f6h1_FeCl3_Live_Day15', 'f6h1_FeCl3_HK_Day15'))

resRaw <- lapply(cond,
                 function(x) {
                   degres %>%
                     results(contrast = c('condition', x)) %T>%
                     summary %>%
                     as_tibble %>%
                     select(pvalue, padj, log2FoldChange) %>%
                     rename_all(.funs = list(~paste0(paste(x, collapse = '_vs_'), '_', .)))
                 }) %>%
  bind_cols

res <- cbind.data.frame(as.matrix(mcols(degres)[, 1:10]), assay(ntd), stringsAsFactors = FALSE) %>%
  rownames_to_column(., var = 'ID') %>%
  as_tibble %>%
  bind_cols(resRaw) %>%
  inner_join(anno, by = 'ID') %>%
  select(ID, Gene : Description, Col0_FeCl3_HK_Day15_Rep1 : f6h1_FeCl3_Live_Day15_vs_f6h1_FeCl3_HK_Day15_log2FoldChange) %>%
  arrange(Col0_FeCl3_HK_Day15_vs_Col0_FeEDTA_HK_Day15_padj)

write_csv(res, 'eachGroup_Day15.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PCA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library('directlabels')
library('ggplot2')
library('RColorBrewer')
library('limma')
library('sva')

rldDay15 <- colData(rld)[, 1] %>%
  as.character %>%
  strsplit(split = '_', fixed = TRUE) %>%
  do.call(rbind, .) %>%
  set_colnames(c('Genotype', 'Iron', 'SynCom', 'Time')) %>%
  as_tibble %>%
  mutate(Conditions = paste(Iron, SynCom, sep = '_') %>%
           factor,
         Genotype = Genotype %<>% factor)

cols <- brewer.pal(4, name = 'Set1')
cols[3:4] <- cols[4:3]

group <- sampleTable$condition
design <- model.matrix(~ group)
rldData <- dat %>%
  removeBatchEffect(covariates = svobj$sv,
                    design = design)

pca <- prcomp(t(rldData))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
percentVar <- round(100 * percentVar)
pca1 <- pca$x[,1]
pca2 <- pca$x[,2]

pcaData <- rldDay15 %>%
  mutate(Colours = factor(Conditions, labels = cols)) %>%
  select(Conditions, Genotype, Colours) %>%
  mutate(PC1 = pca1, PC2 = pca2)

ggplot(pcaData, aes(x = PC1, y = PC2, colour = Conditions)) +
  geom_point(aes(shape = Genotype), size = 4) +
  scale_colour_manual(values = levels(pcaData$Colours),
                      name = 'Experimental\nCondition',
                      labels = expression(FeCl[3]+HK, FeCl[3]+Live, FeEDTA+HK, FeEDTA+Live)) +
  scale_shape_manual(values = c(15, 17)) +
  xlab(paste0('PC1: ',percentVar[1],'% variance')) +
  ylab(paste0('PC2: ',percentVar[2],'% variance')) +
  ggtitle('Day15') +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'))

ggsave('../results/PCA_Day15_sva.pdf', width = 12)
ggsave('../results/PCA_Day15_sva.jpg', width = 12)

write_csv(res, 'eachGroup_Day15.csv')
save(degres, rldData, file = 'eachGroup_Day15.RData')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################
