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
library('tibble')
library('readr')
library('dplyr')
library('stringr')

anno <- read_csv('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/Ensembl_ath_Anno.csv',
                 col_types = cols(Chromosome = col_character())) %>%
  mutate(Gene = Gene %>% {if_else(is.na(.), '', .)}) %>%
  mutate(Description = Description %>% {if_else(is.na(.), '', .)})


##~~~~~~~~~~~~~~~~~~~~load k alignments~~~~~~~~~~~~~~~~~~~~~~~~~~
wd <- '/extDisk1/RESEARCH/MPIPZ_CJ_RNASeq/bamfiles'
setwd(wd)

labelanno <- read_csv('../results/LibraryIDs.csv') %>%
  mutate(Days = rep(c(8, 15, 8, 14), each = 24)) %>%
  filter(Sample %>% str_detect('4206')) %>%
  arrange(Anno) %>%
  arrange(Days)

slabel <- labelanno$Sample %>%
  paste0('_ath_kallisto')

files <- file.path(wd, slabel, 'abundance.h5')
names(files) <- labelanno$Anno

kres <- tximport(files, type = 'kallisto', txOut = TRUE)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~DEG Mock vs. 3 conditions~~~~~~~~~~~~~~~~~
setwd('/extDisk1/RESEARCH/MPIPZ_CJ_RNASeq/results/')

## sampleTable
condi <- labelanno$Anno %>% substring(first = 1, last = nchar(.) - 5) %>% unique
sampleTable <- data.frame(condition = factor(rep(condi, each = 3), levels = condi))
rownames(sampleTable) <- colnames(kres$counts)

degres <- DESeqDataSetFromTximport(kres, sampleTable, ~condition)

## ## Col0_FeCl3_HK_Day8
## degresCol0Day15 <- str_detect(labelanno$Anno, 'Col0.*Day8') %>%
##   degres[, .]

## DEGs
degres %<>%
  estimateSizeFactors %>%
  counts(normalized = TRUE) %>%
  apply(1, checkFe, 1) %>%
  degres[., ]
## degres <- degres[rowSums(counts(degres)) > 1, ]
save(degres, file = 'degres_1stadd.RData')

degres <- DESeq(degres)
## resultsNames(degres)

## count transformation
rld <- rlog(degres)
vst <- varianceStabilizingTransformation(degres)
ntd <- normTransform(degres)

cond <- degres %>%
  resultsNames %>%
  str_extract('(?<=condition_).*') %>%
  .[!is.na(.)]

resRaw <- lapply(cond,
                 function(x) {
                   degres %>%
                     results(name = paste0('condition_', x)) %T>%
                     summary %>%
                     as_tibble %>%
                     select(pvalue, padj, log2FoldChange) %>%
                     rename_all(.funs = list(~paste0(x, '_', .)))
                 }) %>%
  bind_cols


res <- cbind.data.frame(as.matrix(mcols(degres)[, 1:10]), assay(ntd), stringsAsFactors = FALSE) %>%
  rownames_to_column(., var = 'ID') %>%
  as_tibble %>%
  bind_cols(resRaw) %>%
  inner_join(anno, by = 'ID') %>%
  select(ID, Gene : Description, Mock_1 : Flg22_SynCom35_vs_Mock_log2FoldChange) %>%
  arrange(Flg22_vs_Mock_padj)

write_csv(res, 'eachGroup_vs_Mock_k.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


##~~~~~~~~~~~~~~~~~~~~~~~~DEG flg22_SynCom vs flg22~~~~~~~~~~~~~~~~~
setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/')

## sampleTable
sampleTable <- data.frame(condition = factor(rep(c('Mock', 'Flg22', 'Flg22_SynCom33', 'Flg22_SynCom35'), each = 3)), process = factor(rep(c('N', 'Y', 'Y', 'Y'), each = 3)))
sampleTable$condition %<>% relevel(ref = 'Flg22')
rownames(sampleTable) <- colnames(kres$counts)

degres <- DESeqDataSetFromTximport(kres,
                                   sampleTable,
                                   ~condition)

## DEGs role out two zeros in one group
## degres <- degres[rowSums(counts(degres)) > 1, ]
degres %<>%
  estimateSizeFactors %>%
  counts(normalized = TRUE) %>%
  apply(1, checkFlg22, 1) %>%
  degres[., ]
degres <- DESeq(degres)

## count transformation
rld <- rlog(degres)
vst <- varianceStabilizingTransformation(degres)
ntd <- normTransform(degres)

cond <- degres %>%
  resultsNames %>%
  str_extract('(?<=condition_).*') %>%
  .[!is.na(.)] %>%
  `[`(1:2)

resRaw <- lapply(cond,
                 function(x) {
                   degres %>%
                     results(name = paste0('condition_', x)) %T>%
                     summary %>%
                     as_tibble %>%
                     select(pvalue, padj, log2FoldChange) %>%
                     rename_all(.funs = list(~paste0(x, '_', .)))
                 }) %>%
  bind_cols

res <- cbind.data.frame(as.matrix(mcols(degres)[, 1:10]), assay(ntd), stringsAsFactors = FALSE) %>%
  rownames_to_column(., var = 'ID') %>%
  as_tibble %>%
  bind_cols(resRaw) %>%
  inner_join(anno, by = 'ID') %>%
  select(ID, Gene : Description, Flg22_1 : Flg22_SynCom35_vs_Flg22_log2FoldChange) %>%
  arrange(Flg22_SynCom33_vs_Flg22_padj)

write_csv(res, 'SynCom_vs_flg22_k.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~SynCom35 vs. SynCom35~~~~~~~~~~~~~~~~~~
setwd('/extDisk1/RESEARCH/MPIPZ_KaWai_RNASeq/results/')

## sampleTable
sampleTable <- data.frame(condition = factor(rep(c('Mock', 'Flg22', 'Flg22_SynCom33', 'Flg22_SynCom35'), each = 3)), process = factor(rep(c('N', 'Y', 'Y', 'Y'), each = 3)))
sampleTable$condition %<>% relevel(ref = 'Flg22_SynCom33')
rownames(sampleTable) <- colnames(kres$counts)

degres <- DESeqDataSetFromTximport(kres,
                                   sampleTable,
                                   ~condition)

## DEGs role out two zeros in one group
## degres <- degres[rowSums(counts(degres)) > 1, ]
degres %<>%
  estimateSizeFactors %>%
  counts(normalized = TRUE) %>%
  apply(1, checkFlg22, 1) %>%
  degres[., ]
degres <- DESeq(degres)

## count transformation
rld <- rlog(degres)
vst <- varianceStabilizingTransformation(degres)
ntd <- normTransform(degres)

resRaw <- degres %>%
  results(name = 'condition_Flg22_SynCom35_vs_Flg22_SynCom33') %T>%
  summary %>%
  as_tibble %>%
  select(pvalue, padj, log2FoldChange) %>%
  rename_all(.funs = list(~paste0('Flg22_SynCom35_vs_Flg22_SynCom33_', .)))

res <- cbind.data.frame(as.matrix(mcols(degres)[, 1:10]), assay(ntd), stringsAsFactors = FALSE) %>%
  rownames_to_column(., var = 'ID') %>%
  as_tibble %>%
  bind_cols(resRaw) %>%
  inner_join(anno, by = 'ID') %>%
  select(ID, Gene : Description, Flg22_SynCom33_1 : Flg22_SynCom35_vs_Flg22_SynCom33_log2FoldChange) %>%
  arrange(Flg22_SynCom35_vs_Flg22_SynCom33_padj)

write_csv(res, 'SynCom35_vs_SynCom33_k_full.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PCA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library('directlabels')
library('ggplot2')

rldDay8 <- str_detect(colnames(rld), 'Day8') %>% rld[, .]

pca <- prcomp(t(assay(rldDay8)))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
percentVar <- round(100 * percentVar)
pca1 <- pca$x[,1]
pca2 <- pca$x[,2]
pcaData <- data.frame(PC1 = pca1, PC2 = pca2, Group = colData(rldDay8)[, 1], ID = rownames(colData(rldDay8)))

ggplot(pcaData, aes(x = PC1, y = PC2, colour = Group)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_dl(aes(label = ID, color = Group), method = 'smart.grid')

ggsave('../results/PCA_Day8_1stadd.pdf', width = 12)
ggsave('../results/PCA_Day8_1stadd.jpg', width = 12)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################