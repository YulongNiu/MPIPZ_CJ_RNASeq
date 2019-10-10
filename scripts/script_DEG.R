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
  arrange(Days)

slabel <- labelanno$Sample %>%
  paste0('_ath_kallisto')

files <- file.path(wd, slabel, 'abundance.h5')
names(files) <- labelanno$Anno

kres <- tximport(files, type = 'kallisto', txOut = TRUE)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~pair-wise comparison~~~~~~~~~~~~~~~~~
setwd('/extDisk1/RESEARCH/MPIPZ_CJ_RNASeq/results/')

## sampleTable
condi <- labelanno$Anno %>% substring(first = 1, last = nchar(.) - 5) %>% unique
sampleTable <- data.frame(condition = factor(rep(condi, each = 3), levels = condi))
sampleTable$condition %<>% relevel(ref = 'Col0_FeCl3_HK_Day8')
rownames(sampleTable) <- colnames(kres$counts)

degres <- DESeqDataSetFromTximport(kres, sampleTable, ~condition)

## DEGs
degres %<>%
  estimateSizeFactors %>%
  counts(normalized = TRUE) %>%
  apply(1, checkFe, 1) %>%
  degres[., ]
## degres <- degres[rowSums(counts(degres)) > 1, ]
save(degres, file = '../results/degres.RData')

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
  select(ID, Gene : Description, Col0_FeCl3_HK_Day8_Rep1 : f6h1_FeEDTA_Live_Day15_vs_Col0_FeCl3_HK_Day15_log2FoldChange) %>%
  arrange(Col0_FeEDTA_Live_Day15_vs_Col0_FeCl3_HK_Day15_padj)

write_csv(res, 'eachGroup_vs_Col0_FeCl3_HK_Day15_k.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~merge pair-wise compare~~~~~~~~~~~~~~~~~~~~~~
setwd('/extDisk1/RESEARCH/MPIPZ_CJ_RNASeq/results/')

DEGf <- dir(pattern = 'eachGroup_vs')[c(1, 2, 5, 6)]

mergeRes <- foreach(i = seq_along(DEGf), .combine = bind_cols) %do% {

  eachpre <- DEGf[i] %>%
    str_extract('(?<=vs_).*?(?=_k)') %>%
    {
      fpart <- str_replace(., 'FeCl3', 'FeEDTA')
      cond <- paste(fpart, ., sep = '_vs_')
    }

  eachf <- read_csv(DEGf[i]) %>%
    select(starts_with(eachpre))
}

read_csv('eachGroup_vs_Col0_FeCl3_Live_Day8_k.csv',
         col_types = cols(Chromosome = col_character())) %>%
  mutate(Gene = Gene %>% {if_else(is.na(.), '', .)}) %>%
  mutate(Description = Description %>% {if_else(is.na(.), '', .)}) %>%
  select(ID : f6h1_FeEDTA_Live_Day15_Rep3) %>%
  bind_cols(mergeRes) %>%
  write_csv('eachGroup_vs_iron.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PCA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library('directlabels')
library('ggplot2')

rldDay15 <- str_detect(colnames(rld), 'Day15') %>% rld[, .]

pca <- prcomp(t(assay(rldDay15)))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
percentVar <- round(100 * percentVar)
pca1 <- pca$x[,1]
pca2 <- pca$x[,2]
pcaData <- data.frame(PC1 = pca1, PC2 = pca2, Group = colData(rldDay15)[, 1], ID = rownames(colData(rldDay15)))

ggplot(pcaData, aes(x = PC1, y = PC2, colour = Group)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_dl(aes(label = ID, color = Group), method = 'smart.grid')

ggsave('../results/PCA_Day15.pdf', width = 12)
ggsave('../results/PCA_Day15.jpg', width = 12)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################
