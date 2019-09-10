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
  filter(Sample %>% str_detect('4206')) %>%
  arrange(Anno) %>%
  arrange(Days)

slabel <- labelanno$Sample %>%
  paste0('_ath_kallisto')

files <- file.path(wd, slabel, 'abundance.h5')
names(files) <- labelanno$Anno

kres <- tximport(files, type = 'kallisto', txOut = TRUE)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~conditions~~~~~~~~~~~~~~~~~
setwd('/extDisk1/RESEARCH/MPIPZ_CJ_RNASeq/results/')

## sampleTable
condi <- labelanno$Anno %>% substring(first = 1, last = nchar(.) - 5) %>% unique
sampleTable <- data.frame(condition = factor(rep(condi, each = 3), levels = condi))
sampleTable$condition %<>% relevel(ref = 'Col0_FeCl3_Live_Day14')
rownames(sampleTable) <- colnames(kres$counts)

degres <- DESeqDataSetFromTximport(kres, sampleTable, ~condition)

## DEGs
degres %<>%
  estimateSizeFactors %>%
  counts(normalized = TRUE) %>%
  apply(1, checkFe, 1) %>%
  degres[., ]
## degres <- degres[rowSums(counts(degres)) > 1, ]
save(degres, file = 'degres_1stadd_raw.RData')

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
  select(ID, Gene : Description, Col0_FeCl3_HK_Day8_Rep4 : f6h1_FeEDTA_Live_Day14_vs_Col0_FeCl3_Live_Day14_log2FoldChange) %>%
  arrange(Col0_FeEDTA_Live_Day14_vs_Col0_FeCl3_Live_Day14_padj)

write_csv(res, 'eachGroup_vs_Col0_FeCl3_Live_Day14_k.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~merge pair-wise compare~~~~~~~~~~~~~~~~~~~~~~
setwd('/extDisk1/RESEARCH/MPIPZ_CJ_RNASeq/results/')

DEGf <- dir(pattern = 'eachGroup_vs')

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
  select(ID : f6h1_FeEDTA_Live_Day14_Rep6) %>%
  bind_cols(mergeRes) %>%
  write_csv('eachGroup_vs_iron.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PCA~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library('directlabels')
library('ggplot2')
library('RColorBrewer')

rldDay14 <- colData(rld)[, 1] %>%
  as.character %>%
  strsplit(split = '_', fixed = TRUE) %>%
  do.call(rbind, .) %>%
  set_colnames(c('Genotype', 'Iron', 'SynCom', 'Time')) %>%
  as_tibble %>%
  mutate(Conditions = paste(Iron, SynCom, sep = '_') %>%
           factor,
         Genotype = Genotype %<>% factor) %>%
  filter(Time == 'Day14')

cols <- brewer.pal(4, name = 'Set1')
cols[3:4] <- cols[4:3]

pca <- prcomp(t(assay(str_detect(colnames(rld), 'Day14') %>% rld[, .])))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
percentVar <- round(100 * percentVar)
pca1 <- pca$x[,1]
pca2 <- pca$x[,2]

pcaData <- rldDay14 %>%
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
  ggtitle('Day14') +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'))

ggsave('../results/PCA_Day14_1stadd.pdf', width = 12)
ggsave('../results/PCA_Day14_1stadd.jpg', width = 12)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################
