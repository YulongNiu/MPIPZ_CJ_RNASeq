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
  mutate_all(list(~ str_replace(., 'Day14', 'Day15'))) %>%
  mutate(Days = rep(c(8, 15, 8, 15), each = 24)) %>%
  filter(Timepoint %>% str_detect('Day15')) %>%
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
sampleTable <- data.frame(condition = factor(rep(condi, each = 6), levels = condi), batch = factor(rep(rep(1 : 2, each = 3), 8)))
rownames(sampleTable) <- colnames(kres$counts)

degres <- DESeqDataSetFromTximport(kres, sampleTable, ~ condition)

## DEGs
degres %<>%
  estimateSizeFactors %>%
  counts(normalized = TRUE) %>%
  apply(1, checkFe, 1) %>%
  degres[., ]
## degres <- degres[rowSums(counts(degres)) > 1, ]
save(degres, file = 'degres_merge.RData')

degres <- DESeq(degres)
## resultsNames(degres)

## count transformation
rld <- rlog(degres)
vst <- varianceStabilizingTransformation(degres)
ntd <- normTransform(degres)

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
                     rename_all(.funs =  list(~paste0(paste(x, collapse = '_vs_'), '_', .)))
                 }) %>%
  bind_cols

res <- cbind.data.frame(as.matrix(mcols(degres)[, 1:10]), assay(ntd), stringsAsFactors = FALSE) %>%
  rownames_to_column(., var = 'ID') %>%
  as_tibble %>%
  bind_cols(resRaw) %>%
  inner_join(anno, by = 'ID') %>%
  select(ID, Gene : Description, Col0_FeCl3_HK_Day15_Rep1 : f6h1_FeCl3_Live_Day15_vs_f6h1_FeCl3_HK_Day15_log2FoldChange) %>%
  arrange(Col0_FeCl3_HK_Day15_vs_Col0_FeEDTA_HK_Day15_padj)

write_csv(res, 'eachGroup_mergeDay15.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~merge pair-wise comparison~~~~~~~~~~~~~~~~~~~~~~
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
library('limma')
library('sva')

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Day8~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rldDay8 <- colData(rld)[, 1] %>%
  as.character %>%
  strsplit(split = '_', fixed = TRUE) %>%
  do.call(rbind, .) %>%
  set_colnames(c('Genotype', 'Iron', 'SynCom', 'Time')) %>%
  as_tibble %>%
  mutate(Conditions = paste(Iron, SynCom, sep = '_') %>%
           factor,
         Genotype = Genotype %<>% factor,
         Batch = c(rep(rep(1:2, each = 3), 8), rep(1, 24), rep(2, 24))) %>%
  filter(Time == 'Day8')

cols <- brewer.pal(4, name = 'Set1')
cols[3:4] <- cols[4:3]

## raw
rldarray <- assay(str_detect(colnames(rld), 'Day8') %>% rld[, .])

## limma: remove batch effect
rldarray <- assay(str_detect(colnames(rld), 'Day8') %>% rld[, .]) %>%
  removeBatchEffect(rep(rep(1:2, each = 3), 8) %>% factor)

## sva: remove batch effect
modcombat <- model.matrix(~1, data = rldDay8)
rldarray <- assay(str_detect(colnames(rld), 'Day8') %>% rld[, .]) %>%
  ComBat(dat = ., batch = rep(rep(1:2, each = 3), 8) %>% factor, mod = modcombat, par.prior = TRUE, prior.plots = FALSE)

pca <- prcomp(t(rldarray))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
percentVar <- round(100 * percentVar)
pca1 <- pca$x[,1]
pca2 <- pca$x[,2]

pcaData <- rldDay8 %>%
  mutate(Colours = factor(Conditions, labels = cols)) %>%
  select(Conditions, Genotype, Colours, Batch) %>%
  mutate(PC1 = pca1, PC2 = pca2)

ggplot(pcaData, aes(x = PC1, y = PC2, colour = Conditions)) +
  geom_point(aes(shape = Genotype), size = 4) +
  scale_colour_manual(values = levels(pcaData$Colours),
                      name = 'Experimental\nCondition',
                      labels = expression(FeCl[3]+HK, FeCl[3]+Live, FeEDTA+HK, FeEDTA+Live)) +
  scale_shape_manual(values = c(15, 17)) +
  geom_text(aes(label = Batch), hjust = -0.5, vjust = 0) +
  xlab(paste0('PC1: ',percentVar[1],'% variance')) +
  ylab(paste0('PC2: ',percentVar[2],'% variance')) +
  ggtitle('Day8') +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'))

ggsave('../results/PCA_mergeDay8_raw.pdf', width = 12)
ggsave('../results/PCA_mergeDay8_raw.jpg', width = 12)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~day 14/15~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rldDay14 <- colData(rld)[, 1] %>%
  as.character %>%
  strsplit(split = '_', fixed = TRUE) %>%
  do.call(rbind, .) %>%
  set_colnames(c('Genotype', 'Iron', 'SynCom', 'Time')) %>%
  as_tibble %>%
  mutate(Conditions = paste(Iron, SynCom, sep = '_') %>%
           factor,
         Genotype = Genotype %<>% factor,
         Batch = c(rep(rep(1:2, each = 3), 8), rep(1, 24), rep(2, 24))) %>%
  filter(Time %in% c('Day14', 'Day15'))

cols <- brewer.pal(4, name = 'Set1')
cols[3:4] <- cols[4:3]

## raw
rldarray <- assay(str_detect(colnames(rld), 'Day14|Day15') %>% rld[, .])

## limma: remove batch effect
rldarray <- assay(str_detect(colnames(rld), 'Day14|Day15') %>% rld[, .]) %>%
  removeBatchEffect(rep(1:2, each = 24) %>% factor)

## sva: remove batch effect
modcombat <- model.matrix(~1, data = rldDay8)
rldarray <- assay(str_detect(colnames(rld), 'Day14|Day15') %>% rld[, .]) %>%
  ComBat(dat = ., batch = rep(1:2, each = 24) %>% factor, mod = modcombat, par.prior = TRUE, prior.plots = FALSE)

pca <- prcomp(t(rldarray))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
percentVar <- round(100 * percentVar)
pca1 <- pca$x[,1]
pca2 <- pca$x[,2]

pcaData <- rldDay14 %>%
  mutate(Colours = factor(Conditions, labels = cols)) %>%
  select(Conditions, Genotype, Colours, Batch) %>%
  mutate(PC1 = pca1, PC2 = pca2)

ggplot(pcaData, aes(x = PC1, y = PC2, colour = Conditions)) +
  geom_point(aes(shape = Genotype), size = 4) +
  scale_colour_manual(values = levels(pcaData$Colours),
                      name = 'Experimental\nCondition',
                      labels = expression(FeCl[3]+HK, FeCl[3]+Live, FeEDTA+HK, FeEDTA+Live)) +
  scale_shape_manual(values = c(15, 17)) +
  geom_text(aes(label = Batch), hjust = -0.5, vjust = 0) +
  xlab(paste0('PC1: ',percentVar[1],'% variance')) +
  ylab(paste0('PC2: ',percentVar[2],'% variance')) +
  ggtitle('Day14/15') +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'))

ggsave('../results/PCA_mergeDay14_brsva.pdf', width = 12)
ggsave('../results/PCA_mergeDay14_brsva.jpg', width = 12)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################
