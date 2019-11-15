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
sampleTable <- data.frame(condition = factor(rep(condi, each = 6), levels = condi), batch = factor(rep(rep(1 : 2, each = 3), 8)))
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
ggsave('auto_sv_day8.jpg')
ggsave('auto_sv_day8.pdf')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
degres$sv1 <- svobj$sv[, 1]
degres$sv2 <- svobj$sv[, 2]
degres$sv3 <- svobj$sv[, 3]
degres$sv4 <- svobj$sv[, 4]
design(degres) <- ~sv1 + sv2 + sv3 + sv4 + condition

cond <- list(c('Col0_FeCl3_HK_Day8', 'Col0_FeEDTA_HK_Day8'),
             c('Col0_FeCl3_Live_Day8', 'Col0_FeEDTA_Live_Day8'),
             c('f6h1_FeCl3_HK_Day8', 'f6h1_FeEDTA_HK_Day8'),
             c('f6h1_FeCl3_Live_Day8', 'f6h1_FeEDTA_Live_Day8'),
             c('Col0_FeEDTA_Live_Day8', 'Col0_FeEDTA_HK_Day8'),
             c('Col0_FeCl3_Live_Day8', 'Col0_FeCl3_HK_Day8'),
             c('f6h1_FeEDTA_Live_Day8', 'f6h1_FeEDTA_HK_Day8'),
             c('f6h1_FeCl3_Live_Day8', 'f6h1_FeCl3_HK_Day8'))

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
  select(ID, Gene : Description, Col0_FeCl3_HK_Day8_Rep1 : f6h1_FeCl3_Live_Day8_vs_f6h1_FeCl3_HK_Day8_log2FoldChange) %>%
  arrange(Col0_FeCl3_HK_Day8_vs_Col0_FeEDTA_HK_Day8_padj)

write_csv(res, 'eachGroup_mergeDay8.csv')
save(degres, rldData, file = 'eachGroup_mergeDay8.RData')
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
         Batch = c(rep(rep(1:2, each = 3), 8))) %>%
  filter(Time == 'Day8')

cols <- brewer.pal(4, name = 'Set1')
cols[3:4] <- cols[4:3]

dat <- rld %>%
  assay %>%
  {.[rowMeans(.) > 1, ]}

group <- sampleTable$condition
design <- model.matrix(~ group)
rldData <- dat %>%
  removeBatchEffect(covariates = svobj$sv,
                    design = design)

## ## limma: remove batch effect
## rldData <- assay(str_detect(colnames(rld), 'Day8') %>% rld[, .]) %>%
##   removeBatchEffect(rep(rep(1:2, each = 3), 8) %>% factor)

## ## sva: remove batch effect
## modcombat <- model.matrix(~1, data = rldDay8)
## rldData <- assay(str_detect(colnames(rld), 'Day8') %>% rld[, .]) %>%
##   ComBat(dat = ., batch = rep(rep(1:2, each = 3), 8) %>% factor, mod = modcombat, par.prior = TRUE, prior.plots = FALSE)

pca <- prcomp(t(rldData))
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

ggsave('../results/PCA_mergeDay8_sva.pdf', width = 12)
ggsave('../results/PCA_mergeDay8_sva.jpg', width = 12)
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
