######################load K alignment###############################
library('tximport')
library('rhdf5')
library('magrittr')
library('DESeq2')
library('tibble')
library('readr')
library('dplyr')
library('stringr')
library('subSeq')
library('ggplot2')

##~~~~~~~~~~~~~~~~~~~~~~~~~read in sample annotation~~~~~~~~~~~~~~~~~
sampleAnno <- read_delim('/extDisk1/RESEARCH/MPIPZ_CJ_RNASeq/results/sample_anno.txt', delim = '\t') %>%
  rename(ID = `Library number`, SampleAnno = `Library Name`) %>%
  select(ID, SampleAnno) %>%
  mutate(ID = ID %>% str_replace('\\.', '_'), SampleAnno = SampleAnno %>% str_replace('-', '_')) %>%
  arrange(ID)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~read in k~~~~~~~~~~~~~~~~~~~~~~~~~~
wd <- '/home/Yulong/netscratch/CJFe/align_data'
setwd(wd)

slabel <- sampleAnno$ID %>%
  paste0('_ath_kallisto')

files <- file.path(wd, slabel, 'abundance.h5') %>%
  set_names(sampleAnno$SampleAnno)
kres <- tximport(files, type = 'kallisto', txOut = TRUE)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~DESeq2~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sampleTable <- data.frame(condition = sampleAnno$SampleAnno %>% substr(., 1, nchar(.) - 5) %>% factor)
rownames(sampleTable) <- colnames(kres$counts)

## unnormalized counts
degres <- DESeqDataSetFromTximport(kres, sampleTable, ~condition)

## degres %<>%
##   estimateSizeFactors %>%
##   counts(normalized = TRUE)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~subseq~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tvscCondi <- c('Col0_FeCl3_HK_Day8|Col0_FeCl3_Live_Day8', 'f6h1_FeCl3_HK_Day8|f6h1_FeCl3_Live_Day8', 'Col0_FeCl3_HK_Day15|Col0_FeCl3_Live_Day15', 'f6h1_FeCl3_HK_Day15|f6h1_FeCl3_Live_Day15')
savepath <- '/extDisk1/RESEARCH/MPIPZ_CJ_RNASeq/results'

for(i in seq_along(tvscCondi)) {

  eachIdx <- degres %>%
    colnames %>%
    str_detect(tvscCondi[i]) %>%
    which

  testCounts <- degres %>%
    counts %>%
    .[, eachIdx] %>%
    .[rowSums(.) >= 5, ]

  testCondi <- degres$condition %>%
    .[eachIdx] %>%
    droplevels

  proportions <- 10^seq(-2, 0, .1)
  ss <- subsample(testCounts, proportions, method = c('edgeR', 'voomLimma', 'DESeq2'), treatment = testCondi)

  seed <- getSeed(ss)
  ssglm <- subsample(testCounts, proportions, method = c('edgeR.glm'), mm= model.matrix(~ testCondi), seed = seed)
  ssComb <- combineSubsamples(ss, ssglm)

  plot(ssComb)
  ggsave(paste0(str_replace(tvscCondi[i], '\\|', '_vs_'), '.pdf') %>% file.path(savepath, .))
  ggsave(paste0(str_replace(tvscCondi[i], '\\|', '_vs_'), '.jpg') %>% file.path(savepath, .))
}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#####################################################################
