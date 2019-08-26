###########################Raw reads######################
rawpath <- '/netscratch/dep_psl/grp_rgo/yniu/CJFe/raw_data_1stadd'
setwd(rawpath)

library('magrittr')
library('doParallel')
library('foreach')
library('tibble')
library('readr')

ncore <- 45

fqs <- dir(rawpath,
           pattern = 'R1.fq.gz',
           full.names = TRUE)

registerDoParallel(cores = ncore)

rn <- foreach(i = seq_along(fqs), .combine = c) %dopar% {

  eachrn <- paste('zcat', fqs[i], '| awk "END{print NR/4}"') %>%
    system(inter = TRUE) %>%
    as.numeric

  return(eachrn)
}

stopImplicitCluster()

snames <- fqs %>%
  strsplit(split = '/', fixed = TRUE) %>%
  sapply('[[', 8) %>%
  strsplit(split = '.', fixed = TRUE) %>%
  sapply('[[', 1) %>%
  {substr(., 1, nchar(.) - 3)}

tibble(sample = snames,
       rawfq = rn) %>%
  write_csv('raw_seqnumber_1stadd.csv')
##########################################################

#################extract Kallisto and HISAT2 output########
setwd('/extDisk1/RESEARCH/MPIPZ_CJ_RNASeq/results/')

library('dplyr')
library('readr')

KHoutput <- function(op, type = 'PE', org = 'hsa') {

  ## INPUT: 'op' is a character vector. 'type' is 'PE' (pair-end) or 'SE' (single-end). 'org' is the organism name.
  ## OUTPUT: A tibble, 1st column is input reads number, 2nd column is Kallisto aligned read number, and 3rd column is the HISAT2 aligned read number.
  ## USAGE: Extract the number of aligned reads.

  require('stringr')
  require('magrittr')
  require('tibble')

  ## input reads number
  fqnum <- op %>%
    str_detect('\\d+ reads; of these:') %>%
    op[.] %>%
    strsplit(split = ' ', fixed = TRUE) %>%
    sapply('[[', 1) %>%
    as.numeric

  ## HISAT2 aligned
  hmapnum <- op %>%
    str_detect('.* aligned 0 times$') %>%
    op[.] %>%
    strsplit(split = ' ', fixed = TRUE) %>%
    {  if (type == 'PE') {
         sapply(., '[[', 9)
       } else {
         sapply(., '[[', 5)
       }} %>%
    as.numeric %>%
    {  if (type == 'PE') {
         ./2
       } else .} %>%
    {fqnum - .}

  ## Kallisto aligned
  kmapnum <- op %>%
    str_detect('.* reads pseudoaligned') %>%
    op[.] %>%
    strsplit(split = ' ', fixed = TRUE) %>%
    sapply('[[', 5) %>%
    str_replace_all(',', '') %>%
    as.numeric

  ## sample names
  snames <- op %>%
    str_detect('HISAT2 using') %>%
    op[.] %>%
    strsplit(split = ' ', fixed = TRUE) %>%
    sapply('[[', 6) %>%
    {substr(., 1, nchar(.) - 1)}

  res <- tibble(sample = snames,
                trimfq = fqnum,
                hmap = hmapnum,
                kmap = kmapnum,
                org = org)

  return(res)
}

##~~~~~~~~~~~~~~~~~~~~~~~test contamination~~~~~~~~~~~~~~~~~~~~~~~~~~~
sampleAnno <- read_delim('/extDisk1/RESEARCH/MPIPZ_CJ_RNASeq/results/sample_anno.txt', delim = '\t') %>%
  rename(ID = `Library number`, SampleAnno = `Library Name`) %>%
  select(ID, SampleAnno) %>%
  mutate(ID = ID %>% str_replace('\\.', '_'), SampleAnno = SampleAnno %>% str_replace('-', '_')) %>%
  arrange(ID)

athout <- 'alignment_nohup_1stadd.out' %>%
  readLines %>%
  KHoutput(type = 'PE', org = 'ath') %>%
  mutate(H_ath = round(hmap/trimfq, 3), K_ath = round(kmap/trimfq, 3)) %>%
  select(c(-hmap, -kmap, -org))

## raw reads
rawrd <- read_csv('raw_seqnumber_1stadd.csv')

contam <- rawrd %>%
  inner_join(athout)

write_csv(contam, 'ath_alignment_1stadd.csv')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######################################################
