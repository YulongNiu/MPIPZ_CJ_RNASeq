
## originally by Yulong Niu
## yulong.niu@hotmail.com

library('ontologyIndex')
library('ontologyPlot')
library('tidyverse')

data(go)

setwd('/extDisk1/RESEARCH/MPIPZ_CJ_RNASeq/geneset/')

files <- dir() %>%
  .[str_detect(., 'CC.csv')]

for (i in files) {

  prefix <- i %>%
    strsplit(i, split = '.', fixed = TRUE) %>%
    sapply(., '[[', 1)

  clusterGO <- read_csv(i) %>%
    select(-X1)

  ## all inter nodes
  selectedGO <- clusterGO %>%
    ## filter(str_detect(Annotation, 'iron')) %>%
    filter(over_represented_pvalue < 0.05) %>%
    filter(numInCat >= 3)

  selectedTerm <- selectedGO %>%
    .$category %>%
    get_ancestors(go, .) %>%
    remove_links(go, terms = .)

  termCol <- selectedTerm %>%
    {rep('powderblue', length(.))} %>%
    {
      sig <- selectedTerm %in% selectedGO$category
      .[sig] <- 'orange'
      return(.)
    }

  ## plot
  BPplot <- onto_plot(go,
                      terms = selectedTerm,
                      fillcolor = termCol,
                      fontsize = 10)

  write_dot(BPplot, paste0(prefix, '.dot'))

  ## convert to pdf
  command <- str_replace_all('dot -Tpdf PREFIX.dot -o PREFIX.pdf', 'PREFIX', prefix) %>%
    system
}

