args <- commandArgs(trailing = TRUE)

suppressPackageStartupMessages(library(tidyverse))

GET_DENOIS                = args[1] %>% as.character()

#transform
DENOIS_DATA_INPUT=read.delim(GET_DENOIS, stringsAsFactors = F) %>%
  select(chromosome, start, end, log2) %>%
  mutate(start = start + 1)
names(DENOIS_DATA_INPUT) = c("CONTIG", "START", "END", "LOG2_COPY_RATIO")

write_tsv(DENOIS_DATA_INPUT, paste0(GET_DENOIS), col_names = T)
system(paste0("touch status_ok"))

