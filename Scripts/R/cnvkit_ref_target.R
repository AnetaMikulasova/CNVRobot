args <- commandArgs(trailing = TRUE)

suppressPackageStartupMessages(library(tidyverse))

GET_INTERVALS             = args[1] %>% as.character()
GET_OUTPUT_TARGET         = args[2] %>% as.character()

#read GATK intervals and tranform to bed file for CNVkit
INTERVALS = read.delim(GET_INTERVALS, comment.char = "@", head=FALSE, stringsAsFactors = F) %>% 
  select(V1, V2, V3) %>%
  mutate(V2 = V2-1) %>%
  mutate(V4 = "-")

write_tsv(INTERVALS, paste0(GET_OUTPUT_TARGET), col_names = F)
system(paste0("touch status_ok"))
