args <- commandArgs(trailing = TRUE)

suppressPackageStartupMessages(library(tidyverse))

GET_REF_FAI               = args[1] %>% as.character()
GET_OUTPUT_CONTIG         = args[2] %>% as.character()
GET_OUTPUT_CONTIG_SIZES   = args[3] %>% as.character()

#define list of possible contigs
CONTIGS_NUM = c(1:22,"X","Y")
CONTIGS = c(CONTIGS_NUM, as.vector(outer("chr", CONTIGS_NUM, paste, sep="")))

#extract contigs (46,XY) and their length
REF_FAI = read.delim(GET_REF_FAI, header = F, stringsAsFactors = F) %>%
  select(V1, V2) %>%
  filter(V1 %in% CONTIGS)
write_tsv(REF_FAI, paste0(GET_OUTPUT_CONTIG_SIZES), col_names = F) #replace the chromosome size file

REF_FAI = REF_FAI %>%
  mutate(V3 = 0) %>%
  select(V1, V3, V2)
write_tsv(REF_FAI, paste0(GET_OUTPUT_CONTIG), col_names = F) #replace the chromosome size file



