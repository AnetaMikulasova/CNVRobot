args <- commandArgs(trailing = TRUE)

suppressPackageStartupMessages(library(tidyverse))

DIR     = args[1] %>% as.character()

SEX_FILES = list.files(DIR, recursive = T, full.names = T, pattern = ".cov")

DATASET = data.frame(NA,NA,NA,NA,NA,NA)
names(DATASET) = c("NAME","sexA","SRY","CTRL_GENE","sexB","sexC")

for (file in SEX_FILES) {
  NAME = basename(file) %>% gsub("_sextest.cov","",.)
  sexA = strsplit(NAME, "_") %>% as.vector() %>% unlist()
  sexA=tail(sexA, n=1)
  FILE = read.delim(file, header = F, stringsAsFactors = F)
  SRY=FILE[1,5] %>% as.integer()
  CTRL_GENE=FILE[2,5] %>% as.integer()
  sexB=SRY/CTRL_GENE
  sexC=if(sexB < 0.15) {"F"} else {"M"}
  x=cbind(NAME,sexA,SRY,CTRL_GENE,sexB,sexC)
  
  DATASET=rbind(DATASET,x)
  }

DATASET = DATASET %>%
  filter(NAME != "NA") %>%
  mutate(TEST = ifelse(sexA == sexC, "yes", "no"))

names(DATASET) = c("id", "sex_expected_from_metadata", "SRY_coverage", "CTRL_GENE_coverage" ,"SRY/CTRL_GENE_ratio", "sex_detected_by_coverage", "result")

correct = filter(DATASET, result == "yes")
wrong = filter(DATASET, result == "no")

write_tsv(correct, paste0(DIR, "results_correct.tsv"), col_names = T)
if (nrow(wrong)==0) {write_tsv(wrong, paste0(DIR, "results_wrong.tsv"), col_names = F)}
if (nrow(wrong)!=0) {write_tsv(wrong, paste0(DIR, "results_wrong.tsv"), col_names = T)}


