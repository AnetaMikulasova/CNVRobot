args <- commandArgs(trailing = TRUE)

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(expss))  #vlookup

INTERVALS     = args[1] %>% as.character()
INDIR         = args[2] %>% as.character()
PATTERN       = args[3] %>% as.character()
N             = args[4] %>% as.numeric()
OUT           = args[5] %>% as.character()

#get denoised data 
files_denois          <- list.files(INDIR, recursive = T, full.names = T, pattern = PATTERN)

#get intervals
POSITIONS = read.delim(INTERVALS, comment.char = "@", head=FALSE, stringsAsFactors = F) %>% 
  select(V1, V2, V3)
names(POSITIONS) = c("CONTIG", "START", "END")
POSITIONS = mutate(POSITIONS, POSITION = paste0(CONTIG, "_", START, "_", END)) %>% arrange(CONTIG, START, END)

CASE_N = 0

DATASET = POSITIONS %>% select(CONTIG, START, END)
  
for (file in files_denois) {
     CASE_N = CASE_N + 1
     NAME    <- basename(file) %>% gsub(".tsv","",.)
     message("   R ... ", paste0(date(), " - adding ", CASE_N, "/", N, " ", NAME))
     FILE = read.delim(file, comment.char = "@", stringsAsFactors = F) %>%
       mutate(POSITION = paste0(CONTIG, "_", START, "_", END))
     
     WORK = select(POSITIONS,POSITION) %>%
       mutate(x = vlookup(POSITION, dict = FILE, lookup_column = 5, result_column = 4))
     names(WORK) = c("POSITION", NAME) 
     WORK = WORK %>% select(NAME)
    
    DATASET = cbind(DATASET, WORK)
}


message("   R ... ", paste0(date(), " writing output"))
saveRDS(DATASET, file = paste0(OUT))



