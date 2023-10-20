args <- commandArgs(trailing = TRUE)

message(paste0("   R ... ", date(), " - STARTING"))

suppressPackageStartupMessages(library(tidyverse))

GET_CONTIG                      = args[1] %>% as.character()
GET_PRIMARY_DENOIS              = args[2] %>% as.character()
GET_SECONDARY_DENOIS            = args[3] %>% as.character()
GET_OUTPUT_DENOIS               = args[4] %>% as.character()
# GET_PRIMARY_STANDARD            = args[5] %>% as.character()
# GET_SECONDARY_STANDARD          = args[6] %>% as.character()
# GET_OUTPUT_STANDARD             = args[7] %>% as.character()

message(paste0("   R ... ", date(), " - STAGE 1/3 - loading data"))

CONTIGS = read.delim(GET_CONTIG, header = F, stringsAsFactors = F) %>% select(V1)
CONTIGS = CONTIGS$V1 %>% as.vector()
PRIMARY_DENOIS   <- read.delim(paste0(GET_PRIMARY_DENOIS), comment.char = "@", stringsAsFactors = F, colClasses = "character") 
SECONDARY_DENOIS   <- read.delim(paste0(GET_SECONDARY_DENOIS), comment.char = "@", stringsAsFactors = F, colClasses = "character") 
# PRIMARY_STANDARD <- read.delim(paste0(GET_PRIMARY_STANDARD), comment.char = "@", stringsAsFactors = F, colClasses = "character")
# SECONDARY_STANDARD <- read.delim(paste0(GET_SECONDARY_STANDARD), comment.char = "@", stringsAsFactors = F, colClasses = "character") 


message(paste0("   R ... ", date(), " - STAGE 2/3 - process data"))

JOIN_DATA = function(DATA_PRIMARY, DATA_SECONDARY) {
  DATA_PRIMARY_AUTOSOMES = DATA_PRIMARY %>%
    filter(CONTIG != "chrX" & CONTIG != "X" & CONTIG != "chrY" & CONTIG != "Y")
  DATA_SECONDARY_GONOSOMES = DATA_SECONDARY %>%
    filter(CONTIG == "chrX" | CONTIG == "X" | CONTIG == "chrY" | CONTIG == "Y")
  DATA_JOIN = rbind(DATA_PRIMARY_AUTOSOMES, DATA_SECONDARY_GONOSOMES) %>%
    mutate(CONTIG = factor(CONTIG, levels = CONTIGS)) %>%
    arrange(CONTIG) %>%
    mutate(CONTIG = as.character(CONTIG))
  DATA_JOIN = mutate_all(DATA_JOIN, as.character)
}

JOIN_DENOIS_DATA = JOIN_DATA(PRIMARY_DENOIS, SECONDARY_DENOIS)
# JOIN_STANDARD_DATA = JOIN_DATA(PRIMARY_STANDARD, SECONDARY_STANDARD)

message(paste0("   R ... ", date(), " - STAGE 3/3 - output"))
write_tsv(JOIN_DENOIS_DATA, GET_OUTPUT_DENOIS, col_names = T)
# write_tsv(JOIN_STANDARD_DATA, GET_OUTPUT_STANDARD, col_names = T)

message(paste0("   R ... ", date(), " - FINISHED"))
system(paste0("touch status_ok"))


