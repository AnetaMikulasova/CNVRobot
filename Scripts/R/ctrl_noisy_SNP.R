args <- commandArgs(trailing = TRUE)

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(expss))  #vlookup

GNOMAD_FILTERED     = args[1] %>% as.character()
FILES               = args[2] %>% as.character()
GNOMAD_OUT          = args[3] %>% as.character()
AFDIF               = args[4] %>% as.numeric()
AFDEPTH             = args[5] %>% as.numeric()
AFPERC              = args[6] %>% as.numeric()
N                   = args[7] %>% as.numeric()


message(paste0("   R ... ", date(), " - STAGE 1/5 - loading data"))

#get maf data 
files_maf          <- list.files(FILES, recursive = T, full.names = T, pattern = ".rds")

#get filtered (for capture) gnomad
GNOMAD = read.delim(GNOMAD_FILTERED, comment.char = "#", head=FALSE, stringsAsFactors = F)

message(paste0("   R ... ", date(), " - STAGE 2/5 - processing gnomad"))

POSITIONS = GNOMAD %>%
  select(V1, V2) %>%
  mutate(POSITION = paste0(V1, "_", V2))
names(POSITIONS) = c("CONTIG", "POSITION_BP", "POSITION")

DATASET_MAF = POSITIONS %>% select(CONTIG, POSITION_BP)
DATASET_N = POSITIONS %>% select(CONTIG, POSITION_BP)

CASE_N = 0

message(paste0("   R ... ", date(), " - STAGE 3/5 - processing files of controls"))
for (file in files_maf) {
    CASE_N = CASE_N + 1
    NAME    <- basename(file) %>% gsub(".rds","",.)
    message("   R ...... ", paste0(date(), " - adding ", CASE_N, "/", N, " ", NAME))
    FILE = readRDS(file) %>%
      mutate(POSITION = paste0(CONTIG, "_", POSITION)) %>%
      mutate(MAF1 = ALT_COUNT/(ALT_COUNT+REF_COUNT)) %>%
      mutate(N = REF_COUNT + ALT_COUNT)

    WORK_MAF = select(FILE, MAF1)
    names(WORK_MAF) = c(NAME)
    DATASET_MAF = cbind(DATASET_MAF, WORK_MAF)

    WORK_N = select(FILE, N)
    names(WORK_N) = c(NAME)
    DATASET_N = cbind(DATASET_N, WORK_N)

}

message(paste0("   R ... ", date(), " - STAGE 4/5 - identifying noise"))


DATASET_MAF_WORK = DATASET_MAF %>% select(matches("allelicCounts"))
DATASET_N_WORK = DATASET_N %>% select(matches("allelicCounts"))

DATASET = POSITIONS %>% select(CONTIG, POSITION_BP) %>% mutate("NofCases" = N)

DATASET$MAF_N_without_NA = rowSums(DATASET_MAF_WORK >= 0, na.rm=T)
DATASET$MAF_N_out_AFDIF = rowSums((DATASET_MAF_WORK > (0+AFDIF) & DATASET_MAF_WORK < (0.5-AFDIF)) | (DATASET_MAF_WORK > (0.5+AFDIF) & DATASET_MAF_WORK < (1-AFDIF)), na.rm=T)

DATASET$DEPTH_N_without_NA = rowSums(DATASET_N_WORK >= 0, na.rm=T)
DATASET$DEPTH_N_below_cutoff_AFDEPTH = rowSums(DATASET_N_WORK < AFDEPTH)

NOISE = DATASET %>%
  mutate(FofMAF_Noise = MAF_N_out_AFDIF/NofCases) %>%
  mutate(FofDEPTH_Noise = DEPTH_N_below_cutoff_AFDEPTH/NofCases) %>%
  mutate(FILTER = ifelse(FofMAF_Noise > AFPERC | FofDEPTH_Noise > AFPERC, "FAILED", "PASS"))

FILTER = NOISE %>% select(FILTER)

GNOMAD = cbind(GNOMAD, FILTER)
GNOMAD = GNOMAD %>% 
  filter(FILTER == "PASS") %>%
  select(matches("V"))

message(paste0("   R ... ", date(), " - STAGE 5/5 - writing output"))

#output for filtering gnomad 
write_tsv(GNOMAD, path = paste0(GNOMAD_OUT), col_names = F)

