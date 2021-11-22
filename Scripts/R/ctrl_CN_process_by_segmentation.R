args <- commandArgs(trailing = TRUE)

suppressPackageStartupMessages(library(tidyverse))

INPUT    = args[1]
MINLOSS  = args[2] %>% as.numeric()
MINGAIN  = args[3] %>% as.numeric()
OUTPUT   = args[4]

message(paste0("   R ... ", date(), " - STAGE 1/3 - loading data"))
CTRL_CN = readRDS(INPUT)

message(paste0("   R ... ", date(), " - STAGE 2/3 - processing"))

#Get self-denois file
CTRL_CN_POS = select(CTRL_CN, CONTIG, START, END)
CTRL_CN_DATA = select(CTRL_CN, matches("denoisedCR"))

N = ncol(CTRL_CN_DATA)

CTRL_CN_PREP = CTRL_CN_POS %>% mutate("NofALL" = N)
CTRL_CN_PREP$NofCASES = rowSums(CTRL_CN_DATA == 0 | CTRL_CN_DATA > 0 | CTRL_CN_DATA < 0, na.rm=T)
CTRL_CN_PREP$NofLOSS = rowSums(CTRL_CN_DATA < MINLOSS, na.rm=T)
CTRL_CN_PREP$NofGAIN = rowSums(CTRL_CN_DATA > MINGAIN, na.rm=T)
CTRL_CN_PREP = CTRL_CN_PREP %>%
  mutate("FofCASES" = NofCASES/NofALL) %>%
  mutate("FofLOSS" = NofLOSS/NofCASES) %>%
  mutate("FofGAIN" = NofGAIN/NofCASES)

message(paste0("   R ... ", date(), " - STAGE 3/3 - writing output"))

write_tsv(CTRL_CN_PREP, OUTPUT, col_names = T)


#create IGV tracks
#---------------------------
#create IGV file - freq of cases with successful output
IGV_FofCASES = CTRL_CN_PREP %>% select(CONTIG, START, END, FofCASES) %>%
  mutate(START = START - 1) %>%
  arrange(CONTIG, START, END)
IGV_FofCASES_OUTPUT = paste0(gsub(".tsv", "", OUTPUT), "_freqofcases")
write_tsv(IGV_FofCASES, IGV_FofCASES_OUTPUT, col_names = F)
#create IGV file - losses
IGV_FofLOSS = CTRL_CN_PREP %>% select(CONTIG, START, END, FofLOSS) %>% filter(!is.na(FofLOSS)) %>%
  mutate(START = START - 1) %>%
  arrange(CONTIG, START, END)
IGV_FofLOSS_OUTPUT = paste0(gsub(".tsv", "", OUTPUT), "_freqofloss")
write_tsv(IGV_FofLOSS, IGV_FofLOSS_OUTPUT, col_names = F)
#create IGV file - losses
IGV_FofGAIN = CTRL_CN_PREP %>% select(CONTIG, START, END, FofGAIN) %>% filter(!is.na(FofGAIN)) %>%
  mutate(START = START - 1) %>%
  arrange(CONTIG, START, END)
IGV_FofGAIN_OUTPUT = paste0(gsub(".tsv", "", OUTPUT), "_freqofgain")
write_tsv(IGV_FofGAIN, IGV_FofGAIN_OUTPUT, col_names = F)










