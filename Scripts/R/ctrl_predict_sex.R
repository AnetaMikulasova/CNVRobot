args <- commandArgs(trailing = TRUE)

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rhdf5))
suppressPackageStartupMessages(library(expss))



GETRAWCOUNT     = args[1] %>% as.character()

#define autosomes
AUTOSOMES_NUM= c(1:22)
AUTOSOMES = c(AUTOSOMES_NUM, as.vector(outer("chr", AUTOSOMES_NUM, paste, sep="")))


if (grepl(".hdf5$", GETRAWCOUNT)) {
  RAWCOUNT_DATA_INPUT_COUNTS = h5read(GETRAWCOUNT, "/counts/values")
  RAWCOUNT_DATA_INPUT_INTERVALS = h5read(GETRAWCOUNT, "/intervals/transposed_index_start_end")
  RAWCOUNT_DATA_INPUT_CONTIGS = h5read(GETRAWCOUNT, "/intervals/indexed_contig_names") %>% as.data.frame() 
  names(RAWCOUNT_DATA_INPUT_CONTIGS) = c("CONTIG")
  RAWCOUNT_DATA_INPUT_CONTIGS = RAWCOUNT_DATA_INPUT_CONTIGS %>% mutate(INDEX = row_number()-1)
  RAWCOUNT_DATA_INPUT = cbind(RAWCOUNT_DATA_INPUT_INTERVALS, RAWCOUNT_DATA_INPUT_COUNTS) %>% as.data.frame()
  names(RAWCOUNT_DATA_INPUT) = c("INDEX", "START", "END", "COUNT")
  RAWCOUNT_DATA_INPUT = RAWCOUNT_DATA_INPUT %>% 
    mutate(CONTIG = vlookup(INDEX, dict = RAWCOUNT_DATA_INPUT_CONTIGS, lookup_column = which(colnames(RAWCOUNT_DATA_INPUT_CONTIGS)=="INDEX"), result_column = which(colnames(RAWCOUNT_DATA_INPUT_CONTIGS)=="CONTIG"))) %>%
    select(CONTIG, START, END, COUNT)
}

if (grepl(".cnn$", GETRAWCOUNT)) {
  RAWCOUNT_DATA_INPUT=read.delim(GETRAWCOUNT, stringsAsFactors = F) %>%
    select(chromosome, start, end, depth)
  names(RAWCOUNT_DATA_INPUT) = c("CONTIG", "START", "END", "COUNT")
}



RAWCOUNT_DATA_INPUT = RAWCOUNT_DATA_INPUT %>%
  mutate(SIZE_KB = (END-START+1)/1000) %>%
  mutate(COUNT_KB = COUNT/SIZE_KB) %>%
  mutate(CONTIG_CATEG = case_when(
    CONTIG %in% AUTOSOMES ~ "A",
    CONTIG %in% c("X", "chrX") ~ "X",
    CONTIG %in% c("Y", "chrY") ~ "Y",
    TRUE ~ "OTHER"
  ))

MEDIAN_COUNTS <- RAWCOUNT_DATA_INPUT %>%
  group_by(CONTIG_CATEG) %>%
  summarize(MEDIAN_KB = median(COUNT_KB, na.rm = TRUE), N = dplyr::n()) 

# MEDIAN_COUNTS = MEDIAN_COUNTS %>% filter(CONTIG_CATEG == "X")


CHR_X_N = nrow(RAWCOUNT_DATA_INPUT %>% filter(CONTIG_CATEG == "X")) %>% as.numeric()
CHR_Y_N = nrow(RAWCOUNT_DATA_INPUT %>% filter(CONTIG_CATEG == "Y")) %>% as.numeric()
CHR_A_N = nrow(RAWCOUNT_DATA_INPUT %>% filter(CONTIG_CATEG == "A")) %>% as.numeric()




#if autosomal coverage exist
if(CHR_A_N > 0) {

  AUTOSOME = MEDIAN_COUNTS %>% filter(CONTIG_CATEG == "A") %>% select(MEDIAN_KB) %>% as.numeric()
  
  #if there is any gonosomal coverage
  if(CHR_X_N+CHR_Y_N > 0) {
    
    #if X has more intervals (which it will have as it is bigger, so X is default for sex determination - which ok as XXY would be "F" and not having X doubled)
    if(CHR_X_N >= CHR_Y_N) {
      #determine sex based on X/A ratio
      CHR_X = MEDIAN_COUNTS %>% filter(CONTIG_CATEG == "X") %>% select(MEDIAN_KB) %>% as.numeric()
      if(CHR_X/AUTOSOME > 0.75) {SEX = "F"} else {SEX = "M"}
    }
    #if Y has more intervals, then use Y (only in cases where X is "under-covered")
    if(CHR_X_N < CHR_Y_N) {
      #determine sex based on Y/A ratio
      CHR_Y = MEDIAN_COUNTS %>% filter(CONTIG_CATEG == "Y") %>% select(MEDIAN_KB) %>% as.numeric()
      if(CHR_Y/AUTOSOME < 0.15) {SEX = "F"} else {SEX = "M"}
    }
  }

  #if there is no gonosomal coverage, sex stays unknown
  if(CHR_X_N+CHR_Y_N == 0) {SEX = "unk"}
  
}

#if no autosomal coverage, try to use gonosomes
if(CHR_A_N == 0) {
  
  #if both gonosomes have any coverage, predict; if only one gonosome coveraged, make it unknown
  if(CHR_X_N > 0 & CHR_Y_N > 0) {
    CHR_X = MEDIAN_COUNTS %>% filter(CONTIG_CATEG == "X") %>% select(MEDIAN_KB) %>% as.numeric()
    CHR_Y = MEDIAN_COUNTS %>% filter(CONTIG_CATEG == "Y") %>% select(MEDIAN_KB) %>% as.numeric()
    if(CHR_Y/CHR_X < 0.15) {SEX = "F"} else {SEX = "M"}
  } else {
    SEX = "unk"
  }
  
}




system(paste0("touch status_ok"))

cat(SEX)




