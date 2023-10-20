args <- commandArgs(trailing = TRUE)

message(paste0("   R ... ", date(), " - STARTING"))

suppressPackageStartupMessages(library(expss))  #vlookup
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rhdf5))
suppressPackageStartupMessages(library(expss))


GETBANDUCSC          = args[1] %>% as.character()
GETCAPTUREID         = args[2] %>% as.character()
GETPROJECTID         = args[3] %>% as.character()
GETPROJECTTYPE       = args[4] %>% as.character()
GETRAWCOUNT          = args[5] %>% as.character()
GETALLELIC           = args[6] %>% as.character()
# GETSTANDARD          = args[7] %>% as.character() #deleted
GETDENOIS            = args[8] %>% as.character()
GETGROUPID           = args[9] %>% as.character()
GETCASEID            = args[10] %>% as.character()
GETCASETYPE          = args[11] %>% as.character()
GETCASESEX           = args[12] %>% as.character()
GETPONSEX            = args[13] %>% as.character()
NofCTRL_TOTAL        = args[14] %>% as.numeric()
NofCTRL_USED         = args[15] %>% as.numeric()
DIFFERENCE           = args[16] %>% as.numeric()
MINSIZE              = args[17] %>% as.numeric()
MINSIZE_SUB          = args[18] %>% as.numeric()
MINPROBE             = args[19] %>% as.numeric()
MINPROBE_SUB         = args[20] %>% as.numeric()
MINKEEP              = args[21] %>% as.numeric()
MINLOSS              = args[22] %>% as.numeric()
MINLOSS_SUB          = args[23] %>% as.numeric()
MINLOSS_BIAL         = args[24] %>% as.numeric()
MINGAIN              = args[25] %>% as.numeric()
MINGAIN_SUB          = args[26] %>% as.numeric()
GAP                  = args[27] %>% as.numeric()
SMOOTHPERC           = args[28] %>% as.numeric()
AFDIF                = args[29] %>% as.numeric()
AFSIZE               = args[30] %>% as.numeric()
AFPROBES             = args[31] %>% as.numeric()
GETSEGMENTCONDITION  = args[32] %>% as.character()
QCTABLEIN            = args[33] %>% as.character()
AUTOSOME_MODE        = args[34] %>% as.character()
OUTPUT_STANDARD      = args[35] %>% as.character()


message(paste0("   R ... ", date(), " - STAGE 1/5 - loading data"))

#Get data (only for GATK)
if(file.exists(GETRAWCOUNT)==T) {
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

# STANDARD_DATA_INPUT <- read.delim(paste0(GETSTANDARD), comment.char = "@", stringsAsFactors = F) #%>% mutate(CONTIG = ifelse(grepl("^chr*", CONTIG), CONTIG, paste0("chr", CONTIG)))
DENOIS_DATA_INPUT   <- read.delim(paste0(GETDENOIS), comment.char = "@", stringsAsFactors = F) #%>% mutate(CONTIG = ifelse(grepl("^chr*", CONTIG), CONTIG, paste0("chr", CONTIG)))

ALL_DATA_INPUT      <-   readRDS(paste0(GETALLELIC)) %>%
  mutate(N = REF_COUNT + ALT_COUNT) %>%
  filter(N > 0) %>%
  mutate(MAF = ALT_COUNT/(ALT_COUNT+REF_COUNT)) %>%
  mutate(MAFTRANS = abs(0.5-MAF)) %>%
  mutate(GENOTYPE = ifelse(MAFTRANS < 0.5*AFDIF, "NORMAL", "ABNORMAL"))

if(file.exists(QCTABLEIN)==F) {QCTABLE = data.frame()}
if(file.exists(QCTABLEIN)==T) {
  QCTABLE = read.delim(paste0(QCTABLEIN), stringsAsFactors = F) %>%
    mutate(CASE_SEX = ifelse(CASE_SEX == "FALSE", "F", CASE_SEX)) %>%
    mutate(PON_SEX = ifelse(PON_SEX == "FALSE", "F", PON_SEX))
}

#define smooth coefficient for the difference in abnormal segments
#why? Because in germline analysis I expect gain, loss of one or two allele; but in cancer there can be more CN states due to aneuploidies
if(GETPROJECTTYPE == "germline") {
  DIFFERENCE_SMOOTH_COEF = 2
  DIFFERENCE_SMOOTH_COEF_SUB = 1.5
}
if(GETPROJECTTYPE != "germline") {
  DIFFERENCE_SMOOTH_COEF = 1
  DIFFERENCE_SMOOTH_COEF_SUB = 1
}


#Functions
rowShift <- function(x, shiftLen = 1L) {
  r <- (1L + shiftLen):(length(x) + shiftLen)
  r[r<1] <- NA
  return(x[r])
}

fillNAgaps <- function(x, firstBack=FALSE) {
  lvls <- NULL
  if (is.factor(x)) {
    lvls <- levels(x)
    x    <- as.integer(x)
  }
  goodIdx <- !is.na(x)
  if (firstBack)   goodVals <- c(x[goodIdx][1], x[goodIdx])
  else             goodVals <- c(NA,            x[goodIdx])
  fillIdx <- cumsum(goodIdx)+1
  x <- goodVals[fillIdx]
  if (!is.null(lvls)) {
    x <- factor(x, levels=seq_along(lvls), labels=lvls)
  }
  x
}

#introduce sex switch smart segmentation

# divide contigs to autosomes and gonosomes
# AUTOSOMES_NUM = c(1:22)
# AUTOSOMES = c(AUTOSOMES_NUM, as.vector(outer("chr", AUTOSOMES_NUM, paste, sep="")))
# GONOSOMES_NUM = c("X","Y")
# GONOSOMES = c(GONOSOMES_NUM, as.vector(outer("chr", GONOSOMES_NUM, paste, sep="")))



#recalculate loss and gain cut-offs for gonosomes based on given standard cut-offs
#a) when gonosomal counts are NOT doubled in males (only for GATK default analysis)
if(AUTOSOME_MODE == "A") {
    #clonal cut-off
    #------------------------------------------------
    #percentage of the expected
    MINLOSS_PERC = (1-(2^MINLOSS))/0.5
    MINLOSS_BIAL_PERC = (1-(2^MINLOSS_BIAL))/1
    MINGAIN_PERC = ((2^MINGAIN)-1)/0.5
    
    #cut-off for male sample denoised by male controls 
    MINLOSS_MvsPoNM_chrXY = log(1-(MINLOSS_PERC*1),2)
    if(MINLOSS_MvsPoNM_chrXY < -2) {MINLOSS_MvsPoNM_chrXY = -2}
    MINGAIN_MvsPoNM_chrXY = log(1+(MINGAIN_PERC*1),2)    
    
    #cut-off for male sample denoised by female controls (chrY is lost in GATK, but kept haploid in CNVkit)
    MINLOSS_MvsPoNF_chrXY = log(0.5-(MINLOSS_PERC*0.5),2)
    if(MINLOSS_MvsPoNF_chrXY < -2) {MINLOSS_MvsPoNF_chrXY = -2}
    MINGAIN_MvsPoNF_chrXY = log(0.5+(MINGAIN_PERC*0.5),2)
    
    #cut-off for female sample denoised by male controls
    MINLOSS_FvsPoNM_chrX = log(2-(MINLOSS_PERC*1),2)
    MINLOSS_BIAL_FvsPoNM_chrX = log(2-(MINLOSS_BIAL_PERC*2),2)
    if(MINLOSS_BIAL_FvsPoNM_chrX < -1) {MINLOSS_BIAL_FvsPoNM_chrX = -1}
    MINGAIN_FvsPoNM_chrX = log(2+(MINGAIN_PERC*1),2)
    #MINLOSS_FvsPoNM_chrY - cannot happen as there is no chrY expected in female, all is lost
    MINGAIN_FvsPoNM_chrY = log(0+(MINGAIN_PERC*1),2)
    #------------------------------------------------
    
    #sub-clonal cut-off
    #------------------------------------------------
    #percentage of the expected
    MINLOSS_SUB_PERC = (1-(2^MINLOSS_SUB))/0.5
    MINGAIN_SUB_PERC = ((2^MINGAIN_SUB)-1)/0.5
    
    #cut-off for male sample denoised by male controls 
    MINLOSS_SUB_MvsPoNM_chrXY = log(1-(MINLOSS_SUB_PERC*1),2)
    if(MINLOSS_SUB_MvsPoNM_chrXY < -2) {MINLOSS_SUB_MvsPoNM_chrXY = -2}
    MINGAIN_SUB_MvsPoNM_chrXY = log(1+(MINGAIN_SUB_PERC*1),2) 
    
    #cut-off for male sample denoised by female controls (chrY is lost in GATK, but kept haploid in CNVkit)
    MINLOSS_SUB_MvsPoNF_chrXY = log(0.5-(MINLOSS_SUB_PERC*0.5),2)
    if(MINLOSS_SUB_MvsPoNF_chrXY < -2) {MINLOSS_SUB_MvsPoNF_chrXY = -2}
    MINGAIN_SUB_MvsPoNF_chrXY = log(0.5+(MINGAIN_SUB_PERC*0.5),2)
    
    #cut-off for female sample denoised by male controls
    MINLOSS_SUB_FvsPoNM_chrX = log(2-(MINLOSS_SUB_PERC*1),2)
    MINGAIN_SUB_FvsPoNM_chrX = log(2+(MINGAIN_SUB_PERC*1),2)
    #MINLOSS_SUB_FvsPoNM_chrY - cannot happen as there is no chrY expected in female, all is lost
    MINGAIN_SUB_FvsPoNM_chrY = log(0+(MINGAIN_SUB_PERC*1),2)
    #------------------------------------------------
}


#b) when gonosomal counts ARE doubled in males in GATK /// and for CNVkit
if(AUTOSOME_MODE != "A") {
  #clonal cut-off
  #------------------------------------------------
  #percentage of the expected (same)
  MINLOSS_PERC = (1-(2^MINLOSS))/0.5
  MINLOSS_BIAL_PERC = (1-(2^MINLOSS_BIAL))/1
  MINGAIN_PERC = ((2^MINGAIN)-1)/0.5
  
  #cut-off for male sample denoised by male controls (same as "A")
  MINLOSS_MvsPoNM_chrXY = log(1-(MINLOSS_PERC*1),2)
  if(MINLOSS_MvsPoNM_chrXY < -2) {MINLOSS_MvsPoNM_chrXY = -2}
  MINGAIN_MvsPoNM_chrXY = log(1+(MINGAIN_PERC*1),2)    
  
  #cut-off for male sample denoised by female controls (chrY is lost in GATK, but kept haploid in CNVkit)
  MINLOSS_MvsPoNF_chrXY = log(1-(MINLOSS_PERC*1),2)
  if(MINLOSS_MvsPoNF_chrXY < -2) {MINLOSS_MvsPoNF_chrXY = -2}
  MINGAIN_MvsPoNF_chrXY = log(1+(MINGAIN_PERC*1),2)
  
  #cut-off for female sample denoised by male controls (as autosomes)
  MINLOSS_FvsPoNM_chrX = log(1-(MINLOSS_PERC*0.5),2)
  MINLOSS_BIAL_FvsPoNM_chrX = log(1-(MINLOSS_BIAL_PERC*1),2)
  if(MINLOSS_BIAL_FvsPoNM_chrX < -1) {MINLOSS_BIAL_FvsPoNM_chrX = -1}
  MINGAIN_FvsPoNM_chrX = log(1+(MINGAIN_PERC*0.5),2)
  #MINLOSS_FvsPoNM_chrY - cannot happen as there is no chrY expected in female, all is lost
  if(AUTOSOME_MODE != "CNVKIT") {MINGAIN_FvsPoNM_chrY = log(0+(MINGAIN_PERC*0.5),2)}
  #recalculate for CNVkit as chrY is alsways haploid
  if(AUTOSOME_MODE == "CNVKIT") {MINGAIN_FvsPoNM_chrY = log(0+(MINGAIN_PERC*1),2)}
  #------------------------------------------------

  #sub-clonal cut-off
  #------------------------------------------------
  #percentage of the expected
  MINLOSS_SUB_PERC = (1-(2^MINLOSS_SUB))/0.5
  MINGAIN_SUB_PERC = ((2^MINGAIN_SUB)-1)/0.5
  
  #cut-off for male sample denoised by male controls (same) 
  MINLOSS_SUB_MvsPoNM_chrXY = log(1-(MINLOSS_SUB_PERC*1),2)
  if(MINLOSS_SUB_MvsPoNM_chrXY < -2) {MINLOSS_SUB_MvsPoNM_chrXY = -2}
  MINGAIN_SUB_MvsPoNM_chrXY = log(1+(MINGAIN_SUB_PERC*1),2) 
  
  #cut-off for male sample denoised by female controls (chrY is lost in GATK, but kept haploid in CNVkit)
  MINLOSS_SUB_MvsPoNF_chrXY = log(1-(MINLOSS_SUB_PERC*1),2)
  if(MINLOSS_SUB_MvsPoNF_chrXY < -2) {MINLOSS_SUB_MvsPoNF_chrXY = -2}
  MINGAIN_SUB_MvsPoNF_chrXY = log(1+(MINGAIN_SUB_PERC*1),2)
  
  #cut-off for female sample denoised by male controls
  MINLOSS_SUB_FvsPoNM_chrX = log(1-(MINLOSS_SUB_PERC*0.5),2)
  MINGAIN_SUB_FvsPoNM_chrX = log(1+(MINGAIN_SUB_PERC*0.5),2)
  #MINLOSS_SUB_FvsPoNM_chrY - cannot happen as there is no chrY expected in female, all is lost
  if(AUTOSOME_MODE != "CNVKIT") {MINGAIN_SUB_FvsPoNM_chrY = log(0+(MINGAIN_SUB_PERC*0.5),2)}
  #recalculate for CNVkit as chrY is alsways haploid
  if(AUTOSOME_MODE == "CNVKIT") {MINGAIN_SUB_FvsPoNM_chrY = log(0+(MINGAIN_SUB_PERC*1),2)}
  #------------------------------------------------
}


#c) incorporated edited cut-offs into the original table
#this gives general, not sex-switch understanding results with general cut-offs, true for autosomes or chrX in females when PON-F
DENOIS_DATA_INPUT = DENOIS_DATA_INPUT %>%
  mutate(MINLOSS_EDIT = MINLOSS) %>%
  mutate(MINLOSS_SUB_EDIT = MINLOSS_SUB) %>%
  mutate(MINLOSS_BIAL_EDIT = MINLOSS_BIAL) %>%
  mutate(MINGAIN_EDIT = MINGAIN) %>%
  mutate(MINGAIN_SUB_EDIT = MINGAIN_SUB)


#if sample is male and controls are males, chrX and chrY cut-offs need to be edited
if(GETCASESEX == "M" & GETPONSEX == "M") {
  DENOIS_DATA_INPUT = DENOIS_DATA_INPUT %>%
    mutate(MINLOSS_EDIT = ifelse(CONTIG == "chrX" | CONTIG == "X" | CONTIG == "chrY" | CONTIG == "Y", MINLOSS_MvsPoNM_chrXY, MINLOSS_EDIT)) %>%
    mutate(MINLOSS_SUB_EDIT = ifelse(CONTIG == "chrX" | CONTIG == "X" | CONTIG == "chrY" | CONTIG == "Y", MINLOSS_SUB_MvsPoNM_chrXY,  MINLOSS_SUB_EDIT)) %>%
    # mutate(MINLOSS_BIAL_EDIT = MINLOSS_BIAL_EDIT) %>% #bi-all del not possible as only one allele present in males
    mutate(MINGAIN_EDIT = ifelse(CONTIG == "chrX" | CONTIG == "X" | CONTIG == "chrY" | CONTIG == "Y", MINGAIN_MvsPoNM_chrXY,  MINGAIN_EDIT)) %>%
    mutate(MINGAIN_SUB_EDIT = ifelse(CONTIG == "chrX" | CONTIG == "X" | CONTIG == "chrY" | CONTIG == "Y", MINGAIN_SUB_MvsPoNM_chrXY,   MINGAIN_SUB_EDIT))
}

#if sample is male and are controls females, chrX cut-offs need to be edited
#In GATK, chrY is completely taken out from analysis during denoising
#In CNVkit, chrY stays and is haploid as it is always haploid in CNVkit
if (GETCASESEX == "M" & (GETPONSEX == "F" | GETPONSEX == "mixed")) {
  DENOIS_DATA_INPUT = DENOIS_DATA_INPUT %>%
    mutate(MINLOSS_EDIT = ifelse(CONTIG == "chrX" | CONTIG == "X" | CONTIG == "chrY" | CONTIG == "Y", MINLOSS_MvsPoNF_chrXY, MINLOSS_EDIT)) %>%
    mutate(MINLOSS_SUB_EDIT = ifelse(CONTIG == "chrX" | CONTIG == "X" | CONTIG == "chrY" | CONTIG == "Y", MINLOSS_SUB_MvsPoNF_chrXY,  MINLOSS_SUB_EDIT)) %>%
    # mutate(MINLOSS_BIAL_EDIT = MINLOSS_BIAL_EDIT) %>% #bi-all del not possible as only one allele present in males
    mutate(MINGAIN_EDIT = ifelse(CONTIG == "chrX" | CONTIG == "X" | CONTIG == "chrY" | CONTIG == "Y", MINGAIN_MvsPoNF_chrXY,  MINGAIN_EDIT)) %>%
    mutate(MINGAIN_SUB_EDIT = ifelse(CONTIG == "chrX" | CONTIG == "X" | CONTIG == "chrY" | CONTIG == "Y", MINGAIN_SUB_MvsPoNF_chrXY,   MINGAIN_SUB_EDIT))
}

#if sample is female and controls are males, chrX and chrY cut-offs need to be edited
if(GETCASESEX == "F" & (GETPONSEX == "M" | GETPONSEX == "mixed")) {
  DENOIS_DATA_INPUT = DENOIS_DATA_INPUT %>%
    #chrX
    mutate(MINLOSS_EDIT = ifelse(CONTIG == "chrX" | CONTIG == "X", MINLOSS_FvsPoNM_chrX, MINLOSS_EDIT)) %>%
    mutate(MINLOSS_SUB_EDIT = ifelse(CONTIG == "chrX" | CONTIG == "X", MINLOSS_SUB_FvsPoNM_chrX,  MINLOSS_SUB_EDIT)) %>%
    mutate(MINLOSS_BIAL_EDIT = ifelse(CONTIG == "chrX" | CONTIG == "X", MINLOSS_BIAL_FvsPoNM_chrX, MINLOSS_BIAL_EDIT)) %>%
    mutate(MINGAIN_EDIT = ifelse(CONTIG == "chrX" | CONTIG == "X", MINGAIN_FvsPoNM_chrX, MINGAIN_EDIT)) %>%
    mutate(MINGAIN_SUB_EDIT = ifelse(CONTIG == "chrX" | CONTIG == "X", MINGAIN_SUB_FvsPoNM_chrX, MINGAIN_SUB_EDIT)) %>%
    #chrY
    mutate(MINLOSS_EDIT = ifelse(CONTIG == "chrY" | CONTIG == "Y", -Inf, MINLOSS_EDIT)) %>% #-Inf because whole chrY is just lost 
    mutate(MINLOSS_SUB_EDIT = ifelse(CONTIG == "chrY" | CONTIG == "Y", -Inf,  MINLOSS_SUB_EDIT)) %>% #-Inf because whole chrY is just lost 
    # mutate(MINLOSS_BIAL_EDIT = MINLOSS_BIAL_EDIT) %>% #whole chrY is lost, bial value can be still used for smoothing
    mutate(MINGAIN_EDIT = ifelse(CONTIG == "chrY" | CONTIG == "Y", MINGAIN_FvsPoNM_chrY, MINGAIN_EDIT)) %>%
    mutate(MINGAIN_SUB_EDIT = ifelse(CONTIG == "chrY" | CONTIG == "Y", MINGAIN_SUB_FvsPoNM_chrY, MINGAIN_SUB_EDIT))
}


  
  

message(paste0("   R ... ", date(), " - STAGE 2/5 - CN segmentation"))

#SEGMENTATION
#-------------------------------------------------------------------------------------------
#PHASE1 - divide chromosomes to parts with specific distance (GAP), filter out parts with less than n of probes (MINKEEP),
#then divide to segments that differ by DIFFERENCE and filter out segments with less than n of probes (MINKEEP)

PHASE1 = data.frame("NA",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
names(PHASE1) = c("CONTIG", "START", "END", "LOG2_COPY_RATIO", "MINLOSS_EDIT", "MINLOSS_SUB_EDIT", "MINLOSS_BIAL_EDIT", "MINGAIN_EDIT", "MINGAIN_SUB_EDIT",
                  "DISTANCE", "N1", "CHRPART", "CHRPARTN", "NofPinCHRPART", "DIF", "DIFFERENCE_EDIT", "N2", "SEG", 
                  # "LOG2_BEFORE", "LOG2_AFTER", "DIFA", "DIFB", "DIFAB", "DIFAB_TEST", "DIFAB_TEST_BEFORE", "DIFAB_TEST_AFTER",
                  "SEGMENT", "NofP", "NofP_BEFORE", "NofP_AFTER", "KEEPING",
                  "DIF_WAY_A", "DIF_WAY_B", "EXTR_A", "EXTR_B", "EXTR")

LISTofCHR = unique(DENOIS_DATA_INPUT$CONTIG) %>% as.list()

for (CHR in (LISTofCHR)) {
  
 #find segments of the chromosome that are futher than value GAP and exclude segments that are smaller than MINKEEP value
  DENOIS_DATA_PH1 = DENOIS_DATA_INPUT %>%
    filter(CONTIG == CHR) %>%
    arrange(START) %>%
    mutate(DISTANCE = START-rowShift(END, -1)) %>%
    mutate(DISTANCE = ifelse(is.na(DISTANCE), -1, DISTANCE))
  # DENOIS_DATA_PH1[c("DISTANCE")][is.na(DENOIS_DATA_PH1[c("DISTANCE")])] <- -1
  DENOIS_DATA_PH1 = DENOIS_DATA_PH1 %>%
    mutate(N1 = 1:n()) %>%
    mutate(CHRPART = ifelse(DISTANCE > GAP | DISTANCE == -1, "NEWPART", NA))
    
  SUBSET = DENOIS_DATA_PH1 %>%
    filter(CHRPART == "NEWPART") %>%
    mutate(CHRPARTN = 1:n()) %>%
    select(N1, CHRPARTN)
  
  DENOIS_DATA_PH1 = DENOIS_DATA_PH1 %>%
    mutate(CHRPARTN = vlookup(N1, dict = SUBSET, lookup_column = which(colnames(SUBSET)=="N1"), result_column = which(colnames(SUBSET)=="CHRPARTN"))) %>%
    mutate(CHRPARTN = fillNAgaps(CHRPARTN))
  
  FINDNofPinCHRPART = DENOIS_DATA_PH1 %>%
    group_by(CHRPARTN) %>%
    summarise(NofP = n())
  
  DENOIS_DATA_PH1 = DENOIS_DATA_PH1 %>%
    mutate(NofPinCHRPART = vlookup(CHRPARTN, dict = FINDNofPinCHRPART, lookup_column = which(colnames(FINDNofPinCHRPART)=="CHRPARTN"), result_column = which(colnames(FINDNofPinCHRPART)=="NofP"))) %>%
    filter(NofPinCHRPART >= MINKEEP)
  
  LISTofCHRPART = unique(DENOIS_DATA_PH1$CHRPARTN) %>% as.list()

  
      #loop though the chromosome parts and do segmentation
    
      for (PART in (LISTofCHRPART)) {

          DENOIS_DATA_PART_PH1 = DENOIS_DATA_PH1 %>%
            filter(CHRPARTN == PART) %>%
            mutate(DIF = abs(LOG2_COPY_RATIO-rowShift(LOG2_COPY_RATIO, -1))) %>%
            mutate(DIF = ifelse(is.na(DIF), -1, DIF))
          # DENOIS_DATA_PART_PH1[c("DIF")][is.na(DENOIS_DATA_PART_PH1[c("DIF")])] <- -1
          
          #NEW - define difference
          DENOIS_DATA_PART_PH1 = DENOIS_DATA_PART_PH1 %>%
            mutate(DIFFERENCE_EDIT = DIFFERENCE) %>%
            mutate(DIFFERENCE_EDIT = ifelse(LOG2_COPY_RATIO <= MINLOSS_SUB_EDIT & rowShift(LOG2_COPY_RATIO, -1) <= MINLOSS_SUB_EDIT, DIFFERENCE*DIFFERENCE_SMOOTH_COEF_SUB, DIFFERENCE_EDIT)) %>%
            mutate(DIFFERENCE_EDIT = ifelse(LOG2_COPY_RATIO >= MINGAIN_SUB_EDIT & rowShift(LOG2_COPY_RATIO, -1) >= MINGAIN_SUB_EDIT, DIFFERENCE*DIFFERENCE_SMOOTH_COEF_SUB, DIFFERENCE_EDIT)) %>%
            mutate(DIFFERENCE_EDIT = ifelse(LOG2_COPY_RATIO <= MINLOSS_EDIT & rowShift(LOG2_COPY_RATIO, -1) <= MINLOSS_EDIT, DIFFERENCE*DIFFERENCE_SMOOTH_COEF, DIFFERENCE_EDIT)) %>%
            mutate(DIFFERENCE_EDIT = ifelse(LOG2_COPY_RATIO >= MINGAIN_EDIT & rowShift(LOG2_COPY_RATIO, -1) >= MINGAIN_EDIT, DIFFERENCE*DIFFERENCE_SMOOTH_COEF, DIFFERENCE_EDIT)) %>%
            mutate(DIFFERENCE_EDIT = ifelse(is.na(DIFFERENCE_EDIT), DIFFERENCE, DIFFERENCE_EDIT))
            
          DENOIS_DATA_PART_PH1 = DENOIS_DATA_PART_PH1 %>%
            mutate(N2 = 1:n()) %>%
            mutate(SEG = ifelse(DIF > DIFFERENCE_EDIT | DIF == -1, "NEW", NA)) %>%
            
            # # LOG2_BEFORE, LOG2_AFTER, DIFA, DIFB, DIFAB, DIFAB_TEST, DIFAB_TEST_BEFORE, DIFAB_TEST_AFTER
            # 
            # #handling stairs-like abnormalities
            # mutate(LOG2_BEFORE = rowShift(LOG2_COPY_RATIO, -1)) %>%
            # mutate(LOG2_AFTER = rowShift(LOG2_COPY_RATIO, +1)) %>%
            # mutate(DIFA = LOG2_COPY_RATIO-LOG2_BEFORE) %>%
            # mutate(DIFB = LOG2_AFTER-LOG2_COPY_RATIO) %>%
            # mutate(DIFAB = (DIFA+DIFB)) %>%
            # 
            # #1) stairs up to GAIN
            # #(LOG2_BEFORE<MINGAIN_SUB_EDIT)&(LOG2_AFTER>=MINGAIN_SUB_EDIT) - crossing states NORMAL up to (at least) SUB-GAINS
            # #(DIFAB > 0) - log is going up
            # #(LOG2_COPY_RATIO >= 0) - GAIN start
            # #(abs(DIFAB)) - difference between "before" and "after" is more than segmentation difference
            # mutate(DIFAB_TEST = ifelse(  (LOG2_BEFORE<MINGAIN_SUB_EDIT)&(LOG2_AFTER>=MINGAIN_SUB_EDIT)  &  (DIFAB > 0)  &  (LOG2_COPY_RATIO >= 0)  &  (abs(DIFAB) > DIFFERENCE_EDIT)  &  (abs(DIFA) == abs(DIFB))   , "new_here",  NA)) %>%
            # mutate(DIFAB_TEST = ifelse(  (LOG2_BEFORE<MINGAIN_SUB_EDIT)&(LOG2_AFTER>=MINGAIN_SUB_EDIT)  &  (DIFAB > 0)  &  (LOG2_COPY_RATIO >= 0)  &  (abs(DIFAB) > DIFFERENCE_EDIT)  &  (abs(DIFA) >  abs(DIFB))   , "new_here",  DIFAB_TEST)) %>%
            # mutate(DIFAB_TEST = ifelse(  (LOG2_BEFORE<MINGAIN_SUB_EDIT)&(LOG2_AFTER>=MINGAIN_SUB_EDIT)  &  (DIFAB > 0)  &  (LOG2_COPY_RATIO >= 0)  &  (abs(DIFAB) > DIFFERENCE_EDIT)  &  (abs(DIFA) <  abs(DIFB))   , "new_after", DIFAB_TEST)) %>%
            # 
            # #2) stairs up from LOSS
            # #(LOG2_BEFORE<=MINLOSS_SUB_EDIT)&(LOG2_AFTER>MINLOSS_SUB_EDIT) - crossing states (at least) SUB-LOSS to NORMAL
            # #(DIFAB > 0) - log is going up
            # #(LOG2_COPY_RATIO < 0) - LOSS end
            # #(abs(DIFAB)) - difference between "before" and "after" is more than segmentation difference
            # mutate(DIFAB_TEST = ifelse(  (LOG2_BEFORE<=MINLOSS_SUB_EDIT)&(LOG2_AFTER>MINLOSS_SUB_EDIT)  &  (DIFAB > 0)  &  (LOG2_COPY_RATIO < 0)  &  (abs(DIFAB) > DIFFERENCE_EDIT)  &  (abs(DIFA) == abs(DIFB))   , "new_after", DIFAB_TEST)) %>%
            # mutate(DIFAB_TEST = ifelse(  (LOG2_BEFORE<=MINLOSS_SUB_EDIT)&(LOG2_AFTER>MINLOSS_SUB_EDIT)  &  (DIFAB > 0)  &  (LOG2_COPY_RATIO < 0)  &  (abs(DIFAB) > DIFFERENCE_EDIT)  &  (abs(DIFA) > abs(DIFB))    , "new_here",  DIFAB_TEST)) %>%
            # mutate(DIFAB_TEST = ifelse(  (LOG2_BEFORE<=MINLOSS_SUB_EDIT)&(LOG2_AFTER>MINLOSS_SUB_EDIT)  &  (DIFAB > 0)  &  (LOG2_COPY_RATIO < 0)  &  (abs(DIFAB) > DIFFERENCE_EDIT)  &  (abs(DIFA) < abs(DIFB))    , "new_after", DIFAB_TEST)) %>%
            # 
            # #3) stairs down from GAIN
            # #(LOG2_BEFORE>=MINGAIN_SUB_EDIT)&(LOG2_AFTER<MINGAIN_SUB_EDIT) - crossing states (at least) SUB-GAIN to NORMAL
            # #(DIFAB < 0) - log is going down
            # #(LOG2_COPY_RATIO >= 0) - GAIN end
            # #(abs(DIFAB)) - difference between "before" and "after" is more than segmentation difference
            # mutate(DIFAB_TEST = ifelse(  (LOG2_BEFORE>=MINGAIN_SUB_EDIT)&(LOG2_AFTER<MINGAIN_SUB_EDIT)  &  (DIFAB < 0)  &  (LOG2_COPY_RATIO >= 0)  &  (abs(DIFAB) > DIFFERENCE_EDIT)  &  (abs(DIFA) == abs(DIFB))   , "new_after", DIFAB_TEST)) %>%
            # mutate(DIFAB_TEST = ifelse(  (LOG2_BEFORE>=MINGAIN_SUB_EDIT)&(LOG2_AFTER<MINGAIN_SUB_EDIT)  &  (DIFAB < 0)  &  (LOG2_COPY_RATIO >= 0)  &  (abs(DIFAB) > DIFFERENCE_EDIT)  &  (abs(DIFA) > abs(DIFB))    , "new_here",  DIFAB_TEST)) %>%
            # mutate(DIFAB_TEST = ifelse(  (LOG2_BEFORE>=MINGAIN_SUB_EDIT)&(LOG2_AFTER<MINGAIN_SUB_EDIT)  &  (DIFAB < 0)  &  (LOG2_COPY_RATIO >= 0)  &  (abs(DIFAB) > DIFFERENCE_EDIT)  &  (abs(DIFA) < abs(DIFB))    , "new_after", DIFAB_TEST)) %>%
            # 
            # #4) stairs down to LOSS
            # #(LOG2_BEFORE>MINLOSS_SUB_EDIT)&(LOG2_AFTER<=MINLOSS_SUB_EDIT) - crossing states NORMAL to (at least) SUB-LOSS
            # #(DIFAB < 0) - log is going down
            # #(LOG2_COPY_RATIO < 0) - LOSS start
            # #(abs(DIFAB)) - difference between "before" and "after" is more than segmentation difference
            # mutate(DIFAB_TEST = ifelse(  (LOG2_BEFORE>MINLOSS_SUB_EDIT)&(LOG2_AFTER<=MINLOSS_SUB_EDIT)  &  (DIFAB < 0)  &  (LOG2_COPY_RATIO < 0)  &  (abs(DIFAB) > DIFFERENCE_EDIT)  &  (abs(DIFA) == abs(DIFB))    , "new_here", DIFAB_TEST)) %>%
            # mutate(DIFAB_TEST = ifelse(  (LOG2_BEFORE>MINLOSS_SUB_EDIT)&(LOG2_AFTER<=MINLOSS_SUB_EDIT)  &  (DIFAB < 0)  &  (LOG2_COPY_RATIO < 0)  &  (abs(DIFAB) > DIFFERENCE_EDIT)  &  (abs(DIFA) > abs(DIFB))     , "new_here", DIFAB_TEST)) %>%
            # mutate(DIFAB_TEST = ifelse(  (LOG2_BEFORE>MINLOSS_SUB_EDIT)&(LOG2_AFTER<=MINLOSS_SUB_EDIT)  &  (DIFAB < 0)  &  (LOG2_COPY_RATIO < 0)  &  (abs(DIFAB) > DIFFERENCE_EDIT)  &  (abs(DIFA) < abs(DIFB))     , "new_after", DIFAB_TEST)) %>%
            # 
            # mutate(DIFAB_TEST_BEFORE = rowShift(DIFAB_TEST, -1)) %>%
            # mutate(DIFAB_TEST_AFTER = rowShift(DIFAB_TEST, +1)) %>%
            # 
            # mutate(SEG = case_when(DIFAB_TEST == "new_here" ~ "NEW", TRUE ~ as.character(SEG))) %>%
            # mutate(SEG = case_when(DIFAB_TEST_BEFORE == "new_after" ~ "NEW", TRUE ~ as.character(SEG))) %>%
            # mutate(SEG = case_when(DIFAB_TEST_AFTER == "new_before" ~ "NEW", TRUE ~ as.character(SEG))) %>%

            mutate(SEG = ifelse((N2 > 1) & (LOG2_COPY_RATIO <= MINLOSS_BIAL_EDIT) & (rowShift(LOG2_COPY_RATIO, -1) <= MINLOSS_BIAL_EDIT), NA, SEG )) #bi-allelic deletions to be merge


          SUBSET = DENOIS_DATA_PART_PH1 %>%
            filter(SEG == "NEW") %>%
            mutate(SEGMENT = 1:n()) %>%
            select(N2, SEGMENT)
          
          DENOIS_DATA_PART_PH1 = DENOIS_DATA_PART_PH1 %>%
            mutate(SEGMENT = vlookup(N2, dict = SUBSET, lookup_column = which(colnames(SUBSET)=="N2"), result_column = which(colnames(SUBSET)=="SEGMENT"))) %>%
            mutate(SEGMENT = fillNAgaps(SEGMENT))
        
          SKIP = DENOIS_DATA_PART_PH1 %>%
            group_by(SEGMENT) %>%
            summarise(NofP = n())
          
          DENOIS_DATA_PART_PH1 = DENOIS_DATA_PART_PH1 %>%
            mutate(NofP = vlookup(SEGMENT, dict = SKIP, lookup_column = which(colnames(SKIP)=="SEGMENT"), result_column = which(colnames(SKIP)=="NofP"))) %>%
            mutate(NofP_BEFORE = rowShift(NofP, -1)) %>%
            mutate(NofP_AFTER = rowShift(NofP, +1)) %>%
            mutate(KEEPING = ifelse(NofP >= MINKEEP, "keep", ".")) %>%
            mutate(KEEPING = ifelse(NofP < MINKEEP & (NofP_BEFORE < MINKEEP | NofP_AFTER < MINKEEP), "keep", KEEPING)) %>% #keep if several 1-probes segments next to each other
            
            mutate(DIF_WAY_A = LOG2_COPY_RATIO-(rowShift(LOG2_COPY_RATIO, -1))) %>%
            mutate(DIF_WAY_B = LOG2_COPY_RATIO-(rowShift(LOG2_COPY_RATIO, +1))) %>%
            mutate(EXTR_A = ifelse(DIF_WAY_A <= (MINLOSS_EDIT*2) | DIF_WAY_A >= (MINGAIN_EDIT*2), "EXTREME", ".")) %>%
            mutate(EXTR_B = ifelse(DIF_WAY_B <= (MINLOSS_EDIT*2) | DIF_WAY_B >= (MINGAIN_EDIT*2), "EXTREME", ".")) %>%
            mutate(EXTR = ifelse(EXTR_A == "EXTREME" & EXTR_B == "EXTREME", "EXTREME", ".")) %>%
            mutate(KEEPING = ifelse(NofP < MINKEEP & EXTR == "EXTREME", ".", KEEPING)) #but exclude those 1-probes in extremes as they influence mean of segments too much
            
          

          PHASE1 = rbind(PHASE1, DENOIS_DATA_PART_PH1)
          
          }}
          
          PHASE1 = filter(PHASE1, CONTIG != "NA") %>%
                   # filter(NofP >= MINKEEP) %>%
                   filter(KEEPING == "keep")


                 #trial 2
                   # mutate(DIFaround = ifelse(NofP < MINKEEP, abs(rowShift(LOG2_COPY_RATIO, -1)-rowShift(LOG2_COPY_RATIO, +1)), 0)) %>%
                   # mutate(DIFaround = ifelse(NofP < MINKEEP & is.na(DIFaround), 0, DIFaround)) %>%
                   # 
                   # mutate(DIFFERENCE_EDIT2 = ifelse(NofP < MINKEEP, DIFFERENCE, 0)) %>%
                   # mutate(DIFFERENCE_EDIT2 = ifelse(NofP < MINKEEP & rowShift(LOG2_COPY_RATIO, -1) <= MINLOSS & rowShift(LOG2_COPY_RATIO, +1) <= MINLOSS, DIFFERENCE*2, DIFFERENCE_EDIT2)) %>%
                   # mutate(DIFFERENCE_EDIT2 = ifelse(NofP < MINKEEP & rowShift(LOG2_COPY_RATIO, -1) >= MINGAIN & rowShift(LOG2_COPY_RATIO, +1) >= MINGAIN, DIFFERENCE*2, DIFFERENCE_EDIT2)) %>%
                   # mutate(DIFFERENCE_EDIT2 = ifelse(NofP < MINKEEP & is.na(DIFFERENCE_EDIT2), 0, DIFFERENCE_EDIT2)) %>%
                   # 
                   # mutate(FILTERING = ifelse(DIFaround != 0 & DIFaround <= DIFFERENCE_EDIT2, "FILTER_OUT", ".")) #%>%
                   # filter(FILTERING != "FILTER_OUT")
                  #trial 1
                   # mutate(FILTERING = ifelse(NofP < MINKEEP & (rowShift(GENOTYPE_RAW, -1) == rowShift(GENOTYPE_RAW, 1)), "FILTER_OUT", "." ))
          
#-------------------------------------------------------------------------------------------
      
#PHASE2 - for each segment that came from phase 1, create segmentation based on DIFFERENCE and smooth data

PHASE2 = data.frame(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
names(PHASE2)= c("CHRPARTN", "CONTIG", "START", "END", "MEAN_L2R", "MEAN_L2R_EDIT", "SIZE", "PROBES", "TYPE", "MINLOSS_EDIT", "MINLOSS_SUB_EDIT", "MINLOSS_BIAL_EDIT", "MINGAIN_EDIT", "MINGAIN_SUB_EDIT")
      
LISTofCHR = unique(PHASE1$CONTIG) %>% as.list()
      
      
    for (CHR in (LISTofCHR)) { 
      
      DENOIS_DATA_PH2 = PHASE1 %>% 
          filter(CONTIG == CHR)
      LISTofCHRPART = unique(DENOIS_DATA_PH2$CHRPARTN) %>% as.list()
      
      
      
        for (PART in (LISTofCHRPART)) {
          
          DENOIS_DATA_PART_PH2 = DENOIS_DATA_PH2 %>%
              filter(CHRPARTN == PART) %>%
              mutate(DIF2 = abs(LOG2_COPY_RATIO-rowShift(LOG2_COPY_RATIO, -1))) %>%
              mutate(DIF2 = ifelse(is.na(DIF2), -1, DIF2))
              # DENOIS_DATA_PART_PH2[c("DIF2")][is.na(DENOIS_DATA_PART_PH2[c("DIF2")])] <- -1

          #NEW - define difference
          DENOIS_DATA_PART_PH2 = DENOIS_DATA_PART_PH2 %>%
            mutate(DIFFERENCE_EDIT2 = DIFFERENCE) %>%
            mutate(DIFFERENCE_EDIT2 = ifelse(LOG2_COPY_RATIO <= MINLOSS_SUB_EDIT & rowShift(LOG2_COPY_RATIO, -1) <= MINLOSS_SUB_EDIT, DIFFERENCE*DIFFERENCE_SMOOTH_COEF_SUB, DIFFERENCE_EDIT2)) %>%
            mutate(DIFFERENCE_EDIT2 = ifelse(LOG2_COPY_RATIO >= MINGAIN_SUB_EDIT & rowShift(LOG2_COPY_RATIO, -1) >= MINGAIN_SUB_EDIT, DIFFERENCE*DIFFERENCE_SMOOTH_COEF_SUB, DIFFERENCE_EDIT2)) %>%
            mutate(DIFFERENCE_EDIT2 = ifelse(LOG2_COPY_RATIO <= MINLOSS_EDIT & rowShift(LOG2_COPY_RATIO, -1) <= MINLOSS_EDIT, DIFFERENCE*DIFFERENCE_SMOOTH_COEF, DIFFERENCE_EDIT2)) %>%
            mutate(DIFFERENCE_EDIT2 = ifelse(LOG2_COPY_RATIO >= MINGAIN_EDIT & rowShift(LOG2_COPY_RATIO, -1) >= MINGAIN_EDIT, DIFFERENCE*DIFFERENCE_SMOOTH_COEF, DIFFERENCE_EDIT2)) %>%
            mutate(DIFFERENCE_EDIT2 = ifelse(is.na(DIFFERENCE_EDIT2), DIFFERENCE, DIFFERENCE_EDIT2))
          
          DENOIS_DATA_PART_PH2 = DENOIS_DATA_PART_PH2 %>%
              mutate(N3 = 1:n()) %>%
              mutate(SEG2 = ifelse(DIF2 > DIFFERENCE_EDIT2 | DIF2 == -1, "NEW", NA)) %>%
            
            # LOG2_BEFORE, LOG2_AFTER, DIFA, DIFB, DIFAB, DIFAB_TEST, DIFAB_TEST_BEFORE, DIFAB_TEST_AFTER
            # LOG2_BEFORE2, LOG2_AFTER2, DIF2A, DIF2B, DIF2AB, DIF2AB_TEST, DIF2AB_TEST_BEFORE, DIF2AB_TEST_AFTER
            
            #handling stairs-like abnormalities
            mutate(LOG2_BEFORE2 = rowShift(LOG2_COPY_RATIO, -1)) %>%
            mutate(LOG2_AFTER2 = rowShift(LOG2_COPY_RATIO, +1)) %>%
            mutate(DIF2A = LOG2_COPY_RATIO-LOG2_BEFORE2) %>%
            mutate(DIF2B = LOG2_AFTER2-LOG2_COPY_RATIO) %>%
            mutate(DIF2AB = (DIF2A+DIF2B)) %>%
            
            #1) stairs up to GAIN
            #(LOG2_BEFORE2<MINGAIN_SUB_EDIT)&(LOG2_AFTER2>=MINGAIN_SUB_EDIT) - crossing states NORMAL up to (at least) SUB-GAINS
            #(DIF2AB > 0) - log is going up
            #(LOG2_COPY_RATIO >= 0) - GAIN start
            #(abs(DIF2AB)) - difference between "before" and "after" is more than segmentation difference
            mutate(DIF2AB_TEST = ifelse(  (LOG2_BEFORE2<MINGAIN_SUB_EDIT)&(LOG2_AFTER2>=MINGAIN_SUB_EDIT)  &  (DIF2AB > 0)  &  (LOG2_COPY_RATIO >= 0)  &  (abs(DIF2AB) > DIFFERENCE_EDIT2)  &  (abs(DIF2A) == abs(DIF2B))   , "new_here",  NA)) %>%
            mutate(DIF2AB_TEST = ifelse(  (LOG2_BEFORE2<MINGAIN_SUB_EDIT)&(LOG2_AFTER2>=MINGAIN_SUB_EDIT)  &  (DIF2AB > 0)  &  (LOG2_COPY_RATIO >= 0)  &  (abs(DIF2AB) > DIFFERENCE_EDIT2)  &  (abs(DIF2A) >  abs(DIF2B))   , "new_here",  DIF2AB_TEST)) %>%
            mutate(DIF2AB_TEST = ifelse(  (LOG2_BEFORE2<MINGAIN_SUB_EDIT)&(LOG2_AFTER2>=MINGAIN_SUB_EDIT)  &  (DIF2AB > 0)  &  (LOG2_COPY_RATIO >= 0)  &  (abs(DIF2AB) > DIFFERENCE_EDIT2)  &  (abs(DIF2A) <  abs(DIF2B))   , "new_after", DIF2AB_TEST)) %>%
            
            #2) stairs up from LOSS
            #(LOG2_BEFORE2<=MINLOSS_SUB_EDIT)&(LOG2_AFTER2>MINLOSS_SUB_EDIT) - crossing states (at least) SUB-LOSS to NORMAL
            #(DIF2AB > 0) - log is going up
            #(LOG2_COPY_RATIO < 0) - LOSS end
            #(abs(DIF2AB)) - difference between "before" and "after" is more than segmentation difference
            mutate(DIF2AB_TEST = ifelse(  (LOG2_BEFORE2<=MINLOSS_SUB_EDIT)&(LOG2_AFTER2>MINLOSS_SUB_EDIT)  &  (DIF2AB > 0)  &  (LOG2_COPY_RATIO < 0)  &  (abs(DIF2AB) > DIFFERENCE_EDIT2)  &  (abs(DIF2A) == abs(DIF2B))   , "new_after", DIF2AB_TEST)) %>%
            mutate(DIF2AB_TEST = ifelse(  (LOG2_BEFORE2<=MINLOSS_SUB_EDIT)&(LOG2_AFTER2>MINLOSS_SUB_EDIT)  &  (DIF2AB > 0)  &  (LOG2_COPY_RATIO < 0)  &  (abs(DIF2AB) > DIFFERENCE_EDIT2)  &  (abs(DIF2A) > abs(DIF2B))    , "new_here",  DIF2AB_TEST)) %>%
            mutate(DIF2AB_TEST = ifelse(  (LOG2_BEFORE2<=MINLOSS_SUB_EDIT)&(LOG2_AFTER2>MINLOSS_SUB_EDIT)  &  (DIF2AB > 0)  &  (LOG2_COPY_RATIO < 0)  &  (abs(DIF2AB) > DIFFERENCE_EDIT2)  &  (abs(DIF2A) < abs(DIF2B))    , "new_after", DIF2AB_TEST)) %>%
            
            #3) stairs down from GAIN
            #(LOG2_BEFORE2>=MINGAIN_SUB_EDIT)&(LOG2_AFTER2<MINGAIN_SUB_EDIT) - crossing states (at least) SUB-GAIN to NORMAL
            #(DIF2AB < 0) - log is going down
            #(LOG2_COPY_RATIO >= 0) - GAIN end
            #(abs(DIF2AB)) - difference between "before" and "after" is more than segmentation difference
            mutate(DIF2AB_TEST = ifelse(  (LOG2_BEFORE2>=MINGAIN_SUB_EDIT)&(LOG2_AFTER2<MINGAIN_SUB_EDIT)  &  (DIF2AB < 0)  &  (LOG2_COPY_RATIO >= 0)  &  (abs(DIF2AB) > DIFFERENCE_EDIT2)  &  (abs(DIF2A) == abs(DIF2B))   , "new_after", DIF2AB_TEST)) %>%
            mutate(DIF2AB_TEST = ifelse(  (LOG2_BEFORE2>=MINGAIN_SUB_EDIT)&(LOG2_AFTER2<MINGAIN_SUB_EDIT)  &  (DIF2AB < 0)  &  (LOG2_COPY_RATIO >= 0)  &  (abs(DIF2AB) > DIFFERENCE_EDIT2)  &  (abs(DIF2A) > abs(DIF2B))    , "new_here",  DIF2AB_TEST)) %>%
            mutate(DIF2AB_TEST = ifelse(  (LOG2_BEFORE2>=MINGAIN_SUB_EDIT)&(LOG2_AFTER2<MINGAIN_SUB_EDIT)  &  (DIF2AB < 0)  &  (LOG2_COPY_RATIO >= 0)  &  (abs(DIF2AB) > DIFFERENCE_EDIT2)  &  (abs(DIF2A) < abs(DIF2B))    , "new_after", DIF2AB_TEST)) %>%
            
            #4) stairs down to LOSS
            #(LOG2_BEFORE2>MINLOSS_SUB_EDIT)&(LOG2_AFTER2<=MINLOSS_SUB_EDIT) - crossing states NORMAL to (at least) SUB-LOSS
            #(DIF2AB < 0) - log is going down
            #(LOG2_COPY_RATIO < 0) - LOSS start
            #(abs(DIF2AB)) - difference between "before" and "after" is more than segmentation difference
            mutate(DIF2AB_TEST = ifelse(  (LOG2_BEFORE2>MINLOSS_SUB_EDIT)&(LOG2_AFTER2<=MINLOSS_SUB_EDIT)  &  (DIF2AB < 0)  &  (LOG2_COPY_RATIO < 0)  &  (abs(DIF2AB) > DIFFERENCE_EDIT2)  &  (abs(DIF2A) == abs(DIF2B))    , "new_here", DIF2AB_TEST)) %>%
            mutate(DIF2AB_TEST = ifelse(  (LOG2_BEFORE2>MINLOSS_SUB_EDIT)&(LOG2_AFTER2<=MINLOSS_SUB_EDIT)  &  (DIF2AB < 0)  &  (LOG2_COPY_RATIO < 0)  &  (abs(DIF2AB) > DIFFERENCE_EDIT2)  &  (abs(DIF2A) > abs(DIF2B))     , "new_here", DIF2AB_TEST)) %>%
            mutate(DIF2AB_TEST = ifelse(  (LOG2_BEFORE2>MINLOSS_SUB_EDIT)&(LOG2_AFTER2<=MINLOSS_SUB_EDIT)  &  (DIF2AB < 0)  &  (LOG2_COPY_RATIO < 0)  &  (abs(DIF2AB) > DIFFERENCE_EDIT2)  &  (abs(DIF2A) < abs(DIF2B))     , "new_after", DIF2AB_TEST)) %>%
            
            mutate(DIF2AB_TEST_BEFORE = rowShift(DIF2AB_TEST, -1)) %>%
            mutate(DIF2AB_TEST_AFTER = rowShift(DIF2AB_TEST, +1)) %>%

            mutate(SEG2 = case_when(DIF2AB_TEST == "new_here" ~ "NEW", TRUE ~ as.character(SEG2))) %>%
            mutate(SEG2 = case_when(DIF2AB_TEST_BEFORE == "new_after" ~ "NEW", TRUE ~ as.character(SEG2))) %>%
            mutate(SEG2 = case_when(DIF2AB_TEST_AFTER == "new_before" ~ "NEW", TRUE ~ as.character(SEG2))) %>%
            
            mutate(SEG2 = ifelse((N3 > 1) & (LOG2_COPY_RATIO <= MINLOSS_BIAL_EDIT) & (rowShift(LOG2_COPY_RATIO, -1) <= MINLOSS_BIAL_EDIT), NA, SEG2 )) #bi-allelic deletions to be merge
      
          SUBSET = DENOIS_DATA_PART_PH2 %>%
            filter(SEG2 == "NEW") %>%
            mutate(SEGMENT2 = 1:n()) %>%
            select(N3, SEGMENT2)
          
          DENOIS_DATA_PART_PH2 = DENOIS_DATA_PART_PH2 %>%
            mutate(SEGMENT2 = vlookup(N3, dict = SUBSET, lookup_column = which(colnames(SUBSET)=="N3"), result_column = which(colnames(SUBSET)=="SEGMENT2"))) %>%
            mutate(SEGMENT2 = fillNAgaps(SEGMENT2))

          SEGMENTS_1 = DENOIS_DATA_PART_PH2 %>%
            group_by(SEGMENT2) %>%
            summarise(MEAN_L2R = mean(LOG2_COPY_RATIO), START = min(START), END = max(END), PROBES = n())

          #Smooth segments
          SEGMENTS_1 = SEGMENTS_1 %>%
            mutate(SEGDIF = abs(MEAN_L2R-rowShift(MEAN_L2R, -1))) %>%
            mutate(SEGDIF = ifelse(is.na(SEGDIF), -1, SEGDIF))
          # SEGMENTS_1[c("SEGDIF")][is.na(SEGMENTS_1[c("SEGDIF")])] <- -1
          SEGMENTS_1 = SEGMENTS_1 %>%
            mutate(N4 = 1:n()) %>%
            mutate(SEG3 = ifelse(SEGDIF > (DIFFERENCE*SMOOTHPERC) | SEGDIF == -1, "NEW", NA))
          
          #find those created by "stair" algorithm that become less than MINKEEP and join them to the closest segment by log2r difference
          SEGMENTS_1 = SEGMENTS_1 %>%
            mutate(UNDER_MINKEEP = ifelse(PROBES < MINKEEP, "HERE", NA)) %>%
            # mutate(UNDER_MINKEEP_PREV = rowShift(UNDER_MINKEEP, -1)) %>%
            mutate(SEGDIFA = abs(MEAN_L2R-rowShift(MEAN_L2R, -1))) %>%
            mutate(SEGDIFB = abs(rowShift(MEAN_L2R, +1)-MEAN_L2R)) %>%
            mutate(SEGDIFA = ifelse(is.na(SEGDIFA), Inf, SEGDIFA)) %>%
            mutate(SEGDIFB = ifelse(is.na(SEGDIFB), Inf, SEGDIFB)) %>%
            mutate(TEST = ifelse(UNDER_MINKEEP == "HERE" & SEGDIFA < SEGDIFB, "join_previous", NA)) %>%
            mutate(TEST = ifelse(UNDER_MINKEEP == "HERE" & SEGDIFA > SEGDIFB, "join_following", TEST)) %>%
            
            
            mutate(TEST_BEFORE = rowShift(TEST, -1)) %>%
            mutate(TEST_AFTER = rowShift(TEST, +1)) %>%
           
            mutate(SEG3 = case_when(TEST == "join_previous" ~ NA, TRUE ~ as.character(SEG3))) %>%
            mutate(SEG3 = case_when(TEST_BEFORE == "join_following" ~ NA, TRUE ~ as.character(SEG3))) %>%
            
            #add the probe to n of probes to avoid segments with 3 probes to look like 2 probes in the next step
            mutate(PROBES2 = PROBES) %>%
            mutate(PROBES2 = case_when(TEST_BEFORE == "join_following" ~ (PROBES2+1), TRUE ~ PROBES2)) %>%
            mutate(PROBES2 = case_when(TEST_AFTER == "join_previous" ~ (PROBES2+1), TRUE ~ PROBES2))
            

          #find those that became MINKEEP size of probes as they can also often imitate stairs "up and down"
          SEGMENTS_1 = SEGMENTS_1 %>%
            mutate(SEGDIFC = abs(rowShift(MEAN_L2R, +1)-rowShift(MEAN_L2R, -1))) %>%
            mutate(JUST_MINKEEP = ifelse(PROBES2 == MINKEEP, "HERE", NA)) %>% #small segment that is minimaly kept n of probes
            mutate(JUST_MINKEEP2 = ifelse(   SEGDIFA < (DIFFERENCE*(1+(1-SMOOTHPERC)))  |  SEGDIFB < (DIFFERENCE*(1+(1-SMOOTHPERC))),  "HERE", NA )) %>% #at least one log2 difference is below cut-off (a bit more relaxed DIFFERENCE)
            mutate(JUST_MINKEEP3 = ifelse(SEGDIFC < DIFFERENCE*SMOOTHPERC, "HERE", NA)) %>% #previous and next segments are able to be smoothed (~not stair patern where I want high resolution)
            mutate(TEST2 = ifelse(JUST_MINKEEP == "HERE" & JUST_MINKEEP2 == "HERE" & JUST_MINKEEP3 == "HERE" & SEGDIFA < SEGDIFB, "join_previous", NA)) %>% #closer to previous one
            mutate(TEST2 = ifelse(JUST_MINKEEP == "HERE" & JUST_MINKEEP2 == "HERE" & JUST_MINKEEP3 == "HERE" & SEGDIFA > SEGDIFB, "join_following", TEST2)) %>% #closer to next one
            mutate(TEST_BEFORE2 = rowShift(TEST2, -1)) %>%

            mutate(SEG3 = case_when(TEST2 == "join_previous" ~ NA, TRUE ~ as.character(SEG3))) %>%
            mutate(SEG3 = case_when(TEST_BEFORE2 == "join_following" ~ NA, TRUE ~ as.character(SEG3)))
            
 
          SUBSET = SEGMENTS_1 %>%
            filter(SEG3 == "NEW") %>%
            mutate(SEGMENT3 = 1:n()) %>%
            select(N4, SEGMENT3)
          
          SEGMENTS_1 = SEGMENTS_1 %>%
            mutate(SEGMENT3 = vlookup(N4, dict = SUBSET, lookup_column = which(colnames(SUBSET)=="N4"), result_column = which(colnames(SUBSET)=="SEGMENT3"))) %>%
            mutate(SEGMENT3 = fillNAgaps(SEGMENT3))
          
          #Add new segments back to big table - to calculate precize L2R
          DATA_BACK = PHASE1 %>%
            filter(CONTIG == CHR) %>%
            filter(CHRPARTN == PART) %>%
            mutate(SEGMENT3 = vlookup(START, dict = SEGMENTS_1, lookup_column = which(colnames(SEGMENTS_1)=="START"), result_column = which(colnames(SEGMENTS_1)=="SEGMENT3"))) %>% #3vs9
            mutate(SEGMENT3 = fillNAgaps(SEGMENT3))
          
          SEGMENTS_2 = DATA_BACK %>%
            group_by(SEGMENT3) %>%
            summarise(MINLOSS_EDIT=mean(MINLOSS_EDIT), MINLOSS_SUB_EDIT=mean(MINLOSS_SUB_EDIT), MINLOSS_BIAL_EDIT=mean(MINLOSS_BIAL_EDIT), MINGAIN_EDIT=mean(MINGAIN_EDIT), MINGAIN_SUB_EDIT=mean(MINGAIN_SUB_EDIT), 
                      MEAN_L2R = mean(LOG2_COPY_RATIO), START = min(START), END = max(END), PROBES = n()) %>%
            mutate(CONTIG = CHR) %>%
            mutate(CHRPARTN = PART) %>%
            mutate(SIZE = END-START+1) %>%
            mutate(MEAN_L2R_EDIT = ifelse(MEAN_L2R < -2.25, -2.25, MEAN_L2R)) %>%
            mutate(MEAN_L2R_EDIT = ifelse(MEAN_L2R > 2.25, 2.25, MEAN_L2R_EDIT)) %>%
            
            mutate(TYPE = ifelse(MEAN_L2R <= MINLOSS_SUB_EDIT, "subLOSS", "NORMAL")) %>%
            mutate(TYPE = ifelse(MEAN_L2R >= MINGAIN_SUB_EDIT, "subGAIN", TYPE)) %>%
            mutate(TYPE = ifelse(MEAN_L2R <= MINLOSS_EDIT, "LOSS", TYPE)) %>%
            mutate(TYPE = ifelse(MEAN_L2R >= MINGAIN_EDIT, "GAIN", TYPE)) %>%
            
            mutate(FILTERPROBE = 0) %>%
            mutate(FILTERPROBE = ifelse((TYPE == "LOSS" | TYPE == "GAIN") & PROBES >= MINPROBE, 1, FILTERPROBE)) %>%
            mutate(FILTERPROBE = ifelse((TYPE == "subLOSS" | TYPE == "subGAIN") & PROBES >= MINPROBE_SUB, 1, FILTERPROBE)) %>%
            mutate(FILTERSIZE = 0) %>%
            mutate(FILTERSIZE = ifelse((TYPE == "LOSS" | TYPE == "GAIN") & SIZE >= MINSIZE, 1, FILTERSIZE)) %>%
            mutate(FILTERSIZE = ifelse((TYPE == "subLOSS" | TYPE == "subGAIN") & SIZE >= MINSIZE_SUB, 1, FILTERSIZE)) %>%
            
            # mutate(FILTER = ifelse((FILTERPROBE+FILTERSIZE) == 2 | TYPE == "NORMAL", "PASS", ".")) %>%
            mutate(FILTER = ifelse(TYPE == "NORMAL", "na", ".")) %>%
            mutate(FILTER = ifelse(TYPE != "NORMAL" & (FILTERPROBE+FILTERSIZE) == 2, "PASS", FILTER)) %>%
            mutate(FILTER = ifelse(TYPE != "NORMAL" & (FILTERPROBE+FILTERSIZE) < 2, "cnv_small_size", FILTER)) %>%
            
            select(CHRPARTN, CONTIG, START, END, MEAN_L2R, MEAN_L2R_EDIT, SIZE, PROBES, TYPE, MINLOSS_EDIT, MINLOSS_SUB_EDIT, MINLOSS_BIAL_EDIT, MINGAIN_EDIT, MINGAIN_SUB_EDIT)
          
          PHASE2 = rbind(PHASE2, SEGMENTS_2)
    
          }}

          PHASE2 = filter(PHASE2, CONTIG != "NA")
          
#-------------------------------------------------------------------------------------------

#PHASE3 - for each segment that came from phase 2, smooth data even more & merge bi-allelic deletions
# now limited to two rounds of smoothing, but can be done until there is no more segments to be smoothed,
# but it doesn't give too much value and takes time. IF max smoothing, then while(TEST != 0) {
          
PHASE_INPUT = PHASE2
TEST = 1 #for the first loop
N = 0
          
      while(TEST != 0 & N<2) {
            
            N=N+1
            # message(N)
            PHASE_OUTPUT = data.frame(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
            names(PHASE_OUTPUT)=c("CHRPARTN", "SEGMENTATION_CONDITION", "PON_SEX_N", "GROUP", "CASE_ID", "CASE_TYPE", "CASE_SEX", "CONTIG", "START", "END", "MEAN_L2R", "MEAN_L2R_EDIT", "SIZE", "PROBES", "FILTER", "TYPE", "MINLOSS_EDIT", "MINLOSS_SUB_EDIT", "MINLOSS_BIAL_EDIT", "MINGAIN_EDIT", "MINGAIN_SUB_EDIT", "NofNAs")
            
            LISTofCHR = unique(PHASE_INPUT$CONTIG) %>% as.list()
            
            for (CHR in (LISTofCHR)) { 
              
              DENOIS_DATA_PH3 = PHASE_INPUT %>%
                filter(CONTIG == CHR)
              LISTofCHRPART = unique(DENOIS_DATA_PH3$CHRPARTN) %>% as.list()
              
              for (PART in (LISTofCHRPART)) {
                
                DENOIS_DATA_PART_PH3 = DENOIS_DATA_PH3 %>%
                  filter(CHRPARTN == PART) %>%
                  mutate(SEGDIF = abs(MEAN_L2R-rowShift(MEAN_L2R, -1))) %>%
                  mutate(SEGDIF = ifelse(is.na(SEGDIF), -1, SEGDIF))
                # DENOIS_DATA_PART_PH3[c("SEGDIF")][is.na(DENOIS_DATA_PART_PH3[c("SEGDIF")])] <- -1
                
                #NEW - define difference
                DENOIS_DATA_PART_PH3 = DENOIS_DATA_PART_PH3 %>%
                  mutate(DIFFERENCE_EDIT3 = DIFFERENCE) %>%
                  mutate(DIFFERENCE_EDIT3 = ifelse(MEAN_L2R <= MINLOSS_SUB_EDIT & rowShift(MEAN_L2R, -1) <= MINLOSS_SUB_EDIT, DIFFERENCE*DIFFERENCE_SMOOTH_COEF_SUB, DIFFERENCE_EDIT3)) %>%
                  mutate(DIFFERENCE_EDIT3 = ifelse(MEAN_L2R >= MINGAIN_SUB_EDIT & rowShift(MEAN_L2R, -1) >= MINGAIN_SUB_EDIT, DIFFERENCE*DIFFERENCE_SMOOTH_COEF_SUB, DIFFERENCE_EDIT3)) %>%       
                  mutate(DIFFERENCE_EDIT3 = ifelse(MEAN_L2R <= MINLOSS_EDIT & rowShift(MEAN_L2R, -1) <= MINLOSS_EDIT, DIFFERENCE*DIFFERENCE_SMOOTH_COEF, DIFFERENCE_EDIT3)) %>%
                  mutate(DIFFERENCE_EDIT3 = ifelse(MEAN_L2R >= MINGAIN_EDIT & rowShift(MEAN_L2R, -1) >= MINGAIN_EDIT, DIFFERENCE*DIFFERENCE_SMOOTH_COEF, DIFFERENCE_EDIT3)) %>%
                  mutate(DIFFERENCE_EDIT3 = ifelse(is.na(DIFFERENCE_EDIT3), DIFFERENCE, DIFFERENCE_EDIT3))
                
                DENOIS_DATA_PART_PH3 = DENOIS_DATA_PART_PH3 %>%
                  mutate(N5 = 1:n()) %>%
                  mutate(SEG4 = ifelse(SEGDIF > (DIFFERENCE_EDIT3*SMOOTHPERC) | SEGDIF == -1, "NEW", NA)) %>%
                  mutate(SEG4 = ifelse((N5 > 1) & (MEAN_L2R <= MINLOSS_BIAL_EDIT) & (rowShift(MEAN_L2R, -1) <= MINLOSS_BIAL_EDIT), NA, SEG4 )) #bi-allelic deletions to be merge
                
                NofNAsum = nrow(DENOIS_DATA_PART_PH3 %>% filter(is.na(SEG4)))
                
                SUBSET = DENOIS_DATA_PART_PH3 %>%
                  filter(SEG4 == "NEW") %>%
                  mutate(SEGMENT4 = 1:n()) %>%
                  select(N5, SEGMENT4)
                
                DENOIS_DATA_PART_PH3 = DENOIS_DATA_PART_PH3 %>%
                  mutate(SEGMENT4 = vlookup(N5, dict = SUBSET, lookup_column = which(colnames(SUBSET)=="N5"), result_column = which(colnames(SUBSET)=="SEGMENT4"))) %>%
                  mutate(SEGMENT4 = fillNAgaps(SEGMENT4))
                
                #Add new segments back to big table - to calculate precize L2R
                DATA_BACK_2 = PHASE1 %>%
                  filter(CONTIG == CHR) %>%
                  filter(CHRPARTN == PART) %>%
                  mutate(SEGMENT4 = vlookup(START, dict = DENOIS_DATA_PART_PH3, lookup_column = which(colnames(DENOIS_DATA_PART_PH3)=="START"), result_column = which(colnames(DENOIS_DATA_PART_PH3)=="SEGMENT4"))) %>% #3vs14
                  mutate(SEGMENT4 = fillNAgaps(SEGMENT4))
                
                SEGMENTS_3 = DATA_BACK_2 %>%
                  group_by(SEGMENT4) %>%
                  summarise(MINLOSS_EDIT=mean(MINLOSS_EDIT), MINLOSS_SUB_EDIT=mean(MINLOSS_SUB_EDIT), MINLOSS_BIAL_EDIT=mean(MINLOSS_BIAL_EDIT), MINGAIN_EDIT=mean(MINGAIN_EDIT), MINGAIN_SUB_EDIT=mean(MINGAIN_SUB_EDIT), 
                            MEAN_L2R = mean(LOG2_COPY_RATIO), START = min(START), END = max(END), PROBES = n()) %>%
                  mutate(SEGMENTATION_CONDITION = GETSEGMENTCONDITION) %>%
                  mutate(CASE_ID = GETCASEID) %>%
                  mutate(CASE_TYPE = GETCASETYPE) %>%
                  mutate(CASE_SEX = GETCASESEX) %>%
                  mutate(GROUP = GETGROUPID) %>%
                  mutate(PON_SEX_N = paste0(GETPONSEX, "_", NofCTRL_USED)) %>%
                  
                  mutate(CONTIG = CHR) %>%
                  mutate(CHRPARTN = PART) %>%
                  mutate(SIZE = END-START+1) %>%
                  mutate(MEAN_L2R_EDIT = ifelse(MEAN_L2R < -2.25, -2.25, MEAN_L2R)) %>%
                  mutate(MEAN_L2R_EDIT = ifelse(MEAN_L2R > 2.25, 2.25, MEAN_L2R_EDIT)) %>%
                  
                  mutate(TYPE = ifelse(MEAN_L2R <= MINLOSS_SUB_EDIT, "subLOSS", "NORMAL")) %>%
                  mutate(TYPE = ifelse(MEAN_L2R >= MINGAIN_SUB_EDIT, "subGAIN", TYPE)) %>%
                  mutate(TYPE = ifelse(MEAN_L2R <= MINLOSS_EDIT, "LOSS", TYPE)) %>%
                  mutate(TYPE = ifelse(MEAN_L2R >= MINGAIN_EDIT, "GAIN", TYPE)) %>%
                  
                  mutate(FILTERPROBE = 0) %>%
                  mutate(FILTERPROBE = ifelse((TYPE == "LOSS" | TYPE == "GAIN") & PROBES >= MINPROBE, 1, FILTERPROBE)) %>%
                  mutate(FILTERPROBE = ifelse((TYPE == "subLOSS" | TYPE == "subGAIN") & PROBES >= MINPROBE_SUB, 1, FILTERPROBE)) %>%
                  mutate(FILTERSIZE = 0) %>%
                  mutate(FILTERSIZE = ifelse((TYPE == "LOSS" | TYPE == "GAIN") & SIZE >= MINSIZE, 1, FILTERSIZE)) %>%
                  mutate(FILTERSIZE = ifelse((TYPE == "subLOSS" | TYPE == "subGAIN") & SIZE >= MINSIZE_SUB, 1, FILTERSIZE)) %>%
                  
                  # mutate(FILTER = ifelse((FILTERPROBE+FILTERSIZE) == 2 | TYPE == "NORMAL", "PASS", ".")) %>%
                  mutate(FILTER = ifelse(TYPE == "NORMAL", "na", ".")) %>%
                  mutate(FILTER = ifelse(TYPE != "NORMAL" & (FILTERPROBE+FILTERSIZE) == 2, "PASS", FILTER)) %>%
                  mutate(FILTER = ifelse(TYPE != "NORMAL" & (FILTERPROBE+FILTERSIZE) < 2, "cnv_small_size", FILTER)) %>%
                  
                  mutate(NofNAs = NofNAsum) %>%
                  
                  select(CHRPARTN, SEGMENTATION_CONDITION, PON_SEX_N, GROUP, CASE_ID, CASE_TYPE, CASE_SEX, CONTIG, START, END, MEAN_L2R, MEAN_L2R_EDIT, SIZE, PROBES, FILTER, TYPE, MINLOSS_EDIT, MINLOSS_SUB_EDIT, MINLOSS_BIAL_EDIT, MINGAIN_EDIT, MINGAIN_SUB_EDIT, NofNAs)
                
                PHASE_OUTPUT = rbind(PHASE_OUTPUT, SEGMENTS_3)
                
              }}
            
            PHASE_OUTPUT = filter(PHASE_OUTPUT, CONTIG != "NA")
            TEST = sum(PHASE_OUTPUT$NofNAs)
            # message(TEST)
            
            PHASE_INPUT = PHASE_OUTPUT
            
          }
          
          PHASE_OUTPUT = PHASE_OUTPUT %>% select(-CHRPARTN, -MINLOSS_EDIT, -MINLOSS_SUB_EDIT, -MINLOSS_BIAL_EDIT, -MINGAIN_EDIT, -MINGAIN_SUB_EDIT, -NofNAs)
          
#-------------------------------------------------------------------------------------------
            
            
COPYNUMBER_SEGMENTS = PHASE_OUTPUT
#-------------------------------------------------------------------------------------------            

message(paste0("   R ... ", date(), " - STAGE 3/5 - LOH segmentation"))
            
            
#SEGMENTATION - for MAF
#-------------------------------------------------------------------------------------------
#PHASE1 - create part of chromosomes that are futher than "AFSIZE"

PHASE1 = data.frame(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
names(PHASE1) = c("CONTIG", "POSITION", "REF_COUNT", "ALT_COUNT", "REF_NUCLEOTIDE", "ALT_NUCLEOTIDE", "N", "MAF", "MAFTRANS", "GENOTYPE", "DISTANCE", "N1", "CHRPART", "CHRPARTN", "NofPinCHRPART")
            

LISTofCHR = unique(ALL_DATA_INPUT$CONTIG) %>% as.list()

for (CHR in (LISTofCHR)) {
  
  #find segments of the chromosome that are futher than value GAP and exclude segments that are smaller than MINKEEP value
  ALL_DATA_PH1 = ALL_DATA_INPUT %>%
    filter(CONTIG == CHR) %>%
    arrange(POSITION) %>%
    mutate(DISTANCE = POSITION-rowShift(POSITION, -1)) %>%
    mutate(DISTANCE = ifelse(is.na(DISTANCE), -1, DISTANCE))
  # ALL_DATA_PH1[c("DISTANCE")][is.na(ALL_DATA_PH1[c("DISTANCE")])] <- -1
  ALL_DATA_PH1 = ALL_DATA_PH1 %>%
    mutate(N1 = 1:n()) %>%
    # mutate(CHRPART = ifelse(DISTANCE > (AFSIZE/2) | DISTANCE == -1, "NEWPART", NA))
    mutate(CHRPART = ifelse(DISTANCE > (AFSIZE) | DISTANCE == -1, "NEWPART", NA))
  
  SUBSET = ALL_DATA_PH1 %>%
    filter(CHRPART == "NEWPART") %>%
    mutate(CHRPARTN = 1:n()) %>%
    select(N1, CHRPARTN)
  
  ALL_DATA_PH1 = ALL_DATA_PH1 %>%
    mutate(CHRPARTN = vlookup(N1, dict = SUBSET, lookup_column = which(colnames(SUBSET)=="N1"), result_column = which(colnames(SUBSET)=="CHRPARTN"))) %>%
    mutate(CHRPARTN = fillNAgaps(CHRPARTN))
  
  FINDNofPinCHRPART = ALL_DATA_PH1 %>%
    group_by(CHRPARTN) %>%
    summarise(NofP = n())
  
  ALL_DATA_PH1 = ALL_DATA_PH1 %>%
    mutate(NofPinCHRPART = vlookup(CHRPARTN, dict = FINDNofPinCHRPART, lookup_column = which(colnames(FINDNofPinCHRPART)=="CHRPARTN"), result_column = which(colnames(FINDNofPinCHRPART)=="NofP"))) 
  # filter(NofPinCHRPART >= MINKEEP)
  
  PHASE1 = rbind(PHASE1, ALL_DATA_PH1)
}

PHASE1 = filter(PHASE1, CONTIG != "NA")

#-------------------------------------------------------------------------------------------

#PHASE2 - create segments

PHASE2 = data.frame(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
names(PHASE2) = c("SEGMENTATION_CONDITION", "PON_SEX_N", "GROUP", "CASE_ID", "CASE_TYPE", "CASE_SEX", "CONTIG", "START", "END", "MEAN_L2R", "MEAN_L2R_EDIT", "SIZE", "PROBES", "FILTER", "TYPE")

LISTofCHR = unique(PHASE1$CONTIG) %>% as.list()

for (CHR in (LISTofCHR)) {
  
  ALL_DATA_PH2 = PHASE1 %>% 
    filter(CONTIG == CHR)
  LISTofCHRPART = unique(ALL_DATA_PH2$CHRPARTN) %>% as.list()
  
  
  for (PART in (LISTofCHRPART)) {
    
    ALL_DATA_PART_PH2 = ALL_DATA_PH2 %>%
      filter(CHRPARTN == PART) %>%
      mutate(N3 = 1:n()) %>%
      mutate(SEG2 = ifelse(GENOTYPE == rowShift(GENOTYPE, -1), ".", "NEW")) %>%
      mutate(SEG2 = ifelse(is.na(SEG2), "NEW", SEG2))
    # ALL_DATA_PART_PH2[c("SEG2")][is.na(ALL_DATA_PART_PH2[c("SEG2")])] <- "NEW"
    
    SUBSET = ALL_DATA_PART_PH2 %>%
      filter(SEG2 == "NEW") %>%
      mutate(SEGMENT2 = 1:n()) %>%
      select(N3, SEGMENT2)
    
    ALL_DATA_PART_PH2 = ALL_DATA_PART_PH2 %>%
      mutate(SEGMENT2 = vlookup(N3, dict = SUBSET, lookup_column = which(colnames(SUBSET)=="N3"), result_column = which(colnames(SUBSET)=="SEGMENT2"))) %>%
      mutate(SEGMENT2 = fillNAgaps(SEGMENT2))
    
    SEGMENTS_1 = ALL_DATA_PART_PH2 %>%
      group_by(SEGMENT2) %>%
      summarise(SEG_GENOTYPE = (GENOTYPE %>% unique()), START = min(POSITION), END = max(POSITION), PROBES = n())

    #Smooth
    SEGMENT_x = SEGMENTS_1 %>%
      mutate(SIZE = END-START+1) %>%
      mutate(BEFORE_PROBES = rowShift(PROBES, -1), BEFORE_SIZE = rowShift(SIZE, -1), AFTER_PROBES = rowShift(PROBES, +1), AFTER_SIZE = rowShift(SIZE, +1)) %>%
      mutate(SMOOTH = ifelse(PROBES < MINKEEP & BEFORE_SIZE >= (AFSIZE/2) &  AFTER_SIZE >= (AFSIZE/2), "SKIP", ".")) %>%
      mutate(SMOOTH = ifelse(is.na(SMOOTH), ".", SMOOTH))
    # SEGMENT_x[c("SMOOTH")][is.na(SEGMENT_x[c("SMOOTH")])] <- "."
    SEGMENT_x = SEGMENT_x %>%
      filter(SMOOTH != "SKIP") %>%
      mutate(N4 = 1:n()) %>%
      mutate(SEG3 = ifelse(SEG_GENOTYPE == rowShift(SEG_GENOTYPE, -1), ".", "NEW")) %>%
      mutate(SEG3 = ifelse(is.na(SEG3), "NEW", SEG3))
    # SEGMENT_x[c("SEG3")][is.na(SEGMENT_x[c("SEG3")])] <- "NEW"
    
    SUBSET = SEGMENT_x %>%
      filter(SEG3 == "NEW") %>%
      mutate(SEGMENT3 = 1:n()) %>%
      select(N4, SEGMENT3)
    
    SEGMENT_x = SEGMENT_x %>%
      mutate(SEGMENT3 = vlookup(N4, dict = SUBSET, lookup_column = which(colnames(SUBSET)=="N4"), result_column = which(colnames(SUBSET)=="SEGMENT3"))) %>%
      mutate(SEGMENT3 = fillNAgaps(SEGMENT3))
    
    #Add new segments back to big table - to calculate precize MAFTRANS
    DATA_BACK = PHASE1 %>%
      filter(CONTIG == CHR) %>%
      filter(CHRPARTN == PART) %>%
      mutate(SEGMENT3 = vlookup(POSITION, dict = SEGMENT_x, lookup_column = which(colnames(SEGMENT_x)=="START"), result_column = which(colnames(SEGMENT_x)=="SEGMENT3"))) %>% #3vs14
      mutate(SEGMENT3 = fillNAgaps(SEGMENT3))
    
    SEGMENTS_2 = DATA_BACK %>%
      group_by(SEGMENT3) %>%
      summarise(MAF_N1 = sum(MAFTRANS >= 0.495), MAF_N2 = sum(MAFTRANS < 0.495), START = min(POSITION), END = max(POSITION), PROBES = n()) %>%
      
      mutate(SEGMENTATION_CONDITION = GETSEGMENTCONDITION) %>%
      mutate(CASE_ID = GETCASEID) %>%
      mutate(CASE_TYPE = GETCASETYPE) %>%
      mutate(CASE_SEX = GETCASESEX) %>%
      mutate(GROUP = GETGROUPID) %>%
      mutate(PON_SEX_N = paste0(GETPONSEX, "_", NofCTRL_USED)) %>%
      
      mutate(CONTIG = CHR) %>%
      mutate(CHRPARTN = PART) %>%
      mutate(SIZE = END-START+1) %>%
      mutate(MEAN_L2R = "NA") %>%
      mutate(MEAN_L2R_EDIT = "NA") %>%
      mutate(FILTERPROBE = ifelse(PROBES >= AFPROBES, 1, 0)) %>%
      mutate(FILTERSIZE = ifelse(SIZE >= AFSIZE, 1, 0)) %>%
      mutate(FILTER = ifelse((FILTERPROBE+FILTERSIZE) == 2, "PASS", ".")) %>% #filtered below, all have to be PASS to continue
      
      mutate(MAF_R = MAF_N2/MAF_N1) %>%
      
      mutate(TYPE = ifelse(MAF_R >= 0.015, "HOM", ".")) %>%
      mutate(TYPE = ifelse(MAF_R < 0.015, "OTHER", TYPE)) %>%

      mutate(MEAN_L2R = MAF_R) %>%
      
      select(SEGMENTATION_CONDITION, PON_SEX_N, GROUP, CASE_ID, CASE_TYPE, CASE_SEX, CONTIG, START, END, MEAN_L2R, MEAN_L2R_EDIT, SIZE, PROBES, FILTER, TYPE)
    
    
    PHASE2 = rbind(PHASE2, SEGMENTS_2)
    
  }
  
  }

#-------------------------------------------------------------------------------------------
LOH_SEGMENTS = PHASE2 %>% filter(SIZE >= AFSIZE) %>% filter(PROBES >= AFPROBES)
#-------------------------------------------------------------------------------------------            

    
ALLSEGMENTS = rbind(COPYNUMBER_SEGMENTS, LOH_SEGMENTS)


message(paste0("   R ... ", date(), " - STAGE 4/5 - chromosome bands annotation"))

#Get chromosome bands
BANDSUCSC = read.delim(GETBANDUCSC, header = FALSE, stringsAsFactors = F) %>%
  select(V1, V2, V3, V4) %>%
  mutate(N = 1:n())
names(BANDSUCSC) = c("CONTIG", "START", "END", "BANDID", "N")

ALLSEGMENTS_BAND = ALLSEGMENTS %>%
  select(CONTIG, START, END) %>%
  mutate(CHRBAND_START = ".") %>%
  mutate(CHRBAND_END = ".")
            
LISTofBANDS = unique(BANDSUCSC$N) %>% as.list()
            
for (BAND in (LISTofBANDS)) {
  BANDPART = filter(BANDSUCSC, N == BAND)
  BAND_CONTIG = BANDPART$CONTIG %>% as.character
  BAND_START = BANDPART$START %>% as.integer
  BAND_END = BANDPART$END %>% as.integer
  BAND_ID = BANDPART$BANDID %>% as.character
              
ALLSEGMENTS_BAND = ALLSEGMENTS_BAND %>%
  mutate(CHRBAND_START = ifelse(CONTIG == BAND_CONTIG & START >= BAND_START & START <= BAND_END, BAND_ID, CHRBAND_START)) %>%
  mutate(CHRBAND_END = ifelse(CONTIG == BAND_CONTIG & END >= BAND_START & END <= BAND_END, BAND_ID, CHRBAND_END))
}  
            
ALLSEGMENTS_BAND = ALLSEGMENTS_BAND %>%
  mutate(CHRBAND = ifelse(CHRBAND_START == CHRBAND_END, CHRBAND_START, paste0(CHRBAND_START, "-", CHRBAND_END))) %>%
  mutate(CHRBAND = paste0(CONTIG, CHRBAND)) 
            
ALLSEGMENTS_BAND$CHRBAND = gsub("chr", "", ALLSEGMENTS_BAND$CHRBAND)
ALLSEGMENTS_BAND = select(ALLSEGMENTS_BAND, CHRBAND)
ALLSEGMENTS = cbind(ALLSEGMENTS, ALLSEGMENTS_BAND)
ALLSEGMENTS$START = as.integer(as.numeric(ALLSEGMENTS$START))            
ALLSEGMENTS$END = as.integer(as.numeric(ALLSEGMENTS$END))            

names(ALLSEGMENTS) = c("SEGMENTATION_CONDITION", "PON_SEX_N", "SAMPLE_ID1", "SAMPLE_ID2", "SAMPLE_TYPE", "SAMPLE_SEX", "CONTIG", "START", "END", "MEAN_L2R", "MEAN_L2R_EDIT", "SIZE", "PROBES", "FILTER", "TYPE", "CHRBAND")


write_tsv(ALLSEGMENTS, paste0(OUTPUT_STANDARD), col_names = T)





message(paste0("   R ... ", date(), " - STAGE 5/5 - QC metrics"))


if(file.exists(GETRAWCOUNT)==T) {
#QC CALCULATION
#-------------------------------------------------------------------------------------------

#depth data
RAWCOUNT = data.frame(NA,NA,NA,NA,NA)
names(RAWCOUNT) = c("CONTIG", "START", "END", "COUNT", "DIF")

LISTofCHR = unique(RAWCOUNT_DATA_INPUT$CONTIG) %>% as.list()
for (CHR in (LISTofCHR)) {
  RAWCOUNT_DATA = RAWCOUNT_DATA_INPUT %>%
    filter(CONTIG == CHR) %>%
    arrange(START) %>%
    mutate(DIF = abs(COUNT-rowShift(COUNT, -1)))
  RAWCOUNT = rbind(RAWCOUNT, RAWCOUNT_DATA)
}

RAWCOUNT = filter(RAWCOUNT, CONTIG != "NA")
SD_RAWCOUNT=sd(RAWCOUNT$DIF, na.rm=T)
QUANTILES=quantile(RAWCOUNT$COUNT, probs = c(0, 0.25, 0.5, 0.75, 1)) %>% unname() %>% as.data.frame() %>% t() %>% as.data.frame()

# STANDARD = STANDARD_DATA_INPUT %>% mutate(DIF = abs(LOG2_COPY_RATIO-rowShift(LOG2_COPY_RATIO, -1)))
# SD_STANDARD=sd(STANDARD$DIF, na.rm=T)
# STANDART_AUTOSOME = filter(STANDARD, CONTIG != "chrX" & CONTIG != "chrY" & CONTIG != "X" & CONTIG != "Y")
# SD_STANDARD_AUTOSOME=sd(STANDART_AUTOSOME$DIF, na.rm=T)

DENOIS = DENOIS_DATA_INPUT %>% mutate(DIF = abs(LOG2_COPY_RATIO-rowShift(LOG2_COPY_RATIO, -1)))
SD_DENOIS=sd(DENOIS$DIF, na.rm=T)
DENOIS_AUTOSOME = filter(DENOIS, CONTIG != "chrX" & CONTIG != "chrY" & CONTIG != "X" & CONTIG != "Y")
SD_DENOIS_AUTOSOME=sd(DENOIS_AUTOSOME$DIF, na.rm=T)

SEGMENT_COND = paste0(GETSEGMENTCONDITION) %>% as.data.frame()
CAPTURE_ID = paste0(GETCAPTUREID) %>% as.data.frame()
PROJECT_ID = paste0(GETPROJECTID)  %>% as.data.frame()
GROUP = paste0(GETGROUPID) %>% as.data.frame()
ID = paste0(GETCASEID) %>% as.data.frame()
TYPE = paste0(GETCASETYPE) %>% as.data.frame()
SEX = paste0(GETCASESEX) %>% as.data.frame()
PONSEX = paste0(GETPONSEX) %>% as.data.frame()

#segmentation data
N_CN_SEGM = ALLSEGMENTS %>% filter(TYPE == "NORMAL" | TYPE == "LOSS" | TYPE == "GAIN") %>% nrow() %>% as.integer()
N_CN_SEGM_NORM = ALLSEGMENTS %>% filter(TYPE == "NORMAL") %>% nrow() %>% as.integer()
N_CN_SEGM_LOSS = ALLSEGMENTS %>% filter(TYPE == "LOSS") %>% nrow() %>% as.integer()
N_CN_SEGM_GAIN = ALLSEGMENTS %>% filter(TYPE == "GAIN") %>% nrow() %>% as.integer()

N_CN_SEGM_LOSS_FILTER_PASS = ALLSEGMENTS %>% filter(FILTER == "PASS") %>% filter(TYPE == "LOSS") %>% nrow() %>% as.integer()
N_CN_SEGM_GAIN_FILTER_PASS = ALLSEGMENTS %>% filter(FILTER == "PASS") %>% filter(TYPE == "GAIN") %>% nrow() %>% as.integer()

N_LOH_SEGM = ALLSEGMENTS %>% filter(TYPE == "HOM" | TYPE == "OTHER") %>% nrow() %>% as.integer()
N_LOH_SEGM_FILTER_PASS = ALLSEGMENTS %>% filter(FILTER == "PASS") %>% filter(TYPE == "HOM" | TYPE == "OTHER") %>% nrow() %>% as.integer()


#sum together
RESULTS = cbind(SEGMENT_COND, CAPTURE_ID, PROJECT_ID, GROUP, ID, TYPE, SEX, PONSEX, NofCTRL_USED, NofCTRL_TOTAL,
                QUANTILES, SD_RAWCOUNT, 
                # SD_STANDARD, SD_STANDARD_AUTOSOME, 
                SD_DENOIS, SD_DENOIS_AUTOSOME,
                N_CN_SEGM, N_CN_SEGM_NORM, N_CN_SEGM_LOSS, N_CN_SEGM_GAIN,
                N_CN_SEGM_LOSS_FILTER_PASS, N_CN_SEGM_GAIN_FILTER_PASS,
                N_LOH_SEGM, N_LOH_SEGM_FILTER_PASS)
names(RESULTS) = c("SEGMENTATION_CONDITION", "CAPTURE_ID", "PROJECT_ID", "GROUP_ID", "CASE_ID", "CASE_TYPE", "CASE_SEX", "PON_SEX", "PON_N_CTRL_USED", "PON_N_CTRL_TOTAL",
                   "DEPTH.MIN", "DEPTH.Q25", "DEPTH.MEDIAN", "DEPTH.Q75", "DEPTH.MAX", "QC.RAWCOUNT", 
                   # "QC.STANDARD", "QC.STANDARD.AUTOSOME", 
                   "QC.DENOIS", "QC.DENOIS.AUTOSOME",
                   "N_CN_SEGMENTS", "N_CN_NORMAL", "N_LOSS", "N_GAIN",
                   "N_LOSS_FILTER_PASS", "N_GAIN_FILTER_PASS",
                   "N_LOH_SEGMENTS", "N_LOH_SEGMENTS_FILTER_PASS")

RESULTS = rbind(QCTABLE, RESULTS) %>%
  unique()
write_tsv(RESULTS, QCTABLEIN, col_names = T)
#-------------------------------------------------------------------------------------------
}


message(paste0("   R ... ", date(), " - FINISHED"))

system(paste0("touch status_ok"))


