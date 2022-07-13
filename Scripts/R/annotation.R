args <- commandArgs(trailing = TRUE)

message(paste0("   R ... ", date(), " - STARTING"))

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(splitstackshape))
suppressPackageStartupMessages(library(expss))

GET_DTB_DGVDATA         = args[1]%>% as.character()
GET_GENES_DTB_PROCESSED = args[2] %>% as.character()
GET_MAPPABILITY_DTB     = args[3] %>% as.character()
GET_ENCODE_BLACKLIST    = args[4]%>% as.character()
GET_SEGMENTATION        = args[5] %>% as.character()
GET_CTRL_CN             = args[6] %>% as.character()
GET_CTRL_N              = args[7] %>% as.numeric()
GET_CTRL_TYPE           = args[8] %>% as.character()
MINGAIN                 = args[9] %>% as.numeric() * 0.95
MINGAIN_SUB             = args[10] %>% as.numeric()
MINLOSS                 = args[11] %>% as.numeric() * 0.95
MINLOSS_SUB             = args[12]%>% as.numeric()
OUTPUT                  = args[13]%>% as.character()


message(paste0("   R ... ", date(), " - STAGE 1/4 - loading data"))

#Get data
#segmentation table
ORIG_FILE_INPUT = read.delim(paste0(GET_SEGMENTATION), header=TRUE, stringsAsFactors = F)
#get dataset denoised data
if(file.exists(GET_CTRL_CN)==T) {
  CTRL_CN = readRDS(GET_CTRL_CN)
  #exclude relative controls from dataset
  CTRL_CN_POS = select(CTRL_CN, c(1:3))
  CTRL_CN_DATA = select(CTRL_CN, matches("denoisedCR"))
  MAIN_ID = ORIG_FILE_INPUT$SAMPLE_ID1 %>% unique() %>% as.character()
  X = colnames(CTRL_CN_DATA) %>% as.data.frame()
  names(X) = "NAME"
  X = X %>% mutate(NAME_part = NAME)
  X = cSplit(X,'NAME_part', sep='_', type.convert = F) 
  X = X %>% mutate(MATCH = ifelse(NAME_part_1 == MAIN_ID, "yes", "no")) %>% filter(MATCH == "no")
  X = X$NAME %>% as.vector()
  CTRL_CN_DATA = subset(CTRL_CN, select=X)
  CTRL_CN = cbind(CTRL_CN_POS, CTRL_CN_DATA)
}
#load DGV database
DTB_DGVDATA = readRDS(GET_DTB_DGVDATA) %>% mutate(SIZE = END - START + 1) #%>% mutate(CONTIG = ifelse(grepl("^chr*", CONTIG), CONTIG, paste0("chr", CONTIG)))
#load gene dtb
GENES_DTB_PROCESSED = read.delim(paste0(GET_GENES_DTB_PROCESSED), header=TRUE, stringsAsFactors = F)
#load mappability dtb
MAPPABILITY_DTB = readRDS(GET_MAPPABILITY_DTB)
#load ENCODE blacklist
ENCODE_BLACKLIST = read.delim(paste0(GET_ENCODE_BLACKLIST), header=TRUE, stringsAsFactors = F)

#remove "chr" if not present in sample file
CONTIGS_in_SAMPLE=ORIG_FILE_INPUT$CONTIG %>% unique()
if (("TRUE" %in% grepl("chr", CONTIGS_in_SAMPLE))==FALSE) {
  DTB_DGVDATA$CONTIG = gsub("chr", "", DTB_DGVDATA$CONTIG)
  GENES_DTB_PROCESSED$CONTIG = gsub("chr", "", GENES_DTB_PROCESSED$CONTIG)
  MAPPABILITY_DTB$CONTIG = gsub("chr", "", MAPPABILITY_DTB$CONTIG)
  ENCODE_BLACKLIST$CONTIG = gsub("chr", "", ENCODE_BLACKLIST$CONTIG)
}


#Functions
#----------------------------------------------------------------------------------------------------------------------
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
#----------------------------------------------------------------------------------------------------------------------

message(paste0("   R ... ", date(), " - STAGE 2/4 - annotating abnormal segments"))

#----------------------------------------------------------------------------------------------------------------------
#Prepare original table
ORIG_FILE_INPUT$SAMPLE_SEX = as.character(ORIG_FILE_INPUT$SAMPLE_SEX)
ORIG_FILE_INPUT = ORIG_FILE_INPUT %>%
  mutate(SAMPLE_SEX = ifelse(SAMPLE_SEX == "FALSE", "F", SAMPLE_SEX)) %>%
  mutate(N = 1:n())

ORIG_FILE_INPUT_TO_ANNOTATE = ORIG_FILE_INPUT %>% filter(TYPE == "LOSS" | TYPE == "GAIN" | TYPE == "subLOSS" | TYPE == "subGAIN" | TYPE == "HOM" | TYPE == "OTHER" | TYPE == "HET")
ORIG_FILE_INPUT_NOT_TO_ANNOTATE = ORIG_FILE_INPUT %>% filter(TYPE == "NORMAL")

LISTofN = ORIG_FILE_INPUT_TO_ANNOTATE$N %>% as.list()
ANNOTATION = data.frame()

for (Neach in (LISTofN)) {
  
  TEMP = ORIG_FILE_INPUT_TO_ANNOTATE %>% filter(N == Neach)
  TEMP_CONTIG = TEMP$CONTIG %>% as.character()
  TEMP_START = TEMP$START %>% as.integer()
  TEMP_END = TEMP$END %>% as.integer()
  TEMP_TYPE = TEMP$TYPE %>% as.character()

  # BASE = ORIG_FILE_INPUT %>% filter(N == Neach) %>% select(-c(N))
  
  # if(TEMP_TYPE != "NORMAL") {
  #filter gene database for genes that overlap the segment
  GENES_DTB_PROCESSED_TEMP = GENES_DTB_PROCESSED %>%
    filter(CONTIG == TEMP_CONTIG) %>%
    filter(TxSTART <= TEMP_END) %>%
    filter(TxEND >= TEMP_START) 
    if(nrow(GENES_DTB_PROCESSED_TEMP) > 0) {
      GENE = GENES_DTB_PROCESSED_TEMP$GENE_ID %>% as.vector()
      GENE = paste(as.character(GENE), collapse=";")
    } else {GENE = "non-coding"} #when abnormal segments not overlapping any gene
  # } else {GENE = "-"} #when segment is not abnormal
  
  TEMP_LINE = data.frame(TEMP, GENE, stringsAsFactors = FALSE)
  ANNOTATION = rbind(ANNOTATION, TEMP_LINE)
  
}

ORIG_FILE_INPUT_NOT_TO_ANNOTATE = ORIG_FILE_INPUT_NOT_TO_ANNOTATE %>% mutate(GENE = "-")
ANNOTATION = rbind(ANNOTATION, ORIG_FILE_INPUT_NOT_TO_ANNOTATE) %>% arrange(N) #%>% select(-c("N"))



message(paste0("   R ... ", date(), " - STAGE 3/4 - incorporating annotation by controls, DGV, mappability and ENCODE blacklist"))
#----------------------------------------------------------------------------------------------------------------------

#2) INCORPORATE CTRLs, DGV
SEGMENTATION_TO_ANNOTATE = ANNOTATION %>% filter(TYPE == "LOSS" | TYPE == "GAIN" | TYPE == "subLOSS" | TYPE == "subGAIN")
SEGMENTATION_NOT_TO_ANNOTATE = ANNOTATION %>% filter(TYPE == "NORMAL" | TYPE == "HOM" | TYPE == "OTHER" | TYPE == "HET")

#annotate only if there are any abnormal segments
if(nrow(SEGMENTATION_TO_ANNOTATE)>0) {
#CN
#----------------------------------------------------------------------
ADD = data.frame(0.1, "filterout", 0.1,"filterout", 0.1,0.1,"filterout",0.1,0.1,"filterout", 0.1,"filterout", 0.1, "filterout")
names(ADD) = c("CTRL_N", "CTRL_TYPE", "CTRL_NA_FRACTION", "CTRL_NOISE_PREDICTION", "CTRL_FofGAIN", "CTRL_FofLOSS", "CTRL_CNV_PREDICTION", "DGV_NofGAIN_AROUND", "DGV_NofLOSS_AROUND", "DGV_PREDICTION", "MAPPABILITY", "MAPPABILITY_PREDICTION", "ENCODE_BLACKLIST_OVERLAP", "ENCODE_BLACKLIST_PREDICTION")

GENERAL_CUT_OFF=0.25

for(line in 1:nrow(SEGMENTATION_TO_ANNOTATE)) {
  SEGMENT_CONTIG = SEGMENTATION_TO_ANNOTATE[line,7] %>% as.character()
  SEGMENT_START  = SEGMENTATION_TO_ANNOTATE[line,8] %>% as.integer()
  SEGMENT_END    = SEGMENTATION_TO_ANNOTATE[line,9] %>% as.integer()
  SEGMENT_TYPE   = SEGMENTATION_TO_ANNOTATE[line,15] %>% as.character
  
  #CTRLs
    if(file.exists(GET_CTRL_CN)==T) {
      
      # CTRL_N = GET_CTRL_N
      CTRL_N = ncol(CTRL_CN_DATA)
      CTRL_TYPE = GET_CTRL_TYPE

      CTRL_CN_FILTER = CTRL_CN %>% filter((CONTIG == SEGMENT_CONTIG) & (START >= SEGMENT_START) & (END <= SEGMENT_END))
      CTRL_CN_FILTER = select(CTRL_CN_FILTER, matches("denoisedCR"))
      
      #ctrl noise
      CTRL_N_SQUARE = nrow(CTRL_CN_FILTER)*ncol(CTRL_CN_FILTER)
      CTRL_NA_FRACTION=sum(is.na(CTRL_CN_FILTER))
      CTRL_NA_FRACTION = CTRL_NA_FRACTION/CTRL_N_SQUARE
      # if(CTRL_NA_FRACTION == 0) {CTRL_NOISE_PREDICTION = "none"}
      # if(CTRL_NA_FRACTION > 0) {CTRL_NOISE_PREDICTION = "low"}
      # if(CTRL_NA_FRACTION >= 0.25) {CTRL_NOISE_PREDICTION = "intermediate"}
      # if(CTRL_NA_FRACTION >= 0.5) {CTRL_NOISE_PREDICTION = "high"}

      if(CTRL_NA_FRACTION >= GENERAL_CUT_OFF) {CTRL_NOISE_PREDICTION = "high"} else {CTRL_NOISE_PREDICTION = "low"}
      
      
      
      #frequency of CNV in controls
      CTRL_CN_FILTER_MEAN = colMeans(CTRL_CN_FILTER, na.rm = TRUE) %>% as.vector()
      TOTAL = length(CTRL_CN_FILTER_MEAN) %>% as.numeric()
      
      if(SEGMENT_TYPE == "LOSS" | SEGMENT_TYPE == "GAIN") {
        CTRL_FofLOSS_SUM = sum(CTRL_CN_FILTER_MEAN <= MINLOSS, na.rm = TRUE) %>% as.numeric()
        CTRL_FofLOSS     = CTRL_FofLOSS_SUM/TOTAL %>% as.numeric()
        CTRL_FofGAIN_SUM = sum(CTRL_CN_FILTER_MEAN >= MINGAIN, na.rm = TRUE) %>% as.numeric()
        CTRL_FofGAIN     = CTRL_FofGAIN_SUM/TOTAL %>% as.numeric()
      }
      if(SEGMENT_TYPE == "subLOSS" | SEGMENT_TYPE == "subGAIN") {
        CTRL_FofLOSS_SUM = sum(CTRL_CN_FILTER_MEAN <= MINLOSS_SUB, na.rm = TRUE) %>% as.numeric()
        CTRL_FofLOSS     = CTRL_FofLOSS_SUM/TOTAL %>% as.numeric()
        CTRL_FofGAIN_SUM = sum(CTRL_CN_FILTER_MEAN >= MINGAIN_SUB, na.rm = TRUE) %>% as.numeric()
        CTRL_FofGAIN     = CTRL_FofGAIN_SUM/TOTAL %>% as.numeric()
      }
      
      CTRL_FofCNV = (CTRL_FofLOSS_SUM+CTRL_FofGAIN_SUM)/TOTAL
      
      # if(CTRL_FofCNV == 0) {CTRL_CNV_PREDICTION = "none"}
      # if(CTRL_FofCNV > 0) {CTRL_CNV_PREDICTION = "low"}
      # if(CTRL_FofCNV >= 0.05) {CTRL_CNV_PREDICTION = "intermediate"}
      # if(CTRL_FofCNV >= 0.1) {CTRL_CNV_PREDICTION = "high"}
      
      if(CTRL_FofCNV >= GENERAL_CUT_OFF) {CTRL_CNV_PREDICTION = "high"} else {CTRL_CNV_PREDICTION = "low"}
      

    }
    
    if(file.exists(GET_CTRL_CN)==F) {
      CTRL_N = "not_performed"
      CTRL_TYPE = "not_performed"
      CTRL_NA_FRACTION = "not_performed"
      CTRL_NOISE_PREDICTION = "not_performed"
      CTRL_FofGAIN = "not_performed"
      CTRL_FofLOSS = "not_performed"
      CTRL_CNV_PREDICTION = "not_performed"
    }
  
    
  
  # DGV
    #---------------------------------------------------------
    SEGM_SIZE = SEGMENT_END - SEGMENT_START + 1
    
    DTB_DGV_DATA_PROCESSED_AROUND = DTB_DGVDATA %>% 
      filter((CONTIG == SEGMENT_CONTIG) & (START >= (SEGMENT_START - SEGM_SIZE)) & (END <= (SEGMENT_END + SEGM_SIZE))) %>%
      filter((START <= SEGMENT_END) & (END >= SEGMENT_START)) %>%
      filter(SIZE <= (SEGM_SIZE + (SEGM_SIZE*0.5)) & SIZE >= (SEGM_SIZE - (SEGM_SIZE*0.5)))
    
    DGV_NofLOSS_AROUND = filter(DTB_DGV_DATA_PROCESSED_AROUND, TYPE == "LOSS" | TYPE == "COMPLEX")
    DGV_NofLOSS_AROUND = DGV_NofLOSS_AROUND$TYPE %>% as.vector()
    DGV_NofLOSS_AROUND = length(DGV_NofLOSS_AROUND) #/ (SEGM_SIZE/1000)
    DGV_NofGAIN_AROUND = filter(DTB_DGV_DATA_PROCESSED_AROUND, TYPE == "GAIN" | TYPE == "COMPLEX")
    DGV_NofGAIN_AROUND = DGV_NofGAIN_AROUND$TYPE %>% as.vector()
    DGV_NofGAIN_AROUND = length(DGV_NofGAIN_AROUND) #/ (SEGM_SIZE/1000)
    
    if(SEGMENT_TYPE == "GAIN" | SEGMENT_TYPE == "subGAIN") {
      if(DGV_NofGAIN_AROUND >= (GENERAL_CUT_OFF*100)) {DGV_PREDICTION = "high"}
      if(DGV_NofGAIN_AROUND < (GENERAL_CUT_OFF*100)) {DGV_PREDICTION = "low"}
    }
    if(SEGMENT_TYPE == "LOSS" | SEGMENT_TYPE == "subLOSS") {
      if(DGV_NofLOSS_AROUND >= (GENERAL_CUT_OFF*100)) {DGV_PREDICTION = "high"}
      if(DGV_NofLOSS_AROUND < (GENERAL_CUT_OFF*100)) {DGV_PREDICTION = "low"}
    }
  
    
    #mappability
    MAPPABILITY_TEMP = MAPPABILITY_DTB %>%
      filter(CONTIG == SEGMENT_CONTIG) %>% 
      filter(START <= SEGMENT_END) %>% 
      filter(END >= SEGMENT_START) %>%
      mutate(START = ifelse(START < SEGMENT_START, SEGMENT_START, START)) %>%
      mutate(END = ifelse(END > SEGMENT_END, SEGMENT_END, END)) %>%
      mutate(SIZE = END - START + 1) %>%
      mutate(MAPPABILITY_SUM = SIZE*VALUE)
    
    MAPPABILITY_BEST = SEGMENT_END - SEGMENT_START + 1
    MAPPABILITY_REAL = sum(MAPPABILITY_TEMP$MAPPABILITY_SUM)
    MAPPABILITY =  MAPPABILITY_REAL/MAPPABILITY_BEST
    
    # if(MAPPABILITY >= 0.5) {MAPPABILITY_PREDICTION = "high"}
    # if(MAPPABILITY < 0.5) {MAPPABILITY_PREDICTION = "low"}
    
    if(MAPPABILITY >= (1-GENERAL_CUT_OFF)) {MAPPABILITY_PREDICTION = "high"} else {MAPPABILITY_PREDICTION = "low"}


    #ENCODE black list
      SEGM_SIZE = SEGMENT_END - SEGMENT_START + 1
      
      ENCODE_BLACKLIST_TEMP = ENCODE_BLACKLIST %>%
        filter(CONTIG == SEGMENT_CONTIG) %>% 
        filter(START <= SEGMENT_END) %>% 
        filter(END >= SEGMENT_START) 
      
      if(nrow(ENCODE_BLACKLIST_TEMP) == 0) {
        ENCODE_BLACKLIST_OVERLAP=0
        ENCODE_BLACKLIST_PREDICTION="no"
      }
      
      if(nrow(ENCODE_BLACKLIST_TEMP) > 0) {
        ENCODE_BLACKLIST_TEMP = ENCODE_BLACKLIST_TEMP%>%
          mutate(START = ifelse(START < SEGMENT_START, SEGMENT_START, START)) %>%
          mutate(END = ifelse(END > SEGMENT_END, SEGMENT_END, END)) %>%
          mutate(SIZE = END - START + 1)
        
        ENCODE_BLACKLIST_OVERLAP=sum(ENCODE_BLACKLIST_TEMP$SIZE)/SEGM_SIZE
        if(ENCODE_BLACKLIST_OVERLAP >= GENERAL_CUT_OFF) {ENCODE_BLACKLIST_PREDICTION="yes"}
        if(ENCODE_BLACKLIST_OVERLAP < GENERAL_CUT_OFF) {ENCODE_BLACKLIST_PREDICTION="no"}
        
      }

  #---------------------------------------------------------
  
  x = data.frame(CTRL_N, CTRL_TYPE, CTRL_NA_FRACTION, CTRL_NOISE_PREDICTION, CTRL_FofGAIN, CTRL_FofLOSS, CTRL_CNV_PREDICTION, DGV_NofGAIN_AROUND, DGV_NofLOSS_AROUND, DGV_PREDICTION, MAPPABILITY, MAPPABILITY_PREDICTION, ENCODE_BLACKLIST_OVERLAP, ENCODE_BLACKLIST_PREDICTION)
  
  ADD = rbind(ADD, x)
  ADD = ADD %>% filter(CTRL_CNV_PREDICTION != "filterout")
}

#----------------------------------------------------------------------
SEGMENTATION_TO_ANNOTATE = cbind(SEGMENTATION_TO_ANNOTATE, ADD)
}

SEGMENTATION_NOT_TO_ANNOTATE = SEGMENTATION_NOT_TO_ANNOTATE %>%
  mutate(CTRL_N = "NA") %>%
  mutate(CTRL_TYPE = "NA") %>%
  mutate(CTRL_NA_FRACTION = "NA") %>%
  mutate(CTRL_NOISE_PREDICTION = "NA") %>%
  mutate(CTRL_FofGAIN = "NA") %>%
  mutate(CTRL_FofLOSS = "NA") %>%
  mutate(CTRL_CNV_PREDICTION = "NA") %>%
  mutate(DGV_NofGAIN_AROUND = "NA") %>%
  mutate(DGV_NofLOSS_AROUND = "NA") %>%
  mutate(DGV_PREDICTION = "NA") %>%
  mutate(MAPPABILITY = "NA") %>%
  mutate(MAPPABILITY_PREDICTION = "NA") %>%
  mutate(ENCODE_BLACKLIST_OVERLAP = "NA") %>%
  mutate(ENCODE_BLACKLIST_PREDICTION = "NA")

if(nrow(SEGMENTATION_TO_ANNOTATE)>0) {SEGMENTATION = rbind(SEGMENTATION_TO_ANNOTATE, SEGMENTATION_NOT_TO_ANNOTATE) %>% arrange(N)}
if(nrow(SEGMENTATION_TO_ANNOTATE)==0) {SEGMENTATION = SEGMENTATION_NOT_TO_ANNOTATE}

#----------------------------------------------------------------------

#change FILTER
SEGMENTATION = SEGMENTATION %>%
  mutate(FILTER = ifelse(FILTER != "PASS" & CTRL_NOISE_PREDICTION == "high", paste0(FILTER, ";", "ctrl_high_noise"), FILTER)) %>%
  mutate(FILTER = ifelse(FILTER == "PASS" & CTRL_NOISE_PREDICTION == "high", "ctrl_high_noise", FILTER)) %>%
  mutate(FILTER = ifelse(FILTER != "PASS" & CTRL_CNV_PREDICTION == "high", paste0(FILTER, ";", "ctrl_high_cnv"), FILTER)) %>%
  mutate(FILTER = ifelse(FILTER == "PASS" & CTRL_CNV_PREDICTION == "high", "ctrl_high_cnv", FILTER)) %>%
  mutate(FILTER = ifelse(FILTER != "PASS" & DGV_PREDICTION == "high", paste0(FILTER, ";", "dgv_high_cnv"), FILTER)) %>%
  mutate(FILTER = ifelse(FILTER == "PASS" & DGV_PREDICTION == "high", "dgv_high_cnv", FILTER)) %>%
  mutate(FILTER = ifelse(FILTER != "PASS" & MAPPABILITY_PREDICTION == "low", paste0(FILTER, ";", "mappability_low"), FILTER)) %>%
  mutate(FILTER = ifelse(FILTER == "PASS" & MAPPABILITY_PREDICTION == "low", "mappability_low", FILTER)) %>%
  mutate(FILTER = ifelse(FILTER != "PASS" & ENCODE_BLACKLIST_PREDICTION == "yes", paste0(FILTER, ";", "encode_blacklist"), FILTER)) %>%
  mutate(FILTER = ifelse(FILTER == "PASS" & ENCODE_BLACKLIST_PREDICTION == "yes", "encode_blacklist", FILTER))

message(paste0("   R ... ", date(), " - STAGE 4/4 - writing output"))


#improve how the columns are called
#----------------------------------------------------------------------------------------------------------------------
names(SEGMENTATION) = c(colnames(ANNOTATION), "CONTROLS_N", "CONTROLS_DENOIS_STYLE", "CONTROLS_NA_FRACTION", "CONTROLS_NOISE_PREDICTION", 
                        "CONTROLS_GAIN_FREQ", "CONTROLS_LOSS_FREQ", "CONTROLS_CN_FREQ_PREDICTION", "DGV_GAIN_n_per_kb_around", "DGV_LOSS_n_per_kb_around", "DGV_CN_FREQ_PREDICTION",
                        "MAPPABILITY", "MAPPABILITY_PREDICTION", "ENCODE_BLACKLIST_OVERLAP", "ENCODE_BLACKLIST_PREDICTION")
SEGMENTATION = SEGMENTATION %>% select(-c("N"))
#----------------------------------------------------------------------------------------------------------------------


#Write output
#----------------------------------------------------------------------------------------------------------------------
write_tsv(SEGMENTATION, paste0(OUTPUT), col_names = T)

message(paste0("   R ... ", date(), " - FINISHED"))

system(paste0("touch status_ok"))


