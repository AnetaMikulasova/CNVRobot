args <- commandArgs(trailing = TRUE)

message(paste0("   R ... ", date(), " - STARTING"))

options("scipen"=100)

# library(ggplot2)
# library(dplyr)
# library(grid)
suppressPackageStartupMessages(library(expss))  #vlookup
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(karyoploteR))
suppressPackageStartupMessages(library(splitstackshape))
# https://bioconductor.org/packages/devel/bioc/vignettes/karyoploteR/inst/doc/karyoploteR.html
# https://bernatgel.github.io/karyoploter_tutorial/
# https://www.statmethods.net/advgraphs/parameters.html

GETCAPTURE          = args[1] %>% as.character()
ASSEMBLY            = args[2] %>% as.character()
CHROMSIZE           = args[3] %>% as.character()
CENTROMERE          = args[4] %>% as.character()
DGV                 = args[5] %>% as.character()
# GNOMADSV            = args[6] %>% as.character()
GET_MAPPABILITY_DTB = args[6] %>% as.character()
GENES_DTB_PROCESSED = args[7] %>% as.character()
PROJECT_ID          = args[8] %>% as.character()
INDIR               = args[9] %>% as.character()
GNOMAD_PATTERN      = args[10] %>% as.character()
GROUP               = args[11] %>% as.character()
SAMPLE1_ID          = args[12] %>% as.character()
SAMPLE1_TYPE        = args[13] %>% as.character()
SAMPLE1_SEX         = args[14] %>% as.character()
SAMPLE2_ID          = args[15] %>% as.character()
SAMPLE2_TYPE        = args[16] %>% as.character()
SAMPLE2_SEX         = args[17] %>% as.character()
SAMPLE3_ID          = args[18] %>% as.character()
SAMPLE3_TYPE        = args[19] %>% as.character()
SAMPLE3_SEX         = args[20] %>% as.character()
PON_STYLE           = args[21] %>% as.character()
OUTDIR              = args[22] %>% as.character()
OUTDIR_IGV          = args[23] %>% as.character()
GETSEGMENTCONDITION = args[24] %>% as.character()
MINLOSS  = args[25] %>% as.numeric()
MINGAIN  = args[26] %>% as.numeric()
GET_CTRL_CN         = args[27] %>% as.character()
PLOT_PROJECT_SOURCE = args[28] %>% as.character()
  
message(paste0("   R ... ", date(), " - STAGE 1/5 - loading and processing databases"))

#Get ctrl CN file
if(file.exists(GET_CTRL_CN)==T) {
  #get data
  CTRL_CN = readRDS(GET_CTRL_CN)
  #exclude relative controls from dataset
  CTRL_CN_POS = select(CTRL_CN, c(1:3))
  CTRL_CN_DATA = select(CTRL_CN, matches("denoisedCR"))
  X = colnames(CTRL_CN_DATA) %>% as.data.frame()
  names(X) = "NAME"
  X = X %>% mutate(NAME_part = NAME)
  X = cSplit(X,'NAME_part', sep='_', type.convert = F) 
  X = X %>% mutate(MATCH = ifelse(NAME_part_1 == GROUP, "yes", "no")) %>% filter(MATCH == "no")
  X = X$NAME %>% as.vector()
  CTRL_CN_DATA = subset(CTRL_CN, select=X)
  CTRL_CN = cbind(CTRL_CN_POS, CTRL_CN_DATA)
  #treat controls using segmentation
  CTRL_CN_POS = CTRL_CN_POS
  CTRL_CN_DATA = CTRL_CN_DATA
  N = ncol(CTRL_CN_DATA)
  CTRL_CN_PREP = CTRL_CN_POS %>% mutate("NofALL" = N)
  CTRL_CN_PREP$NofCASES = rowSums(CTRL_CN_DATA == 0 | CTRL_CN_DATA > 0 | CTRL_CN_DATA < 0, na.rm=T)
  CTRL_CN_PREP$NofLOSS = rowSums(CTRL_CN_DATA < MINLOSS, na.rm=T)
  CTRL_CN_PREP$NofGAIN = rowSums(CTRL_CN_DATA > MINGAIN, na.rm=T)
  CTRL_CN = CTRL_CN_PREP %>%
    mutate("FofCASES" = NofCASES/NofALL) %>%
    mutate("FofLOSS" = NofLOSS/NofCASES) %>%
    mutate("FofGAIN" = NofGAIN/NofCASES)
  #from original
  CTRL_CN = CTRL_CN %>% mutate(STARText = START-100, ENDext = END+100)
  CTRL_CN_temp = CTRL_CN %>% select(CONTIG, START, END)
  names(CTRL_CN_temp) = c("seqnames", "start", "end")
  CTRL_CN = cbind(CTRL_CN_temp, CTRL_CN) %>% toGRanges()
  remove(CTRL_CN_temp)
}

#Get chromosome sizes to generate size of final figure of the chromosome
DTB_CHRSIZE = read.delim(CHROMSIZE, header = FALSE, stringsAsFactors = F)
names(DTB_CHRSIZE) = c("CONTIG", "SIZE")
DTB_CHRSIZE = DTB_CHRSIZE %>% mutate(CONTIG = ifelse(grepl("^chr*", CONTIG), CONTIG, paste0("chr", CONTIG)))
DTB_CHRSIZE$SIZE= as.numeric(as.integer(DTB_CHRSIZE$SIZE))


#Get capture
if(file.exists(GETCAPTURE)) {
  DTB_CAPTURE = read.delim(GETCAPTURE, header = FALSE, stringsAsFactors = F) %>% select(V1, V2, V3)
  names(DTB_CAPTURE) = c("CONTIG", "START", "END")
  DTB_CAPTURE = DTB_CAPTURE %>%
                mutate(CONTIG = ifelse(grepl("^chr*", CONTIG), CONTIG, paste0("chr", CONTIG))) %>%
                filter(CONTIG != "chrMT") %>%
                mutate(STARText1 = START-10000, ENDext1 = END+10000) %>%
                mutate(STARText2 = START-7000, ENDext2 = END+7000) %>%
                mutate(STARText3 = START-76, ENDext3 = END+76)
  DTB_CAPTURE_temp = DTB_CAPTURE %>% select(CONTIG, START, END)
  names(DTB_CAPTURE_temp) = c("seqnames", "start", "end")
  DTB_CAPTURE = cbind(DTB_CAPTURE_temp, DTB_CAPTURE) %>% toGRanges()
  remove(DTB_CAPTURE_temp)
}


#Get chromosome centromeres
DTB_CENTROMERE = read.delim(CENTROMERE, header = FALSE, stringsAsFactors = F)
names(DTB_CENTROMERE) = c("CONTIG", "POSITION")
DTB_CENTROMERE = DTB_CENTROMERE %>% mutate(CONTIG = ifelse(grepl("^chr*", CONTIG), CONTIG, paste0("chr", CONTIG)))
DTB_CENTROMERE_temp = DTB_CENTROMERE %>% mutate(POSITION2 = POSITION) 
names(DTB_CENTROMERE_temp) = c("seqnames", "start", "end")
DTB_CENTROMERE = cbind(DTB_CENTROMERE_temp, DTB_CENTROMERE) %>% toGRanges()
remove(DTB_CENTROMERE_temp)

#Get DGV
DTB_DGVDATA = readRDS(DGV) %>% mutate(CONTIG = ifelse(grepl("^chr*", CONTIG), CONTIG, paste0("chr", CONTIG)))
DTB_DGVDATA_temp = DTB_DGVDATA %>% select(CONTIG, START, END) 
names(DTB_DGVDATA_temp) = c("seqnames", "start", "end")
DTB_DGVDATA = cbind(DTB_DGVDATA_temp, DTB_DGVDATA)
DTB_DGV_LOSS = filter(DTB_DGVDATA, TYPE == "LOSS") %>% toGRanges()
DTB_DGV_GAIN = filter(DTB_DGVDATA, TYPE == "GAIN") %>% toGRanges()
DTB_DGV_COMPLEX = filter(DTB_DGVDATA, TYPE == "COMPLEX") %>% toGRanges()
remove(DTB_DGVDATA, DTB_DGVDATA_temp)

# #Get gnomAD structure database
# DTB_GNOMADSV = readRDS(GNOMADSV) %>% mutate(CONTIG = ifelse(grepl("^chr*", CONTIG), CONTIG, paste0("chr", CONTIG)))
# DTB_GNOMADSV_temp = DTB_GNOMADSV %>% select(CONTIG, START, END)
# names(DTB_GNOMADSV_temp) = c("seqnames", "start", "end")
# DTB_GNOMADSV = cbind(DTB_GNOMADSV_temp, DTB_GNOMADSV)
# DTB_GNOMADSV_LOSS = filter(DTB_GNOMADSV, TYPE == "DEL") %>% toGRanges()
# DTB_GNOMADSV_GAIN = filter(DTB_GNOMADSV, TYPE == "DUP") %>% toGRanges()
# DTB_GNOMADSV_COMPLEX = filter(DTB_GNOMADSV, TYPE == "CPX") %>% toGRanges()
# remove(DTB_GNOMADSV, DTB_GNOMADSV_temp)

#load mappability dtb
MAPPABILITY = readRDS(GET_MAPPABILITY_DTB)

#Get processed gene database
DTB_GENES_DTB_PROCESSED  = read.delim(GENES_DTB_PROCESSED, header = TRUE, stringsAsFactors = F)


message(paste0("   R ... ", date(), " - STAGE 2/5 - loading samples data"))

#Get data
for (i in list(paste0(SAMPLE1_TYPE, "_", SAMPLE1_ID, "_", SAMPLE1_SEX), paste0(SAMPLE2_TYPE, "_", SAMPLE2_ID, "_", SAMPLE2_SEX), paste0(SAMPLE3_TYPE, "_", SAMPLE3_ID, "_", SAMPLE3_SEX))) {

  if(PON_STYLE == "M") {PON = "PoN-M"}
  if(PON_STYLE == "F") {PON = "PoN-F"}
  if(PON_STYLE == "mixed") {PON = "PoN-mixed"}
  if(PON_STYLE == "none") {PON = "PoN-none"}
  if(PON_STYLE == "matched_main") {PON = paste0("PoN-", SAMPLE1_SEX)}
  if(PON_STYLE == "matched_each" & i == paste0(SAMPLE1_TYPE, "_", SAMPLE1_ID, "_", SAMPLE1_SEX)) {PON = paste0("PoN-", SAMPLE1_SEX)}
  if(PON_STYLE == "matched_each" & i == paste0(SAMPLE2_TYPE, "_", SAMPLE2_ID, "_", SAMPLE2_SEX)) {PON = paste0("PoN-", SAMPLE2_SEX)}
  if(PON_STYLE == "matched_each" & i == paste0(SAMPLE3_TYPE, "_", SAMPLE3_ID, "_", SAMPLE3_SEX)) {PON = paste0("PoN-", SAMPLE3_SEX)}
  
  # #define sex for case and for PoN
  # if(i == paste0(SAMPLE1_TYPE, "_", SAMPLE1_ID, "_", SAMPLE1_SEX)) {SEX_A = paste0("PoN-", SAMPLE1_SEX)}
  # if(i == paste0(SAMPLE2_TYPE, "_", SAMPLE2_ID, "_", SAMPLE2_SEX)) {SEX_A = paste0("PoN-", SAMPLE2_SEX)}
  # if(i == paste0(SAMPLE3_TYPE, "_", SAMPLE3_ID, "_", SAMPLE3_SEX)) {SEX_A = paste0("PoN-", SAMPLE3_SEX)}
  # SEX_B = PON
  
  if(file.exists(paste0(INDIR, i, "_SEGMENTS_", PON, "_", GETSEGMENTCONDITION, "_segmentation.tsv"))) {

    MODEL <- read.delim(paste0(INDIR, i, "_SEGMENTS_", PON, "_", GETSEGMENTCONDITION, "_segmentation.tsv"), stringsAsFactors = F) %>% mutate(CONTIG = ifelse(grepl("^chr*", CONTIG), CONTIG, paste0("chr", CONTIG))) %>%
      # filter(FILTER == "PASS") %>%
      filter(TYPE == "NORMAL" | TYPE == "LOSS" | TYPE == "GAIN" | TYPE == "subLOSS" | TYPE == "subGAIN") %>%
      mutate(COLOR = "#000000") %>%
      mutate(COLOR = ifelse(TYPE == "LOSS", "#C02600", COLOR)) %>%
      mutate(COLOR = ifelse(TYPE == "GAIN", "#006200", COLOR))
    
    # #In many cases, MODEL_LOH is empty, therefore one helping line 'x' is created here
    # # x = data.frame(GETSEGMENTCONDITION, "x", "x", "x", "x", "x", "chrZ", 1, 1, NA, NA, 1, 1, "x", "x", "x", "x",   NA, NA, NA, NA, NA, NA, NA, NA, "x")
    # x = MODEL[1,] %>% mutate(CONTIG = "chrZ")
    # MODEL_LOH <- read.delim(paste0(INDIR, i, "_SEGMENTS_", PON, "_", GETSEGMENTCONDITION, ".tsv"), stringsAsFactors = F)%>% mutate(CONTIG = ifelse(grepl("^chr*", CONTIG), CONTIG, paste0("chr", CONTIG))) %>%
    #   # filter(FILTER == "PASS") %>%
    #   filter(TYPE == "HOM" | TYPE == "HET" | TYPE == "OTHER") %>%
    #   mutate(COLOR = ifelse(TYPE == "HOM", "#B0C0E2", "#DBE3F1"))
    # MODEL_LOH_COLUMNS = colnames(MODEL_LOH) %>% as.vector()
    # names(x) = MODEL_LOH_COLUMNS
    # MODEL_LOH = rbind(MODEL_LOH, x)
    # names(MODEL_LOH) = MODEL_LOH_COLUMNS

    MODEL_LOH <- read.delim(paste0(INDIR, i, "_SEGMENTS_", PON, "_", GETSEGMENTCONDITION, "_segmentation.tsv"), stringsAsFactors = F)%>% mutate(CONTIG = ifelse(grepl("^chr*", CONTIG), CONTIG, paste0("chr", CONTIG))) %>%
      # filter(FILTER == "PASS") %>%
      filter(TYPE == "HOM" | TYPE == "HET" | TYPE == "OTHER") %>%
      mutate(COLOR = ifelse(TYPE == "HOM", "#B0C0E2", "#DBE3F1"))
    
    #bed file for IGV with abnormal segments
    IGV_MODEL <- read.delim(paste0(INDIR, i, "_SEGMENTS_", PON, "_", GETSEGMENTCONDITION, "_segmentation.tsv"), stringsAsFactors = F) %>% #mutate(CONTIG = ifelse(grepl("^chr*", CONTIG), CONTIG, paste0("chr", CONTIG))) %>%
      # filter(FILTER == "PASS") %>%
      select(CONTIG, START, END, TYPE) %>%
      mutate(START = START-1) %>%
      filter(TYPE != "NORMAL") %>%
      mutate(x1 = ".") %>%
      mutate(x2 = ".") %>%
      mutate(START2 = START) %>%
      mutate(END2 = END) %>%
      mutate(COLOR = ifelse(TYPE == "LOSS", "192,38,0", ".")) %>%
      mutate(COLOR = ifelse(TYPE == "GAIN", "0,98,0", COLOR)) %>% 
      mutate(COLOR = ifelse(TYPE == "HOM", "142,169,219", COLOR)) %>%
      mutate(COLOR = ifelse(TYPE == "OTHER", "142,169,219", COLOR)) %>%
      mutate(COLOR = ifelse(TYPE == "subLOSS", "255,172,168", COLOR)) %>%
      mutate(COLOR = ifelse(TYPE == "subGAIN", "168,211,121", COLOR)) %>%
      arrange(CONTIG, START, END)
      
    write_tsv(IGV_MODEL, path = paste0(OUTDIR_IGV, GROUP, "_", i, "_abnormal_segments.bed"), col_names = F)
    
    DENOIS  <-  read.delim(paste0(INDIR, i, "_denoisedCR_", PON, ".tsv"), comment.char = "@", stringsAsFactors = F) %>%
      mutate(CONTIG = ifelse(grepl("^chr*", CONTIG), CONTIG, paste0("chr", CONTIG))) %>%
      mutate(COLOR = "#A6A6A6") %>%
      # mutate(COLOR = ifelse(LOG2_COPY_RATIO < -0.9, "#FF7100", COLOR)) %>%
      # mutate(COLOR = ifelse(LOG2_COPY_RATIO < -1.5, "#FF2600", COLOR)) %>%
      # mutate(COLOR = ifelse(LOG2_COPY_RATIO > 0.65, "#008F00", COLOR)) %>%
      # mutate(COLOR = ifelse(LOG2_COPY_RATIO > 1.00, "#47AFEB", COLOR)) %>%
      mutate(COLOLCONTIG = ".") %>%
      # mutate(COLOLCONTIG = ifelse(CONTIG == "chr1" | CONTIG == "chr3" | CONTIG == "chr5" | CONTIG == "chr7" | CONTIG == "chr9" | CONTIG == "chr11" | CONTIG == "chr13" | CONTIG == "chr15" | CONTIG == "chr17" | CONTIG == "chr19" | CONTIG == "chr21" | CONTIG == "chrX", "#941100", "#011893")) %>%
      mutate(COLOLCONTIG = ifelse(CONTIG == "chr1" | CONTIG == "chr3" | CONTIG == "chr5" | CONTIG == "chr7" | CONTIG == "chr9" | CONTIG == "chr11" | CONTIG == "chr13" | CONTIG == "chr15" | CONTIG == "chr17" | CONTIG == "chr19" | CONTIG == "chr21" | CONTIG == "chrX", "#FEB49B", "#9BC2E6")) %>%
      mutate(L2R = ifelse(LOG2_COPY_RATIO < -2.25, -2.25, LOG2_COPY_RATIO)) %>%
      mutate(L2R = ifelse(LOG2_COPY_RATIO > 2.25, 2.25, L2R)) %>%
      mutate(COLORABN = COLOR)
    
          #re-color denois based on segments from segmentations
          MODELtemp = MODEL %>% select(CONTIG, START, END, TYPE)
          write_tsv(MODELtemp, path = paste0(INDIR, "temp1"), col_names = F)
          MODELtemp = paste0(INDIR, "temp1")
          
          DENOIStemp = DENOIS %>% select(CONTIG, START, END)
          write_tsv(DENOIStemp, path = paste0(INDIR, "temp2"), col_names = F)
          DENOIStemp = paste0(INDIR, "temp2")
          
          #process segments to denois data by bedtools and remove temp files
          x = read.table(text = system(paste0("bedtools intersect -a ", DENOIStemp, " -b ", MODELtemp, " -wb -wa"), intern=T), stringsAsFactors = F) %>%
            mutate(IDENTIF = paste0(V1, ".", V2, ".", V3)) %>%
            select(IDENTIF, V7)
          system(paste0("rm ", DENOIStemp, " ", MODELtemp))
          
          DENOIStemp = DENOIS %>% 
            mutate(IDENTIF = paste0(CONTIG, ".", START, ".", END)) %>%
            select(IDENTIF) %>%
            mutate(TYPE = vlookup(IDENTIF, x, lookup_column = 1, result_column = 2)) %>%
            mutate(TYPE = ifelse(is.na(TYPE), ".", TYPE)) %>%
            select(TYPE)
          
          DENOIS = cbind(DENOIS, DENOIStemp) %>%
            mutate(COLORABN = ifelse(TYPE == "subLOSS",  "#FFACA9", COLORABN)) %>%
            mutate(COLORABN = ifelse(TYPE == "subGAIN",  "#A8D379", COLORABN)) %>%
            mutate(COLORABN = ifelse(TYPE == "LOSS",  "#FF2600", COLORABN)) %>%
            mutate(COLORABN = ifelse(TYPE == "GAIN",  "#008F00", COLORABN))
    
    DENOISABN = DENOIS %>%
      filter(TYPE == "LOSS" | TYPE == "GAIN")
    
    DENOIS_FILTERED = DENOIS %>% sample_frac(0.5) %>% arrange(CONTIG, START)
    
    #bedGraph file with denoised copy-number data
    IGV_DENOIS  <-  read.delim(paste0(INDIR, i, "_denoisedCR_", PON, ".tsv"), comment.char = "@", stringsAsFactors = F) %>%
      mutate(L2R = ifelse(LOG2_COPY_RATIO < -2.25, -2.25, LOG2_COPY_RATIO)) %>%
      mutate(L2R = ifelse(LOG2_COPY_RATIO > 2.25, 2.25, L2R)) %>%
      select(CONTIG, START, END, L2R) %>%
      mutate(START = START-1) %>%
      arrange(CONTIG, START, END)
    write_tsv(IGV_DENOIS, path = paste0(OUTDIR_IGV, GROUP, "_", i, "_CN"), col_names = F)
    
    # ALL    <-   readRDS(paste0(INDIR, i, "_allelicCounts", GNOMAD_PATTERN, "_", PON, ".rds")) %>%
    ALL    <-   readRDS(paste0(INDIR, i, "_allelicCounts", GNOMAD_PATTERN, ".rds")) %>%
      mutate(N = REF_COUNT + ALT_COUNT) %>%
      mutate(CONTIG = ifelse(grepl("^chr*", CONTIG), CONTIG, paste0("chr", CONTIG)))
    ALL = ALL %>%
      mutate(TEST = sample(1:2, nrow(ALL), replace = T)) %>%
      mutate(MAF = ifelse(TEST == 1, ALT_COUNT/N, REF_COUNT/N)) %>%
      mutate(COL = "#000000") %>%
      mutate(COL = ifelse(MAF > 0.375 & MAF < 0.625, "#008F00", COL)) %>%
      mutate(COL = ifelse(MAF > 0.125 & MAF < 0.375, "#FF2600", COL)) %>%
      mutate(COL = ifelse(MAF > 0.625 & MAF < 0.875, "#FF2600", COL)) %>%
      mutate(COLOLCONTIG = ".") %>%
      # mutate(COLOLCONTIG = ifelse(CONTIG == "chr1" | CONTIG == "chr3" | CONTIG == "chr5" | CONTIG == "chr7" | CONTIG == "chr9" | CONTIG == "chr11" | CONTIG == "chr13" | CONTIG == "chr15" | CONTIG == "chr17" | CONTIG == "chr19" | CONTIG == "chr21" | CONTIG == "chrX", "#941100", "#011893"))
      mutate(COLOLCONTIG = ifelse(CONTIG == "chr1" | CONTIG == "chr3" | CONTIG == "chr5" | CONTIG == "chr7" | CONTIG == "chr9" | CONTIG == "chr11" | CONTIG == "chr13" | CONTIG == "chr15" | CONTIG == "chr17" | CONTIG == "chr19" | CONTIG == "chr21" | CONTIG == "chrX", "#FEB49B", "#9BC2E6"))
      
    ALL_FILTERED = ALL %>% sample_frac(0.5) %>% arrange(CONTIG, POSITION)
    
    #bedGraph file with denoised copy-number data
    # ALL_IGV <- readRDS(paste0(INDIR, i, "_allelicCounts", GNOMAD_PATTERN, "_", PON, ".rds"))
    ALL_IGV <- readRDS(paste0(INDIR, i, "_allelicCounts", GNOMAD_PATTERN, ".rds"))
    ALL_IGV = ALL_IGV %>%
      mutate(N = REF_COUNT + ALT_COUNT) %>%
      mutate(TEST = sample(1:2, nrow(ALL_IGV), replace = T)) %>%
      mutate(MAF = ifelse(TEST == 1, ALT_COUNT/N, REF_COUNT/N)) %>%
      mutate(POSITION = POSITION-1) %>%
      mutate(POSITION2 = POSITION) %>%
      select(CONTIG, POSITION, POSITION2, MAF) %>%
      arrange(CONTIG, POSITION)
    write_tsv(ALL_IGV, path = paste0(OUTDIR_IGV, GROUP, "_", i, "_SNP"), col_names = F)
    
    #transform to GRanges

    TRANS_TO_GRange = function(DATA) {
      DATA_temp = DATA %>% select(CONTIG, START, END)
      names(DATA_temp) = c("seqnames", "start", "end")
      DATA = cbind(DATA_temp, DATA) %>% toGRanges()
      }

    DENOIS = TRANS_TO_GRange(DENOIS)
    DENOISABN = TRANS_TO_GRange(DENOISABN)
    DENOIS_FILTERED = TRANS_TO_GRange(DENOIS_FILTERED)
    MODEL = TRANS_TO_GRange(MODEL)
    MODEL_LOH = TRANS_TO_GRange(MODEL_LOH)
    
    TRANS_TO_GRange2 = function(DATA) {
      DATA_temp = DATA %>% mutate(POSITION2 = POSITION) %>% select(CONTIG, POSITION, POSITION2)
      names(DATA_temp) = c("seqnames", "start", "end")
      DATA = cbind(DATA_temp, DATA) %>% toGRanges()
    }
    
    ALL = TRANS_TO_GRange2(ALL)
    ALL_FILTERED = TRANS_TO_GRange2(ALL_FILTERED)
    
    
    if(i == paste0(SAMPLE1_TYPE, "_", SAMPLE1_ID, "_", SAMPLE1_SEX)) {
      assign(paste0("SAMPLE1", "_DENOIS"), DENOIS)
      assign(paste0("SAMPLE1", "_DENOISABN"), DENOISABN)
      assign(paste0("SAMPLE1", "_DENOIS_FILTERED"), DENOIS_FILTERED)
      assign(paste0("SAMPLE1", "_ALL"), ALL)
      assign(paste0("SAMPLE1", "_ALL_FILTERED"), ALL_FILTERED)
      assign(paste0("SAMPLE1", "_MODEL"), MODEL)
      assign(paste0("SAMPLE1", "_MODEL_LOH"), MODEL_LOH)
      
    }
    
    if(i == paste0(SAMPLE2_TYPE, "_", SAMPLE2_ID, "_", SAMPLE2_SEX)) {
      assign(paste0("SAMPLE2", "_DENOIS"), DENOIS)
      assign(paste0("SAMPLE2", "_DENOISABN"), DENOISABN)
      assign(paste0("SAMPLE2", "_DENOIS_FILTERED"), DENOIS_FILTERED)
      assign(paste0("SAMPLE2", "_ALL"), ALL)
      assign(paste0("SAMPLE2", "_ALL_FILTERED"), ALL_FILTERED)
      assign(paste0("SAMPLE2", "_MODEL"), MODEL)
      assign(paste0("SAMPLE2", "_MODEL_LOH"), MODEL_LOH)
    }
    
    if(i == paste0(SAMPLE3_TYPE, "_", SAMPLE3_ID, "_", SAMPLE3_SEX)) {
      assign(paste0("SAMPLE3", "_DENOIS"), DENOIS)
      assign(paste0("SAMPLE3", "_DENOISABN"), DENOISABN)
      assign(paste0("SAMPLE3", "_DENOIS_FILTERED"), DENOIS_FILTERED)
      assign(paste0("SAMPLE3", "_ALL"), ALL)
      assign(paste0("SAMPLE3", "_ALL_FILTERED"), ALL_FILTERED)
      assign(paste0("SAMPLE3", "_MODEL"), MODEL)
      assign(paste0("SAMPLE3", "_MODEL_LOH"), MODEL_LOH)
    }
    

  }}
  

if(PON_STYLE == "M") {
  PON1 = "PoN-M"
  PON2 = PON1
  PON3 = PON1
  }
if(PON_STYLE == "F") {
  PON1 = "PoN-F"
  PON2 = PON1
  PON3 = PON1
  }
if(PON_STYLE == "mixed") {
  PON1 = "PoN-mixed"
  PON2 = PON1
  PON3 = PON1
  }
if(PON_STYLE == "none") {
  PON1 = "PoN-none"
  PON2 = PON1
  PON3 = PON1
  }
if(PON_STYLE == "matched_main") {
  PON1 = paste0("PoN-", SAMPLE1_SEX)
  PON2 = PON1
  PON3 = PON1  
  }
if(PON_STYLE == "matched_each") {
  PON1 = paste0("PoN-", SAMPLE1_SEX)
  PON2 = paste0("PoN-", SAMPLE2_SEX)
  PON3 = paste0("PoN-", SAMPLE3_SEX)
  }

#Define style
if(file.exists(paste0(INDIR, SAMPLE1_TYPE, "_", SAMPLE1_ID, "_", SAMPLE1_SEX, "_SEGMENTS_", PON1, "_", GETSEGMENTCONDITION, "_segmentation.tsv")) & 
   file.exists(paste0(INDIR, SAMPLE2_TYPE, "_", SAMPLE2_ID, "_", SAMPLE2_SEX, "_SEGMENTS_", PON2, "_", GETSEGMENTCONDITION, "_segmentation.tsv")) &
   file.exists(paste0(INDIR, SAMPLE3_TYPE, "_", SAMPLE3_ID, "_", SAMPLE3_SEX, "_SEGMENTS_", PON3, "_", GETSEGMENTCONDITION, "_segmentation.tsv"))) {PLOT_STYLE = "RUN_TRIO"}

if(file.exists(paste0(INDIR, SAMPLE1_TYPE, "_", SAMPLE1_ID, "_", SAMPLE1_SEX, "_SEGMENTS_", PON1, "_", GETSEGMENTCONDITION, "_segmentation.tsv")) & 
   file.exists(paste0(INDIR, SAMPLE2_TYPE, "_", SAMPLE2_ID, "_", SAMPLE2_SEX, "_SEGMENTS_", PON2, "_", GETSEGMENTCONDITION, "_segmentation.tsv")) &
   file.exists(paste0(INDIR, SAMPLE3_TYPE, "_", SAMPLE3_ID, "_", SAMPLE3_SEX, "_SEGMENTS_", PON3, "_", GETSEGMENTCONDITION, "_segmentation.tsv"))==F) {PLOT_STYLE = "RUN_DUO"}

if(file.exists(paste0(INDIR, SAMPLE1_TYPE, "_", SAMPLE1_ID, "_", SAMPLE1_SEX, "_SEGMENTS_", PON1, "_", GETSEGMENTCONDITION, "_segmentation.tsv")) & 
   file.exists(paste0(INDIR, SAMPLE2_TYPE, "_", SAMPLE2_ID, "_", SAMPLE2_SEX, "_SEGMENTS_", PON2, "_", GETSEGMENTCONDITION, "_segmentation.tsv"))==F &
   file.exists(paste0(INDIR, SAMPLE3_TYPE, "_", SAMPLE3_ID, "_", SAMPLE3_SEX, "_SEGMENTS_", PON3, "_", GETSEGMENTCONDITION, "_segmentation.tsv"))) {
  PLOT_STYLE = "RUN_DUO"
  SAMPLE2_ID=SAMPLE3_ID
  SAMPLE2_TYPE=SAMPLE3_TYPE
  SAMPLE2_SEX=SAMPLE3_SEX
  assign("SAMPLE2_DENOIS", SAMPLE3_DENOIS)
  assign("SAMPLE2_DENOISABN", SAMPLE3_DENOISABN)
  assign("SAMPLE2_DENOIS_FILTERED", SAMPLE3_DENOIS_FILTERED)
  assign("SAMPLE2_ALL", SAMPLE3_ALL)
  assign("SAMPLE2_ALL_FILTERED", SAMPLE3_ALL_FILTERED)
  assign("SAMPLE2_MODEL", SAMPLE3_MODEL)
  assign("SAMPLE2_MODEL_LOH", SAMPLE3_MODEL_LOH)
  remove(SAMPLE3_DENOIS, SAMPLE3_DENOISABN, SAMPLE3_ALL, SAMPLE3_ALL_FILTERED, SAMPLE3_MODEL, SAMPLE3_MODEL_LOH)
  }

if(file.exists(paste0(INDIR, SAMPLE1_TYPE, "_", SAMPLE1_ID, "_", SAMPLE1_SEX, "_SEGMENTS_", PON1, "_", GETSEGMENTCONDITION, "_segmentation.tsv")) & 
   file.exists(paste0(INDIR, SAMPLE2_TYPE, "_", SAMPLE2_ID, "_", SAMPLE2_SEX, "_SEGMENTS_", PON2, "_", GETSEGMENTCONDITION, "_segmentation.tsv"))==F &
   file.exists(paste0(INDIR, SAMPLE3_TYPE, "_", SAMPLE3_ID, "_", SAMPLE3_SEX, "_SEGMENTS_", PON3, "_", GETSEGMENTCONDITION, "_segmentation.tsv"))==F) {PLOT_STYLE = "RUN_SINGLE"}


message(paste0("   R ... ", date(), " - STAGE 3/5 - genome plot"))
#-----GENOME-VIEW------------------------------------------------------------------------------------------------------------------------------------
GENOME_CHART = function (ID, DATADENOIS, DATAMODEL, DATAALL, YMIN, YMAX, R0, R1, R0b, R1b) {
  kpAddLabels(kp, labels = c(ID), srt=90, pos=1, label.margin = 0.04, r0=R0, r1=R1b, family="Arial", cex=0.8)
  #COPY-NUMBER
  kpAddLabels(kp, labels = c("Log2R"), srt=90, pos=1, label.margin = 0.025, r0=R0, r1=R1, cex=0.5, family="Arial")
  kpPoints(kp, data=DATADENOIS, chr=DATADENOIS$CONTIG, x=DATADENOIS$START, y=DATADENOIS$L2R, cex=0.2, r0=R0, r1=R1, ymin=YMIN, ymax=YMAX, pch = 20, col = DATADENOIS$COLOLCONTIG)
  kpAxis(kp, ymin=YMIN, ymax=YMAX, r0=R0, r1=R1, col="black", cex=0.5, family="Arial", tick.pos=c(-2, -1, 0, 1, 2))
  kpAbline(kp, h=c(0), col="grey70", ymin=YMIN, ymax=YMAX, r0=R0, r1=R1, lwd=2, lty=1)
  kpSegments(kp, data=DATAMODEL, chr=DATAMODEL$CONTIG, x0=DATAMODEL$START, x1=DATAMODEL$END, y0=DATAMODEL$MEAN_L2R_EDIT, y1=DATAMODEL$MEAN_L2R_EDIT, col="black", r0=R0, r1=R1, ymin=YMIN, ymax=YMAX, lwd=2)
  #MAF
  kpAddLabels(kp, labels = c("AF"), srt=90, pos=1, label.margin = 0.025, r0=R0b, r1=R1b, cex=0.5, family="Arial")
  if(file.exists(GETCAPTURE)) {kpRect(kp, data=DTB_CAPTURE, chr=DTB_CAPTURE$CONTIG, x0=DTB_CAPTURE$STARText1, x1=DTB_CAPTURE$ENDext1, y0=0.45, y1=0.55, col=rgb(250, 240, 0, max=255, alpha=255), r0=R0b, r1=R1b, ymin=0, ymax=1, border=NA)}
  kpPoints(kp, data=DATAALL, chr=DATAALL$CONTIG, x=DATAALL$POSITION, y=DATAALL$MAF, cex=0.1, r0=R0b, r1=R1b, ymax=1, ymin=0, pch = 20, col = DATAALL$COLOLCONTIG)
  kpAxis(kp, ymin=0, ymax=1, r0=R0b, r1=R1b, col="black", cex=0.5, family="Arial", tick.pos=c(0, 0.5, 1))
}

output = paste0(OUTDIR, GROUP, "_genome_profile_", SAMPLE1_TYPE, "_", SAMPLE1_ID, "_", SAMPLE1_SEX)

if(PLOT_STYLE == "RUN_TRIO")   {png(paste0(output, ".png"), width = 4000, height = 2500, res = 300)}
if(PLOT_STYLE == "RUN_DUO")    {png(paste0(output, ".png"), width = 4000, height = 1500, res = 300)}
if(PLOT_STYLE == "RUN_SINGLE") {png(paste0(output, ".png"), width = 4000, height = 1000, res = 300)}

#Create karyotype
plot.params <- getDefaultPlotParams(plot.type=4)
kp <- plotKaryotype(genome=ASSEMBLY,
                    plot.type = 4,
                    plot.params = plot.params,
                    ideogram.plotter = NULL,
                    labels.plotter = NULL,
                    main=GROUP, family = "Arial", cex=1.5, font=2)
kpAddCytobandsAsLine(kp, color.schema = "circos", lwd=7)
kpAddChromosomeNames(kp, srt=45, family="Arial")

if(PLOT_STYLE == "RUN_TRIO") {
  kpDataBackground(kp, r0=0, r1=0.2, color = "grey95")
  kpDataBackground(kp, r0=0.3, r1=0.5, color = "grey95")
  kpDataBackground(kp, r0=0.6, r1=0.8, color = "grey95")
  SAMPLE1_GENOME = GENOME_CHART(paste0(SAMPLE1_TYPE, " ", SAMPLE1_ID), SAMPLE1_DENOIS_FILTERED, SAMPLE1_MODEL, SAMPLE1_ALL_FILTERED, -2.5, 2.5, 0.0, 0.2, 0.21, 0.29)
  SAMPLE2_GENOME = GENOME_CHART(paste0(SAMPLE2_TYPE, " ", SAMPLE2_ID), SAMPLE2_DENOIS_FILTERED, SAMPLE2_MODEL, SAMPLE2_ALL_FILTERED, -2.5, 2.5, 0.3, 0.5, 0.51, 0.59)
  SAMPLE3_GENOME = GENOME_CHART(paste0(SAMPLE3_TYPE, " ", SAMPLE3_ID), SAMPLE3_DENOIS_FILTERED, SAMPLE3_MODEL, SAMPLE3_ALL_FILTERED, -2.5, 2.5, 0.6, 0.8, 0.81, 0.89)}
if(PLOT_STYLE == "RUN_DUO") {
  kpDataBackground(kp, r0=0, r1=0.3, color = "grey95")
  kpDataBackground(kp, r0=0.5, r1=0.8, color = "grey95")
  SAMPLE1_GENOME = GENOME_CHART(paste0(SAMPLE1_TYPE, " ", SAMPLE1_ID), SAMPLE1_DENOIS_FILTERED, SAMPLE1_MODEL, SAMPLE1_ALL_FILTERED, -2.5, 2.5, 0.0, 0.3, 0.33, 0.45)
  SAMPLE2_GENOME = GENOME_CHART(paste0(SAMPLE2_TYPE, " ", SAMPLE2_ID), SAMPLE2_DENOIS_FILTERED, SAMPLE2_MODEL, SAMPLE2_ALL_FILTERED, -2.5, 2.5, 0.5, 0.8, 0.83, 0.95)}
if(PLOT_STYLE == "RUN_SINGLE") {
  kpDataBackground(kp, r0=0, r1=0.6, color = "grey95")
  SAMPLE1_GENOME = GENOME_CHART(paste0(SAMPLE1_TYPE, " ", SAMPLE1_ID), SAMPLE1_DENOIS_FILTERED, SAMPLE1_MODEL, SAMPLE1_ALL_FILTERED, -2.5, 2.5, 0.0, 0.6, 0.65, 0.9)}
invisible(dev.off())
#-----------------------------------------------------------------------------------------------------------------------------------------


message(paste0("   R ... ", date(), " - STAGE 4/5 - chromosome plot"))

#-----CHROMOSOME-VIEW------------------------------------------------------------------------------------------------------------------------------------
# for (CHR in list("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")) {

#chromosomes from capture file
if(file.exists(GETCAPTURE)) {
  LISTofCHR = unique(DTB_CAPTURE$CONTIG) %>% as.list()
}
#chromosomes for WGS
if(file.exists(GETCAPTURE)==F) {
  LISTofCHR = list("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
}

for (CHR in (LISTofCHR)) {
  
message(paste0("   R ... ", date(), " - STAGE 4/5 - chromosome plot: ", CHR))

output = paste0(OUTDIR, GROUP, "_", CHR, "_", SAMPLE1_TYPE, "_", SAMPLE1_ID, "_", SAMPLE1_SEX)

RES = 300
MIN = vlookup("chr21", dict = DTB_CHRSIZE, result_column = 2, lookup_column = 1)
SIZE = vlookup(CHR, dict = DTB_CHRSIZE, result_column = 2, lookup_column = 1)
WIDTH_DATA = (3.5*RES)*(SIZE/MIN) #with of the data for required chromosome
MARGIN_LEFT_PIX = RES/3
MARGIN_RIGHT_PIX = RES/5
WIDTH=WIDTH_DATA+MARGIN_LEFT_PIX+MARGIN_RIGHT_PIX
# margin_pixel = (x*WIDTH_ORIG)+(0.01*RES)
# x = (margin_pixel - (0.01*RES))/WIDTH_ORIG
marginleft = (MARGIN_LEFT_PIX - (0.01*RES))/WIDTH
marginright = (MARGIN_RIGHT_PIX - (0.01*RES))/WIDTH
IDPOS = SIZE/2


#Create chromosome picture
CHROMOSOME_CHART = function (ID, DATADENOIS, DATADENOISABN, DATAMODEL, DATAMODEL_LOH, DATAALL, YMIN, YMAX, R0, R1, R0b, R1b) {
  kpText(kp, chr=CHR, ymin=0, ymax=1, x=IDPOS, y=1, labels=ID, col="black", pos=3, r0=R0b, r1=R1b, cex=0.5, family="Arial")

  #COPY-NUMBER
  if(nrow(DATAMODEL_LOH %>% as.data.frame())>0) {kpRect(kp, data=DATAMODEL_LOH, chr=DATAMODEL_LOH$CONTIG, x0=DATAMODEL_LOH$START, x1=DATAMODEL_LOH$END, y0=-2.1, y1=2.1, col=rgb(142, 169, 219, max=255, alpha=70), r0=R0, r1=R1, ymin=YMIN, ymax=YMAX, border=NA)}
  kpText(kp, chr=CHR, ymin=YMIN, ymax=YMAX, x=0, y=YMAX-0.5, labels="Log2R", col="black", pos=4, r0=R0, r1=R1, cex=0.5, family="Arial")
  kpAbline(kp, h=c(-2, -1, 0, 1, 2), col="grey70", ymin=YMIN, ymax=YMAX, r0=R0, r1=R1, lwd=0.5, lty=2)
  #note that kpPoints for DATADENOIS and DATADENOISABN should be perhaps kpSegments as it is in detailed chart, not probles as long as capture regions are small
  kpPoints(kp, data=DATADENOIS, chr=DATADENOIS$CONTIG, x=DATADENOIS$START, y=DATADENOIS$L2R, cex=0.5, r0=R0, r1=R1, ymin=YMIN, ymax=YMAX, pch=20, col=DATADENOIS$COLORABN)
  if(nrow(DATADENOISABN %>% as.data.frame())>0) {kpPoints(kp, data=DATADENOISABN, chr=DATADENOISABN$CONTIG, x=DATADENOISABN$START, y=DATADENOISABN$L2R, cex=0.5, r0=R0, r1=R1, ymin=YMIN, ymax=YMAX, pch=20, col=DATADENOISABN$COLORABN)}
  kpAxis(kp, ymin=YMIN, ymax=YMAX, r0=R0, r1=R1, col="black", cex=0.5, family="Arial", tick.pos=c(-2, -1, 0, 1, 2), tick.len=0)
  kpSegments(kp, data=DATAMODEL, chr=DATAMODEL$CONTIG, x0=DATAMODEL$START, x1=DATAMODEL$END, y0=DATAMODEL$MEAN_L2R_EDIT, y1=DATAMODEL$MEAN_L2R_EDIT, col="black", r0=R0, r1=R1, ymin=YMIN, ymax=YMAX, lwd=1.7)
  kpPoints(kp, data=DTB_CENTROMERE, chr=DTB_CENTROMERE$CONTIG, x=DTB_CENTROMERE$POSITION, y=0, cex=2, r0=R0, r1=R1, ymin=YMIN, ymax=YMAX, pch=20, col="black")

  #MAF
  kpText(kp, chr=CHR, ymin=0, ymax=1, x=0, y=1.1, labels="AF", col="black", pos=4, r0=R0b, r1=R1b, cex=0.5, family="Arial")
  kpAbline(kp, h=c(0, 0.5, 1), col="grey70", ymin=0, ymax=1, r0=R0b, r1=R1b, lwd=0.5, lty=2)
  if(file.exists(GETCAPTURE)) {kpRect(kp, data=DTB_CAPTURE, chr=DTB_CAPTURE$CONTIG, x0=DTB_CAPTURE$STARText2, x1=DTB_CAPTURE$ENDext2, y0=0, y1=1, col=rgb(255, 235, 120, max=255, alpha=120), r0=R0b, r1=R1b, ymin=0, ymax=1, border=NA)}
  kpPoints(kp, data=DATAALL, chr=DATAALL$CONTIG, x=DATAALL$POSITION, y=DATAALL$MAF, cex=0.4, r0=R0b, r1=R1b, ymax=1, ymin=0, pch=20, col=DATAALL$COL)
  kpAxis(kp, ymin=0, ymax=1, r0=R0b, r1=R1b, col="black", cex=0.5, family="Arial", tick.pos=c(0, 0.5, 1), tick.len=0)
}

if(PLOT_STYLE == "RUN_TRIO")   {png(paste0(output, ".png"), width=WIDTH, height=2500, res=RES)}
if(PLOT_STYLE == "RUN_DUO")    {png(paste0(output, ".png"), width=WIDTH, height=2000, res=RES)}
if(PLOT_STYLE == "RUN_SINGLE") {png(paste0(output, ".png"), width=WIDTH, height=1000, res=RES)}

plot.params <- getDefaultPlotParams(plot.type=4)
plot.params$leftmargin=marginleft
plot.params$rightmargin=marginright
if(PLOT_STYLE == "RUN_SINGLE") {plot.params$ideogramheight=15}
kp <- plotKaryotype(genome=ASSEMBLY,
                    chromosomes=CHR,
                    plot.type=1,
                    plot.params=plot.params,
                    main=GROUP, family="Arial", cex=1.5, font=2,
                    labels.plotter = NULL) # no chromosomes numbers on left, add manually below

kpAddBaseNumbers(kp, tick.dist=10000000, tick.len=5, cex=0.7,
                      minor.tick.dist=1000000, minor.tick.len=2.5, family="Arial")
kpAddCytobandLabels(kp, cex=0.5, family="Arial")
kpAddChromosomeNames(kp, xoffset = -28, yoffset=0.5, srt=0, family="Arial", font=2)


if(PLOT_STYLE == "RUN_TRIO") {
  SAMPLE1_CHROM = CHROMOSOME_CHART(paste0(SAMPLE1_TYPE, " ", SAMPLE1_ID), SAMPLE1_DENOIS, SAMPLE1_DENOISABN, SAMPLE1_MODEL, SAMPLE1_MODEL_LOH, SAMPLE1_ALL, -2.5, 2.5, 0.0, 0.2, 0.21, 0.29)
  SAMPLE2_CHROM = CHROMOSOME_CHART(paste0(SAMPLE2_TYPE, " ", SAMPLE2_ID), SAMPLE2_DENOIS, SAMPLE2_DENOISABN, SAMPLE2_MODEL, SAMPLE2_MODEL_LOH, SAMPLE2_ALL, -2.5, 2.5, 0.3, 0.5, 0.51, 0.59)
  SAMPLE3_CHROM = CHROMOSOME_CHART(paste0(SAMPLE3_TYPE, " ", SAMPLE3_ID), SAMPLE3_DENOIS, SAMPLE3_DENOISABN, SAMPLE3_MODEL, SAMPLE3_MODEL_LOH, SAMPLE3_ALL, -2.5, 2.5, 0.6, 0.8, 0.81, 0.89)}
if(PLOT_STYLE == "RUN_DUO") {
  SAMPLE1_CHROM = CHROMOSOME_CHART(paste0(SAMPLE1_TYPE, " ", SAMPLE1_ID), SAMPLE1_DENOIS, SAMPLE1_DENOISABN, SAMPLE1_MODEL, SAMPLE1_MODEL_LOH, SAMPLE1_ALL, -2.5, 2.5, 0.0, 0.3, 0.33, 0.45)
  SAMPLE2_CHROM = CHROMOSOME_CHART(paste0(SAMPLE2_TYPE, " ", SAMPLE2_ID), SAMPLE2_DENOIS, SAMPLE2_DENOISABN, SAMPLE2_MODEL, SAMPLE2_MODEL_LOH, SAMPLE2_ALL,  -2.5, 2.5, 0.5, 0.8, 0.83, 0.95)}
if(PLOT_STYLE == "RUN_SINGLE") {
  SAMPLE1_CHROM = CHROMOSOME_CHART(paste0(SAMPLE1_TYPE, " ", SAMPLE1_ID), SAMPLE1_DENOIS, SAMPLE1_DENOISABN, SAMPLE1_MODEL, SAMPLE1_MODEL_LOH, SAMPLE1_ALL, -2.5, 2.5, 0.0, 0.6, 0.65, 0.9)}

invisible(dev.off())
}
#-----------------------------------------------------------------------------------------------------------------------------------------



message(paste0("   R ... ", date(), " - STAGE 5/5 - abnormality plot"))

defaultW <- getOption("warn")
options(warn = -1)

#------DETAIL-VIEW-----------------------------------------------------------------------------------------------------------------------------------
DETAILS = as.data.frame(SAMPLE1_MODEL) %>% filter(TYPE == "LOSS" | TYPE == "GAIN") #%>% filter(PROBES > 2)

if(nrow(SAMPLE1_MODEL_LOH %>% as.data.frame())>0) {
DETAILS = rbind(DETAILS, as.data.frame(SAMPLE1_MODEL_LOH)) #%>% filter(CONTIG != "chrZ") #filter(TYPE !="x")
}

if(nrow(DETAILS)>0) {
#print only PASS segments
# DETAILS = DETAILS %>% filter(FILTER == "PASS")

DETAILS = mutate(DETAILS, N = 1:n())

LISTofABNORM = unique(DETAILS$N) %>% as.list()
NofABN = max(unlist(LISTofABNORM))

for (ABN in (LISTofABNORM)) {

  message(paste0("   R ... ", date(), " - STAGE 5/5 - abnormality plot: ", ABN, "/", NofABN))

  DETAILS_PART = DETAILS %>%
    filter(N == ABN)

  FROMABN = DETAILS_PART$START
  TOABN = DETAILS_PART$END
  SIZEABN = TOABN-FROMABN+1
  TYPEABN = DETAILS_PART$TYPE
  PROBE = DETAILS_PART$PROBES %>% as.character()
  CHRBAND = DETAILS_PART$CHRBAND

  CHR = DETAILS_PART$CONTIG %>% as.character()
  FROM = ifelse(TYPEABN == "LOSS" | TYPEABN == "GAIN", (DETAILS_PART$START - DETAILS_PART$SIZE*10), (DETAILS_PART$START - DETAILS_PART$SIZE*0.5))
  if (FROM < 0) {FROM = 0}
  TO = ifelse(TYPEABN == "LOSS" | TYPEABN == "GAIN", (DETAILS_PART$END + DETAILS_PART$SIZE*10), (DETAILS_PART$END + DETAILS_PART$SIZE*0.5))
  SIZE = TO-FROM

  output = paste0(OUTDIR, "detail/", GROUP, "_", SAMPLE1_TYPE, "_", SAMPLE1_ID, "_", SAMPLE1_SEX, "_", CHR, "-", FROMABN, "-", TOABN, "_", TYPEABN, "_size-", SIZEABN, "bp_", PROBE, "-probes")

  ZOOM = toGRanges(data.frame(CHR, FROM, TO))
  IDPOS = FROM+((TO-FROM)/2)

  #Data function
  DETAIL_CHART = function (ID, DATADENOIS, DATADENOISABN, DATAMODEL, DATAMODEL_LOH, DATAALL, YMIN, YMAX, R0, R1, R0b, R1b) {
    kpText(kp, chr=CHR, ymin=0, ymax=1, x=IDPOS, y=0.95, labels=ID, col="black", pos=3, r0=R0b, r1=R1b, cex=0.5, family="Arial")

    #COPY-NUMBER
    kpRect(kp, chr=CHR, x0=FROM, x1=TO, y0=YMIN+0.5, y1=YMAX-0.5, col="grey95", r0=R0, r1=R1, ymin=YMIN, ymax=YMAX, border=NA)
    if(nrow(DATAMODEL_LOH %>% as.data.frame())>0) {kpRect(kp, data=DATAMODEL_LOH, chr=DATAMODEL_LOH$CONTIG, x0=DATAMODEL_LOH$START, x1=DATAMODEL_LOH$END, y0=-2.1, y1=2.1, col=rgb(142, 169, 219, max=255, alpha=70), r0=R0, r1=R1, ymin=YMIN, ymax=YMAX, border=NA)}
    kpText(kp, chr=CHR, ymin=YMIN, ymax=YMAX, x=FROM, y=YMAX, labels="Log2R", col="black", pos=4, r0=R0, r1=R1, cex=0.5, family="Arial")
    kpAbline(kp, h=c(-2, -1, 1, 2), col="grey70", ymin=YMIN, ymax=YMAX, r0=R0, r1=R1, lwd=0.5, lty=2)
    kpAbline(kp, h=c(0), col="grey10", ymin=YMIN, ymax=YMAX, r0=R0, r1=R1, lwd=0.9, lty=2)
    kpSegments(kp, chr=CHR, x0=FROM, x1=FROM, y0=YMIN, y1=YMAX, col="black", r0=R0, r1=R1, ymin=YMIN, ymax=YMAX, lwd=1.7)
    kpText(kp, chr=CHR, ymin=YMIN, ymax=YMAX, x=c(FROM, FROM, FROM, FROM, FROM), y=c(-2, -1, 0, 1, 2), labels=c(-2, -1, 0, 1, 2), col="black", pos=2, r0=R0, r1=R1, cex=0.5, family="Arial", clipping = FALSE)
    kpSegments(kp, data=DATADENOIS, chr=DATADENOIS$CONTIG, x0=DATADENOIS$START, x1=DATADENOIS$END, y0=DATADENOIS$L2R, y1=DATADENOIS$L2R, r0=R0, r1=R1, ymin=YMIN, ymax=YMAX, col=DATADENOIS$COLORABN, lwd=7)
    if(nrow(DATADENOISABN %>% as.data.frame())>0) {kpSegments(kp, data=DATADENOISABN, chr=DATADENOISABN$CONTIG, x0=DATADENOISABN$START, x1=DATADENOISABN$END, y0=DATADENOISABN$L2R, y1=DATADENOISABN$L2R, r0=R0, r1=R1, ymin=YMIN, ymax=YMAX, col=DATADENOISABN$COLORABN, lwd=7)}
    kpSegments(kp, data=DATAMODEL, chr=DATAMODEL$CONTIG, x0=DATAMODEL$START, x1=DATAMODEL$END, y0=DATAMODEL$MEAN_L2R_EDIT, y1=DATAMODEL$MEAN_L2R_EDIT, col="black", r0=R0, r1=R1, ymin=YMIN, ymax=YMAX, lwd=1.7)
    kpPoints(kp, data=DTB_CENTROMERE, chr=DTB_CENTROMERE$CONTIG, x=DTB_CENTROMERE$POSITION, y=0, cex=2, r0=R0, r1=R1, ymin=YMIN, ymax=YMAX, pch=20, col="black")

    #MAF
    kpText(kp, chr=CHR, ymin=0, ymax=1, x=FROM, y=1.1, labels="AF", col="black", pos=4, r0=R0b, r1=R1b, cex=0.5, family="Arial")
    kpAbline(kp, h=c(0, 0.5, 1), col="grey70", ymin=0, ymax=1, r0=R0b, r1=R1b, lwd=0.5, lty=2)
    kpSegments(kp, chr=CHR, x0=FROM, x1=FROM, y0=0, y1=1, col="black", r0=R0b, r1=R1b, ymin=0, ymax=1, lwd=1.7)
    kpText(kp, chr=CHR, ymin=0, ymax=1, x=c(FROM, FROM, FROM), y=c(0, 0.5, 1), labels=c(0, 0.5, 1), col="black", pos=2, r0=R0b, r1=R1b, cex=0.5, family="Arial", clipping = FALSE)
    if(file.exists(GETCAPTURE)) {kpRect(kp, data=DTB_CAPTURE, chr=DTB_CAPTURE$CONTIG, x0=DTB_CAPTURE$STARText3, x1=DTB_CAPTURE$ENDext3, y0=0, y1=1, col=rgb(255, 235, 120, max=255, alpha=100), r0=R0b, r1=R1b, ymin=0, ymax=1, border=NA)}
    kpPoints(kp, data=DATAALL, chr=DATAALL$CONTIG, x=DATAALL$POSITION, y=DATAALL$MAF, cex=0.4, r0=R0b, r1=R1b, ymax=1, ymin=0, pch=20, col=DATAALL$COL)

  }

  #DGV function
  DGV_CHART = function(DATA, POS1, POS2, COLOR) {
    #DATA PANEL 2 - DGV
    kpRect(kp, data.panel=2, chr=DATA$CONTIG, x0=DATA$START, x1=DATA$END, y0=POS1, y1=POS2, r0=0.25, r1=0.35, col=COLOR, border=NA)
  }

  # #gnomAD function
  # GNOMADSV_CHART = function(DATA, POS1, POS2, COLOR) {
  #   #DATA PANEL 2 - GNOMADSV
  #   kpRect(kp, data.panel=2, chr=DATA$CONTIG, x0=DATA$START, x1=DATA$END, y0=POS1, y1=POS2, r0=0.35, r1=0.45, col=COLOR, border=NA)
  # }

  #Prepare gene database
  #------------------------------------------------------------------------------------------------------------------------------------------------
  GENES_DTB_PROCESSED_REGION = DTB_GENES_DTB_PROCESSED %>% filter(CONTIG == CHR) %>% filter(TxEND >= FROM) %>% filter(TxSTART <= TO) %>%
    mutate(TxCenter = (TxSTART+TxEND)/2) %>%
    mutate(WAY = ifelse(STRAND=="+", TxSTART, TxEND))

    if(nrow(GENES_DTB_PROCESSED_REGION)>0) {
      GENES_DTB_PROCESSED_REGION = GENES_DTB_PROCESSED_REGION %>%
      mutate(N = (1:n()))

    #create bins, position for each gene
    GENES_DTB_PROCESSED_BINS = GENES_DTB_PROCESSED_REGION %>%
      mutate(TxSTART = TxSTART - SIZE*0.004) %>% #test if this is needed, what it does
      mutate(TxEND = TxEND + SIZE*0.004) %>% #test if this is needed, what it does
      select(CONTIG, TxSTART, TxEND) %>% toGRanges() %>% disjointBins() %>% as.data.frame()

    names(GENES_DTB_PROCESSED_BINS) = c("BIN")
    GENES_DTB_PROCESSED_NofLayers = max(GENES_DTB_PROCESSED_BINS)
    if(GENES_DTB_PROCESSED_NofLayers <= 20) {GENES_DTB_PROCESSED_LayerMargin = 0.01}
    if(GENES_DTB_PROCESSED_NofLayers > 20) {GENES_DTB_PROCESSED_LayerMargin = (20*0.01)/GENES_DTB_PROCESSED_NofLayers}
    GENES_DTB_PROCESSED_LayerHeight = (1-((GENES_DTB_PROCESSED_NofLayers-1)*GENES_DTB_PROCESSED_LayerMargin))/GENES_DTB_PROCESSED_NofLayers

    y0 = (GENES_DTB_PROCESSED_LayerHeight+GENES_DTB_PROCESSED_LayerMargin)*(GENES_DTB_PROCESSED_BINS-1)
    names(y0) = c("Y0")
    y1 = GENES_DTB_PROCESSED_LayerHeight*GENES_DTB_PROCESSED_BINS+GENES_DTB_PROCESSED_LayerMargin*(GENES_DTB_PROCESSED_BINS-1)
    names(y1) = c("Y1")
    GENES_DTB_PROCESSED_BINS = cbind(GENES_DTB_PROCESSED_BINS, y0, y1)
    GENES_DTB_PROCESSED_BINS = GENES_DTB_PROCESSED_BINS %>%
      mutate(CENTER = Y0+((Y1-Y0)/2)) %>%
      mutate(EXON_Y0 = CENTER-(0.05*GENES_DTB_PROCESSED_LayerHeight)) %>%
      mutate(EXON_Y1 = CENTER+(0.05*GENES_DTB_PROCESSED_LayerHeight)) %>%
      mutate(EXON_CODING_Y0 = CENTER-(0.1*GENES_DTB_PROCESSED_LayerHeight)) %>%
      mutate(EXON_CODING_Y1 = CENTER+(0.1*GENES_DTB_PROCESSED_LayerHeight)) %>%
      mutate(EXON_ID_Y = Y1)

    GENES_DTB_PROCESSED_REGION = cbind(GENES_DTB_PROCESSED_BINS, GENES_DTB_PROCESSED_REGION)

    #extract exons and coding exons
    EXONS = data.frame(NA,NA,NA,NA,NA)
    names(EXONS) = c("CONTIG", "EXON_START", "EXON_END", "EXON_Y0", "EXON_Y1")
    EXONS_CODING = data.frame(NA,NA,NA,NA,NA)
    names(EXONS_CODING) = c("CONTIG", "EXON_CODING_START", "EXON_CODING_END", "EXON_CODING_Y0", "EXON_CODING_Y1")

    for (line in 1:nrow(GENES_DTB_PROCESSED_REGION)) {
      CONTIG = GENES_DTB_PROCESSED_REGION[line,11]
      EXON_Y0 = GENES_DTB_PROCESSED_REGION[line,5]
      EXON_Y1 = GENES_DTB_PROCESSED_REGION[line,6]
      EXON_CODING_Y0 = GENES_DTB_PROCESSED_REGION[line,7]
      EXON_CODING_Y1 = GENES_DTB_PROCESSED_REGION[line,8]

      #overal exons
      STARTS = GENES_DTB_PROCESSED_REGION[line,16]
      STARTS = unlist(strsplit(STARTS, ",")) %>% as.data.frame() %>% mutate_if(is.factor, as.character) %>% mutate_if(is.character, as.integer)
      ENDS = GENES_DTB_PROCESSED_REGION[line,17]
      ENDS = unlist(strsplit(ENDS, ",")) %>% as.data.frame() %>% mutate_if(is.factor, as.character) %>% mutate_if(is.character, as.integer)
      POSITIONS = cbind(STARTS, ENDS)
      names(POSITIONS) = c("EXON_START", "EXON_END")
      POSITIONS = POSITIONS %>%
        mutate(CONTIG = CONTIG, EXON_Y0 = EXON_Y0, EXON_Y1 = EXON_Y1) %>%
        select(CONTIG, EXON_START, EXON_END, EXON_Y0, EXON_Y1)
      EXONS = rbind(EXONS, POSITIONS)

      #coding exons, only if they exist
      NofCodingExons = GENES_DTB_PROCESSED_REGION[line,18] %>% as.integer()
      if (NofCodingExons >0) {
        STARTS = GENES_DTB_PROCESSED_REGION[line,19]
        STARTS = unlist(strsplit(STARTS, ",")) %>% as.data.frame() %>% mutate_if(is.factor, as.character) %>% mutate_if(is.character, as.integer)
        ENDS = GENES_DTB_PROCESSED_REGION[line,20]
        ENDS = unlist(strsplit(ENDS, ",")) %>% as.data.frame() %>% mutate_if(is.factor, as.character) %>% mutate_if(is.character, as.integer)
        POSITIONS = cbind(STARTS, ENDS)
        names(POSITIONS) = c("EXON_CODING_START", "EXON_CODING_END")
        POSITIONS = POSITIONS %>%
          mutate(CONTIG = CONTIG, EXON_CODING_Y0 = EXON_CODING_Y0, EXON_CODING_Y1 = EXON_CODING_Y1) %>%
          select(CONTIG, EXON_CODING_START, EXON_CODING_END, EXON_CODING_Y0, EXON_CODING_Y1)
        EXONS_CODING = rbind(EXONS_CODING, POSITIONS)
      }
    }

    EXONS = EXONS %>% filter(CONTIG != "NA")
    EXONS_CODING = EXONS_CODING %>% filter(CONTIG != "NA")

    #create position for arrows
    SLICE = as.integer(SIZE*2/100)
    GENES_DTB_PROCESSED_ARROWS = GENES_DTB_PROCESSED_REGION %>% select(GENE_ID, STRAND, CENTER, CONTIG, TxSTART, TxEND) %>%
      mutate(SIZE = TxEND - TxSTART) %>%
      mutate(NofARROWS = ceiling(SIZE/SLICE)) %>%
      expandRows("NofARROWS") %>%
      group_by(GENE_ID) %>% mutate(N = row_number()-1) %>%
      mutate(ARROW_START = ifelse(STRAND == "+", TxSTART + (N*SLICE), NA)) %>%
      mutate(ARROW_END = ifelse(STRAND == "+", ARROW_START+SLICE, NA)) %>%
      mutate(ARROW_END = ifelse(STRAND == "+" & ARROW_END > TxEND, TxEND, ARROW_END)) %>%
      mutate(ARROW_SIZE = ifelse(STRAND == "+", ARROW_END - ARROW_START, NA)) %>%
      mutate(ARROW_START = ifelse(STRAND == "-", TxEND - (N*SLICE), ARROW_START)) %>%
      mutate(ARROW_END = ifelse(STRAND == "-", ARROW_START-SLICE, ARROW_END)) %>%
      mutate(ARROW_END = ifelse(STRAND == "-" & ARROW_END < TxSTART, TxSTART, ARROW_END)) %>%
      mutate(ARROW_SIZE = ifelse(STRAND == "-", ARROW_START - ARROW_END, ARROW_SIZE)) %>%
      mutate(FILTER = ifelse(N != 0 & ARROW_SIZE < SLICE, "OUT", ".")) %>%
      filter(FILTER != "OUT")

    } else {
    GENES_DTB_PROCESSED_ARROWS = data.frame()
    EXONS = data.frame()
    EXONS_CODING = data.frame()
  }

  #------------------------------------------------------------------------------------------------------------------------------------------------


  if(PLOT_STYLE == "RUN_TRIO")   {png(paste0(output, ".png"), width = 5000, height = 3000, res = 300)}
  if(PLOT_STYLE == "RUN_DUO")    {png(paste0(output, ".png"), width = 5000, height = 2000, res = 300)}
  if(PLOT_STYLE == "RUN_SINGLE") {png(paste0(output, ".png"), width = 5000, height = 1500, res = 300)}

  plot.params <- getDefaultPlotParams(plot.type=3)
  plot.params$ideogramheight=15
  if(PLOT_STYLE == "RUN_TRIO") {
    plot.params$data1height=270
    plot.params$data2height=130}
  if(PLOT_STYLE == "RUN_DUO") {
    plot.params$data1height=230
    plot.params$data2height=170}
  if(PLOT_STYLE == "RUN_SINGLE") {
    plot.params$data1height=180
    plot.params$data2height=220}

  kp <- plotKaryotype(genome=ASSEMBLY,
                      chromosomes = CHR,
                      plot.type = 2,
                      plot.params = plot.params,
                      # ideogram.plotter = NULL,
                      labels.plotter = NULL,
                      main=paste0(GROUP, "  ●  ", TYPEABN, " ● ", CHRBAND, " ● ", CHR, ":", FROMABN, "-", TOABN, " ● ", SIZEABN, "bp"),
                      family = "Arial", cex=1,
                      zoom=ZOOM)

  # kpDataBackground(kp, data.panel = 2, col="white")

  TICK = 10000000
  TICK = ifelse(SIZE < 10000000, 1000000, TICK)
  TICK = ifelse(SIZE < 1000000, 100000, TICK)
  TICK = ifelse(SIZE < 100000, 10000, TICK)
  TICK = ifelse(SIZE < 10000, 1000, TICK)
  TICK = ifelse(SIZE < 1000, 100, TICK)
  TICK = ifelse(SIZE < 100, 10, TICK)
  TICK = ifelse(SIZE < 10, 1, TICK)

  kpAddBaseNumbers(kp, tick.dist=TICK, tick.len=5, cex=0.7,
                   minor.tick.dist = TICK/10, minor.tick.len=2.5,
                   family="Arial")
  kpAddCytobandLabels(kp, cex=0.7, family="Arial", font=2)

  if(PLOT_STYLE == "RUN_TRIO") {
    SAMPLE1_GENE = DETAIL_CHART(paste0(SAMPLE1_TYPE, " ", SAMPLE1_ID), SAMPLE1_DENOIS, SAMPLE1_DENOISABN, SAMPLE1_MODEL, SAMPLE1_MODEL_LOH, SAMPLE1_ALL, -2.5, 2.5, 0.00, 0.18, 0.20, 0.28)
    SAMPLE2_GENE = DETAIL_CHART(paste0(SAMPLE2_TYPE, " ", SAMPLE2_ID), SAMPLE2_DENOIS, SAMPLE2_DENOISABN, SAMPLE2_MODEL, SAMPLE2_MODEL_LOH, SAMPLE2_ALL, -2.5, 2.5, 0.32, 0.50, 0.52, 0.60)
    SAMPLE3_GENE = DETAIL_CHART(paste0(SAMPLE3_TYPE, " ", SAMPLE3_ID), SAMPLE3_DENOIS, SAMPLE3_DENOISABN, SAMPLE3_MODEL, SAMPLE3_MODEL_LOH, SAMPLE3_ALL, -2.5, 2.5, 0.64, 0.82, 0.84, 0.92)}
  if(PLOT_STYLE == "RUN_DUO") {
    CHILD_GENE   = DETAIL_CHART(paste0(SAMPLE1_TYPE, " ", SAMPLE1_ID), SAMPLE1_DENOIS, SAMPLE1_DENOISABN, SAMPLE1_MODEL, SAMPLE1_MODEL_LOH, SAMPLE1_ALL, -2.5, 2.5, 0.0, 0.3, 0.33, 0.45)
    FATHER_GENE  = DETAIL_CHART(paste0(SAMPLE2_TYPE, " ", SAMPLE2_ID), SAMPLE2_DENOIS, SAMPLE2_DENOISABN, SAMPLE2_MODEL, SAMPLE2_MODEL_LOH, SAMPLE2_ALL, -2.5, 2.5, 0.5, 0.8, 0.83, 0.95)}
  if(PLOT_STYLE == "RUN_SINGLE") {
    CHILD_GENE   = DETAIL_CHART(paste0(SAMPLE1_TYPE, " ", SAMPLE1_ID), SAMPLE1_DENOIS, SAMPLE1_DENOISABN, SAMPLE1_MODEL, SAMPLE1_MODEL_LOH, SAMPLE1_ALL, -2.5, 2.5, 0.00, 0.6, 0.65, 0.9)}


  #BUILD DGV DATABASE
  GAINCHART = DGV_CHART(DTB_DGV_GAIN, 0.01, 0.3, rgb(84, 130, 53, max = 255, alpha = 30))
  LOSSCHART = DGV_CHART(DTB_DGV_LOSS, 0.31, 0.6, rgb(192, 0, 0, max = 255, alpha = 30))
  COMPLEXCHART = DGV_CHART(DTB_DGV_COMPLEX, 0.61, 0.9, rgb(123, 35, 170, max = 255, alpha = 30))

  # #BUILD gnomAD SV DATABASE
  # GAINCHART = GNOMADSV_CHART(DTB_GNOMADSV_GAIN, 0.01, 0.3, rgb(84, 130, 53, max = 255, alpha = 50))
  # LOSSCHART = GNOMADSV_CHART(DTB_GNOMADSV_LOSS, 0.31, 0.6, rgb(192, 0, 0, max = 255, alpha = 50))
  # COMPLEXCHART = GNOMADSV_CHART(DTB_GNOMADSV_COMPLEX, 0.61, 0.9, rgb(123, 35, 170, max = 255, alpha = 50))

  #BUILD MAPPABILITY
  kpBars(kp, data.panel=2, chr=MAPPABILITY$CONTIG, x0=MAPPABILITY$START, x1=MAPPABILITY$END, y1=MAPPABILITY$VALUE, r0=0.44, r1=0.35, col="black")

  if(file.exists(GET_CTRL_CN)==T) {
  #Add ctrl CN data itself
  kpRect(kp, data.panel=2, chr=CHR, x0=FROM, x1=TO, y0=0, y1=1, col=NA, r0=0.12, r1=0.02, border=NA)
  kpRect(kp, data.panel=2, chr=CHR, x0=FROM, x1=TO, y0=0, y1=1, col=NA, r0=0.13, r1=0.23, border=NA)
  kpRect(kp, data.panel=2, chr=CTRL_CN$CONTIG, x0=CTRL_CN$STARText, x1=CTRL_CN$ENDext, y0=0, y1=CTRL_CN$FofCASES, r0=0.12, r1=0.02, col=rgb(255, 235, 120, max=255, alpha=200), border=NA)
  kpRect(kp, data.panel=2, chr=CTRL_CN$CONTIG, x0=CTRL_CN$STARText, x1=CTRL_CN$ENDext, y0=0, y1=CTRL_CN$FofCASES, r0=0.13, r1=0.23, col=rgb(255, 235, 120, max=255, alpha=200), border=NA)
  kpRect(kp, data.panel=2, chr=CTRL_CN$CONTIG, x0=CTRL_CN$START, x1=CTRL_CN$END, y0=0, y1=CTRL_CN$FofGAIN, r0=0.12, r1=0.02, col=rgb(84, 130, 53, max = 255, alpha = 255), border=NA)
  kpRect(kp, data.panel=2, chr=CTRL_CN$CONTIG, x0=CTRL_CN$START, x1=CTRL_CN$END, y0=0, y1=CTRL_CN$FofLOSS, r0=0.13, r1=0.23, col=rgb(192, 0, 0, max = 255, alpha = 255), border=NA)
  }


  #GENES_DTB_PROCESSED
  #transcription
  if(nrow(GENES_DTB_PROCESSED_REGION)>0) {kpArrows(kp, data.panel=2, chr=GENES_DTB_PROCESSED_REGION$CONTIG, x0=GENES_DTB_PROCESSED_REGION$TxSTART, x1=GENES_DTB_PROCESSED_REGION$TxEND, y0=GENES_DTB_PROCESSED_REGION$CENTER, y1=GENES_DTB_PROCESSED_REGION$CENTER, r0=0.5, r1=1, length=0, lwd=0.5*GENES_DTB_PROCESSED_LayerHeight)}
  if(nrow(GENES_DTB_PROCESSED_ARROWS)>0) {kpArrows(kp, data.panel=2, chr=GENES_DTB_PROCESSED_ARROWS$CONTIG, x0=GENES_DTB_PROCESSED_ARROWS$ARROW_START, x1=GENES_DTB_PROCESSED_ARROWS$ARROW_END, y0=GENES_DTB_PROCESSED_ARROWS$CENTER, y1=GENES_DTB_PROCESSED_ARROWS$CENTER, r0=0.5, r1=1, lwd=0.5*GENES_DTB_PROCESSED_LayerHeight, length=0.13*GENES_DTB_PROCESSED_LayerHeight, col="black")}
  #exons
  if(nrow(EXONS)>0) {kpRect(kp, data.panel=2, chr=EXONS$CONTIG, x0=EXONS$EXON_START, x1=EXONS$EXON_END, y0=EXONS$EXON_Y0, y1=EXONS$EXON_Y1, col="black", r0=0.5, r1=1, border=NA)}
  #coding
  if(nrow(EXONS_CODING)>0) {kpRect(kp, data.panel=2, chr=EXONS_CODING$CONTIG, x0=EXONS_CODING$EXON_CODING_START, x1=EXONS_CODING$EXON_CODING_END, y0=EXONS_CODING$EXON_CODING_Y0, y1=EXONS_CODING$EXON_CODING_Y1, col="black", r0=0.5, r1=1, border=NA)}
  #gene symbols
  if(nrow(GENES_DTB_PROCESSED_REGION)>0) {kpText(kp, data.panel=2, chr=GENES_DTB_PROCESSED_REGION$CONTIG, x=GENES_DTB_PROCESSED_REGION$TxCenter, y=GENES_DTB_PROCESSED_REGION$EXON_ID_Y, labels=GENES_DTB_PROCESSED_REGION$GENE_ID, pos=3, font=3, r0=0.5, r1=1, cex=GENES_DTB_PROCESSED_LayerHeight, clipping = FALSE)}


  invisible(dev.off())
}
}

options(warn = defaultW)

#-----------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------
#project specific continue
defaultW <- getOption("warn")
options(warn = -1)
if(file.exists(PLOT_PROJECT_SOURCE)==T) {
  message(paste0("   R ... ", date(), " - STAGE 5+/5 - project specific plot"))
  source(PLOT_PROJECT_SOURCE)
}
options(warn = defaultW)

message(paste0("   R ... ", date(), " - FINISHED"))
  
  
  
  

