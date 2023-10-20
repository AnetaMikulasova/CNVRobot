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
CHROMSIZE           = args[2] %>% as.character()
CHROMBAND           = args[3] %>% as.character()
CENTROMERE          = args[4] %>% as.character()
DGV                 = args[5] %>% as.character()
GET_MAPPABILITY_DTB = args[6] %>% as.character()
GET_ENCODE_BLACKLIST = args[7] %>% as.character()
GENES_DTB_PROCESSED = args[8] %>% as.character()
PROJECT_ID          = args[9] %>% as.character()
INDIR               = args[10] %>% as.character()
GNOMAD_PATTERN      = args[11] %>% as.character()
GROUP               = args[12] %>% as.character()
SAMPLE1_ID          = args[13] %>% as.character()
SAMPLE1_TYPE        = args[14] %>% as.character()
SAMPLE1_SEX         = args[15] %>% as.character()
SAMPLE2_ID          = args[16] %>% as.character()
SAMPLE2_TYPE        = args[17] %>% as.character()
SAMPLE2_SEX         = args[18] %>% as.character()
SAMPLE3_ID          = args[19] %>% as.character()
SAMPLE3_TYPE        = args[20] %>% as.character()
SAMPLE3_SEX         = args[21] %>% as.character()
CTRL_N              = args[22] %>% as.character()
PON_STYLE           = args[23] %>% as.character()
OUTDIR              = args[24] %>% as.character()
# OUTDIR_IGV          = args[25] %>% as.character()
GETSEGMENTCONDITION = args[26] %>% as.character()
MINLOSS  = args[27] %>% as.numeric()
MINGAIN  = args[28] %>% as.numeric()
GET_CTRL_CN              = args[29] %>% as.character()
# PLOT_PROJECT_SOURCE = args[30] %>% as.character()
# DETAIL_PLOT_SKIP = args[31] %>% as.character()
REGIONS_TABLE       = args[32] %>% as.character()
METHOD = args[33] %>% as.character()
  
message(paste0("   R ... ", date(), " - STAGE 1/3 - loading and processing databases"))

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
DTB_CHRSIZE = DTB_CHRSIZE
DTB_CHRSIZE$SIZE= as.numeric(as.integer(DTB_CHRSIZE$SIZE))
#order as chr1, chr2, chr3...Y or 1, 2, 3...Y (if not genome view can get strange orders)
DTB_CHRSIZE= DTB_CHRSIZE %>%
  mutate(contig_temp = gsub("chr", "", DTB_CHRSIZE$CONTIG)) %>%
  mutate(contig_temp = ifelse(contig_temp == "X", 23, contig_temp)) %>%
  mutate(contig_temp = ifelse(contig_temp == "Y", 24, contig_temp)) %>%
  mutate(contig_temp = as.integer(contig_temp)) %>%
  arrange(contig_temp) %>%
  select(CONTIG, SIZE)
CONTIGS_in_SAMPLE=DTB_CHRSIZE$CONTIG %>% unique()

#Get chromosome cytobands
DTB_CHROMBAND = read.delim(CHROMBAND, header = FALSE, stringsAsFactors = F) %>% mutate(V2 = V2 + 1)
names(DTB_CHROMBAND) = c("chr", "start", "end", "name", "gieStain")
DTB_CHROMBAND = DTB_CHROMBAND

#Use chromosome sizes for defining the custom genome
CUSTOM_GENOME = toGRanges(data.frame(chr=DTB_CHRSIZE$CONTIG, start=1, end=DTB_CHRSIZE$SIZE))

#Get capture
if(file.exists(GETCAPTURE)) {
  DTB_CAPTURE = read.delim(GETCAPTURE, header = FALSE, stringsAsFactors = F) %>% select(V1, V2, V3)
  names(DTB_CAPTURE) = c("CONTIG",  "START", "END")
  DTB_CAPTURE = DTB_CAPTURE %>%
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
DTB_CENTROMERE = DTB_CENTROMERE
DTB_CENTROMERE_temp = DTB_CENTROMERE %>% mutate(POSITION2 = POSITION) 
names(DTB_CENTROMERE_temp) = c("seqnames", "start", "end")
DTB_CENTROMERE = cbind(DTB_CENTROMERE_temp, DTB_CENTROMERE) %>% toGRanges()
remove(DTB_CENTROMERE_temp)

#Get DGV
DTB_DGVDATA = readRDS(DGV)
#remove "chr" if not present in sample file
if (("TRUE" %in% grepl("chr", CONTIGS_in_SAMPLE))==FALSE) {
  DTB_DGVDATA$CONTIG = gsub("chr", "", DTB_DGVDATA$CONTIG)
}
DTB_DGVDATA_temp = DTB_DGVDATA %>% select(CONTIG, START, END) 
names(DTB_DGVDATA_temp) = c("seqnames", "start", "end")
DTB_DGVDATA = cbind(DTB_DGVDATA_temp, DTB_DGVDATA)
DTB_DGV_LOSS = filter(DTB_DGVDATA, TYPE == "LOSS") %>% toGRanges()
DTB_DGV_GAIN = filter(DTB_DGVDATA, TYPE == "GAIN") %>% toGRanges()
DTB_DGV_COMPLEX = filter(DTB_DGVDATA, TYPE == "COMPLEX") %>% toGRanges()
remove(DTB_DGVDATA, DTB_DGVDATA_temp)

#load mappability dtb
MAPPABILITY = readRDS(GET_MAPPABILITY_DTB)

#load ENCODE blacklist
ENCODE_BLACKLIST = read.delim(paste0(GET_ENCODE_BLACKLIST), header=TRUE, stringsAsFactors = F)

#Get processed gene database
DTB_GENES_DTB_PROCESSED  = read.delim(GENES_DTB_PROCESSED, header = TRUE, stringsAsFactors = F)

#remove "chr" if not present in sample file
if (("TRUE" %in% grepl("chr", CONTIGS_in_SAMPLE))==FALSE) {
  MAPPABILITY$CONTIG = gsub("chr", "", MAPPABILITY$CONTIG)
  ENCODE_BLACKLIST$CONTIG = gsub("chr", "", ENCODE_BLACKLIST$CONTIG)
  DTB_GENES_DTB_PROCESSED$CONTIG = gsub("chr", "", DTB_GENES_DTB_PROCESSED$CONTIG)
}

#define odd contigs (for coloring contigs in genome plot)
CONTIG_ODD=c(seq(1, 22, 2), "X")
CONTIG_ODD = c(CONTIG_ODD, as.vector(outer("chr", CONTIG_ODD, paste, sep="")))

#if(is.na(DETAIL_PLOT_SKIP)) {DETAIL_PLOT_SKIP="none"}

message(paste0("   R ... ", date(), " - STAGE 2/3 - loading samples data"))

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
  
  
  #define files patterns based on method:
  if(METHOD == "gatk") {
    MODEL_FILE_PATTERN = "_segmentation.tsv"
    DENOIS_FILE_PATTERN = ".tsv"
    ALL_FILE_PATTERN = ".rds"
  }
  
  if(METHOD == "cnvkit") {
    MODEL_FILE_PATTERN = "_segmentation_CNVkit.tsv"
    DENOIS_FILE_PATTERN = ".cnr"
    ALL_FILE_PATTERN = ".rds"
  }
  
  if(file.exists(paste0(INDIR, i, "_SEGMENTS_C-", CTRL_N, "_", PON, "_", GETSEGMENTCONDITION, MODEL_FILE_PATTERN))) {

    MODEL <- read.delim(paste0(INDIR, i, "_SEGMENTS_C-", CTRL_N, "_", PON, "_", GETSEGMENTCONDITION, MODEL_FILE_PATTERN), stringsAsFactors = F) %>% 
      # filter(FILTER == "PASS") %>%
      filter(TYPE == "NORMAL" | TYPE == "LOSS" | TYPE == "GAIN" | TYPE == "subLOSS" | TYPE == "subGAIN") %>%
      mutate(COLOR = "#000000") %>%
      mutate(COLOR = ifelse(TYPE == "LOSS", "#C02600", COLOR)) %>%
      mutate(COLOR = ifelse(TYPE == "GAIN", "#006200", COLOR))
 
    MODEL_LOH <- read.delim(paste0(INDIR, i, "_SEGMENTS_C-", CTRL_N, "_", PON, "_", GETSEGMENTCONDITION, MODEL_FILE_PATTERN), stringsAsFactors = F) %>% 
      # filter(FILTER == "PASS") %>%
      filter(TYPE == "HOM" | TYPE == "HET" | TYPE == "OTHER") %>%
      mutate(COLOR = ifelse(TYPE == "HOM", "#B0C0E2", "#DBE3F1"))

    # #bed file for IGV with abnormal segments
    # IGV_MODEL <- read.delim(paste0(INDIR, i, "_SEGMENTS_C-", CTRL_N, "_", PON, "_", GETSEGMENTCONDITION, MODEL_FILE_PATTERN), stringsAsFactors = F) %>% 
    #   # filter(FILTER == "PASS") %>%
    #   select(CONTIG, START, END, TYPE, MEAN_L2R_EDIT) %>%
    #   mutate(START = START-1) %>%
    #   mutate(x1 = ".") %>%
    #   mutate(x2 = ".") %>%
    #   mutate(START2 = START) %>%
    #   mutate(END2 = END) %>%
    #   mutate(COLOR = ifelse(TYPE == "LOSS", "192,38,0", ".")) %>%
    #   mutate(COLOR = ifelse(TYPE == "GAIN", "0,98,0", COLOR)) %>% 
    #   mutate(COLOR = ifelse(TYPE == "HOM", "142,169,219", COLOR)) %>%
    #   mutate(COLOR = ifelse(TYPE == "OTHER", "142,169,219", COLOR)) %>%
    #   mutate(COLOR = ifelse(TYPE == "subLOSS", "255,172,168", COLOR)) %>%
    #   mutate(COLOR = ifelse(TYPE == "subGAIN", "168,211,121", COLOR)) %>%
    #   arrange(CONTIG, START, END)
    # 
    # IGV_MODEL_BED = IGV_MODEL %>% 
    #   filter(TYPE != "NORMAL") %>%
    #   select(CONTIG, START, END, TYPE, x1, x2, START2, END2, COLOR)
    # write_tsv(IGV_MODEL, paste0(OUTDIR_IGV, GROUP, "_", i, "_abnormal_segments.bed"), col_names = F)
    # 
    # IGV_MODEL_WIG = IGV_MODEL %>% 
    #   filter(TYPE != "HOM") %>%
    #   filter(TYPE != "OTHER") %>%
    #   select(CONTIG, START, END, MEAN_L2R_EDIT)
    # write_tsv(IGV_MODEL_WIG, paste0(OUTDIR_IGV, GROUP, "_", i, "_segments"), col_names = F)
    
    DENOIS  <-  read.delim(paste0(INDIR, i, "_denoisedCR_C-", CTRL_N, "_", PON, DENOIS_FILE_PATTERN), comment.char = "@", stringsAsFactors = F) %>%
      mutate(COLOR = "#A6A6A6") %>%
      # mutate(COLOR = ifelse(LOG2_COPY_RATIO < -0.9, "#FF7100", COLOR)) %>%
      # mutate(COLOR = ifelse(LOG2_COPY_RATIO < -1.5, "#FF2600", COLOR)) %>%
      # mutate(COLOR = ifelse(LOG2_COPY_RATIO > 0.65, "#008F00", COLOR)) %>%
      # mutate(COLOR = ifelse(LOG2_COPY_RATIO > 1.00, "#47AFEB", COLOR)) %>%
      mutate(COLOLCONTIG = ".") %>%
      mutate(COLOLCONTIG = ifelse(CONTIG %in% CONTIG_ODD, "#FEB49B", "#9BC2E6")) %>%
      mutate(L2R = ifelse(LOG2_COPY_RATIO < -2.25, -2.25, LOG2_COPY_RATIO)) %>%
      mutate(L2R = ifelse(LOG2_COPY_RATIO > 2.25, 2.25, L2R)) %>%
      mutate(COLORABN = COLOR)
    
          #re-color denois based on segments from segmentations
          MODELtemp = MODEL %>% select(CONTIG, START, END, TYPE)
          write_tsv(MODELtemp, paste0(INDIR, "temp1"), col_names = F)
          MODELtemp = paste0(INDIR, "temp1")
          
          DENOIStemp = DENOIS %>% select(CONTIG, START, END)
          write_tsv(DENOIStemp, paste0(INDIR, "temp2"), col_names = F)
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
    # IGV_DENOIS  <-  read.delim(paste0(INDIR, i, "_denoisedCR_C-", CTRL_N, "_", PON, DENOIS_FILE_PATTERN), comment.char = "@", stringsAsFactors = F) %>%
    #   mutate(L2R = ifelse(LOG2_COPY_RATIO < -2.25, -2.25, LOG2_COPY_RATIO)) %>%
    #   mutate(L2R = ifelse(LOG2_COPY_RATIO > 2.25, 2.25, L2R)) %>%
    #   select(CONTIG, START, END, L2R) %>%
    #   mutate(START = START-1) %>%
    #   arrange(CONTIG, START, END)
    # write_tsv(IGV_DENOIS, paste0(OUTDIR_IGV, GROUP, "_", i, "_CN"), col_names = F)
    
    ALL    <-   readRDS(paste0(INDIR, i, "_allelicCounts", GNOMAD_PATTERN, ALL_FILE_PATTERN)) %>%
      mutate(N = REF_COUNT + ALT_COUNT) %>%
      filter(N > 2)
    ALL = ALL %>%
      mutate(TEST = sample(1:2, nrow(ALL), replace = T)) %>%
      mutate(MAF = ifelse(TEST == 1, ALT_COUNT/N, REF_COUNT/N)) %>%
      mutate(COL = "#000000") %>%
      mutate(COL = ifelse(MAF > 0.375 & MAF < 0.625, "#008F00", COL)) %>%
      mutate(COL = ifelse(MAF > 0.125 & MAF < 0.375, "#FF2600", COL)) %>%
      mutate(COL = ifelse(MAF > 0.625 & MAF < 0.875, "#FF2600", COL)) %>%
      mutate(COLOLCONTIG = ".") %>%
      mutate(COLOLCONTIG = ifelse(CONTIG %in% CONTIG_ODD, "#FEB49B", "#9BC2E6"))

    ALL_FILTERED = ALL %>% sample_frac(0.5) %>% arrange(CONTIG, POSITION)
    
    #bedGraph file with denoised copy-number data
    # ALL_IGV <- readRDS(paste0(INDIR, i, "_allelicCounts", GNOMAD_PATTERN, ALL_FILE_PATTERN)) %>%
    #   mutate(N = REF_COUNT + ALT_COUNT) %>%
    #   filter(N > 2)
    # ALL_IGV = ALL_IGV %>%
    #   mutate(TEST = sample(1:2, nrow(ALL_IGV), replace = T)) %>%
    #   mutate(MAF = ifelse(TEST == 1, ALT_COUNT/N, REF_COUNT/N)) %>%
    #   mutate(POSITION = POSITION-1) %>%
    #   mutate(POSITION2 = POSITION) %>%
    #   select(CONTIG, POSITION, POSITION2, MAF) %>%
    #   arrange(CONTIG, POSITION)
    # write_tsv(ALL_IGV, paste0(OUTDIR_IGV, GROUP, "_", i, "_SNP"), col_names = F)
    
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
if(file.exists(paste0(INDIR, SAMPLE1_TYPE, "_", SAMPLE1_ID, "_", SAMPLE1_SEX, "_SEGMENTS_C-", CTRL_N, "_", PON1, "_", GETSEGMENTCONDITION, MODEL_FILE_PATTERN)) & 
   file.exists(paste0(INDIR, SAMPLE2_TYPE, "_", SAMPLE2_ID, "_", SAMPLE2_SEX, "_SEGMENTS_C-", CTRL_N, "_", PON2, "_", GETSEGMENTCONDITION, MODEL_FILE_PATTERN)) &
   file.exists(paste0(INDIR, SAMPLE3_TYPE, "_", SAMPLE3_ID, "_", SAMPLE3_SEX, "_SEGMENTS_C-", CTRL_N, "_", PON3, "_", GETSEGMENTCONDITION, MODEL_FILE_PATTERN))) {PLOT_STYLE = "RUN_TRIO"}

if(file.exists(paste0(INDIR, SAMPLE1_TYPE, "_", SAMPLE1_ID, "_", SAMPLE1_SEX, "_SEGMENTS_C-", CTRL_N, "_", PON1, "_", GETSEGMENTCONDITION, MODEL_FILE_PATTERN)) & 
   file.exists(paste0(INDIR, SAMPLE2_TYPE, "_", SAMPLE2_ID, "_", SAMPLE2_SEX, "_SEGMENTS_C-", CTRL_N, "_", PON2, "_", GETSEGMENTCONDITION, MODEL_FILE_PATTERN)) &
   file.exists(paste0(INDIR, SAMPLE3_TYPE, "_", SAMPLE3_ID, "_", SAMPLE3_SEX, "_SEGMENTS_C-", CTRL_N, "_", PON3, "_", GETSEGMENTCONDITION, MODEL_FILE_PATTERN))==F) {PLOT_STYLE = "RUN_DUO"}

if(file.exists(paste0(INDIR, SAMPLE1_TYPE, "_", SAMPLE1_ID, "_", SAMPLE1_SEX, "_SEGMENTS_C-", CTRL_N, "_", PON1, "_", GETSEGMENTCONDITION, MODEL_FILE_PATTERN)) & 
   file.exists(paste0(INDIR, SAMPLE2_TYPE, "_", SAMPLE2_ID, "_", SAMPLE2_SEX, "_SEGMENTS_C-", CTRL_N, "_", PON2, "_", GETSEGMENTCONDITION, MODEL_FILE_PATTERN))==F &
   file.exists(paste0(INDIR, SAMPLE3_TYPE, "_", SAMPLE3_ID, "_", SAMPLE3_SEX, "_SEGMENTS_C-", CTRL_N, "_", PON3, "_", GETSEGMENTCONDITION, MODEL_FILE_PATTERN))) {
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

if(file.exists(paste0(INDIR, SAMPLE1_TYPE, "_", SAMPLE1_ID, "_", SAMPLE1_SEX, "_SEGMENTS_C-", CTRL_N, "_", PON1, "_", GETSEGMENTCONDITION, MODEL_FILE_PATTERN)) & 
   file.exists(paste0(INDIR, SAMPLE2_TYPE, "_", SAMPLE2_ID, "_", SAMPLE2_SEX, "_SEGMENTS_C-", CTRL_N, "_", PON2, "_", GETSEGMENTCONDITION, MODEL_FILE_PATTERN))==F &
   file.exists(paste0(INDIR, SAMPLE3_TYPE, "_", SAMPLE3_ID, "_", SAMPLE3_SEX, "_SEGMENTS_C-", CTRL_N, "_", PON3, "_", GETSEGMENTCONDITION, MODEL_FILE_PATTERN))==F) {PLOT_STYLE = "RUN_SINGLE"}
# 
# 









message(paste0("   R ... ", date(), " - STAGE 3/3 - region plot"))

defaultW <- getOption("warn") 
options(warn = -1)

#------DETAIL-VIEW-----------------------------------------------------------------------------------------------------------------------------------
DETAILS = read.delim(REGIONS_TABLE, header = TRUE, stringsAsFactors = F) 
DETAILS = mutate(DETAILS, N = 1:n())

LISTofABNORM = unique(DETAILS$N) %>% as.list()
NofABN = max(unlist(LISTofABNORM))

for (ABN in (LISTofABNORM)) {

  message(paste0("   R ... ", date(), " - STAGE 3/3 - region plot: ", ABN, "/", NofABN))

  DETAILS_PART = DETAILS %>%
    filter(N == ABN)

  FROMABN = DETAILS_PART$START
  TOABN = DETAILS_PART$END
  SIZEABN = TOABN-FROMABN
  ZOOMCOEF = DETAILS_PART$ZOOM

  CHR = DETAILS_PART$CONTIG %>% as.character()
  FROM = DETAILS_PART$START - SIZEABN*ZOOMCOEF
  if (FROM < 0) {FROM = 0}
  TO = DETAILS_PART$END + SIZEABN*ZOOMCOEF
  SIZE = TO-FROM
  
  output = paste0(OUTDIR, "detail_regions_of_interest/", GROUP, "_", SAMPLE1_TYPE, "_", SAMPLE1_ID, "_", SAMPLE1_SEX, "_", CHR, "-", FROMABN, "-", TOABN, "_size-", SIZEABN, "bp")

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
  
  CUSTOM_CYTOBAND = DTB_CHROMBAND %>% filter(chr == CHR)
  CUSTOM_CYTOBAND  = toGRanges(CUSTOM_CYTOBAND)
  
  kp <- plotKaryotype(genome=CUSTOM_GENOME, cytobands = CUSTOM_CYTOBAND,
                      chromosomes = CHR,
                      plot.type = 2,
                      plot.params = plot.params,
                      # ideogram.plotter = NULL,
                      labels.plotter = NULL,
                      main=paste0(GROUP, "  ●  ", CHR, ":", FROMABN, "-", TOABN, " ● ", SIZEABN, "bp"),
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

  #BUILD MAPPABILITY
  kpBars(kp, data.panel=2, chr=MAPPABILITY$CONTIG, x0=MAPPABILITY$START, x1=MAPPABILITY$END, y1=MAPPABILITY$VALUE, r0=0.44, r1=0.35, col="black")
  
  #BUILD ENCODE BLACKLIST
  kpRect(kp, data.panel=2, chr=ENCODE_BLACKLIST$CONTIG, x0=ENCODE_BLACKLIST$START, x1=ENCODE_BLACKLIST$END, y0=0, y1=1, r0=0.48, r1=0.46, col=rgb(130, 130, 130, max=255), border=NA)

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

options(warn = defaultW)

#-----------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------

message(paste0("   R ... ", date(), " - FINISHED"))
  
system(paste0("touch status_ok"))
  
  

