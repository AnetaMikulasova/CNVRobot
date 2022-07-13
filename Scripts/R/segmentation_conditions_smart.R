args <- commandArgs(trailing = TRUE)

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(splitstackshape))


SEGMENTATION=args[1] %>% as.character()
GETDENOIS= args[2] %>% as.character()
GETINTERVALS= args[3] %>% as.character()
SEQ_TYPE = args[4] %>% as.character()

options("scipen"=999)

#Functions
rowShift <- function(x, shiftLen = 1L) {
  r <- (1L + shiftLen):(length(x) + shiftLen)
  r[r<1] <- NA
  return(x[r])
}


# # if(DETECTION == "smart-sensitive") {COEFICIENT = 0.75}
# # if(DETECTION == "smart-specific") {COEFICIENT = 0.90}
# 
# COEFICIENT=gsub("smart-","", SEGMENTATION) 
# COEFICIENT=gsub("-sub","", COEFICIENT) %>% as.numeric()

LIST=unlist(strsplit(SEGMENTATION, "-|_"))
COEFICIENT=LIST[2] %>% as.numeric()
COEFICIENT_SUB=LIST[4] %>% as.numeric()


if(GETDENOIS != "na") {
#calculate smart qc
DENOIS_DATA_INPUT   <- read.delim(paste0(GETDENOIS), comment.char = "@", stringsAsFactors = F)
DENOIS = DENOIS_DATA_INPUT %>% mutate(DIF = abs(LOG2_COPY_RATIO-rowShift(LOG2_COPY_RATIO, -1)))
DENOIS_AUTOSOME = filter(DENOIS, CONTIG != "chrX" & CONTIG != "chrY" & CONTIG != "X" & CONTIG != "Y")
SD_DENOIS_AUTOSOME=sd(DENOIS_AUTOSOME$DIF, na.rm=T)
}
if(GETDENOIS == "na") {
SD_DENOIS_AUTOSOME=0.5
}

#calculate smart bin
INTERVALS=read.delim(GETINTERVALS, comment.char = "@", stringsAsFactors = F) %>%
  mutate(SIZE = END-START+1)
BIN_MEDIAN=median(INTERVALS$SIZE)

#calculate smart segm conditions
SEGM_DIFFERENCE = COEFICIENT*(2*SD_DENOIS_AUTOSOME^2+0.6*SD_DENOIS_AUTOSOME+0.25)
if(SEGM_DIFFERENCE > 1) {SEGM_DIFFERENCE = 1}
SEGM_MINSIZE = 2*BIN_MEDIAN
if(grepl("-sub", SEGMENTATION)==T) {SEGM_MINSIZE_SUB=5*BIN_MEDIAN} else {SEGM_MINSIZE_SUB = SEGM_MINSIZE}
SEGM_MINPROBE = 2
if(grepl("-sub", SEGMENTATION)==T) {SEGM_MINPROBE_SUB=5} else {SEGM_MINPROBE_SUB = SEGM_MINPROBE}
SEGM_MINKEEP = 2
SEGM_MINLOSS = log(1-0.5*COEFICIENT,2)
if(grepl("-sub", SEGMENTATION)==T) {SEGM_MINLOSS_SUB = log(1-0.5*(COEFICIENT*COEFICIENT_SUB),2)} else {SEGM_MINLOSS_SUB = SEGM_MINLOSS}
SEGM_MINLOSS_BIAL = -1.5
SEGM_MINGAIN = log(1+0.5*COEFICIENT,2)
if(grepl("-sub", SEGMENTATION)==T) {SEGM_MINGAIN_SUB = log(1+0.5*(COEFICIENT*COEFICIENT_SUB),2)} else {SEGM_MINGAIN_SUB = SEGM_MINGAIN}
SEGM_GAP = 500000
SEGM_SMOOTHPERC = 0.9
SEGM_AFDIF = 0.1
if(SEQ_TYPE == "WGS") {SEGM_AFSIZE = 5000000/2} else {SEGM_AFSIZE = 5000000}
if(SEQ_TYPE == "WGS") {SEGM_AFPROBE = 10000/2} else {SEGM_AFPROBE = 10000}

# 
# #round numbers
SEGM_DIFFERENCE = format(round(SEGM_DIFFERENCE, 2), nsmall = 2)
SEGM_MINSIZE = round(SEGM_MINSIZE, 0)
SEGM_MINSIZE_SUB = round(SEGM_MINSIZE_SUB, 0)
SEGM_MINPROBE = round(SEGM_MINPROBE, 0)
SEGM_MINPROBE_SUB = round(SEGM_MINPROBE_SUB, 0)
SEGM_MINKEEP = round(SEGM_MINKEEP, 0)
SEGM_MINLOSS = format(round(SEGM_MINLOSS, 2), nsmall = 2)
SEGM_MINLOSS_SUB = format(round(SEGM_MINLOSS_SUB, 2), nsmall = 2)
SEGM_MINLOSS_BIAL = format(round(SEGM_MINLOSS_BIAL, 2), nsmall = 2)
SEGM_MINGAIN = format(round(SEGM_MINGAIN, 2), nsmall = 2)
SEGM_MINGAIN_SUB = format(round(SEGM_MINGAIN_SUB, 2), nsmall = 2)
SEGM_GAP = round(SEGM_GAP, 0)
SEGM_SMOOTHPERC = format(round(SEGM_SMOOTHPERC, 2), nsmall = 2)
SEGM_AFDIF = format(round(SEGM_AFDIF, 2), nsmall = 2)
SEGM_AFSIZE = round(SEGM_AFSIZE, 0)
SEGM_AFPROBE = round(SEGM_AFPROBE, 0)

# # paste0(SD_DENOIS_AUTOSOME, "_", BIN_MEDIAN)
# c(SEGM_DIFFERENCE, 
#        SEGM_MINSIZE, SEGM_MINSIZE_SUB, SEGM_MINPROBE, SEGM_MINPROBE_SUB, SEGM_MINKEEP,
#        SEGM_MINLOSS, SEGM_MINLOSS_SUB, SEGM_MINLOSS_BIAL, 
#        SEGM_MINGAIN, SEGM_MINGAIN_SUB,
#        SEGM_GAP, SEGM_SMOOTHPERC,
#        SEGM_AFDIF, SEGM_AFSIZE, SEGM_AFPROBE
#        )

paste0(SEGM_DIFFERENCE, "_",
  SEGM_MINSIZE, "_", SEGM_MINSIZE_SUB, "_", SEGM_MINPROBE, "_", SEGM_MINPROBE_SUB, "_", SEGM_MINKEEP, "_",
  SEGM_MINLOSS, "_", SEGM_MINLOSS_SUB, "_", SEGM_MINLOSS_BIAL, "_", 
  SEGM_MINGAIN, "_", SEGM_MINGAIN_SUB, "_",
  SEGM_GAP, "_", SEGM_SMOOTHPERC, "_",
  SEGM_AFDIF, "_", SEGM_AFSIZE, "_", SEGM_AFPROBE
)
