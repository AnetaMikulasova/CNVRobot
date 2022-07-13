args <- commandArgs(trailing = TRUE)

message(paste0("   R ... ", date(), " - STARTING"))

suppressPackageStartupMessages(library(tidyverse))

GET_PROJECT_TYPE     = args[1]%>% as.character()
GET_SAMPLE1_ID       = args[2]%>% as.character()
GET_SAMPLE1_SEX      = args[3]%>% as.character()
GET_SAMPLE1_PON_SEX  = args[4]%>% as.character()
GET_SEGMENTATION     = args[5] %>% as.character()
GET_SAMPLE1_DENOIS   = args[6] %>% as.character()
GET_SAMPLE2_ID       = args[7]%>% as.character()
GET_SAMPLE2_SEX      = args[8]%>% as.character()
GET_SAMPLE2_PON_SEX  = args[9]%>% as.character()
GET_SAMPLE2_DENOIS   = args[10] %>% as.character()
GET_SAMPLE3_ID       = args[11] %>% as.character()
GET_SAMPLE3_SEX      = args[12]%>% as.character()
GET_SAMPLE3_PON_SEX  = args[13]%>% as.character()
GET_SAMPLE3_DENOIS   = args[14] %>% as.character()
MINGAIN              = args[15] %>% as.numeric() #* 0.95
MINGAIN_SUB          = args[16] %>% as.numeric()
MINLOSS              = args[17] %>% as.numeric() #* 0.95
MINLOSS_SUB          = args[18] %>% as.numeric()
OUTPUT               = args[19]%>% as.character()

message(paste0("   R ... ", date(), " - STAGE 1/4 - loading data"))


#project type
PROJECT_TYPE = GET_PROJECT_TYPE

# #rename SAMPLE2 and SAMPLE3 if they don't exist
# if(file.exists(GET_SAMPLE2_DENOIS)==F) {GET_SAMPLE2_ID="x"}
# if(file.exists(GET_SAMPLE3_DENOIS)==F) {GET_SAMPLE3_ID="x"}

#get segmentation table of proband
SEGMENTATION_INPUT = read.delim(GET_SEGMENTATION, stringsAsFactors = F) %>%
  filter(TYPE != "NORMAL") %>%
  # filter(FILTER == "PASS") %>%
  mutate(SAMPLE2_TYPE_ID_SEX = paste0(GET_SAMPLE2_ID, "_", GET_SAMPLE2_SEX)) %>%
  mutate(SAMPLE3_TYPE_ID_SEX = paste0(GET_SAMPLE3_ID, "_", GET_SAMPLE3_SEX)) 
  # select(SEGMENTATION_CONDITION, PON_SEX_N, SAMPLE_ID1, SAMPLE_ID2, SAMPLE_TYPE, SAMPLE_SEX, SAMPLE2_TYPE_ID_SEX, SAMPLE3_TYPE_ID_SEX, CONTIG, START, END, MEAN_L2R, SIZE, PROBES, FILTER, TYPE, CHRBAND, GENE, CTRL_N, CTRL_TYPE, CTRL_FofGAIN, CTRL_FofLOSS, CTRL_PREDICTION, DGV_NofGAIN_AROUND, DGV_NofLOSS_AROUND, DGV_PREDICTION)


# #select what part from SAMPLE2 and SAMPLE3 will go to the inheritance
# if(file.exists(GET_SAMPLE2_DENOIS)) {
#   TYPE = sapply(strsplit(GET_SAMPLE2_ID, "_"), "[", 1)
#   ID = sapply(strsplit(GET_SAMPLE2_ID, "_"), "[", 2)
#   GET_SAMPLE2_ID = paste0(TYPE, "_", ID)
#   }
# if(file.exists(GET_SAMPLE3_DENOIS)) {
#   TYPE = sapply(strsplit(GET_SAMPLE3_ID, "_"), "[", 1)
#   ID = sapply(strsplit(GET_SAMPLE3_ID, "_"), "[", 2)
#   GET_SAMPLE3_ID = paste0(TYPE, "_", ID)
# }


SEGMENTATION_INPUT$SAMPLE_SEX = as.character(SEGMENTATION_INPUT$SAMPLE_SEX)
SEGMENTATION_INPUT = SEGMENTATION_INPUT %>%
  mutate(SAMPLE_SEX = ifelse(SAMPLE_SEX == "FALSE", "F", SAMPLE_SEX))

SEGMENTATION = SEGMENTATION_INPUT %>% filter(TYPE == "LOSS" | TYPE == "GAIN" | TYPE == "subLOSS" | TYPE == "subGAIN")
SEGMENTATION_LOH = SEGMENTATION_INPUT %>% filter(TYPE == "HOM" | TYPE == "OTHER")





#recalculate loss and gain cut-offs for gonosomes based on given standard cut-offs
#clonal cut-off
#------------------------------------------------
#percentage of the expected
MINLOSS_PERC = (1-(2^MINLOSS))/0.5
# MINLOSS_BIAL_PERC = (1-(2^MINLOSS_BIAL))/1
MINGAIN_PERC = ((2^MINGAIN)-1)/0.5

#cut-off for male sample denoised by female controls (chrY is lost)
MINLOSS_MvsPoNF_chrX = log(0.5-(MINLOSS_PERC*0.5),2)
if(MINLOSS_MvsPoNF_chrX < -2) {MINLOSS_MvsPoNF_chrX = -2}
MINGAIN_MvsPoNF_chrX = log(0.5+(MINGAIN_PERC*0.5),2)

#cut-off for female sample denoised by male controls
MINLOSS_FvsPoNM_chrX = log(2-(MINLOSS_PERC*1),2)
# MINLOSS_BIAL_FvsPoNM_chrX = log(2-(MINLOSS_BIAL_PERC*2),2)
# if(MINLOSS_BIAL_FvsPoNM_chrX < -1) {MINLOSS_BIAL_FvsPoNM_chrX = -1}
MINGAIN_FvsPoNM_chrX = log(2+(MINGAIN_PERC*1),2)
#MINLOSS_FvsPoNM_chrY - cannot happen as there is no chrY expected in female, all is lost
MINGAIN_FvsPoNM_chrY = log(0+(MINGAIN_PERC*1),2)
#------------------------------------------------

#sub-clonal cut-off
#------------------------------------------------
#percentage of the expected
MINLOSS_SUB_PERC = (1-(2^MINLOSS_SUB))/0.5
MINGAIN_SUB_PERC = ((2^MINGAIN_SUB)-1)/0.5

#cut-off for male sample denoised by female controls (chrY is lost)
MINLOSS_SUB_MvsPoNF_chrX = log(0.5-(MINLOSS_SUB_PERC*0.5),2)
if(MINLOSS_SUB_MvsPoNF_chrX < -2) {MINLOSS_SUB_MvsPoNF_chrX = -2}
MINGAIN_SUB_MvsPoNF_chrX = log(0.5+(MINGAIN_SUB_PERC*0.5),2)

#cut-off for female sample denoised by male controls
MINLOSS_SUB_FvsPoNM_chrX = log(2-(MINLOSS_SUB_PERC*1),2)
MINGAIN_SUB_FvsPoNM_chrX = log(2+(MINGAIN_SUB_PERC*1),2)
#MINLOSS_SUB_FvsPoNM_chrY - cannot happen as there is no chrY expected in female, all is lost
MINGAIN_SUB_FvsPoNM_chrY = log(0+(MINGAIN_SUB_PERC*1),2)
#------------------------------------------------



# #define sex patterns
# if(GET_SAMPLE1_SEX == GET_SAMPLE1_PON_SEX) {SAMPLE1_SEX_vs_PON = "SAME"} #both M, F, mixed or none
# if(GET_SAMPLE1_SEX == "M" & GET_SAMPLE1_PON_SEX == "F") {SAMPLE1_SEX_vs_PON = "M_PoN-F"}
# if(GET_SAMPLE1_SEX == "F" & GET_SAMPLE1_PON_SEX == "M") {SAMPLE1_SEX_vs_PON = "F_PoN-M"}
# if(file.exists(GET_SAMPLE2_DENOIS)) {
#   if(GET_SAMPLE2_SEX == GET_SAMPLE2_PON_SEX) {SAMPLE2_SEX_vs_PON = "SAME"} #both M, F, mixed or none
#   if(GET_SAMPLE2_SEX == "M" & GET_SAMPLE2_PON_SEX == "F") {SAMPLE2_SEX_vs_PON = "M_PoN-F"}
#   if(GET_SAMPLE2_SEX == "F" & GET_SAMPLE2_PON_SEX == "M") {SAMPLE2_SEX_vs_PON = "F_PoN-M"}
# }
# if(file.exists(GET_SAMPLE3_DENOIS)) {
#   if(GET_SAMPLE3_SEX == GET_SAMPLE3_PON_SEX) {SAMPLE3_SEX_vs_PON = "SAME"} #both M, F, mixed or none
#   if(GET_SAMPLE3_SEX == "M" & GET_SAMPLE3_PON_SEX == "F") {SAMPLE3_SEX_vs_PON = "M_PoN-F"}
#   if(GET_SAMPLE3_SEX == "F" & GET_SAMPLE3_PON_SEX == "M") {SAMPLE3_SEX_vs_PON = "F_PoN-M"}
# }

PROCESS_DENOIS_FCE = function(DATA, SAMPLE_SEX, SAMPLE_PON_SEX) {
  DATA = DATA %>%
    mutate(MINLOSS_EDIT = MINLOSS) %>%
    mutate(MINLOSS_SUB_EDIT = MINLOSS_SUB) %>%
    mutate(MINGAIN_EDIT = MINGAIN) %>%
    mutate(MINGAIN_SUB_EDIT = MINGAIN_SUB)
  if (SAMPLE_SEX == "M" & SAMPLE_PON_SEX == "F") {
    DATA = DATA %>%
      mutate(MINLOSS_EDIT = ifelse(CONTIG == "chrX" | CONTIG == "X", MINLOSS_MvsPoNF_chrX, MINLOSS_EDIT)) %>%
      mutate(MINLOSS_SUB_EDIT = ifelse(CONTIG == "chrX" | CONTIG == "X", MINLOSS_SUB_MvsPoNF_chrX,  MINLOSS_SUB_EDIT)) %>%
      mutate(MINGAIN_EDIT = ifelse(CONTIG == "chrX" | CONTIG == "X", MINGAIN_MvsPoNF_chrX,  MINGAIN_EDIT)) %>%
      mutate(MINGAIN_SUB_EDIT = ifelse(CONTIG == "chrX" | CONTIG == "X", MINGAIN_SUB_MvsPoNF_chrX,   MINGAIN_SUB_EDIT))
  }
  if (SAMPLE_SEX == "F" & SAMPLE_PON_SEX == "M") {
    DATA = DATA %>%
      #chrX
      mutate(MINLOSS_EDIT = ifelse(CONTIG == "chrX" | CONTIG == "X", MINLOSS_FvsPoNM_chrX, MINLOSS_EDIT)) %>%
      mutate(MINLOSS_SUB_EDIT = ifelse(CONTIG == "chrX" | CONTIG == "X", MINLOSS_SUB_FvsPoNM_chrX,  MINLOSS_SUB_EDIT)) %>%
      # mutate(MINLOSS_BIAL_EDIT = ifelse(CONTIG == "chrX" | CONTIG == "X", MINLOSS_BIAL_FvsPoNM_chrX, MINLOSS_BIAL_EDIT)) %>%
      mutate(MINGAIN_EDIT = ifelse(CONTIG == "chrX" | CONTIG == "X", MINGAIN_FvsPoNM_chrX, MINGAIN_EDIT)) %>%
      mutate(MINGAIN_SUB_EDIT = ifelse(CONTIG == "chrX" | CONTIG == "X", MINGAIN_SUB_FvsPoNM_chrX, MINGAIN_SUB_EDIT)) %>%
      #chrY
      mutate(MINLOSS_EDIT = ifelse(CONTIG == "chrY" | CONTIG == "Y", -Inf, MINLOSS_EDIT)) %>% #-Inf because whole chrY is just lost 
      mutate(MINLOSS_SUB_EDIT = ifelse(CONTIG == "chrY" | CONTIG == "Y", -Inf,  MINLOSS_SUB_EDIT)) %>% #-Inf because whole chrY is just lost 
      # mutate(MINLOSS_BIAL_EDIT = MINLOSS_BIAL_EDIT) %>% #whole chrY is lost, bial value can be still used for smoothing
      mutate(MINGAIN_EDIT = ifelse(CONTIG == "chrY" | CONTIG == "Y", MINGAIN_FvsPoNM_chrY, MINGAIN_EDIT)) %>%
      mutate(MINGAIN_SUB_EDIT = ifelse(CONTIG == "chrY" | CONTIG == "Y", MINGAIN_SUB_FvsPoNM_chrY, MINGAIN_SUB_EDIT))
  }
  return(DATA)
}
  
#get denoised data and define sex pattern against PoN based cut-offs
SAMPLE1_DENOIS = read.delim(GET_SAMPLE1_DENOIS, comment.char = "@", head=TRUE, stringsAsFactors = F)
SAMPLE1_DENOIS = PROCESS_DENOIS_FCE(SAMPLE1_DENOIS, GET_SAMPLE1_SEX, GET_SAMPLE1_PON_SEX)

if(file.exists(GET_SAMPLE2_DENOIS)) {
  SAMPLE2_DENOIS = read.delim(GET_SAMPLE2_DENOIS, comment.char = "@", head=TRUE, stringsAsFactors = F)
  SAMPLE2_DENOIS = PROCESS_DENOIS_FCE(SAMPLE2_DENOIS, GET_SAMPLE2_SEX, GET_SAMPLE2_PON_SEX)
  }
if(file.exists(GET_SAMPLE3_DENOIS)) {
  SAMPLE3_DENOIS = read.delim(GET_SAMPLE3_DENOIS, comment.char = "@", head=TRUE, stringsAsFactors = F)
  SAMPLE3_DENOIS = PROCESS_DENOIS_FCE(SAMPLE3_DENOIS, GET_SAMPLE3_SEX, GET_SAMPLE3_PON_SEX)
  }

message(paste0("   R ... ", date(), " - STAGE 2/4 - processing CN"))

#CN
#----------------------------------------------------------------------
if(nrow(SEGMENTATION) != 0) {
  
ADD = data.frame(NA,NA,NA,NA,NA)
names(ADD) = c("SAMPLE1_PREDICTION", "SAMPLE1_L2R", "SAMPLE2_L2R", "SAMPLE3_L2R", "PREDICTION")

for(line in 1:nrow(SEGMENTATION)) {
  SEGMENT_CONTIG = SEGMENTATION[line,7] %>% as.character()
  SEGMENT_START  = SEGMENTATION[line,8] %>% as.integer()
  SEGMENT_END    = SEGMENTATION[line,9] %>% as.integer()
  SEGMENT_TYPE   = SEGMENTATION[line,15] %>% as.character()
  
  #L2R of family - IMPORTANT: this L2R is not smooth as it is during segmentation - think about inter-segmentation???
  SAMPLE1 = SAMPLE1_DENOIS %>% filter((CONTIG == SEGMENT_CONTIG) & (START >= SEGMENT_START) & (END <= SEGMENT_END))
  SAMPLE1_L2R = mean(SAMPLE1$LOG2_COPY_RATIO)
  MINLOSS_TEMP = mean(SAMPLE1$MINLOSS_EDIT)
  MINLOSS_SUB_TEMP = mean(SAMPLE1$MINLOSS_SUB_EDIT)
  MINGAIN_TEMP = mean(SAMPLE1$MINGAIN_EDIT)
  MINGAIN_SUB_TEMP = mean(SAMPLE1$MINGAIN_SUB_EDIT)
  SAMPLE1_PREDICTION = ifelse(SAMPLE1_L2R >= MINGAIN_SUB_TEMP, "subGAIN", "NORMAL")
  SAMPLE1_PREDICTION = ifelse(SAMPLE1_L2R <= MINLOSS_SUB_TEMP, "subLOSS", SAMPLE1_PREDICTION)
  SAMPLE1_PREDICTION = ifelse(SAMPLE1_L2R >= MINGAIN_TEMP, "GAIN", SAMPLE1_PREDICTION)
  SAMPLE1_PREDICTION = ifelse(SAMPLE1_L2R <= MINLOSS_TEMP, "LOSS", SAMPLE1_PREDICTION)
    
  # CHECK = ifelse(SAMPLE1_PREDICTION == SEGMENT_TYPE, "pass", "failed")
    
  if(SEGMENT_TYPE == "LOSS" | SEGMENT_TYPE == "GAIN") {
    
    if(file.exists(GET_SAMPLE2_DENOIS)) {
    SAMPLE2_DATA = SAMPLE2_DENOIS %>% filter((CONTIG == SEGMENT_CONTIG) & (START >= SEGMENT_START) & (END <= SEGMENT_END))
    SAMPLE2_L2R = mean(SAMPLE2_DATA$LOG2_COPY_RATIO)
    MINLOSS_TEMP = mean(SAMPLE2_DATA$MINLOSS_EDIT)
    MINLOSS_SUB_TEMP = mean(SAMPLE2_DATA$MINLOSS_SUB_EDIT)
    MINGAIN_TEMP = mean(SAMPLE2_DATA$MINGAIN_EDIT)
    MINGAIN_SUB_TEMP = mean(SAMPLE2_DATA$MINGAIN_SUB_EDIT)
    SAMPLE2_PREDICTION = ifelse(SAMPLE2_L2R >= MINGAIN_TEMP*0.5, "GAIN_like", "NORMAL")
    SAMPLE2_PREDICTION = ifelse(SAMPLE2_L2R >= MINGAIN_TEMP, "GAIN", SAMPLE2_PREDICTION)
    SAMPLE2_PREDICTION = ifelse(SAMPLE2_L2R <= MINLOSS_TEMP, "LOSS", SAMPLE2_PREDICTION)}
    
    if(file.exists(GET_SAMPLE3_DENOIS)) {
    SAMPLE3_DATA = SAMPLE3_DENOIS %>% filter((CONTIG == SEGMENT_CONTIG) & (START >= SEGMENT_START) & (END <= SEGMENT_END))
    SAMPLE3_L2R = mean(SAMPLE3_DATA$LOG2_COPY_RATIO)
    MINLOSS_TEMP = mean(SAMPLE3_DATA$MINLOSS_EDIT)
    MINLOSS_SUB_TEMP = mean(SAMPLE3_DATA$MINLOSS_SUB_EDIT)
    MINGAIN_TEMP = mean(SAMPLE3_DATA$MINGAIN_EDIT)
    MINGAIN_SUB_TEMP = mean(SAMPLE3_DATA$MINGAIN_SUB_EDIT)
    SAMPLE3_PREDICTION = ifelse(SAMPLE3_L2R >= MINGAIN_TEMP*0.5, "GAIN_like", "NORMAL")
    SAMPLE3_PREDICTION = ifelse(SAMPLE3_L2R >= MINGAIN_TEMP, "GAIN", SAMPLE3_PREDICTION)
    SAMPLE3_PREDICTION = ifelse(SAMPLE3_L2R <= MINLOSS_TEMP, "LOSS", SAMPLE3_PREDICTION)}
    
    #predict segment origin by the type of the project
    #---------------------------------------------------------
    #for family - trio
    if(PROJECT_TYPE == "germline" & file.exists(GET_SAMPLE2_DENOIS) & file.exists(GET_SAMPLE3_DENOIS)) {
      #de novo
      PREDICTION = ifelse((SAMPLE2_PREDICTION == "NORMAL") & (SAMPLE3_PREDICTION == "NORMAL"), "de_novo", "unknown")
      #GAIN
      PREDICTION = ifelse((SEGMENT_TYPE == "GAIN") & (SAMPLE2_PREDICTION == "GAIN")      & (SAMPLE3_PREDICTION == "GAIN"),        paste0("inherited_from_both"), PREDICTION)
      PREDICTION = ifelse((SEGMENT_TYPE == "GAIN") & (SAMPLE2_PREDICTION == "GAIN")      & (SAMPLE3_PREDICTION == "GAIN_like"),   paste0("inherited_from_both-likely"), PREDICTION)
      PREDICTION = ifelse((SEGMENT_TYPE == "GAIN") & (SAMPLE2_PREDICTION == "GAIN_like") & (SAMPLE3_PREDICTION == "GAIN"),        paste0("inherited_from_both-likely"), PREDICTION)
      PREDICTION = ifelse((SEGMENT_TYPE == "GAIN") & (SAMPLE2_PREDICTION == "GAIN_like") & (SAMPLE3_PREDICTION == "GAIN_like"),   paste0("inherited_from_both-likely"), PREDICTION)
      PREDICTION = ifelse((SEGMENT_TYPE == "GAIN") & (SAMPLE2_PREDICTION == "GAIN")      & (SAMPLE3_PREDICTION == "NORMAL"),      paste0("inherited_from_", GET_SAMPLE2_ID), PREDICTION)
      PREDICTION = ifelse((SEGMENT_TYPE == "GAIN") & (SAMPLE2_PREDICTION == "GAIN_like") & (SAMPLE3_PREDICTION == "NORMAL"),      paste0("inherited_from_", GET_SAMPLE2_ID, "-likely"), PREDICTION)
      PREDICTION = ifelse((SEGMENT_TYPE == "GAIN") & (SAMPLE2_PREDICTION == "NORMAL")    & (SAMPLE3_PREDICTION == "GAIN"),        paste0("inherited_from_", GET_SAMPLE3_ID), PREDICTION)
      PREDICTION = ifelse((SEGMENT_TYPE == "GAIN") & (SAMPLE2_PREDICTION == "NORMAL")    & (SAMPLE3_PREDICTION == "GAIN_like"),   paste0("inherited_from_", GET_SAMPLE3_ID, "-likely"), PREDICTION)
      #LOSS
      PREDICTION = ifelse((SEGMENT_TYPE == "LOSS") & (SAMPLE2_PREDICTION == "LOSS")      & (SAMPLE3_PREDICTION == "LOSS"),        paste0("inherited_from_both"), PREDICTION)
      PREDICTION = ifelse((SEGMENT_TYPE == "LOSS") & (SAMPLE2_PREDICTION == "LOSS")      & (SAMPLE3_PREDICTION == "NORMAL"),      paste0("inherited_from_", GET_SAMPLE2_ID), PREDICTION)
      PREDICTION = ifelse((SEGMENT_TYPE == "LOSS") & (SAMPLE2_PREDICTION == "NORMAL")    & (SAMPLE3_PREDICTION == "LOSS"),        paste0("inherited_from_", GET_SAMPLE3_ID), PREDICTION)
      # #unknown for gonozomes, as it depends on used controls
      # PREDICTION = ifelse(SEGMENT_CONTIG == "chrX" | SEGMENT_CONTIG == "chrY" | SEGMENT_CONTIG == "X" | SEGMENT_CONTIG == "Y", "unknown", PREDICTION)
    }
    
    #for family - only sample2 available
    if(PROJECT_TYPE == "germline" & file.exists(GET_SAMPLE2_DENOIS) & file.exists(GET_SAMPLE3_DENOIS)==F) {
      SAMPLE3_L2R = "NA"
      PREDICTION = "unknown"
      #GAIN
      PREDICTION = ifelse((SEGMENT_TYPE == "GAIN") & (SAMPLE2_PREDICTION == "GAIN"),       paste0("inherited_from_", GET_SAMPLE2_ID), PREDICTION)
      PREDICTION = ifelse((SEGMENT_TYPE == "GAIN") & (SAMPLE2_PREDICTION == "GAIN_like"),  paste0("inherited_from_", GET_SAMPLE2_ID, "-likely"), PREDICTION)
      #LOSS
      PREDICTION = ifelse((SEGMENT_TYPE == "LOSS") & (SAMPLE2_PREDICTION == "LOSS"),       paste0("inherited_from_", GET_SAMPLE2_ID), PREDICTION)
    }
    
    #for family - only sample3 available
    if(PROJECT_TYPE == "germline" & file.exists(GET_SAMPLE2_DENOIS)==F & file.exists(GET_SAMPLE3_DENOIS)) {
      SAMPLE2_L2R = "NA"
      PREDICTION = "unknown"
      #GAIN
      PREDICTION = ifelse((SEGMENT_TYPE == "GAIN") & (SAMPLE3_PREDICTION == "GAIN"),        paste0("inherited_from_", GET_SAMPLE3_ID), PREDICTION)
      PREDICTION = ifelse((SEGMENT_TYPE == "GAIN") & (SAMPLE3_PREDICTION == "GAIN_like"),   paste0("inherited_from_", GET_SAMPLE3_ID, "-likely"), PREDICTION)
      #LOSS
      PREDICTION = ifelse((SEGMENT_TYPE == "LOSS") & (SAMPLE3_PREDICTION == "LOSS"),        paste0("inherited_from_", GET_SAMPLE3_ID), PREDICTION)
    }
    
    #for family - singleton
    if(PROJECT_TYPE == "germline" & file.exists(GET_SAMPLE2_DENOIS)==F & file.exists(GET_SAMPLE3_DENOIS)==F) {
      SAMPLE2_L2R = "NA"
      SAMPLE3_L2R = "NA"
      PREDICTION = "unknown"
    }
    
    #for tumor - germline sample available
    if(PROJECT_TYPE == "tumor" & file.exists(GET_SAMPLE2_DENOIS) & file.exists(GET_SAMPLE3_DENOIS)==F) {
      SAMPLE3_L2R = "NA"
      PREDICTION = "unknown"
      #GAIN
      PREDICTION = ifelse((SEGMENT_TYPE == "GAIN") & (SAMPLE2_PREDICTION == "GAIN_like"),   "germline-likely", PREDICTION)
      PREDICTION = ifelse((SEGMENT_TYPE == "GAIN") & (SAMPLE2_PREDICTION == "GAIN"),        "germline", PREDICTION)
      PREDICTION = ifelse((SEGMENT_TYPE == "GAIN") & (SAMPLE2_PREDICTION == "NORMAL"),      "somatic", PREDICTION)
      #LOSS
      PREDICTION = ifelse((SEGMENT_TYPE == "LOSS") & (SAMPLE2_PREDICTION == "LOSS"),        "germline", PREDICTION)
      PREDICTION = ifelse((SEGMENT_TYPE == "LOSS") & (SAMPLE2_PREDICTION == "NORMAL"),      "somatic", PREDICTION)
    }
    
    #for tumor - germline sample NOT available
    if(PROJECT_TYPE == "tumor" & file.exists(GET_SAMPLE2_DENOIS)==F & file.exists(GET_SAMPLE3_DENOIS)==F) {
      SAMPLE2_L2R = "NA"
      SAMPLE3_L2R = "NA"
      PREDICTION = "unknown"
    }
  }
  
    #for subclonal segments
    if(SEGMENT_TYPE == "subLOSS" | SEGMENT_TYPE == "subGAIN") {
      SAMPLE2_L2R = "NA"
      SAMPLE3_L2R = "NA"
      PREDICTION = "unknown"
    }
  
  x = data.frame(SAMPLE1_PREDICTION, SAMPLE1_L2R, SAMPLE2_L2R, SAMPLE3_L2R, PREDICTION)
  
  ADD = rbind(ADD, x)
 }

ADD = filter(ADD, SAMPLE1_L2R != "NA")
#----------------------------------------------------------------------
SEGMENTATION = cbind(SEGMENTATION, ADD)
#----------------------------------------------------------------------
}
message(paste0("   R ... ", date(), " - STAGE 3/4 - processing LOH"))


#LOH
#----------------------------------------------------------------------
if(nrow(SEGMENTATION_LOH) != 0) {

ADD_LOH = data.frame(NA,NA,NA,NA,NA)
names(ADD_LOH) = c("SAMPLE1_PREDICTION", "SAMPLE1_L2R", "SAMPLE2_L2R", "SAMPLE3_L2R", "PREDICTION")

for(line in 1:nrow(SEGMENTATION_LOH)) {
  SEGMENT_CONTIG = SEGMENTATION_LOH[line,7] %>% as.character()
  SEGMENT_START  = SEGMENTATION_LOH[line,8] %>% as.integer()
  SEGMENT_END    = SEGMENTATION_LOH[line,9] %>% as.integer()
  SEGMENT_TYPE   = SEGMENTATION_LOH[line,15] %>% as.character
  
  SAMPLE1 = SAMPLE1_DENOIS %>% filter((CONTIG == SEGMENT_CONTIG) & (START >= SEGMENT_START) & (END <= SEGMENT_END))
  SAMPLE1_L2R = mean(SAMPLE1$LOG2_COPY_RATIO)
  SAMPLE1_PREDICTION = ifelse(SAMPLE1_L2R <= MINGAIN & SAMPLE1_L2R >= MINLOSS, "cnnLOH", "other")
  SAMPLE1_PREDICTION = ifelse(SAMPLE1_L2R <= MINLOSS, "LOH-LOSS", SAMPLE1_PREDICTION)
  
  if(file.exists(GET_SAMPLE2_DENOIS)) {
  SAMPLE2_DATA = SAMPLE2_DENOIS %>% filter((CONTIG == SEGMENT_CONTIG) & (START >= SEGMENT_START) & (END <= SEGMENT_END))
  SAMPLE2_L2R = mean(SAMPLE2_DATA$LOG2_COPY_RATIO)}
  if(file.exists(GET_SAMPLE2_DENOIS)==F) {SAMPLE2_L2R = "NA"}
  
  if(file.exists(GET_SAMPLE3_DENOIS)) {
  SAMPLE3_DATA = SAMPLE3_DENOIS %>% filter((CONTIG == SEGMENT_CONTIG) & (START >= SEGMENT_START) & (END <= SEGMENT_END))
  SAMPLE3_L2R = mean(SAMPLE3_DATA$LOG2_COPY_RATIO)}
  if(file.exists(GET_SAMPLE3_DENOIS)==F) {SAMPLE3_L2R = "NA"}
  
  PREDICTION = "NA"

  x = data.frame(SAMPLE1_PREDICTION, SAMPLE1_L2R, SAMPLE2_L2R, SAMPLE3_L2R, PREDICTION)

  ADD_LOH = rbind(ADD_LOH, x)

}

ADD_LOH = filter(ADD_LOH, SAMPLE1_L2R != "NA")
SEGMENTATION_LOH = cbind(SEGMENTATION_LOH, ADD_LOH)
}

if(nrow(SEGMENTATION) != 0 & nrow(SEGMENTATION_LOH) != 0) {SEGMENTATION = rbind(SEGMENTATION, SEGMENTATION_LOH)}
if(nrow(SEGMENTATION) != 0 & nrow(SEGMENTATION_LOH) == 0) {SEGMENTATION = SEGMENTATION}
if(nrow(SEGMENTATION) == 0 & nrow(SEGMENTATION_LOH) != 0) {SEGMENTATION = SEGMENTATION_LOH}


if(nrow(SEGMENTATION) != 0) {
if(PROJECT_TYPE == "germline") {names(SEGMENTATION) = c(colnames(SEGMENTATION_INPUT), "TYPE_unsmooth", "SAMPLE1_L2R_unsmooth", "SAMPLE2_L2R_unsmooth", "SAMPLE3_L2R_unsmooth", "INHERITANCE_PREDICTION")}
if(PROJECT_TYPE == "tumor")    {names(SEGMENTATION) = c(colnames(SEGMENTATION_INPUT), "TYPE_unsmooth", "TUMOR_L2R_unsmooth", "GERMLINE_L2R_unsmooth", "OTHER_L2R_unsmooth", "ORIGIN_PREDICTION")}
}

if(nrow(SEGMENTATION) == 0 & nrow(SEGMENTATION_LOH) == 0) {SEGMENTATION = data.table("no abnormal segments to report")}



message(paste0("   R ... ", date(), " - STAGE 4/4 - writing output"))

#Write output
#----------------------------------------------------------------------------------------------------------------------
write_tsv(SEGMENTATION, paste0(OUTPUT), col_names = T)
  
message(paste0("   R ... ", date(), " - FINISHED"))

system(paste0("touch status_ok"))
