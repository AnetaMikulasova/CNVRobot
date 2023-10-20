args <- commandArgs(trailing = TRUE)

message(paste0("   R ... ", date(), " - STARTING"))

suppressPackageStartupMessages(library(ASCAT))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(plyranges))
suppressPackageStartupMessages(library(karyoploteR))

FRACTION           = args[1] %>% as.numeric()
GET_ID             = args[2] %>% as.character()
GET_TUM_SEGM       = args[3] %>% as.character()
GET_TUM_DENOIS     = args[4] %>% as.character()
GET_TUM_SNP        = args[5] %>% as.character()
GET_GER_DENOIS     = args[6] %>% as.character()
GET_GER_SNP        = args[7] %>% as.character()
OUTPUT_ABSOLUTE_CN_SEGM = args[8] %>% as.character()
OUTPUT_ABSOLUTE_CN_DENOIS = args[9] %>% as.character()
OUTPUT_ASCAT       = args[10] %>% as.character()
OUTPUT_DIR         = args[11] %>% as.character()
OUTDIR_IGV         = args[12] %>% as.character()

message(paste0("   R ... ", date(), " - STAGE 1/5 - loading and preparing input data"))

TUMOR_SEGMENTATION = read.delim(paste0(GET_TUM_SEGM), stringsAsFactors = F)
TUMOR_DENOIS = read.delim(GET_TUM_DENOIS, comment.char = "@", stringsAsFactors = F)

PROCESS_LOG_TO_SNP = function(DENOIS_DATA, SNP_DATA) {
  DENOIS = read.delim(DENOIS_DATA, comment.char = "@", stringsAsFactors = F)
  GR_DENOIS=DENOIS
  names(GR_DENOIS) = c("seqnames", "start", "end", "log")
  GR_DENOIS=GR_DENOIS %>% toGRanges()

  ALL = readRDS(SNP_DATA) %>%
    dplyr::mutate(N = REF_COUNT + ALT_COUNT)
  ALL = ALL %>%
    dplyr::mutate(TEST = sample(1:2, nrow(ALL), replace = T)) %>%
    dplyr::mutate(MAF = ALT_COUNT/N) %>%
    dplyr::select(CONTIG, POSITION, MAF)
  GR_SNP= ALL %>%
    dplyr::mutate(POSITION2 = POSITION) %>%
    dplyr::select(CONTIG, POSITION, POSITION2, MAF)
  names(GR_SNP) = c("seqnames", "start", "end", "maf")
  GR_SNP=GR_SNP %>% toGRanges()

  DATA = join_overlap_inner(GR_SNP, GR_DENOIS)
  return(DATA)
}

TUMOR_DATA = PROCESS_LOG_TO_SNP(GET_TUM_DENOIS, GET_TUM_SNP)
GERM_DATA = PROCESS_LOG_TO_SNP(GET_GER_DENOIS, GET_GER_SNP)

#subsets
TUMOR_DATA_FILTER=subsetByOverlaps(TUMOR_DATA, GERM_DATA)
GER_DATA_FILTER=subsetByOverlaps(GERM_DATA, TUMOR_DATA)

#just training if too much SNPs to reduce time by analysing random fraction
SUBSET=TUMOR_DATA %>% as.data.frame() %>% sample_frac(FRACTION)
SUBSET = SUBSET %>% toGRanges()
#subset
TUMOR_DATA_FILTER=subsetByOverlaps(TUMOR_DATA_FILTER, SUBSET)
GER_DATA_FILTER=subsetByOverlaps(GER_DATA_FILTER, SUBSET)

GIVE_ID = function(DATA, ID) {
  names(DATA) = c(" ", "chrs", "pos", ID)
  return(DATA)
}

max = nrow(TUMOR_DATA_FILTER %>% as.data.frame)
#run ASCAT
TUM_LOG=TUMOR_DATA_FILTER %>% as.data.frame() %>% mutate(ID = paste0("SNP", 1:max)) %>% select(ID, seqnames, start, log) %>% GIVE_ID(GET_ID)
TUM_SNP=TUMOR_DATA_FILTER %>% as.data.frame() %>% mutate(ID = paste0("SNP", 1:max)) %>% select(ID, seqnames, start, maf) %>% GIVE_ID(GET_ID)
GER_LOG=GER_DATA_FILTER %>% as.data.frame() %>% mutate(ID = paste0("SNP", 1:max)) %>% select(ID, seqnames, start, log) %>% GIVE_ID(GET_ID)
GER_SNP=GER_DATA_FILTER %>% as.data.frame() %>% mutate(ID = paste0("SNP", 1:max)) %>% select(ID, seqnames, start, maf) %>% GIVE_ID(GET_ID)

OUTDIR=paste0(OUTPUT_DIR, "/")
setwd(OUTDIR)

write_tsv(TUM_LOG, paste0(OUTDIR, "tum_log.txt"), col_names = T)
write_tsv(TUM_SNP, paste0(OUTDIR, "tum_snp.txt"), col_names = T)
write_tsv(GER_LOG, paste0(OUTDIR, "ger_log.txt"), col_names = T)
write_tsv(GER_SNP, paste0(OUTDIR, "ger_snp.txt"), col_names = T)

CONTIGS = TUMOR_DATA_FILTER %>% as.data.frame()
CONTIGS=CONTIGS$seqnames %>% as.vector() %>% unique()


message(paste0("   R ... ", date(), " - STAGE 2/5 - running ascat"))

ascat.bc = ascat.loadData("tum_log.txt", "tum_snp.txt", "ger_log.txt", "ger_snp.txt", chrs = CONTIGS)
ascat.plotRawData(ascat.bc)
ascat.bc = ascat.aspcf(ascat.bc)
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc, gamma = 1)

message(paste0("   R ... ", date(), " - STAGE 3/5 - calculating absolute CN"))

# median_normalized_log_ratio=median(TUMOR_DENOIS$LOG2_COPY_RATIO)
purity = ascat.output$aberrantcellfraction %>% as.numeric()
ploidy = ascat.output$ploidy %>% as.numeric()

# TUMOR_SEGMENTATION = TUMOR_SEGMENTATION %>%
#   mutate(purity_adj_log_ratio = MEAN_L2R - median_normalized_log_ratio + log((purity*ploidy + (1-purity)*2)/2,2)) %>%
#   mutate(ABS_CN = (2^(1+purity_adj_log_ratio) - 2* (1-purity))/purity) %>%
#   mutate(ABS_CN_EDIT = ifelse(ABS_CN > 5, 5, ABS_CN))
#
# TUMOR_DENOIS = TUMOR_DENOIS %>%
#   mutate(purity_adj_log_ratio = LOG2_COPY_RATIO - median_normalized_log_ratio + log((purity*ploidy + (1-purity)*2)/2,2)) %>%
#   mutate(ABS_CN = (2^(1+purity_adj_log_ratio) - 2* (1-purity))/purity) %>%
#   mutate(ABS_CN_EDIT = ifelse(ABS_CN > 5, 5, ABS_CN))

TUMOR_SEGMENTATION = TUMOR_SEGMENTATION %>%
  mutate(ABS_CN = ((2^MEAN_L2R - (1-purity))/purity)*ploidy) %>%
  mutate(ABS_CN_R = round(ABS_CN)) %>%
  mutate(ABS_CN_R_EDIT = ABS_CN_R)
  # mutate(ABS_CN_R_EDIT = ifelse(ABS_CN_R > 5, 5, ABS_CN_R))

TUMOR_DENOIS = TUMOR_DENOIS %>%
  mutate(ABS_CN = ((2^LOG2_COPY_RATIO - (1-purity))/purity)*ploidy) %>%
  mutate(ABS_CN_R = round(ABS_CN)) %>%
  mutate(ABS_CN_R_EDIT = ABS_CN_R)
  # mutate(ABS_CN_R_EDIT = ifelse(ABS_CN_R > 5, 5, ABS_CN_R))

message(paste0("   R ... ", date(), " - STAGE 4/5 - writing outpus"))

saveRDS(ascat.output, file = "OUTPUT_ASCAT")
write_tsv(TUMOR_SEGMENTATION, paste0(OUTPUT_ABSOLUTE_CN_SEGM), col_names = T)
write_tsv(TUMOR_DENOIS, paste0(OUTPUT_ABSOLUTE_CN_DENOIS), col_names = T)

#cat(c(purity, ploidy), sep="_")


message(paste0("   R ... ", date(), " - STAGE 5/5 - creating IGV files"))

#1) CNVRobot segments with absolute CN
#bed file for IGV with abnormal segments
IGV_MODEL = TUMOR_SEGMENTATION %>% 
  # filter(FILTER == "PASS") %>%
  select(CONTIG, START, END, TYPE, ABS_CN_R_EDIT) %>%
  mutate(START = START-1) %>%
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

IGV_MODEL_BED_CN = IGV_MODEL %>% 
  filter(TYPE != "NORMAL") %>%
  filter(TYPE != "HOM") %>%
  filter(TYPE != "OTHER") %>%
  mutate(TYPE = paste0(TYPE, "_", ABS_CN_R_EDIT)) %>%
  select(CONTIG, START, END, TYPE, x1, x2, START2, END2, COLOR)
write_tsv(IGV_MODEL_BED_CN, paste0(OUTDIR_IGV, GET_ID, "_abnormal_CN_segments_absCN.bed"), col_names = F)


IGV_MODEL_WIG = IGV_MODEL %>% 
  filter(TYPE != "HOM") %>%
  filter(TYPE != "OTHER") %>%
  select(CONTIG, START, END, ABS_CN_R_EDIT)
write_tsv(IGV_MODEL_WIG, paste0(OUTDIR_IGV, GET_ID, "_segments_absCN"), col_names = F)

# #2) CNVRobot de ois with absolute CN
# #bedGraph file with denoised copy-number data
# IGV_DENOIS  <- TUMOR_DENOIS %>%
#   select(CONTIG, START, END, ABS_CN_R_EDIT) %>%
#   mutate(START = START-1) %>%
#   arrange(CONTIG, START, END)
# write_tsv(IGV_DENOIS, paste0(OUTDIR_IGV, GET_ID, "_CN_absCN"), col_names = F)

#3) ASCAT segments
ASCAT_MODEL = ascat.output$segments %>%
  mutate(total = nMajor + nMinor) %>%
  mutate(nMajor_edit = ifelse(nMajor > 5, 5, nMajor)) %>%
  mutate(nMinor_edit = ifelse(nMinor > 5, 5, nMinor)) %>%
  mutate(total_edit = ifelse(total > 5, 5, total)) %>%
  mutate(startpos = startpos-1) %>%
  arrange(chr, startpos, endpos)
  
MAJOR = ASCAT_MODEL %>% select(chr, startpos, endpos, nMajor)
MINOR = ASCAT_MODEL %>% select(chr, startpos, endpos, nMinor)
TOTAL = ASCAT_MODEL %>% select(chr, startpos, endpos, total)

write_tsv(MAJOR, paste0(OUTDIR_IGV, GET_ID, "_nMajor_allele"), col_names = F)
write_tsv(MINOR, paste0(OUTDIR_IGV, GET_ID, "_nMinor_allele"), col_names = F)
write_tsv(TOTAL, paste0(OUTDIR_IGV, GET_ID, "_nTotal_allele"), col_names = F)








