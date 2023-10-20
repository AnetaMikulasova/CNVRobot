args <- commandArgs(trailing = TRUE)

message(paste0("   R ... ", date(), " - STARTING"))

# # You need to install the Bioconductor BiocManager first
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# # Install rhdf5 package
# BiocManager::install("rhdf5")

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rhdf5))
suppressPackageStartupMessages(library(expss))

GET_FILE          = args[1] %>% as.character()
OUTPUT            = args[2] %>% as.character()

#load data from hdf5 file
#----------------------------------------------------------------
message(paste0("   R ... ", date(), " - STAGE 1/3 - loading data"))

# FILE_LS = h5ls(GET_FILE)
COUNTS = h5read(GET_FILE, "/counts/values")
INTERVALS = h5read(GET_FILE, "/intervals/transposed_index_start_end")
CONTIGS = h5read(GET_FILE, "/intervals/indexed_contig_names") %>% as.data.frame() 
names(CONTIGS) = c("CONTIG")
CONTIGS = CONTIGS %>% mutate(INDEX = row_number()-1)
#----------------------------------------------------------------

#double counts for chrX and chrY
#----------------------------------------------------------------
message(paste0("   R ... ", date(), " - STAGE 2/3 - processing data"))

DATA = cbind(INTERVALS, COUNTS) %>% as.data.frame()
names(DATA) = c("INDEX", "START", "END", "COUNT")
DATA = DATA %>% 
  mutate(CONTIG = vlookup(INDEX, dict = CONTIGS, lookup_column = which(colnames(CONTIGS)=="INDEX"), result_column = which(colnames(CONTIGS)=="CONTIG"))) %>%
  mutate(COUNT_EDIT = ifelse(CONTIG == "chrX" | CONTIG == "X" | CONTIG == "chrY" | CONTIG == "Y", COUNT*2, COUNT))

COUNTS_EDIT = DATA %>% select(COUNT_EDIT) %>% as.matrix()
colnames(COUNTS_EDIT) <- NULL
COUNTS_EDIT = matrix(COUNTS_EDIT, ncol=1)
#----------------------------------------------------------------

#write output
#----------------------------------------------------------------
message(paste0("   R ... ", date(), " - STAGE 3/3 - writing output"))
system(paste0("cp ", GET_FILE, " ", OUTPUT))
h5write(COUNTS_EDIT, OUTPUT, "/counts/values")
#----------------------------------------------------------------
message(paste0("   R ... ", date(), " - FINISHED"))
system(paste0("touch status_ok"))


