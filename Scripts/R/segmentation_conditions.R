args <- commandArgs(trailing = TRUE)

# message(paste0("   R ... ", date(), " - STARTING"))

suppressPackageStartupMessages(library(tidyverse))

# message(paste0("   R ... ", date(), " - STAGE 1/3 - loading segmentation table"))
# GET_SEGMENTATION = "/Users/aneta/Documents/projects/GATK_CN/Masters_workstation/setting_segmentation_temp.txt"
# GET_SEGMENTATION_ID = "SEGM-BrnoA"

GET_SEGMENTATION = args[1] %>% as.character()
GET_SEGMENTATION_ID = args[2] %>% as.character()
GET_SEGMENTATION_ID_OUT = args[3] %>% as.character()

# message(paste0("   R ... ", date(), " - STAGE 2/3 - processing segmentation table"))
SEGMENTATION = read.delim(GET_SEGMENTATION, header = TRUE, stringsAsFactors = F, check.names = F) %>%
  select(FEATURE, matches(paste0('^', GET_SEGMENTATION_ID, '$')))

# message(paste0("   R ... ", date(), " - STAGE 3/3 - writing output"))
write_tsv(SEGMENTATION, path = paste0(GET_SEGMENTATION_ID_OUT), col_names = T)

# message(paste0("   R ... ", date(), " - FINISHED"))
