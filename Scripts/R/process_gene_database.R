args <- commandArgs(trailing = TRUE)

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(splitstackshape))

GET_GENES_DTB            =  args[1] %>% as.character()
GET_GENES_DTB_PROCESSED  =  args[2] %>% as.character()

GENES_DTB = read.delim(GET_GENES_DTB, header = T, stringsAsFactors = F) %>%
  mutate(ID = paste0(gene_id, "_", chrom, "_", strand ))

LISTofID = unique(GENES_DTB$ID) %>% as.list()

GENES_DTB_PROCESSED=data.frame()

for (each in LISTofID) {
  
  TABLE = GENES_DTB %>% filter(ID == each)
  
  #first, tx start and end
  GENE_ID = TABLE$gene_id %>% unique() %>% as.character()
  CONTIG = TABLE$chrom %>% unique() %>% as.character()
  STRAND = TABLE$strand %>% unique() %>% as.character()
  TxSTART = min(TABLE$txStart) %>% as.integer()
  TxEND = max(TABLE$txEnd) %>% as.integer()
  
  #unique overall exons and coding exons
  EXONS = data.frame(NA,NA)
  names(EXONS) = c("start", "end")
  
  EXONS_CODING = data.frame(NA,NA)
  names(EXONS_CODING) = c("start", "end")
  line=5
  for (line in 1:nrow(TABLE)) {
    cdcSTART = TABLE %>% select(cdsStart) %>% slice(line) %>% as.integer()
    cdcEND = TABLE %>% select(cdsEnd) %>% slice(line) %>% as.integer()
    STARTS = TABLE %>% select(exonStarts) %>% slice(line) %>% as.character()
    STARTS = unlist(strsplit(STARTS, ",")) %>% as.data.frame() %>% mutate_if(is.factor, as.character) %>% mutate_if(is.character, as.integer)
    ENDS = TABLE %>% select(exonEnds) %>% slice(line) %>% as.character()
    ENDS = unlist(strsplit(ENDS, ",")) %>% as.data.frame() %>% mutate_if(is.factor, as.character) %>% mutate_if(is.character, as.integer)

    POSITIONS = cbind(STARTS, ENDS)
    names(POSITIONS) = c("start", "end")
    
    EXONS = rbind(EXONS, POSITIONS)
    
    POSITIONS2 = POSITIONS %>% mutate(start = ifelse(start < cdcSTART, cdcSTART, start)) %>%
      mutate(end = ifelse(end > cdcEND, cdcEND, end)) %>%
      mutate(size = end - start) %>%
      filter(size > 0) %>%
      select(start, end)
    
    EXONS_CODING = rbind(EXONS_CODING, POSITIONS2)

  }
  
  EXONS = filter(EXONS, start !="NA") %>% unique() %>% arrange(start, end)
  NofEXONS = nrow(EXONS)
  EXON_STARTS = paste((EXONS$start), collapse=",")
  EXON_ENDS = paste((EXONS$end), collapse=",")

  EXONS_CODING = filter(EXONS_CODING, start !="NA") %>% unique() %>% arrange(start, end)
  NofEXONS_CODING = nrow(EXONS_CODING)
  EXON_STARTS_CODING = if (NofEXONS_CODING > 0) {EXON_STARTS_CODING = paste((EXONS_CODING$start), collapse=",")} else {EXON_STARTS_CODING = "NA"}
  EXON_ENDS_CODING = if (NofEXONS_CODING > 0) {EXON_ENDS_CODING = paste((EXONS_CODING$end), collapse=",")} else {EXON_ENDS_CODING = "NA"}
  
  EXPORTED_TABLE = data.frame(GENE_ID, CONTIG, STRAND, TxSTART, TxEND, NofEXONS, EXON_STARTS, EXON_ENDS, NofEXONS_CODING, EXON_STARTS_CODING, EXON_ENDS_CODING)
  
  GENES_DTB_PROCESSED = rbind(GENES_DTB_PROCESSED, EXPORTED_TABLE)
} 

write_tsv(GENES_DTB_PROCESSED, paste0(GET_GENES_DTB_PROCESSED), col_names = T)


  




