args <- commandArgs(trailing = TRUE)

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(splitstackshape))
suppressPackageStartupMessages(library(expss))
suppressPackageStartupMessages(library(IRanges))

GET_GENES_DTB            =  args[1] %>% as.character()
GET_GENES_DTB_PROCESSED  =  args[2] %>% as.character()

GENES_DTB = read.delim(GET_GENES_DTB, header = T, stringsAsFactors = F) %>%
  mutate(ID = paste0(gene_id, "_", chrom, "_", strand ))

LISTofID = unique(GENES_DTB$ID) %>% as.list()

GENES_DTB_PROCESSED=data.frame()

for (each in LISTofID) {
  
  message(date(), " --- processing gene number ", which(unlist(LISTofID)==each), " out of ", length(LISTofID))

  TABLE = GENES_DTB %>% filter(ID == each)
  
  GROUPS = IRanges(TABLE$txStart, TABLE$txEnd)
  TABLE$PART = subjectHits(findOverlaps(GROUPS, reduce(GROUPS)))

  LISTofID2 = unique(TABLE$PART) %>% as.list()
  
  for (each2 in LISTofID2) {

    TABLE2 = TABLE %>% filter(PART == each2)
    
    #first, tx start and end
    GENE_ID = TABLE2$gene_id %>% unique() %>% as.character()
    CONTIG = TABLE2$chrom %>% unique() %>% as.character()
    STRAND = TABLE2$strand %>% unique() %>% as.character()
    TxSTART = min(TABLE2$txStart) %>% as.integer()
    TxEND = max(TABLE2$txEnd) %>% as.integer()
    
    #unique overall exons and coding exons
    EXONS = data.frame(NA,NA)
    names(EXONS) = c("start", "end")
    
    EXONS_CODING = data.frame(NA,NA)
    names(EXONS_CODING) = c("start", "end")
    line=5
    for (line in 1:nrow(TABLE2)) {
      cdcSTART = TABLE2[line,] %>% select(cdsStart) %>% as.integer()
      cdcEND = TABLE2[line,] %>% select(cdsEnd) %>% as.integer()
      STARTS = TABLE2[line,] %>% select(exonStarts) %>% as.character()
      STARTS = unlist(strsplit(STARTS, ",")) %>% as.data.frame() %>% mutate_if(is.factor, as.character) %>% mutate_if(is.character, as.integer)
      ENDS = TABLE2[line,] %>% select(exonEnds) %>% as.character()
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
}

write_tsv(GENES_DTB_PROCESSED, paste0(GET_GENES_DTB_PROCESSED), col_names = T)

system(paste0("touch status_ok"))
  




