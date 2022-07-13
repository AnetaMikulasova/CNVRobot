# 
# MYC	chr8	124987758	129987754	4999996
# IGH	chr14	105533663	107000000	1466337
# IGK	chr2	88200481	90700481	2500000
# IGL	chr22	21645711	23657813	2012102
# CCND1	chr11	69032532	69785232	752700
# MAF	chr16	77466103	79966103	2500000
# FGFR3NSD2	chr4	1698273	2098273	400000

#specific annotation
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------
#Get IGH annotation
# DTB_IG_ANNOT = read.delim("/Users/aneta/Documents/projects/IGH_project/bed_files/gencode_genes_hg38_V31_IGinColor.bed", header = FALSE, stringsAsFactors = F)
DTB_IG_ANNOT = read.delim("/media/aneta/Data/projects/IGH_project/bed_files/gencode_genes_hg38_V31_IGinColor.bed", header = FALSE, stringsAsFactors = F)
names(DTB_IG_ANNOT) = c("CONTIG", "START1", "END1", "ID", "x", "y", "START2", "END2", "COLORORIG")
DTB_IG_ANNOT = DTB_IG_ANNOT %>%
  filter(ID != "IGLL5") %>% #delete IGLL5 fragment
  filter(ID != "IGLL1") %>%
  filter(ID != "IGKV3OR2-268") %>%
  filter(ID != "IGKV1OR2-108") %>%
  mutate(COLOR = ifelse(COLORORIG == "240,50,80", "#CD1C00", ".")) %>%
  mutate(COLOR = ifelse(COLORORIG == "0,170,220", "#4F94CD", COLOR)) %>%
  mutate(COLOR = ifelse(COLORORIG == "100,180,50", "#548B54", COLOR)) %>%
  mutate(COLOR = ifelse(COLORORIG == "100,70,230", "#8B4789", COLOR)) %>%
  mutate(COLOR = ifelse(COLORORIG == "120,120,120", "#787878", COLOR)) %>%
  mutate(POSITION = (START1+((END1-START1)/2))) %>%
  mutate(ID_SHORT = gsub("IGH|IGL|IGK", "", ID))

# DTB_IG_ANNOT_REG = read.delim("/Users/aneta/Documents/projects/IGH_project/bed_files/gencode_genes_hg38_V31_IG_VDJC.bed", header = FALSE, stringsAsFactors = F)
DTB_IG_ANNOT_REG = read.delim("/media/aneta/Data/projects/IGH_project/bed_files/gencode_genes_hg38_V31_IG_VDJC.bed", header = FALSE, stringsAsFactors = F)
names(DTB_IG_ANNOT_REG) = c("CONTIG", "START1", "END1", "ID", "x", "y", "START2", "END2", "COLORORIG")
DTB_IG_ANNOT_REG = DTB_IG_ANNOT_REG %>%
  filter(ID != "IGLC" & ID != "IGLJ") %>%
  mutate(COLOR = ifelse(COLORORIG == "240,50,80", "#CD1C00", ".")) %>%
  mutate(COLOR = ifelse(COLORORIG == "0,170,220", "#4F94CD", COLOR)) %>%
  mutate(COLOR = ifelse(COLORORIG == "100,180,50", "#548B54", COLOR)) %>%
  mutate(COLOR = ifelse(COLORORIG == "100,70,230", "#8B4789", COLOR)) %>%
  mutate(COLOR = ifelse(COLORORIG == "0,0,0", "#CD1C00", COLOR)) %>%
  mutate(POSITION = (START1+((END1-START1)/2))) %>%
  mutate(IDSHORT = ifelse(ID == "IGHC" | ID == "IGKC" | ID == "IGLC", "C", ".")) %>%
  mutate(IDSHORT = ifelse(ID == "IGHD" | ID == "IGKD" | ID == "IGLD", "D", IDSHORT)) %>%
  mutate(IDSHORT = ifelse(ID == "IGHJ" | ID == "IGKJ" | ID == "IGLJ", "J", IDSHORT)) %>%
  mutate(IDSHORT = ifelse(ID == "IGHV" | ID == "IGKV" | ID == "IGLV", "V", IDSHORT)) %>%
  mutate(IDSHORT = ifelse(ID == "IGLCJ", "J-C", IDSHORT))


#Blueprits regions
# DTB_BLUEPRINT = read.delim("/Users/aneta/Documents/projects/IGH_project/blueprint_ChIPseq_chromatin-state/get_ChIPseq_consensus/regions_results_NEW_NEW_hg38.bed", header = FALSE, stringsAsFactors = F) %>%
DTB_BLUEPRINT = read.delim("/media/aneta/Data/projects/IGH_project/bed_files/regions_results_NEW_NEW_hg38.bed", header = FALSE, stringsAsFactors = F) %>%
  select(V1, V2, V3, V4) %>%
  mutate(V2 = as.integer(V2+1))
names(DTB_BLUEPRINT) = c("CONTIG", "START", "END", "ID")

#Get S-regions
# DTB_SREG = read.delim("/Users/aneta/Documents/projects/IGH_project/bed_files/SupplFigx_AID_200to2500andmore_hg38.bed",  header = FALSE)
DTB_SREG = read.delim("/media/aneta/Data/projects/IGH_project/bed_files/SupplFigx_AID_200to2500andmore_hg38.bed",  header = FALSE)
names(DTB_SREG) = c("CONTIG", "START1", "END1", "ID", "x", "y", "START2", "END2", "COLORORIG")
DTB_SREG = DTB_SREG %>% mutate(POSITION = (START1+((END1-START1)/2)))

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------


#define regions of interest
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------
#Create data frame with regions of interest - regions without extra annotation
REGIONSofINTEREST = data.frame(NA, NA, NA, NA)
MYC = c("MYC-region", "chr8", 124987758, 129987754)
# MYC = c("MYC-region", "chr8", 126000000, 129500000)
CCND1 = c("CCND1-region", "chr11", 69032532, 69785232)
MAF = c("MAF-region", "chr16", 77466103, 79966103)
FGFR3NSD2 = c("FGFR3NSD2-region", "chr4", 1698273, 2098273)
IGH = c("IGH-region", "chr14", 105533663, 107000000)
IGK = c("IGK-region", "chr2", 88200481, 90700481)
IGL = c("IGL-region", "chr22", 21645711, 23657813)


REGIONSofINTEREST = rbind(REGIONSofINTEREST, MYC, CCND1, MAF, FGFR3NSD2, IGH, IGK, IGL)
# REGIONSofINTEREST = rbind(REGIONSofINTEREST, IGH, IGK, IGL)

names(REGIONSofINTEREST) = c("ID", "CONTIG", "START", "END")
REGIONSofINTEREST = REGIONSofINTEREST %>% filter(CONTIG != "NA")

REGIONSofINTEREST = mutate(REGIONSofINTEREST, N = 1:n())
LISTofREGIONS = unique(REGIONSofINTEREST$N) %>% as.list()
NofREG = max(unlist(LISTofREGIONS))
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------

#run plots
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------
#same as detailed plot, only some small changes - no FROMABN, TOABN, SIZEABN, TYPEABN, PROBE, CHRBAND (in output= and main= )
for (REG in (LISTofREGIONS)) {
  
  options("scipen"=100)
  
  REGION = REGIONSofINTEREST %>%
    filter(N == REG)
  
  CHR = REGION$CONTIG %>% as.character()
  FROM = REGION$START %>% as.numeric()
  TO = REGION$END %>% as.numeric()
  SIZE = TO-FROM
  REGID = REGION$ID %>% as.character()
  
  message(paste0("   R ... ", date(), " - STAGE 5+/5 - project specific plot: ", REG, "/", NofREG, " - ", REGID))
  
  output = paste0(OUTDIR, "project_specific_regions/", GROUP, "_", SAMPLE1_TYPE, "_", SAMPLE1_ID, "_", SAMPLE1_SEX, "_", CHR, "-", FROM, "-", TO, "_", REGID, "_size-", SIZE)
  
  ZOOM = toGRanges(data.frame(CHR, FROM, TO))
  IDPOS = FROM+((TO-FROM)/2)
  
  #Data function
  DETAIL_CHART = function (ID, DATADENOIS, DATADENOISABN, DATAMODEL, DATAMODEL_LOH, DATAALL, YMIN, YMAX, R0, R1, R0b, R1b) {
    kpText(kp, chr=CHR, ymin=0, ymax=1, x=IDPOS, y=0.95, labels=ID, col="black", pos=3, r0=R0b, r1=R1b, cex=0.5, family="Arial")
    
    #COPY-NUMBER
    kpRect(kp, chr=CHR, x0=FROM, x1=TO, y0=YMIN+0.5, y1=YMAX-0.5, col="grey95", r0=R0, r1=R1, ymin=YMIN, ymax=YMAX, border=NA)
    kpRect(kp, data=DATAMODEL_LOH, chr=DATAMODEL_LOH$CONTIG, x0=DATAMODEL_LOH$START, x1=DATAMODEL_LOH$END, y0=-2.1, y1=2.1, col=rgb(142, 169, 219, max=255, alpha=70), r0=R0, r1=R1, ymin=YMIN, ymax=YMAX, border=NA)
    # kpRect(kp, data=DATAMODEL_LOH, chr=DATAMODEL_LOH$CONTIG, x0=DATAMODEL_LOH$START, x1=DATAMODEL_LOH$END, y0=-2.1, y1=2.1, col=DATAMODEL_LOH$COLOR, r0=R0, r1=R1, ymin=YMIN, ymax=YMAX, border=NA)
    kpText(kp, chr=CHR, ymin=YMIN, ymax=YMAX, x=FROM, y=YMAX, labels="Log2R", col="black", pos=4, r0=R0, r1=R1, cex=0.5, family="Arial")
    kpAbline(kp, h=c(-2, -1, 1, 2), col="grey70", ymin=YMIN, ymax=YMAX, r0=R0, r1=R1, lwd=0.5, lty=2)
    kpAbline(kp, h=c(0), col="grey10", ymin=YMIN, ymax=YMAX, r0=R0, r1=R1, lwd=0.9, lty=2)
    kpSegments(kp, chr=CHR, x0=FROM, x1=FROM, y0=YMIN, y1=YMAX, col="black", r0=R0, r1=R1, ymin=YMIN, ymax=YMAX, lwd=1.7)
    kpText(kp, chr=CHR, ymin=YMIN, ymax=YMAX, x=c(FROM, FROM, FROM, FROM, FROM), y=c(-2, -1, 0, 1, 2), labels=c(-2, -1, 0, 1, 2), col="black", pos=2, r0=R0, r1=R1, cex=0.5, family="Arial", clipping = FALSE)
    # kpPoints(kp, data=DATADENOIS, chr=DATADENOIS$CONTIG, x=DATADENOIS$START, y=DATADENOIS$L2R, cex=0.5, r0=R0, r1=R1, ymin=YMIN, ymax=YMAX, pch=20, col=DATADENOIS$COLORABN)
    # kpPoints(kp, data=DATADENOISABN, chr=DATADENOISABN$CONTIG, x=DATADENOISABN$START, y=DATADENOISABN$L2R, cex=0.5, r0=R0, r1=R1, ymin=YMIN, ymax=YMAX, pch=20, col=DATADENOISABN$COLORABN)
    kpSegments(kp, data=DATADENOIS, chr=DATADENOIS$CONTIG, x0=DATADENOIS$START, x1=DATADENOIS$END, y0=DATADENOIS$L2R, y1=DATADENOIS$L2R, r0=R0, r1=R1, ymin=YMIN, ymax=YMAX, col=DATADENOIS$COLORABN, lwd=7)
    if(nrow(DATADENOISABN %>% as.data.frame())>0) {kpSegments(kp, data=DATADENOISABN, chr=DATADENOISABN$CONTIG, x0=DATADENOISABN$START, x1=DATADENOISABN$END, y0=DATADENOISABN$L2R, y1=DATADENOISABN$L2R, r0=R0, r1=R1, ymin=YMIN, ymax=YMAX, col=DATADENOISABN$COLORABN, lwd=7)}
    kpSegments(kp, data=DATAMODEL, chr=DATAMODEL$CONTIG, x0=DATAMODEL$START, x1=DATAMODEL$END, y0=DATAMODEL$MEAN_L2R_EDIT, y1=DATAMODEL$MEAN_L2R_EDIT, col="black", r0=R0, r1=R1, ymin=YMIN, ymax=YMAX, lwd=1.7)
    kpPoints(kp, data=DTB_CENTROMERE, chr=DTB_CENTROMERE$CONTIG, x=DTB_CENTROMERE$POSITION, y=0, cex=2, r0=R0, r1=R1, ymin=YMIN, ymax=YMAX, pch=20, col="black")
    
    #MAF
    # kpRect(kp, data=DATAMODEL_LOH, chr=DATAMODEL_LOH$CONTIG, x0=DATAMODEL_LOH$START, x1=DATAMODEL_LOH$END, y0=-0.2, y1=1.2, col=rgb(142, 169, 219, max=255, alpha=70), r0=R0b, r1=R1b, ymin=0, ymax=1, border=NA)
    kpText(kp, chr=CHR, ymin=0, ymax=1, x=FROM, y=1.1, labels="AF", col="black", pos=4, r0=R0b, r1=R1b, cex=0.5, family="Arial")
    kpAbline(kp, h=c(0, 0.5, 1), col="grey70", ymin=0, ymax=1, r0=R0b, r1=R1b, lwd=0.5, lty=2)
    kpSegments(kp, chr=CHR, x0=FROM, x1=FROM, y0=0, y1=1, col="black", r0=R0b, r1=R1b, ymin=0, ymax=1, lwd=1.7)
    kpText(kp, chr=CHR, ymin=0, ymax=1, x=c(FROM, FROM, FROM), y=c(0, 0.5, 1), labels=c(0, 0.5, 1), col="black", pos=2, r0=R0b, r1=R1b, cex=0.5, family="Arial", clipping = FALSE)
    if(file.exists(GETCAPTURE)) {kpRect(kp, data=DTB_CAPTURE, chr=DTB_CAPTURE$CONTIG, x0=DTB_CAPTURE$STARText3, x1=DTB_CAPTURE$ENDext3, y0=0, y1=1, col=rgb(255, 235, 120, max=255, alpha=100), r0=R0b, r1=R1b, ymin=0, ymax=1, border=NA)}
    # kpPoints(kp, data=DATAALL, chr=DATAALL$CONTIG, x=DATAALL$POSITION, y=DATAALL$MAF1, cex=0.4, r0=R0b, r1=R1b, ymax=1, ymin=0, pch=20, col=DATAALL$COL1)
    # kpPoints(kp, data=DATAALL, chr=DATAALL$CONTIG, x=DATAALL$POSITION, y=DATAALL$MAF2, cex=0.4, r0=R0b, r1=R1b, ymax=1, ymin=0, pch=20, col=DATAALL$COL2)
    kpPoints(kp, data=DATAALL, chr=DATAALL$CONTIG, x=DATAALL$POSITION, y=DATAALL$MAF, cex=0.4, r0=R0b, r1=R1b, ymax=1, ymin=0, pch=20, col=DATAALL$COL)
    
  }
  
  #DGV function
  DGV_CHART = function(DATA, POS1, POS2, COLOR) {
    #DATA PANEL 2 - DGV
    kpRect(kp, data.panel=2, chr=DATA$CONTIG, x0=DATA$START, x1=DATA$END, y0=POS1, y1=POS2, r0=0.25, r1=0.35, col=COLOR, border=NA)
  }
  
  #gnomAD function
  GNOMADSV_CHART = function(DATA, POS1, POS2, COLOR) {
    #DATA PANEL 2 - GNOMADSV
    kpRect(kp, data.panel=2, chr=DATA$CONTIG, x0=DATA$START, x1=DATA$END, y0=POS1, y1=POS2, r0=0.35, r1=0.45, col=COLOR, border=NA)
  }
  
  #Prepare RefSeq
  #------------------------------------------------------------------------------------------------------------------------------------------------
  GENES_DTB_PROCESSED_REGION = DTB_GENES_DTB_PROCESSED %>% filter(CONTIG == CHR) %>% filter(TxEND >= FROM) %>% filter(TxSTART <= TO) %>%
    mutate(TxCenter = (TxSTART+TxEND)/2) %>%
    mutate(WAY = ifelse(STRAND=="+", TxSTART, TxEND)) %>%
    filter(!grepl("^LOC", GENE_ID))
  
  if(nrow(GENES_DTB_PROCESSED_REGION)>0) {
    GENES_DTB_PROCESSED_REGION = GENES_DTB_PROCESSED_REGION %>%
      mutate(N = (1:n()))
    
    #create bins, position for each gene
    GENES_DTB_PROCESSED_BINS = GENES_DTB_PROCESSED_REGION %>% 
      # mutate(GENESIZE = TxEND-TxSTART) %>%
      mutate(TxSTART = TxSTART - SIZE*0.004) %>%
      mutate(TxEND = TxEND + SIZE*0.004) %>%
      # mutate(TxSTART = TxSTART - SIZE*0.005) %>%
      # mutate(TxEND = TxEND + SIZE*0.005) %>%
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
    
    if(REGID == "MYC-region") {GENES_DTB_PROCESSED_BINS = GENES_DTB_PROCESSED_BINS %>% mutate(EXON_ID_Y = 1.07*EXON_ID_Y)}
    
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
  # plot.params$leftmargin = MARGIN #if you want to change default params
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
                      main=paste0(GROUP, "  ●  ", REGID, " ● ", CHR, ":", FROM, "-", TO, " ● ", SIZE, "bp"),
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
  
  # #DGV DATABASE - part 1 of data.panel 2
  # GAINCHART = DGV_CHART(DTB_DGV_GAIN, 0.11, 0.28, rgb(84, 130, 53, max = 255, alpha = 50))
  # LOSSCHART = DGV_CHART(DTB_DGV_LOSS, 0.29, 0.46, rgb(192, 0, 0, max = 255, alpha = 50))
  # COMPLEXCHART = DGV_CHART(DTB_DGV_COMPLEX, 0.47, 0.64, rgb(123, 35, 170, max = 255, alpha = 50))
  # INVERSIONCHART = DGV_CHART(DTB_DGV_INVERSION, 0.65, 0.82, rgb(0, 112, 192, max = 255, alpha = 50))
  # OTHERCHART = DGV_CHART(DTB_DGV_OTHER, 0.83, 1.0, rgb(0, 0, 0, max = 255, alpha = 50))
  
  #DGV DATABASE - part 1 of data.panel 2 with selfdenoising data - only gain/loss/clx and different values
  GAINCHART = DGV_CHART(DTB_DGV_GAIN, 0.01, 0.3, rgb(84, 130, 53, max = 255, alpha = 50))
  LOSSCHART = DGV_CHART(DTB_DGV_LOSS, 0.31, 0.6, rgb(192, 0, 0, max = 255, alpha = 50))
  COMPLEXCHART = DGV_CHART(DTB_DGV_COMPLEX, 0.61, 0.9, rgb(123, 35, 170, max = 255, alpha = 50))
  
  
  # #BUIOLD gnomAD SV DATABASE
  # GAINCHART = GNOMADSV_CHART(DTB_GNOMADSV_GAIN, 0.01, 0.3, rgb(84, 130, 53, max = 255, alpha = 50))
  # LOSSCHART = GNOMADSV_CHART(DTB_GNOMADSV_LOSS, 0.31, 0.6, rgb(192, 0, 0, max = 255, alpha = 50))
  # COMPLEXCHART = GNOMADSV_CHART(DTB_GNOMADSV_COMPLEX, 0.61, 0.9, rgb(123, 35, 170, max = 255, alpha = 50))
  
  #BUILD MAPPABILITY
  kpBars(kp, data.panel=2, chr=MAPPABILITY$CONTIG, x0=MAPPABILITY$START, x1=MAPPABILITY$END, y1=MAPPABILITY$VALUE, r0=0.44, r1=0.35, col="black")
  
  if(file.exists(CTRLCN)==T) {
    #Add ctrl CN data itself
    kpRect(kp, data.panel=2, chr=CHR, x0=FROM, x1=TO, y0=0, y1=1, col=NA, r0=0.12, r1=0.02, border=NA)
    kpRect(kp, data.panel=2, chr=CHR, x0=FROM, x1=TO, y0=0, y1=1, col=NA, r0=0.13, r1=0.23, border=NA)
    kpRect(kp, data.panel=2, chr=CTRL_CN$CONTIG, x0=CTRL_CN$STARText, x1=CTRL_CN$ENDext, y0=0, y1=CTRL_CN$FofCASES, r0=0.12, r1=0.02, col=rgb(255, 235, 120, max=255, alpha=200), border=NA)
    kpRect(kp, data.panel=2, chr=CTRL_CN$CONTIG, x0=CTRL_CN$STARText, x1=CTRL_CN$ENDext, y0=0, y1=CTRL_CN$FofCASES, r0=0.13, r1=0.23, col=rgb(255, 235, 120, max=255, alpha=200), border=NA)
    kpRect(kp, data.panel=2, chr=CTRL_CN$CONTIG, x0=CTRL_CN$START, x1=CTRL_CN$END, y0=0, y1=CTRL_CN$FofGAIN, r0=0.12, r1=0.02, col=rgb(84, 130, 53, max = 255, alpha = 255), border=NA)
    kpRect(kp, data.panel=2, chr=CTRL_CN$CONTIG, x0=CTRL_CN$START, x1=CTRL_CN$END, y0=0, y1=CTRL_CN$FofLOSS, r0=0.13, r1=0.23, col=rgb(192, 0, 0, max = 255, alpha = 255), border=NA)
  }
  
  
  #GENES_DTB_PROCESSED
  if(REGID != "IGH-region" & REGID != "IGK-region" & REGID != "IGL-region") {
  #transcription
  if(nrow(GENES_DTB_PROCESSED_REGION)>0) {kpArrows(kp, data.panel=2, chr=GENES_DTB_PROCESSED_REGION$CONTIG, x0=GENES_DTB_PROCESSED_REGION$TxSTART, x1=GENES_DTB_PROCESSED_REGION$TxEND, y0=GENES_DTB_PROCESSED_REGION$CENTER, y1=GENES_DTB_PROCESSED_REGION$CENTER, r0=0.5, r1=1, length=0, lwd=0.5*GENES_DTB_PROCESSED_LayerHeight)}
  if(nrow(GENES_DTB_PROCESSED_ARROWS)>0) {kpArrows(kp, data.panel=2, chr=GENES_DTB_PROCESSED_ARROWS$CONTIG, x0=GENES_DTB_PROCESSED_ARROWS$ARROW_START, x1=GENES_DTB_PROCESSED_ARROWS$ARROW_END, y0=GENES_DTB_PROCESSED_ARROWS$CENTER, y1=GENES_DTB_PROCESSED_ARROWS$CENTER, r0=0.5, r1=1, lwd=0.5*GENES_DTB_PROCESSED_LayerHeight, length=0.13*GENES_DTB_PROCESSED_LayerHeight, col="black")}
  #exons
  if(nrow(EXONS)>0) {kpRect(kp, data.panel=2, chr=EXONS$CONTIG, x0=EXONS$EXON_START, x1=EXONS$EXON_END, y0=EXONS$EXON_Y0, y1=EXONS$EXON_Y1, col="black", r0=0.5, r1=1, border=NA)}
  #coding
  if(nrow(EXONS_CODING)>0) {kpRect(kp, data.panel=2, chr=EXONS_CODING$CONTIG, x0=EXONS_CODING$EXON_CODING_START, x1=EXONS_CODING$EXON_CODING_END, y0=EXONS_CODING$EXON_CODING_Y0, y1=EXONS_CODING$EXON_CODING_Y1, col="black", r0=0.5, r1=1, border=NA)}
  #gene symbols
  if(nrow(GENES_DTB_PROCESSED_REGION)>0) {kpText(kp, data.panel=2, chr=GENES_DTB_PROCESSED_REGION$CONTIG, x=GENES_DTB_PROCESSED_REGION$TxCenter, y=GENES_DTB_PROCESSED_REGION$EXON_ID_Y, labels=GENES_DTB_PROCESSED_REGION$GENE_ID, pos=3, font=3, r0=0.5, r1=1, cex=GENES_DTB_PROCESSED_LayerHeight)}
  }
  
  #Ig genes annotation

  if(REGID == "IGH-region" | REGID == "IGK-region" | REGID == "IGL-region") {
    
    #a) IG parts - VDJC
    #Big enough
    DTB_IG_ANNOT_REG_temp1 = DTB_IG_ANNOT_REG %>%
      filter(ID != "IGHJ" & ID != "IGKC" & ID != "IGKJ")  %>%
      filter(CONTIG == CHR)
    if(nrow(DTB_IG_ANNOT_REG_temp1)>0) {
    kpArrows(kp, data.panel=2, chr=DTB_IG_ANNOT_REG_temp1$CONTIG, x0=DTB_IG_ANNOT_REG_temp1$START1, x1=DTB_IG_ANNOT_REG_temp1$END1, y0=0.1, y1=0.1, r0=0.8, r1=1, lwd=1.5, col=DTB_IG_ANNOT_REG_temp1$COLOR, length=0.05)
    kpArrows(kp, data.panel=2, chr=DTB_IG_ANNOT_REG_temp1$CONTIG, x0=DTB_IG_ANNOT_REG_temp1$END1, x1=DTB_IG_ANNOT_REG_temp1$START1, y0=0.1, y1=0.1, r0=0.8, r1=1, lwd=1.5, col=DTB_IG_ANNOT_REG_temp1$COLOR, length=0.05)
    kpText(kp, data.panel=2, chr=DTB_IG_ANNOT_REG_temp1$CONTIG, x=DTB_IG_ANNOT_REG_temp1$POSITION, labels=DTB_IG_ANNOT_REG_temp1$IDSHORT, y=0, r0=0.8, r1=1, pos=1, font=2, cex=0.8, col=DTB_IG_ANNOT_REG_temp1$COLOR)
    }
    #Small
    DTB_IG_ANNOT_REG_temp2 = DTB_IG_ANNOT_REG %>%
      filter(ID == "IGHJ") %>%
      filter(CONTIG == CHR)
    if(nrow(DTB_IG_ANNOT_REG_temp2)>0) {
    kpArrows(kp, data.panel=2, chr=DTB_IG_ANNOT_REG_temp2$CONTIG, x0=DTB_IG_ANNOT_REG_temp2$START1, x1=DTB_IG_ANNOT_REG_temp2$END1, y0=0.1, y1=0.1, r0=0.8, r1=1, lwd=1.5, col=DTB_IG_ANNOT_REG_temp2$COLOR, length=0)
    # kpArrows(kp, data.panel=2, chr=DTB_IG_ANNOT_REG_temp2$CONTIG, x0=DTB_IG_ANNOT_REG_temp2$POSITION, x1=DTB_IG_ANNOT_REG_temp2$START1, y0=0.1, y1=0.4, r0=0.8, r1=1, lwd=0.5, lty=1, col=DTB_IG_ANNOT_REG_temp2$COLOR, length=0)
    kpArrows(kp, data.panel=2, chr=DTB_IG_ANNOT_REG_temp2$CONTIG, x0=DTB_IG_ANNOT_REG_temp2$POSITION, x1=DTB_IG_ANNOT_REG_temp2$POSITION, y0=0.1, y1=0.17, r0=0.8, r1=1, lwd=0.5, lty=1, col=DTB_IG_ANNOT_REG_temp2$COLOR, length=0)
    kpText(kp, data.panel=2, chr=DTB_IG_ANNOT_REG_temp2$CONTIG, x=DTB_IG_ANNOT_REG_temp2$START1, labels=DTB_IG_ANNOT_REG_temp2$IDSHORT, y=0, r0=0.8, r1=1, pos=1, font=2, cex=0.8, col=DTB_IG_ANNOT_REG_temp2$COLOR)
    }
    DTB_IG_ANNOT_REG_temp2 = DTB_IG_ANNOT_REG %>%
      filter(ID == "IGKC") %>%
      filter(CONTIG == CHR)
    if(nrow(DTB_IG_ANNOT_REG_temp2)>0) {
      kpArrows(kp, data.panel=2, chr=DTB_IG_ANNOT_REG_temp2$CONTIG, x0=DTB_IG_ANNOT_REG_temp2$START1, x1=DTB_IG_ANNOT_REG_temp2$END1, y0=0.1, y1=0.1, r0=0.8, r1=1, lwd=1.5, col=DTB_IG_ANNOT_REG_temp2$COLOR, length=0)
      kpArrows(kp, data.panel=2, chr=DTB_IG_ANNOT_REG_temp2$CONTIG, x0=DTB_IG_ANNOT_REG_temp2$POSITION, x1=DTB_IG_ANNOT_REG_temp2$START1-(SIZE*0.003), y0=0.1, y1=0.17, r0=0.8, r1=1, lwd=0.5, lty=1, col=DTB_IG_ANNOT_REG_temp2$COLOR, length=0)
      kpText(kp, data.panel=2, chr=DTB_IG_ANNOT_REG_temp2$CONTIG, x=DTB_IG_ANNOT_REG_temp2$START1-(SIZE*0.003), labels=DTB_IG_ANNOT_REG_temp2$IDSHORT, y=0, r0=0.8, r1=1, pos=1, font=2, cex=0.8, col=DTB_IG_ANNOT_REG_temp2$COLOR)
    }
    DTB_IG_ANNOT_REG_temp2 = DTB_IG_ANNOT_REG %>%
      filter(ID == "IGKJ") %>%
      filter(CONTIG == CHR)
    if(nrow(DTB_IG_ANNOT_REG_temp2)>0) {
      kpArrows(kp, data.panel=2, chr=DTB_IG_ANNOT_REG_temp2$CONTIG, x0=DTB_IG_ANNOT_REG_temp2$START1, x1=DTB_IG_ANNOT_REG_temp2$END1, y0=0.1, y1=0.1, r0=0.8, r1=1, lwd=1.5, col=DTB_IG_ANNOT_REG_temp2$COLOR, length=0)
      kpArrows(kp, data.panel=2, chr=DTB_IG_ANNOT_REG_temp2$CONTIG, x0=DTB_IG_ANNOT_REG_temp2$POSITION, x1=DTB_IG_ANNOT_REG_temp2$START1+(SIZE*0.003), y0=0.1, y1=0.17, r0=0.8, r1=1, lwd=0.5, lty=1, col=DTB_IG_ANNOT_REG_temp2$COLOR, length=0)
      kpText(kp, data.panel=2, chr=DTB_IG_ANNOT_REG_temp2$CONTIG, x=DTB_IG_ANNOT_REG_temp2$START1+(SIZE*0.003), labels=DTB_IG_ANNOT_REG_temp2$IDSHORT, y=0, r0=0.8, r1=1, pos=1, font=2, cex=0.8, col=DTB_IG_ANNOT_REG_temp2$COLOR)   
      }
    
    
  #b) IG detail - IGHC - with description
  DTB_IG_ANNOT_temp1 = DTB_IG_ANNOT %>% filter(grepl("IGHM|IGHG|IGHA|IGHE|IGHD$", ID)) %>%
    filter(CONTIG == CHR)
  if(nrow(DTB_IG_ANNOT_temp1) > 0) {
  kpRect(kp, data.panel=2, chr=DTB_IG_ANNOT_temp1$CONTIG, x0=DTB_IG_ANNOT_temp1$START1, x1=DTB_IG_ANNOT_temp1$END1, y0=0.15, y1=0.75, r0=0.67, r1=0.75, col=DTB_IG_ANNOT_temp1$COLOR, border=NA)
  kpRect(kp, data.panel=2, chr=DTB_IG_ANNOT_temp1$CONTIG, x0=DTB_IG_ANNOT_temp1$START2, x1=DTB_IG_ANNOT_temp1$END2, y0=0, y1=1, r0=0.67, r1=0.75, col=DTB_IG_ANNOT_temp1$COLOR, border=NA)
  kpText(kp, data.panel=2, chr=DTB_IG_ANNOT_temp1$CONTIG, x=DTB_IG_ANNOT_temp1$POSITION, y=0.6, labels=DTB_IG_ANNOT_temp1$ID_SHORT, col=DTB_IG_ANNOT_temp1$COLOR, pos=1, r0=0.67, r1=0.75, font=3, cex=0.5)
  }

  #b) IG detail - IGHVDJ & IGK, IGL 
  DTB_IG_ANNOT_temp2 = DTB_IG_ANNOT %>% filter(!grepl("IGHM|IGHG|IGHA|IGHE|IGHD$", ID)) %>%
    filter(CONTIG == CHR)
  if(nrow(DTB_IG_ANNOT_temp2) > 0) {
  kpRect(kp, data.panel=2, chr=DTB_IG_ANNOT_temp2$CONTIG, x0=DTB_IG_ANNOT_temp2$START1-(SIZE/5000), x1=DTB_IG_ANNOT_temp2$END1+(SIZE/5000), y0=0.15, y1=0.75, r0=0.67, r1=0.75, col=DTB_IG_ANNOT_temp2$COLOR, border=NA)
  kpRect(kp, data.panel=2, chr=DTB_IG_ANNOT_temp2$CONTIG, x0=DTB_IG_ANNOT_temp2$START2-(SIZE/5000), x1=DTB_IG_ANNOT_temp2$END2+(SIZE/5000), y0=0, y1=1, r0=0.67, r1=0.75, col=DTB_IG_ANNOT_temp2$COLOR, border=NA)
  }
  }
  

  if(REGID == "IGH-region") {
    
    #Promoters
    AN_PROMOTER = DTB_BLUEPRINT %>% filter(grepl("IGH_PROMOTER", ID)) %>% mutate(PR_SIZE = END-START) %>% filter(PR_SIZE > 10000)
    kpRect(kp, data.panel=2, chr=AN_PROMOTER$CONTIG, x0=AN_PROMOTER$START, x1=AN_PROMOTER$END, y0=-0.07, y1=1.07, col="#E899B1",  r0=0.63, r1=0.65, border=NA)
    kpText(kp, data.panel=2, chr=CHR, ymin=0, ymax=1, x=AN_PROMOTER$END-(SIZE*0.005), y=0.6, labels="P", col="#E899B1", pos=4, r0=0.63, r1=0.65, cex=0.5, family="Arial", font=2, srt=0, clipping=F)
    
    #Switch regions (AID rich)
    kpPoints(kp, data.panel=2, chr=DTB_SREG$CONTIG, x=DTB_SREG$POSITION, y=0.5, cex=0.3, r0=0.63, r1=0.65, pch=4, col="black")
  
    #ENHANCERS
    CHART_ENH_FCE = function(SE_ID, YPOS, SE_NAME) {
      AN_SE = DTB_BLUEPRINT %>% filter(ID == SE_ID)
      kpRect(kp, data.panel=2, chr=AN_SE$CONTIG, x0=AN_SE$START, x1=AN_SE$END, y0=0, y1=1, r0=0.55, r1=0.6, col="#EF6F2E", border=NA)
      kpText(kp, data.panel=2, chr=AN_SE$CONTIG, x=(AN_SE$START+AN_SE$END)/2, y=YPOS, labels=SE_NAME, col="#EF6F2E", pos=3, r0=0.55, r1=0.6, font=2, cex=0.5)
    }
    Ealfa2 = CHART_ENH_FCE("IGH_ENHANCER_A2", 0.7, expression(bold(paste("E", alpha, "2"))))
    Ealfa1 = CHART_ENH_FCE("IGH_ENHANCER_A1", 0.7, expression(bold(paste("E", alpha, "1"))))
    # Emu1 = CHART_ENH_FCE("IGH_ENHANCER_M", -0.85, expression(bold(paste("E", mu))))
    # Edelta = CHART_ENH_FCE("IGH_ENHANCER_D", -0.77, expression(bold(paste("E", delta))))
    Emu = CHART_ENH_FCE("IGH_ENHANCER_D", 0.75, expression(bold(paste("E", mu))))
    Emu = CHART_ENH_FCE("IGH_ENHANCER_M", 0.75, expression(bold(paste("E", mu))))
    
    #Broad domain
    AN_BD = DTB_BLUEPRINT %>% filter(grepl("IGH_BD", ID))
    kpRect(kp, data.panel=2, chr="chr14", x0=AN_BD$START, x1=AN_BD$END, y0=0, y1=1, col="#E899B1", r0=0.55, r1=0.6, border=NA)
    kpText(kp, data.panel=2, chr="chr14", x=(AN_BD$START+AN_BD$END)/2, y=0.6, labels="BD", pos=3, font=2, r0=0.55, r1=0.6, cex=0.5, col="#E899B1")
    
    
    
    
    
    
    }

  
  
  
  
  
  invisible(dev.off())
}

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------

