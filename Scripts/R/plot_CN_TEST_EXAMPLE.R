

#ANNOTATION FILE
ANNOTATION = data.frame(CONTIG=c("chr21", "chr21"), START=c(15200001, 15700001), END=c(15300000, 15800000), TEXT=c("simulated loss", "simulated gain"), stringsAsFactors = F)
ANNOTATION = ANNOTATION %>%
  mutate(MID_POSITION = START+(END-START)/2)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------


#define regions of interest
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------
#Create data frame with regions of interest

REGIONSofINTEREST=data.frame(ID = c("my_region_loss", "my_region_gain"), 
                             CONTIG = c("chr21", "chr21"), 
                             START = c(15200001-100000, 15700001-100000), 
                             END = c(15300000+100000,15800000+100000),
                             stringsAsFactors = F)
# 
# REGIONSofINTEREST = data.frame()
# MY_REGION_LOSS = c("my_region_loss", "chr21", 15200001-100000, 15300000+100000)
# MY_REGION_GAIN = c("my_region_gain", "chr21", 15700001-100000, 15800000+100000)
# REGIONSofINTEREST = rbind(MY_REGION_LOSS, MY_REGION_GAIN)
# names(REGIONSofINTEREST) = c("ID", "CONTIG", "START", "END")
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
  SIZE = TO-FROM+1
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
  #DGV DATABASE - part 1 of data.panel 2 with selfdenoising data - only gain/loss/clx and different values
  GAINCHART = DGV_CHART(DTB_DGV_GAIN, 0.01, 0.3, rgb(84, 130, 53, max = 255, alpha = 50))
  LOSSCHART = DGV_CHART(DTB_DGV_LOSS, 0.31, 0.6, rgb(192, 0, 0, max = 255, alpha = 50))
  COMPLEXCHART = DGV_CHART(DTB_DGV_COMPLEX, 0.61, 0.9, rgb(123, 35, 170, max = 255, alpha = 50))
  
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
  
  
  # #GENES_DTB_PROCESSED
  # #transcription
  # if(nrow(GENES_DTB_PROCESSED_REGION)>0) {kpArrows(kp, data.panel=2, chr=GENES_DTB_PROCESSED_REGION$CONTIG, x0=GENES_DTB_PROCESSED_REGION$TxSTART, x1=GENES_DTB_PROCESSED_REGION$TxEND, y0=GENES_DTB_PROCESSED_REGION$CENTER, y1=GENES_DTB_PROCESSED_REGION$CENTER, r0=0.5, r1=1, length=0, lwd=0.5*GENES_DTB_PROCESSED_LayerHeight)}
  # if(nrow(GENES_DTB_PROCESSED_ARROWS)>0) {kpArrows(kp, data.panel=2, chr=GENES_DTB_PROCESSED_ARROWS$CONTIG, x0=GENES_DTB_PROCESSED_ARROWS$ARROW_START, x1=GENES_DTB_PROCESSED_ARROWS$ARROW_END, y0=GENES_DTB_PROCESSED_ARROWS$CENTER, y1=GENES_DTB_PROCESSED_ARROWS$CENTER, r0=0.5, r1=1, lwd=0.5*GENES_DTB_PROCESSED_LayerHeight, length=0.13*GENES_DTB_PROCESSED_LayerHeight, col="black")}
  # #exons
  # if(nrow(EXONS)>0) {kpRect(kp, data.panel=2, chr=EXONS$CONTIG, x0=EXONS$EXON_START, x1=EXONS$EXON_END, y0=EXONS$EXON_Y0, y1=EXONS$EXON_Y1, col="black", r0=0.5, r1=1, border=NA)}
  # #coding
  # if(nrow(EXONS_CODING)>0) {kpRect(kp, data.panel=2, chr=EXONS_CODING$CONTIG, x0=EXONS_CODING$EXON_CODING_START, x1=EXONS_CODING$EXON_CODING_END, y0=EXONS_CODING$EXON_CODING_Y0, y1=EXONS_CODING$EXON_CODING_Y1, col="black", r0=0.5, r1=1, border=NA)}
  # #gene symbols
  # if(nrow(GENES_DTB_PROCESSED_REGION)>0) {kpText(kp, data.panel=2, chr=GENES_DTB_PROCESSED_REGION$CONTIG, x=GENES_DTB_PROCESSED_REGION$TxCenter, y=GENES_DTB_PROCESSED_REGION$EXON_ID_Y, labels=GENES_DTB_PROCESSED_REGION$GENE_ID, pos=3, font=3, r0=0.5, r1=1, cex=GENES_DTB_PROCESSED_LayerHeight)}


  #Extra genes annotation instead of refseq
  kpRect(kp, data.panel=2, chr=ANNOTATION$CONTIG, x0=ANNOTATION$START, x1=ANNOTATION$END, y0=0, y1=1, col="black", r0=0.6, r1=0.7, border=NA)
  kpText(kp, data.panel=2, chr=ANNOTATION$CONTIG, x=ANNOTATION$MID_POSITION, y=-0.1, labels=ANNOTATION$TEXT, pos=1, font=2, r0=0.6, r1=0.7, cex=1, col="white")

  invisible(dev.off())
}

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------

