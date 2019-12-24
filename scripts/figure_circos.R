read_tsv("inteneded_error")


# load ---------------------------------------------------------------------------------------
library(readr)
library(dplyr)
library(stringr)
library(circlize)

# snp_pal <- c("C>A"="#15A0EC",# blue
#              "C>G"="#0D0C1B",#black
#              "C>T"="#F23A29",#red
#              "T>A"="#A1A1A1",#grey
#              "T>C"="#5AB440",#green
#              "T>G"="#F2BBC5")#pink
# line_pal = c("5>5"="purple",
#              "3>3"="green",
#              "3>5"="red",
#              "5>3"="blue",
#              "interchromosomal"="purple")
# /home/users/sypark/00_Project/01_thymoma/10_Final_data/Rscripts/WGS_various_plots.R

filepaths <- c("~sypark/00_Project/01_thymoma/03_WGS/15_Timing/SNU_19_C.079_21.timing/SNU_19_C_wgs.snv.timing_input.cninfo.scF.sorted.kat.proba",
               "~sypark/00_Project/01_thymoma/03_WGS/15_Timing/SNU_26_C.026_31.timing/SNU_26_C_wgs.snv.timing_input.cninfo.scF.sorted.kat.proba",
               "~sypark/00_Project/01_thymoma/03_WGS/15_Timing/SNU_19_C.079_21.timing/SNU_19_C_wgs_segments.txt",
               "~sypark/00_Project/01_thymoma/03_WGS/15_Timing/SNU_26_C.026_31.timing/SNU_26_C_wgs.026_31.XY.sequenza3.withBP2_segments.txt",
               "~sypark/00_Project/01_thymoma/03_WGS/03_DELLY/SNU_19_C_WGS_delly.vcf.rough_fi.igv_fi.simple_edit",
               "~sypark/00_Project/01_thymoma/03_WGS/03_DELLY/SNU_26_C_WGS_delly.vcf.rough_fi.igv_fi.simple_edit.man_edit",
               "~sypark/00_Project/01_thymoma/03_WGS/07_Sequenza/SNU_19_C_wgs.079_21.XY.sequenza/SNU_19_C_wgs_sequenza_extract.RData",
               "~sypark/00_Project/01_thymoma/03_WGS/07_Sequenza/SNU_26_C_wgs.026_31.XY.sequenza3.withBP2/sequenza.extract")


SNU_19_C.SNV <- read_tsv(filepaths[1],col_types = cols(`#CHROM`="c"))
SNU_26_C.SNV <- read_tsv(filepaths[2],col_types = cols(`#CHROM`="c"))
SNU_19_C.SNV$chr <-paste0("chr",str_replace(SNU_19_C.SNV$`#CHROM`, "MT", "M"))
SNU_26_C.SNV$chr <-paste0("chr",str_replace(SNU_26_C.SNV$`#CHROM`, "MT", "M"))
SNU_19_C.SNV$start <- SNU_19_C.SNV$POS
SNU_26_C.SNV$start <- SNU_26_C.SNV$POS
SNU_19_C.SNV$end <- SNU_19_C.SNV$POS+1
SNU_26_C.SNV$end <- SNU_26_C.SNV$POS+1
SNU_19_C.SNV$vaf <- SNU_19_C.SNV$var_readN/SNU_19_C.SNV$ref_readN
SNU_26_C.SNV$vaf <- SNU_26_C.SNV$var_readN/SNU_26_C.SNV$ref_readN
SNU_19_C.SNV$type <- paste0(SNU_19_C.SNV$REF,">",SNU_19_C.SNV$ALT)
SNU_26_C.SNV$type <- paste0(SNU_26_C.SNV$REF,">",SNU_26_C.SNV$ALT)
SNU_19_C.SNV <- dplyr::select(SNU_19_C.SNV,chr,start,end,vaf,type,mutCN,totCN)
SNU_26_C.SNV <- dplyr::select(SNU_26_C.SNV,chr,start,end,vaf,type,mutCN,totCN)


SNU_19_C.CNV <- read_tsv(filepaths[3],col_types = cols(chromosome="c"))
SNU_26_C.CNV <- read_tsv(filepaths[4],col_types = cols(chromosome="c"))

SNU_19_C.CNV$chr <- paste0("chr",str_replace(SNU_19_C.CNV$chromosome,"MT","M"))
SNU_19_C.CNV$start <- SNU_19_C.CNV$start.pos
SNU_19_C.CNV$end <- SNU_19_C.CNV$end.pos
# SNU_19_C.CNV <- select(SNU_19_C.CNV,chr,start,end,CNt,A,B)

SNU_26_C.CNV$chr <- paste0("chr",str_replace(SNU_26_C.CNV$chromosome,"MT","M"))
SNU_26_C.CNV$start <- SNU_26_C.CNV$start.pos
SNU_26_C.CNV$end <- SNU_26_C.CNV$end.pos
# SNU_26_C.CNV <- select(SNU_26_C.CNV,chr,start,end,CNt,A,B)

load(filepaths[7])
dt_19 <- `SNU_19_C_wgs_sequenza_extract`
rm(`SNU_19_C_wgs_sequenza_extract`)

Cellularity=0.79
Ploidy=2.1

tail(SNU_19_C.CNV)
SNU_19_C.CNV$CNt0 <- SNU_19_C.CNV$CNt
SNU_19_C.CNV$CNt <- (SNU_19_C.CNV$depth.ratio/dt_19$avg.depth.ratio - (1-Cellularity))*Ploidy/Cellularity
SNU_19_C.CNV$A[is.na(SNU_19_C.CNV$A)] = SNU_19_C.CNV$CNt0[is.na(SNU_19_C.CNV$A)]
SNU_19_C.CNV$B[is.na(SNU_19_C.CNV$B)] = 0
# SNU_19_C.CNV$CNt[SNU_26_C.CNV$chromosome == "Y"] = SNU_19_C.CNV$CNt[SNU_26_C.CNV$chromosome == "Y"]/2
tail(SNU_19_C.CNV)
# dt_19$ratio$`4`[dt_19$ratio$`4`$N>20000,] %>%
# {(.$mean/(dt_19$avg.depth.tumor/dt_19$avg.depth.normal) - (1-Cellularity))*Ploidy/Cellularity -> xx;plot(xx)}
# dt_19$ratio$`4` %>%
# {(.$mean/(dt_19$avg.depth.tumor/dt_19$avg.depth.normal) - (1-Cellularity))*Ploidy/Cellularity -> xx;plot(xx)}
# (dt_19$ratio$`4`$mean/(dt_19$avg.depth.ratio) - (1-Cellularity))*Ploidy/Cellularity -> xx;plot(xx)
# SNU_19_C.CNV %>% filter(chromosome == "4") %>% .$CNt %>% plot


dt_26 <- read_rds(filepaths[8])
# rm(SNU_26_C_wgs.026_31.XY.sequenza3.withBP2_sequenza_extract)

Cellularity=0.26
Ploidy=3.1
SNU_26_C.CNV$CNt0 <- SNU_26_C.CNV$CNt
SNU_26_C.CNV$CNt <- (SNU_26_C.CNV$depth.ratio/dt_26$avg.depth.ratio - (1-Cellularity))*Ploidy/Cellularity
SNU_26_C.CNV$A[is.na(SNU_26_C.CNV$A)] = SNU_26_C.CNV$CNt0[is.na(SNU_26_C.CNV$A)]
SNU_26_C.CNV$B[is.na(SNU_26_C.CNV$B)] = 0
tail(SNU_26_C.CNV)
#

SNU_19_C.StV <- read_tsv(filepaths[5])
SNU_26_C.StV <- read_tsv(filepaths[6])

SNU_19_C.StV$chr1 = paste0("chr",str_replace(SNU_19_C.StV$chr1,"MT","M"))
SNU_19_C.StV$chr2 = paste0("chr",str_replace(SNU_19_C.StV$chr2,"MT","M"))
SNU_26_C.StV$chr1 = paste0("chr",str_replace(SNU_26_C.StV$chr1,"MT","M"))
SNU_26_C.StV$chr2 = paste0("chr",str_replace(SNU_26_C.StV$chr2,"MT","M"))

SNU_19_C.StV$ori = c("3to5"="3'-5'",
                     "3to3"="3'-3'",
                     "5to5"="5'-5'",
                     "5to3"="5'-3'")[SNU_19_C.StV$ori]
SNU_26_C.StV$ori = c("3to5"="3'-5'",
                     "3to3"="3'-3'",
                     "5to5"="5'-5'",
                     "5to3"="5'-3'")[SNU_26_C.StV$ori]
# circos function  ------------------------------------------------------------------------------

if(F){
  SNV=SNU_26_C.SNV;CNV=SNU_26_C.CNV;Stv=SNU_26_C.StV;seqz=dt_26;purity=0.12;ploidy=3.8;ymax1=7;ymax2=10;StVstyle=1;chrlist=paste0("chr",c(1:22,"X","Y"))
  snp_pal = c("C>A"="#15A0EC","G>T"="#15A0EC",# blue
              "C>G"="#0D0C1B","G>C"="#0D0C1B",#black
              "C>T"="#F23A29","G>A"="#F23A29",#red
              "T>A"="#A1A1A1","A>T"="#A1A1A1",#grey
              "T>C"="#5AB440","A>G"="#5AB440",#green
              "T>G"="#F2BBC5","A>C"="#F2BBC5")
}

circos <- function(SNV=SNU_19_C.SNV, # chr("chr1"...), start, end, type("C>A","C>T"....)
                   CNV=SNU_19_C.CNV, # chr("chr1"...), start, end, tCN(total copy number)
                   StV=SNU_19_C.StV, # chr, start, end, chr_, start_, end_, type("5'-5'","3'-3'","3'-5'","5'-3'","Interchromosomal")
                   seqz=dt_19,       # sequenza_extract.RData
                   purity=0.79,
                   ploidy=2.1,
                   ymax1=NULL,
                   ymax2=NULL,
                   snp_pal = c("C>A"="#15A0EC","G>T"="#15A0EC",# blue
                               "C>G"="#0D0C1B","G>C"="#0D0C1B",#black
                               "C>T"="#F23A29","G>A"="#F23A29",#red
                               "T>A"="#A1A1A1","A>T"="#A1A1A1",#grey
                               "T>C"="#5AB440","A>G"="#5AB440",#green
                               "T>G"="#F2BBC5","A>C"="#F2BBC5"),
                   line_pal = c("5'-5'"="#F23A29",
                                "3'-3'"="#F2B705",
                                "3'-5'"="#15A0EC",
                                "5'-3'"="#39A649",
                                "translocation"="#DD4DF0"),
                   chrlist=paste0("chr",c(1:22,"X","Y"))){
  circos.clear()
  circos.par("start.degree" = 90, # start from 12 oclock direction
             cell.padding   =c(0.00, 1.00, 0.00, 1.00), # exact ylim will be applied
             track.margin = c(0.015, 0.015), # increased gap between track
             gap.degree = c(rep(1,length(chrlist)-1),4))
             # cell.padding   =c(0.02, 1.00, 0.02, 1.00)) # original margin
  # Ideogram
  # circos.initializeWithIdeogram()
  
  circos.initializeWithIdeogram(plotType = NULL,chromosome.index = chrlist)
  circos.genomicIdeogram(track.height = 0.05)
  circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
    chr = get.cell.meta.data("sector.index") %>% str_replace("chr","")
    xcenter = get.cell.meta.data("xcenter")
    ycenter = get.cell.meta.data("ylim")[2]
    circos.text(xcenter, ycenter, chr,adj = c(0.5,-1.0)) # adj = position of  chromosome label
    # circos.genomicAxis(h = "top") # for ticks and axis (coordinates)
  })
  
  # snp tract
  # ymax=ceiling(max(SNV$mutCN))
  if(is.null(ymax1)){
    ymax = ceiling(max(c(CNV$A[CNV$end-CNV$start > 1E6])))
  } else {
    ymax=ymax1
  }
  circos.track(factors=chrlist, ylim=c(0,ymax),track.height =0.22,bg.border="#00000030") # initialize empty track
  for(CHR in chrlist) {
    # grey lines in background
    for(lvl in 1:(ymax)){ 
      circos.lines(get.cell.meta.data(name="xlim",sector.index = CHR), c(lvl,lvl),col="grey", sector.index=CHR)
    }
    
    # major copy number segments
    adj_w=0.1 # prevent overlapping of major and minor cn
    CNV$A_ = ifelse(CNV$A == CNV$B, CNV$A+adj_w,CNV$A)
    CNV$B_ = ifelse(CNV$A == CNV$B, CNV$B-adj_w,CNV$B)
    
    CNV[CNV$chr==CHR & CNV$end-CNV$start > 1E6,] %>% 
      with(circos.segments(x0=start,x1=end,y0=A_,y1=A_,sector.index=CHR,
                           col="#EC651A", # red>orange
                           lwd=3,lend=1))
    # minor copy number segments
    if(!CHR %in% c("chrX","chrY")){
      CNV[CNV$chr==CHR & CNV$end-CNV$start > 1E6,] %>% 
        with(circos.segments(x0=start,x1=end,y0=B_,y1=B_,sector.index=CHR,
                             col="#0D33A6", # blue
                             lwd=3,lend=1))
    }
    # black horizontal line at total copy number==2
    # circos.lines(get.cell.meta.data(name="xlim",sector.index = CHR), c(2,2),col="black", sector.index=CHR)
    
    # snv points
    if(CHR %in% c("chrX","chrY")){
      circos.points(x=SNV[SNV$chr==CHR,]$start,
                    y=SNV[SNV$chr==CHR,]$mutCN/2,
                    sector.index=CHR,
                    col=snp_pal[SNV[SNV$chr==CHR,]$type],
                    pch=20,cex=0.5)
    } else {
      circos.points(x=SNV[SNV$chr==CHR,]$start,
                    y=SNV[SNV$chr==CHR,]$mutCN,
                    sector.index=CHR,
                    col=snp_pal[SNV[SNV$chr==CHR,]$type],
                    pch=20,cex=0.5)
    }
  }
  circos.yaxis(side="left",at=(0:ymax),sector.index="chr1",labels.cex=0.6,tick.length = 0.1)

  # circos.track(factors=chrlist, ylim=c(0,cnv_ymax),track.height =0.0001,bg.border="#00000000")
  
  # copy number tract
  if(is.null(ymax2)){
    cnv_ymax = ceiling(max(CNV$CNt[CNV$end-CNV$start > 1E6]))+1
  } else {
    cnv_ymax=ymax2
  }
  circos.track(factors=chrlist,
               ylim=c(0,cnv_ymax),
               track.height =0.22,bg.border="#00000030",bg.col="#90efa810")
  # circos.yaxis(side="left",at=(0:cnv_ymax),sector.index="chr1",labels.cex=0.6,tick.length = 0.1)
  # circos.yaxis(side="right",at=(0:cnv_ymax),sector.index="chrX",labels.cex=0.5,lwd = 0)
  cnt <- function(dr, pu, pl, adr, cnn) {
    (dr/adr*(pu*pl/2+1-pu) - (1-pu))/pu*cnn
  }
  for(CHR in chrlist) {
    # grey lines in background
    if(T){
      seqz.window <- seqz$ratio[[str_replace(CHR,"chr","")]]
      seqz.window <- seqz.window[seqz.window$N >= 1, ]
      # seqz.window$mean <- (seqz.window$mean/(seqz$avg.depth.ratio) - (1-purity))*ploidy/purity
      # seqz.window$q0 <- (seqz.window$q0/(seqz$avg.depth.ratio) - (1-purity))*ploidy/purity
      # seqz.window$q1 <- (seqz.window$q1/(seqz$avg.depth.ratio) - (1-purity))*ploidy/purity
      seqz.window$mean <- cnt(seqz.window$mean, purity, ploidy, seqz$avg.depth.ratio,2)
      seqz.window$q0 <- cnt(seqz.window$q0, purity, ploidy, seqz$avg.depth.ratio,2)
      seqz.window$q1 <- cnt(seqz.window$q1, purity, ploidy, seqz$avg.depth.ratio,2)
      if(CHR %in% c("chrX","chrY")){
        seqz.window$mean <- seqz.window$mean/2
        seqz.window$q0 <- seqz.window$q0/2
        seqz.window$q1 <- seqz.window$q1/2
      }
      seqz.window$q1 <- pmin(seqz.window$q1,cnv_ymax)
      seqz.window$mean <- pmin(seqz.window$mean,cnv_ymax)
      seqz.window$q0 <- pmax(seqz.window$q0,0)
      
      circos.rect(xleft = seqz.window$start, ybottom = seqz.window$q0,
                  xright = seqz.window$end, ytop = seqz.window$q1,
                  col="grey90",border=NA,sector.index=CHR)
      for(lvl in 1:(cnv_ymax-1)){
        circos.lines(get.cell.meta.data(name="xlim",sector.index = CHR), c(lvl,lvl),col="grey70", sector.index=CHR)
      }
      # circos.segments(y0 = seqz.window$mean,y1 = seqz.window$mean, x0 = seqz.window$start,
      #          x1 = seqz.window$end, lty = 1, lwd = 1, col = "black",sector.index=CHR)
    }
  }
  # 
  
  if(T){
    for (i in 1:nrow(StV)){
      if(StV$chr1[i] == StV$chr2[i]){
        # h=ifelse(StV$ori[i] %in% c("5'-5'","3'-3'"), height1,height2)
        Col=line_pal[StV$ori[i]]
      } else {
        # h=NULL
        Col=line_pal["translocation"]
      }
      circos.segments(sector.index = StV$chr1[i],
                      x0 = StV$pos1[i], y0 = 0,
                      x1 = StV$pos1[i], y1 = cnv_ymax,
                      col = Col,
                      lwd=0.7)
      circos.segments(sector.index = StV$chr2[i],
                      x0 = StV$pos2[i], y0 = 0,
                      x1 = StV$pos2[i], y1 = cnv_ymax,
                      col = Col,
                      lwd=0.7)
    }
  }
  for(CHR in chrlist) {
    {
      seqz.window <- seqz$ratio[[str_replace(CHR,"chr","")]]
      seqz.window <- seqz.window[seqz.window$N >= 1, ]
      # seqz.window$q1 <- pmin(seqz.window$q1,cnv_ymax)
      # seqz.window$mean <- pmax(pmin(seqz.window$mean,cnv_ymax),0)
      # seqz.window$q0 <- pmax(seqz.window$q0,0)
      seqz.window$mean <- cnt(seqz.window$mean, purity, ploidy, seqz$avg.depth.ratio,2)
      seqz.window$q0 <- cnt(seqz.window$q0, purity, ploidy, seqz$avg.depth.ratio,2)
      seqz.window$q1 <- cnt(seqz.window$q1, purity, ploidy, seqz$avg.depth.ratio,2)
      if(CHR %in% c("chrX","chrY")){
        seqz.window$mean <- seqz.window$mean/2
        seqz.window$q0 <- seqz.window$q0/2
        seqz.window$q1 <- seqz.window$q1/2
      }
      seqz.window$q1 <- pmin(seqz.window$q1,cnv_ymax)
      seqz.window$mean <- pmin(seqz.window$mean,cnv_ymax)
      seqz.window$q0 <- pmax(seqz.window$q0,0)
      circos.segments(y0 = seqz.window$mean,y1 = seqz.window$mean, x0 = seqz.window$start,
                      x1 = seqz.window$end, lty = 1, lwd = 1, col = "black",sector.index=CHR)
    }
  }
  circos.yaxis(side="left",at=(0:cnv_ymax),sector.index="chr1",labels.cex=0.6,tick.length = 0.1)
  
  
  baseline1=get.cell.meta.data("cell.bottom.radius")
  baseline2=get.cell.meta.data("cell.bottom.radius")
  height1=0.1
  height2=0.1
  
  
  for (i in 1:nrow(StV)){
    if(StV$chr1[i] == StV$chr2[i]){
      h=ifelse(StV$ori[i] %in% c("5'-5'","3'-3'"), height1,height2)
      Col=line_pal[StV$ori[i]]
    } else {
      h=NULL
      Col=line_pal["translocation"]
    }
    circos.link(sector.index1 = StV$chr1[i], StV$pos1[i], 
                sector.index2 = StV$chr2[i], StV$pos2[i],
                col = Col,
                h=h,
                lwd=0.7,
                rou=ifelse(StV$type[i] %in% c("5'-5'","3'-3'"), baseline1,baseline2))
  }
}

# draw --------------------------------------------------------------------------------------------------------------

cairo_pdf("figures/circos.pdf", 
          height = 12/2.54, 
          width=24/2.54, 
          pointsize = 12*0.7)
snp_pal2 <- c("C>A"="#15A0EC",# blue
              "C>G"="#0D0C1B",#black
              "C>T"="#F23A29",#red
              "T>A"="#A1A1A1",#grey
              "T>C"="#5AB440",#green
              "T>G"="#F2BBC5")#pink
line_pal = c("5'-5'"="#F23A29",
             "3'-3'"="#F2B705",
             "3'-5'"="#15A0EC",
             "5'-3'"="#39A649",
             "translocation"="#DD4DF0")
par(mfrow=c(1,2))
circos(SNU_19_C.SNV,SNU_19_C.CNV,SNU_19_C.StV,dt_19,0.79,2.1,4,6)
mtext(side=3,"SNU_19_C",outer = F,line=-1,adj = 0.00,cex=1.1,font=2)
mtext(side=3,"purity: 79%",outer = F,line=-2.2,adj = 0.00,cex=1,font=1)
mtext(side=3,"ploidy: 2.1",outer = F,line=-3.2,adj = 0.00,cex=1,font=1)
lgd1=legend("bottomleft",lty=1, legend=names(line_pal),col=line_pal,bty='n')
circos(SNU_26_C.SNV,SNU_26_C.CNV,SNU_26_C.StV,dt_26,0.26,3.1,4,5)
mtext(side=3,"SNU_26_C",outer = F,line=-1,adj = 0.00,cex=1.1,font=2)
mtext(side=3,"purity: 26%",outer = F,line=-2.2,adj = 0.00,cex=1,font=1)
mtext(side=3,"ploidy: 3.1",outer = F,line=-3.2,adj = 0.00,cex=1,font=1)
lgd2=legend("bottomright",legend=c("major copy number","minor copy number"),
            bty='n',col=c("#EC651A","#0D33A6"),lty=1,lwd=3)
lgd3=legend(lgd2$rect$left+lgd2$rect$w,lgd2$rect$top-0.05,xjust = 1,yjust=0,
            pch=21, legend=names(snp_pal2),pt.bg=snp_pal2,bty='n')
dev.off()



load("~sypark/00_Project/01_thymoma/03_WGS/07_Sequenza/SNU_26_C_wgs.015_56.XY.sequenza3.withBP2/SNU_26_C_wgs.015_56.XY.sequenza3.withBP2_sequenza_cp_table.RData")
mycptable <- SNU_26_C_wgs.015_56.XY.sequenza3.withBP2_sequenza_cp_table
rm(SNU_26_C_wgs.015_56.XY.sequenza3.withBP2_sequenza_cp_table)
mycptable
par(mfrow=c(2,2))
image(pmax(log(mycptable$lpp),-10000))
image(pmax(log(mycptable$lpp),-500))
image(pmax(log(mycptable$lpp),-50))
image(pmax(log(mycptable$lpp),-100))
abline(h=0.5)

max(mycptable$lpp)
image(matrix(1:6,ncol=2))
Heatmap(mycptable$lpp,col = )
