"------------------------------------------------------------------------------"
"                                                                              "
"                           Fig. 4h snu19c yilongplot                          "
"                                                                              "
"------------------------------------------------------------------------------"



"------------------------------------------------------------------------------"
"                                 load files                                   "
"------------------------------------------------------------------------------"


library(tidyverse)
filepaths <- c("~sypark/00_Project/01_thymoma/03_WGS/15_Timing/SNU_19_C.079_21.timing/SNU_19_C_wgs.snv.timing_input.cninfo.scF.sorted.kat.proba",
               "~sypark/00_Project/01_thymoma/03_WGS/15_Timing/SNU_19_C.079_21.timing/SNU_19_C_wgs_segments.txt",
               "~sypark/00_Project/01_thymoma/03_WGS/03_DELLY/SNU_19_C_WGS_delly.vcf.rough_fi.igv_fi.simple_edit",
               "~sypark/00_Project/01_thymoma/03_WGS/07_Sequenza/SNU_19_C_wgs.079_21.XY.sequenza/SNU_19_C_wgs_sequenza_extract.RData")
SNU_19_C.SNV <- read_tsv(filepaths[1],col_types = cols(`#CHROM`="c"))
SNU_19_C.SNV$chr <-paste0("chr",str_replace(SNU_19_C.SNV$`#CHROM`, "MT", "M"))
SNU_19_C.SNV$start <- SNU_19_C.SNV$POS
SNU_19_C.SNV$end <- SNU_19_C.SNV$POS+1
SNU_19_C.SNV$vaf <- SNU_19_C.SNV$var_readN/(SNU_19_C.SNV$var_readN+SNU_19_C.SNV$ref_readN)
SNU_19_C.SNV$type <- paste0(SNU_19_C.SNV$REF,">",SNU_19_C.SNV$ALT)
SNU_19_C.SNV <- dplyr::select(SNU_19_C.SNV,chr,start,end,vaf,type,mutCN,totCN)
SNU_19_C.CNV <- read_tsv(filepaths[2],col_types = cols(chromosome="c"))
SNU_19_C.CNV$chr <- paste0("chr",str_replace(SNU_19_C.CNV$chromosome,"MT","M"))
SNU_19_C.CNV$start <- SNU_19_C.CNV$start.pos
SNU_19_C.CNV$end <- SNU_19_C.CNV$end.pos
load(filepaths[4])
dt_19 <- `SNU_19_C_wgs_sequenza_extract`
rm(`SNU_19_C_wgs_sequenza_extract`)
Cellularity=0.79
purity = Cellularity
Ploidy=2.1
ploidy=Ploidy
tail(SNU_19_C.CNV)
SNU_19_C.CNV$CNt0 <- SNU_19_C.CNV$CNt
SNU_19_C.CNV$CNt <- (SNU_19_C.CNV$depth.ratio/dt_19$avg.depth.ratio - (1-Cellularity))*Ploidy/Cellularity
SNU_19_C.CNV$A[is.na(SNU_19_C.CNV$A)] = SNU_19_C.CNV$CNt0[is.na(SNU_19_C.CNV$A)]
SNU_19_C.CNV$B[is.na(SNU_19_C.CNV$B)] = 0
tail(SNU_19_C.CNV)
SNU_19_C.StV <- read_tsv(filepaths[3])
SNU_19_C.StV$chr1 = paste0("chr",str_replace(SNU_19_C.StV$chr1,"MT","M"))
SNU_19_C.StV$chr2 = paste0("chr",str_replace(SNU_19_C.StV$chr2,"MT","M"))
SNU_19_C.StV$ori = c("3to5"="3'-5'",
                     "3to3"="3'-3'",
                     "5to5"="5'-5'",
                     "5to3"="5'-3'")[SNU_19_C.StV$ori]

"------------------------------------------------------------------------------"
"                               plot function                                  "
"------------------------------------------------------------------------------"

cn_dt <- read_tsv('~sypark/00_Project/01_thymoma/03_WGS/10_Smoothened_CN/SNU_26_C_wgs.q20Q20.mpileup.100kbcov.026_330.absCN', comment = '#', col_names = c('chr','pos','tdp','ndp','abscn'), col_types = "cdddd")
sv_dt <- read_tsv('~sypark/00_Project/01_thymoma/03_WGS/03_DELLY/SNU_26_C_WGS_delly.vcf.rough_fi.igv_fi.simple_edit.man_edit')

cn_dt <- read_tsv('~sypark/00_Project/01_thymoma/03_WGS/10_Smoothened_CN/SNU_19_C_wgs.q20Q20.mpileup.100kbcov.079_210.absCN', comment = '#', col_names = c('chr','pos','tdp','ndp','abscn'), col_types = "cdddd")
sv_dt <- read_tsv('~sypark/00_Project/01_thymoma/03_WGS/03_DELLY/SNU_19_C_WGS_delly.vcf.rough_fi.igv_fi.simple_edit')


cairo_pdf("figures/Fig5h.yilong_plot.pdf",height = 102/25.4,width=115/25.4,pointsize = 12*0.7*0.8)
# par(mfcol=c(3,2))
layout(matrix(c(1,4,
                2,5,
                3,6),ncol=2,byrow=T),widths=c(3.3,2))
filter(cn_dt, chr=="1") %>% .$pos %>% range
filter(cn_dt, chr=="8") %>% .$pos %>% range

# draw_svsketch(cn_dt, sv_dt, 1, 0, 249200001)
# axis(1, at=c(0:100)*50000000, labels = F)
# 
# draw_svsketch(cn_dt, sv_dt, 8, 0, 146300001)
# axis(1, at=c(0:100)*50000000, labels = F)

par(mar=c(1,2,1,1),oma=c(5.5,3,3,0))

SNV=SNU_19_C.SNV
CNV=SNU_19_C.CNV
seqz=dt_19
StV=SNU_19_C.StV
snp_pal = c("C>A"="#15A0EC","G>T"="#15A0EC",# blue
            "C>G"="#0D0C1B","G>C"="#0D0C1B",#black
            "C>T"="#F23A29","G>A"="#F23A29",#red
            "T>A"="#A1A1A1","A>T"="#A1A1A1",#grey
            "T>C"="#5AB440","A>G"="#5AB440",#green
            "T>G"="#F2BBC5","A>C"="#F2BBC5")
snp_pal = c("coamp" = "#F23A29",
            "onecopy" = "0D0C1B")


cnt <- function(dr, pu, pl, adr, cnn) {
  (dr/adr*(pu*pl/2+1-pu) - (1-pu))/pu*cnn
}
CHR='chr1'
chr_range= c(0, 249200001)
cnv_ymax = 7 # 4


plot(-1,-1,xlim=chr_range, ylim = c(0,cnv_ymax), xlab="",ylab="",bty='n',xaxt='n')
axis(1, at=c(0:100)*50000000, labels = F)
mtext("Absolute CN",2,2.5)

seqz.window <- seqz$ratio[[str_replace(CHR,"chr","")]]
seqz.window <- seqz.window[seqz.window$N >= 1, ]

seqz.window$mean <- cnt(seqz.window$mean, purity, ploidy, seqz$avg.depth.ratio,2)
seqz.window$q0 <- cnt(seqz.window$q0, purity, ploidy, seqz$avg.depth.ratio,2)
seqz.window$q1 <- cnt(seqz.window$q1, purity, ploidy, seqz$avg.depth.ratio,2)

seqz.window$q1 <- pmin(seqz.window$q1,cnv_ymax)
seqz.window$mean <- pmin(seqz.window$mean,cnv_ymax)
seqz.window$q0 <- pmax(seqz.window$q0,0)

rect(xleft = seqz.window$start, ybottom = seqz.window$q0,
            xright = seqz.window$end, ytop = seqz.window$q1,
            col="grey90",border=NA)
abline(h = 0:7, col="grey")

for (i in 1:nrow(StV)){
  if(StV$chr1[i] == CHR){
    abline(v= StV$pos1[i],
           col = "#DD4DF0",
           lwd=0.7)
  }
  if(StV$chr2[i] == CHR){
    abline(v= StV$pos2[i],
           col = "#DD4DF0",
           lwd=0.7)
  }
}
segments(y0 = seqz.window$mean,y1 = seqz.window$mean, x0 = seqz.window$start,
         x1 = seqz.window$end, lty = 1, lwd = 1, col = "black")

# --------------
snp_pal_fun=circlize::colorRamp2(c(2,3),c("grey","red"))
plot(x=SNV[SNV$chr==CHR,]$start,
     y=SNV[SNV$chr==CHR,]$vaf,
     # col=snp_pal[SNV[SNV$chr==CHR,]$type],
     col=snp_pal_fun(SNV[SNV$chr==CHR,]$mutCN),
     ylim=c(0,1),bty='n',xaxt='n',
     pch=20,cex=1)
mtext("Variant-allele frequency",2,2.5)
axis(1, at=c(0:100)*50000000, labels = F)
# --------------
plot(-1,-1,xlim=chr_range, ylim = c(0,4), xlab="",ylab="",bty='n',xaxt='n')
mtext(paste0("Positions on chr", 1, " (Mbp)"),1,2.5)
mtext("Mutated CN",2,2.5)
axis(1, at=c(0:100)*50000000, labels = c(0:100)*50)
abline(h = 0:3, col="grey")
# major copy number
adj_w=0.1 # prevent overlapping of major and minor cn
CNV$A_ = ifelse(CNV$A == CNV$B, CNV$A+adj_w,CNV$A)
CNV$B_ = ifelse(CNV$A == CNV$B, CNV$B-adj_w,CNV$B)
CNV[CNV$chr==CHR & CNV$end-CNV$start > 1E6,] %>% 
  with(segments(x0=start,x1=end,y0=A_,y1=A_,
                col="#EC651A", # red>orange
                lwd=3,lend=1))
CNV[CNV$chr==CHR & CNV$end-CNV$start > 1E6,] %>% 
  with(segments(x0=start,x1=end,y0=B_,y1=B_,
                col="#0D33A6", # blue
                lwd=3,lend=1))
# snv points
points(x=SNV[SNV$chr==CHR,]$start,
       y=SNV[SNV$chr==CHR,]$mutCN,
       # col=snp_pal[SNV[SNV$chr==CHR,]$type],
       col=snp_pal_fun(SNV[SNV$chr==CHR,]$mutCN),
       pch=20,cex=0.5)


# ------------------------------------------------------------------------------

CHR='chr8'
chr_range= range(0, 146300001)
cnv_ymax = 4

plot(-1,-1,xlim=chr_range, ylim = c(0,cnv_ymax), xlab="",ylab="",bty='n',xaxt='n')
axis(1, at=c(0:100)*50000000, labels = F)


seqz.window <- seqz$ratio[[str_replace(CHR,"chr","")]]
seqz.window <- seqz.window[seqz.window$N >= 1, ]

seqz.window$mean <- cnt(seqz.window$mean, purity, ploidy, seqz$avg.depth.ratio,2)
seqz.window$q0 <- cnt(seqz.window$q0, purity, ploidy, seqz$avg.depth.ratio,2)
seqz.window$q1 <- cnt(seqz.window$q1, purity, ploidy, seqz$avg.depth.ratio,2)

seqz.window$q1 <- pmin(seqz.window$q1,cnv_ymax)
seqz.window$mean <- pmin(seqz.window$mean,cnv_ymax)
seqz.window$q0 <- pmax(seqz.window$q0,0)

rect(xleft = seqz.window$start, ybottom = seqz.window$q0,
     xright = seqz.window$end, ytop = seqz.window$q1,
     col="grey90",border=NA)
abline(h = 0:7, col="grey")

for (i in 1:nrow(StV)){
  if(StV$chr1[i] == CHR){
    abline(v= StV$pos1[i],
           col = "#DD4DF0",
           lwd=0.7)
  }
  if(StV$chr2[i] == CHR){
    abline(v= StV$pos2[i],
           col = "#DD4DF0",
           lwd=0.7)
  }
}
segments(y0 = seqz.window$mean,y1 = seqz.window$mean, x0 = seqz.window$start,
         x1 = seqz.window$end, lty = 1, lwd = 1, col = "black")
# -----------------------
snp_pal_fun=circlize::colorRamp2(c(1,2),c("grey","red"))
plot(x=SNV[SNV$chr==CHR,]$start,
     y=SNV[SNV$chr==CHR,]$vaf,
     # col=snp_pal[SNV[SNV$chr==CHR,]$type],
     col=snp_pal_fun(SNV[SNV$chr==CHR,]$mutCN),
     ylim=c(0,1),bty='n',xaxt='n',
     pch=20,cex=1)
axis(1, at=c(0:100)*50000000, labels = F)
# ---------------------------
plot(-1,-1,xlim=chr_range, ylim = c(0,4), xlab="",ylab="",bty='n',xaxt='n')
axis(1, at=c(0:100)*50000000, labels = c(0:100)*50)
mtext(paste0("Positions on chr", 8, " (Mbp)"),1,2.5)
abline(h = 0:3, col="grey")
# major copy number
adj_w=0.1 # prevent overlapping of major and minor cn
CNV$A_ = ifelse(CNV$A == CNV$B, CNV$A+adj_w,CNV$A)
CNV$B_ = ifelse(CNV$A == CNV$B, CNV$B-adj_w,CNV$B)
CNV[CNV$chr==CHR & CNV$end-CNV$start > 1E6,] %>% 
  with(segments(x0=start,x1=end,y0=A_,y1=A_,
                col="#EC651A", # red>orange
                lwd=3,lend=1))
CNV[CNV$chr==CHR & CNV$end-CNV$start > 1E6,] %>% 
  with(segments(x0=start,x1=end,y0=B_,y1=B_,
                col="#0D33A6", # blue
                lwd=3,lend=1))
# snv points
points(x=SNV[SNV$chr==CHR,]$start,
       y=SNV[SNV$chr==CHR,]$mutCN,
       col=snp_pal_fun(SNV[SNV$chr==CHR,]$mutCN),
       pch=20,cex=0.5)
par(mfrow=c(1,1),new=T,mar=rep(0,4),oma=rep(0,4))
plot.new()
legend("bottom",legend = c("early (pre-amplification) base substitution", "post-amplification or unamplified base substitution"),col=c("red","grey"),pch=16,xpd=T,horiz = F,bty='n',bg="#00000000")

dev.off()















