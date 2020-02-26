"------------------------------------------------------------------------------"
"                           Supplementary Fig.2                                "
"                               CYLD in TC                                     "
"                       Diploid proportion in 3 group                          "
"                         Chromothripsis yilong plot                           "
"------------------------------------------------------------------------------"




"------------------------------------------------------------------------------"
"                               CYLD in TC                                     "
"------------------------------------------------------------------------------"

source("~kjyi/Projects/Thymoma_single_cell/final2/scripts/figures_big_complex_heatmap.R")

TC_meta_dt <- meta_dt %>% filter(GTF2I_status2 == 'c',histologic_type!= "NE") %>% 
  arrange(CYLD,`16q`)

TC_onco_dt <- oncodt[c("CYLD","TP53","SMARCA4","NRAS"),TC_meta_dt$id] 
TC_cn_dt <- cndt %>% dplyr::select(chr_arm, TC_meta_dt$id) %>% as.data.frame() %>% column_to_rownames('chr_arm') %>% as.matrix()
TC_cn_dt <- TC_cn_dt[c("1q","16q"),,drop=F]
TC_occn_dt <- rbind(TC_onco_dt, TC_cn_dt)
class(TC_occn_dt) <- "numeric"
onco_pct <- c(paste0(round(rowSums(is.na(TC_onco_dt)==F)*100/ncol(TC_onco_dt),1),'%'), rep('',nrow(TC_cn_dt)))
TC_right= rowAnnotation(pct = anno_text(onco_pct, gp = gpar(fontsize =9)))
TC_body <- Heatmap(TC_occn_dt, cluster_rows = F, cluster_columns = F, col=occn_pal, show_heatmap_legend = F, na_col = "gray90", 
                   row_gap = unit(1,'mm'), row_title = NULL, row_names_side = "left", show_column_names = F,
                   column_gap = unit(1,'mm'), column_title = NULL)

cairo_pdf("figures/sup.fig2a.tc.oncogrid.pdf",height = 5/2.54,width=9/2.54,pointsize = 12*0.7)
draw(TC_body)
dev.off()


"------------------------------------------------------------------------------"
"                       Diploid proportion in 3 group                          "
"------------------------------------------------------------------------------"
meta_dt=read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/02_metadata/Thymoma_summary_191227_1stSheet.txt')
gtf2i_pal = c(m="#3C5488FF",w="#E64B35FF",c="#00A087FF")

cairo_pdf("figures/sup.fig2b.diploid_proportion.pdf",height = 6/2.54,width=5/2.54,pointsize = 12*0.7)
par(mar=c(5,4,1,1))
plt = boxplot(meta_dt$diploid_proportion~factor(meta_dt$GTF2I_status2, levels = c("m","w","c")),
       col=gtf2i_pal, ylab = "Diploid proportion of genome",xlab="",xaxt='n')
text(1:3, par("usr")[3],
     labels = c(expression(GTF2I^mut),expression(GTF2I^WT),"Thymic carcinoma"), 
     srt = 35, adj = c(1.1,1.1), xpd = TRUE, cex=1) 
dev.off()


"------------------------------------------------------------------------------"
"                         Chromothripsis yilong plot                           "
"------------------------------------------------------------------------------"

# input
#colnames(cninfo) = c("chr", "pos", "abscn")
#colnames(svinfo) = c("chr1", "pos1", "chr2", "pos2", "ori", "type")
#ori example: 3to3, 3to5, ...
#type example: DEL, DUP, TRA, INV
draw_svsketch <- function(cninfo, svinfo, chri, inipos, endpos){
  cninfo$id = paste(cninfo$chr, cninfo$pos, sep=":")
  ### Load and parse the SV information for driver chain (from dchain file)
  isv <- svinfo
  isv$group[isv$type=="DEL"] = 1
  isv$group[isv$type=="DUP"] = 2
  isv$group[isv$type=="TRA"] = 5
  isv$group[isv$type=="unknown"] = 6
  isv$group[isv$type=="INV" & isv$ori == "5to5"] = 3
  isv$group[isv$type=="INV" & isv$ori == "3to3"] = 4
  isv$svid = paste(isv$chr1, isv$pos1, isv$group, sep=":")
  isv$svid2 = paste(isv$chr2, isv$pos2, isv$group, sep=":")
  isv$chrpos1 = paste(isv$chr1, isv$pos1, sep=":")
  isv$chrpos2 = paste(isv$chr2, isv$pos2, sep=":")
  ## Plotting SVs for chr4
  isvintra <- subset(isv, isv$chr1 == chri & isv$chr2 == chri)
  isvinter1 <- subset(isv, isv$chr1 == chri & isv$chr2 != chri)
  isvinter2 <- subset(isv, isv$chr1 != chri & isv$chr2 == chri)
  isvsingle <- subset(isv, (isv$chr1 == chri & is.na(pos2) == T) | (isv$chr2 == chri & is.na(pos1) == T))
  par(mar=c(5.1,4.1,10.1,2.1))
  k <<- max(cninfo$abscn[cninfo$chr==chri], na.rm=T)
  plot(cninfo$abscn[cninfo$chr==chri & (cninfo$pos > inipos & cninfo$pos < endpos)] ~ cninfo$pos[cninfo$chr==chri & (cninfo$pos > inipos & cninfo$pos < endpos)], pch=20, col=rgb(0,0,0,.7), ylim=c(0,1.1*k), xlim = c(inipos, endpos), ylab="Absolute CN", xlab=paste0("Positions on chromosome ", chri), frame=F, xaxt = "n")
  axis(1, at = cninfo$pos[cninfo$chr==chri & (cninfo$pos > inipos & cninfo$pos < endpos)]%/%10000000*10000000)
  ## Breakpoints
  ### DELETIONS (blue)
  if (length(isvintra$svid[isvintra$group==1]) != 0){
    segments(isvintra$pos1[isvintra$group==1], 0, isvintra$pos1[isvintra$group==1], k*1.6, lty=1, xpd=TRUE, col=rgb(0,0,255/255,.4))
    segments(isvintra$pos2[isvintra$group==1], 0, isvintra$pos2[isvintra$group==1], k*1.6, lty=1, xpd=TRUE, col=rgb(0,0,255/255,.4))
  }
  ### TANDEM DUPLICATIONS (green)
  if (length(isvintra$svid[isvintra$group==2]) != 0){
    segments(isvintra$pos1[isvintra$group==2], 0, isvintra$pos1[isvintra$group==2], k*1.6, lty=1, xpd=TRUE, col=rgb(110/255,139/255,61/255,.4))
    segments(isvintra$pos2[isvintra$group==2], 0, isvintra$pos2[isvintra$group==2], k*1.6, lty=1, xpd=TRUE, col=rgb(110/255,139/255,61/255,.4))
  }
  ### HEAD-TO-HEAD INVERSIONS (red)
  if (length(isvintra$svid[isvintra$group==3]) != 0){
    segments(isvintra$pos1[isvintra$group==3], 0, isvintra$pos1[isvintra$group==3], k*1.3, lty=1, xpd=TRUE, col=rgb(255/255,0,0,.4))
    segments(isvintra$pos2[isvintra$group==3], 0, isvintra$pos2[isvintra$group==3], k*1.3, lty=1, xpd=TRUE, col=rgb(255/255,0,0,.4))
  }
  ### TAIL-TO-TAIL INVERSIONS (orange)
  if (length(isvintra$svid[isvintra$group==4]) != 0){
    segments(isvintra$pos1[isvintra$group==4], 0, isvintra$pos1[isvintra$group==4], k*1.3, lty=1, xpd=TRUE, col=rgb(255/255,128/255,0,.4))
    segments(isvintra$pos2[isvintra$group==4], 0, isvintra$pos2[isvintra$group==4], k*1.3, lty=1, xpd=TRUE, col=rgb(255/255,128/255,0,.4))
  }
  ### INTERCHROMOSOMAL TRANSLOCATIONS (purple)
  if (length(isvinter1$svid[isvinter1$group==5]) != 0){
    segments(isvinter1$pos1[isvinter1$group==5], 0, isvinter1$pos1[isvinter1$group==5], k*1.9, lty=1, xpd=TRUE, col=rgb(75/255,0/255,130/255,.4))
    for (i in 1:length(isvinter1$pos1)){
      text(x = isvinter1$pos1[isvinter1$group==5][i], y = k*1.1, labels = isvinter1$chrpos2[i], cex =0.5)
    }
  }
  if (length(isvinter2$svid[isvinter2$group==5]) != 0){
    segments(isvinter2$pos2[isvinter2$group==5], 0, isvinter2$pos2[isvinter2$group==5], k*1.9, lty=1, xpd=TRUE, col=rgb(75/255,0/255,130/255,.4))
    for (i in 1:length(isvinter2$pos2)){
      text(x = isvinter2$pos2[isvinter2$group==5][i], y = k*1.1, labels = isvinter2$chrpos1[i], cex =0.5)
    } 
  }
  ### Unknown mate(single) (gray)
  if (length(isvsingle$svid) != 0){
    segments(isvsingle$pos1, 0, isvsingle$pos1, k*1.9, lty=1, xpd=TRUE, col=rgb(96/255,96/255,96/255,.4))
    segments(isvsingle$pos2, 0, isvsingle$pos2, k*1.9, lty=1, xpd=TRUE, col=rgb(96/255,96/255,96/255,.4))
  }
  ## Horizontal Lines for Formatting
  segments(inipos,k*1.3,endpos,k*1.3, col=rgb(0,0,0,.4), lty=1, xpd=TRUE)
  segments(inipos,k*1.6,endpos,k*1.6, col=rgb(0,0,0,.4), lty=1, xpd=TRUE)
  segments(inipos,k*1.9,endpos,k*1.9, col=rgb(0,0,0,.4), lty=1, xpd=TRUE)
  
  theta=seq(0,pi, len=100)
  if (nrow(isvintra) !=0 ){
    isvintra$rad=abs(isvintra$pos2-isvintra$pos1)/2
    deldf <- isvintra[isvintra$group==1,]
    dupdf <- isvintra[isvintra$group==2,]
    hhidf <- isvintra[isvintra$group==3,]
    ttidf <- isvintra[isvintra$group==4,]
  }
  ## DELETIONS
  if (length(isvintra$svid[isvintra$group==1])!=0){
    for (j in seq(1,length(unique(deldf$svid)),1)){
      x = deldf$rad[j]*cos(theta)+((deldf$pos1[j]+deldf$pos2[j])/2)
      y = k*0.1*sin(theta)+k*1.6
      lines(x,y,col=rgb(0,0,255/255,.4), xpd=TRUE)
    }
  }
  
  ## DUPLICATIONS
  if (length(isvintra$svid[isvintra$group==2])!=0){
    for (j in seq(1,length(unique(dupdf$svid)),1)){
      x = dupdf$rad[j]*cos(theta)+((dupdf$pos1[j]+dupdf$pos2[j])/2)
      y = -k*0.1*sin(theta)+k*1.6
      lines(x,y,col=rgb(110/255,139/255,61/255,.4), xpd=TRUE)
    }
  }
  
  ## HEAD-TO-HEAD INVERSIONS
  if (length(isvintra$svid[isvintra$group==3])!=0){
    for (j in seq(1,length(unique(hhidf$svid)),1)){
      x = hhidf$rad[j]*cos(theta)+((hhidf$pos1[j]+hhidf$pos2[j])/2)
      y = k*0.1*sin(theta)+k*1.3
      lines(x,y,col=rgb(255/255,0,0,.4), xpd=TRUE)
    }
  }
  
  ## TAIL-TO-TAIL INVERSIONS
  if (length(isvintra$svid[isvintra$group==4])!=0){
    for (j in seq(1,length(unique(ttidf$svid)),1)){
      x = ttidf$rad[j]*cos(theta)+((ttidf$pos1[j]+ttidf$pos2[j])/2)
      y = -k*0.1*sin(theta)+k*1.3
      lines(x,y,col=rgb(255/255,128/255,0,.4), xpd=TRUE)
    }
  }
}

library(tidyverse)

cn_dt <- read_tsv('~sypark/00_Project/01_thymoma/03_WGS/10_Smoothened_CN/SNU_26_C_wgs.q20Q20.mpileup.100kbcov.026_330.absCN', comment = '#', col_names = c('chr','pos','tdp','ndp','abscn'), col_types = "cdddd")
sv_dt <- read_tsv('~sypark/00_Project/01_thymoma/03_WGS/03_DELLY/SNU_26_C_WGS_delly.vcf.rough_fi.igv_fi.simple_edit.man_edit')

cairo_pdf("figures/sup.fig2c.yilong_plot.pdf",height = 7/2.54,width=18/2.54,pointsize = 12*0.7)
par(mfrow=c(1,2))
draw_svsketch(cn_dt, sv_dt, 5, 50000000, 180000000)
draw_svsketch(cn_dt, sv_dt, 16, 70000000, 95000000)
dev.off()


