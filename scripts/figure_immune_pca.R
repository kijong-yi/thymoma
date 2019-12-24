library(tidyverse)
library(GSVA)
source("~kjyi/Projects/thymus_single_cell/final2/scripts/biplotfun.R")

color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  # dev.new(width=1.75, height=5)
  plot(c(min,max),c(0,10), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title,yaxs='i')
  axis(1, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(y,0,y+1/scale,2, col=lut[i], border=NA)
  }
  rect(min,0,max,2)
  mtext(expression(log[10](CPM~"x"~100)),side = 1,line=2.5)
}

color.bar2 <- function(lut, min, max=-min,ylim=c(min,max), nticks=11, ticks=seq(min, max, len=nticks), title='', thickness=c(0,2)) {
  scale = (length(lut)-1)/(max-min)
  # dev.new(width=1.75, height=5)
  plot(c(0,10),c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title,xaxs='i',ylim=ylim)
  # axis(2, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(thickness[1],y,thickness[2],y+1/scale, col=lut[i], border=NA)
  }
  rect(thickness[1],min,thickness[2],max)
  # mtext(expression(log[10](CPM~"x"~100)),side = 1,line=2.5)
}




histo_pal = c("#631879FF", "#3B4992FF", "#5F559BFF", "#A20056FF", "#BB0021FF", "#EE0000FF", "#008B45FF", "#008280FF")
names(histo_pal) = c("A","AB","MN-T","B1","B2","B3","NE","TC")
gtf2i_pal = c(m="#3C5488FF",w="#E64B35FF",c="#00A087FF")

meta_dt <- read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/02_metadata/Thymoma_summary_191205_1stSheet.txt')
exp_dt <- read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/01_expression/IRS4_corrected_v2/thym_137s_tpm_geneSymbol_tumorOnly_IRS4cor.tsv')


markers <- list("B cell" = c("BLK","CD19","FCRL2","MS4A1","TNFRSF17","TCL1A","SPIB","PNOC"),
                "Cytotoxic" = c("PRF1","GZMA","GZMB","NKG7","GZMH","KLRK1","KLRB1","KLRD1","CTSW","GNLY"),
                "Dendritic cell" = c("CCL13","CD209","HSD11B1"),
                "Exhausted" = c("LAG3","TIGIT","HAVCR2","CTLA4"),
                "Macrophage" = c("CD68","CD163","MS4A4A", "CD84"),
                "Neutrophil" = c("FPR1","SIGLEC5","CSF3R","FCAR","FCGR3B"),
                "NK cell" = c("KIR2DL3","KIR3DL1","KIR3DL2","XCL1","XCL2","NCR1"),
                "proT" = c("ADGRG1","ATP6AP1L","CDK6","CEP70","ETS2","FXYD2","GUCY1A1","GUCY1B1","GXYLT2","HIVEP3","HOXA9","JCHAIN","MEST","NDN","RAB13","RGPD1","RIMS3","RRAS2","TLR7","TNFSF4"),
                "Doublepolar" = c("ACSBG1","AL357060.1","BCL6","C3orf52","CALN1","CD1C","CD1D","CD8A","CD8B","CIB2","CPLX1","CYP2U1","DNMBP","EFNB2","ELOVL4","HRK","LHFPL2","LYST","MCTP1","MIR646HG","RASD1","RIPK4","RMND5A","RORC","SH2D1A","SLAMF1","SLC7A3","TBC1D19"),
                "Singlepositive" = c("CLDN1","CTSL","DUSP4","EGR3","HTR2B","IER3","IFI44L","IRAK2","NR4A2","P2RY1","PDE4D","RSAD2","SERPINE2","SPRY2","TPRG1"),
                "naive" = c("ADTRP","AK5","ATP10A","C1orf162","EPPK1","GIMAP5","GIMAP8","IL6R","INPP4B","LDLRAP1","NOG","PASK","PCED1B","PLEKHA1","RAB30","RPS4Y1","SHISAL2A","TAGAP","UPP1","VSIG1"),
                "Late thymocytes" = c("CD1A","CD1B","DNTT"))
markers[["Thymopoiesis"]] <- c(markers[["proT"]],markers[["Doublepolar"]])
# markers[["Cytotoxic/exhausted"]] <- c(markers[["Exhausted"]],markers[["Cytotoxic"]])
# markers[["TNFa signaling"]] <- c("ABCA1","AC129492.1","ACKR3","AREG","ATF3","ATP2B1","B4GALT1","B4GALT5","BCL2A1","BCL3","BCL6","BHLHE40","BIRC2","BIRC3","BMP2","BTG1","BTG2","BTG3","CCL2","CCL20","CCL4","CCL5","CCN1","CCND1","CCNL1","CCRL2","CD44","CD69","CD80","CD83","CDKN1A","CEBPB","CEBPD","CFLAR","CLCF1","CSF1","CSF2","CXCL1","CXCL10","CXCL11","CXCL2","CXCL3","CXCL6","DDX58","DENND5A","DNAJB4","DRAM1","DUSP1","DUSP2","DUSP4","DUSP5","EDN1","EFNA1","EGR1","EGR2","EGR3","EHD1","EIF1","ETS2","F2RL1","F3","FJX1","FOS","FOSB","FOSL1","FOSL2","FUT4","G0S2","GADD45A","GADD45B","GCH1","GEM","GFPT2","GPR183","HBEGF","HES1","ICAM1","ICOSLG","ID2","IER2","IER3","IER5","IFIH1","IFIT2","IFNGR2","IL12B","IL15RA","IL18","IL1A","IL1B","IL23A","IL6","IL6ST","IL7R","INHBA","IRF1","IRS2","JAG1","JUN","JUNB","KDM6B","KLF10","KLF2","KLF4","KLF6","KLF9","KYNU","LAMB3","LDLR","LIF","LITAF","MAFF","MAP2K3","MAP3K8","MARCKS","MCL1","MSC","MXD1","MYC","NAMPT","NFAT5","NFE2L2","NFIL3","NFKB1","NFKB2","NFKBIA","NFKBIE","NINJ1","NR4A1","NR4A2","NR4A3","OLR1","PANX1","PDE4B","PDLIM5","PFKFB3","PHLDA1","PHLDA2","PLAU","PLAUR","PLEK","PLK2","PLPP3","PMEPA1","PNRC1","PPP1R15A","PTGER4","PTGS2","PTPRE","PTX3","RCAN1","REL","RELA","RELB","RHOB","RIPK2","RNF19B","SAT1","SDC4","SERPINB2","SERPINB8","SERPINE1","SGK1","SIK1","SLC16A6","SLC2A3","SLC2A6","SMAD3","SNN","SOCS3","SOD2","SPHK1","SPSB1","SQSTM1","STAT5A","TANK","TAP1","TGIF1","TIPARP","TLR2","TNC","TNF","TNFAIP2","TNFAIP3","TNFAIP6","TNFAIP8","TNFRSF9","TNFSF9","TNIP1","TNIP2","TRAF1","TRIB1","TRIP10","TSC22D1","TUBB2A","VEGFA","YRDC","ZBTB10","ZC3H12A","ZFP36")
# markers[["TNFa signaling"]] <- markers[["TNFa signaling"]][markers[["TNFa signaling"]] %in% rownames(l10_exp_dt2)]
# markers[["Neutrophil/macrophage"]] <- c(markers[["Neutrophil"]],markers[["Macrophage"]])



l10_exp_dt2 <- exp_dt %>% mutate_at(vars(-gene), list(~log10(.+0.01))) %>% column_to_rownames("gene")
exp_dt2 <- exp_dt %>% column_to_rownames("gene")
df22 <- gsva(as.matrix(l10_exp_dt2), markers[c("Thymopoiesis","Cytotoxic","Exhausted", "Neutrophil","NK cell", "Macrophage","Dendritic cell", "B cell")], method = "ssgsea")

write_rds(df22,"~kjyi/Projects/thymus_single_cell/final2/data/immune_ssgsea_all.Rds")

l10_exp_dt2 <- l10_exp_dt2[,!colnames(l10_exp_dt2) %in% c("TCGA-5U-AB0D", "SNU_09_C")]
exp_dt2 <- exp_dt2[,!colnames(exp_dt2) %in% c("TCGA-5U-AB0D", "SNU_09_C")]
# l10_exp_dt3 <- l10_exp_dt2[,meta_dt$id[meta_dt$final_cellularity <= 0.5]]
# df <- gsva(as.matrix(l10_exp_dt3), markers[c("Thymopoiesis","Cytotoxic","Exhausted", "Neutrophil","NK cell", "Macrophage","Dendritic cell", "B cell")], method = "ssgsea")
df2 <- gsva(as.matrix(l10_exp_dt2), markers[c("Thymopoiesis","Cytotoxic","Exhausted", "Neutrophil","NK cell", "Macrophage","Dendritic cell", "B cell")], method = "ssgsea")

write_rds(df2,"~kjyi/Projects/thymus_single_cell/final2/data/immune_ssgsea.Rds")

IGH_gene_ids  <- rownames(l10_exp_dt2)[grep("^IGH",rownames(l10_exp_dt2))]
IGH_sum <- colSums(exp_dt2[IGH_gene_ids,]) %>% {log10(.+0.01)}


pca2 <- prcomp(t(df2))
s2 <- summary(pca2)


# 912&height=496
# 496/72
cairo_pdf("~kjyi/Projects/thymus_single_cell/final2/figures/pca_immune_score.pdf",height = 10/2.54,width=19/2.54,pointsize = 12*0.7)

layout(matrix(c(1,1,2,3,
                1,1,4,5),byrow=T,ncol=4), widths=c(1,1,1,1),heights=c(1,1))
par(mar=c(0,0,0,2),oma=c(5,5,1,1))
biplot.prcomp2(pca2, bg=gtf2i_pal[meta_dt$GTF2I_status2[match(colnames(l10_exp_dt2),meta_dt$id)]],pt.cex=2.5,
               xlim=c(-0.35,0.25),ylim=c(-0.2,0.25),axis.text.cex = 1.5,
               xlab="",ylab="",
               pch=21)
mtext(paste("PCA 1 (", round(s2$importance[2]*100, 1), "%)", sep = ""),side=1,line=2.5)
mtext(paste("PCA 2 (", round(s2$importance[5]*100, 1), "%)", sep = ""),side=2,line=2.5)
abline(v=0,h=0,lty=2,col="grey40")

# 
par(mar=c(0,0,0,0))
biplot.prcomp2(pca2, bg=circlize::colorRamp2(seq(0,2.5,length.out=6),
                                             c("#F2F2F2","#BFBFBF","#8C8C8C","#595959","#242424","#BF0404")
)(l10_exp_dt2["CD1A",]),pt.cex=1.7,bty='n',xaxt='n',yaxt='n',xlab='',ylab='',textaxis=F,
xlim=c(-0.35,0.25),ylim=c(-0.2,0.25),
pch=21);mtext("CD1A",cex=1,font=2, side=3,line=-2,adj=c(0.03))

par(new=T,pty='m',mar=c(0,3,0,0))
color.bar2(colorRampPalette(c("#F2F2F2","#BFBFBF","#8C8C8C","#595959","#242424","#BF0404"))(100),0,2.5,ylim=c(0,2.5*2.5),thickness=c(0,0.5))
axis(side=2,at = log10(c(0.99,10,100,300)+0.01),labels = c(0,10,100,300), las=1)
par(mar=c(0,0,0,0))


biplot.prcomp2(pca2, bg=circlize::colorRamp2(seq(0,6,length.out=6),
                                             c("#F2F2F2","#BFBFBF","#8C8C8C","#595959","#242424","#BF0404")
)(exp_dt2["IFNG",]),pt.cex=1.7,bty='n',xaxt='n',yaxt='n',xlab='',ylab='',textaxis=F,
xlim=c(-0.35,0.25),ylim=c(-0.2,0.25),
pch=21);mtext("IFNG",cex=1,font=2, side=3,line=-2,adj=c(0.03))

par(new=T,pty='m',mar=c(0,3,0,0))
color.bar2(colorRampPalette(c("#F2F2F2","#BFBFBF","#8C8C8C","#595959","#242424","#BF0404"))(100),0,6,ylim=c(0,6*2.5),thickness=c(0,0.5))
axis(side=2,at = seq(0,6,length.out=4),labels = seq(0,6,length.out=4), las=1)
par(mar=c(0,0,0,0))




biplot.prcomp2(pca2, bg=circlize::colorRamp2(seq(-0.5,1.3,length.out=6),
                                             c("#F2F2F2","#BFBFBF","#8C8C8C","#595959","#242424","#BF0404")
)(l10_exp_dt2["CTLA4",]),pt.cex=1.7,bty='n',xaxt='n',yaxt='n',xlab='',ylab='',textaxis=F,
xlim=c(-0.35,0.25),ylim=c(-0.2,0.25),
pch=21);mtext("CTLA4",cex=1,font=2, side=3,line=-2,adj=c(0.03))
par(new=T,pty='m',mar=c(0,3,0,0))
color.bar2(colorRampPalette(c("#F2F2F2","#BFBFBF","#8C8C8C","#595959","#242424","#BF0404"))(100),0,1.3,ylim=c(0,1.3*2.5),thickness=c(0,0.5))
axis(side=2,at = log10(c(0.99,5,10,20)+0.01),labels = c(0,5,10,20), las=1)
par(mar=c(0,0,0,0))


biplot.prcomp2(pca2, bg=circlize::colorRamp2(seq(1,4,length.out=6),
                                             c("#F2F2F2","#BFBFBF","#8C8C8C","#595959","#242424","#BF0404")
)(IGH_sum),pt.cex=1.7,bty='n',xaxt='n',yaxt='n',xlab='',ylab='',textaxis=F,
xlim=c(-0.35,0.25),ylim=c(-0.2,0.25),
pch=21);mtext("IGH (sum)",cex=c(1),font=2, side=3,line=-2,adj=c(0.03))
par(new=T,pty='m',mar=c(0,3,0,0))
color.bar2(colorRampPalette(c("#F2F2F2","#BFBFBF","#8C8C8C","#595959","#242424","#BF0404"))(100),1,4,ylim=c(1,4*2.5),thickness=c(0,0.5))
axis(side=2,at = log10(c(10,100,1000,10000)+0.01),labels = c(10,expression(10^2),expression(10^3),expression(10^4)), las=1)
par(mar=c(0,0,0,0))
dev.off()




biplot.prcomp2(pca2, bg=circlize::colorRamp2(seq(1,4,length.out=6),
                                             c("#F2F2F2","#BFBFBF","#8C8C8C","#595959","#242424","#BF0404")
)(IGH_sum),pt.cex=2,bty='n',xaxt='n',yaxt='n',xlab='',ylab='',textaxis=F,
xlim=c(-0.35,0.25),ylim=c(-0.2,0.25),draw.box=F,
pch=21);mtext("IGH*",cex=c(1),font=2, side=3,line=-2,adj=c(0.03))

par(new=T,pty='m',mar=c(0,3,0,0))
color.bar2(colorRampPalette(c("#29394c","#3288bd","#66c2a5","#abdda4","#e6f598",
                             "#ffff33","#fdae61","#f76234","#e84053","#ee165d"))(100),0,2,ylim=c(0,4),thickness=c(0,0.5))
axis(side=2,at = 0:2,labels = c(0,10,100), las=1)
par(mar=c(0,0,0,0))










