library(tidyverse)
meta_dt <- read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/02_metadata/Thymoma_summary_191204_1stSheet.txt')
cairo_pdf("figures/vaf_wes_rna_gtf2i.pdf",height = 16/2.54,width=10/2.54,pointsize = 12*0.7)
par(mfrow=c(2,1),oma=c(4,0,3,0),mar=c(0,0,0,1))

gr.pal = c("c"="#4CBCAB","m"="#7687AB","w"="#EE8071")
meta_dt$is.rescued = meta_dt$TCGA_paper_GTF2Imt == "w" & meta_dt$GTF2I_status2 == "m"
sum(meta_dt$is.rescued,na.rm=T)
meta_dt$final_cellularity
purity.pal=circlize::colorRamp2(0:10/10,RColorBrewer::brewer.pal(11,"Spectral"))

par(pty='s',xpd=F)
plot(100,100,
     xaxt='n',yaxt='n',xlim=c(0,log10(101)),ylim=c(0,log10(101)),
     xlab="",ylab="")
axis(side = 3,at = log10(c(0,0.01,0.05,0.1,0.5,1)*100+1), labels = c(0,0.01,0.05,0.1,0.5,1))
axis(side = 2,at = log10(c(0,0.01,0.05,0.1,0.5,1)*100+1), labels = c(0,0.01,0.05,0.1,0.5,1),las=2)
par(xpd=F)
abline(0,1,lty=2,col="grey")
# abline(h=0,v=0,lty=2,col="grey")
points(log10(meta_dt$GTF2I_WES_VAF+1),
       log10(meta_dt$GTF2I_RNAseq_VAF+1),
       pch = 21, bg = gr.pal[meta_dt$GTF2I_status2],
       col = ifelse(meta_dt$is.rescued, "red","black"),
       cex=ifelse(meta_dt$is.rescued, 1.5,1.2),
       lwd=ifelse(meta_dt$is.rescued, 1.5,1))

legend("bottomright",legend = c("Thymic carcinoma",expression("GTF2I"^mut),expression("GTF2I"^WT),"Rescued (n=16)"),pch=21,
       pt.bg=c(gr.pal,"white"),col=c("black","black","black","red"),pt.cex=c(1.2,1.2,1.2,1.5),pt.lwd=c(1,1,1,1.5))

# text(log10(meta_dt$GTF2I_WES_VAF+1),
#        log10(meta_dt$GTF2I_RNAseq_VAF+1),labels = meta_dt$id
#        )


mtext("Variant frequency of RNA-seq",2,line=2.6)

table(log10(meta_dt$GTF2I_WES_VAF+1)==0 & log10(meta_dt$GTF2I_RNAseq_VAF+1)==0)

# x=WESvaf, y=purity, col=type
plot(100,100,
     xaxt='n',yaxt='n',
     xlim=c(0,log10(101)),ylim=c(0.1,1),
     xlab="",ylab="")
mtext("Tumor cell fraction",2,line=2.5)
axis(side = 2, las=2)
axis(side = 1,at = log10(c(0,0.01,0.05,0.1,0.5,1)*100+1), labels = c(0,0.01,0.05,0.1,0.5,1))
# axis(side = 2,at = log10(c(0,0.01,0.05,0.1,0.5,1)*100+1), labels = rep("",6),las=2)
points(log10(meta_dt$GTF2I_WES_VAF+1),
       meta_dt$final_cellularity,
       pch = 21,
       col = ifelse(meta_dt$is.rescued, "red","black"),
       bg = gr.pal[meta_dt$GTF2I_status2],
       cex=ifelse(meta_dt$is.rescued, 1.5,1.2),
       lwd=ifelse(meta_dt$is.rescued, 1.5,1))

mtext("Variant frequency of WES",1,line=2)
dev.off()
# 
# # x=WESvaf, y=RNAvaf, col=purity
# plot(100,100,
#      xaxt='n',yaxt='n',
#      xlim=c(0,log10(101)),ylim=c(0,log10(101)),
#      xlab="",ylab="")
# axis(side = 1,at = log10(c(0,0.01,0.05,0.1,0.5,1)*100+1), labels = c(0,0.01,0.05,0.1,0.5,1))
# axis(side = 2,at = log10(c(0,0.01,0.05,0.1,0.5,1)*100+1), labels = rep("",6),las=2)
# par(xpd=F)
# abline(0,1,lty=2,col="grey")
# abline(h=0,v=0,lty=2,col="grey")
# points(log10(meta_dt$GTF2I_WES_VAF+1),
#        log10(meta_dt$GTF2I_RNAseq_VAF+1),
#        pch = 21,
#        bg = purity.pal(meta_dt$final_cellularity),
#        cex=ifelse(meta_dt$is.rescued, 1.5,1.2),
#        lwd=ifelse(meta_dt$is.rescued, 1.5,1))
# # mtext("Variant frequency of RNA-seq",2,line=2.6)
# mtext("Variant frequency of WES",1,line=2)
# legend_image <- as.raster(matrix(purity.pal(10:0/10), ncol=1))
# rasterImage(legend_image, 1.9, 0.0, 2,0.7)
# rect(1.9, 0.0, 2,0.7)
# segments(1.85,seq(0,0.7,length.out = 6),1.9)
# text(1.85,seq(0,0.7,length.out = 6),labels = seq(0,1,length.out=6),adj=1)
# text(2,0.7,"Purity",adj=c(1,-1),font=2)
# # text(x=1.5, y = seq(0,1,l=5), labels = seq(0,1,l=5))
# 
# 
# par(xpd=T)
# segments(0.25,0.15,0,0,col="grey")
# text(0.25,0.15,"(45 cases at 0,0)",adj=c(0,0.5),col="grey30")
# 
# na.omit(meta_dt$final_cellularity[meta_dt$is.rescued]) %>% summary
# mtext("Average purity of rescued cases = 15.3%",3,-1)
# 
# dev.off()





