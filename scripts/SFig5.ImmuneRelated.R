"------------------------------------------------------------------------------"
"                           Supplementary Fig.5                                "
"                     immune scores with histologic type                       "
"                     immune scores with purity                                "
"------------------------------------------------------------------------------"

library(tidyverse)

histo_pal = c("#631879FF", "#3B4992FF", "#5F559BFF", "#A20056FF", "#BB0021FF", "#EE0000FF", "#008B45FF", "#008280FF")
names(histo_pal) = c("A","AB","MN-T","B1","B2","B3","NE","TC")
gtf2i_pal = c(m="#3C5488FF",w="#E64B35FF",c="#00A087FF")
meta_dt <- read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/02_metadata/Thymoma_summary_191205_1stSheet.txt')

meta_dt <- meta_dt[!meta_dt$histologic_type %in% c("NE","MN-T"),]
df2 <- read_rds("~kjyi/Projects/thymus_single_cell/final2/data/immune_ssgsea.Rds")
df22 <-read_rds("~kjyi/Projects/thymus_single_cell/final2/data/immune_ssgsea_all.Rds")
dim(df22)
dim(df2) # two NE excluded
df1  <- df22[,meta_dt$id] 
dim(df1) # two NE excluded and one MNT


cairo_pdf("~kjyi/Projects/thymus_single_cell/final2/figures/SFig5.immune_histol.pdf",
          width=150/25.4,height = 100/25.4,pointsize = 12*0.7*0.9)
# s <- svglite::svgstring(width=(200)/25.4,height = 45/25.4,pointsize = 12*0.7*0.9)
# 
# layout(matrix(c(1,2,3,4,5),byrow=T,ncol=5), widths=c(1.2,1.6,1,1,1),heights=c(1,1))
par(mfrow=c(2,3))

"Boxplot: Thymopoiesis ~ histology"
par(mar=c(2.5,2.5,1,1),bty="L", pty='s',oma = c(1,2,0,0))
boxplot(df22["Thymopoiesis",match(meta_dt$id,colnames(df22))]~meta_dt$histologic_type,
        xlab="",ylab="",xaxt='n',yaxt='n',
        col=histo_pal[levels(factor(meta_dt$histologic_type))],cex.axis=0.8)
axis(2,cex.axis=0.9,padj=0.5)
mapply(axis, side = 1, at = axTicks(1), labels = levels(factor(meta_dt$histologic_type)),cex.axis=0.9)
mtext(text = "Thymopoiesis score", side = 2, line=2,cex=0.65)

"Boxplot: Cytotoxic ~ histology"
boxplot(df22["Cytotoxic",match(meta_dt$id,colnames(df22))]~meta_dt$histologic_type,
        xlab="",ylab="",xaxt='n',yaxt='n',
        col=histo_pal[levels(factor(meta_dt$histologic_type))])
axis(2,cex.axis=0.9,padj=0.5)
mapply(axis, side = 1, at = axTicks(1), labels = levels(factor(meta_dt$histologic_type)),cex.axis=0.9)
mtext(text = "Cytotoxic score", side = 2, line=2,cex=0.65)

"Boxplot: B cell ~ histology"
boxplot(df22["B cell",match(meta_dt$id,colnames(df22))]~meta_dt$histologic_type,
        xlab="",ylab="",xaxt='n',yaxt='n',
        col=histo_pal[levels(factor(meta_dt$histologic_type))],cex.axis=0.8)
axis(2,cex.axis=0.9,padj=0.5)
mapply(axis, side = 1, at = axTicks(1), labels = levels(factor(meta_dt$histologic_type)),cex.axis=0.9)
mtext(text = "B cell score", side = 2, line=2,cex=0.65)

"Scatterplot: cytotoxic ~ purity"
plot(1-meta_dt$final_cellularity, df22["Cytotoxic",match(meta_dt$id,colnames(df22))],
     xlab="",ylab="", pch=21,bg=histo_pal[meta_dt$histologic_type],cex=1.7,lwd=0.7,yaxt='n')
axis(2,cex.axis=0.9,padj=0.5)
mtext(text = "Cytotoxic score", side=2,line=2,cex=0.65)
mtext(text = "1-purity", side=1,line=2,cex=0.65)
"Scatterplot: B cell  ~ purity"
plot(1-meta_dt$final_cellularity, df22["B cell",match(meta_dt$id,colnames(df22))],
     xlab="",ylab="", pch=21,bg=histo_pal[meta_dt$histologic_type],cex=1.7,lwd=0.7,yaxt='n')
axis(2,cex.axis=0.9,padj=0.5)
mtext(text = "B cell score", side=2,line=2,cex=0.65)
mtext(text = "1-purity", side=1,line=2,cex=0.65)

# "Biplot"
# source("~kjyi/Projects/thymus_single_cell/final2/scripts/biplotfun.R")
# 
# pca2 <- prcomp(t(df1))
# biplot.prcomp2(pca2, bg=histo_pal[meta_dt$histologic_type],
#                pt.cex=2,
#                xlim=c(-0.37,0.25),ylim=c(-0.2,0.25),axis.text.cex = 1.5,arrow.len=0.05,
#                xlab="",ylab="",
#                pch=21,lwd=0.8,
#                textaxisbg="white",textaxisborder="black")

dev.off()
