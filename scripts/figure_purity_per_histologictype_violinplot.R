meta_dt <- read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/02_metadata/Thymoma_summary_191204_1stSheet.txt')
meta_dt$histologic_type2 = meta_dt$histologic_type
meta_dt$histologic_type2[meta_dt$histologic_type=="MN-T"] = "A"
meta_dt$histologic_type2[meta_dt$histologic_type=="NE"] = "C"
meta_dt$histologic_type2[meta_dt$histologic_type=="TC"] = "C"
meta_dt$histologic_type2 = factor(meta_dt$histologic_type2, levels=rev(c("A","AB","B1","B2","B3","C")))
meta_dt$`Histologic type` = factor(meta_dt$histologic_type2, levels=(c("A","AB","B1","B2","B3","C")))



cairo_pdf("figures/purity_per_histology.pdf",height = 10/2.54,width=13/2.54,pointsize = 12*0.7)
library(ggplot2)
histo_pal = c("#631879FF", "#3B4992FF", "#5F559BFF", "#A20056FF", "#BB0021FF", "#EE0000FF", "#008B45FF", "#008280FF")
names(histo_pal) = c("A","AB","MN-T","B1","B2","B3","NE","C")
library(ggridges)
ggplot(meta_dt, aes(x = final_cellularity, y = histologic_type2)) +
  geom_density_ridges(aes(fill = `Histologic type`)) +
  scale_fill_manual(values = histo_pal)+theme_classic2() + scale_x_continuous("Tumor cell fraction",breaks = c(0:5/5), limits = c(-0.1,1.1)) + scale_y_discrete("Histologic type")
dev.off()

# x <- density(meta_dt$final_cellularity)
# polygon(x$x,x$y)
# par(mfrow=c(2,1),mar=c(0,4,0,1),oma=c(4,0,1,1))
# plot(density(meta_dt$final_cellularity),xlab="",ylab="",xaxt='n',yaxt='n',xlim=c(0.1,1),main="")
# plot(100,xlab="",ylab="",ylim=c(0,8),xlim=c(0.1,1),main="")
# density(meta_dt$final_cellularity) %>% {polygon(.$x,.$y,col="#00000010")}
# density(meta_dt$final_cellularity[meta_dt$histologic_type2=="C"]) %>% {polygon(.$x,.$y+1,col="#00000010")}
# density(meta_dt$final_cellularity[meta_dt$histologic_type2=="B3"]) %>% {polygon(.$x,.$y+2,col="#00000010")}
# density(meta_dt$final_cellularity[meta_dt$histologic_type2=="B2"]) %>% {polygon(.$x,.$y+3,col="#00000010")}
# density(meta_dt$final_cellularity[meta_dt$histologic_type2=="B1"]) %>% {polygon(.$x,.$y+4,col="#00000010")}
# density(meta_dt$final_cellularity[meta_dt$histologic_type2=="AB"]) %>% {polygon(.$x,.$y+5,col="#00000010")}
# density(meta_dt$final_cellularity[meta_dt$histologic_type2=="A"]) %>% {polygon(.$x,.$y+6,col="#00000010")}
# 
# polygon(density(meta_dt$final_cellularitydensity(meta_dt$final_cellularity) %>% {polygon(.$x,.$y,col="#00000010")}),col="#01000010")
# 
# boxplot(meta_dt$final_cellularity~meta_dt$histologic_type2,horizontal=T,xlab="",ylab="",xaxt='n',yaxt='n',ylim=c(0.1,1))
# axis(1,at=c(0.1,0:5/5))
# axis(2,at=1:6,levels(meta_dt$histologic_type2),las=2)
# mtext("Tumor cell fraction",1,2)
# mtext("Histologic type",2,2.5)
