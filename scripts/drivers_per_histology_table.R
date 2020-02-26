
library(tidyverse)
meta_dt <- read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/02_metadata/Thymoma_summary_191204_1stSheet.txt')  

meta_dt$g = ifelse(
  !is.na(meta_dt$CYLD),"CYLD",
  ifelse(meta_dt$GTF2I_status2 == "m", "GTF2I",
         ifelse(meta_dt$corrected_IRS4_TPM > 30, "IRS4>30",
                ifelse(meta_dt$corrected_IRS4_TPM > 3, "IRS4>3",
                       "?"))))

table(meta_dt$corrected_IRS4_TPM > 3, meta_dt$GTF2I_status2 == "m")

table(meta_dt$corrected_IRS4_TPM > 5, meta_dt$GTF2I_status2 == "m")
table(meta_dt$g, meta_dt$GTF2I_status2)

max(meta_dt$corrected_IRS4_TPM[meta_dt$GTF2I_status2 == "m"])

#

table(meta_dt$g, meta_dt$histologic_type)[c(3,4,2,1),rev(c(6,1,2,3,4,5,8))] %>% {t(t(.)/colSums(.))} %>% barplot(horiz=T,las=1)
y <- table(meta_dt$g, meta_dt$histologic_type)[c(3,4,2,1),rev(c(6,1,2,3,4,5,8))]
x <- barplot(y,horiz=T,las=1)

meta_dt2 = dplyr::filter(meta_dt, GTF2I_status2 != "c")
table(meta_dt2$corrected_IRS4_TPM > 3,  meta_dt2$GTF2I_status2 == "m")
chisq.test(table(meta_dt2$corrected_IRS4_TPM > 3,  meta_dt2$GTF2I_status2 == "m"))
fisher.test(table(meta_dt2$corrected_IRS4_TPM > 3,  meta_dt2$GTF2I_status2 == "m"))


cairo_pdf("figures/driver_table.pdf",width=670/254,height = 620/254,pointsize = 12*0.7)
# s <- svglite::svgstring(width=670/254,height = 620/254,pointsize = 12*0.7)
par(mar=c(3,3.1,2,2.7),xpd=T)
y <- table(meta_dt$g, meta_dt$histologic_type)[c(3,5,4,2,1),rev(c(6,1,2,3,4,5,8))]

x <- barplot(y,horiz=T,las=1,
          main="",
          col=c("#395081","#f53c30","#f29381","#59B580","#AAAAAA"))
title("Demographics of driver mutations",cex=1.05)
legend("bottomright", c("GTF2I","IRS4 (TPM > 30)","IRS4 (30 > TPM > 3)","CYLD"),pt.bg=c("#395081","#f53c30","#f29381","#59B580","#AAAAAA"),pch=22,bty='n',pt.cex=2)
# })
  
text(colSums(y)[1],x[1],paste0(" ",round(y["CYLD",1]/colSums(y)[1] * 100,1),"%"),adj=0, cex=0.7)
text(colSums(y)[5],x[5],paste0(" ",round(y["GTF2I",5]/colSums(y)[5] * 100,1),"%\n ",
                               round(y["IRS4>30",5]/colSums(y)[5] * 100,1),"%"),adj=0, cex=0.7)
for(i in 2:4){
  text(colSums(y)[i],x[i],paste0(" ",round(y["GTF2I",i]/colSums(y)[i] * 100,1),"%\n ",
                                 round(y["IRS4>30",i]/colSums(y)[i] * 100,1),"%\n ",
                                 round(y["IRS4>3",i]/colSums(y)[i] * 100,1),"%"),adj=0, cex=0.7)
}
dev.off()
htmltools::browsable(htmltools::HTML(s()))


c("#56BEEC", "#050607", "#D53C32","#CBCACB","#ABCC72","#E7C9C6")




