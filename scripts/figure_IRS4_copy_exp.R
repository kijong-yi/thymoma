

meta_dt <- read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/02_metadata/Thymoma_summary_191204_1stSheet.txt')
meta_dt$IRS4_tCN[meta_dt$id=="TCGA-X7-A8M3"] = 4
cairo_pdf("figures/boxplot_xpcopy_irs4tpm.pdf",height = 7/2.54,width=6/2.54,pointsize = 12*0.7)
par(mfrow=c(1,1),oma=c(0,0,0,0),mar=c(5,4,1,1))
meta_dt[meta_dt$GTF2I_status2=="w",] %>% 
  boxplot(log10(corrected_IRS4_TPM+1)~pmin(5,IRS4_tCN),data=.,xlab="",ylab="",xaxt="n",yaxt="n",
          col=c("#595045","#A69886","#BF2626","#8C1B1B","#400101"),bty="n")
axis(side = 1,at = 1:5,labels = c(1:4,"≥ 5"),line=0)
axis(side = 1,at = c(2,4,5), labels = c(1,2,"≥ 3"),tick=F,line = 1)
mtext(text = "Female",side = 1,line = 1,at = 0)
mtext(text = "Male",side = 1,line = 2,at = 0)
mtext(text = "Copy number (Xp)",side = 1,line = 3)
mtext(text = "IRS4 expression (TPM)",side = 2,line=2.5)
axis(side = 2,at = log10(c(0,10,50,100,200,400)+1), labels = c(0,10,50,100,200,400),las=2)
dev.off()
