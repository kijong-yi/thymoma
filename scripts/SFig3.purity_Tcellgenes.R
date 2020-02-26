"------------------------------------------------------------------------------"
"                           Supplementary Fig.3                                "
"                       purity ~ T cell expression genes                       "
"                          GTEx plot (illustrator)                             "
"------------------------------------------------------------------------------"

gtf2i_pal = c(m="#3C5488FF",w="#E64B35FF",c="#00A087FF")

meta_dt=read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/02_metadata/Thymoma_summary_191227_1stSheet.txt')

exp_dt <- read_tsv(paste0("~sypark/00_Project/01_thymoma/10_Final_data/01_expression/",
                          'IRS4_corrected_v2/thym_137s_tpm_geneSymbol_tumorOnly_IRS4cor.tsv')) %>%
  as.data.frame %>% column_to_rownames("gene") %>% .[,meta_dt$id]

plot_fx <-function(gene, label=gene){
  y = log10(exp_dt[gene,]+0.01)
  x = 1-meta_dt$final_cellularity
  fit = lm(unlist(y)~x)
  CI <- predict(fit, newdata=data.frame(x=x), interval="prediction",level=0.9) %>% as.data.frame
  plot(x,y,col="#00000000", xlab="",ylab="",yaxt="n")
  lines(x=x, y=CI$fit, lwd=1.2,col="grey50")
  polygon(c(sort(x),rev(sort(x))),c(CI$lwr[order(x)],rev(CI$upr[order(x)])),
          col="#00000010",lty = 0)
  points(x,y, pch=21, bg = gtf2i_pal[meta_dt$GTF2I_status2],cex=2.2,lwd=0.7)
  axis(side = 2,at = log10(c(1,10,100,1000)+0.01),labels = c(1,10,100,1000))
  mtext(text = bquote("  "~r~"="~.(round(summary(fit)$r.squared^0.5,3))),
        line=-2,side = 3,adj=c(0))
  mtext(text = bquote("  "~italic(p)~"="~.(format(anova(fit)$'Pr(>F)'[1],scientific = T,digits=3))),
        line=-3.2,side = 3,adj=c(0))
  mtext(paste0(label," expression (TPM)"),2,2.0)
  mtext("1 - purity",1,2.0)
}



cairo_pdf("figures/sup.fig3a.purity.Tcellgenes.pdf",height = 4.5/2.54,width=18/2.54,pointsize = 12*0.7*0.8)
par(mfrow=c(1,4),mar=c(3,4,.5,.5))
# plot_fx("PTPRC")
plot_fx("CD3E")
plot_fx("CD4")
plot_fx("CD8A")
plot_fx("CD1A")
# plot_fx("DNTT")
dev.off()

