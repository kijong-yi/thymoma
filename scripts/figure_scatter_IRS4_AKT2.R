
#load libraries
library(tidyverse)
library(ComplexHeatmap)
library(ggsci)
library("scales")
library(circlize)

# load meta data
meta_dt <- read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/02_metadata/Thymoma_summary_191106_1stSheet.txt')
MT_ids <- meta_dt$id[meta_dt$GTF2I_status2 == 'm']
WT_ids <- meta_dt$id[meta_dt$GTF2I_status2 == 'w']

# load biomart data
bm_dt <- read_tsv('~sypark/02_Reference/13_biomart/biomart_human_geneID_transcriptID_hgncsymbol_genetype_ensemblgenename_190522.txt')
bm_dt <- bm_dt %>% dplyr::rename(gene = `Gene name`, gene_type = `Gene type`) %>% dplyr::select(gene, gene_type) %>% unique()

# load expression data
exp_dt <- read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/01_expression/IRS4_corrected_v2/thym_137s_tpm_geneSymbol_tumorOnly_IRS4cor.tsv')
exp_dt_pcg <- left_join(exp_dt, bm_dt) %>% dplyr::filter(gene_type == 'protein_coding') %>% dplyr::select(-gene_type) # filter protein coding only
l10_exp_dt <- exp_dt %>% mutate_at(vars(-gene), funs(log10(.+0.01)))

#GSVA
l10_exp_mx <- l10_exp_dt %>% as.data.frame() %>% column_to_rownames('gene') %>% as.matrix()
require(GSVA)
c2cp_dt <- read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/04_GSVA/genesets/c2.cp.v6.2.symbols.gmt', col_names = F)
c2cp_list <- c2cp_dt %>% as.data.frame() %>% column_to_rownames('X1') %>% select(-X2) %>% t %>% as.data.frame() %>% as.list
c2cp_res <- gsva(l10_exp_mx, c2cp_list)
tmp_dt <- c2cp_res %>% as.data.frame() %>% rownames_to_column('pathway') %>% as.tibble()
tmp_dt <- tmp_dt %>% filter(pathway == 'REACTOME_PI3K_CASCADE') %>% gather(-pathway, key=id, value=REACTOME_PI3K_CASCADE) %>% dplyr::select(-pathway)

c5_list <- getGmt('~kjyi/ref/msigdb/c5.all.v6.2.symbols.gmt') %>% geneIds()
c2_list <- getGmt('~kjyi/ref/msigdb/c2.all.v6.2.symbols.gmt') %>% geneIds()
all_list <- append(c2_list,c5_list)

all_list <- all_list[c("GO_INSULIN_LIKE_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY","REACTOME_PI3K_CASCADE")]

all_gs_res <- gsva(l10_exp_mx,all_list)
all_gs_res[1:2,1:2]
meta_dt[1:2,1:2]

l10_exp_mx[1:2,1:2]
plot(l10_exp_mx["IRS4",WT_ids],all_gs_res[2,WT_ids])
plot(log10(as.data.frame(exp_dt[,-1])[exp_dt$gene == "IRS4",WT_ids]+0.01) %>% as.numeric,all_gs_res[2,WT_ids])
plot(log10(as.data.frame(exp_dt[,-1])[exp_dt$gene == "IRS4",WT_ids]+0.1) %>% as.numeric,all_gs_res[2,WT_ids])
plot(log10(meta_dt$corrected_IRS4_TPM[match(WT_ids,meta_dt$id)]+0.1),all_gs_res[2,WT_ids])


gsva_dt <- left_join(meta_dt, all_gs_res)
gsva_dt <- left_join(meta_dt, tmp_dt)
meta_dt$corrected_IRS4_TPM
ggplot(subset(gsva_dt, GTF2I_status2 == 'w'), aes(x=log10(corrected_IRS4_TPM+0.1), y=REACTOME_PI3K_CASCADE))+
  geom_point(aes(color = histologic_type), size=4, alpha=0.7)+
  scale_color_manual(values = histo_pal)+
  scale_x_continuous(breaks=log10(c(0,1,10,100,1000)+0.1), labels=c(0,1,10,100,1000))+
  scale_y_continuous(breaks = seq(-0.4,0.4,0.2), limits = c(-0.4, 0.4))+
  theme_bw()+ theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line())+
  ggtitle("GTF2I wild-type")+xlab("IRS4 expression(TPM)")



cairo_pdf("figures/IRS4_pathway_cor.1.pdf",height = 12/2.54,width=7/2.54,pointsize = 12*0.7)
par(mar=c(0,0,0,0),oma=c(4.5,4.5,0.7,0.7),mfrow=c(2,1))
tmp_df <- data.frame(IRS4=pmax(0,geneexpresionmatrix_WT["IRS4",]),
                     IGFsig=all_gs_res["GO_INSULIN_LIKE_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY",],
                     PI3K=all_gs_res["REACTOME_PI3K_CASCADE",],
                     histol=histo_pal[meta_dt$histologic_type[match(WT_samplenames,meta_dt$id)]])
model1 <- lm(IGFsig ~ IRS4, data=tmp_df)
model2 <- lm(PI3K ~ IRS4, data=tmp_df)
CI <- predict(model1, newdata=tmp_df, interval="prediction",level=0.9) %>% as.data.frame
CI2 <- predict(model2, newdata=tmp_df, interval="prediction",level=0.9) %>% as.data.frame
plot(1000,xlim=c(0,max(tmp_df$IRS4)),ylim=c(min(tmp_df$IGFsig),max(tmp_df$IGFsig)),
     xlab="",ylab="",xaxt="n",las=2)
mtexti("IGF receptor signaling score",2,0.45)
lines(x=tmp_df$IRS4, y=CI$fit, lwd=1.2,col="grey50")
polygon(c(sort(tmp_df$IRS4),rev(sort(tmp_df$IRS4))),c(CI$lwr[order(tmp_df$IRS4)],rev(CI$upr[order(tmp_df$IRS4)])),
        col="#00000010",lty = 0)
points(tmp_df$IRS4,tmp_df$IGFsig,bg=histo_pal[meta_dt$histologic_type[match(WT_samplenames,meta_dt$id)]],pch=21,cex=1.5)
text(0,max(tmp_df$IGFsig),
     bquote(R^2~"="~.(round(summary(model1)$r.squared,3))
     ),adj=c(0,1))
text(0,max(tmp_df$IGFsig),
     bquote(
       italic(p)~"="~.(format(anova(model1)$'Pr(>F)'[1],scientific = T,digits=3))
     ),adj=c(0,2.5))


plot(1000,xlim=c(0,max(tmp_df$IRS4)),ylim=c(min(tmp_df$PI3K),max(tmp_df$PI3K)),
     xlab=expression(log[10]~TPM~of~IRS4),ylab="PI3K cascade score",las=2,xaxt="n")
axis(1,at = c(0,1,2,log10(400)),labels = c(0,10,100,400))
mtexti("PI3K cascade score",2,0.45)
lines(x=tmp_df$IRS4, y=CI2$fit, lwd=1.2,col="grey50")
polygon(c(sort(tmp_df$IRS4),rev(sort(tmp_df$IRS4))),c(CI2$lwr[order(tmp_df$IRS4)],rev(CI2$upr[order(tmp_df$IRS4)])),
        col="#00000010",lty = 0)
points(tmp_df$IRS4,tmp_df$PI3K,bg=histo_pal[meta_dt$histologic_type[match(WT_samplenames,meta_dt$id)]],pch=21,cex=1.5)

text(0,max(tmp_df$PI3K),
     bquote(R^2~"="~.(round(summary(model2)$r.squared,3))
     ),adj=c(0,1))
text(0,max(tmp_df$PI3K),
     bquote(
       italic(p)~"="~.(format(anova(model2)$'Pr(>F)'[1],scientific = T,digits=3))
     ),adj=c(0,2.5))
# mtexti(expression(log[10]~TPM~of~IRS4),1,0.5)
mtexti(expression(TPM~of~IRS4),1,0.4)
legend("bottomright",legend=names(histo_pal[-c(3,7)]),pch=21,pt.bg=histo_pal[-c(3,7)],bty="n",pt.cex=1.5)
dev.off()

# ?predict.lm
# The prediction intervals are for a single observation at each case
# in ‘newdata’ (or by default, the data used for the fit) with error
# variance(s) ‘pred.var’.  This can be a multiple of ‘res.var’, the
# estimated value of sigma^2: the default is to assume that future
# observations have the same error variance as those used for
# fitting.  If ‘weights’ is supplied, the inverse of this is used as
# a scale factor.  For a weighted fit, if the prediction is for the
# original data frame, ‘weights’ defaults to the weights used for
# the model fit, with a warning since it might not be the intended
# result.  If the fit was weighted and ‘newdata’ is given, the
# default is to assume constant prediction variance, with a warning.



