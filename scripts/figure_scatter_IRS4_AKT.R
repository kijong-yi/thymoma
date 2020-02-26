"------------------------------------------------------------------------------"
"                                 Fig.4a                                       "
"                                 IRS4~pathway                                 "
"------------------------------------------------------------------------------" 

library(tidyverse)
library(scales)
library(circlize)
library(fgsea)
require(GSVA)
# pals
histo_pal = c("#631879FF", "#3B4992FF", "#5F559BFF", "#A20056FF", "#BB0021FF", "#EE0000FF", "#008B45FF", "#008280FF")
names(histo_pal) = c("A","AB","MN-T","B1","B2","B3","NE","TC")
gtf2i_pal = c(m="#3C5488FF",w="#E64B35FF",c="#00A087FF")

# load exp data
meta_dt <- read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/02_metadata/Thymoma_summary_191227_1stSheet.txt')
bm_dt <- read_tsv(paste0('~sypark/02_Reference/13_biomart/',
                         'biomart_human_geneID_transcriptID_hgncsymbol',
                         '_genetype_ensemblgenename_190522.txt')) %>%
  dplyr::rename(gene = `Gene name`, gene_type = `Gene type`) %>% dplyr::select(gene, gene_type) %>% unique()
exp_dt <- read_tsv(paste0("~sypark/00_Project/01_thymoma/10_Final_data/01_expression/",
                          'IRS4_corrected_v2/',
                          'thym_137s_tpm_geneSymbol_tumorOnly_IRS4cor.tsv'))
exp_dt_pcg <- left_join(exp_dt, bm_dt, by="gene") %>%
  dplyr::filter(gene_type == 'protein_coding') %>% dplyr::select(-gene_type) # filter protein coding only
l10_exp_dt <- exp_dt %>% mutate_at(vars(-gene), list(~log10(.+0.01)))
l10_exp_dt <- left_join(l10_exp_dt, bm_dt, by="gene") %>%
  dplyr::filter(gene_type == 'protein_coding') %>% dplyr::select(-gene_type) # filter protein coding only
rm(bm_dt,exp_dt_pcg,exp_dt)
# load gene set
h_list <- GSEABase::getGmt('~kjyi/ref/msigdb/h.all.v6.2.symbols.gmt') %>% GSEABase::geneIds()
c5_list <- GSEABase::getGmt('~kjyi/ref/msigdb/c5.all.v6.2.symbols.gmt') %>% GSEABase::geneIds()
c2_list <- GSEABase::getGmt('~kjyi/ref/msigdb/c2.all.v6.2.symbols.gmt') %>% GSEABase::geneIds()
all_list <- append(append(h_list,c2_list),c5_list);rm(h_list,c2_list,c5_list)
all_list <- all_list[str_replace(names(all_list),"_.*","") %in% c("GO","HALLMARK","KEGG","REACTOME")]

# run gsva with log10(x+0.01)
l10_exp_dt_mat <- l10_exp_dt %>% as.data.frame %>% column_to_rownames("gene") %>% as.matrix
l10_exp_dt_mat <- l10_exp_dt_mat[,meta_dt$id]
all_gs_res <- l10_exp_dt_mat %>% gsva(all_list)
# all_gs_res <- all_gs_res[,meta_dt$id]

l10_exp_dt_mat[1:3,1:3]
all_gs_res[1:3,1:3]
colnames(l10_exp_dt_mat) == colnames(all_gs_res)
lm_p1 = foreach(i = 1:nrow(all_gs_res),.combine=rbind) %do% {
  res = cor.test(y = l10_exp_dt_mat["IRS4",meta_dt$GTF2I_status2=="w"],
                 x = all_gs_res[i,meta_dt$GTF2I_status2=="w"], method = "pearson")
  res2 = cor.test(y = l10_exp_dt_mat["IRS4",meta_dt$GTF2I_status2=="w"],
                  x = all_gs_res[i,meta_dt$GTF2I_status2=="w"], method = "spearman")
  c(res$estimate,res$p.value,res2$estimate,res2$p.value)
}
rownames(lm_p1) = rownames(all_gs_res)
colnames(lm_p1) = c("Pearson coeficient", "Pearson p-value", "Spearman coeficient","Spearman p-value")
View(lm_p1)

lm_p2 = foreach(i = 1:nrow(all_gs_res),.combine=rbind) %do% {
  summary(lm(all_gs_res[i,meta_dt$GTF2I_status2=="w"] ~ 
               meta_dt$final_cellularity[meta_dt$GTF2I_status2=="w"] + 
               l10_exp_dt_mat["IRS4",meta_dt$GTF2I_status2=="w"]))$coefficients -> x
  c(x[2:3,4],x[2:3,1])
}
rownames(lm_p2) = rownames(all_gs_res)
colnames(lm_p2) = c("lm purity p-value", "lm IRS4 p-value", "lm purity coeficient", "lm IRS4 coeficient")
lm_bind = cbind(lm_p1,lm_p2)
if(F) View(lm_bind)
library(rJava,lib.loc="/home/users/kjyi/R/x86_64-redhat-linux-gnu-library/3.6/")
library(xlsx,lib.loc="/home/users/kjyi/R/x86_64-redhat-linux-gnu-library/3.6/")
lm_bind %>% xlsx::write.xlsx(file = "tables/S.Table3.xlsx")



# plot -------------------------------------------------------------------------
all_gs_res <- all_gs_res[,WT_samplenames]
geneexpresionmatrix_WT <- l10_exp_dt_mat[,WT_samplenames]
cor(t(all_gs_res),geneexpresionmatrix_WT["IRS4",]) %>% sort(decreasing=T)

cor.test(y = geneexpresionmatrix_WT["IRS4",],
               x = all_gs_res["GO_INTRACELLULAR_LIPID_TRANSPORT",WT_samplenames], method = "pearson")



tmp_df <- data.frame(IRS4=pmax(0,geneexpresionmatrix_WT["IRS4",]),
                     IRS4_orig=geneexpresionmatrix_WT["IRS4",],
                     IGFsig=all_gs_res["GO_INSULIN_LIKE_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY",],
                     lipid=all_gs_res["GO_INTRACELLULAR_LIPID_TRANSPORT",],
                     PI3Kcascade=all_gs_res["REACTOME_PI3K_CASCADE",],
                     mito=all_gs_res["GO_MITOCHONDRION_DISTRIBUTION",],
                     histol=histo_pal[meta_dt$histologic_type[match(WT_samplenames,meta_dt$id)]],
                     purity=meta_dt$final_cellularity[match(WT_samplenames,meta_dt$id)])
summary(lm(IGFsig ~ IRS4_orig+purity, data=tmp_df))
summary(lm(IGFsig ~ purity, data=tmp_df))

model1 <- lm(IGFsig ~ IRS4_orig, data=tmp_df)
model2 <- lm(lipid ~ IRS4_orig, data=tmp_df)
model3 <- lm(PI3Kcascade ~ IRS4_orig, data=tmp_df)
model4 <- lm(mito ~ IRS4_orig, data=tmp_df)

cor(tmp_df$IRS4_orig,tmp_df$IGFsig)
cor(tmp_df$IRS4_orig,tmp_df$lipid)
cor(tmp_df$IRS4_orig,tmp_df$PI3Kcascade)
cor(tmp_df$IRS4_orig,tmp_df$mito)
CI <- predict(model1, newdata=tmp_df, interval="prediction",level=0.9) %>% as.data.frame
CI2 <- predict(model2, newdata=tmp_df, interval="prediction",level=0.9) %>% as.data.frame
CI3 <- predict(model3, newdata=tmp_df, interval="prediction",level=0.9) %>% as.data.frame
CI4 <- predict(model4, newdata=tmp_df, interval="prediction",level=0.9) %>% as.data.frame


# svglite::htmlSVG(width=1150/254,height=1150/254,pointsize = 12*0.7, {
cairo_pdf("figures/IRS4_pathway_cor.1.pdf",width=1150/254,height=1050/254,pointsize = 12*0.7)

par(mar=c(4.5,4.5,0,0),oma=c(0,0,1,1),mfrow=c(2,2))



# left bottom  -> left top
plot(1000,xlim=c(0,max(tmp_df$IRS4)),ylim=c(min(tmp_df$PI3Kcascade),max(tmp_df$PI3Kcascade)),
     ylab="",las=2,xaxt="n",xlab="")
axis(1,at = c(0,1,2,log10(400)),labels = c(0,10,100,400))
mtext("PI3K cascade score",2,2.45)
lines(x=tmp_df$IRS4_orig, y=CI3$fit, lwd=1.2,col="grey50")
polygon(c(sort(tmp_df$IRS4_orig),rev(sort(tmp_df$IRS4_orig))),c(CI3$lwr[order(tmp_df$IRS4_orig)],rev(CI3$upr[order(tmp_df$IRS4_orig)])),
        col="#00000010",lty = 0)
points(tmp_df$IRS4,tmp_df$PI3Kcascade,bg=histo_pal[meta_dt$histologic_type[match(WT_samplenames,meta_dt$id)]],pch=21,cex=1.5)
text(0,max(tmp_df$PI3Kcascade),
     bquote(r~"="~.(round(summary(model3)$r.squared^0.5,3))
     ),adj=c(0,1))
text(0,max(tmp_df$PI3Kcascade),
     bquote(
       italic(p)~"="~.(format(anova(model3)$'Pr(>F)'[1],scientific = T,digits=3))
     ),adj=c(0,2.5))
mtext("IRS4 expression (TPM)",1,2.0)
# legend("bottomright",legend=names(histo_pal[-c(1,3,7,8)]),pch=21,pt.bg=histo_pal[-c(1,3,7,8)],bty="n",pt.cex=1.5)


# left top -> right top
plot(1000,xlim=c(0,max(tmp_df$IRS4)),ylim=c(min(tmp_df$IGFsig),max(tmp_df$IGFsig)),
     xlab="",ylab="",xaxt="n",las=2)
mtext("IGF receptor signaling score",2,2.45)
lines(x=tmp_df$IRS4_orig, y=CI$fit, lwd=1.2,col="grey50")
polygon(c(sort(tmp_df$IRS4_orig),rev(sort(tmp_df$IRS4_orig))),c(CI$lwr[order(tmp_df$IRS4_orig)],rev(CI$upr[order(tmp_df$IRS4_orig)])),
        col="#00000010",lty = 0)
points(tmp_df$IRS4,tmp_df$IGFsig,bg=histo_pal[meta_dt$histologic_type[match(WT_samplenames,meta_dt$id)]],pch=21,cex=1.5)
axis(1,at = c(0,1,2,log10(400)),labels = c(0,10,100,400))
text(0,max(tmp_df$IGFsig),
     bquote(r~"="~.(round(summary(model1)$r.squared^0.5,3))
     ),adj=c(0,1))
text(0,max(tmp_df$IGFsig),
     bquote(
       italic(p)~"="~.(format(anova(model1)$'Pr(>F)'[1],scientific = T,digits=3))
     ),adj=c(0,2.5))
mtext("IRS4 expression (TPM)",1,2.0)
# mtext(bquote(paste(bold("GTF2I"^WT)~"(n"~"="~.(nrow(tmp_df)),")")))

# right top -> left bottom
plot(1000,xlim=c(0,max(tmp_df$IRS4)),ylim=c(min(tmp_df$lipid),max(tmp_df$lipid)),
     xlab="",ylab="",xaxt="n",las=2,yaxt='n')
axis(1,at = c(0,1,2,log10(400)),labels = c(0,10,100,400))
axis(2,las=2)
mtext("Intracellular lipid transport score",2,2.45)
lines(x=tmp_df$IRS4_orig, y=CI2$fit, lwd=1.2,col="grey50")
polygon(c(sort(tmp_df$IRS4_orig),rev(sort(tmp_df$IRS4_orig))),c(CI2$lwr[order(tmp_df$IRS4_orig)],rev(CI2$upr[order(tmp_df$IRS4_orig)])),
        col="#00000010",lty = 0)
points(tmp_df$IRS4,tmp_df$lipid,bg=histo_pal[meta_dt$histologic_type[match(WT_samplenames,meta_dt$id)]],pch=21,cex=1.5)
text(0,max(tmp_df$lipid),
     bquote(r~"="~.(round(summary(model2)$r.squared^0.5,3))
     ),adj=c(0,1))
text(0,max(tmp_df$lipid),
     bquote(
       italic(p)~"="~.(format(anova(model2)$'Pr(>F)'[1],scientific = T,digits=3))
     ),adj=c(0,2.5))
# mtext(bquote(paste(bold("GTF2I"^WT)~"(n"~"="~.(nrow(tmp_df)),")")),outer=T)
axis(1,at = c(0,1,2,log10(400)),labels = c(0,10,100,400))
mtext("IRS4 expression (TPM)",1,2.0)



# right bottom

plot(1000,xlim=c(0,max(tmp_df$IRS4)),ylim=c(min(tmp_df$mito),max(tmp_df$mito)),
     ylab="",las=2,xaxt="n",yaxt="n",xlab="")
axis(2,las=2)
axis(1,at = c(0,1,2,log10(400)),labels = c(0,10,100,400))
mtext("Mitochondrial distribution score",2,2.45)
lines(x=tmp_df$IRS4_orig, y=CI4$fit, lwd=1.2,col="grey50")
polygon(c(sort(tmp_df$IRS4_orig),rev(sort(tmp_df$IRS4_orig))),c(CI4$lwr[order(tmp_df$IRS4_orig)],rev(CI4$upr[order(tmp_df$IRS4_orig)])),
        col="#00000010",lty = 0)
points(tmp_df$IRS4,tmp_df$mito,bg=histo_pal[meta_dt$histologic_type[match(WT_samplenames,meta_dt$id)]],pch=21,cex=1.5)
text(0,max(tmp_df$mito),
     bquote(r~"="~.(round(summary(model4)$r.squared^0.5,3))
     ),adj=c(0,1))
text(0,max(tmp_df$mito),
     bquote(
       italic(p)~"="~.(format(anova(model4)$'Pr(>F)'[1],scientific = T,digits=3))
     ),adj=c(0,2.5))
mtext("IRS4 expression (TPM)",1,2.0)
legend("bottomright",legend=names(histo_pal[-c(1,3,7,8)]),pch=21,pt.bg=histo_pal[-c(1,3,7,8)],bty="n",pt.cex=1.5)
# })
dev.off()
















# ==============================================================================


# /home/users/sypark/00_Project/01_thymoma/10_Final_data/17_Figures_for_publication
library(tidyverse)
library(ggsci)
library(ggsci)
library(scales)
library(circlize)
# library(fgsea)
require(GSVA)
# require(GSEABase)

# meta_dt <- read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/02_metadata/Thymoma_summary_191115_1stSheet.txt')
# volc_dt <- read_tsv(paste0("~sypark/00_Project/01_thymoma/10_Final_data/01_expression/",
#                            'IRS4_corrected_v2/thym_137s_tpm_geneSymbol_tumorOnly_IRS4cor.tsv')) %>%
#   mutate_at(vars(-gene), list(~log10(.+0.1)))

# load exp data
meta_dt <- read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/02_metadata/Thymoma_summary_191227_1stSheet.txt')
bm_dt <- read_tsv(paste0('~sypark/02_Reference/13_biomart/',
                         'biomart_human_geneID_transcriptID_hgncsymbol',
                         '_genetype_ensemblgenename_190522.txt')) %>%
  dplyr::rename(gene = `Gene name`, gene_type = `Gene type`) %>% dplyr::select(gene, gene_type) %>% unique()
exp_dt <- read_tsv(paste0("~sypark/00_Project/01_thymoma/10_Final_data/01_expression/",
                          'IRS4_corrected_v2/',
                          'thym_137s_tpm_geneSymbol_tumorOnly_IRS4cor.tsv'))
exp_dt_pcg <- left_join(exp_dt, bm_dt, by="gene") %>%
  dplyr::filter(gene_type == 'protein_coding') %>% dplyr::select(-gene_type) # filter protein coding only
l10_exp_dt <- exp_dt %>% mutate_at(vars(-gene), list(~log10(.+0.01)))
l10_exp_dt <- left_join(l10_exp_dt, bm_dt, by="gene") %>%
  dplyr::filter(gene_type == 'protein_coding') %>% dplyr::select(-gene_type) # filter protein coding only
rm(bm_dt,exp_dt_pcg,exp_dt)


volc_dt <- l10_exp_dt


histo_pal = pal_aaas("default")(10)[c(4,1,7,8,6,2,3,5)]
names(histo_pal) = c("A","AB","MN-T","B1","B2","B3","NE","TC")

meanexp = list(
  m=rowMeans(volc_dt[,with(meta_dt,id[GTF2I_status2=='m'&final_cellularity>0.3])]),
  c=rowMeans(volc_dt[,with(meta_dt,id[GTF2I_status2=='c'&final_cellularity>0.3])]),
  w=rowMeans(volc_dt[,with(meta_dt,id[GTF2I_status2=='w'&final_cellularity>0.3])]),
  C=rowMeans(volc_dt[,with(meta_dt,id[GTF2I_status2!='c'&final_cellularity>0.3])])#,
)

stats= list(
  C=structure(with(meanexp,c-C),names=volc_dt$gene),
  mw=structure(with(meanexp,m-w),names=volc_dt$gene)#,
)


h_list <- GSEABase::getGmt('~kjyi/ref/msigdb/h.all.v6.2.symbols.gmt') %>% GSEABase::geneIds()
c5_list <- GSEABase::getGmt('~kjyi/ref/msigdb/c5.all.v6.2.symbols.gmt') %>% GSEABase::geneIds()
c2_list <- GSEABase::getGmt('~kjyi/ref/msigdb/c2.all.v6.2.symbols.gmt') %>% GSEABase::geneIds()
all_list <- append(append(h_list,c2_list),c5_list);rm(h_list,c2_list,c5_list)
all_list <- all_list[str_replace(names(all_list),"_.*","") %in% c("GO","HALLMARK","KEGG","REACTOME")]

volc_dt[1:3,1:3]

# purityfiltered_WT_samplenames=meta_dt$id[meta_dt$GTF2I_status2=="w" & meta_dt$final_cellularity>0.2]
WT_samplenames=meta_dt$id[meta_dt$GTF2I_status2=="w"]

# geneexpresionmatrix_WT_purityfiltered <- volc_dt[,c("gene",purityfiltered_WT_samplenames)] %>% as.data.frame %>% column_to_rownames("gene") %>% as.matrix

if(F){
  geneexpresionmatrix_WT <- volc_dt[,c("gene",WT_samplenames)] %>% as.data.frame %>% column_to_rownames("gene") %>% as.matrix
  all_gs_res <- gsva(geneexpresionmatrix_WT,all_list)
  
} else {
  gexpmat <- volc_dt %>% as.data.frame %>% column_to_rownames("gene") %>% as.matrix
  all_gs_res <- gsva(gexpmat,all_list)
  all_gs_res_WT <- all_gs_res[,WT_samplenames]
  gexpmat_WT <- gexpmat[,WT_samplenames]
}




# write_rds(all_gs_res,"data/gsva_WT_all_go_hm_kegg_reactome.Rds")

correlation_with_IRS4 = cor(geneexpresionmatrix_WT["IRS4",],t(all_gs_res))
correlation_with_IRS4=t(correlation_with_IRS4)%>%{colnames(.)="IRS4";.} %>% as.data.frame
correlation_with_IRS4=correlation_with_IRS4[order(correlation_with_IRS4$IRS4,decreasing=T),,drop=F]
head(correlation_with_IRS4,100)
"GO_INTRACELLULAR_LIPID_TRANSPORT"
"GO_MITOCHONDRION_DISTRIBUTION"
"GO_INSULIN_LIKE_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY"
"cor with purity"
# correlation_with_STAR = cor(geneexpresionmatrix_WT["STAR",],t(all_gs_res))
# correlation_with_STAR=t(correlation_with_STAR)%>%{colnames(.)="STAR";.} %>% as.data.frame
# correlation_with_STAR=correlation_with_STAR[order(correlation_with_STAR$STAR,decreasing=T),,drop=F]
# head(correlation_with_STAR,100)

# if(F){
#   ##############################################################################
#   "Supplementary table 2" deprecated
#   ##############################################################################
#   library(rJava,lib.loc="/home/users/kjyi/R/x86_64-redhat-linux-gnu-library/3.6/")
#   library(xlsx,lib.loc="/home/users/kjyi/R/x86_64-redhat-linux-gnu-library/3.6/")
#   degs <- c("IRS4","STAR","GPR87","C4orf50","ATP10B","NEFL")
#   degs %in% rownames(geneexpresionmatrix_WT)
#   correlation_with_degs <- cor(t(all_gs_res),t(geneexpresionmatrix_WT[degs,]))
#   correlation_with_degs <- correlation_with_degs[order(correlation_with_degs[,"IRS4"],decreasing=T),]
#   write.xlsx(correlation_with_degs, 
#              "tables/sup.table2.correlation_with_degs.xlsx",
#                           sheetName = "Sup.Table2.Correlation_with_ES_and_gene_expression_in_GTF2I_WT")
#   # for(deg in degs){
#   #   correlation_with_deg <- cor(geneexpresionmatrix_WT[deg,],t(all_gs_res))
#   #   correlation_with_deg <- t(correlation_with_deg) %>%
#   #     {colnames(.)=deg;.} %>% as.data.frame
#   #   correlation_with_deg <- correlation_with_deg[order(correlation_with_deg[[deg]],decreasing=T),,drop=F]
#   #   write.xlsx(correlation_with_deg, 
#   #              "tables/sup.table2.correlation_with_degs.xlsx",
#   #              sheetName = deg,
#   #              append = deg!="IRS4")
#   # }
#   # write.csv(correlation_with_IRS4,"tables/correlation_with_IRS4.csv")
#   # write.csv(correlation_with_STAR,"tables/correlation_with_STAR.csv")
# }

# correlation_with_IRS4["REACTOME_PI3K_CASCADE",]
# # GO_INSULIN_LIKE_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY
# GO_INTRACELLULAR_LIPID_TRANSPORT
# plot(pmax(0,geneexpresionmatrix_WT["IRS4",]),all_gs_res["GO_MAP_KINASE_ACTIVITY",],
#      bg=circlize::colorRamp2(0:1,c("red","black"))(meta_dt$final_cellularity[match(WT_samplenames,meta_dt$id)]),pch=21,cex=2)
# 
# plot(meta_dt$final_cellularity[match(WT_samplenames,meta_dt$id)],all_gs_res["GO_INTRACELLULAR_LIPID_TRANSPORT",],
#      bg=circlize::colorRamp2(c(0,0.3,0.7),c("blue","grey","red"))(meta_dt$final_cellularity[match(WT_samplenames,meta_dt$id)]),pch=21,cex=2)
# 
# 
# plot(pmax(0,geneexpresionmatrix_WT["IRS4",]),all_gs_res["GO_INTRACELLULAR_LIPID_TRANSPORT",],
#      bg=circlize::colorRamp2(c(0,0.3,0.7),c("blue","grey","red"))(meta_dt$final_cellularity[match(WT_samplenames,meta_dt$id)]),pch=21,cex=2)
# plot(pmax(0,geneexpresionmatrix_WT["IRS4",]),all_gs_res["GO_WNT_SIGNALOSOME",],
#      bg=circlize::colorRamp2(0:1,c("red","black"))(meta_dt$final_cellularity[match(WT_samplenames,meta_dt$id)]),pch=21,cex=2)
# plot(pmax(0,geneexpresionmatrix_WT["IRS4",]),all_gs_res["REACTOME_ENERGY_DEPENDENT_REGULATION_OF_MTOR_BY_LKB1_AMPK",],
#      bg=circlize::colorRamp2(0:1,c("red","black"))(meta_dt$final_cellularity[match(WT_samplenames,meta_dt$id)]),pch=21,cex=2)
# plot(pmax(0,geneexpresionmatrix_WT["IRS4",]),all_gs_res["GO_INSULIN_LIKE_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY",],
#      bg=circlize::colorRamp2(0:1,c("red","black"))(meta_dt$final_cellularity[match(WT_samplenames,meta_dt$id)]),pch=21,cex=2)
# plot(pmax(0,geneexpresionmatrix_WT["IRS4",]),all_gs_res["REACTOME_PI3K_CASCADE",],
#      bg=circlize::colorRamp2(0:1,c("red","black"))(meta_dt$final_cellularity[match(WT_samplenames,meta_dt$id)]),pch=21,cex=2)
# 
# 
# all_list["GO_INTRACELLULAR_LIPID_TRANSPORT"]
# all_list["REACTOME_ENERGY_DEPENDENT_REGULATION_OF_MTOR_BY_LKB1_AMPK"]
# all_list["GO_INSULIN_LIKE_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY"]
# all_list["GO_WNT_SIGNALOSOME"]
# all_list["GO_MITOCHONDRION_DISTRIBUTION"]
# 
# # rSquared <- summary(model1)$r.squared
# # pVal <- anova(model1)$'Pr(>F)'[1]
# 
# # cor(tmp_df$IRS4, tmp_df$IGFsig)
# # cor(tmp_df$IGFsig,tmp_df$IRS4)
# # cor(tmp_df$IGFsig,tmp_df$IRS4)^2

all_gs_res <- all_gs_res[,WT_samplenames]
geneexpresionmatrix_WT <- gexpmat[,WT_samplenames]
cor(t(all_gs_res),geneexpresionmatrix_WT["IRS4",]) %>% sort(decreasing=T)

tmp_df <- data.frame(IRS4=pmax(0,geneexpresionmatrix_WT["IRS4",]),
                     IRS4_orig=geneexpresionmatrix_WT["IRS4",],
                     IGFsig=all_gs_res["GO_INSULIN_LIKE_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY",],
                     lipid=all_gs_res["GO_INTRACELLULAR_LIPID_TRANSPORT",],
                     PI3Kcascade=all_gs_res["REACTOME_PI3K_CASCADE",],
                     mito=all_gs_res["GO_MITOCHONDRION_DISTRIBUTION",],
                     histol=histo_pal[meta_dt$histologic_type[match(WT_samplenames,meta_dt$id)]],
                     purity=meta_dt$final_cellularity[match(WT_samplenames,meta_dt$id)])
summary(lm(IGFsig ~ IRS4_orig+purity, data=tmp_df))
summary(lm(IGFsig ~ purity, data=tmp_df))

model1 <- lm(IGFsig ~ IRS4_orig, data=tmp_df)
model2 <- lm(lipid ~ IRS4_orig, data=tmp_df)
model3 <- lm(PI3Kcascade ~ IRS4_orig, data=tmp_df)
model4 <- lm(mito ~ IRS4_orig, data=tmp_df)

cor(tmp_df$IRS4_orig,tmp_df$IGFsig)
cor(tmp_df$IRS4_orig,tmp_df$lipid)
cor(tmp_df$IRS4_orig,tmp_df$PI3Kcascade)
cor(tmp_df$IRS4_orig,tmp_df$mito)
CI <- predict(model1, newdata=tmp_df, interval="prediction",level=0.9) %>% as.data.frame
CI2 <- predict(model2, newdata=tmp_df, interval="prediction",level=0.9) %>% as.data.frame
CI3 <- predict(model3, newdata=tmp_df, interval="prediction",level=0.9) %>% as.data.frame
CI4 <- predict(model4, newdata=tmp_df, interval="prediction",level=0.9) %>% as.data.frame


# svglite::htmlSVG(width=1150/254,height=1150/254,pointsize = 12*0.7, {
cairo_pdf("figures/IRS4_pathway_cor.1.pdf",width=1150/254,height=1050/254,pointsize = 12*0.7)
par(mar=c(4.5,4.5,0,0),oma=c(0,0,1,1),mfrow=c(2,2))

# left top
plot(1000,xlim=c(0,max(tmp_df$IRS4)),ylim=c(min(tmp_df$IGFsig),max(tmp_df$IGFsig)),
     xlab="",ylab="",xaxt="n",las=2)
mtext("IGF receptor signaling score",2,2.45)
lines(x=tmp_df$IRS4_orig, y=CI$fit, lwd=1.2,col="grey50")
polygon(c(sort(tmp_df$IRS4_orig),rev(sort(tmp_df$IRS4_orig))),c(CI$lwr[order(tmp_df$IRS4_orig)],rev(CI$upr[order(tmp_df$IRS4_orig)])),
        col="#00000010",lty = 0)
points(tmp_df$IRS4,tmp_df$IGFsig,bg=histo_pal[meta_dt$histologic_type[match(WT_samplenames,meta_dt$id)]],pch=21,cex=1.5)
axis(1,at = c(0,1,2,log10(400)),labels = c(0,10,100,400))
text(0,max(tmp_df$IGFsig),
     bquote(r~"="~.(round(summary(model1)$r.squared^0.5,3))
     ),adj=c(0,1))
text(0,max(tmp_df$IGFsig),
     bquote(
       italic(p)~"="~.(format(anova(model1)$'Pr(>F)'[1],scientific = T,digits=3))
     ),adj=c(0,2.5))
mtext("IRS4 expression (TPM)",1,2.0)
# mtext(bquote(paste(bold("GTF2I"^WT)~"(n"~"="~.(nrow(tmp_df)),")")))

# right top
plot(1000,xlim=c(0,max(tmp_df$IRS4)),ylim=c(min(tmp_df$lipid),max(tmp_df$lipid)),
     xlab="",ylab="",xaxt="n",las=2,yaxt='n')
axis(1,at = c(0,1,2,log10(400)),labels = c(0,10,100,400))
axis(2,las=2)
mtext("Intracellular lipid transport score",2,2.45)
lines(x=tmp_df$IRS4_orig, y=CI2$fit, lwd=1.2,col="grey50")
polygon(c(sort(tmp_df$IRS4_orig),rev(sort(tmp_df$IRS4_orig))),c(CI2$lwr[order(tmp_df$IRS4_orig)],rev(CI2$upr[order(tmp_df$IRS4_orig)])),
        col="#00000010",lty = 0)
points(tmp_df$IRS4,tmp_df$lipid,bg=histo_pal[meta_dt$histologic_type[match(WT_samplenames,meta_dt$id)]],pch=21,cex=1.5)
text(0,max(tmp_df$lipid),
     bquote(r~"="~.(round(summary(model2)$r.squared^0.5,3))
     ),adj=c(0,1))
text(0,max(tmp_df$lipid),
     bquote(
       italic(p)~"="~.(format(anova(model2)$'Pr(>F)'[1],scientific = T,digits=3))
     ),adj=c(0,2.5))
# mtext(bquote(paste(bold("GTF2I"^WT)~"(n"~"="~.(nrow(tmp_df)),")")),outer=T)
axis(1,at = c(0,1,2,log10(400)),labels = c(0,10,100,400))
mtext("IRS4 expression (TPM)",1,2.0)


# left bottom 
plot(1000,xlim=c(0,max(tmp_df$IRS4)),ylim=c(min(tmp_df$PI3Kcascade),max(tmp_df$PI3Kcascade)),
     ylab="",las=2,xaxt="n",xlab="")
axis(1,at = c(0,1,2,log10(400)),labels = c(0,10,100,400))
mtext("PI3K cascade score",2,2.45)
lines(x=tmp_df$IRS4_orig, y=CI3$fit, lwd=1.2,col="grey50")
polygon(c(sort(tmp_df$IRS4_orig),rev(sort(tmp_df$IRS4_orig))),c(CI3$lwr[order(tmp_df$IRS4_orig)],rev(CI3$upr[order(tmp_df$IRS4_orig)])),
        col="#00000010",lty = 0)
points(tmp_df$IRS4,tmp_df$PI3Kcascade,bg=histo_pal[meta_dt$histologic_type[match(WT_samplenames,meta_dt$id)]],pch=21,cex=1.5)
text(0,max(tmp_df$PI3Kcascade),
     bquote(r~"="~.(round(summary(model3)$r.squared^0.5,3))
     ),adj=c(0,1))
text(0,max(tmp_df$PI3Kcascade),
     bquote(
       italic(p)~"="~.(format(anova(model3)$'Pr(>F)'[1],scientific = T,digits=3))
     ),adj=c(0,2.5))
mtext("IRS4 expression (TPM)",1,2.0)
# legend("bottomright",legend=names(histo_pal[-c(1,3,7,8)]),pch=21,pt.bg=histo_pal[-c(1,3,7,8)],bty="n",pt.cex=1.5)


# right bottom

plot(1000,xlim=c(0,max(tmp_df$IRS4)),ylim=c(min(tmp_df$mito),max(tmp_df$mito)),
     ylab="",las=2,xaxt="n",yaxt="n",xlab="")
axis(2,las=2)
axis(1,at = c(0,1,2,log10(400)),labels = c(0,10,100,400))
mtext("Mitochondrial distribution score",2,2.45)
lines(x=tmp_df$IRS4_orig, y=CI4$fit, lwd=1.2,col="grey50")
polygon(c(sort(tmp_df$IRS4_orig),rev(sort(tmp_df$IRS4_orig))),c(CI4$lwr[order(tmp_df$IRS4_orig)],rev(CI4$upr[order(tmp_df$IRS4_orig)])),
        col="#00000010",lty = 0)
points(tmp_df$IRS4,tmp_df$mito,bg=histo_pal[meta_dt$histologic_type[match(WT_samplenames,meta_dt$id)]],pch=21,cex=1.5)
text(0,max(tmp_df$mito),
     bquote(r~"="~.(round(summary(model4)$r.squared^0.5,3))
     ),adj=c(0,1))
text(0,max(tmp_df$mito),
     bquote(
       italic(p)~"="~.(format(anova(model4)$'Pr(>F)'[1],scientific = T,digits=3))
     ),adj=c(0,2.5))
mtext("IRS4 expression (TPM)",1,2.0)
legend("bottomright",legend=names(histo_pal[-c(1,3,7,8)]),pch=21,pt.bg=histo_pal[-c(1,3,7,8)],bty="n",pt.cex=1.5)
# })
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



meta_dt$GTF2I_status2 %>% table()
cor_gexp_with_gsva <- cor(t(all_gs_res),t(geneexpresionmatrix_WT[c("STAR","IRS4","GPR87","C4orf50","ATP10B","NEFL"),]))
cor_gexp_with_gsva%>% as.data.frame()   %>% rownames_to_column("gs") %>% arrange(desc(IRS4))
cor_gexp_with_gsva%>% as.data.frame()   %>% rownames_to_column("gs") %>% arrange(desc(STAR))

write_csv(cor_gexp_with_gsva %>% as.data.frame %>% rownames_to_column("gs"),"data/cor_gexp_with_gsva.csv")




