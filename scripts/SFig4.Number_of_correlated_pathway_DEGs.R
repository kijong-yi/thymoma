"------------------------------------------------------------------------------"
"                           Supplementary Fig.4                                "
"                 IRS4, IGFBPL1, others cor.with IGFR sig.pathway              "
"------------------------------------------------------------------------------"
gtf2i_pal = c(m="#3C5488FF",w="#E64B35FF",c="#00A087FF")

fast_cor <- function(xt,yt=NULL){
  if(is.null(yt)){x <- t(xt) - colMeans(xt);return(tcrossprod(x / sqrt(rowSums(x ^ 2))))} else {
    x <- t(xt) - colMeans(xt);y <- t(yt) - colMeans(yt);
    return(tcrossprod(x / sqrt(rowSums(x ^ 2)),y / sqrt(rowSums(y ^ 2))))}}

# correlation with IRS4
cairo_pdf("~kjyi/Projects/thymus_single_cell/final2/figures/SFig4.IRS4_genes_scatter.pdf",
          width=65/25.4,height = 65/25.4,pointsize = 12*0.7*0.9)

# metadt
meta_dt <- read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/02_metadata/Thymoma_summary_191227_1stSheet.txt')
WT_samplenames=meta_dt$id[meta_dt$GTF2I_status2=="w"] 
MT_samplenames=meta_dt$id[meta_dt$GTF2I_status2=="m"] 
# exp
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
l10_exp_dt <- l10_exp_dt %>% as.data.frame() %>% column_to_rownames("gene") %>% .[,meta_dt$id]
rm(exp_dt_pcg, exp_dt, bm_dt)


# par(mfcol=c(1,3), mar=c(4,4,1,1))
par(mfcol=c(1,1), mar=c(4,4,1,1))
MM <- l10_exp_dt %>% as.matrix
# plot(MM["IGFBPL1",],MM["IRS4",],pch=21, cex=1.5, bg=gtf2i_pal[meta_dt$GTF2I_status2],
#      xlab="IGFBPL1 (log10(tpm+0.01))", ylab="IRS4 (log10(tpm+0.01))")
plot(10^MM["IGFBPL1",meta_dt$id]-0.01,10^MM["IRS4",meta_dt$id]-0.01,pch=21,cex=1.5,bg=gtf2i_pal[meta_dt$GTF2I_status2],
     xlab="IGFBPL1 (TPM)",ylab="IRS4 (TPM)")
# plot(MM["GPR87",],MM["IRS4",],pch=21, cex=1.5, bg=gtf2i_pal[meta_dt$GTF2I_status2],
#      xlab="GPR87 (log10(tpm+0.01))", ylab="IRS4 (log10(tpm+0.01))")
# plot(10^MM["GPR87",meta_dt$id]-0.01,10^MM["IRS4",meta_dt$id]-0.01,pch=21,cex=1.5,bg=gtf2i_pal[meta_dt$GTF2I_status2],
#      xlab="GPR87",ylab="IRS4 (log10(tpm+0.01))")
# plot(MM["CHST4",],MM["IRS4",],pch=21, cex=1.5, bg=gtf2i_pal[meta_dt$GTF2I_status2],
#      xlab="CHST4 (log10(tpm+0.01))", ylab="IRS4 (log10(tpm+0.01))")
# plot(10^MM["CHST4",meta_dt$id]-0.01,10^MM["IRS4",meta_dt$id]-0.01,pch=21,cex=1.5,bg=gtf2i_pal[meta_dt$GTF2I_status2],
#      xlab="CHST4",ylab="IRS4 (log10(tpm+0.01))")
dev.off()


cairo_pdf("~kjyi/Projects/thymus_single_cell/final2/figures/SFig4.IGFBPL1_cor.pdf",
          width=115/25.4,height = 65/25.4,pointsize = 12*0.7*0.9)


  histo_pal = c("#631879FF", "#3B4992FF", "#5F559BFF", "#A20056FF", "#BB0021FF", "#EE0000FF", "#008B45FF", "#008280FF")
  names(histo_pal) = c("A","AB","MN-T","B1","B2","B3","NE","TC")
  gtf2i_pal = c(m="#3C5488FF",w="#E64B35FF",c="#00A087FF")
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
  all_gs_res <- all_gs_res[,WT_samplenames]
  geneexpresionmatrix_WT <- l10_exp_dt_mat[,WT_samplenames]
   
  
  tmp_df <- data.frame(IGFBPL1=pmax(0,geneexpresionmatrix_WT["IGFBPL1",]),
                       IGFBPL1_orig=geneexpresionmatrix_WT["IGFBPL1",],
                       IGFsig=all_gs_res["GO_INSULIN_LIKE_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY",],
                       lipid=all_gs_res["GO_INTRACELLULAR_LIPID_TRANSPORT",],
                       PI3Kcascade=all_gs_res["REACTOME_PI3K_CASCADE",],
                       mito=all_gs_res["GO_MITOCHONDRION_DISTRIBUTION",],
                       histol=histo_pal[meta_dt$histologic_type[match(WT_samplenames,meta_dt$id)]],
                       purity=meta_dt$final_cellularity[match(WT_samplenames,meta_dt$id)])
  
  tmp_df <- data.frame(IGFBPL1=10^geneexpresionmatrix_WT["IGFBPL1",] + 0.01,
                       IGFBPL1_orig=10^geneexpresionmatrix_WT["IGFBPL1",] + 0.01,
                       IGFsig=all_gs_res["GO_INSULIN_LIKE_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY",],
                       lipid=all_gs_res["GO_INTRACELLULAR_LIPID_TRANSPORT",],
                       PI3Kcascade=all_gs_res["REACTOME_PI3K_CASCADE",],
                       mito=all_gs_res["GO_MITOCHONDRION_DISTRIBUTION",],
                       histol=histo_pal[meta_dt$histologic_type[match(WT_samplenames,meta_dt$id)]],
                       purity=meta_dt$final_cellularity[match(WT_samplenames,meta_dt$id)])
  
  summary(lm(IGFsig ~ IRS4_orig+purity, data=tmp_df))
  summary(lm(IGFsig ~ purity, data=tmp_df))
  
  model1 <- lm(IGFsig ~ IGFBPL1_orig, data=tmp_df)
  model2 <- lm(lipid ~ IGFBPL1_orig, data=tmp_df)
  model3 <- lm(PI3Kcascade ~ IGFBPL1_orig, data=tmp_df)
  model4 <- lm(mito ~ IGFBPL1_orig, data=tmp_df)
  
  cor(tmp_df$IGFBPL1_orig,tmp_df$IGFsig)
  # cor(tmp_df$IRS4_orig,tmp_df$lipid)
  # cor(tmp_df$IRS4_orig,tmp_df$PI3Kcascade)
  # cor(tmp_df$IRS4_orig,tmp_df$mito)
  CI <- predict(model1, newdata=tmp_df, interval="prediction",level=0.9) %>% as.data.frame
  CI2 <- predict(model2, newdata=tmp_df, interval="prediction",level=0.9) %>% as.data.frame
  CI3 <- predict(model3, newdata=tmp_df, interval="prediction",level=0.9) %>% as.data.frame
  CI4 <- predict(model4, newdata=tmp_df, interval="prediction",level=0.9) %>% as.data.frame
  # 
  
  # svglite::htmlSVG(width=1150/254,height=1150/254,pointsize = 12*0.7, {
  # cairo_pdf("figures/IRS4_pathway_cor.1.pdf",width=1150/254,height=1050/254,pointsize = 12*0.7)
  # 
  par(mar=c(4.5,4.5,0,0),oma=c(0,0,1,1),mfrow=c(1,2))
  
  
  
  # left bottom  -> left top
  plot(1000,xlim=c(0,max(tmp_df$IGFBPL1)),ylim=c(min(tmp_df$PI3Kcascade),max(tmp_df$PI3Kcascade)),
       ylab="",las=2,xlab="")
  # axis(1,at = c(log10(c(0.99,2,4,6,8,10)+0.01)),labels = c(0,2,4,6,8,10))
  # axis(1,at = c(0,2,4,6,8,10),labels = c(0,2,4,6,8,10))
  mtext("PI3K cascade score",2,2.45)
  lines(x=tmp_df$IGFBPL1_orig, y=CI3$fit, lwd=1.2,col="grey50")
  polygon(c(sort(tmp_df$IGFBPL1_orig),rev(sort(tmp_df$IGFBPL1_orig))),c(CI3$lwr[order(tmp_df$IGFBPL1_orig)],rev(CI3$upr[order(tmp_df$IGFBPL1_orig)])),
          col="#00000010",lty = 0)
  points(tmp_df$IGFBPL1,tmp_df$PI3Kcascade,bg=histo_pal[meta_dt$histologic_type[match(WT_samplenames,meta_dt$id)]],pch=21,cex=1.5)
  text(0,max(tmp_df$PI3Kcascade),
       bquote(r~"="~.(round(summary(model3)$r.squared^0.5,3))
       ),adj=c(0,1))
  text(0,max(tmp_df$PI3Kcascade),
       bquote(
         italic(p)~"="~.(format(anova(model3)$'Pr(>F)'[1],scientific = T,digits=3))
       ),adj=c(0,2.5))
  mtext("IGFBPL1 expression (TPM)",1,2.0)
  # legend("bottomright",legend=names(histo_pal[-c(1,3,7,8)]),pch=21,pt.bg=histo_pal[-c(1,3,7,8)],bty="n",pt.cex=1.5)
  
  
  # left top -> right top
  plot(1000,xlim=c(0,max(tmp_df$IGFBPL1)),ylim=c(min(tmp_df$IGFsig),max(tmp_df$IGFsig)),
       xlab="",ylab="",las=2)
  mtext("IGF receptor signaling score",2,2.45)
  lines(x=tmp_df$IGFBPL1_orig, y=CI$fit, lwd=1.2,col="grey50")
  polygon(c(sort(tmp_df$IGFBPL1_orig),rev(sort(tmp_df$IGFBPL1_orig))),c(CI$lwr[order(tmp_df$IGFBPL1_orig)],rev(CI$upr[order(tmp_df$IGFBPL1_orig)])),
          col="#00000010",lty = 0)
  points(tmp_df$IGFBPL1,tmp_df$IGFsig,bg=histo_pal[meta_dt$histologic_type[match(WT_samplenames,meta_dt$id)]],pch=21,cex=1.5)
  # axis(1,at = c(log10(c(0.99,2,4,6,8,10)+0.01)),labels = c(0,2,4,6,8,10))
  # axis(1,at = c(0,2,4,6,8,10),labels = c(0,2,4,6,8,10))
  text(0,max(tmp_df$IGFsig),
       bquote(r~"="~.(round(summary(model1)$r.squared^0.5,3))
       ),adj=c(0,1))
  text(0,max(tmp_df$IGFsig),
       bquote(
         italic(p)~"="~.(format(anova(model1)$'Pr(>F)'[1],scientific = T,digits=3))
       ),adj=c(0,2.5))
  mtext("IGFBPL1 expression (TPM)",1,2.0)
  # mtext(bquote(paste(bold("GTF2I"^WT)~"(n"~"="~.(nrow(tmp_df)),")")))
  # 
  # # right top -> left bottom
  # plot(1000,xlim=c(0,max(tmp_df$IRS4)),ylim=c(min(tmp_df$lipid),max(tmp_df$lipid)),
  #      xlab="",ylab="",xaxt="n",las=2,yaxt='n')
  # axis(1,at = c(0,1,2,log10(400)),labels = c(0,10,100,400))
  # axis(2,las=2)
  # mtext("Intracellular lipid transport score",2,2.45)
  # lines(x=tmp_df$IRS4_orig, y=CI2$fit, lwd=1.2,col="grey50")
  # polygon(c(sort(tmp_df$IRS4_orig),rev(sort(tmp_df$IRS4_orig))),c(CI2$lwr[order(tmp_df$IRS4_orig)],rev(CI2$upr[order(tmp_df$IRS4_orig)])),
  #         col="#00000010",lty = 0)
  # points(tmp_df$IRS4,tmp_df$lipid,bg=histo_pal[meta_dt$histologic_type[match(WT_samplenames,meta_dt$id)]],pch=21,cex=1.5)
  # text(0,max(tmp_df$lipid),
  #      bquote(r~"="~.(round(summary(model2)$r.squared^0.5,3))
  #      ),adj=c(0,1))
  # text(0,max(tmp_df$lipid),
  #      bquote(
  #        italic(p)~"="~.(format(anova(model2)$'Pr(>F)'[1],scientific = T,digits=3))
  #      ),adj=c(0,2.5))
  # # mtext(bquote(paste(bold("GTF2I"^WT)~"(n"~"="~.(nrow(tmp_df)),")")),outer=T)
  # axis(1,at = c(0,1,2,log10(400)),labels = c(0,10,100,400))
  # mtext("IRS4 expression (TPM)",1,2.0)
  # 
  # 
  # 
  # # right bottom
  # 
  # plot(1000,xlim=c(0,max(tmp_df$IRS4)),ylim=c(min(tmp_df$mito),max(tmp_df$mito)),
  #      ylab="",las=2,xaxt="n",yaxt="n",xlab="")
  # axis(2,las=2)
  # axis(1,at = c(0,1,2,log10(400)),labels = c(0,10,100,400))
  # mtext("Mitochondrial distribution score",2,2.45)
  # lines(x=tmp_df$IRS4_orig, y=CI4$fit, lwd=1.2,col="grey50")
  # polygon(c(sort(tmp_df$IRS4_orig),rev(sort(tmp_df$IRS4_orig))),c(CI4$lwr[order(tmp_df$IRS4_orig)],rev(CI4$upr[order(tmp_df$IRS4_orig)])),
  #         col="#00000010",lty = 0)
  # points(tmp_df$IRS4,tmp_df$mito,bg=histo_pal[meta_dt$histologic_type[match(WT_samplenames,meta_dt$id)]],pch=21,cex=1.5)
  # text(0,max(tmp_df$mito),
  #      bquote(r~"="~.(round(summary(model4)$r.squared^0.5,3))
  #      ),adj=c(0,1))
  # text(0,max(tmp_df$mito),
  #      bquote(
  #        italic(p)~"="~.(format(anova(model4)$'Pr(>F)'[1],scientific = T,digits=3))
  #      ),adj=c(0,2.5))
  # mtext("IRS4 expression (TPM)",1,2.0)
  # legend("bottomright",legend=names(histo_pal[-c(1,3,7,8)]),pch=21,pt.bg=histo_pal[-c(1,3,7,8)],bty="n",pt.cex=1.5)
  # # })
  # dev.off()
  # 

dev.off()

# corelation with pathway
cairo_pdf("~kjyi/Projects/thymus_single_cell/final2/figures/SFig4.pathway_genes_scatter.pdf",
          width=150/25.4,height = 80/25.4,pointsize = 12*0.7*0.9)
labels =c("GO_INSULIN_LIKE_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY" = "IGFR signaling pathway",
          "REACTOME_PI3K_CASCADE" = "PI3K cascade")
par(mfcol=c(2,4), mar=c(4,4,1,1))
for(g in c("IRS4","IGFBPL1","GPR87","CHST4")){
  for(gs in c("GO_INSULIN_LIKE_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY",
              "REACTOME_PI3K_CASCADE")){
    plot(MM[g,],all_gs_res[gs,],
         pch=21, cex=1.7, bg=gtf2i_pal[meta_dt$GTF2I_status2],
         xlab=g, ylab=labels[gs],lwd=0.7)
  }
}
dev.off()

plot(MM["IGFBPL1",],all_gs_res["GO_INSULIN_LIKE_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY",],
     pch=21, cex=1.5, bg=gtf2i_pal[meta_dt$GTF2I_status2],
     xlab="IGFBPL1", ylab="IGFRsig.path")
plot(MM["IGFBPL1",],all_gs_res["REACTOME_PI3K_CASCADE",],
     pch=21, cex=1.5, bg=gtf2i_pal[meta_dt$GTF2I_status2],
     xlab="IGFBPL1", ylab="PI3Kcascade")
plot(MM["IGFBPL1",],all_gs_res["GO_INTRACELLULAR_LIPID_TRANSPORT",],
     pch=21, cex=1.5, bg=gtf2i_pal[meta_dt$GTF2I_status2],
     xlab="IGFBPL1", ylab="Intracellualr lipid transport")
plot(MM["IGFBPL1",],all_gs_res["GO_MITOCHONDRION_DISTRIBUTION",],
     pch=21, cex=1.5, bg=gtf2i_pal[meta_dt$GTF2I_status2],
     xlab="IGFBPL1", ylab="mitochonria dist")


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# degs <- read_rds("data/degs.Rds")
WT_high_genes <- degs$WT_high
WT_high_genes <- volc_dt %>% filter(WT_mean - MT_mean >= 1 & p.adjust < 10e-20) %>% .$gene
WT_high_genes <- volc_dt %>% arrange(MT_mean-WT_mean) %>% .$gene %>% .[1:100]
which(WT_high_genes=="IGFBPL1")
which(WT_high_genes=="PAH")

WT_high_genes <- volc_dt %>% filter(WT_mean - MT_mean >= 0) %>% arrange(p.adjust) %>% .$gene %>% .[1:100]
which(WT_high_genes=="IGFBPL1")

volc_dt[volc_dt$gene == "IGFBPL1",] %>% with(p.adjust)

WT_high_genes



cormat <- fast_cor(t(all_gs_res5[,WT_samplenames]),t(volc_dt2[,WT_samplenames]))


# cormat <- fast_cor(t(all_gs_res2),(log_exp_dt2[WT_samplenames,]))
logfc = colMeans(log_exp_dt2[WT_samplenames,]) - colMeans(log_exp_dt2[MT_samplenames,])
logfc2 = colMeans(log_exp_dt4[WT_samplenames,]) - colMeans(log_exp_dt4[MT_samplenames,])


logfc5 = rowMeans(volc_dt2[,WT_samplenames]) - rowMeans(volc_dt2[,MT_samplenames])
logfc5["IRS4"]
plot(logfc5, cormat["GO_INSULIN_LIKE_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY",],
     col=ifelse(names(logfc5)=="IRS4","red",ifelse(names(logfc5)=="IGFBPL1","green",ifelse(names(logfc5)=="PAH","red","black"))))
abline(v=1,h=0.7)


cormat <- fast_cor(t(all_gs_res5[,WT_samplenames]),t(volc_dt2[,WT_samplenames]))



mydf=data.frame(logfc=logfc5,
           cor_igfrsigpath=cormat["GO_INSULIN_LIKE_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY",],
           gene=names(logfc5))



plotly::plot_ly(mydf,x=~logfc,y=~cor_igfrsigpath,text=~gene)



names(logfc2)[logfc2 > 1 & cormat["GO_INSULIN_LIKE_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY",] > 0.7]



names(logfc)

cormat4[1:3,1:3]

logfc[1:3]
plot(logfc, cormat4["WNT_UP.V1_UP",])





rownames(cormat4)


cormat2 <- cor(t(all_gs_res2),(log_exp_dt2[WT_samplenames,WT_high_genes]),method = "spearman")
# cormat3 <- fast_cor(t(all_gs_res3),(log_exp_dt2[WT_samplenames,WT_high_genes]))

# mat = cormat[grepl("^HALLMARK",rownames(cormat)),]
# mat = cormat[grepl("^GO",rownames(cormat)),]
# mat = cormat[grepl("^KEGG",rownames(cormat)),]
mat = cormat
mat = cormat2
# mat = cormat3

sort(colSums(mat>0.7),decreasing=T)[1:10]
sort(colSums(mat>0.7),decreasing=T)["IRS4"]

mat[,"IRS4"] %>% sort(decreasing=T) %>% .[1:15]
mat[,"CHST4"] %>% sort(decreasing=T) %>% .[1:10]
mat[,"GPR87"] %>% sort(decreasing=T) %>% .[1:10]

which(names(sort(colSums(mat>0.7),decreasing=T))=="IRS4")

sort(mat["GO_INSULIN_LIKE_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY",],decreasing=T)[1:10]
sort(mat["REACTOME_PI3K_CASCADE",],decreasing=T)[1:10]
which(names(sort(mat["GO_INSULIN_LIKE_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY",],decreasing=T))=="IRS4")
which(names(sort(mat["REACTOME_PI3K_CASCADE",],decreasing=T))=="IRS4")


sort(mat["GO_INTRACELLULAR_LIPID_TRANSPORT",],decreasing=T)[1:10]
sort(mat["REACTOME_PI3K_CASCADE",],decreasing=T)[1:10]
length(WT_high_genes)


which(names(sort(mat["GO_INTRACELLULAR_LIPID_TRANSPORT",],decreasing=T))=="IRS4")

sort(mat["GO_INSULIN_LIKE_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY",],decreasing=T)[1:10]
sort(mat["REACTOME_PI3K_CASCADE",],decreasing=T)[1:10]

which(WT_high_genes == "PAH")
which(WT_high_genes == "IGFBPL1")
which(WT_high_genes == "TXLNB")

which(names(sort(mat["GO_INSULIN_LIKE_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY",],decreasing=T))=="IRS4")

log_exp_dt2
plot(log_exp_dt2[,"IGFBPL1"],log_exp_dt2[,"IRS4"])

plot(log_exp_dt3[,"IGFBPL1"],log_exp_dt3[,"IRS4"])

plot(log_exp_dt2[,"IGFBPL1"],log_exp_dt2[,"IRS1"])
plot(log_exp_dt2[,"IGFBPL1"],log_exp_dt2[,"IRS2"])
plot(log_exp_dt2[,"IGFBPL1"])
table(gtf2i = meta_dt$GTF2I_status2=="m", ab= meta_dt$histologic_type %in% c("A","AB"))
table(gtf2i = meta_dt$GTF2I_status2=="m", b123 = meta_dt$histologic_type %in% c("B1","B2","B3"))


fit = lm(all_gs_res2["GO_INSULIN_LIKE_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY",WT_samplenames]
           scale(log_exp_dt2[WT_samplenames,WT_high_genes[1:50]]))
summary(fit)

