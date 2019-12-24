# /home/users/sypark/00_Project/01_thymoma/10_Final_data/17_Figures_for_publication
library(tidyverse)
library(ggsci)
library(ggsci)
library(scales)
library(circlize)
library(fgsea)
require(GSVA)
require(GSEABase)

meta_dt <- read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/02_metadata/Thymoma_summary_191115_1stSheet.txt')
volc_dt <- read_tsv(paste0("~sypark/00_Project/01_thymoma/10_Final_data/01_expression/",
                           'IRS4_corrected_v2/thym_137s_tpm_geneSymbol_tumorOnly_IRS4cor.tsv')) %>%
  mutate_at(vars(-gene), list(~log10(.+0.1)))

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


h_list <- getGmt('~kjyi/ref/msigdb/h.all.v6.2.symbols.gmt') %>% geneIds()
c5_list <- getGmt('~kjyi/ref/msigdb/c5.all.v6.2.symbols.gmt') %>% geneIds()
c2_list <- getGmt('~kjyi/ref/msigdb/c2.all.v6.2.symbols.gmt') %>% geneIds()
all_list <- append(append(h_list,c2_list),c5_list);rm(h_list,c2_list,c5_list)
all_list <- all_list[str_replace(names(all_list),"_.*","") %in% c("GO","HALLMARK","KEGG","REACTOME")]

volc_dt[1:3,1:3]

# purityfiltered_WT_samplenames=meta_dt$id[meta_dt$GTF2I_status2=="w" & meta_dt$final_cellularity>0.2]
WT_samplenames=meta_dt$id[meta_dt$GTF2I_status2=="w"]

# geneexpresionmatrix_WT_purityfiltered <- volc_dt[,c("gene",purityfiltered_WT_samplenames)] %>% as.data.frame %>% column_to_rownames("gene") %>% as.matrix
geneexpresionmatrix_WT <- volc_dt[,c("gene",WT_samplenames)] %>% as.data.frame %>% column_to_rownames("gene") %>% as.matrix

all_gs_res <- gsva(geneexpresionmatrix_WT,all_list)
write_rds(all_gs_res,"data/gsva_WT_all_go_hm_kegg_reactome.Rds")

correlation_with_IRS4 = cor(geneexpresionmatrix_WT["IRS4",],t(all_gs_res))
correlation_with_IRS4=t(correlation_with_IRS4)%>%{colnames(.)="IRS4";.} %>% as.data.frame
correlation_with_IRS4=correlation_with_IRS4[order(correlation_with_IRS4$IRS4,decreasing=T),,drop=F]
head(correlation_with_IRS4,100)

correlation_with_STAR = cor(geneexpresionmatrix_WT["STAR",],t(all_gs_res))
correlation_with_STAR=t(correlation_with_STAR)%>%{colnames(.)="STAR";.} %>% as.data.frame
correlation_with_STAR=correlation_with_STAR[order(correlation_with_STAR$STAR,decreasing=T),,drop=F]
head(correlation_with_STAR,100)

if(F){
  write.csv(correlation_with_IRS4,"tables/correlation_with_IRS4.csv")
  write.csv(correlation_with_STAR,"tables/correlation_with_STAR.csv")
}

# correlation_with_IRS4["REACTOME_PI3K_CASCADE",]
# # GO_INSULIN_LIKE_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY
# GO_INTRACELLULAR_LIPID_TRANSPORT
plot(pmax(0,geneexpresionmatrix_WT["IRS4",]),all_gs_res["GO_MAP_KINASE_ACTIVITY",],
     bg=circlize::colorRamp2(0:1,c("red","black"))(meta_dt$final_cellularity[match(WT_samplenames,meta_dt$id)]),pch=21,cex=2)

plot(meta_dt$final_cellularity[match(WT_samplenames,meta_dt$id)],all_gs_res["GO_INTRACELLULAR_LIPID_TRANSPORT",],
     bg=circlize::colorRamp2(c(0,0.3,0.7),c("blue","grey","red"))(meta_dt$final_cellularity[match(WT_samplenames,meta_dt$id)]),pch=21,cex=2)


plot(pmax(0,geneexpresionmatrix_WT["IRS4",]),all_gs_res["GO_INTRACELLULAR_LIPID_TRANSPORT",],
     bg=circlize::colorRamp2(c(0,0.3,0.7),c("blue","grey","red"))(meta_dt$final_cellularity[match(WT_samplenames,meta_dt$id)]),pch=21,cex=2)
plot(pmax(0,geneexpresionmatrix_WT["IRS4",]),all_gs_res["GO_WNT_SIGNALOSOME",],
     bg=circlize::colorRamp2(0:1,c("red","black"))(meta_dt$final_cellularity[match(WT_samplenames,meta_dt$id)]),pch=21,cex=2)
plot(pmax(0,geneexpresionmatrix_WT["IRS4",]),all_gs_res["REACTOME_ENERGY_DEPENDENT_REGULATION_OF_MTOR_BY_LKB1_AMPK",],
     bg=circlize::colorRamp2(0:1,c("red","black"))(meta_dt$final_cellularity[match(WT_samplenames,meta_dt$id)]),pch=21,cex=2)
plot(pmax(0,geneexpresionmatrix_WT["IRS4",]),all_gs_res["GO_INSULIN_LIKE_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY",],
     bg=circlize::colorRamp2(0:1,c("red","black"))(meta_dt$final_cellularity[match(WT_samplenames,meta_dt$id)]),pch=21,cex=2)
plot(pmax(0,geneexpresionmatrix_WT["IRS4",]),all_gs_res["REACTOME_PI3K_CASCADE",],
     bg=circlize::colorRamp2(0:1,c("red","black"))(meta_dt$final_cellularity[match(WT_samplenames,meta_dt$id)]),pch=21,cex=2)


all_list["GO_INTRACELLULAR_LIPID_TRANSPORT"]
all_list["REACTOME_ENERGY_DEPENDENT_REGULATION_OF_MTOR_BY_LKB1_AMPK"]
all_list["GO_INSULIN_LIKE_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY"]
all_list["GO_WNT_SIGNALOSOME"]
all_list["GO_MITOCHONDRION_DISTRIBUTION"]
mtexti <- function(text, side, off = 0.25,
                   srt = if(side == 2) 90  else
                     if(side == 4) 270 else 0, ...) {
  # dimensions of plotting region in user units
  usr <- par('usr')
  # dimensions of plotting region in inches
  pin <- par('pin')
  # user units per inch
  upi <- c(usr[2]-usr[1],
           usr[4]-usr[3]) / pin
  # default x and y positions
  xpos <- (usr[1] + usr[2])/2
  ypos <- (usr[3] + usr[4])/2
  if(1 == side)
    ypos <- usr[3] - upi[2] * off
  if(2 == side)
    xpos <- usr[1] - upi[1] * off
  if(3 == side)
    ypos <- usr[4] + upi[2] * off
  if(4 == side)
    xpos <- usr[2] + upi[1] * off
  text(x=xpos, y=ypos, text, xpd=NA, srt=srt, ...)
}

# rSquared <- summary(model1)$r.squared
# pVal <- anova(model1)$'Pr(>F)'[1]

# cor(tmp_df$IRS4, tmp_df$IGFsig)
# cor(tmp_df$IGFsig,tmp_df$IRS4)
# cor(tmp_df$IGFsig,tmp_df$IRS4)^2

cairo_pdf("figures/IRS4_pathway_cor.1.pdf",height = 12/2.54,width=7/2.54,pointsize = 12*0.7)
par(mar=c(0,0,0,0),oma=c(4.5,4.5,1.7,0.7),mfrow=c(2,1))
tmp_df <- data.frame(IRS4=pmax(0,geneexpresionmatrix_WT["IRS4",]),
                     IRS4_orig=geneexpresionmatrix_WT["IRS4",],
                     IGFsig=all_gs_res["GO_INSULIN_LIKE_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY",],
                     # lipid=all_gs_res["REACTOME_PI3K_CASCADE",],
                     lipid=all_gs_res["GO_INTRACELLULAR_LIPID_TRANSPORT",],
                     histol=histo_pal[meta_dt$histologic_type[match(WT_samplenames,meta_dt$id)]],
                     purity=meta_dt$final_cellularity[match(WT_samplenames,meta_dt$id)])
summary(lm(IGFsig ~ IRS4_orig+purity, data=tmp_df))
summary(model0)
summary(lm(IGFsig ~ purity, data=tmp_df))

model1 <- lm(IGFsig ~ IRS4_orig, data=tmp_df)
summary(model1)


model2 <- lm(lipid ~ IRS4_orig, data=tmp_df)
cor(tmp_df$IRS4_orig,tmp_df$IGFsig)
cor(tmp_df$IRS4_orig,tmp_df$IGFsig)
CI <- predict(model1, newdata=tmp_df, interval="prediction",level=0.9) %>% as.data.frame
CI2 <- predict(model2, newdata=tmp_df, interval="prediction",level=0.9) %>% as.data.frame
plot(1000,xlim=c(0,max(tmp_df$IRS4)),ylim=c(min(tmp_df$IGFsig),max(tmp_df$IGFsig)),
     xlab="",ylab="",xaxt="n",las=2)
mtext("IGF receptor signaling score",2,2.45)
lines(x=tmp_df$IRS4_orig, y=CI$fit, lwd=1.2,col="grey50")
polygon(c(sort(tmp_df$IRS4_orig),rev(sort(tmp_df$IRS4_orig))),c(CI$lwr[order(tmp_df$IRS4_orig)],rev(CI$upr[order(tmp_df$IRS4_orig)])),
        col="#00000010",lty = 0)
points(tmp_df$IRS4,tmp_df$IGFsig,bg=histo_pal[meta_dt$histologic_type[match(WT_samplenames,meta_dt$id)]],pch=21,cex=1.5)
text(0,max(tmp_df$IGFsig),
     bquote(r~"="~.(round(summary(model1)$r.squared^0.5,3))
     ),adj=c(0,1))
text(0,max(tmp_df$IGFsig),
     bquote(
       italic(p)~"="~.(format(anova(model1)$'Pr(>F)'[1],scientific = T,digits=3))
     ),adj=c(0,2.5))
# mtext(bquote(GTF2I^WT~"(n"~"="~.(nrow(tmp_df))+")"),font=2)
mtext(bquote(paste(bold("GTF2I"^WT)~"(n"~"="~.(nrow(tmp_df)),")")))
plot(1000,xlim=c(0,max(tmp_df$IRS4)),ylim=c(min(tmp_df$lipid),max(tmp_df$lipid)),
     xlab=expression(log[10]~TPM~of~IRS4),ylab="PI3K cascade score",las=2,xaxt="n")
axis(1,at = c(0,1,2,log10(400)),labels = c(0,10,100,400))
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
# mtexti(expression(log[10]~TPM~of~IRS4),1,0.5)
mtext("IRS4 expression (TPM)",1,2.0)
legend("bottomright",legend=names(histo_pal[-c(1,3,7,8)]),pch=21,pt.bg=histo_pal[-c(1,3,7,8)],bty="n",pt.cex=1.5)

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




