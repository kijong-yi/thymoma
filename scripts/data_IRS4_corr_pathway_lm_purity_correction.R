read_tsv("abc")
# /home/users/sypark/00_Project/01_thymoma/10_Final_data/17_Figures_for_publication
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




# try many genes
Y=all_gs_res[,meta_dt$id[meta_dt$GTF2I_status2=="w"]]
x1=meta_dt$final_cellularity[meta_dt$GTF2I_status2=="w"]
X2=l10_exp_dt_mat[,meta_dt$id[meta_dt$GTF2I_status2=="w"]]

try_one <- function(g){
  lm_p = 
    lapply(1:nrow(all_gs_res), function(i){
      summary(lm(Y[i,]~x1+X2[g,]))$coefficients -> b
      c(b[3,1]) # c(x[2:3,4],x[2:3,1])
    })
  unlist(lm_p)
}


degs <- read_rds("data/degs.Rds")
degs$WT_high
registerDoMC(20)
start=Sys.time()

coefs <- foreach(g = degs$WT_high,.combine=cbind) %dopar% {
  try_one(g)
}
dim(coefs)
colnames(coefs) = degs$WT_high
rownames(coefs) = rownames(all_gs_res)
coefs[1:3,1:3]
write_csv(rownames_to_column(as.data.frame(coefs), "geneset"), "tables/coefs.lm.wt_high_DEGS_purity.csv")

# ---------------------------------------------------------------------------- #
"                        regress out strategy - not working                    "
# ------------------------------------------------------------------------------
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

X <- l10_exp_dt %>% as.data.frame() %>% column_to_rownames("gene") %>% t
y <- meta_dt$final_cellularity
dim(X)
X[1:3,1:3]
Xr = X
for(i in 1:ncol(X)){
  qr <- lm(y~X[,i],qr = T)$qr
  Xr[,i] <- qr.resid(qr = qr, y = X[,i])
}

#
h_list <- GSEABase::getGmt('~kjyi/ref/msigdb/h.all.v6.2.symbols.gmt') %>% GSEABase::geneIds()
c5_list <- GSEABase::getGmt('~kjyi/ref/msigdb/c5.all.v6.2.symbols.gmt') %>% GSEABase::geneIds()
c2_list <- GSEABase::getGmt('~kjyi/ref/msigdb/c2.all.v6.2.symbols.gmt') %>% GSEABase::geneIds()
all_list <- append(append(h_list,c2_list),c5_list);rm(h_list,c2_list,c5_list)
all_list <- all_list[str_replace(names(all_list),"_.*","") %in% c("GO","HALLMARK","KEGG","REACTOME")]

all_gs_res <- gsva(t(X), all_list)
dim(all_gs_res)

#
WT_samplenames=meta_dt$id[meta_dt$GTF2I_status2=="w"]
MT_samplenames=meta_dt$id[meta_dt$GTF2I_status2=="m"]

fast_cor <- function(xt,yt=NULL){
  if(is.null(yt)){
    x <- t(xt) - colMeans(xt)
    return(tcrossprod(x / sqrt(rowSums(x ^ 2))))
  } else {
    x <- t(xt) - colMeans(xt)
    y <- t(yt) - colMeans(yt)
    return(tcrossprod(x / sqrt(rowSums(x ^ 2)),y / sqrt(rowSums(y ^ 2))))  
  }
}

cormat <- fast_cor(X[WT_samplenames,],t(all_gs_res[,WT_samplenames]))
cormat_r <- fast_cor(Xr[WT_samplenames,],t(all_gs_res[,WT_samplenames]))

cormat["IRS4",] %>% sort(decreasing=T) %>% .[1:10]
cormat_r["IRS4",] %>% sort(decreasing=T) %>% .[1:10]

#
"mutation burden"
(meta_dt$n_pointmt/meta_dt$bait_size)[meta_dt$histologic_type %in% c("A","AB","B1","B2","B3")] %>% summary
(meta_dt$n_pointmt/meta_dt$bait_size)[meta_dt$histologic_type %in% c("TC")] %>% summary
x1 = (meta_dt$n_pointmt/meta_dt$bait_size)[meta_dt$histologic_type %in% c("A","AB","B1","B2","B3")]
x2 = (meta_dt$n_pointmt/meta_dt$bait_size)[meta_dt$histologic_type %in% c("TC")]
x2[-which.max(x2)] %>% summary
t.test(x1,x2)
t.test(x1,x2[-which.max(x2)])
wilcox.test(x1,x2[-which.max(x2)])
wilcox.test(x1,x2)
meta_dt$histologic_type %in% c("A","AB","B1","B2","B3")
meta_dt$histologic_type %in% "TC"


# ---------------------------------------------------------------------------- #





coefs[1:3,1:3]


meta_dt$GTF2I_status2=="w"
rownames(lm_p) = rownames(all_gs_res)
View(lm_p)

colnames(lm_p) = c("p.purity", "p.logIRS4", "coef.purity", "coef.logIRS4")

lm_p2 = foreach(i = 1:nrow(all_gs_res),.combine=rbind) %do% {
  summary(lm(all_gs_res[i,meta_dt$GTF2I_status2=="w"] ~ 
               meta_dt$final_cellularity[meta_dt$GTF2I_status2=="w"] + 
               log_exp_dt["STAR",meta_dt$GTF2I_status2=="w"]))$coefficients -> x
  c(x[2:3,4],x[2:3,1])
}
rownames(lm_p2) = rownames(all_gs_res)
View(lm_p2)


colnames(lm_p2) = c("p.purity", "p.logSTAR", "coef.purity", "coef.logSTAR")


lm_p3 = foreach(i = 1:nrow(all_gs_res),.combine=rbind) %do% {
  summary(lm(all_gs_res[i,meta_dt$GTF2I_status2=="w"] ~ 
               meta_dt$final_cellularity[meta_dt$GTF2I_status2=="w"] + 
               log_exp_dt["NEFL",meta_dt$GTF2I_status2=="w"]))$coefficients -> x
  c(x[2:3,4],x[2:3,1])
}
rownames(lm_p3) = rownames(all_gs_res)
View(lm_p3)

colnames(lm_p2) = c("p.purity", "p.logNEFL", "coef.purity", "coef.logNEFL")

lm.summaries = list(IRS4 = lm_p, STAR = lm_p2, NEFL = lm_p3)

write_csv(rownames_to_column(as.data.frame(lm_p), "geneset"), "tables/summary.lm.IRS4+purity.wt.csv")
write_csv(rownames_to_column(as.data.frame(lm_p2), "geneset"), "tables/summary.lm.STAR+purity.wt.csv")
write_csv(rownames_to_column(as.data.frame(lm_p3), "geneset"), "tables/summary.lm.NEFL+purity.wt.csv")

lm_p[1:3,1:3]

lm_p[lm_p[,4]>0,] %>% .[order(.[,2]),] %>% head(20)

lm_p2[lm_p2[,4]>0,] %>% .[order(.[,2]),] %>% head(20)
                               

gsoi = c("GO_GAMMA_CATENIN_BINDING",
         "GO_INTRACELLULAR_LIPID_TRANSPORT", 
         "GO_MITOCHONDRION_DISTRIBUTION",
         "GO_PHOSPHATIDYLINOSITOL_3_PHOSPHATE_BINDING",
         "GO_SCF_UBIQUITIN_LIGASE_COMPLEX")

p1 = prcomp(t(all_gs_res[gsoi,]))
s1 = summary(p1)

colnames(all_gs_res) == meta_dt$id

par(mfrow=c(1,2))
biplot.prcomp2(p1, bg=gtf2i_pal[meta_dt$GTF2I_status2],pt.cex=2.5,
               axis.text.cex = 1.5,
               xlab="",ylab="",
               pch=21)

biplot.prcomp2(p1, bg=circlize::colorRamp2(seq(0,1,length.out=3),
                                           c("red","green","blue")
)(meta_dt$final_cellularity),
pt.cex=2.5,
axis.text.cex = 1.5,
xlab="",ylab="",
pch=21)
