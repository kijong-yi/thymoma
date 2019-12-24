
library(tidyverse)
library(Seurat)


# load data ---------------------------------------------------------------------------------------
merged <- read_rds("data/singlecell/merged.Rds")
gene.use <- read_rds("data/gene.use.Rds")
thymoma_tpm <- read_tsv(paste0("~sypark/00_Project/01_thymoma/10_Final_data/01_expression/",
                               'IRS4_corrected_v2/',
                               'thym_137s_tpm_geneSymbol_tumorOnly_IRS4cor.tsv')) %>% 
  as.data.frame() %>%
  column_to_rownames("gene")


# prep thymoma data --------------------------------------------------------------------------------
# tpm -> log10(.+0.1) -> gene-wise scale(mean=0,sd=1) -> done
thymoma_tpm_mggene <- thymoma_tpm %>% 
  rownames_to_column("hgname") %>% 
  dplyr::left_join(gene.use, by = "hgname") %>%
  na.omit() %>%
  dplyr::select(-hgname) %>% 
  group_by(mgname) %>% 
  summarize_all(sum) %>%
  column_to_rownames("mgname") # %>%
# > thymoma_tpm_mggene[1:5,1:7]
#               SNU_03_C SNU_04_C SNU_05_C SNU_06_C SNU_09_C SNU_10_C TCGA-4X-A9FD
# 0610040J01Rik     0.37     0.93     1.04     1.05     0.22     0.68         0.64
# 1520401A03Rik     0.04     5.73     1.27     4.14     3.03     3.73         9.24
thymoma_mggene_scaled <- scale(log10(t(thymoma_tpm_mggene) + .1)) %>% {.[is.nan(.)]=0;.}


# prep thymus single cell data -------------------------------------------------------------------
# thymus: CPM -> knnpool(5) -> log10(.+0.01) -> scale(mean=0,sd=1)
boxplot(merged$nCount_RNA~merged$simplified)
boxplot(merged$nCount_RNA~merged$simplified2)
boxplot(merged$nCount_RNA~merged$age)

CalculateCpm <- function(object) {
  cpm <- GetAssayData(object, slot = "counts") %>%
  {./Matrix::Matrix(rep(Matrix::colSums(.),nrow(.)), ncol = ncol(.),byrow = T)} * 1E6
  cpm[is.nan(as.matrix(cpm))] <- 0
  SetAssayData(object, assay.type = "RNA", slot = "data", new.data = cpm)
}

merged_cpm <- CalculateCpm(merged)

thymus_cluster_cpm <- GetAssayData(merged_cpm, slot = "data")[gene.use$mgname,] %>%
  Matrix::t() %>%
  as.data.frame() %>%
  {rbind(
    cbind(y = "progenitor",.[merged_cpm$simplified == "progenitor",]),
    cbind(y = "cTEC",      .[merged_cpm$simplified == "cTEC",]),
    cbind(y = "mTEC",      .[merged_cpm$simplified == "mTEC",]),
    cbind(y = "Tuft",      .[merged_cpm$simplified == "Tuft",]),
    cbind(y = "jTEC",      .[merged_cpm$simplified2 == "jTEC",]),
    cbind(y = "mTEClo",    .[merged_cpm$simplified2 == "mTEClo",]),
    cbind(y = "mTEChi",    .[merged_cpm$simplified2 == "mTEChi",]),
    cbind(y = "jTEC_inner",.[Idents(merged_cpm) == "jTEC (11wo)",]),
    cbind(y = "mTEClo_old",.[Idents(merged_cpm) == "mTEClo (11wo)",])
  )}
# thymus_cluster_lognormcount[1:10,1:7]
#                     y A830018L16Rik     Sulf1 Slco5a1      Eya1 Msc Sbspon
# foregut_1  progenitor     0.0000000 0.0000000       0 0.0000000   0      0
# foregut_2  progenitor     0.0000000 0.5642949       0 0.0000000   0      0
# foregut_3  progenitor     0.0000000 0.0000000       0 0.0000000   0      0

data.frame(y=thymus_cluster_cpm$y,
           rowSums=thymus_cluster_cpm[,-1] %>% rowSums) %>%
  boxplot(rowSums~y,data=., main="sum(CPM) of gene.use")

data.frame(y=thymus_cluster_lognormcount$y,
           rowSums=thymus_cluster_lognormcount[,-1] %>% rowSums) %>%
  boxplot(rowSums~y,data=., main="sum(logcount) of gene.use")


knn <- function(A,k) { # column-wisely similar columns were summed
  i = nrow(A)
  j = ncol(A)
  B = matrix(0, nrow = i, ncol = j)
  D = HiClimR::fastCor(A, nSplit = 1,optBLAS=T)
  for(m in 1:j){
    r = rank(-D[,m],ties.method = "first") <= k + 1
    B[,m] = rowSums(A[,r,drop = F])
  }
  colnames(B) <- colnames(A)
  rownames(B) <- rownames(A)
  B
}

thymus_cluster_scaled <- scale(log10(t(knn(t(thymus_cluster_cpm[,-1]),5))+0.1)) %>%
                                 {.[is.nan(.)]=0;.}
thymus_cluster_label <- thymus_cluster_cpm$y

# clean up-----------------------------------------------------------------------------------

thymus_cluster_scaled[1:10,1:10]
thymoma_mggene_scaled[1:10,1:10]
mg_genes <- intersect(colnames(thymus_cluster_scaled),colnames(thymoma_mggene_scaled))

thymus_cluster_scaled <- thymus_cluster_scaled[,mg_genes]
thymoma_mggene_scaled <- thymoma_mggene_scaled[,mg_genes]
dim(thymus_cluster_scaled)
dim(thymoma_mggene_scaled)

hg_genes <- structure(gene.use$hgname,names=gene.use$mgname)[mg_genes]

gene.use.matched <- cbind(mg_genes,hg_genes)

# save --------------------------------------------------------------------------------------
if(F){
  save(thymus_cluster_scaled,
       thymoma_mggene_scaled, 
       gene.use.matched,
       thymus_cluster_label, 
       file="data/data.scaled.for_comparison.RData")
  # save.image("data/data.scaled.for_comparison.extended.RData",compress = "gzip")
}

# end of script -----------------------------------------------------------------------------