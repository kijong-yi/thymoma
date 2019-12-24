
library(statmod)
library(Seurat)
library(tidyverse)
library(kjyi, lib.loc= "/home/users/kjyi/R/x86_64-redhat-linux-gnu-library/3.6")

textplot2 <- function(x, y, words, cex = 1, new = TRUE, show.lines = TRUE, xadj=0,yadj=0, ...){
  if (new) 
    plot(x, y, type = "n", ...)
  lay <- wordcloud::wordlayout(c(x+xadj,x), c(y+yadj,y), rep(words,2), cex, ...)
  if (show.lines) {
    for (i in 1:length(x)) {
      xl <- lay[i, 1]
      yl <- lay[i, 2]
      w <- lay[i, 3]
      h <- lay[i, 4]
      if (x[i] < xl || x[i] > xl + w || y[i] < yl || y[i] > 
          yl + h) {
        points(x[i], y[i], pch = 16, col = "red", cex = 0.5)
        nx <- xl + 0.5 * w
        ny <- yl + 0.5 * h
        lines(c(x[i], nx), c(y[i], ny), col = "grey")
      }
    }
  }
  lay <- lay[1:length(x),]
  text(lay[, 1] + 0.5 * lay[, 3], lay[, 2] + 0.5 * lay[, 4], 
       words, cex = cex, ...)
}


# gene to use in nearest k-neighbor mapping
merged <- read_rds("data/singlecell/merged.Rds")
# intersect genes that are not expressed in stroma and variably expressed in thymoma 

bm <- kjyi::load_biomart('<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
             <Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
             
             <Dataset name = "mmusculus_gene_ensembl" interface = "default" >
             <Attribute name = "ensembl_gene_id" />
             <Attribute name = "external_gene_name" />
             </Dataset>
             
             <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
             <Attribute name = "ensembl_gene_id" />
             <Attribute name = "external_gene_name" />
             </Dataset>
             </Query>',col_names=c("mgid","mgname","hgid","hgname"))

thymoma_tpm <- read_tsv(paste0("~sypark/00_Project/01_thymoma/10_Final_data/01_expression/",
                          'IRS4_corrected_v2/',
                          'thym_137s_tpm_geneSymbol_tumorOnly_IRS4cor.tsv')) %>% 
  as.data.frame() %>%
  column_to_rownames("gene")

# thymoma_tpm <- read_tsv("~kjyi/Projects/thymus_single_cell/final/expression/tpm_geneSymbol_tumorOnly.tsv", 
#                         col_types = cols(.default = "d", gene = "c")) %>%
#   as.data.frame() %>%
#   column_to_rownames("gene")
thymoma_tpm[1:3,1:3]

# plot(unlist(thymoma_tpm["SPOCK1",]), unlist(thymoma_tpm["IRS4",]))
# pct of count in stroma, total UMI count
gene_table <- data.frame(mgname = rownames(merged)) %>% tbl_df
mat <- GetAssayData(subset(merged), slot = "counts")
gene_table$stromal_occupacy <- Matrix::rowSums(mat[,merged$simplified == "stroma"])/Matrix::rowSums(mat)
gene_table$stromal_occupacy[is.nan(gene_table$stromal_occupacy)] <- 0
gene_table$total_count_in_single_cell <- Matrix::rowSums(mat)
gene_table <- gene_table %>% left_join(bm) %>% na.omit() %>%
  left_join(data.frame(hgname = rownames(thymoma_tpm), 
                       thymoma_mean_tpm = apply(thymoma_tpm, 1, mean),
                       thymoma_var_tpm = apply(thymoma_tpm, 1, var),
                       stringsAsFactors=F) %>%
              mutate(thymoma_cv2 = thymoma_var_tpm/thymoma_mean_tpm^2)) %>%
  na.omit()

pdf("figures/supplimentary.gene.selection.nonstromal.highvarinthymoma.pdf", 13,7)
par(mfrow = c(1,2))
with(gene_table, smoothScatter(stromal_occupacy, log(total_count_in_single_cell), 
                               xlab = "Stromal transcript occupacy", ylab="log(total UMI count)"))
title("Selection of epithelial genes")
rect(xleft = -0.037, ybottom = log(10), xright = 0.15, ytop = 15.41 ,border = "#00FF00")
# gene_table %>% filter(hgname %in% c("AIRE", "PSMB11", "LY75", "LGALS7", "DMKN", "KRT15", "KRT14", "TBX1",
#                                     "CD40","CLDN4", "PAX9", "PAX1", "FOXN1")) %>%
#   with(points(stromal_occupacy, log(total_count_in_single_cell), col = "#FF0000"))


gene_table %>% filter(hgname %in% c("AIRE", "PSMB11", "LY75", "LGALS7", "DMKN", "KRT15", "KRT14", "TBX1",
                                    "CD40","CLDN4", "PAX9", "PAX1","FOXN1")) %>%
  with(textplot2(stromal_occupacy,log(total_count_in_single_cell), c("AIRE", "PSMB11", "LY75", "LGALS7", "DMKN", "KRT15", "KRT14", "TBX1",
                                                                     "CD40","CLDN4", "PAX9", "PAX1","FOXN1"),new=F,
                 yadj=(log(total_count_in_single_cell)-9)*0.2,
                 xadj=0.15,
                 xlim=c(0.1,15),ylim=c(4,15)))



# gene_table %>% filter(total_count_in_single_cell > 10, stromal_occupacy < 0.15) %>%
#   with(points(log(thymoma_mean_tpm),log(thymoma_cv2), col = "#00FF0030"))
non_stromal_genes <- with(gene_table, total_count_in_single_cell > 2 & stromal_occupacy < 0.15)

with(gene_table, smoothScatter(log(thymoma_mean_tpm),log(thymoma_cv2),
                               xlab = "Mean expression (thymoma)\n(logTPM)",
                               ylab = ""))
mtext(expression(log(CV^2)),2,2)


# gene_table$total_count_in_single_cell[9637]
# gene_table$stromal_occupacy[9637]

# variable non-stromal genes in thymoma dataset
params=list(a=.2,
            b=.6)
abline(h = log(params$a), col = "#00000030")
# fit a regression line
minMeanForFit <- unname(quantile(gene_table$thymoma_mean_tpm[which(gene_table$thymoma_cv2>params$a)], 
                                 params$b) )
abline(v = log(minMeanForFit), col = "#00000030")
useForFit <- gene_table$thymoma_mean_tpm >= minMeanForFit
# useForFit <- useForFit & non_stromal_genes
# gene_table[useForFit,] %>%
#   with(points(log(thymoma_mean_tpm),log(thymoma_cv2), col = "#00FF0030"))


fit <- glmgam.fit(cbind( a0 = 2, a1tilde = 1/gene_table$thymoma_mean_tpm[useForFit]), 
                  gene_table$thymoma_cv2[useForFit] )
a0 <- unname( fit$coefficients["a0"] )
a1 <- unname( fit$coefficients["a1tilde"])
fit$coefficients
# Now add the fit and the 95% confidence interval to our plot:
xg <- exp(seq( min(log(gene_table$thymoma_mean_tpm[gene_table$thymoma_mean_tpm>0])), 
               max(log(gene_table$thymoma_mean_tpm)), length.out=1000 ))
vfit <- a1/xg + a0
# add fit line
lines( log(xg), log(vfit), col="black", lwd=3 )
df <- ncol(thymoma_tpm) - 1
# add confidence interval
# lines(log(xg),log(vfit * qchisq(0.975,df)/df),lty=2,col="black")
# lines(log(xg),log(vfit * qchisq(0.025,df)/df),lty=2,col="black")

# Rank genes by the significance of deviation from the fit
afit <- a1/gene_table$thymoma_mean_tpm+a0
varFitRatio <- gene_table$thymoma_var_tpm/(afit*gene_table$thymoma_mean_tpm^2)
varorder <- order(varFitRatio,decreasing=T)
over_fit_line <- varFitRatio > 1
over_lower_fit_line <- intersect(varorder[1:5000], which(non_stromal_genes))

points(log(gene_table$thymoma_mean_tpm[over_fit_line&non_stromal_genes]),
       log(gene_table$thymoma_cv2[over_fit_line&non_stromal_genes]),
       col="#FFEB3380")
# points(log(gene_table$thymoma_mean_tpm[over_lower_fit_line]),log(gene_table$thymoma_cv2[over_lower_fit_line]),col="#FFEB3380")
# monitoring markers
# gene_table %>% filter(hgname %in% c("AIRE", "PSMB11", "LY75", "LGALS7", "DMKN", "KRT15", "KRT14", "TBX1",
#                                     "CD40","CLDN4", "PAX9", "PAX1","FOXN1")) %>%
#   with(points(log(thymoma_mean_tpm),log(thymoma_cv2), col = "#FF0000"))

gene_table %>% filter(hgname %in% c("AIRE", "PSMB11", "LY75", "LGALS7", "DMKN", "KRT15", "KRT14", "TBX1",
                                    "CD40","CLDN4", "PAX9", "PAX1","FOXN1","UGT1A8")) %>%
  with(textplot2(log(thymoma_mean_tpm),log(thymoma_cv2), c("AIRE", "PSMB11", "LY75", "LGALS7", "DMKN", "KRT15", "KRT14", "TBX1",
                                                           "CD40","CLDN4", "PAX9", "PAX1","FOXN1","UGT1A8"),new=F,
                 yadj=(log(thymoma_cv2)+2)*0.1+0.3,
                 xadj=(log(thymoma_mean_tpm)+4)*0.1+2,
                 xlim=c(-10,8),ylim=c(-2,5)))

cbind(gene_table,non_stromal_genes)[gene_table$hgname=="UGT1A8",]
cbind(gene_table,non_stromal_genes)[gene_table$mgname %in% c("Ugt1a8","Ugt1a10","Ugt1a9"),]

gene.use

lines( log(xg), log(vfit), col="black", lwd=3 )
lines(log(xg),log(vfit * qchisq(0.995,df)/df),lty=2,col="black")
lines(log(xg),log(vfit * qchisq(0.005,df)/df),lty=2,col="black")
legend("topright", col = "#FFEB33", pch = 1, legend = gene_table[over_fit_line&non_stromal_genes,] %>% .$mgname %>% unique %>% length)
title("Selection of high variance gene (thymoma)")
dev.off()

which(gene_table$hgname == "LGALS7")
which(gene_table$hgname == "AIRE")
gene_table[9637,]
non_stromal_genes[9637]
non_stromal_genes[11359]
# final gene.use
gene_table %>% filter(over_fit_line & non_stromal_genes)
gene_table %>% filter(hgname=="FOXN1")

gene.use <- gene_table[over_fit_line & non_stromal_genes,] %>% dplyr::select(mgname,hgname)
c("AIRE", "PSMB11", "LY75", "LGALS7", "DMKN", "KRT15", "KRT14", "TBX1",
  "CD40","CLDN4", "PAX9", "PAX1","FOXN1")  %>% .[!.%in% gene.use$hgname]
gene_table$hgname

if(F){
  write_rds(gene.use, "data/gene.use.Rds")
  rm(mat,a0,a1,afit,df,minMeanForFit,non_stromal_genes,over_fit_line,useForFit,varFitRatio,varorder,vfit,xg, fit, gene_table, over_lower_fit_line,params,g)
  rm(thymoma_tpm,bm)
}
