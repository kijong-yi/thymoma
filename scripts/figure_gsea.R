library(tidyverse)
library(gridExtra)
library(fgsea)
require(GSVA)
require(GSEABase)
library(grid)

table(meta_dt$IRS4_tCN>2, meta_dt$GTF2I_status2)
table(meta_dt$corrected_IRS4_TPM>5, meta_dt$GTF2I_status2)

# prep ----------------------------------------------------------------------------------------------------------------
meta_dt <- read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/02_metadata/Thymoma_summary_191115_1stSheet.txt')
volc_dt <- read_tsv(paste0("~sypark/00_Project/01_thymoma/10_Final_data/01_expression/",
                           'IRS4_corrected_v2/thym_137s_tpm_geneSymbol_tumorOnly_IRS4cor.tsv')) %>%
  mutate_at(vars(-gene), list(~log10(.+0.01)))

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

all_list<-all_list[str_replace(names(all_list),"_.*","") %in% c("GO","HALLMARK","KEGG")]

# calculate -----------------------------------------------------------------------------------------------------------
set.seed(1)
fgseaRes <- fgseaMultilevel(pathways = all_list, 
                            stats = stats$mw,
                            minSize=15,
                            maxSize=500,
                            nproc = 10)
set.seed(1)
collapsedPathways_all_mw <- collapsePathways(fgseaRes[order(pval)][padj < 0.05], all_list, stats$mw)

topPathwaysUp <- fgseaRes[pathway %in% collapsedPathways_all_mw$mainPathways,][ES > 0][head(order(pval), n=10),][order(-NES), pathway]
topPathwaysDown <- fgseaRes[pathway %in% collapsedPathways_all_mw$mainPathways,][ES < 0][head(order(pval), n=10),][order(NES), pathway]
(topPathways <- c(topPathwaysUp, rev(topPathwaysDown)))

resmat <- fgseaRes[match(topPathways,pathway),c("pathway","padj","NES")]
top_pathway_info <-fgseaRes[match(topPathways,pathway),] # for save as table

try(dev.off())

# draw -----------------------------------------------------------------------------------------------------------------

cairo_pdf("figures/gsea.1.pdf",height = 10/2.54,width=12/2.54,pointsize = 12*0.7)
par(mar=c(5,22,3,1))
with(resmat,{
  firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  }
  plot(-structure(NES,names=pathway),structure(1:length(pathway),names=pathway),
       las=1,yaxt="n",ylab="",type='n',xlab="normalized enrichment score\nGTF2I-mutant ↔ wild-type",
       main="Top enriched gene sets");
  segments(x0 = -NES,y0 = 1:length(pathway),x1 = 0,lty = 1,lwd = 5,col = "grey50",lend=1)
  points(-structure(NES,names=pathway),structure(1:length(pathway),names=pathway),
         bg=circlize::colorRamp2(c(10:0),RColorBrewer::brewer.pal(11,"Spectral"))(-log10(padj)),
         pch=21,cex=2)
  axis(2,at = 1:length(pathway),labels = pathway %>% tolower %>%
         str_replace("[^_]*_","") %>%
         str_replace("rna","RNA") %>%
         str_replace("dna","DNA") %>%
         str_replace("_i_","_I_") %>%
         str_replace("Calcium","Ca") %>%
         str_replace("dependent","dep.") %>%
         str_replace("_iii_","_III_") %>%
         str_replace_all("_"," ") %>%
         firstup ,las=2)
  abline(v=0,lty=2,col="grey")
})
text(1.5,20,expression(-log[10]~q))
par(new=T,fig=c(.94,.96,.2,.35),mar=c(0,0,0,0))
plot(100,xlim=c(0,1),ylim=c(0,20),xaxt='n',bty='n',xlab="",las=2)


rasterImage(as.raster(RColorBrewer::brewer.pal(11,"Spectral"), ncol=1),
            0,0,1,20)
rect(0,0,1,20)
dev.off()


# save as table --------------------------------------------------------------------------------------------------------


topPathwaysUp2 <- fgseaRes[pathway %in% collapsedPathways_all_mw$mainPathways,][ES > 0][head(order(pval), n=100),][order(-NES), pathway]
topPathwaysDown2 <- fgseaRes[pathway %in% collapsedPathways_all_mw$mainPathways,][ES < 0][head(order(pval), n=100),][order(NES), pathway]
(topPathways2 <- c(topPathwaysUp2, rev(topPathwaysDown2)))


top_pathway_info <-fgseaRes[match(topPathways2,pathway),] # for save as table

x <- top_pathway_info %>% as.data.frame %>% as.data.frame %>% tbl_df()

x$leadingEdge = lapply(x$leadingEdge,function(x){paste(x,collapse=",")}) %>% unlist()

x$pathwaygroup=str_replace(x$pathway,"_.*","")

write_csv(x,"tables/gsea_top_pathway_info_mt_vs_wt.2.csv")

write_rds(fgseaRes,"data/fgseaRds_mw_all_go_hm_kegg.Rds")


# end of script --------------------------------------------------------------------------------------------------------
