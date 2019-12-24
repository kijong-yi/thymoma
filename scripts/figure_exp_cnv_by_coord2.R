library(tidyverse)
library(ComplexHeatmap)
library(ggsci)
library(scales)
library(circlize)
library(ggrepel)
library(wordcloud)
require(GSEABase)
library(fgsea)

load("data/figure_exp_cnv_by_coord.Rda")


parlayout3 <- function(nrow=3){
  par(mfrow = c(nrow,23),mar=c(0,0,0,0),oma=c(1,26,1.8,0),xpd=NA);
  mat=matrix(1:(23*nrow),nrow,byrow=T)
  chrsize = foreach(chr=chrlist,.combine=c) %do%{
    max(as.numeric(stringr::str_replace(rownames(plogsum[[chr]]),".*_","")))}
  layout(mat, widths = chrsize,
         heights = c(1,2), respect = FALSE)
}

parlayout3(2)
plot_sqz_cohort()
c("a"=1,"b"=5,"c"=4) %>%
{plot(x=.,y=1:length(.),yaxt='n',ylab="");axis(2,at=1:length(.),labels = names(.))}

gs = list("linoleic" = c("MIR4420","RN7SKP91","IRS4"),
          "devel" = c("D","E","F"))
gene.chr=="1"
gene.coord["IRS4"]
c5_list <- getGmt('~kjyi/ref/msigdb/c5.all.v6.2.symbols.gmt') %>% geneIds()
c2_list <- getGmt('~kjyi/ref/msigdb/c2.all.v6.2.symbols.gmt') %>% geneIds()
h_list <- getGmt('~kjyi/ref/msigdb/h.all.v6.2.symbols.gmt') %>% geneIds()
all_list <- append(c5_list,c2_list) %>% append(h_list)
names(all_list)[grepl("EMBRYONIC_ORGAN_DEVEL",names(all_list))]
gs <- all_list[c("GO_CORNIFIED_ENVELOPE",
                 "GO_KERATINIZATION",
                 "GO_KERATINOCYTE_DIFFERENTIATION",
                 "GO_STEROID_BIOSYNTHETIC_PROCESS",
                 "GO_CELLULAR_RESPIRATION",
                 "KEGG_LINOLEIC_ACID_METABOLISM",
                 "GO_DRUG_METABOLIC_PROCESS",
                 "GO_ANDROGEN_METABOLIC_PROCESS",
                 "GO_CORNIFIED_ENVELOPE",
                 "GO_TRICARBOXYLIC_ACID_METABOLIC_PROCESS",
                 "GO_EXTRACELLULAR_STRUCTURE_ORGANIZATION",
                 "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                 "GO_HEART_DEVELOPMENT",
                 "GO_EXTRACELLULAR_MATRIX",
                 "GO_TUBE_MORPHOGENESIS",
                 "GO_VASCULATURE_DEVELOPMENT",
                 "GO_EMBRYONIC_ORGAN_DEVELOPMENT")]

volc_dt <- read_tsv(paste0("~sypark/00_Project/01_thymoma/10_Final_data/01_expression/",
                           'IRS4_corrected_v2/thym_137s_tpm_geneSymbol_tumorOnly_IRS4cor.tsv')) %>%
  mutate_at(vars(-gene), list(~log10(.+0.1)))
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

set.seed(1)
fgseaRes <- fgseaMultilevel(pathways = gs, 
                            stats = stats$mw,
                            minSize=15,
                            maxSize=500,
                            nproc = 10)
leadingEdgeSubset <- fgseaRes$leadingEdge
names(leadingEdgeSubset) <- fgseaRes$pathway


plot_gene_coord <- function(gs, ...){
  for (chr in chrlist){
    # xmin=c(13000000)
    xmin=c(1)
    xmax=max(as.numeric(stringr::str_replace(rownames(plogsum[[chr]]),".*_","")))

    plot(0,xlim=c(xmin,xmax),ylim=c(1,length(gs)),
         xaxt='n',yaxt='n',bty="n",ylab="",
         col="#00000000")
    box(col="grey")
    if(chr=="1"){axis(side = 2,at = 1:length(gs),labels = names(gs),las=2,tick=F)}
    axis(1)
    for(i in 1:length(gs)){
      genes_in_this_chromosome = gs[[i]][gene.chr[gs[[i]]]==chr]
      their_coord = gene.coord[genes_in_this_chromosome]
      col_LES = ifelse(genes_in_this_chromosome %in% leadingEdgeSubset[[names(gs)[i]]],"red","#00000080")
      points(their_coord,rep(i,length(their_coord)),pch="|",col=col_LES)
    }
  }
}

parlayout3(3)
plot_sqz_cohort()
plot_gene_coord(gs)


# get gsea p values ......................
gtf=read_tsv(paste0("/home/users/kjyi/Projects/thymus_single_cell/final/",
                    "expression/IRS4/merged_gtf/Homo_sapiens.GRCh38.95.mod3.gtf"),skip=5,
             col_names = c(
               "chr","spircename","feature","start","end",
               "score","strand","frame","attribute"),
             col_types="ccciicccc") %>%
  dplyr::filter(feature == "gene") %>%
  dplyr::select(chr,start,end,attribute) %>%
  mutate(attribute = str_replace(attribute, '.*gene_name \"','') %>%
           str_replace('\".*','')) %>%
  as.data.frame()
gene.chr = structure(gtf$chr,names=gtf$attribute)
gene.coord = structure((gtf$start+gtf$end)/2,names=gtf$attribute)

h_list <- getGmt('~kjyi/ref/msigdb/h.all.v6.2.symbols.gmt') %>% geneIds()
c5_list <- getGmt('~kjyi/ref/msigdb/c5.all.v6.2.symbols.gmt') %>% geneIds()
c2_list <- getGmt('~kjyi/ref/msigdb/c2.all.v6.2.symbols.gmt') %>% geneIds()
all_list <- append(append(h_list,c2_list),c5_list);rm(h_list,c2_list,c5_list)
all_list <- all_list[str_replace(names(all_list),"_.*","") %in% c("GO","HALLMARK","KEGG")]

# 
ncasecoordampdel <- read_rds("data/ncasecoordampdel.Rds")
regions_amp_0.15 <- rownames(ncasecoordampdel)[ncasecoordampdel$w_amp > 0.15]
regions_amp_0.15.chr = str_replace(regions_amp_0.15,":.*","")
regions_amp_0.15.start = str_replace(regions_amp_0.15,".*:","")%>%str_replace("-.*","")  %>% as.numeric
regions_amp_0.15.end = str_replace(regions_amp_0.15,".*-","")  %>% as.numeric

check_overlap1 <- function(chr,coord,regions.chr, region.start, region.end){
  any(coord>regions.start[regions.chr==chr] & coord<regions.end[regions.chr==chr])
}

check_overlap2 <- function(genes,region.chr, region.start, region.end){
  lapply(genes,function(gene){
    chr=gene.chr[gene] %>% paste0("chr",.)
    coord=gene.coord[gene]
    any(coord>region.start[region.chr==chr] & coord<region.end[region.chr==chr])
  }) %>% unlist()
}

u = names(gene.chr)
y = check_overlap2(genes = u,
                   region.chr = regions_amp_0.15.chr,
                   region.start = regions_amp_0.15.start,
                   region.end = regions_amp_0.15.end) %>% factor
table(u %in% all_list[[1]])

fit <- fisher.test(factor(u %in% all_list[[1]]),factor(y))
fit$p.value
fit$estimate
{
  gs_overlap_amp_0.15.fit <- lapply(all_list,function(gs){
    x = factor(u %in% gs, levels=c(F,T))
    fisher.test(x,y)
  })
  gs_overlap_amp_0.15.p.value <- lapply(gs_overlap_amp_0.15.fit,function(x){x$p.value}) %>% unlist()
  gs_overlap_amp_0.15.oddsratio <- lapply(gs_overlap_amp_0.15.fit,function(x){x$estimate}) %>% unlist()
  names(gs_overlap_amp_0.15.oddsratio) = names(gs_overlap_amp_0.15.p.value)
}


x = factor(u %in% all_list[["GO_CORNIFIED_ENVELOPE"]], levels=c(F,T))
table(gs=x,ol=y)
y[x==T]

all_list[["GO_STEROID_BIOSYNTHETIC_PROCESS"]]
cbind(all_list[["GO_STEROID_BIOSYNTHETIC_PROCESS"]],
      gene.chr[all_list[["GO_STEROID_BIOSYNTHETIC_PROCESS"]]],
      gene.coord[all_list[["GO_STEROID_BIOSYNTHETIC_PROCESS"]]],
      check_overlap2(genes = all_list[["GO_STEROID_BIOSYNTHETIC_PROCESS"]],
                     region.chr = regions_amp_0.15.chr,
                     region.start = regions_amp_0.15.start,
                     region.end = regions_amp_0.15.end))


plot(gs_overlap_amp_0.15.p.value, gs_overlap_amp_0.15.oddsratio)



fgseaRes <- read_rds("data/fgseaRds_mw_all_go_hm_kegg.Rds")
fgseaRes.df <- as.data.frame(fgseaRes)


#factors
gs.df <- data.frame(ol.p = gs_overlap_amp_0.15.p.value[fgseaRes.df$pathway],
                    ol.odds = gs_overlap_amp_0.15.oddsratio[fgseaRes.df$pathway],
                    gsea.p = fgseaRes.df$pval,
                    gsea.nes = fgseaRes.df$NES)


smoothScatter(gs_overlap_amp_0.15.p.value[fgseaRes.df$pathway],fgseaRes.df$pval)
plot(gs_overlap_amp_0.15.p.value[fgseaRes.df$pathway] %>% log,fgseaRes.df$pval %>% log)

gs.df["KEGG_LINOLEIC_ACID_METABOLISM",]

gs.df[gs.df$gsea.nes<0 & gs.df$ol.odds>1,] %>%
  plot_ly(x = ~log10(ol.p), y = ~log10(gsea.p),
          text = rownames(.),
          hoverinfo = 'text')

# gs.df[gs.df$gsea.nes<0 & gs.df$ol.odds>1,] %>%
#   plot_ly(x = ~ol.odds, y = ~gsea.nes,
#         text = rownames(.),
#         hoverinfo = 'text')


mark <- c("GO_CORNIFIED_ENVELOPE",
  "GO_KERATINIZATION",
  "GO_KERATINOCYTE_DIFFERENTIATION",
  "GO_STEROID_BIOSYNTHETIC_PROCESS",
  "GO_CELLULAR_RESPIRATION") # %in% names(all_list)
label <- c("Cornified envelope","Keratinization", "Keratinocyte\ndifferentiation", "Steroid biosynthetic\nprocess", "Cellular respiration")

gs.df.dir <- gs.df[gs.df$gsea.nes<0 & gs.df$ol.odds>1,]
with(gs.df.dir,plot(x = log10(ol.p), y = log10(gsea.p)))

cairo_pdf("figures/scatter_pp_gsea_coord.pdf",height = 8/2.54,width=8/2.54,pointsize = 12*0.7)

par(mfrow=c(1,1),oma=c(0,0,0,0),mar=c(5,5.2,1,1))
with(gs.df.dir,plot(x = log10(ol.p), y = log10(gsea.p),xaxt='n',yaxt='n',xlab="",ylab=""))
axis(side = 2,at = 0:-9, labels = c(0,0.1, 0.01, expression(10^-3),expression(10^-4),
                                     expression(10^-5),expression(10^-6),expression(10^-7),
                                     expression(10^-8),expression(10^-9)),las=2)
axis(side = 1,at = c(0,-1,-2,-4,-6,-8), labels = c(0,0.1, 0.01, 
                                                   expression(10^-4),
                                                   expression(10^-6),
                                                   expression(10^-8)))
abline(v=log10(0.05),h=log10(0.05),col="grey",lty=2)
# with(gs.df.dir,plot(x = log10(ol.p), y = log10(gsea.p),xaxt='n',yaxt='n',xlab="",ylab=""))
with(gs.df.dir[mark,],
     textplot2(x = log10(ol.p),y = log10(gsea.p),words = label,new = F,
               xlim=c(-9,0),xadj = -0.8,yadj=-0.5,
               pt.pch=21,pt.col="black",pt.bg="red",pt.cex=1))
mtext(text = "Gene sets in recurrently amplified regions",side = 1,line=2.1)
mtext(text = expression("("*italic(p)*"-value of Fisher's exact test)"),side = 1,line=3.1)
mtext(text = expression("Enriched in "~"GTF2I"^"WT"~"(compare to "~"GTF2I"^"mut"*")"),side = 2,line=3.5)
mtext(text = expression("("*italic(p)*"-value of GSEA)"), side=2, line=2.5)
text(-9,log10(0.05),"0.05",col="grey",adj=c(0,1))
dev.off()
# save.image("~/Projects/thymus_single_cell/final2/data/snapshop_for_pp_scatter_gsea_overlapfisher.RData")

fgseaRes.df[fgseaRes.df$NES<0&fgseaRes.df$pval<0.01,1:6] %>% arrange(NES)
fgseaRes.df[fgseaRes.df$NES<0&fgseaRes.df$pval<0.01,1:6] %>% arrange(pval) %>% .$pathway %>% paste0('"',.,'"') %>%
  cat(sep=",")

gs_union <- list(union1=unname(unlist(all_list[c(
  # "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "KEGG_LINOLEIC_ACID_METABOLISM",
  # "GO_CELLULAR_RESPIRATION",
  # "GO_DNA_PACKAGING_COMPLEX",
  # "GO_OXIDATIVE_PHOSPHORYLATION",
  "GO_DRUG_METABOLIC_PROCESS",
  "GO_RESPONSE_TO_XENOBIOTIC_STIMULUS",
  "GO_STEROID_METABOLIC_PROCESS",
  "GO_ANDROGEN_METABOLIC_PROCESS",
  "GO_RESPONSE_TO_GROWTH_HORMONE",
  "KEGG_STEROID_HORMONE_BIOSYNTHESIS",
  "GO_CORNIFIED_ENVELOPE",
  # "GO_ORGANELLE_INNER_MEMBRANE",
  "GO_STEROID_BIOSYNTHETIC_PROCESS",
  # "GO_INTERMEDIATE_FILAMENT",
  # "GO_RESPIRATORY_CHAIN",
  # "GO_OXYGEN_BINDING",
  # "GO_AEROBIC_RESPIRATION",
  "GO_CELLULAR_HORMONE_METABOLIC_PROCESS",
  "KEGG_METABOLISM_OF_XENOBIOTICS_BY_CYTOCHROME_P450",
  # "GO_ELECTRON_TRANSPORT_CHAIN",
  # "GO_INTERMEDIATE_FILAMENT_CYTOSKELETON",
  "GO_TRICARBOXYLIC_ACID_METABOLIC_PROCESS",
  "KEGG_STARCH_AND_SUCROSE_METABOLISM",
  # "GO_CHROMATIN_SILENCING",
  "GO_KERATINIZATION"
  # "GO_MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX_I_BIOGENESIS"
)])),
union2=unname(unlist(all_list[c(
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  # "KEGG_LINOLEIC_ACID_METABOLISM",
  "GO_CELLULAR_RESPIRATION",
  # "GO_DNA_PACKAGING_COMPLEX",
  "GO_OXIDATIVE_PHOSPHORYLATION",
  # "GO_DRUG_METABOLIC_PROCESS",
  # "GO_RESPONSE_TO_XENOBIOTIC_STIMULUS",
  "GO_STEROID_METABOLIC_PROCESS",
  # "GO_ANDROGEN_METABOLIC_PROCESS",
  # "GO_RESPONSE_TO_GROWTH_HORMONE",
  # "KEGG_STEROID_HORMONE_BIOSYNTHESIS",
  # "GO_CORNIFIED_ENVELOPE",
  # "GO_ORGANELLE_INNER_MEMBRANE",
  "GO_STEROID_BIOSYNTHETIC_PROCESS",
  # "GO_INTERMEDIATE_FILAMENT",
  # "GO_RESPIRATORY_CHAIN",
  # "GO_OXYGEN_BINDING",
  # "GO_AEROBIC_RESPIRATION",
  # "GO_CELLULAR_HORMONE_METABOLIC_PROCESS",
  # "KEGG_METABOLISM_OF_XENOBIOTICS_BY_CYTOCHROME_P450",
  "GO_ELECTRON_TRANSPORT_CHAIN"
  # "GO_INTERMEDIATE_FILAMENT_CYTOSKELETON",
  # "GO_TRICARBOXYLIC_ACID_METABOLIC_PROCESS",
  # "KEGG_STARCH_AND_SUCROSE_METABOLISM"
  # "GO_CHROMATIN_SILENCING",
  # "GO_KERATINIZATION"
  # "GO_MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX_I_BIOGENESIS"
)])),
union3=unname(unlist(all_list[c(
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "KEGG_LINOLEIC_ACID_METABOLISM",
  "GO_CELLULAR_RESPIRATION",
  "GO_DNA_PACKAGING_COMPLEX",
  "GO_OXIDATIVE_PHOSPHORYLATION",
  "GO_DRUG_METABOLIC_PROCESS",
  "GO_RESPONSE_TO_XENOBIOTIC_STIMULUS",
  "GO_STEROID_METABOLIC_PROCESS",
  "GO_ANDROGEN_METABOLIC_PROCESS",
  "GO_RESPONSE_TO_GROWTH_HORMONE",
  "KEGG_STEROID_HORMONE_BIOSYNTHESIS",
  "GO_CORNIFIED_ENVELOPE",
  "GO_ORGANELLE_INNER_MEMBRANE",
  "GO_STEROID_BIOSYNTHETIC_PROCESS",
  "GO_INTERMEDIATE_FILAMENT",
  "GO_RESPIRATORY_CHAIN",
  "GO_OXYGEN_BINDING",
  "GO_AEROBIC_RESPIRATION",
  "GO_CELLULAR_HORMONE_METABOLIC_PROCESS",
  "KEGG_METABOLISM_OF_XENOBIOTICS_BY_CYTOCHROME_P450",
  "GO_ELECTRON_TRANSPORT_CHAIN",
  "GO_INTERMEDIATE_FILAMENT_CYTOSKELETON",
  "GO_TRICARBOXYLIC_ACID_METABOLIC_PROCESS",
  "KEGG_STARCH_AND_SUCROSE_METABOLISM",
  "GO_CHROMATIN_SILENCING",
  "GO_KERATINIZATION",
  "GO_MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX_I_BIOGENESIS"
)])))

table(gene.set=factor(u %in% gs_union$union3, levels=c(F,T)),coord.overlap=y)
fisher.test(x=factor(u %in% gs_union$union3, levels=c(F,T)),y=y)

table(gene.set=factor(u %in% gs_union$union1, levels=c(F,T)),coord.overlap=y)
fisher.test(x=factor(u %in% gs_union$union1, levels=c(F,T)),y=y)

table(gene.set=factor(u %in% gs_union$union2, levels=c(F,T)),coord.overlap=y)
fisher.test(x=factor(u %in% gs_union$union2, levels=c(F,T)),y=y)

meta_dt$IRS4_tCN
meta_dt$corrected_IRS4_TPM
meta_dt$final_cellularity
meta_dt$final_ploidy
meta_dt$GTF2I_status2
meta_dt$histologic_type
meta_dt[meta_dt$GTF2I_status2=="w",] %>% boxplot(log10(corrected_IRS4_TPM+1)~pmin(5,IRS4_tCN),data=.)




meta_dt <- read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/02_metadata/Thymoma_summary_191204_1stSheet.txt')
cairo_pdf("figures/boxplot_xpcopy_irs4tpm.pdf",height = 7/2.54,width=6/2.54,pointsize = 12*0.7)
par(mfrow=c(1,1),oma=c(0,0,0,0),mar=c(5,4,1,1))
meta_dt[meta_dt$GTF2I_status2=="w",] %>% 
  boxplot(log10(corrected_IRS4_TPM+1)~pmin(5,IRS4_tCN),data=.,xlab="",ylab="",xaxt="n",yaxt="n",
          col=c("#595045","#A69886","#BF2626","#8C1B1B","#400101"))
axis(side = 1,at = 1:5,labels = c(1:4,"≥ 5"),line=0)
axis(side = 1,at = c(2,4,5), labels = c(1,2,"≥ 3"),tick=F,line = 1)
mtext(text = "Female",side = 1,line = 1,at = 0)
mtext(text = "Male",side = 1,line = 2,at = 0)
mtext(text = "Copy number (Xp)",side = 1,line = 3)
mtext(text = "IRS4 expression (TPM)",side = 2,line=2.5)
axis(side = 2,at = log10(c(0,10,50,100,200,400)+1), labels = c(0,10,50,100,200,400),las=2)
dev.off()
