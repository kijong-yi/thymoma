library(tidyverse)
library(ComplexHeatmap, lib.loc = "/home/users/kjyi/R/x86_64-redhat-linux-gnu-library/3.6")
library(ggsci, lib.loc = "/home/users/kjyi/R/x86_64-redhat-linux-gnu-library/3.6")
library(circlize, lib.loc = "/home/users/kjyi/R/x86_64-redhat-linux-gnu-library/3.6")
library(ggrepel, lib.loc = "/home/users/kjyi/R/x86_64-redhat-linux-gnu-library/3.6")
require(GSEABase, lib.loc = "/home/users/kjyi/R/x86_64-redhat-linux-gnu-library/3.6")
library(seriation, lib.loc = "/home/users/kjyi/R/x86_64-redhat-linux-gnu-library/3.6")

# sypark's code --------------------------------------------------------------------------------------------------------
meta_dt <- read_tsv(paste0('~sypark/00_Project/01_thymoma/10_Final_data/02_metadata/Thymoma_summary_191205_1stSheet.txt'))
MT_ids <- meta_dt$id[meta_dt$GTF2I_status2 == 'm']
WT_ids <- meta_dt$id[meta_dt$GTF2I_status2 == 'w']
CA_ids <- meta_dt$id[meta_dt$GTF2I_status2 == 'c']
bm_dt <- read_tsv(paste0('~sypark/02_Reference/13_biomart/',
                         'biomart_human_geneID_transcriptID_hgncsymbol',
                         '_genetype_ensemblgenename_190522.txt')) %>% 
  dplyr::rename(gene = `Gene name`, gene_type = `Gene type`) %>% dplyr::select(gene, gene_type) %>% unique()
# load expression data
exp_dt <- read_tsv(paste0("~sypark/00_Project/01_thymoma/10_Final_data/01_expression/",
                          'IRS4_corrected_v2/',
                          'thym_137s_tpm_geneSymbol_tumorOnly_IRS4cor.tsv'))
# exp_dt_pcg <- left_join(exp_dt, bm_dt, by="gene") %>%
#   filter(gene_type == 'protein_coding') %>% select(-gene_type) # filter protein coding only
l10_exp_dt <- exp_dt %>% mutate_at(vars(-gene), list(~log10(.+0.01)))
l10_exp_dt <- left_join(l10_exp_dt, bm_dt, by="gene") %>%
  dplyr::filter(gene_type == 'protein_coding') %>% dplyr::select(-gene_type) # filter protein coding only
my_pal = pal_npg("nrc")(10)
meta_dt$histologic_type %>% unique()
histo_pal = pal_aaas("default")(10)[c(4,1,7,8,6,2,3,5)]
names(histo_pal) = c("A","AB","MN-T","B1","B2","B3","NE","TC")
stage_pal=c('#ffffd4','#fed98e','#fe9929','#d95f0e','#993404', 'white')
names(stage_pal) <- c('I','II','III','IVa','IVb', '0')
cohort_pal = pal_jama("default", alpha = 0.8)(7)[c(1,7)]
names(cohort_pal) = c("SNUH","TCGA_CancerCell")
gtf2i_pal = pal_npg("nrc")(10)[c(4,1,3)]
names(gtf2i_pal) = c('m','w','c')
exp_pal = colorRamp2(c(-3,0,3), c('#253494',"gray90",'#f03b20'))
exp_pal2 = colorRamp2(c(-1,0,1), c('#253494',"gray90",'#f03b20'))
purity_pal = colorRamp2(c(0,0.3,1), c("red","#272822","#999488"))
myred = pal_npg("nrc")(10)[8]
myblue = pal_npg("nrc")(10)[4]
# unsupervised clustering
hm_dt <- l10_exp_dt
hm_dt <- hm_dt %>% as.data.frame() %>% column_to_rownames('gene') 
hm_dt$var <- apply(hm_dt, 1, function(x) var(x))
hm_dt <- hm_dt %>% rownames_to_column('gene')
var_genes <- hm_dt %>% arrange(desc(var))%>%.$gene%>%.[1:2500] # top 2500 high variance genes, variance in log10 scale
hm_dt <- hm_dt %>% dplyr::filter(gene %in% var_genes) %>% dplyr::select(-var) %>% column_to_rownames('gene') %>%
  as.matrix()
hm_annot_dt <- meta_dt %>% dplyr::select(id, cohort, GTF2I_status2, histologic_type, Stage, ImmuneScore,Purity=final_cellularity) %>% 
  as.data.frame() %>% column_to_rownames('id')
hm_annot_dt <- hm_annot_dt[colnames(hm_dt),]
cohort_pal2 = cohort_pal
names(cohort_pal2)[2] = "TCGA"



hm_dt <- t(scale(t(hm_dt))) # row scaling
dim(hm_dt)
# Heatmap(hm_dt, clustering_distance_columns = 'pearson',
#         clustering_method_columns =  "average",
#         clustering_method_rows = "ward.D",
#         top_annotation = top_anno, col= exp_pal, show_row_names = F)



# gene set prep --------------------------------------------------------------------------------------------------------
# h_list <- getGmt('~kjyi/ref/msigdb/h.all.v6.2.symbols.gmt') %>% geneIds()
# c2_list <- getGmt('~sypark/00_Project/01_thymoma/10_Final_data/04_GSVA/genesets/c2.cp.v6.2.symbols.gmt') %>% geneIds()
# c5_list <- getGmt('~kjyi/ref/msigdb/c5.all.v6.2.symbols.gmt') %>% geneIds()
# all_list <- append(append(h_list,c2_list),c5_list);rm(h_list,c2_list,c5_list)


# names(all_list)[grepl("DEVELOP",names(all_list))&grepl("THYM",names(all_list))]
# 
# names(all_list)[grepl("GO",names(all_list))&grepl("_T_CELL",names(all_list))]

all_list <- list("T cell development" = c("CD3D","CD3E","GPAP2","LCK"),
                 "Metabolic process"=c("CYP2C9","CYP3A5","CYP3A4","STAR","UGT2B7","PPARGC1A","AKR1D1","SRD5A1","WNT4","SLC34A1","ASS1"),
                 "ECM organization" = c("ADAMTS20","COL9A3","MMP3","THSD4","COL11A1"),
                 "TGFβ signaling"=c("GDF5","GDF1","BMP2","BMP4","LEFTY2","BMP8B","LRG1","SMAD7","GDF6","BMP6","TGFB2"),
                 "Development"=c("TBX1","MYH6","WNT2","HMGA2","WNT5A","BMP4","CX3CR1","IRX2","IRX4","SALL1"),
                 "TNFα signaling/inflammation"=c("CCL20","CXCL13"))
selected_names <- names(all_list)
#
# selected_names <- c(
#   # "GO_OXIDATIVE_PHOSPHORYLATION",
#   # "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
#   # "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
#   # "GO_THYMUS_DEVELOPMENT",
#   # "GO_EMBRYONIC_ORGAN_MORPHOGENESIS",
#   "GO_KERATINIZATION",
#   "GO_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
#   # "GO_CHROMATIN_SILENCING", #hist...
#   # "REACTOME_TCR_SIGNALING", # good
#   # "GO_T_CELL_DIFFERENTIATION",
#   # "KEGG_LINOLEIC_ACID_METABOLISM", 
#   # "KEGG_STEROID_HORMONE_BIOSYNTHESIS", # good??
#   "KEGG_TGF_BETA_SIGNALING_PATHWAY")#, # good
# # "HALLMARK_INTERFERON_GAMMA_RESPONSE")#,
# # "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION") #good?)




#
# library(piano)
# gslist2gsc <- function(gslist){
#   lapply(names(gslist),function(gsname){cbind(gslist[[gsname]],gsname)}) %>% do.call(rbind,.) %>% loadGSC}
# gsares <- runGSAhyper(genes=all_names,
#                       gsc = gslist2gsc(all_list))
if(F){
  gsares$pvalues %>% sort %>% names %>% head(100) %>% paste0('"',.,'"') %>% cat(sep=",\n")
  
  gsares$pvalues %>% sort %>% names %>% head(200) %>% .[grepl("EXTRACELLULAR",.)] -> x
  gsares$pvalues %>% sort %>% names %>% head(200) %>% .[grepl("EMBRY",.)] -> x
  structure(lapply(x,function(x){table(all_list[[x]]%in%rownames(hm_dt))}),names=x)
}

# o1 = seriate(dist(hm_dt), method = "TSP") # about genes
o1p = seriate(as.dist(1-cor(t(hm_dt))), method = "TSP") # about genes

# o1pgw = seriate(as.dist(1-cor(t(hm_dt))), method = "GW") # about genes
# hm_dt <- hm_dt[nrow(hm_dt):1,]


length(selected_names)
myninecolor=RColorBrewer::brewer.pal(9,"Set1")[c(1,5,2,4,3,6)]
myninecolor=c("#E41A1C",
              "#FF7F00",
              "#377EB8",
              "#984EA3",
              "#4DAF4A",
              "#008280")
all_names = rownames(hm_dt)
all_color=rep("black",length(all_names))
for(i in 1:length(selected_names)){
  all_color[all_names %in% all_list[[selected_names[i]]]] = myninecolor[i]
}
table(all_color)


o1pgw = hclust(as.dist(1-cor(t(hm_dt)))) # about genes

o2pgw = seriate(as.dist(1-cor(hm_dt)), method = "GW")

# o2 = seriate(dist(t(hm_dt)), method = "TSP") # about samples 
# o2p = seriate(as.dist(1-cor(hm_dt)), method = "TSP")

top_anno <- HeatmapAnnotation("Group"= hm_annot_dt$GTF2I_status2,
                              "Histologic type"= hm_annot_dt$histologic_type,
                              "Purity" = hm_annot_dt$Purity,
                              # "Cohort" = hm_annot_dt$cohort %>% str_replace("_CancerCell",""),
                              col = list("Group" = gtf2i_pal,
                                         "Histologic type" = histo_pal,
                                         "Purity" = purity_pal),
                              # gp = gpar(cex=0.75),
                              annotation_legend_param = list(
                                "Histologic type" = list(
                                  title="Histologic type",
                                  at = c("A","AB","MN-T","B1","B2","B3","NE","TC"),
                                  labels = c("Type A thymoma",
                                             "Type AB thymoma",
                                             "Micronodular thymoma with lymphoid stroma",
                                             "Type B1 thymoma",
                                             "Type B2 thymoma",
                                             "Type B3 thymoma",
                                             "Neuroendocrine carcinoma",
                                             "Thymic carcinoma (Squamous cell carcinoma,Undifferentiated ca)")
                                ),
                                "Group" = list(
                                  title = "Group",
                                  at = c("m", "w","c"),
                                  labels = c("GTF2I-mutant", "Wild-type","Thymic carcinoma"))))
ha = rowAnnotation(foo = anno_mark(at = which(all_color != "black"),
                                   labels = all_names[all_color != "black"], 
                                   lines_gp = gpar(col=all_color[which(all_color != "black")]),
                                   labels_gp = gpar(col=all_color[which(all_color != "black")],cex=1)))
# lgd = Legend(labels = c("Keratinization", "TCR signaling", "TGFβ signaling"), title = "Gene ontology", legend_gp = gpar(fill = myninecolor))
lgd = Legend(labels = selected_names, title = "Gene ontology", legend_gp = gpar(col = myninecolor,bg="white"),
             labels_gp = gpar(col = myninecolor),
             type = "lines")

lgd2 = Legend(col_fun = exp_pal, title = "Normalized gene expression")

lgd12=packLegend(lgd,lgd2)
draw(lgd12)

h1 <- Heatmap(hm_dt,
              row_order = get_order(o1p),width = unit(9,"cm"),height = unit(11,"cm"),
              show_row_dend = F,
              show_column_names = F,
              cluster_columns = as.dendrogram(o2pgw[[1]]),
              cluster_rows = rev(as.dendrogram(o1pgw)),
              top_annotation = top_anno, col= exp_pal, show_row_names = F,
              right_annotation = ha,show_heatmap_legend = F,
              row_split = 4,
              left_annotation = rowAnnotation( width = unit(4.5, "mm"),d=anno_block(
                gp = gpar(fill = c("#FF7F00","#008280","#EE0000","#377EB8")),
                labels = c("group1", "group2", "group3","group4"),
                labels_gp = gpar(col = "white", fontsize = 10)
              )),
              row_title_gp = gpar(col="#00000000"),
              row_gap = unit(0.5, "mm"))

# draw(lgd12)
draw(h1, annotation_legend_list = lgd12)

cairo_pdf("figures/heatmap_rna_expression.1.pdf",height = 15/2.54,width=25/2.54,pointsize = 12*0.7)
draw(h1, annotation_legend_list = lgd12)
dev.off()

