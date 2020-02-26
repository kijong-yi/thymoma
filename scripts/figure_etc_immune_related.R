load("data/snapshot_etc_immune_related.RData")
#load libraries
library(tidyverse)
library(ComplexHeatmap)
library(ggsci)
library(scales)
library(circlize)

# load meta data
meta_dt <- read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/02_metadata/Thymoma_summary_191106_1stSheet.txt')
MT_ids <- meta_dt$id[meta_dt$GTF2I_status2 == 'm']
WT_ids <- meta_dt$id[meta_dt$GTF2I_status2 == 'w']

# load biomart data
bm_dt <- read_tsv('~sypark/02_Reference/13_biomart/biomart_human_geneID_transcriptID_hgncsymbol_genetype_ensemblgenename_190522.txt')
bm_dt <- bm_dt %>% rename(gene = `Gene name`, gene_type = `Gene type`) %>% select(gene, gene_type) %>% unique()

# load expression data
exp_dt <- read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/01_expression/IRS4_corrected_v2/thym_137s_tpm_geneSymbol_tumorOnly_IRS4cor.tsv')
exp_dt_pcg <- left_join(exp_dt, bm_dt) %>% filter(gene_type == 'protein_coding') %>% select(-gene_type) # filter protein coding only
l10_exp_dt <- exp_dt %>% mutate_at(vars(-gene), funs(log10(.+0.01)))
nrow(l10_exp_dt)
l10_exp_dt <- left_join(l10_exp_dt, bm_dt) %>% filter(gene_type == 'protein_coding') %>% select(-gene_type) # filter protein coding only
nrow(l10_exp_dt)

# load cancer gene sensus
cgs_dt <- read_tsv('~sypark/02_Reference/14_cosmic/Cancer_gene_census_GRCh37_v89.tsv')
onc_dt <- cgs_dt %>% rename(gene = `Gene Symbol`, role_in_cancer = `Role in Cancer`) %>% select(gene, role_in_cancer) %>% filter(grepl('oncogene',role_in_cancer)== T)
#onc_dt <- cgs_dt %>% rename(gene = `Gene Symbol`, role_in_cancer = `Role in Cancer`) %>% select(gene, role_in_cancer) %>% filter(role_in_cancer != 'fusion')
onc_dt$role_in_cancer <- 'oncogene'

# load immune score
im_dt <- read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/13_Danaher_immune_cell/thym_159s_Danaher_score_190923.tsv')
im_dt <- im_dt %>% select(-ends_with("_N")) 
colnames(im_dt) <- gsub("\\.","-",colnames(im_dt))

#load mixCR result
mixcr_dt <- read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/15_mixCR/Clone3_Fraction0.01_summary_tbl.tsv')

# load GTEx data
gtex_dt <- read_tsv('~sypark/02_Reference/27_GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.txt')
colnames(gtex_dt)



# color setting
my_pal = pal_npg("nrc")(10)
show_col(my_pal)
meta_dt$histologic_type %>% unique()
histo_pal = pal_aaas("default")(10)[c(4,1,7,8,6,2,3,5)]
show_col(histo_pal)
names(histo_pal) = c("A","AB","MN-T","B1","B2","B3","NE","TC")

stage_pal=c('#ffffd4','#fed98e','#fe9929','#d95f0e','#993404', 'white')
show_col(stage_pal)
names(stage_pal) <- c('I','II','III','IVa','IVb', '0')

cohort_pal = pal_jama("default", alpha = 0.8)(7)[c(1,7)]
names(cohort_pal) = c("SNUH","TCGA_CancerCell")
show_col(cohort_pal)

gtf2i_pal = pal_npg("nrc")(10)[c(4,1,3)]
show_col(gtf2i_pal)
names(gtf2i_pal) = c('m','w','c')

exp_pal = colorRamp2(c(-3,0,3), c('#253494',"gray90",'#f03b20'))

myred = pal_npg("nrc")(10)[8]
myblue = pal_npg("nrc")(10)[4]



# unsupervised clustering
hm_dt <- l10_exp_dt
hm_dt <- hm_dt %>% as.data.frame() %>% column_to_rownames('gene') 
hm_dt$var <- apply(hm_dt, 1, function(x) var(x))
hm_dt <- hm_dt %>% rownames_to_column('gene')
var_genes <- hm_dt %>% arrange(desc(var)) %>% .$gene %>% .[1:2500]
hm_dt <- hm_dt %>% filter(gene %in% var_genes) %>% select(-var) %>% column_to_rownames('gene') %>% as.matrix()
#hm_dt <- hm_dt %>% filter(var >= 0.5)  %>% select(-var) %>% as.matrix()
dim(hm_dt)
hm_annot_dt <- meta_dt %>% select(id, cohort, GTF2I_status2, histologic_type, Stage, ImmuneScore) %>% as.data.frame() %>% column_to_rownames('id')
hm_annot_dt <- hm_annot_dt[colnames(hm_dt),]
top_anno <- HeatmapAnnotation(gtf2i= hm_annot_dt$GTF2I_status2, hist= hm_annot_dt$histologic_type, cohort = hm_annot_dt$cohort,
                              col = list(gtf2i = gtf2i_pal, hist = histo_pal,cohort = cohort_pal))
hm_dt <- t(scale(t(hm_dt))) # row scaling
Heatmap(hm_dt, clustering_distance_columns = 'pearson', clustering_method_columns =  "average",top_annotation = top_anno, col= exp_pal, show_row_names = F)

# unsupervised clustering with immune related genes (deprecated)
import_dt <- read_tsv('~sypark/02_Reference/25_ImmPort/ImmuneGene_GeneSummary_from_immport.rename.txt')
import_dt
imr_genes <- unique(import_dt$new_name)
l10_exp_dt <- exp_dt %>% mutate_at(vars(-gene), funs(log10(.+0.01)))
hm_dt2 <- l10_exp_dt
hm_dt2 <- hm_dt2 %>% as.data.frame() %>% column_to_rownames('gene') 
length(imr_genes)
intersect(rownames(hm_dt2), imr_genes) %>% length()
hm_dt2 <- hm_dt2[intersect(rownames(hm_dt2), imr_genes),]
hm_dt2$var <- apply(hm_dt2, 1, function(x) var(x))
hm_dt2 <- hm_dt2 %>% rownames_to_column('gene')
var_genes <- hm_dt2 %>% arrange(desc(var)) %>% .$gene %>% .[1:50]
hm_dt2 <- hm_dt2 %>% filter(gene %in% var_genes) %>% select(-var) %>% column_to_rownames('gene') %>% as.matrix()
hm_annot_dt <- meta_dt %>% select(id, cohort, GTF2I_status2, histologic_type, Stage, ImmuneScore) %>% as.data.frame() %>% column_to_rownames('id')
hm_annot_dt <- hm_annot_dt[colnames(hm_dt2),]
top_anno <- HeatmapAnnotation(gtf2i= hm_annot_dt$GTF2I_status2, hist= hm_annot_dt$histologic_type, cohort = hm_annot_dt$cohort,
                              col = list(gtf2i = gtf2i_pal, hist = histo_pal,cohort = cohort_pal))
Heatmap(hm_dt2, clustering_distance_columns = 'euclidean', clustering_method_columns =  "complete",top_annotation = top_anno, col= exp_pal, row_names_gp=gpar(fontsize = 8), column_names_gp=gpar(fontsize=8))

# immune score corrplot
dim(im_dt)
fi_im_dt <- im_dt %>% filter(!cell_type %in% c('CD4_T_cells','CD8_T_cells','T_cells','CD45','Early_thymocytes'))
fi_im_dt <- im_dt %>% filter(!cell_type %in% c('Early_thymocytes'))
im_mx <- fi_im_dt %>% as.data.frame() %>% column_to_rownames('cell_type') %>% as.matrix()
im_cor <- cor(t(im_mx))
library(corrplot)
library(RColorBrewer)
mycol <- colorRampPalette(c(myblue, "white", myred)) 

cairo_pdf("figures/immune_sig_correlogram.pdf",height = 12/2.54,width=13/2.54,pointsize = 12*0.7)
x <-corrplot(im_cor, method="color", type="upper", order="hclust", col=mycol(50), tl.col="black", tl.srt=45)
x
dev.off()

# unsupervised clusteirng with immune scores
dim(im_dt)
fi_im_dt <- im_dt %>% filter(!cell_type %in% c('CD4_T_cells','CD8_T_cells','T_cells','CD45','Early_thymocytes'))
fi_im_dt <- im_dt # ????????????????????????????????????????????????????????????????????????????????????????????
im_mx <- fi_im_dt %>% as.data.frame() %>% column_to_rownames('cell_type') %>% as.matrix()
#meta_dt <- left_join(meta_dt, mixcr_dt)
meta_dt[is.na(meta_dt)] <- 0
colnames(meta_dt)
hm_annot_dt <- meta_dt %>% select(id, cohort, GTF2I_status2, histologic_type, Stage, ImmuneScore, history_myasthenia_gravis,sum_frac_adj) %>% as.data.frame() %>% column_to_rownames('id')
hm_annot_dt <- hm_annot_dt[colnames(im_mx),]
hm_annot_dt$Stage[which(hm_annot_dt$Stage==0)] = NA
top_anno <- HeatmapAnnotation("Group"= hm_annot_dt$GTF2I_status2, 
                              "Histology"= hm_annot_dt$histologic_type, 
                              "Cohort" = hm_annot_dt$cohort, 
                              "Masaoka stage" = hm_annot_dt$Stage,
                              # mg = hm_annot_dt$history_myasthenia_gravis,
                              col = list(Group = gtf2i_pal, Histology = histo_pal,Cohort = cohort_pal, "Masaoka stage" = stage_pal),
                              annotation_legend_param = list(
                                Group = list(title = "Group",
                                  at = c("m", "w","c"),
                                  labels = c("GTF2I-mutant", "Wild-type","Thymic carcinoma"))))
im_mx_scale <- t(scale(t(im_mx)))


cairo_pdf("figures/immune_score_heatmap.pdf",height = 14/2.54,width=28/2.54,pointsize = 12*0.7)

Heatmap(im_mx_scale[rownames(x),],
        clustering_distance_columns = 'pearson', clustering_method_columns =  "average",
        # clustering_distance_rows = 'pearson', clustering_method_rows =  "complete", 
        cluster_rows=F,name="Normalized enrichment score",show_column_names=F,
        top_annotation = top_anno, col= exp_pal, row_names_gp=gpar(fontsize = 8),
        width=unit(16,"cm"),height=unit(10,"cm"))
dev.off()




# immune score vs gene
mx1 <- im_dt %>% as.data.frame() %>% column_to_rownames('cell_type') %>% as.matrix()
mx2 <- l10_exp_dt %>% as.data.frame() %>% column_to_rownames('gene') %>% as.matrix()
dim(t(mx1))
dim(t(mx2))
rownames(t(mx1))
rownames(t(mx2))
cor_res <- cor(t(mx1), t(mx2)[rownames(t(mx1)),])
cor_res <- t(cor_res) %>% as.data.frame() %>% rownames_to_column('gene') %>% as.tibble()
cor_res %>% arrange(desc(Cytotoxic_cells))
dim(gtex_dt)
blood_genes <- gtex_dt %>% filter(`Whole Blood` >= 1) %>% .$Description
cor_res %>% filter(!gene %in% blood_genes) %>% arrange(desc(Cytotoxic_cells)) %>% filter(Cytotoxic_cells >= 0.6 | Exhausted_CD8 >= 0.6) %>% select(gene, Cytotoxic_cells, Exhausted_CD8)
cor_res %>% filter(!gene %in% blood_genes) %>% arrange(desc(Cytotoxic_cells)) %>% filter(Cytotoxic_cells >= 0.4) %>% .$gene

gene_dt <- l10_exp_dt %>% filter(gene %in% c('CXCL9','CD70','IDO1','P2RY6','DTHD1')) %>% as.data.frame() %>% column_to_rownames('gene') %>% t() %>% as.data.frame() %>% rownames_to_column('id') %>% as.tibble()
cyto_dt <- im_dt %>% filter(cell_type == 'Cytotoxic_cells') %>% as.data.frame() %>% column_to_rownames('cell_type') %>% t() %>% as.data.frame() %>% rownames_to_column('id') %>% as.tibble()
gene_cyto_dt <- left_join(cyto_dt, gene_dt)
gene_cyto_dt <- left_join(gene_cyto_dt, meta_dt)
dim(gene_cyto_dt)
colnames(gene_cyto_dt)
ggplot(gene_cyto_dt, aes(x=Cytotoxic_cells, y=DTHD1))+
  geom_point(aes(color=histologic_type))+
  scale_color_manual(values = histo_pal)+
  theme_bw()+ theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line())+
  ylab('log10(TPM+0.01)')+xlab("Score of cytotoxic cells")+ggtitle("DTHD1")


target_genes1 <- c('CD274','PDCD1LG2','CD80','CD86','CD276','VTCN1','TNFRSF14','LGALS9','TNFSF18','PVR', 'IDO1','NECTIN2')
label_names1 <- c('PD-L1','PD-L2','CD80','CD86','B7-H3','B7-H4','HVEM','GAL9','GITRL','CD155','IDO1','CD112')
names(label_names1) = target_genes1
# length(target_genes)
cor_res %>% filter(gene %in% target_genes1) %>% nrow()
library(ggrepel)
p1 <- ggplot(subset(cor_res, !cor_res$gene %in% blood_genes), aes(x=Cytotoxic_cells, y=Exhausted_CD8))+
  geom_point(alpha=0.3)+
  geom_point(data=subset(cor_res, cor_res$gene %in% target_genes1), aes(x=Cytotoxic_cells, y=Exhausted_CD8),color = "red")+ 
  geom_text_repel(data=subset(cor_res, cor_res$gene %in% target_genes1), aes(x=Cytotoxic_cells, y=Exhausted_CD8, label=label_names1[gene]), color = "red")+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_vline(xintercept = 0, linetype="dashed")+
  scale_x_continuous(breaks = seq(-0.4, 0.8, 0.2), limits = c(-0.6, 1))+
  scale_y_continuous(breaks = seq(-0.4, 0.8, 0.2), limits = c(-0.6, 1))+
  theme_bw()+theme(panel.grid=element_blank(), panel.border = element_blank(), axis.line = element_line())+
  xlab("Correlation with Cytotoxic cell score")+ylab("Correlation with Exhausted CD8 score")+ggtitle("Non-immune cells")
p1

target_genes <- c('PDCD1','CTLA4','BTLA','LAG3','CD96','TIGIT','PVRIG','HAVCR2',"TNFRSF18")
label_names <- c('PD1','CTLA-4','BTLA','LAG-3','CD96','TIGIT','CD112R','TIM-3',"GITR!!!!!!!!!!!")
names(label_names) = target_genes
length(target_genes)
cor_res %>% filter(gene %in% target_genes) %>% nrow()
library(ggrepel)
p2 <- ggplot(subset(cor_res, cor_res$gene %in% blood_genes), aes(x=Cytotoxic_cells, y=Exhausted_CD8))+
  geom_point(alpha=0.3)+
  geom_point(data=subset(cor_res, cor_res$gene %in% target_genes), aes(x=Cytotoxic_cells, y=Exhausted_CD8),color = "red")+ 
  geom_text_repel(data=subset(cor_res, cor_res$gene %in% target_genes), aes(x=Cytotoxic_cells, y=Exhausted_CD8, label=label_names[gene]), color = "red")+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_vline(xintercept = 0, linetype="dashed")+
  scale_x_continuous(breaks = seq(-0.4, 0.8, 0.2), limits = c(-0.6, 1))+
  scale_y_continuous(breaks = seq(-0.4, 0.8, 0.2), limits = c(-0.6, 1))+
  theme_bw()+theme(panel.grid=element_blank(), panel.border = element_blank(), axis.line = element_line())+
  xlab("Correlation with Cytotoxic cell score")+ylab("Correlation with Exhausted CD8 score")+ggtitle("Immune cells")
print(p1)
print(p2)

tmp_dat <- subset(cor_res, !cor_res$gene %in% blood_genes) %>% dplyr::select(Cytotoxic_cells, Exhausted_CD8) %>% unique
seed=41
p1 <- ggplot(tmp_dat, 
       aes(x=Cytotoxic_cells, y=Exhausted_CD8))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_vline(xintercept = 0, linetype="dashed")+
  stat_density2d(aes(alpha=..level..), geom="polygon",show.legend = FALSE,bins = 30) +
  ggtitle(paste0("Other genes (n=",nrow(tmp_dat),")"))+ 
  # geom_point(colour="black",alpha=0.02)+
  geom_point(data=subset(cor_res, cor_res$gene %in% target_genes1), aes(x=Cytotoxic_cells, y=Exhausted_CD8),color = "red")+ 
  scale_alpha_continuous(limits=c(0,16),breaks=seq(0,16,by=2))+
  geom_label_repel(data=subset(cor_res, cor_res$gene %in% target_genes1), 
                   aes(x=Cytotoxic_cells, y=Exhausted_CD8, label=label_names1[gene]), 
                   color = "black",fill="white",force=3,
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.5, "lines"))+
  xlab("Cytotoxic cell score") + ylab("Exhausted CD8 score") +
  theme_bw()

tmp_dat2 <- subset(cor_res, cor_res$gene %in% blood_genes) %>% dplyr::select(Cytotoxic_cells, Exhausted_CD8) %>% unique

p2 <-  ggplot(tmp_dat2,
       aes(Cytotoxic_cells,Exhausted_CD8))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_vline(xintercept = 0, linetype="dashed")+
  stat_density2d(aes(alpha=..level..), geom="polygon",show.legend = FALSE,bins = 30) +
  ggtitle(paste0("GTEx blood exp. genes (n=",nrow(tmp_dat2),")"))+ 
  # scale_alpha_continuous(limits=c(0,5),breaks=seq(0,5,by=1))+
  # geom_point(colour="black",alpha=0.02)+
  geom_point(data=subset(cor_res, cor_res$gene %in% target_genes), aes(x=Cytotoxic_cells, y=Exhausted_CD8),color = "red")+ 
  geom_label_repel(data=subset(cor_res, cor_res$gene %in% target_genes),
                   aes(x=Cytotoxic_cells, y=Exhausted_CD8, label=label_names[gene]),
                   color = "black",fill="white",force=3,
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.5, "lines"),seed=seed)+ 
  xlab("Correlation with cytotoxic cell score") + ylab("Correlation with exhausted CD8 score") +
  theme(legend.position = "none")+
  theme_bw()
cairo_pdf("figures/scatter_corr_exhaustion_related.pdf",height =10/2.54,width=20/2.54,pointsize = 12*0.7)
cowplot::plot_grid(p2, p1)


dev.off()

#mixCR - Cytotoxic_cells
merged_dt <- left_join(meta_dt, mixcr_dt)
t_im_dt <- im_dt %>% as.data.frame() %>% column_to_rownames('cell_type') %>% t() %>% as.data.frame() %>% rownames_to_column('id') %>% as.tibble()
merged_dt <- left_join(merged_dt, t_im_dt)
colnames(merged_dt)
merged_dt[is.na(merged_dt)] <- 0
library(ggrepel)

cairo_pdf("figures/scatter_clone_cytotoxic_corr.pdf",height = 9/2.54,width=13/2.54,pointsize = 12*0.7)
ggplot(merged_dt, aes(x=Cytotoxic_cells, y=sum_frac_adj))+
  geom_point(aes(color=histologic_type, size=num_clones), alpha=0.7)+
  #geom_text_repel(data=subset(merged_dt, sum_frac_adj >= 0.1), aes(label=id))+
  scale_color_manual(values = histo_pal)+
  scale_size(range=c(3,9), breaks = c(0,2,4,6))+
  theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line())+
  xlab("Cytotoxic cell score")+ylab("Fraction of clonally expanded T cell") + labs(color="Histologic type", size="Number of clones")
dev.off()

colnames(meta_dt)


if(F){save.image("data/snapshot_etc_immune_related.RData",compress="gzip")}

