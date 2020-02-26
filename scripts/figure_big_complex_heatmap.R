# read_tsv('abcd')

#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(tidyverse)
library(ggsignif)
library(circlize)
# Load meta data
meta_dt <- read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/02_metadata/Thymoma_summary_200217_1stSheet.txt')
meta_dt <- meta_dt[is.na(meta_dt$id) == F, ]
meta_dt$age_at_diagnosis %>% summary()
meta_dt$TCGA_paper_GTF2Imt %>% unique()
meta_dt <- meta_dt %>% mutate(gtf2i_rescued = ifelse(TCGA_paper_GTF2Imt == 'w' & GTF2I_status2 == 'm', 'Y','N'))

# Load IM data and merged into meta_dt
im_dt <- read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/13_Danaher_immune_cell/thym_137s_Danaher_score_191017.tsv')
im_dt <- im_dt %>% as.data.frame() %>% column_to_rownames('cell_type') %>% t() %>% as.data.frame() %>% rownames_to_column('id') %>% as.tibble()
nrow(meta_dt)
meta_dt <- left_join(meta_dt, im_dt)
nrow(meta_dt)

# Load oncogrid data and manipulation
dt <- read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/16_point_mutation_final_call/05_recurrent_gene/recurrently_mutated_genes.sort.2.txt.selected')
dt <- dt %>% dplyr::rename(gene_name = `#gene_name`) %>% dplyr::select(gene_name, sample_id, mt_type)
dt$mt_type %>% unique()
#assign number to type of mutations
#101:missense mutation, 104:in-frame indel, 102:nonsense mutation, 105:frame-shift indel, 103:splice site mutation
dt[dt == 'nonsynonymous SNV'] <- 101
dt[dt == 'nonframeshift insertion'] <- 104
dt[dt == 'stopgain'] <- 102
dt[dt == 'frameshift deletion'] <- 105
dt[dt == 'splicing'] <- 103
dt[dt == 'frameshift insertion'] <- 105
dt$mt_type %>% unique()
dt <- unique(dt)
rescued_samples <- meta_dt %>% filter(gtf2i_rescued == 'Y') %>% .$id
dt$mt_type[dt$gene_name == 'GTF2I' & dt$sample_id %in% rescued_samples] <- 101.1

mdt <- as.matrix(spread(dt, sample_id, mt_type) %>% as.data.frame() %>% column_to_rownames('gene_name'))
na_dt <- matrix(ncol=length(setdiff(meta_dt$id,colnames(mdt))), nrow=nrow(mdt))
colnames(na_dt) <- setdiff(meta_dt$id,colnames(mdt))
oncodt <- cbind(mdt, na_dt)
epi_genes <- c('BAZ1A','BCOR','CHD2','SMARCA4')
pwy_genes <- c('HRAS','KIT','KRAS','NF2','NRAS','PIK3CA','PIK3R1')
cycle_genes <- c('CDKN2A','RPL22','TP53')
other_genes <- c('CYLD','IDH1','MSH2','SF3B1')
gene_dt <- rowSums(is.na(oncodt)== F) %>% as.data.frame() %>% rownames_to_column('gene') %>% rename(count='.')  %>% mutate(group= ifelse(gene == 'GTF2I', 'gtf2i',ifelse(gene %in% epi_genes, 'epi', ifelse(gene %in% pwy_genes, 'pwy', ifelse(gene %in% cycle_genes,'cycle', 'other'))))) %>% arrange(group, desc(count))
gene_dt$group <- factor(gene_dt$group, levels = c('gtf2i','pwy','epi','cycle','other'))
onco_gene_order <- gene_dt %>% arrange(group, desc(count)) %>% .$gene
#onco_gene_order <- rowSums(is.na(oncodt)== F) %>% as.data.frame() %>% rownames_to_column('gene') %>% rename(count='.') %>% arrange(desc(count)) %>% .$gene
oncodt <- oncodt[onco_gene_order,]
gene_num <- nrow(oncodt)

# Load CN data
# setwd("")
cndt <- read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/18_sequenza/01_segments/TETs_137s_meanArmCN_191113.txt')
arm_num <- nrow(cndt)

# Load expression data
# setwd("")
exp_dt <- read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/01_expression/IRS4_corrected_v2/thym_137s_tpm_geneSymbol_tumorOnly_IRS4cor.tsv')

# Loaad scRNA score
thymoma_scores <- readRDS("/home/users/kjyi/Projects/thymus_single_cell/final/models/thymoma_scores.Rds")
thymoma_scores <- thymoma_scores %>% as.data.frame() %>% rownames_to_column('id') %>% as.tibble()
meta_dt <- left_join(meta_dt, thymoma_scores)

# color setting
library(ggsci)
library("scales")
my_pal = pal_jama("default")(10)
my_pal = pal_jco("default")(10)
my_pal = pal_aaas("default")(10)
show_col(my_pal)


meta_dt$histologic_type %>% unique()
histo_pal = pal_aaas("default")(10)[c(4,1,7,8,6,2,3,5)]
show_col(histo_pal)
histo_pal= c(histo_pal,"black")
names(histo_pal) = c("A","AB","MN-T","B1","B2","B3","NE","CA-SqCC","CA-UN")

stage_pal=c('#ffffd4','#fed98e','#fe9929','#d95f0e','#993404')
show_col(stage_pal)
names(stage_pal) <- c('I','II','III','IVa','IVb')
cohort_pal = pal_jama("default", alpha = 0.8)(7)[c(1,7)]
names(cohort_pal) = c("SNUH","TCGA_CancerCell")
show_col(cohort_pal)
#101:missense mutation, 104:in-frame indel, 102:nonsense mutation, 105:frame-shift indel, 103:splice site mutation
onco_pal = pal_lancet("lanonc")(5)
names(onco_pal) = c("missense mutation","nonsense mutation","splice site mutation", "in-frame indel","frame-shift indel")
show_col(onco_pal)
cn_pal = colorRamp2(c(0,2,5), c('#253494',"gray90",'#f03b20'))

gtf2i_pal = pal_npg("nrc")(10)[c(4,1,3)]
show_col(gtf2i_pal)
names(gtf2i_pal) = c('m','w','c')

sig_pal=c("#DC143C", "#FF69B4", "#9C9C9C", "#36648B", "#9C661F", "#000000", "#D2691E", "#FFAEB9", "#00E5EE")
names(sig_pal)=c("Signature 1","Signature 2","Signature 4","Signature 5","Signature 6","Signature 13","Signature 15","Signature 17","Signature 18")


library(circlize)
age_pal = colorRamp2(c(20,70), c("white","darkgrey"))
purity_pal = colorRamp2(c(0,1), c("white", pal_jama("default", alph = 1)(7)[1]))
occn_pal = colorRamp2(c(0, 2, 5, 101, 101.1, 102, 103, 104, 105), c('#253494',"gray90",'#f03b20', pal_lancet("lanonc")(5)[1], pal_lancet("lanonc")(5)))
cytotoxic_pal = colorRamp2(c(-0,1), c("white",pal_jama("default", alph = 1)(7)[4]))
thymocyte_pal = colorRamp2(c(-1,2.5), c("white",pal_jama("default", alph = 1)(7)[3]))
cTEC_pal = colorRamp2(c(-10,7),c("white",pal_jama("default", alph = 1)(7)[6]))
mTEC_pal = colorRamp2(c(-3,7), c("white",pal_aaas("default")(10)[3]))
tuft_pal = colorRamp2(c(-2,3),c("white",pal_aaas("default")(10)[5]))
progenitor_pal = colorRamp2(c(-1,1.5),c("white",pal_aaas("default")(10)[8]))

meta_dt$histologic_type <- meta_dt$histologic_type2
meta_dt$histologic_type2 %>% table
# MT plot
meta_dt$histologic_type <- factor(meta_dt$histologic_type, levels = c("A","AB","MN-T","B1","B2","B3","CA-SqCC","CA-UN","NE"))
MT_meta_dt <- meta_dt %>% filter(GTF2I_status2 == 'm') %>% arrange(histologic_type,desc(n_pointmt/bait_size))
histo_num <- MT_meta_dt %>% group_by(histologic_type) %>% count() %>% .$n
names(histo_num) <- MT_meta_dt %>% group_by(histologic_type) %>% count() %>% .$histologic_type
MT_col_split <- c(rep('A', histo_num['A']), rep('AB', histo_num['AB']), rep('MN-T',histo_num['MN-T']), rep('B1', histo_num['B1']), rep('B2', histo_num['B2']), rep('B3', histo_num['B3']))
MT_col_split = factor(MT_col_split, levels = c('A','AB','MN-T','B1','B2','B3'))

MT_top= HeatmapAnnotation(n_pointmt = anno_barplot(MT_meta_dt$n_pointmt/MT_meta_dt$bait_size,bar_width=1, ylim=c(0,3),gp=gpar(fill = "darkgrey", col = "darkgrey")),
                          hist = MT_meta_dt$histologic_type,
                          cohort = MT_meta_dt$cohort,
                          age = MT_meta_dt$age_at_diagnosis,
                          purity= MT_meta_dt$final_cellularity,
                          stage = MT_meta_dt$Stage,
                          irs4 = anno_barplot(log10(MT_meta_dt$corrected_IRS4_TPM+0.1)+1, ylim= c(0,4), axis_param = list(at=log10(c(0, 1, 10, 100)+0.1)+1, labels=c(0,1,10,100))),
                          show_annotation_name = T, 
                          col = list(hist = histo_pal, stage = stage_pal, cohort = cohort_pal, age = age_pal, purity = purity_pal, thymocyte = thymocyte_pal, cytotoxic = cytotoxic_pal,
                                     cTEC=cTEC_pal, mTEC=mTEC_pal, tuft=tuft_pal, progenitor=progenitor_pal),
                          show_legend = c(FALSE)
)

MT_onco_dt <- oncodt[,MT_meta_dt$id] 
MT_cn_dt <- cndt %>% dplyr::select(chr_arm, MT_meta_dt$id) %>% as.data.frame() %>% column_to_rownames('chr_arm') %>% as.matrix()
colnames(MT_onco_dt) == colnames(MT_cn_dt)
MT_occn_dt <- rbind(MT_onco_dt, MT_cn_dt)
class(MT_occn_dt) <- "numeric"
chrom_list <- as.numeric(gsub('X','23',gsub('q','',gsub('p','',cndt$chr_arm))))
univ_row_split = c(rep(0.1,1), rep(0.2,7), rep(0.3,4), rep(0.4,3), rep(0.5,4), chrom_list)
onco_pct <- c(paste0(round(rowSums(is.na(MT_onco_dt)==F)*100/ncol(MT_onco_dt),1),'%'), rep('',nrow(MT_cn_dt)))
MT_right= rowAnnotation(pct = anno_text(onco_pct, gp = gpar(fontsize = 9)))
MT_body <- Heatmap(MT_occn_dt, top_annotation = MT_top, right_annotation = MT_right, cluster_rows = F, cluster_columns = F, col=occn_pal, show_heatmap_legend = F, na_col = "gray90", 
                   row_split = univ_row_split, row_gap = unit(1,'mm'), row_title = NULL,
                   column_split = MT_col_split, column_gap = unit(1,'mm'), column_title = NULL, 
                   row_names_side = "left", row_names_gp=gpar(fontsize=9), show_column_names = F, 
                   cell_fun = function(j, i, x, y, width, height, fill){
                     if(is.na(MT_occn_dt[i,j]) == F){if (MT_occn_dt[i,j] == 101.1){
                       grid.points(x=x, y=y, pch=16, size=unit(2,'mm'), gp = gpar(col = "white"))}}
                   })
if(F){
  draw(MT_body)
}


# WT plot
WT_meta_dt <- meta_dt %>% filter(GTF2I_status2 == 'w') %>% arrange(histologic_type,desc(n_pointmt/bait_size))
histo_num <- WT_meta_dt %>% group_by(histologic_type) %>% count() %>% .$n
names(histo_num) <- WT_meta_dt %>% group_by(histologic_type) %>% count() %>% .$histologic_type
WT_col_split <- c(rep('AB', histo_num['AB']), rep('B1', histo_num['B1']), rep('B2', histo_num['B2']), rep('B3', histo_num['B3']))
WT_top= HeatmapAnnotation(n_pointmt = anno_barplot(WT_meta_dt$n_pointmt/WT_meta_dt$bait_size,bar_width=1, 
                                                   ylim=c(0,3),gp=gpar(fill = "darkgrey", col = "darkgrey")),
                          hist = WT_meta_dt$histologic_type,
                          cohort = WT_meta_dt$cohort,
                          age = WT_meta_dt$age_at_diagnosis,
                          purity= WT_meta_dt$final_cellularity,
                          stage = WT_meta_dt$Stage,
                          irs4 = anno_barplot(log10(WT_meta_dt$corrected_IRS4_TPM+0.1)+1, ylim= c(0,4),
                                              axis_param = list(at=log10(c(0, 1, 10, 100)+0.1)+1, labels=c(0,1,10,100))),
                          show_annotation_name = F, 
                          col = list(hist = histo_pal, cohort = cohort_pal, age = age_pal, purity = purity_pal, stage = stage_pal,
                                     cytotoxic = cytotoxic_pal, thymocyte = thymocyte_pal,
                                     cTEC=cTEC_pal, mTEC=mTEC_pal, tuft=tuft_pal, progenitor=progenitor_pal),
                          show_legend = c(FALSE)
)
WT_onco_dt <- oncodt[,WT_meta_dt$id] 
WT_cn_dt <- cndt %>% dplyr::select(chr_arm, WT_meta_dt$id) %>% as.data.frame() %>% column_to_rownames('chr_arm') %>% as.matrix()
colnames(WT_onco_dt) == colnames(WT_cn_dt)
WT_occn_dt <- rbind(WT_onco_dt, WT_cn_dt)
class(WT_occn_dt) <- "numeric"
onco_pct <- c(paste0(round(rowSums(is.na(WT_onco_dt)==F)*100/ncol(WT_onco_dt),1),'%'), rep('',nrow(WT_cn_dt)))
WT_right= rowAnnotation(pct = anno_text(onco_pct, gp = gpar(fontsize =9)))
WT_body <- Heatmap(WT_occn_dt, top_annotation = WT_top, right_annotation = WT_right, cluster_rows = F, cluster_columns = F, col=occn_pal, show_heatmap_legend = F, na_col = "gray90", 
                   row_split = univ_row_split, row_gap = unit(1,'mm'), row_title = NULL, 
                   column_split = WT_col_split, column_gap = unit(1,'mm'), column_title = NULL, 
                   row_names_side = "left", show_column_names = F)

# draw(WT_body)

# TC plot
TC_meta_dt <- meta_dt %>% filter(GTF2I_status2 == 'c') %>% arrange(histologic_type=="NE",desc(n_pointmt/bait_size))
TC_meta_dt$n_pointmt/TC_meta_dt$bait_size
histo_num <- TC_meta_dt %>% group_by(histologic_type) %>% count() %>% as.data.frame %>% column_to_rownames("histologic_type") %>% t
TC_col_split <- c(rep('TC', histo_num[,'TC']), rep("NE",histo_num[,'NE']))
TC_col_split = factor(TC_col_split, levels = c('TC','NE'))

TC_top= HeatmapAnnotation(n_pointmt = anno_barplot(TC_meta_dt$n_pointmt/TC_meta_dt$bait_size,bar_width=1, ylim=c(0,3),gp=gpar(fill = "darkgrey", col = "darkgrey")),
                          hist = TC_meta_dt$histologic_type,
                          cohort = TC_meta_dt$cohort,
                          age = TC_meta_dt$age_at_diagnosis,
                          purity= TC_meta_dt$final_cellularity,
                          stage = TC_meta_dt$Stage,
                          # cytotoxic = TC_meta_dt$Cytotoxic_cells,
                          # thymocyte = TC_meta_dt$Late_thymocytes,
                          # cTEC = TC_meta_dt$cTEC,
                          # mTEC = TC_meta_dt$mTEC,
                          # tuft = TC_meta_dt$Tuft,
                          # progenitor = TC_meta_dt$Progenitor,
                          irs4 = anno_barplot(log10(TC_meta_dt$corrected_IRS4_TPM+0.1)+1, ylim= c(0,4), axis_param = list(at=log10(c(0, 1, 10, 100)+0.1)+1, labels=c(0,1,10,100))),
                          show_annotation_name = F, 
                          col = list(hist = histo_pal, cohort = cohort_pal, age = age_pal, purity = purity_pal, stage = stage_pal, cytotoxic = cytotoxic_pal, thymocyte = thymocyte_pal,
                                     cTEC=cTEC_pal, mTEC=mTEC_pal, tuft=tuft_pal, progenitor=progenitor_pal),
                          show_legend = c(FALSE)
)
TC_onco_dt <- oncodt[,TC_meta_dt$id] 
TC_cn_dt <- cndt %>% dplyr::select(chr_arm, TC_meta_dt$id) %>% as.data.frame() %>% column_to_rownames('chr_arm') %>% as.matrix()
colnames(TC_onco_dt) == colnames(TC_cn_dt)
TC_occn_dt <- rbind(TC_onco_dt, TC_cn_dt)
class(TC_occn_dt) <- "numeric"
onco_pct <- c(paste0(round(rowSums(is.na(TC_onco_dt)==F)*100/ncol(TC_onco_dt),1),'%'), rep('',nrow(TC_cn_dt)))
TC_right= rowAnnotation(pct = anno_text(onco_pct, gp = gpar(fontsize =9)))
TC_body <- Heatmap(TC_occn_dt, top_annotation = TC_top, right_annotation = TC_right, cluster_rows = F, cluster_columns = F, col=occn_pal, show_heatmap_legend = F, na_col = "gray90", 
                   row_split = univ_row_split, row_gap = unit(1,'mm'), row_title = NULL, row_names_side = "left", show_column_names = F,
                   column_split = TC_col_split, column_gap = unit(1,'mm'), column_title = NULL)

if(F){
  MT_body + WT_body + TC_body
  draw(TC_body)
 }
# draw(TC_body)

# Draw merged plot
plot_list <- MT_body + WT_body + TC_body
# setwd("~/00_Project/01_thymoma/10_Final_data/17_Figures_for_publication")
# pdf('./mutation_CN_heatmap/MT_CN_heatmap_191118.pdf',width=17, height=10)

# dev.off()

# Draw packed legend
lgd_hist = Legend (title = "Histologic types", legend_gp=gpar(fill = histo_pal[c(1,2,3,4,5,6,8,7)]), labels = names(histo_pal)[c(1,2,3,4,5,6,8,7)])
lgd_cohort = Legend (title = "Cohort", legend_gp=gpar(fill = cohort_pal), labels = names(cohort_pal))
lgd_age = Legend (title = "Age", col_fun = age_pal)
lgd_purity = Legend (title = "Tumor cell fraction", col_fun = purity_pal)
lgd_stage = Legend (title = "Stage", legend_gp=gpar(fill = stage_pal), labels = names(stage_pal))
# lgd_cytotoxic = Legend (title = "Cytotoxic cell score", col_fun = cytotoxic_pal)
# lgd_thymocyte = Legend (title = "Thymocyte score", col_fun = thymocyte_pal)
# lgd_cTEC = Legend (title = "cTEC score", col_fun = cTEC_pal)
# lgd_mTEC = Legend (title = "mTEC score", col_fun = mTEC_pal)
# lgd_tuft = Legend (title = "Tuft score", col_fun = tuft_pal)
# lgd_progenitor = Legend (title = "Progenitor score", col_fun = progenitor_pal)
lgd_onco = Legend (title = "Types of mutations", legend_gp=gpar(fill = onco_pal), labels = names(onco_pal))
# lgd_onco2 = Legend (title = "Types of mutations", legend_gp=gpar(fill = onco_pal), labels = names(onco_pal))
lgd_mCN = Legend (title = "Mean copy number", col_fun = cn_pal)

# plgd1 = packLegend(lgd_hist, lgd_cohort, lgd_age, lgd_purity, lgd_stage, max_width = unit(5,'cm'))
# plgd2 = packLegend(lgd_cytotoxic, lgd_thymocyte, lgd_cTEC, lgd_mTEC, lgd_tuft, lgd_progenitor, lgd_onco, lgd_mCN, max_width = unit(5,'cm'))

plgd1 = packLegend(lgd_hist, lgd_cohort, lgd_age, lgd_purity, lgd_stage, max_width = unit(5,'cm'))
plgd2 = packLegend(lgd_onco, lgd_mCN, max_width = unit(5,'cm'))


plgd3 = packLegend(lgd_hist, lgd_cohort, lgd_age, lgd_purity, lgd_stage, max_width = unit(5,'cm'),lgd_onco, lgd_mCN)

# cairo_pdf("figures/bigbigcomplexheatmap.pdf",width=3500/254,height = 2700/254,pointsize = 12*0.7)
if(F){
  cairo_pdf("figures/bigbigcomplexheatmap.pdf",width=1.7*1900/254,height = 1.7*1600/254,pointsize = 12*0.7/1.8)
  draw(plot_list, ht_gap = unit(0.5, "cm"),annotation_legend_list = plgd3)
  dev.off()
}


