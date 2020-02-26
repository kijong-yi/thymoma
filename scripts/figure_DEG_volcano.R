# /home/users/sypark/00_Project/01_thymoma/10_Final_data/17_Figures_for_publication
library(tidyverse)
library(ComplexHeatmap)
library(ggsci)
library(scales)
library(circlize)
library(ggrepel)
library(wordcloud)


# sypark's code --------------------------------------------------------------------------------------------------------
meta_dt <- read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/02_metadata/Thymoma_summary_191115_1stSheet.txt')
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
exp_dt_pcg <- left_join(exp_dt, bm_dt, by="gene") %>%
  dplyr::filter(gene_type == 'protein_coding') %>% dplyr::select(-gene_type) # filter protein coding only
l10_exp_dt <- exp_dt %>% mutate_at(vars(-gene), list(~log10(.+0.01)))
l10_exp_dt <- left_join(l10_exp_dt, bm_dt, by="gene") %>%
  dplyr::filter(gene_type == 'protein_coding') %>% dplyr::select(-gene_type) # filter protein coding only

# load cancer gene sensus
cgs_dt <- read_tsv("~sypark/02_Reference/14_cosmic/Cancer_gene_census_GRCh37_v89.tsv")
onc_dt <- cgs_dt %>% 
  dplyr::rename(gene = `Gene Symbol`, role_in_cancer = `Role in Cancer`) %>%
  dplyr::select(gene, role_in_cancer) %>% filter(grepl('oncogene',role_in_cancer)== T)
onc_dt$role_in_cancer <- 'oncogene'

# color setting
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
myred = pal_npg("nrc")(10)[8]
myblue = pal_npg("nrc")(10)[4]

# volcano plot - total sample
volc_dt <- l10_exp_dt
volc_dt$MT_mean <- rowMeans(volc_dt[,MT_ids])
volc_dt$WT_mean <- rowMeans(volc_dt[,WT_ids])
volc_dt$CA_mean <- rowMeans(volc_dt[,CA_ids])
volc_dt$WTCA_mean <- rowMeans(volc_dt[,c(CA_ids,WT_ids)])
volc_dt$WTMT_mean <- rowMeans(volc_dt[,c(MT_ids,WT_ids)])


my.t.test.p.value <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}
volc_dt$t.test.pvalue <- apply(volc_dt[,c(MT_ids, WT_ids)],1,function(x) my.t.test.p.value(x[MT_ids], x[WT_ids]))

volc_dt <- left_join(volc_dt, onc_dt)

ggplot(volc_dt, aes(x=WT_mean - MT_mean, y=-log10(t.test.pvalue)))+
  geom_point(aes(color = role_in_cancer), size=3, alpha=0.7)+
  geom_text_repel(data = subset(volc_dt, abs(WT_mean - MT_mean) >=2 | -log10(t.test.pvalue) >=40), aes(label = gene))+
  geom_hline(yintercept = -log10(0.05/nrow(volc_dt)), linetype="longdash")+
  scale_x_continuous(limits = c(-3,3))+
  theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line())+
  ggtitle('Total samples')

WT_high_genes <- volc_dt %>% filter(WT_mean - MT_mean >=1 & t.test.pvalue < 0.05/nrow(volc_dt)) %>% .$gene


# kjyi's code-------------------------------------------------------------------------------------------

volc_dt$t.test.pvalue
volc_dt$p.adjust = p.adjust(volc_dt$t.test.pvalue, method="BH")

WT_high_genes <- volc_dt %>% filter(WT_mean - MT_mean >=0 & p.adjust < 0.01) %>% .$gene
MT_high_genes <- volc_dt %>% filter(WT_mean - MT_mean <0 & p.adjust < 0.01) %>% .$gene
length(WT_high_genes)
length(MT_high_genes)
length(MT_high_genes) + length(WT_high_genes)
degs = list(WT_high = WT_high_genes,
            MT_high = MT_high_genes)
if(F) write_rds(degs,"data/degs.Rds")

(p.adjust(volc_dt$t.test.pvalue, method = "BH"))[volc_dt$gene=="IRS4"]

cairo_pdf("figures/volcano.1.pdf",width=1000/254,height = 800/254,pointsize = 12*0.7)
par(mar=c(3,3.5,0.5,0.5))

  x=volc_dt$WT_mean-volc_dt$MT_mean
  y=-log10(p.adjust(volc_dt$t.test.pvalue, method = "BH"))
  print(table(y>2))
  # y = -log10(t.test.pvalue)
  y[is.na(y)] = 0
  labelat = (abs(x)>2) | (y > 40)
  col=ifelse(is.na(volc_dt$role_in_cancer), "grey","red")
  greys=col=="grey"
  reds=col=="red"
  plot(x[greys],y[greys],bg="#00000020",col="#00000030",pch=21,
       xlim=c(-max(abs(x)),max(abs(x)))*1.25,
       ylim=c(0,max(y))*1.15,
       # ylab=expression(-log[10]~q),
       ylab="",
       xlab=""
       # xlab=expression((GTF2Imutant - WT))
       # xlab=expression(log[10]~FC)
  )
  title(ylab=expression(-log[10]~Adj.p), line=2)
  title(xlab=expression(log[10]~FC), line=2)
  # mtext(expression(log[10]~FC),side = 1,line=2)
  points(x[reds],y[reds],bg="#FF000060",col="#FF000070",pch=21)
  # wordcloud::textplot(x=x[labelat],y=y[labelat],words=gene[labelat],cex=1,new=F)
  set.seed(42)
  maptools::pointLabel(x=x[labelat],y=y[labelat],labels=volc_dt$gene[labelat],method = "GA")
  if(F){
    points(x=x[volc_dt$gene=="IGFBPL1"],y=y[volc_dt$gene=="IGFBPL1"],pch=21,col="#0000FFBB")
    text(x=x[volc_dt$gene=="IGFBPL1"],y=y[volc_dt$gene=="IGFBPL1"],labels = "  IGFBPL1",adj=0)
  }
  
  legend(-1.25*max(abs(x)),1.15*max(y),legend="in GTF2I-mutant",pch="↑",xjust=0,yjust=1,bty='n',text.font = 2)
  legend(1.25*max(abs(x)),1.15*max(y),legend="in wild-type",pch="↑",xjust=1,yjust=1,bty='n',text.font = 2)
  legend(1.25*max(abs(x)),1.15*max(y),
         legend = c("","known oncogene"),
         pt.bg=c("#00000000","#FF000060"),
         col=c("#00000000","#FF000070"),
         pch=21,
         xjust=1,yjust=1,bty='n'
  )
  abline(h=-log10(0.01),lty=2)
  text(par()$usr[1],-log10(0.01),"Adj.p = 0.01",adj=c(-0.1,-0.1))
  # points(par()$usr[1],-log10(0.05),col="red")

dev.off()

volc_dt$role_in_cancer %>% table

if(F){save.image("~/Projects/thymus_single_cell/final2/data/snapshop_for_volcanoplot.RData")}

load("~/Projects/thymus_single_cell/final2/data/snapshop_for_volcanoplot.RData")

MM <- volc_dt%>% as.data.frame %>% column_to_rownames("gene") %>% .[,meta_dt$id]%>% as.matrix
colnames(MM)
par(pty="s")

plot(MM["IGFBPL1",],MM["IRS4",],pch=21,cex=1.5,bg=gtf2i_pal[meta_dt$GTF2I_status2],
     xlab="IGFBPL1",ylab="IRS4 (log10(tpm+0.01))")

plot(10^MM["IGFBPL1",meta_dt$id]-0.01,10^MM["IRS4",meta_dt$id]-0.01,pch=21,cex=1.5,bg=gtf2i_pal[meta_dt$GTF2I_status2],
     xlab="IGFBPL1",ylab="IRS4 (log10(tpm+0.01))")





