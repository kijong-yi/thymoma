library(tidyverse)
# library(plotly)
library(plot3D)
library(ggsci)


meta_dt <- read_tsv(paste0('~sypark/00_Project/01_thymoma/10_Final_data/02_metadata',
                           '/Thymoma_summary_190822_1stSheet.txt'))

bm_dt <- read_tsv(paste0('~sypark/02_Reference/13_biomart/',
                         'biomart_human_geneID_transcriptID_hgncsymbol',
                         '_genetype_ensemblgenename_190522.txt')) %>%
  dplyr::rename(gene = `Gene name`, gene_type = `Gene type`) %>% dplyr::select(gene, gene_type) %>% unique()

exp_dt <- read_tsv(paste0("~sypark/00_Project/01_thymoma/10_Final_data/01_expression/",
                          'IRS4_corrected_v2/',
                          'thym_137s_tpm_geneSymbol_tumorOnly_IRS4cor.tsv'))

l10_exp_dt <- exp_dt %>% mutate_at(vars(-gene), list(~log10(.+0.01)))
l10_exp_dt <- left_join(l10_exp_dt, bm_dt, by="gene") %>%
  dplyr::filter(gene_type == 'protein_coding') %>% dplyr::select(-gene_type) # filter protein coding only

hm_dt <- l10_exp_dt
hm_dt <- hm_dt %>% as.data.frame() %>% column_to_rownames('gene') 
hm_dt$var <- apply(hm_dt, 1, function(x) var(x))
hm_dt <- hm_dt %>% rownames_to_column('gene')
var_genes <- hm_dt %>% arrange(desc(var)) %>% .$gene %>% .[1:2500] # using top 2500 high variance genes, variance in log10 scale
hm_dt <- hm_dt %>% dplyr::filter(gene %in% var_genes) %>% dplyr::select(-var) %>% column_to_rownames('gene') %>% as.matrix()
#hm_dt <- hm_dt %>% filter(var >= 0.5)  %>% select(-var) %>% as.matrix()
# dim(hm_dt)


pca_exp <- prcomp(t(hm_dt),scale=T)

# plot_ly(as.data.frame(pca_exp$x[,1:3]), x = ~PC1, y = ~PC2, z = ~PC3,
#         color = ~gtf2i_pal[meta_dt$GTF2I_status2])

gtf2i_pal = pal_npg("nrc")(10)[c(4,1,3)]
names(gtf2i_pal) = c('m','w','c')


# cairo_pdf("figures/pca_gexp.1.pdf",height = 5.5/2.54,width=5.5/2.54,pointsize = 12*0.6)
cairo_pdf("figures/pca_gexp.1.pdf",height = 10/2.54,width=10/2.54,pointsize = 12*0.7)
par(mar=c(1,1,1,1))
scatter3D(-pca_exp$x[,1], -pca_exp$x[,2], pca_exp$x[,3],
          colvar = as.integer(factor(meta_dt$GTF2I_status2,levels=names(gtf2i_pal))),
          col="grey70",
          lwd=0.5,lty=1,
          type = "h", xlim=c(-35,50),ylim=c(-45,25),zlim=c(-50,35),
          phi = 40, theta = 50,
          xlab="PC1",ylab="PC2",zlab="PC3",
          bty='u',col.axis="grey",
          ticktype = "detailed",cex=1.5,cex.axis=1)#,
          # sub = "Principle component space of 2500 high variance gene expression data")
points3D(-pca_exp$x[,1], -pca_exp$x[,2], pca_exp$x[,3],
        # colvar = as.integer(factor(meta_dt$GTF2I_status2,levels=names(gtf2i_pal))),
        bg=gtf2i_pal[meta_dt$GTF2I_status2],add=T,colkey=F,cex=2,pch=21,col="black",bgkey=F)

legend("topright",c("GTF2I-mutant", "GTF2I-wildtype", "Carcinoma"),col="black",pch=21,cex=1,pt.bg=gtf2i_pal,pt.cex=1.5)
dev.off()


