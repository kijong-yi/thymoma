read_tsv('abc')

library(tidyverse)

meta_dt <- read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/02_metadata/Thymoma_summary_191205_1stSheet.txt')
tcga_files <- list.files('/home/users/sypark/00_Project/01_thymoma/04_TCGA/08_RNAseq_hg19/mixCR/', pattern = 'ALL.txt$', full.names = T)
snuh_files <- list.files('/home/users/sypark/00_Project/01_thymoma/01_RNAseq/04_mixCR/', pattern = 'ALL.txt$', full.names = T)
file_list <- c(tcga_files, snuh_files); rm(tcga_files, snuh_files)
gtf2i_pal = c(m="#3C5488FF",w="#E64B35FF",c="#00A087FF")

result_tbl2 <- lapply(1:length(file_list),function(n){
  dt <- read_tsv(file_list[n])
  data.frame(
    sampleid = ifelse(grepl('TCGA',file_list[n]), unlist(strsplit(unlist(strsplit(file_list[n],'//'))[2],'-01A-'))[1], unlist(strsplit(unlist(strsplit(file_list[n],'//'))[2],'\\.'))[1]),
    Bcell_tot_clone_sum = dt %>% filter(substr(allVHitsWithScore,1,2)=='IG') %>% .$cloneCount %>% sum(),
    Bcell_expan_clone_sum = dt %>% filter(substr(allVHitsWithScore,1,2)=='IG' & cloneCount >=3) %>% .$cloneCount %>% sum(),
    Bcell_expan_clone_num = dt %>% filter(substr(allVHitsWithScore,1,2)=='IG' & cloneCount >=3) %>% .$cloneCount %>% length() ,
    Tcell_tot_clone_sum = dt %>% filter(substr(allVHitsWithScore,1,2)=='TR') %>% .$cloneCount %>% sum(),
    Tcell_expan_clone_sum = dt %>% filter(substr(allVHitsWithScore,1,2)=='TR' & cloneCount >=3) %>% .$cloneCount %>% sum(),
    Tcell_expan_clone_num = dt %>% filter(substr(allVHitsWithScore,1,2)=='TR' & cloneCount >=3) %>% .$cloneCount %>% length() ,
    stringsAsFactors = F)
}) %>% do.call(rbind,.)

# (meta_dt$id %in% result_tbl2$sampleid) %>% table()
result_tbl2 <- result_tbl2[match(meta_dt$id,result_tbl2$sampleid),]

write_rds(result_tbl2, "data/mixcr.rds")
mixcr_meta_dt <- read_rds("~kjyi/Projects/thymus_single_cell/final2/data/mixcr.rds")


# table(result_tbl2$sampleid==meta_dt$id)
m_meta_dt <- cbind(meta_dt, result_tbl2)

# rownames(ssgsea_dt)
ssgsea_dt <- t(read_rds('~kjyi/Projects/thymus_single_cell/final2/data/immune_ssgsea_all.Rds'))
# table(meta_dt$id %in% rownames(ssgsea_dt))
ssgsea_dt <- ssgsea_dt[m_meta_dt$id,]

m_meta_dt <- cbind(m_meta_dt, ssgsea_dt)

# library size factor
rsem_count <- read_tsv("/home/users/kjyi/Projects/thymus_single_cell/final/expression/IRS4/genes.count.txt")
library_read_counts <- colSums(rsem_count[,-1])
library_read_counts <- data.frame(name=names(library_read_counts),readcount=library_read_counts)
library_read_counts$name = str_replace(library_read_counts$name, "RSEM/","")
library_read_counts$name = str_replace(library_read_counts$name, "-01A-.*","")
library_read_counts$name = str_replace(library_read_counts$name, ".genes.results","")
rownames(library_read_counts) = library_read_counts$name
m_meta_dt$library_size_factor <- library_read_counts[m_meta_dt$id,"readcount"]/1E6
par(pty="s")
boxplot(m_meta_dt$library_size_factor~m_meta_dt$cohort)

ggplot(m_meta_dt,aes(x=`B cell`, y=Bcell_expan_clone_sum))+
  geom_point(aes(color=cohort, size=Bcell_expan_clone_num))

m_meta_dt %>%
ggplot(aes(x=`B cell`, y=Bcell_expan_clone_sum/library_size_factor))+
  geom_point(aes(color=cohort, size=Bcell_expan_clone_num))

ggplot(m_meta_dt, aes(x=Cytotoxic, y=Tcell_expan_clone_sum))+
  geom_point(aes(color=GTF2I_status2, size=Tcell_expan_clone_num))

ggplot(m_meta_dt, aes(x=Cytotoxic, y=Tcell_expan_clone_sum/Tcell_tot_clone_sum))+
  geom_point(aes(color=GTF2I_status2, size=Tcell_expan_clone_num))+
  theme_bw()+theme(panel.grid= element_blank())


ggplot(m_meta_dt,aes(x=`B cell`, y=Bcell_expan_clone_sum))+
  geom_point(aes(color=cohort, size=Bcell_expan_clone_num))

hist(m_meta_dt$Bcell_expan_clone_num,breaks=30)
hist(m_meta_dt$Tcell_expan_clone_num,breaks=30)

m_meta_dt$GTF2I_status2[m_meta_dt$histologic_type %in% "NE"] = NA

m_meta_dt <- m_meta_dt[!m_meta_dt$histologic_type %in% "NE",]



cairo_pdf("~kjyi/Projects/thymus_single_cell/final2/figures/mixcr.pdf",height = 10/2.54,width=19/2.54,pointsize = 12*0.7)
par(mfrow=c(1,2),pty='s', mar=c(4,4,1,1),oma=c(0,1,0,0))
plot(m_meta_dt$Cytotoxic, m_meta_dt$Tcell_expan_clone_sum,
     cex=as.numeric(Hmisc::cut2(m_meta_dt$Tcell_expan_clone_num,cuts = c(0,2,5,7))),
     pch=21,xlab="",ylab="",
     bg=gtf2i_pal[m_meta_dt$GTF2I_status2] %>% adjustcolor(alpha = 0.6))
mtext("Read count of expanded clones\nfrom TCR assembly (RPM)", side = 2,line = 2)
mtext("Cytotoxic T cell enrichment score", side = 1,line = 2)
par(new=T)
plot(100,100,ylim=c(0,1),xlim=c(0,1),xaxt='n',yaxt='n',xlab="",ylab="")
points(rep(0.03,4),c(0.9,0.85,0.79,0.70)-0.02,pch=21,cex=1:4)
text(rep(0.09,4),c(0.9,0.85,0.79,0.70)-0.02,c("0-1","2-4","5-6","7-10"),adj=0)
text(0.09,0.97,"Number of\nclones")
# legend("left",c("0-1","2-4","5-6","7-10"),pch=21,pt.cex=c(1,2,3,4),y.intersp = c(1,1,1.1,1.3),x.intersp = 1.5)

plot(m_meta_dt$`B cell`, m_meta_dt$Bcell_expan_clone_sum,
     cex=as.numeric(Hmisc::cut2(m_meta_dt$Bcell_expan_clone_num,cuts = c(0,100,200,500,1000))),
     pch=21,xlab="",ylab="",
     bg=gtf2i_pal[m_meta_dt$GTF2I_status2] %>% adjustcolor(alpha = 0.6))
text(m_meta_dt$`B cell`, m_meta_dt$Bcell_expan_clone_sum,m_meta_dt$id)
mtext("Read counts of expanded clones\nfrom Ig genes assembly (RPM)", side = 2,line = 2)
mtext("B cell enrichment score", side = 1,line = 2)
par(new=T)
plot(100,100,ylim=c(0,1),xlim=c(0,1),xaxt='n',yaxt='n',xlab="",ylab="")
points(rep(0.03,4),c(0.9,0.85,0.79,0.70)-0.02,pch=21,cex=1:4)
text(rep(0.09,4),c(0.9,0.85,0.79,0.70)-0.02,c("0-99","100-199","200-499",">500"),adj=0)
text(0.09,0.97,"Number of\nclones")


# col=ifelse(m_meta_dt$history_myasthenia_gravis == "YES", "red","black"),

dev.off()



plotly::plot_ly(m_meta_dt, ~`B cell`, ~Bcell_expan_clone_sum,
             cex=~as.numeric(Hmisc::cut2(m_meta_dt$Bcell_expan_clone_num,cuts = c(0,100,200,500,1000))),
             pch=21,xlab="",ylab="",
             bg=~gtf2i_pal[m_meta_dt$GTF2I_status2] %>% adjustcolor(alpha = 0.6))





m_meta_dt %>% arrange(-Bcell_expan_clone_sum) %>% View()


plot(m_meta_dt$`B cell`, m_meta_dt$Bcell_expan_clone_sum,
     cex=as.numeric(Hmisc::cut2(m_meta_dt$Bcell_expan_clone_num,cuts = c(0,100,200,500,1000))),
     pch=21,xlab="",ylab="",
     bg=ifelse(m_meta_dt$history_myasthenia_gravis == "YES", gtf2i_pal[m_meta_dt$GTF2I_status2],"grey40") %>% adjustcolor(alpha = 0.6))
mtext("Read counts of expanded clones\nfrom Ig genes assembly (RPM)", side = 2,line = 2)
mtext("B cell enrichment score", side = 1,line = 2)
par(new=T)
plot(100,100,ylim=c(0,1),xlim=c(0,1),xaxt='n',yaxt='n',xlab="",ylab="")
points(rep(0.03,4),c(0.9,0.85,0.79,0.70)-0.02,pch=21,cex=1:4)
text(rep(0.09,4),c(0.9,0.85,0.79,0.70)-0.02,c("0-99","100-199","200-499",">500"),adj=0)
text(0.09,0.97,"Number of\nclones")




plot(log10(m_meta_dt$corrected_IRS4_TPM+1), m_meta_dt$Tcell_expan_clone_sum,
     cex=as.numeric(Hmisc::cut2(m_meta_dt$Tcell_expan_clone_num,cuts = c(0,2,5,7))),
     pch=21,xlab="IRS4TPM",ylab="T cell clone expansion",
     bg=gtf2i_pal[m_meta_dt$GTF2I_status2] %>% adjustcolor(alpha = 0.6))
plot(log10(m_meta_dt$corrected_IRS4_TPM+1), m_meta_dt$Bcell_expan_clone_sum,
     cex=as.numeric(Hmisc::cut2(m_meta_dt$Bcell_expan_clone_num,cuts = c(0,100,200,500,1000))),
     pch=21,xlab="IRS4TPM",ylab="B cell clone expansion",
     bg=gtf2i_pal[m_meta_dt$GTF2I_status2] %>% adjustcolor(alpha = 0.6),ylim=c(0,20000))



plot(log10(m_meta_dt$corrected_IRS4_TPM+1), m_meta_dt$`B cell`,
     cex=as.numeric(Hmisc::cut2(m_meta_dt$Bcell_expan_clone_num,cuts = c(0,100,200,500,1000))),
     pch=21,xlab="IRS4TPM",ylab="B cell score",
     bg=gtf2i_pal[m_meta_dt$GTF2I_status2] %>% adjustcolor(alpha = 0.6))

plot(m_meta_dt$final_cellularity, m_meta_dt$`B cell`,
     cex=as.numeric(Hmisc::cut2(m_meta_dt$Bcell_expan_clone_num,cuts = c(0,100,200,500,1000))),
     pch=21,xlab="purity",ylab="B cell score",
     bg=gtf2i_pal[m_meta_dt$GTF2I_status2] %>% adjustcolor(alpha = 0.6))

plot(m_meta_dt$final_cellularity, m_meta_dt$`T cell`,
     cex=as.numeric(Hmisc::cut2(m_meta_dt$Tcell_expan_clone_num,cuts = c(0,2,5,7))),
     pch=21,xlab="purity",ylab="B cell score",
     bg=gtf2i_pal[m_meta_dt$GTF2I_status2] %>% adjustcolor(alpha = 0.6))

m_meta_dt[m_meta_dt$GTF2I_status2=="w",] %>% {
  summary(lm(`B cell`~final_cellularity+corrected_IRS4_TPM,data=.))
}


# export for review
plot(m_meta_dt$`B cell`, m_meta_dt$Bcell_expan_clone_sum,
           cex=as.numeric(Hmisc::cut2(m_meta_dt$Bcell_expan_clone_num,cuts = c(0,100,200,500,1000))),
           pch=21,xlab="",ylab="",
           bg=gtf2i_pal[m_meta_dt$GTF2I_status2] %>% adjustcolor(alpha = 0.6))
text(m_meta_dt$`B cell`, m_meta_dt$Bcell_expan_clone_sum,m_meta_dt$id)


library(plotly)

plot_ly(m_meta_dt, x=~`B cell`, y=~Bcell_expan_clone_sum,type="scatter",
        hovertext=~id, size=2*as.numeric(Hmisc::cut2(m_meta_dt$Bcell_expan_clone_num,cuts = c(0,100,200,500,1000))),
        color=gtf2i_pal[m_meta_dt$GTF2I_status2]
     # cex=as.numeric(Hmisc::cut2(m_meta_dt$Bcell_expan_clone_num,cuts = c(0,100,200,500,1000))),
     # pch=21,xlab="",ylab="",
     # bg=gtf2i_pal[m_meta_dt$GTF2I_status2] %>% adjustcolor(alpha = 0.6))
)
htmlwidgets::saveWidget(p, "test.html")
m_meta_dt %>%
  arrange(desc(Bcell_expan_clone_sum)) %>%
  dplyr::select(id,cohort,histologic_type, GTF2I_status2,Bcell_expan_clone_sum,`B cell`) %>%
  write_csv("data/review_template.csv")
