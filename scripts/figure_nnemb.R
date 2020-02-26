# nearest neighbor embedding
# param: w
# k=20
# 

library(tidyverse)
library(Seurat)


gtf2i_pal = c("#3C5488FF", "#E64B35FF", "#00A087FF")
names(gtf2i_pal) = c('m','w','c')
histo_pal = c("#631879FF", "#3B4992FF", "#5F559BFF", "#A20056FF", "#BB0021FF", "#EE0000FF", "#008B45FF", "#008280FF")
names(histo_pal) = c("A","AB","MN-T","B1","B2","B3","NE","TC")

merged <- read_rds("data/singlecell/merged.Rds")
gene.use <- read_rds("data/gene.use.Rds")
load("data/data.scaled.for_comparison.RData")

# merged@reductions$bbknn@assay.used = "RNA"

T0 <- t(thymus_cluster_scaled[thymus_cluster_label %in% c("progenitor", "cTEC", "mTEC", "Tuft", "jTEC"),])
Tx <- t(thymoma_mggene_scaled)

if(T){
  # gene weight --------------------------------------------------------------------------------------
  mg_markers <- merged[gene.use$mgname,] %>% subset(simplified %in% c("cTEC", "progenitor","mTEC","jTEC","Tuft")) %>%
    FindAllMarkers(test.use = "wilcox",logfc.threshold = 0.1, min.pct = 0.1, only.pos = T)
  score1 <- mg_markers %>% 
    group_by(cluster) %>%
    mutate(logq = log10(p_val_adj), combined = -logq*avg_logFC) %>%
    ungroup() %>%
    dplyr::select(gene, combined) %>%
    group_by(gene) %>% summarize(combined = min(combined)) %>%
    as.data.frame %>% column_to_rownames("gene") %>%
    {.[rownames(Tx),]} %>%
    {.[is.na(.)] <- min(.,na.rm = T);.} %>%
    {.[.==Inf] <- NA;.} %>%
    {.[is.na(.)] <- max(.,na.rm = T);.} %>%
    {.^0.1} %>%
    {(. - min(.))/(max(.)-min(.))}
  hist(score1)
  table(score1>0.2) # using only 1768 genes
  score1 = score1/sum(score1)
  
  # non-vectorized, two column only function
  wmean=function(x,w){sum(x*w)/sum(w)}
  wcov=function(x,y,w){sum(w*(x-wmean(x,w))*(y-wmean(y,w)))/sum(w)}
  wcor=function(x,y,w){wcov(x,y,w)/sqrt(wcov(x,x,w)*wcov(y,y,w))}
  wcor2 <- function(A,B=NULL,w){
    ws=sum(w)
    x=t(A)-colSums(A*w)/ws
    y=t(B)-colSums(B*w)/ws
    r=tcrossprod(sweep(x,2,w,'*')/sqrt(rowSums(sweep(x^2,2,w,'*'))),y/sqrt(rowSums(sweep(y^2,2,w,'*'))))
    r
  }
}
# -------------------------------------------------------------------------------------------------

fast_cor <- function(xt,yt=NULL){
  if(is.null(yt)){
    x <- t(xt) - colMeans(xt)
    return(tcrossprod(x / sqrt(rowSums(x ^ 2))))
  } else {
    x <- t(xt) - colMeans(xt)
    y <- t(yt) - colMeans(yt)
    return(tcrossprod(x / sqrt(rowSums(x ^ 2)),y / sqrt(rowSums(y ^ 2))))  
  }
}


# calculate nn using 2453 non-stromal genes
# HiClimR::fastCor(A, nSplit = 1,optBLAS=T)
nn <- fast_cor(T0, Tx) %>% apply(2, function(x){names(sort(x,decreasing = T)[1:20])})

if(T){
  nnw <- wcor2(T0, Tx, score1) %>% apply(2, function(x){names(sort(x,decreasing = T)[1:20])})
}
colnames(Tx)[1]
meta_dt$GTF2I_status2[meta_dt$id==colnames(Tx)[1]]
thymoma_label <- meta_dt$GTF2I_status2[match(colnames(Tx),meta_dt$id)]
thymoma_histol <- meta_dt$histologic_type[match(colnames(Tx),meta_dt$id)]


mycol <- c(w = "#4dd816", m = "#1c54ac", c = "#f64a1d")[thymoma_label]
mycol <- gtf2i_pal[thymoma_label]

# myhis <- thymoma_histol
mycol2 <- c(A  = "#B60A27",
            AB = "#F4756B",
            B1 = "#C0DDEB",
            B2 = "#6F94C4",
            B3 = "#356199",
            TC = "#936650",
            NE = "#484848",
            "MN-T" = "black")[thymoma_histol]
mycol2 <- histo_pal[thymoma_histol]
mypch <- c(w = 21, m = 23, c = 25)[thymoma_label]

background_plot <- function(){
  my_color_palette <- c("#F78981", "#CE425A", "#9D0721", "#0BE2A1", "#20A27B", "#00A1FF", "#0B7DC0", "#AB07FF", "#624B92",
                        "#5100FF", "#002EFC", "#1F30BF", "#282C4D", "#1C2362", "#E38900", "#8E766B", "#715757", "#926650",
                        "#BCBCC2", "#84848C", "#74E74C", "#6FA75A", "#102607", "#F766BF")
  # par(mar = c(2,2,2,2),pty="s")
  
  plot(merged@reductions$bbknn@cell.embeddings, 
       pch = 20,
       cex=1,
       col = my_color_palette[as.numeric(Idents(merged))],
       # bty = 'l',
       bty='n',
       xaxt='n',yaxt='n',xlab="",ylab="",
       xlim=c(2.5,11.5),ylim=c(-9,-0.57),
       asp = 1)
  
  # mtext(side =1 ,"UMAP (batch-balanced)",cex=0.8)
}

# show all points

cairo_pdf("figures/nearestneighborembedding.pdf",height = 53.461/25.4,width=53.461/25.4,pointsize = 12*0.7)
par(mar=rep(0,4),oma=c(0,0,0,0),pty='s')
# background_plot()
# rect(1.5,-10,12.5,1,col="#FFFFFF60")
# box()
# emb=matrix(0,ncol=2,nrow=ncol(nn))
# for(i in 1:ncol(nn)){
#   position <- nn[,i] %>%
#     merged@reductions$bbknn@cell.embeddings[.,] %>%
#     apply(2,median)
#     emb[i,]=position
#     position %>%
#     {points(.[1],.[2], pch = mypch[i], col = "black", bg=mycol[i], cex = 2)}}

# weighted pearson
emb=matrix(0,ncol=2,nrow=ncol(nnw))
for(i in 1:ncol(nnw)){
  position <- nnw[,i] %>%
    merged@reductions$bbknn@cell.embeddings[.,] %>%
    apply(2,median)
  emb[i,]=position
}
rownames(emb) = colnames(nnw)
colnames(emb)=c("x","y")
background_plot()

rect(1.5,-10,12.5,1,col="#FFFFFF50")
box()
points(emb, pch = mypch, col = "black", bg=mycol, cex = 1.8)
legend("bottomright", legend = c(expression(GTF2I^mut), expression(GTF2I^WT), "TC"),
       pch = c(21,23,25), col = "black", pt.bg=gtf2i_pal,pt.cex=1.7,bty="n",horiz = F,cex=0.9)
dev.off()

if("for marker plotting, case selection" == 1){
  background_plot()
  points(emb, pch = mypch, col = "black", bg=mycol, cex = 2)
  gatepoints::fhs(as.data.frame(emb)) %>% paste0('"',.,'"') %>% cat(sep=",")
  handgated <- list(
    progenitors=c("TCGA-4V-A9QU","TCGA-4V-A9QW","TCGA-XU-AAY0","TCGA-ZB-A96B","TCGA-ZC-AAAF"),
    immaturecortical=c("TCGA-4X-A9FC","TCGA-5V-A9RR","TCGA-X7-A8D6","TCGA-X7-A8M8","TCGA-XU-A92U","TCGA-YT-A95E","TCGA-YT-A95H"),
    maturecortical=c("TCGA-ZB-A961","TCGA-XM-A8RD","TCGA-4V-A9QR","TCGA-ZB-A96F","SNU_19_C"),
    immaturemedullary=c("SNU_10_C","TCGA-ZB-A96H","TCGA-ZB-A96Q","TCGA-ZL-A9V6"),
    maturemedullary=c("TCGA-X7-A8DG","TCGA-ZC-AAAH","TCGA-4X-A9FB","TCGA-ZB-A96R","TCGA-3S-AAYX","TCGA-4V-A9QJ"),
    tuft=c("SNU_04_C","TCGA-5U-AB0D","TCGA-X7-A8DD","TCGA-XU-A936","TCGA-ZB-A96A")
  )
}



cairo_pdf("figures/supplimentary.nearestneighborembedding.histol.pdf",height = 10/2.54,width=10/2.54,pointsize = 12*0.7)
background_plot()
box()
rect(1.5,-10,12.5,1,col="#FFFFFF50")
for(i in 1:ncol(nnw)){nnw[,i] %>%
    merged@reductions$bbknn@cell.embeddings[.,] %>%
    apply(2,median) %>%
    {points(.[1],.[2], pch = mypch[i], col = "black", bg=mycol2[i], cex = 2)}}
legend("bottomright",  legend = names(histo_pal),
       pch = 23, pt.bg = histo_pal,ncol = 2,pt.cex=1.2)

# background_plot()
# rect(1.5,-10,12.5,1,col="#FFFFFF60")
# box()
# for(i in 1:ncol(nn)){nn[,i] %>%
#     merged@reductions$bbknn@cell.embeddings[.,] %>%
#     apply(2,median) %>%
#     {points(.[1],.[2], pch = mypch[i], col = "black", bg=mycol2[i], cex = 2)}}
# legend("bottomright",  legend = names(histo_pal),
#        pch = 23, pt.bg = histo_pal,ncol = 2,pt.cex=1.2)
dev.off()




# per-case bootstrapping
# pdf("figures/fig6f_extended.pdf",5.5,5.5)
for(i in 1:ncol(nn)){
  # i = 2
  background_plot()
  nn[,i] %>%{
    merged@reductions$bbknn@cell.embeddings[.,] %>% points()
    merged@reductions$bbknn@cell.embeddings[.,] %>% apply(2,median) %>%
    {points(.[1],.[2], pch = mypch[i], col = "black",bg=mycol2[i], cex = 2)}
  }
  library(doMC)
  x <- foreach(1:100, .combine = rbind) %do% {
    boot <- sample(nrow(T0),nrow(T0),replace = T) #2452
    cor(T0[boot,], Tx[boot,i])[,1] %>% sort(decreasing = T) %>% .[1:25] %>% names %>%
      merged@reductions$bbknn@cell.embeddings[.,] %>% apply(2,median)
  } 
  myrank <- ((x[,1] - median(x[,1])) ^ 2 + (x[,2] - median(x[,2])) ^ 2 ) ^ 0.5 %>% rank
  ch <- chull(x[which(myrank <= 90),])
  ch <- c(ch, ch[1])
  lines(x[which(myrank <= 90),][ch,], col = mycol2[i])
  title(colnames(nn)[i])
}
dev.off()
