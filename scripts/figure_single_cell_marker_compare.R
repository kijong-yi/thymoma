library(ComplexHeatmap)
load("data/snapshop_for_heatmap_comparison.RData")

T0 <- t(thymus_cluster_scaled[thymus_cluster_label %in% c("progenitor", "cTEC", "mTEC", "Tuft", "jTEC"),])
Tx <- t(thymoma_mggene_scaled)


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

if(F){
  wcor2 <- function(A,B=NULL,w){ws=sum(w);x=t(A)-colSums(A*w)/ws;y=t(B)-colSums(B*w)/ws;
  tcrossprod(sweep(x,2,w,'*')/sqrt(rowSums(sweep(x^2,2,w,'*'))),y/sqrt(rowSums(sweep(y^2,2,w,'*'))))}
  nn <- fast_cor(T0, Tx) %>% apply(2, function(x){names(sort(x,decreasing = T)[1:20])})
  nnw <- wcor2(T0, Tx, score1) %>% apply(2, function(x){names(sort(x,decreasing = T)[1:20])})
}

Tp = Tx
for(i in 1:ncol(Tp)){
  Tp[,i] = rowSums(T0[,nn[,i]])
}
dim(Tp)
dim(Tx)
Tp == Tx
r=lapply(1:nrow(Tp),function(i){cor(Tp[i,],Tx[i,])}) %>% do.call(c,.)
r2=lapply(1:nrow(Tp),function(i){cor(Tp[i,],Tx[i,],method = "spearman")}) %>% do.call(c,.)
which(apply(Tp,1,var)==0)
which(apply(Tp,1,mean)==0)


which(is.na(r))
Tp[which(is.na(r)),] %>% apply(1,sd) # sd==0 in pseudobulk..
which(is.na(r))
which(is.na(r2))
names(r) = rownames(T0)
names(r2) = rownames(T0)
sort(r,decreasing=T)
sort(r2,decreasing=T)
# # A tibble: 3 x 2
#   mgname  hgname
#   <chr>   <chr> 
# 1 Ugt1a10 UGT1A8
# 2 Ugt1a9  UGT1A8
# 3 Ugt1a8  UGT1A8

handgated <- list(progenitors=c("TCGA-4V-A9QU","TCGA-4V-A9QW","TCGA-X7-A8DF","TCGA-XU-AAY0","TCGA-ZB-A96B","TCGA-ZC-AAAF"),
                  maturecortical=c("TCGA-4X-A9FD","TCGA-ZB-A961","TCGA-ZB-A96F","SNU_17_C","SNU_21_C","SNU_19_C"),
                  maturemedullary=c("TCGA-X7-A8DG","TCGA-XU-AAXY","TCGA-3Q-A9WF","TCGA-ZC-AAAH","TCGA-4X-A9FB","TCGA-ZB-A96O","TCGA-3T-AA9L","TCGA-ZB-A96R","TCGA-XM-A8RL","TCGA-3G-AB0O","TCGA-X7-A8DB","TCGA-XM-A8RB","TCGA-5U-AB0F","TCGA-X7-A8M6","TCGA-4V-A9QT","TCGA-3S-AAYX","TCGA-XM-A8RE","SNU_11_C","TCGA-4V-A9QJ"),
                  tuft=c("SNU_04_C","TCGA-XU-A936","TCGA-ZB-A96A"),
                  immaturecortical=c("SNU_27_C","TCGA-4V-A9QM","TCGA-4V-A9QX","TCGA-X7-A8D6","TCGA-X7-A8M8","TCGA-XM-A8RH","TCGA-XU-A92U","TCGA-YT-A95E"),
                  medullarylow=c("TCGA-XU-AAXZ","SNU_22_C","TCGA-3S-A8YW","TCGA-5G-A9ZZ","TCGA-X7-A8DD","TCGA-4V-A9QN","TCGA-ZT-A8OM"),
                  immaturecortical2=c("SNU_28_C","TCGA-X7-A8D8","TCGA-X7-A8DJ","TCGA-X7-A8M1","TCGA-XM-AAZ2","SNU_14_C","TCGA-ZB-A96C"),
                  innerjunctional=c("TCGA-XM-A8RI","TCGA-XM-A8RD","TCGA-X7-A8M0","TCGA-X7-A8M5","TCGA-XU-A92X","TCGA-XM-A8RC","TCGA-ZB-A963"),
                  outerjunctional=c("TCGA-ZB-A96I","TCGA-XU-A92Q","TCGA-X7-A8M3"))
Tp[1:3,1:3]
Tx[1:3,1:3]
handgated$progenitors
handgated[1:4]
hist(r2)
table(r2>0.4)


which(Tx[names(sort(r2,decreasing=T)[1:500]),unlist(handgated)]%>% apply(1,var) == 0)
which(Tp[names(sort(r2,decreasing=T)[1:500]),unlist(handgated)]%>% apply(1,var) == 0)

r2 = r2[names(r2) != "Gm24682"]

Heatmap(Tx[names(sort(r2,decreasing=T)[100:500]),unlist(handgated)],cluster_columns = F, cluster_rows = T,
        clustering_method_rows = "ward.D",clustering_distance_rows = "pearson") + 
  Heatmap(Tp[names(sort(r2,decreasing=T)[100:500]),unlist(handgated)],cluster_columns = F, cluster_rows = T,
          clustering_method_rows = "ward.D",clustering_distance_rows = "pearson")

Tp[names(sort(r2,decreasing=T)[1:500]),unlist(handgated[1:4])]
rank(c(100,20,500,1))
getmarkers <- function(Tx,Tp,query,backgroud=handgated[1:4],n=50){
  genes = rownames(Tx)
  Rx = t(apply(-Tx,1,rank))
  Rp = t(apply(-Tp,1,rank))
  Sx=rowMeans(Rx[,query])
  Sp=rowMeans(Rp[,query])
  mS=sqrt((Sx^2)*Sp)
  x=Heatmap(Tx[names(sort(mS)[1:n]),unlist(backgroud)],
          cluster_columns = F, cluster_rows = T,
          clustering_method_rows = "ward.D",
          clustering_distance_rows = "pearson") + 
    Heatmap(Tp[names(sort(mS)[1:n]),unlist(backgroud)],
            cluster_columns = F, cluster_rows = T,
            clustering_method_rows = "ward.D",
            clustering_distance_rows = "pearson")
  draw(x)
  sort(mS)[1:n]
}
handgated
handgated <- list(
  progenitors=c("TCGA-4V-A9QU","TCGA-4V-A9QW","TCGA-XU-AAY0","TCGA-ZB-A96B","TCGA-ZC-AAAF"),
  immaturecortical=c("TCGA-4X-A9FC","TCGA-5V-A9RR","TCGA-X7-A8D6","TCGA-X7-A8M8","TCGA-XU-A92U","TCGA-YT-A95E","TCGA-YT-A95H"),
  immaturemedullary=c("SNU_10_C","TCGA-ZB-A96H","TCGA-ZB-A96Q","TCGA-ZL-A9V6"),
  maturemedullary=c("TCGA-X7-A8DG","TCGA-ZC-AAAH","TCGA-4X-A9FB","TCGA-ZB-A96R","TCGA-3S-AAYX","TCGA-4V-A9QJ"),
  maturecortical=c("TCGA-ZB-A961","TCGA-XM-A8RD","TCGA-4V-A9QR","TCGA-ZB-A96F","SNU_19_C"),
  tuft=c("SNU_04_C","TCGA-X7-A8DD","TCGA-ZB-A96A","TCGA-5U-AB0D","TCGA-XU-A936")
)

cairo_pdf("tmp.pdf",height = 21/2.54,width=29.7/2.54,pointsize = 12*0.7,onefile = T)
getmarkers(Tx,Tp,handgated$progenitors,handgated,30) %>% names %>% cat(sep="\n")
getmarkers(Tx,Tp,handgated$maturecortical,handgated,50) %>% names %>% cat(sep="\n")
getmarkers(Tx,Tp,handgated$maturemedullary,handgated,50) %>% names %>% cat(sep="\n")
getmarkers(Tx,Tp,handgated$tuft,handgated,50) %>% names %>% cat(sep="\n")
getmarkers(Tx,Tp,handgated$immaturecortical,handgated,50) %>% names %>% cat(sep="\n")
getmarkers(Tx,Tp,handgated$immaturemedullary,handgated,50) %>% names %>% cat(sep="\n")
dev.off()

system("rm tmp.pdf")

epi <- !(merged$simplified %in% c("stroma","etc"))


plotfx <- function(x=handgated$progenitors){
  plot(merged@reductions$bbknn@cell.embeddings[epi,], 
       pch = 20,
       cex=1,
       col = "grey",
       bty='n',
       xaxt='n',yaxt='n',xlab="",ylab="",
       xlim=c(2.5,11.5),ylim=c(-9,-0.57),
       asp = 1)
  points(emb[!rownames(emb) %in% x,], pch = 23, col = "#888888C0",
         bg="#999999C0", cex = 2)
  points(emb[rownames(emb) %in% x,], pch = 23, col = "black",
         bg="red", cex = 2)
}
nnw[,handgated$progenitors]
plotfx2 <- function(x=handgated$progenitors){
  plot(merged@reductions$bbknn@cell.embeddings[epi,], 
       pch = 20,
       cex=1,
       col = "grey",
       bty='n',
       xaxt='n',yaxt='n',xlab="",ylab="",
       xlim=c(2.5,11.5),ylim=c(-9,-0.57),
       asp = 1)
  points(merged@reductions$bbknn@cell.embeddings[unique(unlist(nnw[,x])),],col="#FF000050",pch=20)
  # points(emb[!rownames(emb) %in% x,], pch = 23, col = "#888888C0",
  #        bg="#999999C0", cex = 2)
  # points(emb[rownames(emb) %in% x,], pch = 23, col = "black",
  #        bg="red", cex = 2)
}
par(mar=c(0,0,0,0),mfrow=c(5,10))
plotfx(handgated$progenitors)
plotfx(handgated$immaturecortical)
plotfx(handgated$maturecortical)
plotfx(handgated$maturemedullary)
plotfx(handgated$tuft)


plotfx2(handgated$progenitors)
plotfx2(handgated$immaturecortical)
plotfx2(handgated$maturecortical)
plotfx2(handgated$maturemedullary)
plotfx2(handgated$tuft)

cairo_pdf("figures/figure_umap_above_heatmap.pdf",height = 7/2.54,width=30/2.54,pointsize = 12*0.7)

par(mar=c(0,0,0,0),mfrow=c(1,2), pty='s')
plot(merged@reductions$bbknn@cell.embeddings[epi,],pch = 20,cex=1,
      col = "grey",bty='n',xaxt='n',yaxt='n',xlab="",ylab="",
      xlim=c(2.5,11.5),ylim=c(-9,-0.57),asp = 1)
# points(emb[!rownames(emb) %in% handgated$progenitors,], pch = 23, col = "#888888C0",bg="#999999C0", cex = 2)
points(emb[!rownames(emb) %in% unlist(handgated),], pch = 23, col = "#888888C0",bg="#999999C0", cex = 2)
points(emb[rownames(emb) %in% handgated$progenitors,], pch = 23, col = "black",bg="#0bc187", cex = 2)
points(emb[rownames(emb) %in% handgated$immaturecortical,], pch = 23, col = "black",bg="#ff91ad", cex = 2)
points(emb[rownames(emb) %in% handgated$maturecortical,], pch = 23, col = "black",bg="#d70c24", cex = 2)
points(emb[rownames(emb) %in% handgated$immaturemedullary,], pch = 23, col = "black",bg="#33a7ee", cex = 2)
points(emb[rownames(emb) %in% handgated$maturemedullary,], pch = 23, col = "black",bg="#053061", cex = 2)
points(emb[rownames(emb) %in% handgated$tuft,], pch = 23, col = "black",bg="#f18821", cex = 2)


plot(merged@reductions$bbknn@cell.embeddings[epi,],pch = 20,cex=1,
     col = "grey",bty='n',xaxt='n',yaxt='n',xlab="",ylab="",
     xlim=c(2.5,11.5),ylim=c(-9,-0.57),asp = 1)
points(merged@reductions$bbknn@cell.embeddings[unique(unlist(nnw[,handgated$progenitors])),],col="#0bc187",pch=20)
points(merged@reductions$bbknn@cell.embeddings[unique(unlist(nnw[,handgated$immaturecortical])),],col="#ff91ad",pch=20)
points(merged@reductions$bbknn@cell.embeddings[unique(unlist(nnw[,handgated$maturecortical])),],col="#d70c24",pch=20)
points(merged@reductions$bbknn@cell.embeddings[unique(unlist(nnw[,handgated$immaturemedullary])),],col="#33a7ee",pch=20)
points(merged@reductions$bbknn@cell.embeddings[unique(unlist(nnw[,handgated$maturemedullary])),],col="#053061",pch=20)
points(merged@reductions$bbknn@cell.embeddings[unique(unlist(nnw[,handgated$tuft])),],col="#f18821",pch=20)
dev.off()


markers <- list(
  progenitors = c("Elob","C1qtnf12","Tmem59l"),
  immaturecortical=c("Irx1","Id1","Id3","Spock3"),
  maturecortical = c("Rab7b","Lgals7","Gpt2","Aldh1l1","Tsku","Parp3","Zbp1"),
  immaturemedullary=c("Tbx1","Erbb2","Tnfaip2","Jag1","Zfp750","Esrp1","Cdh1","Cldn8","F3","Myo5b","App"),
  maturemedullary=c("Aire","Aif1","Cytl1","Marco","Tff3","Gzmf","Hmcn2"),
  tuft=c("Neurl1a","Ovol3","Trpm5","Sh2d6","Gnb3","Hepacam2","Pou2f3","Avil"),
  other=c("Irs4")
)

if(F){save(markers,handgated,Tx,T0,nnw,emb,merged,file = "data/forfigure_single_cell_marker_comparison.RData",compress = "gzip")}

unlist(markers) %>% length
gene.use[match(unlist(markers),gene.use$mgname),]$hgname %>% length


gene.group <- names(unlist(markers)) %>% str_replace("[0-9]*$","") %>% 
{c("tuft"="Tuft-like",
   "maturemedullary"="Mature medullary",
   "immaturecortical" = "Developing cortical",
   "progenitors" = "Progenitor",
   "immaturemedullary" = "Developing medullary",
   "maturecortical" = "Mature cortical",
   "other"="Others")[.]} %>% unname() %>% 
  factor(levels=c("Progenitor","Developing cortical","Mature cortical","Developing medullary","Mature medullary","Tuft-like","Others"))

tumor.group <- rep(names(handgated),unlist(lapply(handgated,length))) %>%
{c("tuft"="Tuft-like",
   "maturemedullary"="Mature medullary",
   "immaturecortical" = "Developing cortical",
   "progenitors" = "Progenitor",
   "immaturemedullary" = "Developing medullary",
   "maturecortical" = "Mature cortical")[.]} %>% unname() %>% 
  factor(levels=c("Progenitor","Developing cortical","Mature cortical","Developing medullary","Mature medullary","Tuft-like"))
tumor.genotype <- meta_dt$GTF2I_status2[match(unlist(handgated),meta_dt$id)]
gtf2i_pal = c(m="#3C5488FF",w="#E64B35FF",c="#00A087FF")
tumor.genotype.ha <- HeatmapAnnotation("Group"= tumor.genotype,
                                       col = list("Group" = gtf2i_pal))



cell.group <- rev(unique(rev(c((nnw[,unlist(handgated)])))))
names(cell.group) = cell.group
cell.group[c(nnw[,unlist(handgated$progenitors)])] <- "Progenitor"
cell.group[c(nnw[,unlist(handgated$immaturecortical)])] <- "Developing cortical"
cell.group[c(nnw[,unlist(handgated$maturecortical)])] <- "Mature cortical"
cell.group[c(nnw[,unlist(handgated$immaturemedullary)])] <- "Developing medullary"
cell.group[c(nnw[,unlist(handgated$maturemedullary)])] <- "Mature medullary"
cell.group[c(nnw[,unlist(handgated$tuft)])] <- "Tuft-like"
table(cell.group)
cell.group <- factor(cell.group,levels=c("Progenitor","Developing cortical","Mature cortical","Developing medullary","Mature medullary","Tuft-like"))

mycolors=c("#0bc187",#progenitor
           "#ff91ad",#ic
           "#d70c24",#mc
           "#33a7ee",#im
           "#053061",#mm
           "#f18821",#tuft
           "grey50")#other
exp_pal = circlize::colorRamp2(c(-4,0,4), c('#253494',"gray90",'#f03b20'))

tmpmat <- Tx[unlist(markers),unlist(handgated)]
rownames(tmpmat) <- gene.use[match(unlist(markers),gene.use$mgname),]$hgname
Heatmap(tmpmat, cluster_rows = F,cluster_columns = F,show_heatmap_legend = F,col=exp_pal,
        width=unit(10,"cm"),height=unit(18,"cm"),row_names_side = "left",show_column_names = F,
        row_split=gene.group,column_split=tumor.group,row_title =" ",column_title=" ",
        bottom_annotation = HeatmapAnnotation(Group=tumor.genotype,col=list(Group=gtf2i_pal),
                                              annotation_legend_param = list(
                                                "Group" = list(
                                                  title = "Group",
                                                  at = c("m", "w","c"),
                                                  labels = c("GTF2I-mutant", "Wild-type","Thymic carcinoma")))),
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = mycolors)))) +
T0[unlist(markers),names(cell.group)] %>% t() %>% scale %>% t() %>%
Heatmap(cluster_rows = F,cluster_columns = F,column_title=" ",col=exp_pal,
        width=unit(10,"cm"),height=unit(18,"cm"),show_column_names = F,
        column_split = cell.group,name="Normalized gene expression",
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = mycolors))),
        heatmap_legend_param = list(direction = "horizontal")) -> hml
draw(hml,merge_legend = TRUE,heatmap_legend_side = "bottom")

cairo_pdf("figures/figure_heatmap_tumor_sc_comparison.pdf",height = 25/2.54,width=30/2.54,pointsize = 12*0.7)
draw(hml,merge_legend = TRUE,heatmap_legend_side = "bottom")
dev.off()
if(F){save.image(file="data/snapshop_for_heatmap_comparison.RData",compress = "gzip")}
load("data/snapshop_for_heatmap_comparison.RData")
