library(velocyto.R)
library(tidyverse)
library(Seurat)

# loom files from three dataset

# old mouse
oldmouse.loomfile <- "/home/users/team_projects/thymus_single_cell/mouse_tec/cellranger/TEC/velocyto/TEC.loom"

# kernfeld et al.
kernfeld.loomfile <- list(E12_5_wholeThy_venus_1 = "~kjyi/Projects/thymus_single_cell/Kernfeld/raw/processed/E12_5_wholeThy_venus_1.velocyto_result/E12_5_wholeThy_venus_1_QRYEI.loom",
                          E12_5_wholeThy_venus_2 = "~kjyi/Projects/thymus_single_cell/Kernfeld/raw/processed/E12_5_wholeThy_venus_2.velocyto_result/E12_5_wholeThy_venus_2_HIJL0.loom",
                          E12_5_wholeThy_venus_3 = "~kjyi/Projects/thymus_single_cell/Kernfeld/raw/processed/E12_5_wholeThy_venus_3.velocyto_result/E12_5_wholeThy_venus_3_3J3ZX.loom",
                          E13_5_wholeThy_1 = "~kjyi/Projects/thymus_single_cell/Kernfeld/raw/processed/E13_5_wholeThy_1.velocyto_result/E13_5_wholeThy_1_EEHDV.loom",
                          E13_5_wholeThy_2 = "~kjyi/Projects/thymus_single_cell/Kernfeld/raw/processed/E13_5_wholeThy_2.velocyto_result/E13_5_wholeThy_2_WCE9O.loom",
                          E13_5_wholeThy_3 = "~kjyi/Projects/thymus_single_cell/Kernfeld/raw/processed/E13_5_wholeThy_3.velocyto_result/E13_5_wholeThy_3_I0BEX.loom",
                          E14_5_wholeThy_1 = "~kjyi/Projects/thymus_single_cell/Kernfeld/raw/processed/E14_5_wholeThy_1.velocyto_result/E14_5_wholeThy_1_MXDLN.loom",
                          E14_5_wholeThy_2 = "~kjyi/Projects/thymus_single_cell/Kernfeld/raw/processed/E14_5_wholeThy_2.velocyto_result/E14_5_wholeThy_2_FH00N.loom",
                          E15_5_wholeThy_1 = "~kjyi/Projects/thymus_single_cell/Kernfeld/raw/processed/E15_5_wholeThy_1.velocyto_result/E15_5_wholeThy_1_7H61D.loom",
                          E15_5_wholeThy_2 = "~kjyi/Projects/thymus_single_cell/Kernfeld/raw/processed/E15_5_wholeThy_2.velocyto_result/E15_5_wholeThy_2_5ZHKQ.loom",
                          E16_5_wholeThy_1 = "~kjyi/Projects/thymus_single_cell/Kernfeld/raw/processed/E16_5_wholeThy_1.velocyto_result/E16_5_wholeThy_1_3EY5J.loom",
                          E16_5_wholeThy_2 = "~kjyi/Projects/thymus_single_cell/Kernfeld/raw/processed/E16_5_wholeThy_2.velocyto_result/E16_5_wholeThy_2_XX59A.loom",
                          E16_5_wholeThy_3 = "~kjyi/Projects/thymus_single_cell/Kernfeld/raw/processed/E16_5_wholeThy_3.velocyto_result/E16_5_wholeThy_3_40NCT.loom",
                          E17_5_wholeThy_1 = "~kjyi/Projects/thymus_single_cell/Kernfeld/raw/processed/E17_5_wholeThy_1.velocyto_result/E17_5_wholeThy_1_VIWT3.loom",
                          E17_5_wholeThy_2 = "~kjyi/Projects/thymus_single_cell/Kernfeld/raw/processed/E17_5_wholeThy_2.velocyto_result/E17_5_wholeThy_2_S71AU.loom",
                          E18_5_wholeThy_1 = "~kjyi/Projects/thymus_single_cell/Kernfeld/raw/processed/E18_5_wholeThy_1.velocyto_result/E18_5_wholeThy_1_MCVP9.loom",
                          E18_5_wholeThy_2 = "~kjyi/Projects/thymus_single_cell/Kernfeld/raw/processed/E18_5_wholeThy_2.velocyto_result/E18_5_wholeThy_2_YC96D.loom",
                          P0_wholeThy_1 = "~kjyi/Projects/thymus_single_cell/Kernfeld/raw/processed/P0_wholeThy_1.velocyto_result/P0_wholeThy_1_2P28W.loom",
                          P0_wholeThy_2 = "~kjyi/Projects/thymus_single_cell/Kernfeld/raw/processed/P0_wholeThy_2.velocyto_result/P0_wholeThy_2_KPUU6.loom")

# ibarra 
ibarra.loomfile <- list(E_MTAB_6153_rep1_1 = "/home/users/kjyi/Projects/thymus_single_cell/final/raw/ibarra/E-MTAB-6153/E_MTAB_6153_rep1_2_cr3/velocyto/E_MTAB_6153_rep1_1_cr3.loom",
                        E_MTAB_6153_rep1_2 = "/home/users/kjyi/Projects/thymus_single_cell/final/raw/ibarra/E-MTAB-6153/E_MTAB_6153_rep1_2_cr3/velocyto/E_MTAB_6153_rep1_2_cr3.loom",
                        E_MTAB_6153_rep1 = "/home/users/kjyi/Projects/thymus_single_cell/final/raw/ibarra/E-MTAB-6153/E_MTAB_6153_rep1_cr3/velocyto/E_MTAB_6153_rep1_cr3.loom",
                        E_MTAB_6153_rep2 = "/home/users/kjyi/Projects/thymus_single_cell/final/raw/ibarra/E-MTAB-6153/E_MTAB_6153_rep2_cr3/velocyto/E_MTAB_6153_rep2_cr3.loom")

# oldmouse data can be processed directly

oldmouse.emb <- read_rds("/home/users/kjyi/Projects/thymus_single_cell/final2/data/singlecell/10x.Rds") %>% Embeddings("tsne")
oldmouse.loom <- read.loom.matrices(file = oldmouse.loomfile,engine='h5')
# str(oldmouse.loom) # cellnames : "TEC:ACTATATATATAx" ...
# str(tmp.embedding) # cellnames : "ACTATATATATA"
oldmouse.emb <- read_rds("/home/users/kjyi/Projects/thymus_single_cell/final2/data/singlecell/10x.Rds") %>% Embeddings("tsne") %>% {rownames(.) <- paste0("TEC:",rownames(.),"x");.}


oldmouse.E <- oldmouse.loom$spliced[,rownames(oldmouse.emb)]
oldmouse.N <- oldmouse.loom$unspliced[,rownames(oldmouse.emb)]

# will take time
oldmouse.velocityestimates <- gene.relative.velocity.estimates(oldmouse.E,oldmouse.N,deltaT=1,kCells=10,n.cores=10)
show.velocity.on.embedding.cor(oldmouse.emb,oldmouse.velocityestimates,n=200,scale='sqrt',cex=0.8,arrow.scale=3,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=F,cell.border.alpha = 0.1)
colnames(oldmouse.E) <- colnames(oldmouse.E) %>% stringr::str_replace(".*:","") %>% stringr::str_replace("x","")
colnames(oldmouse.N) <- colnames(oldmouse.N) %>% stringr::str_replace(".*:","") %>% stringr::str_replace("x","")
write_rds(list(oldmouse.E = oldmouse.E, oldmouse.N = oldmouse.N), "data/singlecell/velocyto.oldmouse.matched.loom.Rds", compress = "gz")
x=read_rds("data/singlecell/velocyto.oldmouse.matched.loom.Rds")
oldmouse.E = x[[1]]
oldmouse.N = x[[2]]
rm(x)

# Kernfeld.... .....

# tmp <- read_rds("single_cell/data/kernfeld.Rds") %>% Embeddings("tsne")
# tmp2 <- read_rds("single_cell/data/kernfeld.Rds")
kernfeld.emb <- read_rds("data/singlecell/kernfeld.Rds") %>% Embeddings("tsne") # for cell name
# rownames(kernfeld.emb) %in% colnames(kernfeld.E) %>% table
# check rownames match, colnames not duplicated, then merge
# tmp_rownames <- read.loom.matrices(kernfeld.loomfile[[1]])[[1]] %>% rownames
# for(i in kernfeld.loomfile){print(all(tmp_rownames == rownames(read.loom.matrices(i)[[1]])))}
kernfeld.loom <- list() # initiation
for(i in 1:length(kernfeld.loomfile)){kernfeld.loom[[names(kernfeld.loomfile)[i]]] <- read.loom.matrices(kernfeld.loomfile[[i]],engine='h5')}
rm(i)
# tmp_colnames <- c()
# for(i in kernfeld.loom){ tmp_colnames = c(tmp_colnames,colnames(i[[1]]))}
# tmp_colnames %>% stringr::str_replace(".*:", "") -> x
# table(duplicated(x))
# x %>% make.unique %>% stringr::str_replace("x", "") %>% {.[stringr::str_length(.)==14]} # good
# rownames(kernfeld.emb) %in% colnames(kernfeld.E) %>% table
# tmp2[[]][!rownames(kernfeld.emb) %in% colnames(kernfeld.E),]


kernfeld.E <- cbind(kernfeld.loom[[1]][[1]], kernfeld.loom[[2]][[1]], kernfeld.loom[[3]][[1]], kernfeld.loom[[4]][[1]], kernfeld.loom[[5]][[1]], kernfeld.loom[[6]][[1]], kernfeld.loom[[7]][[1]], kernfeld.loom[[8]][[1]], kernfeld.loom[[9]][[1]], kernfeld.loom[[10]][[1]], kernfeld.loom[[11]][[1]], kernfeld.loom[[12]][[1]], kernfeld.loom[[13]][[1]], kernfeld.loom[[14]][[1]], kernfeld.loom[[15]][[1]], kernfeld.loom[[16]][[1]], kernfeld.loom[[17]][[1]], kernfeld.loom[[18]][[1]], kernfeld.loom[[19]][[1]]) %>% 
{colnames(.) <- colnames(.) %>% stringr::str_replace(".*:", "") %>% stringr::str_replace("x", "") %>% make.unique();.} 
kernfeld.N <- cbind(kernfeld.loom[[1]][[2]], kernfeld.loom[[2]][[2]], kernfeld.loom[[3]][[2]], kernfeld.loom[[4]][[2]], kernfeld.loom[[5]][[2]], kernfeld.loom[[6]][[2]], kernfeld.loom[[7]][[2]], kernfeld.loom[[8]][[2]], kernfeld.loom[[9]][[2]], kernfeld.loom[[10]][[2]], kernfeld.loom[[11]][[2]], kernfeld.loom[[12]][[2]], kernfeld.loom[[13]][[2]], kernfeld.loom[[14]][[2]], kernfeld.loom[[15]][[2]], kernfeld.loom[[16]][[2]], kernfeld.loom[[17]][[2]], kernfeld.loom[[18]][[2]], kernfeld.loom[[19]][[2]]) %>% 
{colnames(.) <- colnames(.) %>% stringr::str_replace(".*:", "") %>% stringr::str_replace("x", "") %>% make.unique();.}
# some cells in embeding are missing in velocyto loom files. draw arrows without these ....
# FALSE  TRUE 
# 2373 10488 
kernfeld.intersection_cellnames <- intersect(rownames(kernfeld.emb), colnames(kernfeld.E))
kernfeld.vel_not_estimated_cellnames <- rownames(kernfeld.emb)[!rownames(kernfeld.emb) %in% colnames(kernfeld.E)]

kernfeld.E <- kernfeld.E[,kernfeld.intersection_cellnames]
kernfeld.N <- kernfeld.N[,kernfeld.intersection_cellnames]
kernfeld.emb.intersect <- kernfeld.emb[kernfeld.intersection_cellnames,]

write_rds(list(kernfeld.E = kernfeld.E, kernfeld.N = kernfeld.N), "data/singlecell/velocyto.kernfeld.matched.loom.Rds", compress = "gz")
x = read_rds("data/singlecell/velocyto.kernfeld.matched.loom.Rds")
kernfeld.E=x[[1]]
kernfeld.N=x[[2]]
rm(x)
# kernfeld.velocityestimates <- gene.relative.velocity.estimates(kernfeld.E,kernfeld.N,deltaT=1,kCells=10,n.cores=11)
# show.velocity.on.embedding.cor(kernfeld.emb.intersect,kernfeld.velocityestimates,n=200,scale='sqrt',cex=0.8,arrow.scale=3,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=F,cell.border.alpha = 0.1)


# oldmouse = match, kernfeld = almost match ....
# but in case of ibarra ... , we have to find nearest cells ... !!!!!!!!!!!!!!!!!!!!!!!!!
ibarra.fp.seurat <- read_rds("data/singlecell/ibarra_foregut_and_pharyngeal_mesoderm.Rds")
ibarra.all.seurat <- read_rds("data/singlecell/ibarra.Rds") # cell names are kind but horrible to match with raw data
ibarra.count.seurat <- GetAssayData(ibarra.all.seurat, slot = "counts")
ibarra.rep1.loom <- read.loom.matrices(ibarra.loomfile[[3]],engine='h5')
ibarra.rep2.loom <- read.loom.matrices(ibarra.loomfile[[4]],engine='h5')


ibarra.count.loom <- cbind(ibarra.rep1.loom[[1]] + ibarra.rep1.loom[[2]] + ibarra.rep1.loom[[3]],
                           ibarra.rep2.loom[[1]] + ibarra.rep2.loom[[2]] + ibarra.rep2.loom[[3]]) 
ibarra.intersect.genes <- intersect(rownames(ibarra.count.loom),rownames(ibarra.count.seurat))
ibarra.count.loom <- ibarra.count.loom[ibarra.intersect.genes, ]
ibarra.count.seurat <- ibarra.count.seurat[ibarra.intersect.genes, ]
write_rds(ibarra.count.loom, "data/singlecell/velocyto.ibarra.count.loom.Rds")
write_rds(ibarra.count.seurat, "data/singlecell/velocyto.ibarra.count.seurat.Rds")

ibarra.count.loom   <- read_rds("data/singlecell/velocyto.ibarra.count.loom.Rds")
ibarra.count.seurat <- read_rds("data/singlecell/velocyto.ibarra.count.seurat.Rds")
# nn_idx <- read_tsv("nn_idx.corr.txt", col_names = "dist")$dist # calculate corr with julia.. length = 20819
# table(duplicated(nn_idx))

dim(ibarra.count.loom) #26347
dim(ibarra.count.seurat) #20819
ibarra.count.seurat[1:3,1:3]
ibarra.count.loom[1:3,1:3]
# for fast calculation of correlation, I used only 2000 genes
topgenes <- rowSums(ibarra.count.seurat) %>% {rownames(ibarra.count.seurat)[rank(-.) < 2000]}
# sum(ibarra.count.seurat["Mrpl15",]) # good
cormat <- qlcMatrix::corSparse(ibarra.count.seurat[topgenes,],ibarra.count.loom[topgenes,])
write_rds(cormat,"data/singlecell/velocyto.ibarra.cormat.Rds")
cormat[1:3,1:3]
dim(cormat)
rownames(cormat) = colnames(ibarra.count.seurat)
colnames(cormat) = colnames(ibarra.count.loom)
cormat[1:3,1:3]
nn_idx=apply(cormat,1,which.max)
length(nn_idx)
# ibarra.count.match <- ibarra.count.loom[,nn_idx]
# ibarra.count.match[1:3,1:3]
# colnames(ibarra.count.match) <- colnames(ibarra.count.seurat)[nn_idx]
ibarra.E <- cbind(ibarra.rep1.loom[[1]], ibarra.rep2.loom[[1]])[,nn_idx]
ibarra.N <- cbind(ibarra.rep1.loom[[2]], ibarra.rep2.loom[[2]])[,nn_idx]
colnames(ibarra.E) <- colnames(ibarra.count.seurat)
colnames(ibarra.N) <- colnames(ibarra.count.seurat)
dim(ibarra.E)

ibarra.cellname <- intersect(colnames(ibarra.fp.seurat),colnames(ibarra.E))
ibarra.E <- ibarra.E[,ibarra.cellname]
ibarra.N <- ibarra.N[,ibarra.cellname]



write_rds(list(ibarra.E = ibarra.E, ibarra.N = ibarra.N), "data/singlecell/ibarra.matched.loom.Rds", compress = "gz")
x=read_rds("data/singlecell/ibarra.matched.loom.Rds")
ibarra.E = x[[1]]
ibarra.N = x[[2]]


oldmouse.E[1:3,1:3]
kernfeld.E[1:3,1:3]
ibarra.E[1:3,1:3]

intersect_genes <- intersect(rownames(oldmouse.E),rownames(kernfeld.E)) %>% intersect(rownames(ibarra.E))

merged.E <- cbind(oldmouse.E[intersect_genes,],kernfeld.E[intersect_genes,],ibarra.E[intersect_genes,])
merged.N <- cbind(oldmouse.N[intersect_genes,],kernfeld.N[intersect_genes,],ibarra.N[intersect_genes,])

write_rds(list(merged.E=merged.E,merged.N=merged.N), "data/singlecell/velocyto.merged.loom.Rds", compress = "gz")

merged.seurat <- read_rds("data/singlecell/merged.Rds")

merged.emb <- Embeddings(merged.seurat,"bbknn")

table(colnames(merged.E) %in% rownames(merged.emb))
table(rownames(merged.emb) %in% colnames(merged.E))
merged.emb.int <- merged.emb[rownames(merged.emb) %in% colnames(merged.E), ]
plot(merged.emb.int)

merged.emb.tec <- merged.emb[!merged.seurat$simplified %in% c("stroma", "etc"),]
intersect_cells <- intersect(rownames(merged.emb.tec), colnames(merged.E))

merged.emb.tec.int <- merged.emb.tec[intersect_cells, ]
merged.emb.tec.exclude <- merged.emb.tec[!rownames(merged.emb.tec) %in% intersect_cells, ]
merged.E.tec <- merged.E[,intersect_cells]
merged.N.tec <- merged.N[,intersect_cells]

plot(merged.emb.tec.int)
points(merged.emb.tec.exclude,col="red")

##
# save.image("~/Projects/thymus_single_cell/final/save_2019_4_18_velocyto_merged.RData")
# 
# load("~/Projects/thymus_single_cell/final/save_2019_4_18_velocyto_merged.RData")

save(merged.E.tec,
     merged.N.tec,
     merged.emb.tec.int,
     merged.seurat,
     file = "data/velocyto.prepall.RData")

merged.velocityestimates <- gene.relative.velocity.estimates(merged.E.tec,merged.N.tec,deltaT=1,kCells=10,n.cores=20)

velocity.on.emb <- show.velocity.on.embedding.cor(merged.emb.tec.int,merged.velocityestimates,
                                                  n=40,scale='sqrt',cex=0.8,arrow.scale=1.7,
                                                  show.grid.flow=TRUE,min.grid.cell.mass=0.5,
                                                  grid.n=60,arrow.lwd=1,do.par=F,cell.border.alpha = 0.1,
                                                  return.details=T,n.cores=35)
# grid.sd= 0.2081906  min.arrow.size= 0.004163811  max.grid.arrow.length= 0.09709767
velocity.on.emb$garrows


save(merged.E.tec,
     merged.N.tec,
     merged.emb.tec.int,
     merged.seurat,
     merged.velocityestimates,
     velocity.on.emb,
     file = "data/velocyto.prepall.RData",
     compress = "gzip")



if("fogi due to long time and high memor require" == 1){
  #cleanup memory
  rm(list=ls()[!ls() %in% c("merged.E.tec","merged.N.tec","merged.emb.tec.int","merged.seurat",
                            "merged.velocityestimates","velocity.on.emb")])
  set.seed(42)
  merged.sampled.cellnames <- sample(colnames(merged.E.tec),400)
  merged.E.tec.sampled <- merged.E.tec[,merged.sampled.cellnames]
  merged.N.tec.sampled <- merged.N.tec[,merged.sampled.cellnames]
  merged.emb.tec.int.sampled <- merged.emb.tec.int[merged.sampled.cellnames,]
  merged.sampled.velocityestimates <- gene.relative.velocity.estimates(merged.E.tec.sampled,merged.N.tec.sampled,deltaT=1,kCells=5,n.cores=20)
  show.velocity.on.embedding.cor(merged.emb.tec.int.sampled,merged.sampled.velocityestimates,n=40,
                                 scale='sqrt',cex=0.8,arrow.scale=3,
                                 show.grid.flow=TRUE,min.grid.cell.mass=0.5,
                                 grid.n=40,arrow.lwd=1,do.par=F,cell.border.alpha = 0.1)
  # trajectory
  trajectory.sampled.default <- show.velocity.on.embedding.eu(
    merged.emb.tec.int.sampled,merged.sampled.velocityestimates,
    n=40,scale='sqrt',cex=0.8,nPcs=30,sigma=2.5,show.trajectories=TRUE,
    diffusion.steps=400,n.trajectory.clusters=10,ntop.trajectories=1,
    embedding.knn=T,control.for.neighborhood.density=TRUE,n.cores=15
  )
  # fogi...
}



velocity.on.emb

merged <- merged.seurat
cairo_pdf("figures/velocyto.pdf",height = 7/2.54,width=7/2.54,pointsize = 12*0.7)
{
  my_color_palette <- c("#F78981", "#CE425A", "#9D0721", "#0BE2A1", "#20A27B", 
                        "#00A1FF", "#0B7DC0", "#AB07FF", "#624B92", "#5100FF", 
                        "#002EFC", "#1F30BF", "#282C4D", "#1C2362", "#E38900", 
                        "#8E766B", "#715757", "#926650", "#BCBCC2", "#84848C", 
                        "#74E74C", "#6FA75A", "#102607", "#F766BF")
  par(mar = c(0.5,0.5,0.5,0.5),pty='s')
  plot(merged@reductions$bbknn@cell.embeddings, 
       pch = 20,
       cex=1,
       col = my_color_palette[as.numeric(Idents(merged))],
       # bty = 'l',
       xaxt='n',yaxt='n',xlab="",ylab="",
       xlim=c(1.9,11.9),ylim=c(-9.5,-0.5),
       asp = 1)
  garrows <- velocity.on.emb$garrows
  # max.grid.arrow.length <- sqrt(sum((par('pin')/c(length(gx),length(gy)))^2))*0.25
  max.grid.arrow.length <- 0.3
  alen <- pmin(max.grid.arrow.length,sqrt( ((garrows[,3]-garrows[,1]) * par('pin')[1] / diff(par('usr')[c(1,2)]) )^2 + ((garrows[,4]-garrows[,2])*par('pin')[2] / diff(par('usr')[c(3,4)]) )^2))
  lapply(1:nrow(garrows),function(i) arrows(garrows[i,1],garrows[i,2],garrows[i,3],garrows[i,4],length=0.3*alen[i],lwd=1,cex=0.05))
}
dev.off()
# velocity.on.emb%>% str


par(mar = c(1, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom",legend=levels(Idents(merged)), col = my_color_palette,pch=19, ncol = 2,bty="n",cex = 0.7)


dev.off()




