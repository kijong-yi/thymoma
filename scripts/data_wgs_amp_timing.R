library(tidyverse)

if(F){
  gmmp <- function(mutCN, k){
    library("mixtools")
    plot_mix_comps <- function(x, mu, sigma, lam) {
      lam * dnorm(x, mu, sigma)
    }
    # set.seed(1)
    wait = mutCN
    mixmdl <- normalmixEM(wait, k = k,mean.constr = c(1:k))
    mixmdl$x %>% density %>% plot()
    pal=rainbow(k)
    for(i in 1:k){
      plot_mix_comps(x=1:(100*(k+1))/100,mu=mixmdl$mu[i], sigma=mixmdl$sigma[i], lam = mixmdl$lambda[i]) %>% lines(1:(100*(k+1))/100,.,col=pal[i])
    }
    mixmdl
  }
  lambda = gmm(tmp$mutCN[1:20],k=A)
}

if(F){
  proba=snu19proba
  seg=snu19seg
  n.iter=10
  n.iter.inner=5
}
# function
main <- function(proba=snu19proba, seg=snu19seg, n.iter=100, n.iter.inner=50){
  require(mixtools)
  chr_arms = unique(seg$broadAMP)
  chr_arms = chr_arms[chr_arms!="."]
  gmm <- function(mutCN, k, iter = n.iter.inner){
    o = rep(NA,iter)
    for(it in 1:iter){
      try({
        invisible(capture.output({m <- normalmixEM(mutCN, k=k,mean.constr = c(1:k),ECM=T)}))
        # m <- suppressWarnings(normalmixEM(mutCN, k = k,mean.constr = c(1:k),ECM=T,verb=F))
        o[it] <- m$lambda[k]
      })
    }
    median(o,na.rm=T)
  }
  require(doMC)
  outmat = foreach(I = 1:length(chr_arms)) %do% {
    chr_arm = chr_arms[I]
    seg %>% filter(broadAMP==chr_arm) %>% {c(min(c(.$start_pos,.$end_pos)),max(c(.$start_pos,.$end_pos)))} -> rrr
    positionstart =  rrr[1]
    positionend = rrr[2]
    A = seg %>% filter(broadAMP==chr_arm) %>% filter(end_pos - start_pos == max(end_pos - start_pos)) %>% .$majCN
    B = seg %>% filter(broadAMP==chr_arm) %>% filter(end_pos - start_pos == max(end_pos - start_pos)) %>% .$minCN
    CHR= seg %>% filter(broadAMP==chr_arm) %>% filter(end_pos - start_pos == max(end_pos - start_pos)) %>% .$`#CHROM`
    tmp <- proba %>% filter(`#CHROM`==CHR & POS >= positionstart & POS <= positionend)
    num.mut = nrow(tmp)
    lambda = gmm(tmp$mutCN,k=A)
    est.early.mut.count = lambda * (ifelse(A==B,1,2)) * num.mut
    permutation = rep(NA,n.iter)
    for(i in 1:n.iter){
      bs <- sample(tmp$mutCN,size = num.mut,replace = T)
      # gmm(bs,k=A) * (ifelse(A==B,1,2)) * num.mut
      try({
        # hist(bs)
        # sum(gmmp(bs,k=A)$lambda[2:A]) * (ifelse(A==B,1,2)) * num.mut
        permutation[i] = gmm(bs,k=A) * (ifelse(A==B,1,2)) * num.mut
      })
    }
    est.early.mut.count.CIL=quantile(permutation, 0.025,na.rm=T)
    est.early.mut.count.CIU=quantile(permutation, 0.975,na.rm=T)
    list(CHR,positionstart,positionend,chr_arm,A,B,num.mut,lambda,est.early.mut.count,est.early.mut.count.CIL,est.early.mut.count.CIU)
  }
  outmat <- outmat %>% do.call(rbind,.) %>% as.data.frame
  colnames(outmat) = c("chr","start","end","chr_arm","A","B","num.mut","lambda","est.early.mut.count",
                       "est.early.mut.count.CIL","est.early.mut.count.CIU")
  outmat
}



# data
snu19proba <- read_tsv('~sypark/00_Project/01_thymoma/03_WGS/15_Timing/SNU_19_C.079_21.timing/SNU_19_C_wgs.snv.timing_input.cninfo.scF.sorted.kat.proba',col_types=cols(`#CHROM`="c"))
snu19seg <- read_tsv('~sypark/00_Project/01_thymoma/03_WGS/15_Timing/SNU_19_C.079_21.timing/SNU_19_C_wgs_segments.txt.edit.broadAMP',col_types=cols(`#CHROM`="c"))
snu26proba <- read_tsv('~sypark/00_Project/01_thymoma/03_WGS/15_Timing/SNU_26_C.026_31.timing/SNU_26_C_wgs.snv.timing_input.cninfo.scF.sorted.kat.proba',col_types=cols(`#CHROM`="c"))
snu26seg <- read_tsv('~sypark/00_Project/01_thymoma/03_WGS/15_Timing/SNU_26_C.026_31.timing/SNU_26_C_wgs.026_31.XY.sequenza3.withBP2_segments.txt.edit.broadAMP',col_types=cols(`#CHROM`="c"))
tail(snu26seg)
tail(snu26proba,20)
# run
snu19.early.mut.info <- main(proba = snu19proba,seg = snu19seg,n.iter = 1000, n.iter.inner = 30)
if(F){
  write_rds(snu19.early.mut.info,"data/snu19.early.mut.info.Rds")
}
snu26.early.mut.info <- main(proba = filter(snu26proba,`#CHROM`!="22"), # only 2 snv in chr22 yield error in gmm
                             seg = filter(snu26seg,`#CHROM`!="22"),
                             n.iter = 1000, n.iter.inner = 30)
snu26.early.mut.infoX <- main(proba = filter(snu26proba,`#CHROM`=="X") %>% mutate(mutCN=mutCN/2), # X,Y chromsome matter
                             seg = filter(snu26seg,`#CHROM`=="X"),
                             n.iter = 1000, n.iter.inner = 30)
snu26.early.mut.info <- bind_rows(snu26.early.mut.info[1:(nrow(snu26.early.mut.info)-1),],snu26.early.mut.infoX)
if(F){
  write_rds(snu26.early.mut.info,"data/snu26.early.mut.info.Rds")
}


snu19.early.mut.info<- read_rds("data/snu19.early.mut.info.Rds")
2.654002e-11:0.2708998

# additional figure - decomposing
# snu19.early.mut.info <- main(proba = snu19proba,seg = snu19seg,n.iter = 1000, n.iter.inner = 30)
gmmp <- function(mutCN, k,ymax=1.5,ylab="density"){
  library("mixtools")
  plot_mix_comps <- function(x, mu, sigma, lam) {
    lam * dnorm(x, mu, sigma)
  }
  set.seed(1)
  wait = mutCN
  mixmdl <- normalmixEM(wait, k = k,mean.constr = c(1:k))
  
  mixmdl$x %>% density %>% plot(main="",lwd=1.5,xlim=c(0,3.6),ylim=c(0,ymax),xaxt='n',las=2,ylab="",bty='L')
  mtext(ylab,side=2,cex=0.8,line=2.5)
  # pal=rainbow(k)
  pal=c("blue","red","red")
  for(i in 1:k){
    plot_mix_comps(
      x=1:(100*(k+1))/100,
      mu=mixmdl$mu[i], 
      sigma=mixmdl$sigma[i], 
      lam = mixmdl$lambda[i]) %>%
      lines(1:(100*(k+1))/100,.,col=pal[i],lty=2)
    
  }
  for(i in 1:3) abline(v=i,col="grey",lty=2)
  mixmdl
}
# lambda = gmm(tmp$mutCN[1:20],k=A)
proba=snu19proba
seg=snu19seg
chr_arms = unique(seg$broadAMP)
chr_arms = chr_arms[chr_arms!="."]





cairo_pdf("figures/gmm.pdf",width=48.1/25.4,height = 50/25.4,pointsize = 12*0.7*0.7)


par(mfrow=c(3,1),mar=c(0,4.5,0,1),oma=c(4.5,0,4.5,0))
chr_arms = c("7pAMP","1qAMP","8qAMP")
for (I in 1:length(chr_arms)) {
  chr_arm = chr_arms[I]
  seg %>% filter(broadAMP==chr_arm) %>% {c(min(c(.$start_pos,.$end_pos)),max(c(.$start_pos,.$end_pos)))} -> rrr
  positionstart =  rrr[1]
  positionend = rrr[2]
  A = seg %>% filter(broadAMP==chr_arm) %>% filter(end_pos - start_pos == max(end_pos - start_pos)) %>% .$majCN
  B = seg %>% filter(broadAMP==chr_arm) %>% filter(end_pos - start_pos == max(end_pos - start_pos)) %>% .$minCN
  CHR= seg %>% filter(broadAMP==chr_arm) %>% filter(end_pos - start_pos == max(end_pos - start_pos)) %>% .$`#CHROM`
  tmp <- proba %>% filter(`#CHROM`==CHR & POS >= positionstart & POS <= positionend)
  num.mut = nrow(tmp)
  lambda = gmmp(tmp$mutCN,k=A,ymax=c(1.6,1.6,0.9)[I],ylab=chr_arms[I])
  if(I==1){legend("topright",bty='n',legend = c("Actual distribution of mutations","Estimated # of post-amp. mutations","Estimated # of co-amp. mutations"),
                  lty=c(1,2,2),col=c("black","blue","red"),lwd=c(1.5,1,1))
  }
  abline(v=A,col="#AA000020",lwd="10")
  abline(v=B,col="#0000AA20",lwd="10")
  # legend("bottomright",legend = c(lambda[1],sum(lambda[-1])))
}
axis(1,at=c(0,1,2,3),labels = c(0,1,2,3))
mtext("Absolute copy number",1,line=2.5,cex=0.8)
mtext("Density distribution of \nestimated point mutation copy number",3,outer=T)
dev.off()

