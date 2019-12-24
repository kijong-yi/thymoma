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
write_rds(snu19.early.mut.info,"data/snu19.early.mut.info.Rds")
snu26.early.mut.info <- main(proba = filter(snu26proba,`#CHROM`!="22"), # only 2 snv in chr22 yield error in gmm
                             seg = filter(snu26seg,`#CHROM`!="22"),
                             n.iter = 1000, n.iter.inner = 30)
snu26.early.mut.infoX <- main(proba = filter(snu26proba,`#CHROM`=="X") %>% mutate(mutCN=mutCN/2), # X,Y chromsome matter
                             seg = filter(snu26seg,`#CHROM`=="X"),
                             n.iter = 1000, n.iter.inner = 30)
snu26.early.mut.info <- bind_rows(snu26.early.mut.info[1:(nrow(snu26.early.mut.info)-1),],snu26.early.mut.infoX)
write_rds(snu26.early.mut.info,"data/snu26.early.mut.info.Rds")





