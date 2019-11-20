library(tidyverse)
library(stringr)
library(doMC)

# prep -------------------------------------------------------------------------------------------------------

plogsum <- read_rds("data/expbycoord/plogtpm_w2x10^7_o_0.2_sum.Rds")
plogavg <- read_rds("data/expbycoord/plogtpm_w2x10^7_o_0.2_average.Rds")
zlogsum <- read_rds("data/expbycoord/zlogtpm_w2x10^7_o_0.2_sum.Rds")
zlogavg <- read_rds("data/expbycoord/zlogtpm_w2x10^7_o_0.2_average.Rds")
# logsum <- read_rds("data/expbycoord/logtpm_w2x10^7_o_0.2_sum.Rds")
# logavg <- read_rds("data/expbycoord/logtpm_w2x10^7_o_0.2_average.Rds")
# tpmsum <- read_rds("data/expbycoord/tpm_w2x10^7_o_0.2_sum.Rds")
# tpmavg <- read_rds("data/expbycoord/tpm_w2x10^7_o_0.2_average.Rds")
# ngenes <- read_rds("data/expbycoord/logtpm_w2x10^7_o_0.2_ngenes.Rds")

metadata <-read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/02_metadata/Thymoma_summary_191115_1stSheet.txt') %>% 
  select(id,GTF2I_status2) %>%
  {.[match(colnames(plogsum[[1]]),.$id),]}
chrlist=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X")
# centromere and telomere

blacklist = circlize::read.cytoband()$df %>% filter(V5 %in% c("gvar","stalk","acen")) %>% .[,1:3] %>% {names(.) = c("ch","start","end");.} %>% mutate(ch=str_replace(ch,"chr",""))
blacklist=rbind(blacklist,data.frame(ch="3",start=200000001,end=252000001))
blacklist=rbind(blacklist,data.frame(ch="X",start=156000001,end=252000001))


# sequenza prep
# load sequenza segment tables
sqzseg <- foreach(i=colnames(plogsum[[1]]))%do%{
  file <- dir("/home/users/sypark/00_Project/01_thymoma/10_Final_data/18_sequenza/01_segments",
              paste0(i,".*segments.txt$"),full.names = T)
  if(length(file) > 1) stop("check file please")
  read_tsv(file, col_types = "ciidddddddddd") %>%
    filter(end.pos-start.pos > 10000000)
}
names(sqzseg) <- colnames(plogsum[[1]]); rm(i,file)


intervals = lapply(chrlist,
                   function(i){
                     data.frame(
                       chr=i,
                       start.pos=as.numeric(stringr::str_replace(rownames(plogsum[[i]]),"_.*","")),
                       end.pos=as.numeric(stringr::str_replace(rownames(plogsum[[i]]),".*_","")))
                   }) %>% do.call(rbind,.)

cnmat <- matrix(0,ncol=length(sqzseg),nrow=nrow(intervals))
colnames(cnmat) = names(sqzseg)
rownames(cnmat) = paste0("chr",intervals$chr,":",intervals$start.pos,"-",intervals$end.pos)
for(sample in names(sqzseg)){
  cnmat[,sample] = foreach(i = 1:nrow(intervals),.combine=c) %do%{
    sqzseg[[sample]] %>% filter(chromosome == intervals$chr[i] & start.pos <= intervals$end.pos[i] & end.pos >= intervals$start.pos[i]) %>%
      mutate(overlap=pmin(intervals$end.pos[i], end.pos) - pmax(start.pos,intervals$start.pos[i]),
             overlapprop = overlap/sum(overlap), wCNt = overlapprop*CNt) %>% .$wCNt %>% sum
  }
}

# rownames(cnmat)

# sqzseg[[7]] %>% filter(chromosome=="1") %>% {
#   plot(0,xlim=c(min(.$start.pos),max(.$end.pos)),ylim=c(0,max(.$CNt)))
#   segments(x0 = .$start.pos,x1 = .$end.pos,y0 = .$CNt, col=ifelse(.$end.pos-.$start.pos>5000000,"black","red"))
# }


sqzseg[[11]] %>% filter(chromosome=="3") %>% {
  plot(0,xlim=c(min(.$start.pos),max(.$end.pos)),ylim=c(0,max(.$CNt)))
  segments(x0 = .$start.pos,x1 = .$end.pos,y0 = .$CNt, col=ifelse(.$end.pos-.$start.pos>5000000,"black","red"))
}



rowProp  = function(x){rowSums(x)/ncol(x)}

colnames(cnmat) == metadata$id
ncasecoordampdel=data.frame(
  chr=rownames(cnmat) %>% str_replace(":.*","") %>% str_replace("chr",""),
  m_amp=rowProp(cnmat[,metadata$GTF2I_status2 == "m"] > 2.5),
  m_del=rowProp(cnmat[,metadata$GTF2I_status2 == "m"] < 1.5),
  w_amp=rowProp(cnmat[,metadata$GTF2I_status2 == "w"] > 2.5),
  w_del=rowProp(cnmat[,metadata$GTF2I_status2 == "w"] < 1.5),
  c_amp=rowProp(cnmat[,metadata$GTF2I_status2 == "c"] > 2.5),
  c_del=rowProp(cnmat[,metadata$GTF2I_status2 == "c"] < 1.5)
)

# cnmat[grepl("chrX",rownames(cnmat)),4:8]

gtf=read_tsv(paste0("/home/users/kjyi/Projects/thymus_single_cell/final/",
                    "expression/IRS4/merged_gtf/Homo_sapiens.GRCh38.95.mod3.gtf"),skip=5,
             col_names = c(
               "chr","spircename","feature","start","end",
               "score","strand","frame","attribute"),
             col_types="ccciicccc") %>%
  dplyr::filter(feature == "gene") %>%
  dplyr::select(chr,start,end,attribute) %>%
  mutate(attribute = str_replace(attribute, '.*gene_name \"','') %>%
           str_replace('\".*','')) %>%
  as.data.frame()
gene.chr = structure(gtf$chr,names=gtf$attribute)
gene.coord = structure((gtf$start+gtf$end)/2,names=gtf$attribute)



# functions --------------------------------------------------------------------------------------------------

plotorignal <- function(data,col=NULL){
  ymax=lapply(data[chrlist],max) %>% unlist %>% max
  ymin=lapply(data[chrlist],min) %>% unlist %>% min
  for (chr in chrlist){
    plot(0,0,xlim=c(0,max(as.numeric(stringr::str_replace(rownames(data[[chr]]),".*_","")))),
         ylim = c(ymin,ymax),col="white",xaxt="n",yaxt="n",xlab="",ylab="")
    mtext(chr)
    for(i in 1:ncol(data[[chr]])){
      if(is.null(col)){
        thiscol=ifelse(metadata$GTF2I_status2[i] == "c","#FF000050",
                       ifelse(metadata$GTF2I_status2[i] == "w", "#00FF0090","#0000FF50"))
      } else {
        thiscol = col[i]
      }
      y=data[[chr]][,i]
      x=(as.numeric(stringr::str_replace(names(data[[chr]][,i]),"_.*","")) + 
           as.numeric(stringr::str_replace(names(data[[chr]][,i]),".*_","")))/2
      mask=filter(blacklist, ch==chr)
      masked=unlist(lapply(x,function(X){ifelse(any(X>mask$start & X<mask$end), T,F)}))
      y[masked]=0
      lines(y = y,
            x = x,
            col=thiscol)}
    if(chr == "X") text(108719949,ymax*0.6,"*",cex=2,col="black")
  }
  legend("topright",legend=c("w","m","c"),col=c("#00FF00", "#0000FF", "#FF0000"),lty=1,bty='n')
}

plot_median_subtracted <- function(data,col=NULL){
  mediansubtracted <- data[chrlist]
  for(i in chrlist){
    mediansubtracted[[i]] = sweep(data[[i]], 1, apply(data[[i]],1,median),`-`)
  }
  ymax=lapply(mediansubtracted,max) %>% unlist %>% max
  ymin=lapply(mediansubtracted,min) %>% unlist %>% min
  for (chr in chrlist){
    plot(0,0,xlim=c(0,max(as.numeric(stringr::str_replace(rownames(mediansubtracted[[chr]]),".*_","")))),
         ylim = c(ymin,ymax),col="white",xaxt="n",yaxt="n",xlab="",ylab="",bty='n')
    box(col="grey")
    # mtext(chr)
    for(i in 1:ncol(mediansubtracted[[chr]])){
      if(is.null(col)){
        thiscol=ifelse(metadata$GTF2I_status2[i] == "c","#FF000050",
                       ifelse(metadata$GTF2I_status2[i] == "w", "#00FF0090","#0000FF50"))
      } else {
        thiscol = col[i]
      }
      
      y = mediansubtracted[[chr]][,i]
      x = (as.numeric(stringr::str_replace(names(mediansubtracted[[chr]][,i]),"_.*","")) + 
             as.numeric(stringr::str_replace(names(mediansubtracted[[chr]][,i]),".*_","")))/2
      mask=filter(blacklist, ch==chr)
      masked=unlist(lapply(x,function(X){ifelse(any(X>mask$start & X<mask$end), T,F)}))
      y[masked]=0
      
      lines(y=y,
            x=x,
            col=thiscol)
    }
    if(chr == "X") text(108719949,12000,"*",cex=2,col="black")
  }; legend("topright",legend=c("w","m","c"),col=c("#00FF00", "#0000FF", "#FF0000"),lty=1,bty='n')
}

# cohort difference
plotwithsd <-function(x,y,minuslogp,sdy,logscale=F, ...) {
  if(logscale){
    mylog <- function(x){
      sign(x)*log10(abs(x)+1)
    }
    plot(x,mylog(y),type = 'n',col="#00000000",...)
    # abline(h=0,col="#00000060",lwd=1.5) 
    segments(0,0,max(x),col="#00000060",lwd=1.5)
    polygon(x = c(x,rev(x)),y = c(mylog(y+sdy),mylog(rev(y)-rev(sdy))),col = "#00000020",lty = 0)
    lineplot_with_p(x,mylog(y),minuslogp)  
  }else{
    plot(x,y,type = 'n',...)
    # abline(h=0,col="#00000060",lwd=1.5)  
    segments(0,0,max(x),col="#00000060",lwd=1.5)
    polygon(x = c(x,rev(x)),y = c(y+sdy,rev(y)-rev(sdy)),col = "#00000020",lty = 0)
    lineplot_with_p(x,y,minuslogp)  
  }
}
# median - median

plot_cohort_diff <- function(data, col=NULL,genes=NULL, gene_mark_y=NULL,
                             logscale=F, ylim=NULL,ylab1=expression(bar(E)[w]-bar(E)[m]),ylab2="TPM/20Mb"){
  if(is.null(ylim)){
    ally<-lapply(chrlist, function(chr){
      apply(data[[chr]][,metadata$GTF2I_status2 == "w"],1,median) - apply(data[[chr]][,metadata$GTF2I_status2 == "m"],1,median)
    }) %>% unlist
    allysd <- lapply(chrlist, function(chr){
      sdy = sqrt(abs(apply(data[[chr]][,metadata$GTF2I_status2 == "w"],1,var) - apply(data[[chr]][,metadata$GTF2I_status2 == "m"],1,var)))
    }) %>% unlist
    ymax=max(ally+allysd)
    ymin=min(ally-allysd)
    ymax=max(c(abs(ymax),abs(ymin)))
    ymin=-max(c(abs(ymax),abs(ymin)))
    if(logscale){
      ymax=log10(ymax+1)
      ymin=-ymax
    }
  }else{
    ymin = ylim[1]
    ymax = ylim[2]
  }
  for (chr in chrlist){
    x=(as.numeric(stringr::str_replace(names(data[[chr]][,1]),"_.*","")) + 
         as.numeric(stringr::str_replace(names(data[[chr]][,1]),".*_","")))/2
    y=apply(data[[chr]][,metadata$GTF2I_status2 == "w"],1,median) - apply(data[[chr]][,metadata$GTF2I_status2 == "m"],1,median)
    p=apply(data[[chr]],1,function(x){t.test(x[metadata$GTF2I_status2 == "w"],
                                             x[metadata$GTF2I_status2 == "m"])$p.value})
    minuslogp=-log10(p)
    sdy = sqrt(abs(apply(data[[chr]][,metadata$GTF2I_status2 == "w"],1,var) - apply(data[[chr]][,metadata$GTF2I_status2 == "m"],1,var)))
    
    
    mask=filter(blacklist, ch==chr)
    masked=unlist(lapply(x,function(X){ifelse(any(X>mask$start & X<mask$end), T,F)}))
    y[masked]=0
    sdy[masked]=0
    minuslogp[masked]=0
    
    plotwithsd(x, y, minuslogp, sdy, xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(ymin,ymax),bty="n",logscale=logscale)
    # box(col="grey")
    # abline(v=max(x),col="grey",lty=2)
    segments(x0 = max(x),y0 = ymin,x1 = max(x),y1 = ymax,col="#00000060",lty=2)
    # abline(h=ymin,col="grey",lty=2)
    # if(chr == "X") text(108719949,2500,"*",cex=2,col="black")
    gidx=which(gene.chr == chr & names(gene.chr) %in% genes)
    this_chr_genes = genes[gene.chr[genes] == chr]
    this_chr_y = gene_mark_y[gene.chr[genes] == chr]
    # par(xpd=F)
    if(length(this_chr_genes)>0){
      if(length(this_chr_genes)==1){
        # text(gene.coord[this_chr_genes],this_chr_y,
             # "*",cex=2,col="black",adj=c(0.5,1))
        text(gene.coord[this_chr_genes],this_chr_y,
             this_chr_genes,adj=c(0.5,0.5),cex=1,col="black")
      }else{
        wordcloud::textplot(gene.coord[this_chr_genes],this_chr_y,this_chr_genes,new = F)
      }
    }
    if(chr == "1"){
      axis(2,col="grey")
      mtext(ylab1, side=2, line=4, adj=0.5)
      mtext(ylab2, side=2, line=2, adj=0.5)
      legend("bottomleft",c("p<0.001","±1 SD"), col=c("red","grey"),lty=c(1,1),lwd=c(1,7),bty='n')
      
    }
  }
}

if(F){
  parlayout(2)
  plot_sqz_cohort()
  plot_cohort_diff(zlogsum, ylim=c(-50,50),genes = concord_genes, gene_mark_y = marky,
                   ylab1="Expression gap",ylab2=expression(GTF2I^WT-GTF2I^mut))
}


lineplot_with_p <- function(x,y,p){ # deprecated
  P=(c(p[1],p) + c(p,p[length(p)]))/2  
  for(i in 1:(length(x)-1)){
    segments(x0 = x[i],y0 = y[i],x1 = x[i+1],y1 = y[i+1],col = circlize::colorRamp2(breaks = c(0,3),
                                                                                    colors=c('black','red'))(P[i]))
  }
}

lineplot_with_p <- function(x,y,p){
  P=(c(p[1],p) + c(p,p[length(p)]))/2  
  for(i in 1:(length(x)-1)){
    segments(x0 = x[i],y0 = y[i],x1 = x[i+1],y1 = y[i+1],col = ifelse(P[i]>3,"red","grey20"))
  }
}


mypolygon <- function(x,y,col=rgb(0.2,0.1,0.5,0.2)) {
  polygon( 
    c(min(x), x, max(x)), 
    c(0, y, 0), 
    col=col, border=F
  )
}


#plot sequenza cohort summary
plot_sqz_cohort <-function(){
  for (chr in chrlist){
    xmin=c(13000000)
    xmax=max(as.numeric(stringr::str_replace(rownames(plogsum[[chr]]),".*_","")))
    plot(0,0,xlim=c(xmin,xmax),
         ylim = c(-1,1),col="#00000000",xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
    # box(col="grey")
    # abline(v=xmax,col="grey",lty=2)
    # abline(h=-1,col="grey")
    segments(x0 = xmax,y0 = -1.2,x1 = xmax,y1 = 1,col="#00000060",lty=2)
    segments(x0 = xmin,y0 = -1.2,x1 = xmax,y1 = -1.2,col="#00000060")
    
    if(chr == "1"){
      axis(2,col="grey")
      mtext("% of cases", side=2, line=4, adj=0.5)
      mtext("gain", side=2, line=2, adj=0.75)
      mtext("loss", side=2, line=2, adj=0.25)
      # mtext(c(1,0,1), side=2, line=0, adj=c(0,0.5,1))
      legend("bottomleft",legend=c(expression(GTF2I^WT),expression(GTF2I^mut)),fill=c(rgb(0.1,0.7,0.1,0.3), rgb(0.1,0.1,0.7,0.3)),bty='n')
    }
    if(chr %in% c("1","2","3","4","5","6","7","8","9","X")){
      mtext(paste0("chr",chr))
    }else{
      mtext(chr)
    }
    x =(intervals$start.pos[intervals$chr==chr] + intervals$end.pos[intervals$chr==chr])/2
    m_amp = ncasecoordampdel$m_amp[ncasecoordampdel$chr == chr]
    w_amp = ncasecoordampdel$w_amp[ncasecoordampdel$chr == chr]
    m_del = ncasecoordampdel$m_del[ncasecoordampdel$chr == chr]
    w_del = ncasecoordampdel$w_del[ncasecoordampdel$chr == chr]
    # c_amp = ncasecoordampdel$c_amp[ncasecoordampdel$chr == chr]
    # c_del = ncasecoordampdel$c_del[ncasecoordampdel$chr == chr]
    
    mask=filter(blacklist, ch==chr)
    masked=unlist(lapply(x,function(X){ifelse(any(X>mask$start & X<mask$end), T,F)}))
    m_amp[masked]=0
    w_amp[masked]=0
    m_del[masked]=0
    w_del[masked]=0
    # c_amp[masked]=0
    # c_del[masked]=0
    
    mypolygon(x,m_amp,col=rgb(0.1,0.1,0.7,0.3))
    mypolygon(x,w_amp,col=rgb(0.1,0.7,0.1,0.3))
    # mypolygon(x,c_amp,col=rgb(0.7,0.1,0.1,0.2))
    # abline(h=0,col="black",lwd=1.5)
    segments(x0 = 0,y0 = 0,x1 = xmax+1000,y1 = 0,col="#000000C0")
    
    mypolygon(x,-m_del,col=rgb(0.1,0.1,0.7,0.3))
    mypolygon(x,-w_del,col=rgb(0.1,0.7,0.1,0.3))
    # mypolygon(x,-c_del,col=rgb(0.7,0.1,0.1,0.2))
    # }; legend("topright",legend=c("w","m","c"),fill=c(rgb(0.1,0.7,0.1,0.3), rgb(0.1,0.1,0.7,0.3), rgb(0.7,0.1,0.1,0.2)),bty='n')
  }
  
}

parlayout <- function(nrow=3){
  par(mfrow = c(nrow,23),mar=c(0.2,0,0.3,0),oma=c(0,6,1.8,0),xpd=NA);
  mat=matrix(1:(23*nrow),nrow,byrow=T)
  chrsize = foreach(chr=chrlist,.combine=c) %do%{
    max(as.numeric(stringr::str_replace(rownames(plogsum[[chr]]),".*_","")))}
  layout(mat, widths = chrsize,
         heights = rep.int(1,nrow(mat)), respect = FALSE)
}

parlayout2 <- function(){
  nrow=2
  par(mfrow = c(nrow,23),mar=c(0,0,0,0),oma=c(1,6,1.8,0),xpd=NA);
  mat=matrix(1:(23*nrow),nrow,byrow=T)
  chrsize = foreach(chr=chrlist,.combine=c) %do%{
    max(as.numeric(stringr::str_replace(rownames(plogsum[[chr]]),".*_","")))}
  layout(mat, widths = chrsize,
         heights = c(1,2), respect = FALSE)
}



# figures ------------------------------------------------------------------------------------

# plot_gather <- function(data,logscale=F){
#   parlayout(3)
#   plot_sqz_cohort()
#   plot_median_subtracted(data)
#   plot_cohort_diff(data,logscale=logscale)
# }

# plot_gather(plogsum)
# plot_gather(plogavg)
# plot_gather(zlogsum)
# plot_gather(zlogavg)
# plot_gather(logsum)
# plot_gather(logavg)
# plot_gather(tpmsum)
# plot_gather(tpmsum,logscale=T)
# 
# parlayout(3)
# plot_sqz_cohort()
# plot_cohort_diff(tpmsum,logscale=F,ylim=c(-3000,3000))
# plot_cohort_diff(tpmsum,logscale=T)

cairo_pdf("figures/coordexp.1.pdf",height = 6.2/2.54,width=21/2.54,pointsize = 12*0.7)
parlayout(2)
plot_sqz_cohort()
plot_cohort_diff(zlogsum, ylim=c(-50,50),genes = c("IRS4","TBX1"), gene_mark_y = c(40,40),
                 ylab1="Expression gap",ylab2=expression(GTF2I^WT-GTF2I^mut))
text(108719949,35,"*",cex=2,col="black")
dev.off()


tmp_zlogsum <- zlogsum
for(i in names(tmp_zlogsum)){
  rownames(tmp_zlogsum[[i]]) = paste(i,rownames(tmp_zlogsum[[i]]),sep = ":")
}
tmp_zlogsum <- tmp_zlogsum %>% do.call(rbind, .)


avg_diff = apply(tmp_zlogsum[,metadata$GTF2I_status2 == "w"],1,median) - apply(tmp_zlogsum[,metadata$GTF2I_status2 == "m"],1,median)
hist(avg_diff)
table(abs(avg_diff)<5)
WT_high_region <- names(avg_diff)[avg_diff > 10]
WT_low_region <- names(avg_diff)[avg_diff < -10]

WT_high_genes <- volc_dt %>% filter(WT_mean - MT_mean >= 1.5 & t.test.pvalue < 0.001/nrow(volc_dt)) %>% .$gene
WT_low_genes <- volc_dt %>% filter(WT_mean - MT_mean <= -1.5 & t.test.pvalue < 0.001/nrow(volc_dt)) %>% .$gene

gene.chr[WT_high_genes]
gene.coord[WT_high_genes]

library(GenomicRanges)
gr.wthighregion <- GRanges(str_replace(WT_high_region,":.*",""),
                           ranges=IRanges(as.numeric(str_extract(WT_high_region,"(?<=:)[0-9]*(?=_)")),
                                          end=as.numeric(str_extract(WT_high_region,"(?<=_)[0-9]*"))))
gr.wthighgenes <- GRanges(gene.chr[WT_high_genes],
                          ranges=IRanges(gene.coord[WT_high_genes],
                                         end=gene.coord[WT_high_genes]+1),
                          id = WT_high_genes)
gr.wtlowregion <- GRanges(str_replace(WT_low_region,":.*",""),
                           ranges=IRanges(as.numeric(str_extract(WT_low_region,"(?<=:)[0-9]*(?=_)")),
                                          end=as.numeric(str_extract(WT_low_region,"(?<=_)[0-9]*"))))
gr.wtlowgenes <- GRanges(gene.chr[WT_low_genes],
                          ranges=IRanges(gene.coord[WT_low_genes],
                                         end=gene.coord[WT_low_genes]+1),
                          id = WT_low_genes)


concord_wt_high_genes <- subsetByOverlaps(gr.wthighgenes,gr.wthighregion)$id
concord_wt_low_genes <- subsetByOverlaps(gr.wtlowgenes,gr.wtlowregion)$id
concord_genes = c(concord_wt_high_genes,concord_wt_low_genes)
# names(concord_genes) = c(rep("wt_high",length(concord_wt_high_genes)),
#                          rep("wt_low",length(concord_wt_low_genes)))
# 22   19770148-19770149      * |        TBX1
marky = (volc_dt$WT_mean - volc_dt$MT_mean)[match(concord_genes,volc_dt$gene)] * 15




cairo_pdf("figures/coordexp.2.pdf",height = 6.2/2.54,width=24/2.54,pointsize = 12*0.7)
parlayout(2)
plot_sqz_cohort()
plot_cohort_diff(zlogsum, ylim=c(-50,50),genes = concord_genes, gene_mark_y = marky,
                 ylab1="Expression gap",ylab2=expression(GTF2I^WT-GTF2I^mut))
# text(108719949,35,"*",cex=2,col="black")
dev.off()

save.image("data/figure_exp_cnv_by_coord.Rda")
