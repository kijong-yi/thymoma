"------------------------------------------------------------------------------"
"                         Supplementary Fig.1                                  "
"                       APOBEC signature in TC                                 "
"                           Hypermutator signature                             "
"------------------------------------------------------------------------------"


dirpath="/home/users/sypark/00_Project/01_thymoma/10_Final_data/16_point_mutation_final_call/01_annotated_file"
meta_dt=read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/02_metadata/Thymoma_summary_191106_1stSheet.txt')
meta_dt$id[meta_dt$GTF2I_status2 == "c"]
cosmic="/home/users/jolim/Projects/S05_Jongsoo_Yoon/code/db/signatures_probabilities.txt" %>% read_tsv()
cosmic = arrange(cosmic,`Substitution Type`, Trinucleotide)
sig_plot <- function(mutant_signature,axis=F){
  x=barplot(mutant_signature$mutDF$Original/sum(mutant_signature$mutDF$Original),
            col = rep(c("#56BEEC", "#050607", "#D53C32","#CBCACB","#ABCC72","#E7C9C6"),each = 16),
            las=2,ylab="Fraction of mutations",border = F,cex.names=0.7,xaxt="n")
  if(axis){axis(side = 1,at=x,line=1,labels = paste0(substr(mutant_signature$mutDF$Trinucleotide,1,1),
                                                     mutant_signature$mutDF$Substitution,
                                                     substr(mutant_signature$mutDF$Trinucleotide,3,3)),
                las=2,cex.axis=0.2,tick = F,hadj = -0.5)
    mtext(line = 2.5,side = 1,
          at = (x[((1:6)-0.5)*16] + x[1+((1:6)-0.5)*16])/2,
          text =  c("C>A","C>G","C>T","T>A","T>C","T>G"),
          col = c("#56BEEC", "#050607", "#D53C32","#CBCACB","#ABCC72","#E7C9C6"))
  }
}
# cosmic$`Substitution Type` == carcinoma_each_signatures[[1]]$mutDF$Substitution
# cosmic$Trinucleotide == carcinoma_each_signatures[[1]]$mutDF$Trinucleotide
sig_plot2 <- function(mutant_signature,axis=F){
  x=barplot(mutant_signature/sum(mutant_signature),
            col = rep(c("#56BEEC", "#050607", "#D53C32","#CBCACB","#ABCC72","#E7C9C6"),each = 16),
            las=2,ylab="Fraction of mutations",border = F,cex.names=0.7,xaxt="n")
  if(axis){
    axis(side = 1,at=x,line=1,labels = paste0(substr(cosmic$Trinucleotide,1,1),
                                              cosmic$`Substitution Type`,
                                              substr(cosmic$Trinucleotide,3,3)),
                las=2,cex.axis=0.2,tick = F,hadj = -0.5)
    mtext(line = 2.5,side = 1,
          at = (x[((1:6)-0.5)*16] + x[1+((1:6)-0.5)*16])/2,
          text =  c("C>A","C>G","C>T","T>A","T>C","T>G"),
          col = c("#56BEEC", "#050607", "#D53C32","#CBCACB","#ABCC72","#E7C9C6"))
  }
}

filetable <- data.frame(id=str_extract(dir(path = dirpath, pattern = "snv"),"[^.]*"),
                        path=dir(path = dirpath, pattern = "snv",full.names = T),stringsAsFactors = F)
source("~kjyi/src/mutationalsignature/mutationalsignature.R")
filetable <- filetable[filetable$id %in% meta_dt$id[meta_dt$GTF2I_status2 == "c"],]
filetable$num_mut <- filetable[filetable$id %in% meta_dt$id[meta_dt$GTF2I_status2 == "c"],] %>%
  .$path %>%
  lapply(read_tsv, col_types=cols(.default = "c",POS = "i")) %>%
  lapply(dplyr::select,1:6) %>%
  lapply(nrow) %>% unlist

carcinoma_each_spectra <- filetable[filetable$num_mut >= 50 & filetable$num_mut < 1000,] %>%
  .$path %>%
  lapply(read_tsv, col_types=cols(.default = "c",POS = "i")) %>%
  lapply(dplyr::select,1:6)

names(carcinoma_each_spectra) = filetable$id[filetable$num_mut >= 50 & filetable$num_mut < 1000]
carcinoma_each_signatures <- lapply(carcinoma_each_spectra,function(x){main(unique(x),signature_number = "1,2,5,13",samplename = "mutant",plot=F)})
names(carcinoma_each_signatures) = filetable$id[filetable$num_mut >= 50 & filetable$num_mut < 1000]

carcinomas.decomp = data.frame(
  id=names(carcinoma_each_signatures),
  s1=lapply(carcinoma_each_signatures,function(x){x$exp[x$exp[,1]=="Signature 1","Proportion"]}) %>% lapply(function(x){ifelse(length(x)==0,0,x)}) %>% unlist,
  s2=lapply(carcinoma_each_signatures,function(x){x$exp[x$exp[,1]=="Signature 2","Proportion"]}) %>% lapply(function(x){ifelse(length(x)==0,0,x)}) %>% unlist,
  s5=lapply(carcinoma_each_signatures,function(x){x$exp[x$exp[,1]=="Signature 5","Proportion"]}) %>% lapply(function(x){ifelse(length(x)==0,0,x)}) %>% unlist,
  # s6=lapply(carcinoma_each_signatures,function(x){x$exp[x$exp[,1]=="Signature 6","Proportion"]}) %>% lapply(function(x){ifelse(length(x)==0,0,x)}) %>% unlist,
  s13=lapply(carcinoma_each_signatures,function(x){x$exp[x$exp[,1]=="Signature 13","Proportion"]}) %>% lapply(function(x){ifelse(length(x)==0,0,x)}) %>% unlist,
  cosSim=unlist(lapply(carcinoma_each_signatures,function(x){x$cosSim[[1]]})),
  stringsAsFactors=F)
left_join(carcinomas.decomp,
filetable[,-2])

"------------------------------------------------------------------------------"
"                         Decomposition plot                                   "
"------------------------------------------------------------------------------"

cairo_pdf("figures/sup.fig1a.tc.signature.pdf",height = 9/2.54,width=9/2.54,pointsize = 12*0.7)
# s <- svglite::svgstring(height = 9/2.54,width=9/2.54,pointsize = 12*0.7)
sigpal <- c("#0BADBF","#ADD4D9","#F2CD5E","#F2AC57","#D94A3D") %>% rev

par(mar=c(2,8,12,1),xpd=T)
decomp=t(carcinomas.decomp[,2:5])
colnames(decomp)=carcinomas.decomp$id
x <- barplot(decomp, border="white", xlab="",horiz=T,las=1,xlim=c(0,100),col=sigpal)
for(i in 1:ncol(decomp)){
  text(cumsum(decomp[,i]) - decomp[,i]/2,x[i],
       # labels = c(5,1,6,2,13),
       labels = c(1,2,5,13),
       cex=0.7,
       # col=ifelse(decomp[,i]==0,"#00000000",c("white","black","black","black","black")))
       col=ifelse(decomp[,i]==0,"#00000000","black"))
}

mtext("Signature decomposition",3,cex=0.8)
par(mar=c(0,0,4,0),new=T)
plot(1,col="#00000000",xlab='',ylab='',xaxt='n',yaxt='n',bty='n')
x2=legend("top",legend = c("Spontaneous 5-methylcytosine deamination",
                           "AID/APOBEC",
                           "Endogenous mutational process",
                           # "Defective DNA mismatch repair",
                           "AID/APOBEC"),title = "COSMIC mutational signatures",
          pt.bg = sigpal,pch=22,pt.cex=2.5,xpd = T,horiz = F,col="#00000080")
text(x2$rect$left+0.035,y = x2$text$y,labels = c("1","2","5","13"),adj=0.5,
     cex=0.7)

dev.off()
htmltools::browsable(htmltools::HTML(s()))

"------------------------------------------------------------------------------"
"                           96 spectra plot                                    "
"------------------------------------------------------------------------------"

cairo_pdf("figures/sup.fig1b.tc.signature2.pdf",height = 9/2.54,width=9/2.54,pointsize = 12*0.7)

par(mfrow=c(5,1),mar=c(0,3,2,0),oma=c(5,0,0,0))

"Signature 2"
sig_plot2(cosmic$`Signature 2`,axis = F)
mtext("COSMIC signature 2",at=0,adj = 0,line=0)

sig_plot2(cosmic$`Signature 13`,axis = F)
mtext("COSMIC signature 13",at=0,adj = 0,line=0)

sig_plot(carcinoma_each_signatures[["TCGA-ZC-AAA7"]], axis=F)
mtext("TCGA-ZC-AAA7",at=0,adj = 0,line=0)

sig_plot(carcinoma_each_signatures[["TCGA-ZB-A96V"]], axis=F)
mtext("TCGA-ZB-A96V",at=0,adj = 0,line=0)

sig_plot(carcinoma_each_signatures[["TCGA-3S-A8YW"]], axis=T)
mtext("TCGA-3S-A8YW",at=0,adj = 0,line=0)

text(10,0.04,"*",cex=2)
text(11,0.04,"*",cex=2,col="red")

dev.off()


"------------------------------------------------------------------------------"
"                           96 spectra plot, hypermutator                      "
"------------------------------------------------------------------------------"

cairo_pdf("figures/sup.fig1b.tc.signature3.pdf",height = 6/2.54,width=9/2.54,pointsize = 12*0.7)
par(mfrow=c(2,1), mar=c(0,3,2,0),oma=c(5,0,0,0))
hypermutator_vcf <- filetable$path[filetable$id %in% "TCGA-ZB-A966"] %>%
  lapply(read_tsv, col_types=cols(.default = "c",POS = "i")) %>%
  lapply(dplyr::select,1:6) %>% do.call(rbind,.)
hypermutator_signature <- main(hypermutator_vcf,signature_number = "1,5,6",samplename = "hypermutator")

sig_plot2(cosmic$`Signature 6`, axis=F)
mtext("COSMIC signature 6",at=0,adj = 0,line=0)
sig_plot(hypermutator_signature, axis=T)
mtext("TCGA-ZB-A966",at=0,adj = 0,line=0)
dev.off()

