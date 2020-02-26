

"
signature plot is updated and can find in 'figure_purity_per_histologictype_vi....R'

"



dirpath="/home/users/sypark/00_Project/01_thymoma/10_Final_data/16_point_mutation_final_call/01_annotated_file"
meta_dt=read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/02_metadata/Thymoma_summary_191106_1stSheet.txt')

meta_dt$id[meta_dt$GTF2I_status2 == "m"]

filetable <- data.frame(id=str_extract(dir(path = dirpath, pattern = "snv"),"[^.]*"),
                        path=dir(path = dirpath, pattern = "snv",full.names = T),stringsAsFactors = F)
wgs_filetable <- c("/home/users/sypark/00_Project/01_thymoma/03_WGS/15_Timing/SNU_19_C_wgs.snv.timing_input",
                   "/home/users/sypark/00_Project/01_thymoma/03_WGS/15_Timing/SNU_26_C_wgs.snv.timing_input")

# TCGA-ZB-A966: hypermutator

mutant_vcf <- filetable$path[filetable$id %in% meta_dt$id[meta_dt$GTF2I_status2 == "m"]] %>%
  lapply(read_tsv, col_types=cols(.default = "c",POS = "i")) %>%
  lapply(select,1:6) %>% do.call(rbind,.)

wildtype_vcf <- filetable$path[filetable$id %in% meta_dt$id[meta_dt$GTF2I_status2 == "w"]] %>%
  lapply(read_tsv, col_types=cols(.default = "c",POS = "i")) %>%
  lapply(select,1:6) %>% do.call(rbind,.)

carcinoma_vcf <- filetable$path[filetable$id %in% meta_dt$id[meta_dt$GTF2I_status2 == "c"] & filetable$id != "TCGA-ZB-A966"] %>%
  lapply(read_tsv, col_types=cols(.default = "c",POS = "i")) %>%
  lapply(select,1:6) %>% do.call(rbind,.)

hypermutator_vcf <- filetable$path[filetable$id %in% "TCGA-ZB-A966"] %>%
  lapply(read_tsv, col_types=cols(.default = "c",POS = "i")) %>%
  lapply(select,1:6) %>% do.call(rbind,.)

SNU19_vcf <- read_tsv(wgs_filetable[1], col_types = cols(.default="c", POS="i")) %>%
  select(`#CHROM`,POS,ID,REF,ALT,QUAL)
SNU26_vcf <- read_tsv(wgs_filetable[2], col_types = cols(.default="c", POS="i")) %>%
  select(`#CHROM`,POS,ID,REF,ALT,QUAL)

source("~kjyi/src/mutationalsignature/mutationalsignature.R")

mutant_signature <- main(unique(mutant_vcf),signature_number = "1,5",samplename = "mutant")
wildtype_signature <- main(unique(wildtype_vcf),signature_number = "1,5",samplename = "wildtype")
carcinoma_signature <- main(unique(carcinoma_vcf),signature_number = "1,2,5,13",samplename = "carcinoma")
hypermutator_signature <- main(hypermutator_vcf,signature_number = "1,5,6",samplename = "hypermutator")
SNU19_signature <- main(SNU19_vcf, signature_number = "1,5",samplename = "SNU_19_C")
SNU26_signature <- main(SNU26_vcf, signature_number = "1,5",samplename = "SNU_26_C")

mutant_signature.d <- main((mutant_vcf),signature_number = "1,5",samplename = "mutant")
wildtype_signature.d <- main((wildtype_vcf),signature_number = "1,5",samplename = "wildtype")
carcinoma_signature.d <- main((carcinoma_vcf),signature_number = "1,2,5,13",samplename = "carcinoma")





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
nrow(mutant_vcf) # 1953
table(meta_dt$GTF2I_status2)
meta_dt$GTF2I_status2 == "m" # 71
nrow(wildtype_vcf) # 1205
nrow(carcinoma_vcf)
nrow(hypermutator_vcf) #1078

mutant_signature$exp$case="GTF2I-mutant"
wildtype_signature$exp$case="GTF2I-wildtype"
carcinoma_signature$exp$case="Thymic carcinoma"
hypermutator_signature$exp$case="TCGA-ZB-A966 (WES)"
SNU19_signature$exp$case="SNU_19_C (WGS)"
SNU26_signature$exp$case="SNU_26_C (WGS)"

decomp <- rbind(mutant_signature$exp,
                wildtype_signature$exp,
                carcinoma_signature$exp,
                hypermutator_signature$exp,
                SNU19_signature$exp,
                SNU26_signature$exp) %>%
  select(-Exposure) %>% spread(key = case,value = Proportion,fill = 0) %>%
  as.data.frame %>% column_to_rownames("Signature #") %>% 
  as.matrix() %>% {./100} %>% .[c(2,3,1,5,4),c(4,3,5,6,2,1)]


cairo_pdf("figures/signature1.pdf",height = 9/2.54,width=20/2.54,pointsize = 12*0.7)
layout(matrix(c(1,4,7,
                2,5,7,
                3,6,7),byrow=T,ncol=3),heights=c(1,1,1),widths=c(7,7,5))
par(mar=c(0,3,4,0),oma=c(5,0,0,0))
#1mut
sig_plot(mutant_signature.d)
mtext(expression(bold("GTF2I"^mut)~"(WES merged,n=71)"),at=0,adj = 0,line=0)
mtext("1,953 substitution",at=0,adj = 0,line=-1,font=1)
text(66.7,0.05,"*",cex=2)
#2wt
sig_plot(wildtype_signature.d)
mtext(expression(bold("GTF2I"^WT)~"(WES merged,n=53)"),at=0,adj = 0,line=0)
mtext("1,205 substitutions",at=0,adj = 0,line=-1,font=1)
#3.ca.merge
sig_plot(carcinoma_signature.d,axis=T)
mtext(expression(bold("Thymic carcinoma") ~ "(WES merged,n=13)"),at=0,adj = 0,line=0)
mtext("2,013 substitutions       (TCGA-ZB-A966 excluded)",at=0,adj = 0,line=-1)
#
sig_plot(hypermutator_signature)
mtext(expression(bold("TCGA-ZB-A966")~"(WES)"),at=0,adj = 0,line=0)
mtext("1,078 substitutions",at=0,adj = 0,line=-1)
#
sig_plot(SNU19_signature)
mtext(expression(bold("SNU_19_C")~"(WGS)"),at=0,adj = 0,line=0)
mtext("1,138 substitutions",at=0,adj = 0,line=-1)
#
sig_plot(SNU26_signature,axis = T)
mtext(expression(bold("SNU_26_C")~"(WGS)"),at=0,adj = 0,line=0)
mtext("364 substitutions",at=0,adj = 0,line=-1)
# dev.off()
par(mar=c(2,9,15,1),xpd=T)

# plot(1);plot(1);plot(1);plot(1);
sigpal <- c("#0BADBF","#ADD4D9","#F2CD5E","#F2AC57","#D94A3D") %>% rev
x <- barplot(decomp, border="white", xlab="",horiz=T,las=2,xlim=c(0,1),
             # xaxt="n",
             col=sigpal)

for(i in 1:ncol(decomp)){
  text(cumsum(decomp[,i]) - decomp[,i]/2,x[i],
       # labels = c(5,1,6,2,13),
       labels = c(1,2,5,6,13),
       cex=c(1,1,1,1,0.6),
       # col=ifelse(decomp[,i]==0,"#00000000",c("white","black","black","black","black")))
       col=ifelse(decomp[,i]==0,"#00000000","black"))
}
mtext("Signature decomposition",3,cex=0.8)
par(mar=c(0,0,4,0),new=T)
plot(1,col="#00000000",xlab='',ylab='',xaxt='n',yaxt='n',bty='n')
x2=legend("top",legend = c("Spontaneous 5-methylcytosine deamination",
                       "AID/APOBEC",
                       "Endogenous mutational process",
                       "Defective DNA mismatch repair",
                       "AID/APOBEC"),title = "COSMIC mutational signatures",
       pt.bg = sigpal,pch=22,pt.cex=2.5,xpd = T,horiz = F,col="#00000080")
text(x2$rect$left+0.035,y = x2$text$y,labels = c("1","2","5","6","13"),adj=0.5,
     cex=0.7)
dev.off()

