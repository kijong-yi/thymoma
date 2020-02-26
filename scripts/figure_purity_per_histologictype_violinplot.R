

"burden is low"
"purity distribution"
"signature"

================================================================================


# load metadata
meta_dt <- read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/02_metadata/Thymoma_summary_191204_1stSheet.txt')  


"Fig 1a. burden is low"

cairo_pdf("figures/mutation_burden_tmtc.pdf",height = 6.5/2.54,width=9.25/2.54,pointsize = 12*0.7)
p = boxplot(list(
  Thymomas = (meta_dt$n_pointmt/meta_dt$bait_size)[meta_dt$histologic_type %in% c("A","AB","B1","B2","B3","MN-T")],
  `Thymic carcinomas` = (meta_dt$n_pointmt/meta_dt$bait_size)[meta_dt$histologic_type %in% c("TC")]
),
ylim=c(0,3),ylab = "N. of point mutations / Mbps"
)
dev.off()

cairo_pdf("figures/mutation_detection_and_purity.pdf",height = 9/2.54,width=9.25/2.54,pointsize = 12*0.7)
library(plotrix)
histo_pal = c("#631879FF", "#3B4992FF", "#5F559BFF", "#A20056FF", "#BB0021FF", "#EE0000FF", "#008B45FF", "#008280FF")
names(histo_pal) = c("A","AB","MN-T","B1","B2","B3","NE","TC")
par(mfrow=c(1,1),mar=c(4,4,1,1))
gap.plot(
  meta_dt$final_cellularity,
  meta_dt$n_pointmt/meta_dt$bait_size,
  gap.axis = "y",
  gap = c(2.7,30.11),
  xlab="Tumor cell fraction",
  ylab="N. point mutation/Mbps",
  xlim=c(0,1),yaxt='n',
  bg=histo_pal[meta_dt$histologic_type],"bgcol"="white",
  pch=ifelse(meta_dt$GTF2I_status2=="c",23,21),
  cex=1.5)

axis.break(2, 2.7, breakcol="snow", style="gap")
axis.break(2, 2.7*(1+0.02), breakcol="black", style="slash")
axis.break(4, 2.7*(1+0.02), breakcol="black", style="slash")
axis(2, at=c(0,1,2,2.5),labels = c(0,1,2,2.5),las=2)
axis(2, at=c(2.7+30.35-30.11),labels = 30.35,las=2)


tmp_df <- data.frame(
  id=meta_dt$id[-which.max(meta_dt$n_pointmt)],
  x=meta_dt$final_cellularity[-which.max(meta_dt$n_pointmt)],
  y=(meta_dt$n_pointmt/meta_dt$bait_size)[-which.max(meta_dt$n_pointmt)],
  h=meta_dt$histologic_type[-which.max(meta_dt$n_pointmt)])

fit <- lm(y~x,data=tmp_df)
CI <- predict(fit, newdata=tmp_df, interval="prediction",level=0.9) %>% as.data.frame
lines(x=tmp_df$x, y=CI$fit, lwd=1.2,col="grey50")
polygon(c(sort(tmp_df$x),rev(sort(tmp_df$x))),c(CI$lwr[order(tmp_df$x)],rev(CI$upr[order(tmp_df$x)])),
        col="#00000010",lty = 0)

# points(tmp_df$x,tmp_df$y,
#        bg=histo_pal[meta_dt$histologic_type[match(tmp_df$id,meta_dt$id)]],
#        pch=21,
#        cex=1.5)

text(0.1,1.5,
     bquote(r~"="~.(round(summary(fit)$r.squared^0.5,3))
     ),adj=c(0,1))
text(0.1,1.5,
     bquote(
       italic(p)~"="~.(format(anova(fit)$'Pr(>F)'[1],scientific = T,digits=3))
     ),adj=c(0,2.5))
dev.off()




meta_dt$n_pointmt/meta_dt$bait_size


if(F){
  
  t.test((meta_dt$n_pointmt/meta_dt$bait_size)[meta_dt$histologic_type %in% c("A","AB","B1","B2","B3","MN-T")],(meta_dt$n_pointmt/meta_dt$bait_size)[meta_dt$histologic_type %in% c("TC") & meta_dt$id != "TCGA-ZB-A966"])
  
  
  # other possible plot
  p = plot(meta_dt$final_cellularity,meta_dt$n_pointmt/meta_dt$bait_size, ylim=c(0,3))
  p = boxplot((meta_dt$n_pointmt/meta_dt$bait_size)~meta_dt$histologic_type,ylim=c(0,3))
}


"Fig 1b. purity distribution"
library(ggplot2)
library(ggridges)

meta_dt2 <- meta_dt %>% filter(!histologic_type %in% c("MN-T", "NE"))
table(meta_dt2$histologic_type)
meta_dt2$histologic_type = factor(meta_dt2$histologic_type, levels=rev(c("A","AB","B1","B2","B3","TC","NE")))

table(is.na(meta_dt2$histologic_type))

cairo_pdf("figures/purity_per_histology.pdf",height = 7/2.54,width=9.575/2.54,pointsize = 12*0.7)

histo_pal = c("#631879FF", "#3B4992FF", "#5F559BFF", "#A20056FF", "#BB0021FF", "#EE0000FF", "#008B45FF", "#008280FF")
names(histo_pal) = c("A","AB","MN-T","B1","B2","B3","NE","TC")
ggplot(meta_dt2, aes(x = final_cellularity, y = histologic_type)) +
  geom_density_ridges(aes(fill = `histologic_type`)) +
  scale_fill_manual(values = histo_pal) + 
  theme_classic() + 
  scale_x_continuous("Tumor cell fraction",breaks = c(0:5/5), limits = c(-0.1,1.1)) + 
  scale_y_discrete("Histologic type") + 
  theme(axis.text.x = element_text(colour = "black")) +
  theme(axis.text.y = element_text(colour = "black")) + 
  theme(legend.position = "none")
dev.off()



"Fig 1c. signature is 1 and 5"
dirpath="/home/users/sypark/00_Project/01_thymoma/10_Final_data/16_point_mutation_final_call/01_annotated_file"
meta_dt <- read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/02_metadata/Thymoma_summary_191204_1stSheet.txt')  
# meta_dt$id[meta_dt$GTF2I_status2 == "m"]
filetable <- data.frame(id=str_extract(dir(path = dirpath, pattern = "snv"),"[^.]*"),
                        path=dir(path = dirpath, pattern = "snv",full.names = T),stringsAsFactors = F)
# wgs_filetable <- c("/home/users/sypark/00_Project/01_thymoma/03_WGS/15_Timing/SNU_19_C_wgs.snv.timing_input",
#                    "/home/users/sypark/00_Project/01_thymoma/03_WGS/15_Timing/SNU_26_C_wgs.snv.timing_input")


thymoma_vcf <- filetable$path[filetable$id %in% meta_dt$id[meta_dt$GTF2I_status2 %in% c("m","w")]] %>%
  lapply(read_tsv, col_types=cols(.default = "c",POS = "i")) %>%
  lapply(dplyr::select,1:6) %>% do.call(rbind,.)

carcinoma_vcf <- filetable$path[filetable$id %in% meta_dt$id[meta_dt$GTF2I_status2 == "c"] & filetable$id != "TCGA-ZB-A966"] %>%
  lapply(read_tsv, col_types=cols(.default = "c",POS = "i")) %>%
  lapply(dplyr::select,1:6) %>% do.call(rbind,.)

source("~kjyi/src/mutationalsignature/mutationalsignature.R")

thymoma_signature <- main(unique(thymoma_vcf),signature_number = "1,2,5,13",samplename = "carcinoma")
carcinoma_signature <- main(unique(carcinoma_vcf),signature_number = "1,2,5,13",samplename = "carcinoma")

sig_plot <- function(mutant_signature,axis=F){
  x=barplot(mutant_signature$mutDF$Original/sum(mutant_signature$mutDF$Original),
            col = rep(c("#56BEEC", "#050607", "#D53C32","#CBCACB","#ABCC72","#E7C9C6"),each = 16),
            las=2,ylab="Fraction of mutations",border = F,cex.names=0.7,xaxt="n")
  if(axis){
    # axis(side = 1,at=x,line=0.5,labels = paste0(substr(mutant_signature$mutDF$Trinucleotide,1,1),
    #                                             mutant_signature$mutDF$Substitution,
    #                                             substr(mutant_signature$mutDF$Trinucleotide,3,3)),
    #             las=2,cex.axis=0.2,tick = F,hadj = -0.5)
    lab = paste0(substr(mutant_signature$mutDF$Trinucleotide,1,1),
                 mutant_signature$mutDF$Substitution,
                 substr(mutant_signature$mutDF$Trinucleotide,3,3))
    mapply(FUN=graphics::axis, side = 1, at = x, labels = lab, cex.axis=0.3,las=2,tick=F,hadj=-0.5,line=0.5)
    # mapply(axis, labels = lab,side = 1, at = x)
    
    mtext(line = 1.5,side = 1,
          at = (x[((1:6)-0.5)*16] + x[1+((1:6)-0.5)*16])/2,
          text =  c("C>A","C>G","C>T","T>A","T>C","T>G"),
          col = c("#56BEEC", "#050607", "#D53C32","#CBCACB","#ABCC72","#E7C9C6"))
  }
}
sig_plot(carcinoma_signature,axis=T)
mutant_signature=carcinoma_signature

thymoma_signature$exp$case = "Thymoma"
carcinoma_signature$exp$case = "Thymic carcinoma"

table(meta_dt$GTF2I_status2 %in% c("m","w"))
dim(thymoma_vcf)
dim(carcinoma_vcf)

decomp <- rbind(thymoma_signature$exp,
                carcinoma_signature$exp) %>%
  dplyr::select(-Exposure) %>% spread(key = case,value = Proportion,fill = 0) %>%
  as.data.frame %>% column_to_rownames("Signature #") %>% 
  as.matrix() %>% {./100} %>% .[c(2,3,1,4),c(2,1)]


cairo_pdf("figures/signature1.pdf",width=925/254,height = 600/254,pointsize = 12*0.7)

# svglite::htmlSVG( 
#   width=925/254,height = 600/254,pointsize = 12*0.7,
#   {
      # layout(matrix(c(1,2,3),byrow=T,ncol=1),heights=c(1,1,1),widths=c(1))
par(mar=c(0,5,2,0),oma=c(4,0,0,0),mfrow=c(2,1))
sig_plot(thymoma_signature)
mtext("Thymoma (WES merged, n = 124)",at=0,adj = 0,line=0.5)
mtext("3,158 substitution",at=0,adj = 0,line=-0.5,font=1)
sig_plot(carcinoma_signature,axis=T)
mtext("Thymic carcinomas(WES merged, n=13)",at=0,adj = 0,line=0.5)
mtext("935 substitutions",at=0,adj = 0,line=-0.5)
  # })
dev.off()


cairo_pdf("figures/signature2.pdf",width=925/254,height = 500/254,pointsize = 12*0.7)

# svglite::htmlSVG(
#   width=925/254,height = 500/254,pointsize = 12*0.7,bg="grey",
#   {
    layout(matrix(c(1,2),byrow=T,ncol=1),heights=c(1,1),widths=c(1))
sigpal <- c("#0BADBF","#ADD4D9","#F2CD5E",          "#D94A3D") %>% rev
par(mar=c(1,8,1,1),oma=c(0,0,0,0))
x <- barplot(decomp, border="white", xlab="",horiz=T,las=2,xlim=c(0,1),col=sigpal,xaxt='n')
axis(1)
for(i in 1:ncol(decomp)){
  text(cumsum(decomp[,i]) - decomp[,i]/2,x[i],
       labels = c(1,2,5,13),
       cex=c(1,1,1,0.7),
       col=ifelse(decomp[,i]<0.05,"#00000000","black"))
}
# mtext("Signature decomposition",3,cex=1.2,line=1.5)
par(xpd=T)
par(mar=c(0.1,0,0,0))

# par(mar=c(0,0,0,0),)
plot(1,col="#00000000",xlab='',ylab='',xaxt='n',yaxt='n',bty='n')
x2=legend("bottom",legend = c("Spontaneous 5-methylcytosine deamination",
                           "AID/APOBEC",
                           "Endogenous mutational process",
                           # "Defective DNA mismatch repair",
                           "AID/APOBEC"),
          # title = "COSMIC mutational signatures",
          pt.bg = sigpal,pch=22,pt.cex=2.5,xpd = T,horiz = F,col="#00000080",cex=1)
text(x2$rect$left+0.026,y = x2$text$y,labels = c("1","2","5","13"),adj=0.5,
     cex=0.7) 
# })
# 


dev.off()


