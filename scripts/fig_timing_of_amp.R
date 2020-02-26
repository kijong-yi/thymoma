
meta_dt <- read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/02_metadata/Thymoma_summary_191205_1stSheet.txt')
# meta_dt$mut_rate <- meta_dt$n_pointmt/75  # panel size of SureSelectXT HumanAllExon V5+UTRs, 75MB (75000 kb), alignable genome 2800mbps
meta_dt$mut_rate <- meta_dt$n_pointmt/meta_dt$bait_size


histo_pal = c("#631879FF", "#3B4992FF", "#5F559BFF", "#A20056FF", "#BB0021FF", "#EE0000FF", "#008B45FF", "#008280FF")
names(histo_pal) = c("A","AB","MN-T","B1","B2","B3","NE","TC")
stage_pal=c('#ffffd4','#fed98e','#fe9929','#d95f0e','#993404', 'white')
stage_pal=c("#4352F7", "#DE75E6", "#F52075", "#780000", "#800000", "#4CE8A5")
names(stage_pal) <- c('I','II','III','IVa','IVb', '0')



tmp_df = meta_dt[meta_dt$GTF2I_status2 %in% c("w","m")&meta_dt$final_cellularity>=0.85,]

tmp_df2 = meta_dt[meta_dt$GTF2I_status2 %in% c("w","m")&meta_dt$final_cellularity>=0.5,]
tmp_df3 = meta_dt[meta_dt$GTF2I_status2 %in% c("w","m")&meta_dt$final_cellularity>=0.75,]
tmp_df4 = meta_dt[meta_dt$GTF2I_status2 %in% c("w","m")&meta_dt$final_cellularity>=0.8,]
tmp_df %>%
  with(plot(age_at_diagnosis,mut_rate,ylim=c(0,max(mut_rate)),xlim=c(0,90),
            col=histo_pal[histologic_type]))

model1 <- lm(mut_rate ~ age_at_diagnosis, data=tmp_df)
anova(model1)
summary(model1)

model2 <- lm(mut_rate ~ age_at_diagnosis, data=tmp_df2)
anova(model2)
summary(model2)

model3 <- lm(mut_rate ~ age_at_diagnosis + final_cellularity, data=tmp_df)
summary(model3)


model4 <- lm(mut_rate ~ age_at_diagnosis, data=tmp_df4)
summary(model4)


meta_dt$diploid_proportion %>% summary
tmp_df5 <- meta_dt %>% filter(GTF2I_status2 %in% c('w','m') & final_cellularity >= 0.5 & diploid_proportion >= 0.85)
model5 <- lm(mut_rate ~ age_at_diagnosis, data = tmp_df5)
summary(model5)
# slope=0.005787, p-value=0.03531
# circle size = tumor cell fraction
# circle fill = diploid proportion
# selected 26 cases = border thicker
n=26


# plot --------------------------------------------------------------------
cairo_pdf("figures/mutation_rate_age.pdf",height = 9/2.54,width=10/2.54,pointsize = 12*0.7)
par(mar=c(4.1,4.1,2.1,1.1))
ymax=max(meta_dt$mut_rate[meta_dt$GTF2I_status2 %in% c("w","m")])
plot(-100,-100,ylim=c(0,ymax),xlim=c(0,90),ylab = "Number of point mutation/Gbps",xlab="",yaxt='n')
mtext("Age at diagnosis",side=1,line=2)
axis(2,at=0:14/10,labels=0:14*100,las=2)
# model1 <- lm(mut_rate ~ age_at_diagnosis, data=tmp_df)
# model2 <- lm(mut_rate ~ age_at_diagnosis, data=tmp_df2)
# anova(model1)
# anova(model2)
model1=model5
CI <- predict(model1, newdata=tmp_df, interval="prediction",level=0.9) %>% as.data.frame
# CI2 <- predict(model2, newdata=tmp_df2, interval="prediction",level=0.9) %>% as.data.frame

abline(model1$coefficients[1],model1$coefficients[2],lwd=1.2,col="grey50",lty=2) # pink
polygon(c(sort(tmp_df$age_at_diagnosis),rev(sort(tmp_df$age_at_diagnosis))),
        c(CI$lwr[order(tmp_df$age_at_diagnosis)],rev(CI$upr[order(tmp_df$age_at_diagnosis)])),
        col="#F0000020",lty = 0)

meta_dt[meta_dt$GTF2I_status2 %in% c("w","m"),] %>%
  with(points(age_at_diagnosis,mut_rate,
            pch=21, 
            cex=ifelse(final_cellularity>=0.85,2.5,ifelse(final_cellularity>=0.5,2,1)),
            bg=circlize::colorRamp2(c(0.1,0.84,0.85,1),c("white","grey","grey","red"))(diploid_proportion),
            lwd=ifelse(final_cellularity>=0.5&diploid_proportion>=0.85,1.5,0.5)))
            # bg=ifelse(diploid_proportion>=0.85, "red","grey")))


ymax=max(meta_dt$mut_rate[meta_dt$GTF2I_status2 %in% c("w","m")])
mtext(side = 3,line = -1.4,at = 0,adj = 0,text = "Purity ≥ 0.5 & Diploid % ≥ 0.85 (n=26)",font=2)
mtext(side = 3,line = -2.8,at = 0,adj = 0,text = bquote(R^2~"="~.(round(summary(model1)$r.squared,3))~","~italic(p)~"="~.(round(anova(model1)$'Pr(>F)'[1],3))))
mtext(side = 3,line = -3.7,at = 0,adj = 0,text = paste0("slope = ",round(1000*model1$coefficients["age_at_diagnosis"],2),"/Gbps/year"))

library(ComplexHeatmap)
lgd = Legend(col_fun = circlize::colorRamp2(c(0.1,0.84,0.85,1),c("white","grey","grey","red")),
             title = "",grid_height = unit(0.25,"cm"),labels_gp=gpar(fontsize=7))
ComplexHeatmap::draw(lgd, x = unit(0.17, "npc"), y = unit(0.3, "npc"), just = c("left", "bottom"))
legend(0,0.78,xjust = 0,yjust=1,bty="n",
       legend=c("≥ 0.85","≥ 0.5","< 0.5"),pch=21,pt.cex=c(3,2,1),y.intersp = c(0.7,1,1))
text(0,0.78,"Purity",adj=c(0,0),cex=1,font=2)
dev.off()


# timing -----------------------------------------------------------------------------------------
tm_dt0 <- read_tsv('~sypark/00_Project/01_thymoma/03_WGS/15_Timing/Early_snv_count_191025.txt')
tm_dt1 <- read_tsv('~sypark/00_Project/01_thymoma/03_WGS/15_Timing/SNU_26_C.026_31.timing/SNU_26_C.026_31.broad_earlysnv')
tm_dt2 <- read_tsv('~sypark/00_Project/01_thymoma/03_WGS/15_Timing/SNU_19_C.079_21.timing/SNU_19_C.079_21.broad_earlysnv')
tm_dt1$early_SNV/tm_dt1$dist*1E6
tm_dt2$case = "SNU_19_C"
tm_dt1$case = "SNU_26_C"

snu19.early.mut.info <- read_rds("data/snu19.early.mut.info.Rds")
snu26.early.mut.info <- read_rds("data/snu26.early.mut.info.Rds")
snu19.early.mut.info$case = "SNU_19_C"
snu26.early.mut.info$case = "SNU_26_C"
tm_dt_new <- rbind(snu19.early.mut.info,snu26.early.mut.info)
for(col in colnames(tm_dt_new)){tm_dt_new[[col]] <- unlist(tm_dt_new[[col]])}
tm_dt_new <- tm_dt_new %>% 
  mutate(mr = est.early.mut.count/(end-start)*1E6,
         mr.upper = est.early.mut.count.CIU/(end-start)*1E6,
         mr.lower = est.early.mut.count.CIL/(end-start)*1E6)


tm_dt <- bind_rows(tm_dt1,tm_dt2)
tm_dt$id = tm_dt$case
# tm_dt$`early_SNVrate/Mb`
# tm_dt$early_SNVrate_range_low
# tm_dt$early_SNVrate_range_high

age_at_diagnosis2 = c("SNU_19_C"=45,
                     "SNU_26_C"=19)
genome_size_in_mb=2.7*10^9/10^6
# total_mut_rate = c("SNU_19_C_wgs"=1138,"SNU_26_C_wgs"=364) / genome_size_in_mb
total_mut_count = c("SNU_19_C"=1138,"SNU_26_C"=364)

# mutation_rate = 0.00615193 # /Mbps/year 6.15193 /Gbps/year
mutation_rate = model1$coefficients[2]
# mutation_rate = 0.00615

#A function to add arrows on the chart
cairo_pdf("figures/timing_amp.pdf",height = 12/2.54,width=16/2.54,pointsize = 12*0.7)
par(mar=c(10,10,4,4))
tm_dt <- arrange(tm_dt,desc(case),desc(`early_SNVrate/Mb`))
tm_dt_new <- arrange(tm_dt_new,desc(case),desc(mr))
tm_dt$total_clonal_SNV_rate
tm_dt$`early_SNVrate/Mb`
total_mut_count/genome_size_in_mb
tm_barplot <- barplot(tm_dt_new$mr,
                      names.arg = tm_dt_new$chr_arm,
                      horiz=T,
                      las=2,
                      xlim=c(-0.005,mutation_rate*(100)),
                      xaxt='n')


axis(side = 3,at=0:5/10,labels = 0:5/10*1000)
axis(side = 1,at=0:5/10,labels = 0:5/10*1000)
# axis(side=1)
axis(side = 1,line=3,at = mutation_rate*(0:8*10),labels = (0:8*10))

mtext(text = paste0("Estimated age of amplification\nbased on mutation rate ",round(mutation_rate*1000,2),"/Gbps/year"),line = 6,side = 1)
mtext("Number of amplified (early) base substitutions / Gbps",line = 2.5)
axis(side = 2,at = tm_barplot,labels = rep("",length(tm_barplot)))
arrows(x0 = pmax(0,tm_dt_new$mr.lower),
       x1 = pmax(0,tm_dt_new$mr.upper),
       y0 = tm_barplot,
       angle=90, code=3, length=0.025)
# tm_dt$id = c(rep("SNU_26_C",18),rep("SNU_19_C",3))
for(i in unique(tm_dt$case)){
  axis(side = 2,at = tm_barplot[c(range(which(tm_dt$case==i)))],labels = c("",""),line = 5.1,tck=0.02)
  mtext(text = str_replace(i,"_wgs",""),side = 2,outer = F,line = 5.6,
        at = median(tm_barplot[c(range(which(tm_dt$case==i)))]),
        col="black")
  segments(x0 = age_at_diagnosis2[i]*mutation_rate,
           y0 = tm_barplot[min(which(tm_dt$case==i))]-diff(tm_barplot)[1]/2,
           y1 = tm_barplot[max(which(tm_dt$case==i))]+diff(tm_barplot)[1]/2,
           lty=2)
  segments(x0 = total_mut_count[i]/genome_size_in_mb,
           y0 = tm_barplot[min(which(tm_dt$case==i))]-diff(tm_barplot)[1]/2,
           y1 = tm_barplot[max(which(tm_dt$case==i))]+diff(tm_barplot)[1]/2,
           lty=1,col="red")
}
# legend("topright",legend = "age at diagnosis",lty=2,bty="n")
legend("bottomright",legend = c("Age at diagnosis","total number of co-amplified\n subsitutions / Gbps"),lty=c(2,1),col=c("black","red"),bty="n")

dev.off()

