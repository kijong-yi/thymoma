# Survival analysis
#KMcurve
library(survminer)
library(survival)
meta_dt <- read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/02_metadata/Thymoma_summary_191115_1stSheet.txt')

meta_dt$group = ifelse(meta_dt$GTF2I_status2 %in% c("m","c"),c(m="mutant",c="carcinoma")[meta_dt$GTF2I_status2],
                       ifelse(meta_dt$corrected_IRS4_TPM >= 30,"wildtypeIRShigh","wildtypeIRSlow"))
meta_dt$stagegroup = ifelse(meta_dt$Stage %in% c('I','II'), 
                                         'I,II',
                                         ifelse(meta_dt$Stage %in% c("IVa","IVb"), "IVa/IVb","III"))

fit1 <- survfit(Surv(RecurFreeSurvival/365, Recur_event) ~ group, data = meta_dt,conf.type = "log-log")
names(fit1$strata) = c("Carcinoma","GTF2I-mutant","GTF2I-WT,IRS4 high","GTF2I-WT,IRS4 low")
fit1 %>% summary

fit2 <- survfit(Surv(RecurFreeSurvival/365, Recur_event) ~ group, data = meta_dt[meta_dt$group %in% c("wildtypeIRShigh","wildtypeIRSlow"),],conf.type = "log-log")
survdiff(Surv(RecurFreeSurvival/365, Recur_event) ~ group, 
         data = meta_dt[meta_dt$group %in% c("wildtypeIRShigh","wildtypeIRSlow"),])

survdiff(Surv(RecurFreeSurvival/365, Recur_event) ~ group, 
         data = meta_dt[meta_dt$group %in% c("wildtypeIRShigh","mutant"),])
"mutant"
"carcinoma"
"wildtypeIRShigh"
"wildtypeIRSlow"
survdiff(Surv(RecurFreeSurvival/365, Recur_event) ~ group, 
         data = meta_dt[meta_dt$group %in% c("wildtypeIRShigh","carcinoma"),])


as.data.frame(meta_dt[meta_dt$group %in% c("wildtypeIRShigh","wildtypeIRSlow"),]) %>% 
{coin::logrank_test(Surv(.$RecurFreeSurvival/365, .$Recur_event) ~ factor(.$group),
                    type = "Tarone-Ware", ties.method="average-scores")}


plot(fit2)
{
  lwd=2
  cairo_pdf("figures/KM_survival.1.pdf",height = 7/2.54,width=7/2.54,pointsize = 12*0.7)
  par(mar=c(4,4,0.7,0.7))
  plot(fit1,
       lwd=2,
       col=c("#00A087FF","#3C5488FF","#DC0000FF","#F39B7FFF"),
       ylab="Recurrence free survival probability",
       xlab="Time (years)",mark.time=T)
  legend("bottomright",legend=c(expression(GTF2I^mut),expression(GTF2I^WT~IRS4^low),expression(GTF2I^WT~IRS4^high),"Thymic carcinoma"),
         bty="n",lty=1,lwd=lwd,col=c("#3C5488FF","#F39B7FFF","#DC0000FF","#00A087FF"))
  dev.off()
  getOption("viewer")("figures/KM_survival.1.pdf")
}
if(F){
  median(meta_dt$final_cellularity[meta_dt$GTF2I_status2 == "w"])
  meta_dt$group3 = ifelse(meta_dt$GTF2I_status2 %in% c("m","c"),c(m="mutant",c="carcinoma")[meta_dt$GTF2I_status2],
           ifelse(meta_dt$final_cellularity>0.18,"wildtypePurityhigh","wildtypePuritylow"))
  fit3 <- survfit(Surv(RecurFreeSurvival/365, Recur_event) ~ group3, data = meta_dt,conf.type = "log-log")
  
  plot(fit3,
           lwd=2,
           col=c("#00A087FF","#3C5488FF","#DC0000FF","#F39B7FFF"),
           ylab="Recurrence free survival probability",
           xlab="Time (years)",mark.time=T)
}



#Coxph
thymoma_scores <- readRDS("~kjyi/Projects/thymus_single_cell/final/models/thymoma_scores.Rds")
thymoma_scores<-thymoma_scores %>% rownames_to_column("id")
m_meta_dt <- meta_dt %>% mutate(his_group = ifelse(histologic_type %in% c('A','AB','MN-T'),'A',
                                                   ifelse(histologic_type %in% c('TC','NE'),'C','B'))) %>%
  mutate(stage_group = ifelse(Stage %in% c('I','II'), 'I_II', 'III_IV'))
m_meta_dt$GTF2I_status2 <- factor(m_meta_dt$GTF2I_status2, levels = c('m','w','c'))

m_meta_dt <- left_join(m_meta_dt,thymoma_scores)

m_meta_dt <- filter(m_meta_dt,GTF2I_status2 != "c")

m_meta_dt$cTEC
m_meta_dt$Tuft
m_meta_dt$Progenitor %>% hist
m_meta_dt$his_group
# m_meta_dt$stage_group
m_meta_dt$"Masaoka stage" = ifelse(m_meta_dt$Stage %in% c('I','II'), 'I,II', 'III,IV')
m_meta_dt$"Masaoka stage group" = ifelse(m_meta_dt$Stage %in% c('I','II'), 'I,II', ifelse(m_meta_dt$Stage %in% c("IVa","IVb"), "IVa/IVb","III"))
m_meta_dt$GTF2I_status2
m_meta_dt$Gender = m_meta_dt$gender
mean(m_meta_dt$age_at_diagnosis)
m_meta_dt$"Age" = m_meta_dt$age_at_diagnosis
m_meta_dt$"Age at diagnosis" = ifelse(m_meta_dt$age_at_diagnosis > 60,"≥ 60","< 60") %>% factor(level=c("< 60","≥ 60"))
m_meta_dt$"Myasthenia gravis"= c("YES"="Present","NO"="Absent")[m_meta_dt$history_myasthenia_gravis]
m_meta_dt$ImmuneScore
m_meta_dt$"Progenitor Index" = m_meta_dt$Progenitor
m_meta_dt$"Progenitor Index group" = ifelse(m_meta_dt$Progenitor > 0, "high", "low") %>% factor(level=c("low","high"))

m_meta_dt$"GTF2I L404H" = ifelse(m_meta_dt$GTF2I_status2 %in% "m","Mutant","Wild-type")

m_meta_dt$"IRS4" = ifelse(m_meta_dt$corrected_IRS4_TPM>30,"TPM > 30", "TPM < 30")
m_meta_dt$"IRS4 overexp." = ifelse(m_meta_dt$corrected_IRS4_TPM>30,"TPM > 30", "TPM < 30")

m_meta_dt$"GTF2I/IRS4 genotype" = ifelse(m_meta_dt$GTF2I_status2 %in% "m","Mutant",
                                         ifelse(m_meta_dt$corrected_IRS4_TPM>30,"WT/IRS4-high",
                                                "WT/IRS4-low"))
m_meta_dt$"GTF2I/IRS4 genotype" = factor(m_meta_dt$"GTF2I/IRS4 genotype",levels = c("Mutant","WT/IRS4-low","WT/IRS4-high"))



bigcoxmodel <-coxph(Surv(RecurFreeSurvival, Recur_event) ~ 
                      Gender +
                      `Age at diagnosis` +
                      # Age +
                      `Myasthenia gravis` +
                      `Masaoka stage` +
                      # `Masaoka stage group` +
                      IRS4 + `GTF2I L404H` + `Progenitor Index`,
                    data=m_meta_dt); summary(bigcoxmodel); ggforest(bigcoxmodel,main = "Hazard ratio of Cox proportional model",fontsize =0.9,cpositions = c(0.02, 0.22, 0.4))

bigcoxmodel2 <-coxph(Surv(RecurFreeSurvival, Recur_event) ~ 
                       Gender +
                       `Age at diagnosis` +
                       # Age +
                       # `Myasthenia gravis` +
                       `Masaoka stage` +
                       # `Masaoka stage group` +
                       `IRS4 overexp.` + `GTF2I L404H` + `Progenitor Index`,
                     data=m_meta_dt); summary(bigcoxmodel2); ggforest(bigcoxmodel2,main = "Hazard ratio of Cox proportional model",fontsize =0.9,cpositions = c(0.02, 0.22, 0.4))

bigcoxmodel3 <-coxph(Surv(RecurFreeSurvival, Recur_event) ~ 
                       Gender +
                       # `Age at diagnosis` +
                       # Age +
                       # `Myasthenia gravis` +
                       # `Masaoka stage` +
                       `Masaoka stage group` +
                       # `IRS4 overexp.` + 
                       # `GTF2I L404H` + 
                       # his_group+
                       `GTF2I/IRS4 genotype`,
                     data=m_meta_dt); summary(bigcoxmodel3); ggforest(bigcoxmodel3,main = "Hazard ratio of Cox proportional model",fontsize =0.9,
                                                                      cpositions = c(0.02, 0.22, 0.4))


model=bigcoxmodel3
data=NULL
# function--------------------------------------------------------------------------------------------------------------
ggforest2 <- function (model, data = NULL, main = "",
                       cpositions = c(0.02, 0.18, 0.3,0.4,0.95),whisker_start_position=0.6,
                       fontsize = 0.7, 
                       refLabel = "reference", 
                       noDigits = 2, second_column_labels=NULL,first_column_labels=NULL) 
{
  
  conf.high <- conf.low <- estimate <- NULL
  stopifnot(class(model) == "coxph")
  data <- eval(model$call$data)
  terms <- attr(model$terms, "dataClasses")[-1]
  coef <- as.data.frame(broom::tidy(model))
  gmodel <- broom::glance(model)
  allTerms <- lapply(seq_along(terms), function(i) {
    var <- names(terms)[i]
    if (terms[i] %in% c("factor", "character")) {
      adf <- as.data.frame(table(data[, var]))
      cbind(var = var, adf, pos = 1:nrow(adf))
    }
    else if (terms[i] == "numeric") {
      data.frame(var = var, Var1 = "", Freq = nrow(data), 
                 pos = 1)
    }
    else {
      vars = grep(paste0("^", var, "*."), coef$term, value = TRUE)
      data.frame(var = vars, Var1 = "", Freq = nrow(data), 
                 pos = seq_along(vars))
    }
  })
  allTermsDF <- do.call(rbind, allTerms)
  colnames(allTermsDF) <- c("var", "level", "N", "pos")
  inds <- apply(allTermsDF[, 1:2], 1, paste0, collapse = "")
  rownames(coef) <- gsub(coef$term, pattern = "`", replacement = "")
  toShow <- cbind(allTermsDF, coef[inds, ])[, c("var", "level", 
                                                "N", "p.value", "estimate", "conf.low", "conf.high", 
                                                "pos")]
  toShowExp <- toShow[, 5:7]
  toShowExp[is.na(toShowExp)] <- 0
  toShowExp <- format(exp(toShowExp), digits = noDigits)
  toShowExpClean <- data.frame(toShow, pvalue = signif(toShow[, 
                                                              4], noDigits + 1), toShowExp)
  toShowExpClean$stars <- paste0(round(toShowExpClean$p.value, 
                                       noDigits + 1), " ", ifelse(toShowExpClean$p.value < 0.05, 
                                                                  "*", ""), ifelse(toShowExpClean$p.value < 0.01, "*", 
                                                                                   ""), ifelse(toShowExpClean$p.value < 0.001, "*", ""))
  toShowExpClean$ci <- paste0("(", toShowExpClean[, "conf.low.1"], 
                              "-", toShowExpClean[, "conf.high.1"], ")")
  toShowExpClean$estimate.1[is.na(toShowExpClean$estimate)] = refLabel
  toShowExpClean$stars[which(toShowExpClean$p.value < 0.001)] = "<0.001 ***"
  toShowExpClean$stars[is.na(toShowExpClean$estimate)] = ""
  toShowExpClean$ci[is.na(toShowExpClean$estimate)] = ""
  toShowExpClean$estimate[is.na(toShowExpClean$estimate)] = 0
  toShowExpClean$var = as.character(toShowExpClean$var)
  toShowExpClean$var[duplicated(toShowExpClean$var)] = ""
  # toShowExpClean$N <- paste0("(N=", toShowExpClean$N, ")")
  toShowExpClean <- toShowExpClean[nrow(toShowExpClean):1,]
  
  if(!is.null(second_column_labels)){toShowExpClean$level = second_column_labels}
  if(!is.null(first_column_labels)){toShowExpClean$var = first_column_labels}
  
  
  rangeb <- range(toShowExpClean$conf.low, toShowExpClean$conf.high, 
                  na.rm = TRUE)
  breaks <- axisTicks(rangeb/2, log = TRUE, nint = 7)
  rangeplot <- rangeb
  # whisker_start_position == 0.25 -> 1/3
  # whisker_start_position == 0.5 --> 1
  # whisker_start_position == 0.75 --> 3
  
  rangeplot[1] <- rangeplot[1] - diff(rangeb)*whisker_start_position/(1-whisker_start_position)
  rangeplot[2] <- rangeplot[2] + 0.17 * diff(rangeb)
  # rangeplot[2] <- rangeplot[2] + 0.1 * diff(rangeb)
  width <- diff(rangeplot)
  y_variable <- rangeplot[1] + cpositions[1] * width
  y_nlevel <- rangeplot[1] + cpositions[2] * width
  y_nn <- rangeplot[1] + cpositions[3] * width
  y_cistring <- rangeplot[1] + cpositions[4] * width
  y_stars <- rangeplot[1] + cpositions[5] * width
  x_annotate <- seq_len(nrow(toShowExpClean))
  annot_size_mm <- fontsize * as.numeric(grid::convertX(unit(theme_get()$text$size, "pt"), "mm"))
  
  p <- ggplot(toShowExpClean, aes(seq_along(var), exp(estimate))) +
    geom_rect(aes(xmin = seq_along(var) - 0.5, xmax = seq_along(var) +
                    0.5, ymin = exp(rangeplot[1]), ymax = exp(rangeplot[2]),
                  fill = ordered(seq_along(var)%%2 + 1))) +
    geom_rect(aes(xmin = length(var) + 1 - 0.5, xmax = length(var) + 1 + 0.5,
                  ymin = exp(rangeplot[1]), ymax = exp(rangeplot[2])),fill = "grey70") +
    scale_fill_manual(values = c("white", "grey90"), guide = "none") + 
    # geom_rect(xmin = 12 - 0.5, xmax = 13 +0.5, ymin = exp(rangeplot[1]), ymax = exp(rangeplot[2]),
    # fill = "red") +
    # geom_rect(xmin=0.5,xmax=1.5,ymin = exp(rangeplot[1]),ymax=exp(rangeplot[2]),fill="red")+
    geom_point(pch = 15, size = 4) +
    geom_errorbar(aes(ymin = exp(conf.low), ymax = exp(conf.high)), width = 0.15) +
    geom_segment(aes(x = 0.5, y = 1, xend = 11.5, yend = 1),linetype=3,size=0.1,col="black") +
    # geom_hline(yintercept = 1, linetype = 3) +
    coord_flip(ylim = exp(rangeplot)) + 
    # ggtitle(main) +
    scale_y_log10(name = "",
                  labels = sprintf("%g", breaks), expand = c(0.02, 0.02),
                  breaks = breaks) +
    theme_light() + 
    theme(panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(),
          legend.position = "none", panel.border = element_blank(),
          axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
    xlab("") +
    annotate(geom = "text", x = x_annotate, y = exp(y_variable),
             label = toShowExpClean$var, fontface = "bold", hjust = 0,
             size = annot_size_mm) +
    annotate(geom = "text", x = x_annotate, parse=T,
             y = exp(y_nlevel), hjust = 0, 
             # label = ifelse(toShowExpClean$level=="","(continuous)",as.character(toShowExpClean$level)),
             label = toShowExpClean$level,
             # vjust = -0.1,
             size = annot_size_mm) +
    annotate(geom = "text",
             x = x_annotate, y = exp(y_nn), label = toShowExpClean$N,
             fontface = "italic", hjust = 0,
             # vjust = ifelse(toShowExpClean$level == "", 0.5, 1.1),
             size = annot_size_mm) +
    annotate(geom = "text",
             x = x_annotate, y = exp(y_cistring), label = paste(toShowExpClean$estimate.1,toShowExpClean$ci),
             size = annot_size_mm, hjust = 0
             # vjust = ifelse(toShowExpClean$estimate.1 == "reference", 0.5, -0.1)
             ) +
    # annotate(geom = "text",
    #          x = x_annotate, y = exp(y_cistring), label = toShowExpClean$ci,
    #          size = annot_size_mm, vjust = 1.1, fontface = "italic") +
    annotate(geom = "text", x = x_annotate, y = exp(y_stars),
             label = toShowExpClean$stars, size = annot_size_mm,
             hjust = 0, fontface = "italic") +
    annotate(geom = "text",
             x = 0.5, y = exp(y_variable),
             label = paste0("# Events:^2 ", gmodel$nevent, "; Global p-value (Log-Rank): ", format.pval(gmodel$p.value.log, eps = ".001"),
                            " \nAIC: ", round(gmodel$AIC, 2), "; Concordance Index: ", round(gmodel$concordance, 2)), 
             size = annot_size_mm, hjust = 0, vjust = 1.2, fontface = "italic") +
    annotate(geom = "text",
             x = max(x_annotate)+1, y = c(exp(y_variable),exp(y_nlevel),exp(y_nn),exp(y_cistring),exp((y_stars+y_cistring)/2),exp(y_stars)),
             label = c("Factors","","n","Hazard ratio (95% CI)","","p"),
             size = annot_size_mm, hjust = 0, fontface = "bold");p
  gt <- ggplot_gtable(ggplot_build(p))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  ggpubr::as_ggplot(gt)
}
# ---------------------------------------------------------------------------------------------------------------------

p <- ggforest2(bigcoxmodel3,main = "",fontsize =1,cpositions = c(0.01, 0.21, 0.40,0.45,0.94),
               whisker_start_position=0.64,
          second_column_labels = c("Female",
                                   "Male",
                                   expression(paste('I, II')),
                                   "III",
                                   expression(paste('IVa, IVb')),
                                   expression(GTF2I^mut),
                                   expression(GTF2I^WT~IRS4^low),
                                   expression(GTF2I^WT~IRS4^high)) %>% rev,
          first_column_labels = c("","","GTF2I/IRS4 status","","","Masaoka stage","","Gender"))


cairo_pdf("figures/coxph.1.pdf",height = 950/254,width=1900/254,pointsize = 12*0.7)
print(p)
dev.off()
12*0.7

getOption("viewer")("figures/coxph.1.pdf")
# -------------------------------------------------------

fit1 <- survfit(Surv(RecurFreeSurvival/365, Recur_event) ~ group, data = meta_dt,conf.type = "log-log")
names(fit1$strata) = c("Carcinoma","GTF2I-mutant","GTF2I-WT,IRS4 high","GTF2I-WT,IRS4 low")

fit2 <- survfit(Surv(RecurFreeSurvival/365, Recur_event) ~ stagegroup, data = m_meta_dt,conf.type = "log-log")
names(fit2$strata) = c("I,II","III","IVa/IVb")

{
  cairo_pdf("figures/KM_survival.2.pdf",height = 7/2.54,width=14/2.54,pointsize = 12*0.7)
  par(mar=c(4,0,0.7,0),mfrow=c(1,2),oma=c(0,4,0,0.7))
  
  plot(fit2,
       lwd=2,
       col=c("grey60","grey40","black"),
       ylab="",
       xlab="",
       # yaxt='n',
       mark.time=T)
  mtext("Time (years)",1,2)
  mtext("Recurrence free survival probability",2,2)
  
  legend("bottomright",legend=c("I,II","III","IVa,IVb"),
         bty="n",lty=1,lwd=2,col=c("grey60","grey40","black"),title = "Masaoka stage")
  # mtext(text = "(Carcinoma excluded)",side = 1,at = 0,line = -2,adj = 0)
  plot(fit1,
       lwd=2,
       col=c("#00A087FF","#3C5488FF","#DC0000FF","#F39B7FFF"),
       ylab="",
       yaxt='n',
       xlab="",mark.time=T)
  legend("bottomright",legend=c(expression(GTF2I^mut),expression(GTF2I^WT~IRS4^low),expression(GTF2I^WT~IRS4^high),"Thymic carcinoma"),
         bty="n",lty=1,lwd=2,col=c("#3C5488FF","#F39B7FFF","#DC0000FF","#00A087FF"))
  mtext("Time (years)",1,2)
  # mtext("Recurrence free survival probability",2,2)
  
  dev.off()
}




