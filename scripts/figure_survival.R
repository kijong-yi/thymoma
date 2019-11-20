# Survival analysis
#KMcurve
library(survminer)
library(survival)
meta_dt <- read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/02_metadata/Thymoma_summary_191115_1stSheet.txt')

meta_dt$group = ifelse(meta_dt$GTF2I_status2 %in% c("m","c"),c(m="mutant",c="carcinoma")[meta_dt$GTF2I_status2],
                       ifelse(meta_dt$corrected_IRS4_TPM >= 30,"wildtypeIRShigh","wildtypeIRSlow"))

fit1 <- survfit(Surv(RecurFreeSurvival/365, Recur_event) ~ group, data = meta_dt,conf.type = "log-log")
names(fit1$strata) = c("Carcinoma","GTF2I-mutant","GTF2I-WT,IRS4 high","GTF2I-WT,IRS4 low")
{
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
}



#Coxph
thymoma_scores <- readRDS("~/Projects/thymus_single_cell/final/models/thymoma_scores.Rds")
thymoma_scores<-thymoma_scores %>% rownames_to_column("id")
m_meta_dt <- meta_dt %>% mutate(his_group = ifelse(histologic_type %in% c('A','AB','MN-T'),'A',
                                                   ifelse(histologic_type %in% c('TC','NE'),'C','B'))) %>%
  mutate(stage_group = ifelse(Stage %in% c('I','II'), 'I_II', 'III_IV'))
m_meta_dt$GTF2I_status2 <- factor(m_meta_dt$GTF2I_status2, levels = c('m','w','c'))

m_meta_dt <- left_join(m_meta_dt,thymoma_scores)

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

bigcoxmodel <-coxph(Surv(RecurFreeSurvival, Recur_event) ~ 
                      Gender +
                      `Age at diagnosis` +
                      # Age +
                      `Myasthenia gravis` +
                      `Masaoka stage` +
                      # `Masaoka stage group` +
                      IRS4 + `GTF2I L404H` + `Progenitor Index`,
                    data=m_meta_dt); summary(bigcoxmodel); ggforest(bigcoxmodel,main = "Hazard ratio of Cox proportional model",fontsize =0.9,cpositions = c(0.02, 0.22, 0.4))

# GenderMALE                   0.4550    1.5762   0.6129  0.742   0.4578  
# Age at diagnosis≥60         -0.4706    0.6246   0.5829 -0.807   0.4195  
# `Myasthenia gravis`Present  -0.9377    0.3915   0.7631 -1.229   0.2191  
# `Masaoka stage`III,IV        1.4059    4.0790   0.5521  2.546   0.0109 *
# IRS4TPM > 30                 1.4299    4.1781   0.6776  2.110   0.0348 *
# `GTF2I L404H`Wild-type       1.4180    4.1289   0.8547  1.659   0.0971 .
# `Progenitor Index`           0.7911    2.2057   0.4382  1.805   0.0710 .


bigcoxmodel2 <-coxph(Surv(RecurFreeSurvival, Recur_event) ~ 
                      Gender +
                      `Age at diagnosis` +
                      # Age +
                      # `Myasthenia gravis` +
                      `Masaoka stage` +
                      # `Masaoka stage group` +
                      `IRS4 overexp.` + `GTF2I L404H` + `Progenitor Index`,
                    data=m_meta_dt); summary(bigcoxmodel2); ggforest(bigcoxmodel2,main = "Hazard ratio of Cox proportional model",fontsize =0.9,cpositions = c(0.02, 0.22, 0.4))

# GenderMALE                   0.5741    1.7755   0.6207  0.925   0.3551  
# Age at diagnosis≥60         -0.2836    0.7530   0.5635 -0.503   0.6147  
# `Masaoka stage`III,IV        1.3201    3.7439   0.5387  2.451   0.0143 *
# IRS4TPM > 30                 1.2808    3.5997   0.6695  1.913   0.0557 .
# `GTF2I L404H`Wild-type       1.4188    4.1324   0.8451  1.679   0.0932 .
# `Progenitor Index`           0.9473    2.5788   0.4091  2.315   0.0206 *



# function--------------------------------------------------------------------------------------------------------------
ggforest2 <- function (model, data = NULL, main = "", cpositions = c(0.02, 0.18, 0.3,0.4,0.95), fontsize = 0.7, refLabel = "reference", noDigits = 2) 
{
  conf.high <- conf.low <- estimate <- NULL
  stopifnot(class(model) == "coxph")
  data <- eval(model$call$data)
  terms <- attr(model$terms, "dataClasses")[-1]
  coef <- as.data.frame(tidy(model))
  gmodel <- glance(model)
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
  toShowExpClean <- toShowExpClean[nrow(toShowExpClean):1, 
                                   ]
  rangeb <- range(toShowExpClean$conf.low, toShowExpClean$conf.high, 
                  na.rm = TRUE)
  breaks <- axisTicks(rangeb/2, log = TRUE, nint = 7)
  rangeplot <- rangeb
  rangeplot[1] <- rangeplot[1] - diff(rangeb)
  rangeplot[2] <- rangeplot[2] + 0.15 * diff(rangeb)
  width <- diff(rangeplot)
  y_variable <- rangeplot[1] + cpositions[1] * width
  y_nlevel <- rangeplot[1] + cpositions[2] * width
  y_nn <- rangeplot[1] + cpositions[3] * width
  y_cistring <- rangeplot[1] + cpositions[4] * width
  y_stars <- rangeplot[1] + cpositions[5] * width
  x_annotate <- seq_len(nrow(toShowExpClean))
  annot_size_mm <- fontsize * as.numeric(convertX(unit(theme_get()$text$size, "pt"), "mm"))
  
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
                  breaks = breaks) + theme_light() + theme(panel.grid.minor.y = element_blank(),
                                                           panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(),
                                                           legend.position = "none", panel.border = element_blank(),
                                                           axis.title.y = element_blank(), axis.text.y = element_blank(),
                                                           axis.ticks.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
    xlab("") +
    annotate(geom = "text", x = x_annotate, y = exp(y_variable),
             label = toShowExpClean$var, fontface = "bold", hjust = 0,
             size = annot_size_mm) +
    annotate(geom = "text", x = x_annotate,
             y = exp(y_nlevel), hjust = 0, label = ifelse(toShowExpClean$level=="","(continuous)",as.character(toShowExpClean$level)),
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
             label = paste0("# Events: ", gmodel$nevent, "; Global p-value (Log-Rank): ", format.pval(gmodel$p.value.log, eps = ".001"),
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

cairo_pdf("figures/coxph.1.pdf",height = 10/2.54,width=23/2.54,pointsize = 12*0.7)
ggforest2(bigcoxmodel2,main = "",fontsize =1,cpositions = c(0.02, 0.18, 0.28,0.33,0.94))
dev.off()




