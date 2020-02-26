--------------------------------------------------------------------------------
"supplementary_01.R"

# utility
library(tidyverse)
histo_pal = c("#631879FF", "#3B4992FF", "#5F559BFF", "#A20056FF", "#BB0021FF", "#EE0000FF", "#008B45FF", "#008280FF")
gtf2i_pal = c(m="#3C5488FF",w="#E64B35FF",c="#00A087FF")

################################################################################
  "Supplementary table 1"
################################################################################
meta_dt <- read_tsv('~sypark/00_Project/01_thymoma/10_Final_data/02_metadata/Thymoma_summary_191204_1stSheet.txt')
meta_dt_clean <- meta_dt %>% 
  transmute(ID=id, Histologic_type=histologic_type, Cohort=cohort, Sex=gender, 
            Age_at_diagnosis=age_at_diagnosis, Massaoka_koga_stage = Stage,
            `Recurrence(0=no,1=yes)` = Recur_event,
            `Recur_free_survival_time(day)`=RecurFreeSurvival,
            `Overall_survival(0=alive,1=death)`=Death_event,
            `Overall_survival_time(day)`=DeathSurvival,
            History_of_myasthenia_gravis=history_myasthenia_gravis,
            GTF2I_L424H_final=ifelse(GTF2I_status2=="m","m","w"),
            GTF2I_L424H_TCGA_reported=TCGA_paper_GTF2Imt,
            n_sub_prv,
            n_sub,
            n_indel,
            n_pointmt,
            bait_size,
            final_cellularity,
            final_ploidy)
library(rJava,lib.loc="/home/users/kjyi/R/x86_64-redhat-linux-gnu-library/3.6/")
library(xlsx,lib.loc="/home/users/kjyi/R/x86_64-redhat-linux-gnu-library/3.6/")
meta_dt_clean %>% xlsx::write.xlsx("tables/metadata_sup_table_1.xlsx")




################################################################################
  "Supplementary Fig1a"
################################################################################
"correlation between tumor cell fraction and the number of point mutation"

meta_dt$histologic_type
meta_dt$GTF2I_status2

# cairo_pdf("figures/purity_nummutation.pdf",height = 10/2.54,width=10/2.54,pointsize = 12*0.7)
# par(pty='s')
# with(meta_dt,
#      plot(final_cellularity, n_pointmt/bait_size,
#           ylim=c(0,2),cex=1.5,
#           pch=21,
#           bg=gtf2i_pal[GTF2I_status2],
#           xlab="Tumor cell fraction",
#           ylab="N. of  point mutation/whole exome-seq bait size")
#      )
# dev.off()

################################################################################
"Supplementary Fig1a"
################################################################################
load("~/Projects/thymus_single_cell/final2/data/snapshop_for_volcanoplot.RData")

cairo_pdf("figures/vocano_per_cohort.pdf",height = 10/2.54,width=20/2.54,pointsize = 12*0.7)
par(mfrow = c(1,2), pty="s")
for (cohort in c("TCGA","SNU")) {
  volc_dt %>%
    select(gene,starts_with(cohort),role_in_cancer) %>%
  {
    MT_ids_0 = MT_ids[grepl(cohort, MT_ids)];
    WT_ids_0 = WT_ids[grepl(cohort, WT_ids)];
    .$MT_mean = rowMeans(.[,MT_ids_0]);
    .$WT_mean = rowMeans(.[,WT_ids_0]);
    .$t.test.pvalue <- apply(.[,c(MT_ids_0, WT_ids_0)],
                             1,
                             function(x) my.t.test.p.value(x[MT_ids_0], x[WT_ids_0]));
    WT_high_genes <- filter(., WT_mean - MT_mean >=1 & t.test.pvalue < 0.05/nrow(.)) %>% .$gene;
    with(.,{
      x=WT_mean-MT_mean
      y=-log10(p.adjust(t.test.pvalue, method = "BH"))
      y[is.na(y)] = 0
      labelat = (abs(x)>2) | (y > 40)
      col=ifelse(is.na(.$role_in_cancer), "grey","red")
      greys=col=="grey"
      reds=col=="red"
      plot(x[greys], y[greys],
           xlim=c(-max(abs(x)),max(abs(x)))*1.25,
           ylim=c(0,max(y))*1.15,
           bg="#00000020", col="#00000030", pch=21, ylab="", xlab="")
      points(x[reds],y[reds],bg="#FF000060",col="#FF000070",pch=pch)
      title(cohort)
      title(ylab=expression(-log[10]~Adj.p), line=2)
      title(xlab=expression(log[10]~FC), line=2)
      set.seed(42)
      maptools::pointLabel(x=x[labelat],y=y[labelat],labels=gene[labelat],method = "GA")
      # legend(-1.25*max(abs(x)),1.15*max(y),legend="in GTF2I-mutant",pch="↑",xjust=0,yjust=1,bty='n',text.font = 2)
      # legend(1.25*max(abs(x)),1.15*max(y),legend="in wild-type",pch="↑",xjust=1,yjust=1,bty='n',text.font = 2)
      legend(1.25*max(abs(x)),1.15*max(y),
             legend = c("known oncogene"),
             pt.bg=c("#FF000060"),
             col=c("#FF000070"),
             pch=21,
             xjust=1,yjust=1,bty='n'
      )
      abline(h=-log10(0.05),lty=2)
      text(par()$usr[1],-log10(0.05/19900),"Adj.p > 0.05",adj=c(-0.1,-0.1))
    }
    )
  }
}
dev.off()


################################################################################
"Supplementary Fig2e"
################################################################################
"Sequenza-like plot, example 2+2"

library(tidyverse,lib.loc="/home/users/kjyi/R/x86_64-redhat-linux-gnu-library/3.6/")
library(stringr,lib.loc="/home/users/kjyi/R/x86_64-redhat-linux-gnu-library/3.6/")
# library(sequenza,lib.loc="/home/users/kjyi/R/x86_64-redhat-linux-gnu-library/3.6/") 
library(sequenza2,lib.loc="/home/users/sypark/R/x86_64-redhat-linux-gnu-library/3.6/") 

#chromsome X ----
hg19Bands <- as.data.frame(matrix(c("chrX", 0,  4300000,  "p22.33", "gneg",
                                    "chrX", 4300000,  6000000,  "p22.32", "gpos50",
                                    "chrX", 6000000,  9500000,  "p22.31", "gneg",
                                    "chrX", 9500000,  17100000, "p22.2",  "gpos50",
                                    "chrX", 17100000, 19300000, "p22.13", "gneg",
                                    "chrX", 19300000, 21900000, "p22.12", "gpos50",
                                    "chrX", 21900000, 24900000, "p22.11", "gneg",
                                    "chrX", 24900000, 29300000, "p21.3",  "gpos100",
                                    "chrX", 29300000, 31500000, "p21.2",  "gneg",
                                    "chrX", 31500000, 37600000, "p21.1",  "gpos100",
                                    "chrX", 37600000, 42400000, "p11.4",  "gneg",
                                    "chrX", 42400000, 46400000, "p11.3",  "gpos75",
                                    "chrX", 46400000, 49800000, "p11.23", "gneg",
                                    "chrX", 49800000, 54800000, "p11.22", "gpos25",
                                    "chrX", 54800000, 58100000, "p11.21", "gneg",
                                    "chrX", 58100000, 60600000, "p11.1",  "acen",
                                    "chrX", 60600000, 63000000, "q11.1",  "acen",
                                    "chrX", 63000000, 64600000, "q11.2",  "gneg",
                                    "chrX", 64600000, 67800000, "q12",  "gpos50",
                                    "chrX", 67800000, 71800000, "q13.1",  "gneg",
                                    "chrX", 71800000, 73900000, "q13.2",  "gpos50",
                                    "chrX", 73900000, 76000000, "q13.3",  "gneg",
                                    "chrX", 76000000, 84600000, "q21.1",  "gpos100",
                                    "chrX", 84600000, 86200000, "q21.2",  "gneg",
                                    "chrX", 86200000, 91800000, "q21.31", "gpos100",
                                    "chrX", 91800000, 93500000, "q21.32", "gneg",
                                    "chrX", 93500000, 98300000, "q21.33", "gpos75",
                                    "chrX", 98300000, 102600000,  "q22.1",  "gneg",
                                    "chrX", 102600000,  103700000,  "q22.2",  "gpos50",
                                    "chrX", 103700000,  108700000,  "q22.3",  "gneg",
                                    "chrX", 108700000,  116500000,  "q23",  "gpos75",
                                    "chrX", 116500000,  120900000,  "q24",  "gneg",
                                    "chrX", 120900000,  128700000,  "q25",  "gpos100",
                                    "chrX", 128700000,  130400000,  "q26.1",  "gneg",
                                    "chrX", 130400000,  133600000,  "q26.2",  "gpos25",
                                    "chrX", 133600000,  138000000,  "q26.3",  "gneg",
                                    "chrX", 138000000,  140300000,  "q27.1",  "gpos75",
                                    "chrX", 140300000,  142100000,  "q27.2",  "gneg",
                                    "chrX", 142100000,  147100000,  "q27.3",  "gpos100",
                                    "chrX", 147100000,  155270560,  "q28",  "gneg"), byrow=T, ncol=5), stringsAsFactors=F)
hg19Bands$V2=as.numeric(hg19Bands$V2)
hg19Bands$V3=as.numeric(hg19Bands$V3)
hg19Bands$V6=unlist(plyr::revalue(hg19Bands$V5, list(gneg = "white", gpos25 = "grey25", gpos50 = "grey50", gpos75 = "grey75", gpos100 = "black", acen = NA)))

# functions ----

DrawChromosome <- function(x,y, width, chrBands, horiz = F) {
  maxY <- max(chrBands$V3)
  aceni <- which(chrBands$V5 == "acen")
  if (length(aceni) == 2) {
    if(horiz) {
      allX <- c(0, 0, (width/2),0,0,0+width, 0+width,0+(width/2),0+width,0+width, 0) +y
      allY <- c(0, chrBands$V2[aceni[1]], chrBands$V3[aceni[1]],  chrBands$V3[aceni[2]], maxY, maxY, chrBands$V3[aceni[2]], chrBands$V3[aceni[1]], chrBands$V2[aceni[1]], 0,0)+x
      nPoints <- length(allX)
      
      polygon(x = c(x+chrBands$V2[aceni[1]],x+chrBands$V2[aceni[1]],x+chrBands$V3[aceni[1]]),
              y = c(y, y+width, y + width/2),
              col = "#8B2323", lty = 0)
      polygon(x = c(x+chrBands$V2[aceni[2]],x+chrBands$V3[aceni[2]],x+chrBands$V3[aceni[2]]),
              y = c(y + width/2, y, y+width),
              col = "#8B2323", lty = 0)
      rect(xleft = x+chrBands$V2,
           xright = x+chrBands$V3,
           ybottom = y,
           ytop = y + width, lwd = 0, col = chrBands$V6)
      segments(x0 = allY[1:nPoints-1], 
               y0 = allX[1:nPoints-1],
               x1 = allY[2:nPoints],
               y1 = allX[2:nPoints])
      
    } else {
      allX <- c(0, 0, (width/2),0,0,0+width, 0+width,0+(width/2),0+width,0+width, 0) +x
      allY <- c(0, chrBands$V2[aceni[1]], chrBands$V3[aceni[1]],  chrBands$V3[aceni[2]], maxY, maxY, chrBands$V3[aceni[2]], chrBands$V3[aceni[1]], chrBands$V2[aceni[1]], 0,0)+y
      nPoints <- length(allX)
      
      polygon(y = c(y+chrBands$V2[aceni[1]],y+chrBands$V2[aceni[1]],y+chrBands$V3[aceni[1]]),
              x = c(x, x+width, x + width/2),
              col = "#8B2323")
      polygon(y = c(y+chrBands$V2[aceni[2]],y+chrBands$V3[aceni[2]],y+chrBands$V3[aceni[2]]),
              x = c(x + width/2, x, x+width),
              col = "#8B2323")
      rect(ybottom = y + chrBands$V2,
           ytop = y + chrBands$V3,
           xleft = x,
           xright = x + width, lwd = 0, col = chrBands$V6)
      segments(x0 = allX[1:nPoints-1], 
               y0 = allY[1:nPoints-1],
               x1 = allX[2:nPoints],
               y1 = allY[2:nPoints])
    }
  }
}

chromosome.view2 <- function (baf.windows, ratio.windows, mut.tab = NULL, segments = NULL, 
                              min.N.baf = 1, min.N.ratio = 10000, main = "", vlines = FALSE, 
                              legend.inset = c(-20 * strwidth("a", units = "figure"), 0), 
                              BAF.style = "lines", CNn = 2, cellularity = NULL, ploidy = NULL, 
                              avg.depth.ratio = NULL, model.lwd = 1, model.lty = "24", 
                              model.col = "#00000040", x.chr.space = 10, ylim = c(0,2.5), xlab = "Position (Mb)", par = F, xaxt = NULL) {
  
  dt2 <- function(x, df, ncp, log = FALSE, mean, sd) {
    x2 <- (x - mean) / sd
    dt(x2, df = df, ncp = ncp, log = log)
  }
  expected.baf <- function(sd, ...) {
    baf      <- theoretical.baf(...)
    baf.t2 <- function(BAF, sd, by = 0.001){
      bafs   <- seq(0,1,0.001)
      b.b    <- dt2(x = bafs, mean = BAF, sd = sd, df = 5)
      b.a    <- dt2(x = bafs, mean = 1-BAF, sd = sd, df = 5)
      half.b <- bafs[bafs <= 0.5]
      b <- (b.b+b.a)[bafs <= 0.5]
      weighted.mean(half.b,b)
    }
    BAF <- mapply(FUN = baf.t2, baf$BAF,
                  MoreArgs = list(sd = sd))
    wgh <- dt2(x = baf$BAF, mean = 0.5, sd = sd, df = 5)
    wgh <- wgh/max(wgh)
    mean.bf <- function(x) {
      weighted.mean(x=c(x["BAF"], x["eBAF"]), w = c((1-x["wgh"]), x["wgh"]))
    }
    baf$BAF <- apply(cbind(BAF = baf$BAF, eBAF = BAF, wgh = wgh), 1, FUN = mean.bf)
    baf
  }
  make.polygons <- function(segments, model.baf) {
    max.B <- max(model.baf$B[model.baf$CNt =t= max(segments$CNt)])
    mat.polygs <- matrix(ncol = max.B + 1, nrow = nrow(segments))
    colnames(mat.polygs) <- 0:max.B
    get.B <- function(CNt, B) model.baf$BAF[model.baf$CNt == CNt & model.baf$B == B]
    polyg.coords <- sapply(0:max.B, FUN = function(k) as.numeric(sapply(segments$CNt, FUN = function(i) get.B(i, k))))
    polyg.coords[is.na(polyg.coords)] <- 1
    polyg.coords <- cbind(polyg.coords, 1)
    polyg.pos <- segments[, c("start.pos", "end.pos")]
    edge1 <- polyg.pos$end.pos[-nrow(polyg.pos)]
    edge2 <- polyg.pos$start.pos[-1]
    no.dat <- c(1:nrow(polyg.pos))[edge2 - edge1 >= 1000000]
    v.gaps <- apply(rbind(edge1[no.dat], edge2[no.dat]), 2, mean)
    v.gaps <- cbind(start.pos = v.gaps, end.pos = v.gaps)
    polyg.p.new <- rbind(polyg.pos, v.gaps)
    polyg.c.new <- polyg.coords
    for (i in no.dat) {
      polyg.c.new <- rbind(polyg.c.new, 0)
    }
    new_idx <- c(seq_along(polyg.pos$start.pos), no.dat + 0.5)
    polyg.pos <- polyg.p.new[order(new_idx), ]
    polyg.coords <- polyg.c.new[order(new_idx), ]
    Xs <- unlist(lapply(1:nrow(polyg.coords), function(k) polyg.pos[k, ]))
    color <- gray.colors((max.B + 1), start = 0.5, end = 0.9, alpha = 0.3)
    extra.x <- c(max(Xs), min(Xs))
    extra.y = c(0, 0)
    for (k in 1:ncol(polyg.coords)) {
      Ys <- unlist(lapply(polyg.coords[, k], function(i) rep(i, 2)))
      polygon(x = c(Xs, extra.x), y = c(Ys, extra.y), col = color[k], border = NA)
      extra.y <- rev(Ys)
      extra.x <- rev(Xs)
    }
  }
  if (is.null(segments)) {
    data.model <- NULL
  }
  else {
    if ("CNt" %in% colnames(segments)) {
      if (length(c(cellularity, ploidy, avg.depth.ratio)) != 3) {
        data.model <- NULL
      }
      else {
        data.model <- list()
        CNt.max <- max(segments$CNt, na.rm = TRUE) + 1
        CNt.min <- 0
        data.model$baf <- expected.baf(sd = mean(segments$sd.BAF, na.rm = TRUE), CNn = CNn, CNt = CNt.max, cellularity = cellularity)
        if (CNn == 2) {
          data.model$baf <- rbind(c(0, 0, max(data.model$baf$BAF), 0), data.model$baf)
        }
        else {
          data.model$baf <- rbind(c(0, 0, 1, 0), data.model$baf)
        }
        types <- types.matrix(CNt.min = CNt.min, CNt.max = CNt.max, CNn = CNn)
        data.model$muf <- cbind(types, model.points(cellularity = cellularity, ploidy = ploidy, types = types, avg.depth.ratio = avg.depth.ratio))
      }
    }
    else {
      data.model <- NULL
    }
  }
  if(par==T){par(mar = c(5,5,5,5), oma = c(0, 0, 0, 0), xaxt = "n")}
  min.x <- min(c(min(baf.windows$start), min(ratio.windows$start)))
  max.x <- max(c(max(baf.windows$end), max(ratio.windows$end)))
  xlim <- c(min.x, max.x)
  plotWindows(ratio.windows, ylab = "", las = 1, 
              n.min = min.N.ratio, ylim = ylim, yaxt = 'n', xlab = '', bty = 'l',xaxt = "n",m.col = "grey60",q.bg = "grey90")
  axis(side = 4, line = 0, las = 1)
  mtext(text = "Depth ratio", side = 4, line = 3, cex = par("cex.lab") * 
          par("cex"))
  if (!is.null(segments)) {
    if (vlines) {
      abline(v = segments$end.pos, lwd = 1, lty = 2)
    }
    if (!is.null(data.model)) {
      ratios.theoric <- unique(data.model$muf[, c("CNt", 
                                                  "depth.ratio")])
      segments(x0 = rep(min(segments$start.pos, na.rm = TRUE),times = nrow(ratios.theoric)),
               x1 = rep(max(segments$end.pos,na.rm = TRUE), times = nrow(ratios.theoric)),
               y0 = ratios.theoric$depth.ratio, lwd = model.lwd,
               lty = model.lty, col = model.col) # grey dot line background grid
      segments(x0 = rep(min(segments$start.pos, na.rm = TRUE),times = nrow(ratios.theoric)), 
               x1 = rep(max(segments$end.pos,na.rm = TRUE), times = nrow(ratios.theoric)), 
               y0 = ratios.theoric$depth.ratio[ratios.theoric$CNt==2], lwd = model.lwd, 
               lty = model.lty, col = model.col) # CN == 2 line 
      axis(labels = as.character(ratios.theoric$CNt), side = 2, 
           line = 0, las = 1, at = ratios.theoric$depth.ratio)
      mtext(text = "Copy number", side = 2, line = 3, cex = par("cex.lab") * 
              par("cex"))
    }
    twoN = ratios.theoric$depth.ratio[ratios.theoric$CNt==2]
    segments(x0 = segments$start.pos, y0 = segments$depth.ratio,
             x1 = segments$end.pos, y1 = segments$depth.ratio,
             col = ifelse(segments$depth.ratio>1,"red","black"), lwd = 3) # segment line
  }
  par(xaxt = "s")
  if(!xaxt == 'n') {
    axis(labels = as.character(round(seq(xlim[1]/1000000, xlim[2]/1000000, by = x.chr.space), 0)), side = 1, line = 0, at = seq(xlim[1], xlim[2], by = 1000000 * x.chr.space), outer = FALSE, cex = par("cex.axis") * par("cex"))
  }
  mtext(xlab, side = 1, line = 3, outer = FALSE, cex = par("cex.lab") * par("cex"))
  # mtext(main, 3, outer = TRUE, cex = par("cex.main") * par("cex"), line = 2)
  # title(main, cex = par("cex.main") * par("cex"))
}

plot_sqz_chr_view <- function(sqz_result_path, ylim = c(0,2.5), xaxt = "s", xlab = "") {
  tmpenv <- new.env()
  # tmpenv2 <- new.env()
  load(dir(sqz_result_path, "sequenza_extract.RData",full.names = T)[1], envir=tmpenv)
  # load(dir(sqz_result_path, "sequenza_cp_table.RData",full.names = T)[1], envir=tmpenv2)
  seg.tab <- read_tsv(dir(sqz_result_path, "segments.txt",full.names = T)[1])
  solution <- read_tsv(dir(sqz_result_path, "confints_CP.txt",full.names = T)[1])[2,]
  test <- tmpenv[[ls(envir=tmpenv)[1]]]
  avg.depth.ratio <-test$avg.depth.ratio
  if(is.null(avg.depth.ratio)){avg.depth.ratio <- mean(test$gc$adj[, 2])}
  chrXname= ifelse("X" %in% names(test$mutations), "X","chrX")
  
  # test$mutations[["X"]] <- na.omit(test$mutations[["X"]])
  # test$BAF[["X"]] <- na.omit(test$BAF[["X"]])
  # test$ratio[["X"]] <- na.omit(test$ratio[["X"]])
  chromosome.view2(mut.tab = test$mutations[[chrXname]], baf.windows = test$BAF[[chrXname]],
                   ratio.windows = test$ratio[[chrXname]], min.N.ratio = 1,
                   segments = seg.tab[seg.tab$chromosome == chrXname,],
                   main = "X",
                   cellularity = solution$cellularity, ploidy = solution$ploidy.estimate,
                   avg.depth.ratio = avg.depth.ratio, ylim = ylim, xlab = xlab,
                   model.col = "#00000060", xaxt = xaxt,model.lwd = 0.6)
}

plot_sqz_chr_view("data/sequenza_wes/SNU/SNU_14_C_sequenza_0.45_2.3/")
plot_sqz_chr_view("data/sequenza_wes/TCGA/TCGA-X7-A8M3")
sqz_result_path=
# sequenza plot ----
# cairo_pdf("figures/sup.fig2e.sequenza.pdf",height = 400/72,width=400/72,pointsize = 12*0.7)
# pdf("figures/sup.fig2e.sequenza.pdf",height = 400/72,width=400/72,pointsize = 12*0.7)
pdf("figures/sup.fig2e.sequenza.pdf",height = 75/25.4,width=75/25.4,pointsize = 12*0.7)
# plot(1, xlim = c(1,155270560),ylim = c(0,1), bty = "n", xaxt = "n", yaxt = "n", col = "white", ylab = "", main = "Chromosome X")
# DrawChromosome(0,0.4,0.6, hg19Bands, horiz = T)
# dev.off()

par(mar = c(0,5,1,5),oma=c(3,0,0,0))
layout(matrix(c(1,2,3,4), 4, 1, byrow = TRUE), widths=c(1), heights=c(1,2,2,2))
# par(mfrow=c(1,3))
plot(1, xlim = c(1,155270560),ylim = c(0,1), bty = "n", xaxt = "n", yaxt = "n", col = "white", ylab = "", main = "Chromosome X")
DrawChromosome(0,0.4,0.6, hg19Bands, horiz = T)
text(106200000,0.1, "IRS4")
segments(x0 = 106200000, y0 = 0.4, y1 = 1, col = "red", lty = 1, lwd = 2)
par(mar = c(1,5,1,5))
# plot_sqz_chr_view(sqz_result_path="data/sequenza_wes/SNU/SNU_14_C_sequenza_0.45_2.3/",
#                   ylim = c(0.8,1.6), xaxt = 'n')
# title("SNU_14_C", adj = 0)
# segments(x0 = 106200000, y0 = 0.45, y1 = 2.4, col = "red", lty = 3, lwd = 2)
# text(x = 106200000+1000000, y = 1.4, "7 copies", adj = c(0,0), font = 2)
plot_sqz_chr_view(sqz_result_path="data/sequenza_wes/SNU/SNU_17_C_sequenza_0.5_2.3/",
                  ylim = c(0.7,1.3), xaxt = 'n')
abline(v = 106200000, col = "red", lty = 3, lwd = 2)
title("SNU_17_C", adj = 0)
text(x = 106200000+1000000, y = 1.13, "5 copies", adj = c(0,0), font = 2)
plot_sqz_chr_view(sqz_result_path="data/sequenza_wes/SNU/SNU_18_C_sequenza/",
                  ylim = c(0.0,1.8), xaxt = 'n')
title("SNU_18_C", adj = 0)
abline(v = 106200000, col = "red", lty = 3, lwd = 2)
text(x = 106200000+1000000, y = 1.39, "3 copies", adj = c(0,0), font = 2)
# dev.off()
# TCGA-XM-A8RC
# plot_sqz_chr_view(sqz_result_path="data/sequenza_wes/TCGA/TCGA-X7-A8M3", ylim = c(0,2.5), xaxt = 'n')
# segments(x0 = 10000,y0  = 1.39, x1 = 20000)
# title("TCGA-X7-A8M3", adj = 0)
# abline(v = 106200000, col = "red", lty = 3, lwd = 2)
# text(x = 106200000+1000000, y = 1.43, "4 copies", adj = c(0,0), font = 2)
# plot_sqz_chr_view(sqz_result_path="data/sequenza_wes/TCGA/TCGA-XU-A931", xaxt = 'n')
plot_sqz_chr_view(sqz_result_path="data/sequenza_wes/TCGA/TCGA-XU-A931", ylim = c(0.3,2.3), xaxt = 'n')
segments(x0 = 10000,y0  = 1.39, x1 = 20000)
title("TCGA-XU-A931", adj = 0)
abline(v = 106200000, col = "red", lty = 3, lwd = 2)
text(x = 106200000+1000000, y = 1.43, "4 copies", adj = c(0,0), font = 2)

dev.off()



