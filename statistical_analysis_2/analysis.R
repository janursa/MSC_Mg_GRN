

library(MSqb) # in-house
library(RNAqb) # in-house
library(ggplot2)
library(ggsci)
library(RColorBrewer)
library(ComplexHeatmap)
library(FactoMineR)
library(factoextra)
library(corrplot)
library(ggrepel)
library(imputeLCMD)


## read data
MSqb:::.mk.dir(path = getwd(),
               sub.dir = c("QCplots", "PCA", "DiffExp_tables", "ProfilePlots", "VennDiagram", "ProteinAbundance_tables"))

library(readr)
dtw <- read_csv("Data/primary_omics.csv" ) %>%
   setDT() %>% .[,-1]



## remove protein duplicates
dp <- c("Q8WZ42", "Q9NVA2", "Q9P2E9", "P41252", "P08727")
dtw %<>% .[!Entry %in% dp]



## convert to long format
dt <- melt(dtw, id.vars = "Entry", value.name = "Abundance", variable.name = "Sample") %>%
   .[, c("Condition", "Day") := tstrsplit(Sample, "_")] %>%
   .[, Day := sprintf("day%02s", Day)] %>%
   setnames(., "Entry", "Protein") %>%
   .[, Time := gsub("day", "", Day) %>% as.numeric] %>%
   char2fact()



## generate phenodata, categorical variables and colors
pdt <- dt[, .(Sample, Condition, Day)] %>% unique
Categorical.vars <- pdt[, -c("Sample")]
manual.color.palette <- list(Condition = "npg", Day = "Set1")
global.colors <- setColors(dt = Categorical.vars,
                           manual.color.palette = manual.color.palette)


## intensity dist. before normalisation
# genHM_denc(dt,
#            val = "Abundance",
#            row.p = "Protein",
#            col.p = names(manual.color.palette),
#            log.trans = FALSE,
#            # bar.palettes = manual.color.palette,
#            hm.palette = "Spectral",
#            fill.box.by = "Condition",
#            sp = "_",
#            biorep = NULL,
#            techrep = "Sample",
#            box.outline = FALSE,
#            box.height = unit(3, "cm"),
#            # heatmap_legend_side = "bottom",
#            annotation_legend_side = "right",
#            column_title = "Density heatmap",
#            filename = "IntensityDensity_PreNormalisation",
#            save.plot = F,
#            suffix = NULL,
#            save.format = "png",
#            hm.path = QCplots.path,
#            width = 7, height = 6,
#            column_names_fontsize = 7)


## PCA before normalisation
# pca.p <- plot_PCA(dt,
#                   sampleid = "Sample",
#                   val = "Abundance",
#                   row.p = "Protein",
#                   col.p = names(manual.color.palette),
#                   topNperc = 100,
#                   sep.col = "_._",
#                   biorep = NULL,
#                   techrep = "Animal",
#                   label = c("Condition", "Day"),
#                   labelsize = 2,
#                   pca.dims = c(1, 2),
#                   complete.rows = TRUE,
#                   repel = TRUE,
#                   manual.color.palette = manual.color.palette,
#                   legend.position = "right",
#                   legend.box = "vertical",
#                   legend.direction = "vertical",
#                   prefix = "",
#                   suffix = "preNorm",
#                   width = 7,
#                   height = 7,
#                   saveplot = TRUE,
#                   plot.seperate = FALSE,
#                   exclude.row.p =  "Sample",
#                   preserve.aspect.ratio = TRUE,
#                   save.format = "PDF",
#                   plot.extra.PCA = TRUE,
#                   path = file.path(PCA.path, "preNorm") )




pdt[, Time := as.numeric(gsub("day", "", Day))] %>% .[, Day := NULL]
mt <- dcast(dt, Protein ~ Sample, value.var = "Abundance") %>%
   data.frame(., row.names = "Protein" ) %>% as.matrix()

nmt <- limma::normalizeMedianValues(mt)
dt <- nmt %>% data.table(., keep.rownames = "Protein") %>%
   melt(., id.vars = "Protein", variable.name = "Sample", value.name = "Norm.Abundance") %>%
   merge(dt, .)


## intensity dist. after normalisation
# genHM_denc(dt,
#            val = "Norm.Abundance",
#            row.p = "Protein",
#            col.p = names(manual.color.palette),
#            log.trans = FALSE,
#            # bar.palettes = manual.color.palette,
#            hm.palette = "Spectral",
#            fill.box.by = "Condition",
#            sp = "_",
#            biorep = NULL,
#            techrep = "Sample",
#            box.outline = FALSE,
#            box.height = unit(3, "cm"),
#            # heatmap_legend_side = "bottom",
#            annotation_legend_side = "right",
#            column_title = "Density heatmap",
#            filename = "IntensityDensity_PostNormalisation",
#            save.plot = TRUE,
#            suffix = NULL,
#            save.format = "png",
#            hm.path = QCplots.path,
#            width = 7, height = 6,
#            column_names_fontsize = 7)


## PCA before normalisation
# pca.p <- plot_PCA(dt,
#                   sampleid = "Sample",
#                   val = "Norm.Abundance",
#                   row.p = "Protein",
#                   col.p = names(manual.color.palette),
#                   topNperc = 100,
#                   sep.col = "_._",
#                   biorep = NULL,
#                   techrep = "Animal",
#                   label = c("Condition", "Day"),
#                   labelsize = 2,
#                   pca.dims = c(1, 2),
#                   complete.rows = TRUE,
#                   repel = TRUE,
#                   manual.color.palette = manual.color.palette,
#                   legend.position = "right",
#                   legend.box = "vertical",
#                   legend.direction = "vertical",
#                   prefix = "",
#                   suffix = "preNorm",
#                   width = 7,
#                   height = 7,
#                   saveplot = TRUE,
#                   plot.seperate = FALSE,
#                   exclude.row.p =  "Sample",
#                   preserve.aspect.ratio = TRUE,
#                   save.format = "PDF",
#                   plot.extra.PCA = TRUE,
#                   path = file.path(PCA.path, "postNorm") )



## check missingness
# genHM_NA(dt,
#          value = "Abundance",
#          variable = "Protein",
#          sampleID = "Sample",
#          annotation.bar.para = names(manual.color.palette),
#          annotation.bar.palettes = manual.color.palette,
#          body.colors = c("grey", "darkred"),
#
#          bar.height = unit(5, "cm"),
#          heatmap_legend_side = "right",
#          annotation_legend_side = "right",
#          column_title = "Heatmap of missingness",
#          prefix = NULL,
#          suffix = NULL,
#          save.format = "png",
#          hm.path = QCplots.path,
#          width = 5, height = 5,
#          save.plot = TRUE,
#          filename = "Global_missingness",
#          show.plot = TRUE)

# check missingness per protein
dtw2 <- dcast(dt, Protein ~ Sample, value.var = "Abundance")
mtw <- data.frame(dtw2, row.names = "Protein") %>% as.matrix()
prot.miss <- apply(mtw, 1, \(x) 100 * round( sum(is.na(x)) / ncol(mtw), 2) )


# png("QCplots/Missingness_per_protein.png")
# prot.miss %>%
#    hist(., main = "missingness per protein", xlab = "missingness [%]", ylab = "Proteins count")
# dev.off()

# proteins with more than 50% missingness
NAprot <- prot.miss[prot.miss > 50]


# remove proteins with more than 50% missingness
dt %<>% .[!Protein %in% names(NAprot)]
nmt <- nmt[which(!rownames(nmt) %in% names(NAprot)),]




# imoute the rest of the missing values
topList <- list()
dtList <- list()
topN <- 20
imp.met <- c("noImputation", "MinProb", "MLE", "SVD", "KNN", "QRILC")
imp.form <- c(expr(imp.mat0),
              expr(imputeLCMD::impute.MinProb(imp.mat0, q = 0.1)),
              expr(imputeLCMD::impute.wrapper.MLE(imp.mat0) ),
              expr(imputeLCMD::impute.wrapper.SVD(imp.mat0, K = 2) ),
              expr(imputeLCMD::impute.wrapper.KNN(imp.mat0, K = 10) ),
              expr(imputeLCMD::impute.QRILC(imp.mat0) %>% .[[1]] ) )

for (im in seq_along(imp.met)) {


   # make matrix from imputed data
   imp.mat0 <- l2w(
      dl = dt,
      col.p = c("Condition", "Day"),
      row.p = "Protein",
      val = "Norm.Abundance", sp = "&.&.&",
      asMat = TRUE,
      Mat.rownm = "Protein"
   )

   imp.mat <- eval(imp.form[[im]])

   dtl <- w2l(dw = imp.mat, col.p = c("Condition", "Day"), row.p = "Protein", val = "Norm.Abundance", sp = "&.&.&")

   dt.imp <-  merge(dt[, -"Norm.Abundance"], dtl, by = c("Protein", "Condition", "Day"))
   dt.imp[, Imputed := ifelse(is.na(Abundance), TRUE, FALSE)]
   dtList[[ imp.met[im] ]] <- dt.imp

   inmt <- l2w(
      dl = dt.imp,
      col.p = "Sample",
      row.p = "Protein",
      val = "Norm.Abundance",
      asMat = TRUE,
      Mat.rownm = "Protein"
   )

   write.csv(inmt,
             file = file.path("ProteinAbundance_tables",
                              paste0("ProtAbundance__Norm_n_Imp_", imp.met[im],".csv")))



   ## stat test
   df <- 3
   i <- 2 + df
   x = splines::ns(pdt$Time, df = df)
   gr <-  factor(pdt$Condition)


   # method1: Limma guide: https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf
   design1 <- model.matrix(~ gr * x) # sae result as 0 + gr * x
   fit01 <- lmFit(inmt, design1)
   fit1 <- eBayes(fit01, )
   f <- ncol(fit1)
   toptb1 <- topTable(fit1, coef = (i+1):f, number = Inf) %>% data.table(., keep.rownames = "Protein" )
   topList[[ imp.met[im] ]] <- toptb1[, method := imp.met[im]]

   ## save time course results
   write.csv(toptb1[, .(Protein, P.Value, adj.P.Val)],
             file = file.path("DiffExp_tables", paste0("TopTable_TimeCourse_", imp.met[im],".csv")))





   ## profile plots
   inD <- dt.imp[Protein %in% toptb1$Protein[1:topN]]
   inDimp.ctr <- inD[which(Imputed)] %>% .[Condition == "ctr"]
   inDimp.mg <- inD[which(Imputed)] %>% .[Condition == "mg"]


   pl <- ggplot(
      inD,
      aes(x = Day, y = Norm.Abundance, col = Condition)) +
      geom_point() +
      geom_point(data = inDimp.mg,
                 mapping = aes(x = Day, y = Norm.Abundance),
                 shape = 21, col = "brown", size = 2.5) +
      geom_point(data = inDimp.ctr,
                 mapping = aes(x = Day, y = Norm.Abundance),
                 shape = 21, col = "grey20", size = 2.5) +

      # geom_line(alpha = 0.9, aes(col = Condition, group = Condition)) + #
      # facet_grid(rows = vars(cluster), cols = vars(Condition)) +
      facet_wrap(facets = "Protein") +
      theme_bw() +

      scale_colour_manual(
         values = c("#999999", "orange2"),
         breaks = c("ctr", "mg")
      ) +
      geom_spline(aes(x = Day, y = Norm.Abundance, group = Condition, col = Condition),
                  df = df, size = 1.
      ) +
      guides(alpha = "none") +
      theme(axis.text.x = element_text(angle = 90, size = 8))

   ggsave(pl, filename = file.path("ProfilePlots", paste0("SplineFit_df3_top",topN,"_", imp.met[im], ".pdf")),
          width = 9, height = 8)

}


# topProt <- lapply(topList, \(x) x[P.Value < 0.05, Protein])

topl <- rbindlist(topList, use.names = T)

vpath <- file.path("VennDiagram")
if (!dir.exists(vpath)) dir.create(vpath, recursive = T)


plot_Venn_rna(topl,
              set.column = "method",
              input.subset = "all",
              venn.group.size = 6,
              p.adj = "P.Value",
              p.adj.cut = 0.05,
              variable = "Protein",
              guide.palette = "Spectral",
              output.plot = "both", # either of "venn", "upset", "both"

              set_size = 0,
              edge_size = 3,
              label = "count",
              label_geom = "label",
              label_alpha = 0.4,
              label_color = "darkred",

              venn.low.col = "grey90",
              venn.high.col = "grey30",

              width.venn = 9,
              height.venn = 6,
              width.upset = 7,
              height.upset = 4,

              path = vpath,
              filename.xls = paste0("VennIntersect___", imp.met[im],".xlsx") )


# # method2: Aron Lun: https://support.bioconductor.org/p/71005/#117540
# design2 <- model.matrix(~ gr + gr:x) # sae result as 0 + gr * x
# fit02 <- lmFit(nmt, design2)
# colnames(design2) %<>% gsub(":", "_", .)
# # cont <- c("grmg_x1 - grctr_x1", "grmg_x2 - grctr_x2", "grmg_x3 - grctr_x3")
# con <- makeContrasts(grmg_x1 - grctr_x1, grmg_x2 - grctr_x2, grmg_x3 - grctr_x3,
#                      levels = design2) # sub ":" with "_"
# fit2 <- contrasts.fit(fit02, contrasts = con) %>% eBayes()
# toptb2 <- topTable(fit2, number = Inf) %>% data.table(., keep.rownames = "Protein" )




