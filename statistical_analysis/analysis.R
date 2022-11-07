

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


## read data
MSqb:::.mk.dir(path = getwd(),
               sub.dir = c("QCplots", "Volcanos", "PCA", "DiffExp_tables",
                           "TopSig_Heatmap", "ProfilePlots", "VennDiagram"))

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
genHM_denc(dt,
           val = "Abundance",
           row.p = "Protein",
           col.p = names(manual.color.palette),
           log.trans = FALSE,
           # bar.palettes = manual.color.palette,
           hm.palette = "Spectral",
           fill.box.by = "Condition",
           sp = "_",
           biorep = NULL,
           techrep = "Sample",
           box.outline = FALSE,
           box.height = unit(3, "cm"),
           # heatmap_legend_side = "bottom",
           annotation_legend_side = "right",
           column_title = "Density heatmap",
           filename = "IntensityDensity_PreNormalisation",
           save.plot = F,
           suffix = NULL,
           save.format = "png",
           hm.path = QCplots.path,
           width = 7, height = 6,
           column_names_fontsize = 7)


## PCA before normalisation
pca.p <- plot_PCA(dt,
                  sampleid = "Sample",
                  val = "Abundance",
                  row.p = "Protein",
                  col.p = names(manual.color.palette),
                  topNperc = 100,
                  sep.col = "_._",
                  biorep = NULL,
                  techrep = "Animal",
                  label = c("Condition", "Day"),
                  labelsize = 2,
                  pca.dims = c(1, 2),
                  complete.rows = TRUE,
                  repel = TRUE,
                  manual.color.palette = manual.color.palette,
                  legend.position = "right",
                  legend.box = "vertical",
                  legend.direction = "vertical",
                  prefix = "",
                  suffix = "preNorm",
                  width = 7,
                  height = 7,
                  saveplot = TRUE,
                  plot.seperate = FALSE,
                  exclude.row.p =  "Sample",
                  preserve.aspect.ratio = TRUE,
                  save.format = "PDF",
                  plot.extra.PCA = TRUE,
                  path = file.path(PCA.path, "preNorm") )




pdt[, Time := as.numeric(gsub("day", "", Day))] %>% .[, Day := NULL]
mt <- dcast(dt, Protein ~ Sample, value.var = "Abundance") %>%
   data.frame(., row.names = "Protein" ) %>% as.matrix()

nmt <- limma::normalizeMedianValues(mt)
dt <- nmt %>% data.table(., keep.rownames = "Protein") %>%
   melt(., id.vars = "Protein", variable.name = "Sample", value.name = "Norm.Abundance") %>%
   merge(dt, .)


## intensity dist. after normalisation
genHM_denc(dt,
           val = "Norm.Abundance",
           row.p = "Protein",
           col.p = names(manual.color.palette),
           log.trans = FALSE,
           # bar.palettes = manual.color.palette,
           hm.palette = "Spectral",
           fill.box.by = "Condition",
           sp = "_",
           biorep = NULL,
           techrep = "Sample",
           box.outline = FALSE,
           box.height = unit(3, "cm"),
           # heatmap_legend_side = "bottom",
           annotation_legend_side = "right",
           column_title = "Density heatmap",
           filename = "IntensityDensity_PostNormalisation",
           save.plot = TRUE,
           suffix = NULL,
           save.format = "png",
           hm.path = QCplots.path,
           width = 7, height = 6,
           column_names_fontsize = 7)


## PCA before normalisation
pca.p <- plot_PCA(dt,
                  sampleid = "Sample",
                  val = "Norm.Abundance",
                  row.p = "Protein",
                  col.p = names(manual.color.palette),
                  topNperc = 100,
                  sep.col = "_._",
                  biorep = NULL,
                  techrep = "Animal",
                  label = c("Condition", "Day"),
                  labelsize = 2,
                  pca.dims = c(1, 2),
                  complete.rows = TRUE,
                  repel = TRUE,
                  manual.color.palette = manual.color.palette,
                  legend.position = "right",
                  legend.box = "vertical",
                  legend.direction = "vertical",
                  prefix = "",
                  suffix = "preNorm",
                  width = 7,
                  height = 7,
                  saveplot = TRUE,
                  plot.seperate = FALSE,
                  exclude.row.p =  "Sample",
                  preserve.aspect.ratio = TRUE,
                  save.format = "PDF",
                  plot.extra.PCA = TRUE,
                  path = file.path(PCA.path, "postNorm") )



## check missingness
genHM_NA(dt,
         value = "Abundance",
         variable = "Protein",
         sampleID = "Sample",
         annotation.bar.para = names(manual.color.palette),
         annotation.bar.palettes = manual.color.palette,
         body.colors = c("grey", "darkred"),

         bar.height = unit(5, "cm"),
         heatmap_legend_side = "right",
         annotation_legend_side = "right",
         column_title = "Heatmap of missingness",
         prefix = NULL,
         suffix = NULL,
         save.format = "png",
         hm.path = QCplots.path,
         width = 5, height = 5,
         save.plot = TRUE,
         filename = "Global_missingness",
         show.plot = TRUE)

# check missingness per protein
dtw2 <- dcast(dt, Protein ~ Sample, value.var = "Abundance")
mtw <- data.frame(dtw2, row.names = "Protein") %>% as.matrix()
prot.miss <- apply(mtw, 1, \(x) 100 * round( sum(is.na(x)) / ncol(mtw), 2) )


png("QCplots/Missingness_per_protein.png")
prot.miss %>%
   hist(., main = "missingness per protein", xlab = "missingness [%]", ylab = "Proteins count")
dev.off()

# proteins with more than 50% missingness
NAprot <- prot.miss[prot.miss > 50]


# remove proteins with more than 50% missingness
dt %<>% .[!Protein %in% names(NAprot)]
nmt <- nmt[which(!rownames(nmt) %in% names(NAprot)),]

# imoute the rest of the missing values
# library(imputeLCMD)
# imp.mat <-  l2w(
#    dl = dt.Liver,
#    col.p = c("Condition", "Animal"),
#    row.p = "Protein",
#    val = "Intensity", sp = "&.&.&",
#    asMat = TRUE,
#    Mat.rownm = "Protein"
# ) %>%
#    imputeLCMD::impute.wrapper.KNN(., K = 15) %>% #.[[1]] %>%
#    w2l(dw = ., col.p = c("Condition", "Animal"), row.p = "Protein", val = "Intensity", sp = "&.&.&")
# dt.Liver <-  merge(dt.Liver[,-"Intensity"], imp.mat, by = c("Protein", "Condition", "Animal"))



library(splines)
## stat test
df <- 3
i <- 2 + df
f <- ncol(fit)
x = splines::ns(pdt$Time, df = df)
gr <-  factor(pdt$Condition)


# method1: Limma guide: https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf
design1 <- model.matrix(~ gr * x) # sae result as 0 + gr * x
fit01 <- lmFit(nmt, design1)
fit1 <- eBayes(fit01, )
toptb1 <- topTable(fit1, coef = (i+1):f, number = Inf) %>% data.table(., keep.rownames = "Protein" )

## save time course results
write.csv(toptb1[, .(Protein, P.Value, adj.P.Val)], file = "TopTable_TimeCourse_final.csv")

# # method2: Aron Lun: https://support.bioconductor.org/p/71005/#117540
# design2 <- model.matrix(~ gr + gr:x) # sae result as 0 + gr * x
# fit02 <- lmFit(nmt, design2)
# colnames(design2) %<>% gsub(":", "_", .)
# # cont <- c("grmg_x1 - grctr_x1", "grmg_x2 - grctr_x2", "grmg_x3 - grctr_x3")
# con <- makeContrasts(grmg_x1 - grctr_x1, grmg_x2 - grctr_x2, grmg_x3 - grctr_x3,
#                      levels = design2) # sub ":" with "_"
# fit2 <- contrasts.fit(fit02, contrasts = con) %>% eBayes()
# toptb2 <- topTable(fit2, number = Inf) %>% data.table(., keep.rownames = "Protein" )



## profile plot
library("ggformula")
{ggplot(
   dt[Protein %in% toptb1$Protein[1:12]],
   aes(x = Day, y = Norm.Abundance, col = Condition)) +
   geom_point() +
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
   theme(axis.text.x = element_text(angle = 90, size = 8))} %>%
   ggsave(filename = "profile_plots_withSpline_df3_top12.pdf", width = 8, height = 7)

