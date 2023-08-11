

library(readr)
library(data.table)
library(magrittr)
library(limma)
library(dplyr)
library(ggformula) # for plotting splines in ggplot




## defince top dir
topdir <- file.path(getwd(), "subsetTime_knn_minprob_20230423")
if (!dir.exists(topdir)) dir.create(topdir)

## subset of the days to be analysed in each itteration
# days.sub <- c("day1_04", "day1_07", "day1_09", "day1_11", "day1_14", "day1_21")
days.sub <- c("day1_11", "day1_21")


## read data
dtw <- readr::read_csv("Data/primary_omics.csv" ) %>% setDT() %>% .[,-1]

## build phenodata
pdt <- data.table(Sample = names(dtw[, -"Entry"])) %>%
   .[, c("Condition", "Time") := tstrsplit(Sample, "_")] %>%
   .[, Time := as.numeric(Time)] %>%
   .[, Day := factor(sprintf("day%02s", Time))]


## remove protein duplicates
dp <- c("Q8WZ42", "Q9NVA2", "Q9P2E9", "P41252", "P08727")
dtw %<>% .[!Entry %in% dp]

## generate a matrix
mt <- data.frame(dtw, row.names = "Entry" ) %>% as.matrix()

## median normalisation
mt <- limma::normalizeMedianValues(mt)

## check missingness per protein
prot.miss <- apply(mt, 1, \(x) 100 * round( sum(is.na(x)) / ncol(mt), 2) )
# hist(prot.miss) $histogram of the missingness

## remove proteins with more than 50% missingness
NAprot <- prot.miss[prot.miss > 50] # this eliminates 1270 proteins
mt %<>% .[which(!rownames(mt) %in% names(NAprot)),]


## generate imputations calls to be applied for NA imputation
imp.form <- c(
   "MinProb" = dplyr::expr(imputeLCMD::impute.MinProb(mt, q = 0.1)),
   "KNN" = dplyr::expr(imputeLCMD::impute.wrapper.KNN(mt, K = 15) )
)



# imoute the rest of the missing values
topList <- list()
dtList <- list()

## empty datatable to store data in
MasterList <- data.table()



## loop over short. long.term
for (analysis in days.sub) {

   if (analysis == "day1_11") {
      days <-  sprintf("day%02s", c(1:4, 7:11))
      period <- "short.term"
   } else if (analysis == "day1_21") {
      days <-  sprintf("day%02s", c(1:4, 7:11, 14, 21))
      period <- "long.term"
   }

   ## loop over imputation method
   for (im in names(imp.form)) {

      ## apply imputation
      imp.mat <- eval(imp.form[[im]]) # imp.form defined above




      ##__________________________ time-course with limma
      # method1: Limma guide: https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf

      df <- degree.of.freedom <- 3
      i <- 2 + df
      x = splines::ns(pdt$Time, df = df)
      gr <- factor(pdt$Condition)

      design <- model.matrix(~ gr * x)
      fit0 <- lmFit(imp.mat, design)
      fit <- eBayes(fit0, trend = T)
      f <- ncol(fit)
      toptb <- topTable(fit, coef = (i+1):f, number = Inf) %>%
         data.table(., keep.rownames = "Protein") %>%
         .[order(P.Value)]
      topList[[im]] <- toptb[, method := im]

      ## save time course results
      # write.csv(toptb1[, .(Protein, P.Value, adj.P.Val)],
      #           file = file.path(DiffExp_tables.path, paste0("TopTable_TimeCourse_", im,".csv")))
      ##__________________________ end



      ## fill in master list
      MasterList <- lapply(topList, \(x) x[, c("Protein", "P.Value", "adj.P.Val", "method")]) %>%
         rbindlist() %>%
         .[, days := analysis] %>%
         .[, effect := period] %>%
         rbind(MasterList, .)


      ## ___________________profile plots
      # make plot data (ng format)
      # topN = 10
      # dt.plot <- imp.mat %>%
      #    data.table(., keep.rownames = "Protein") %>%
      #    melt(., id.vars = "Protein", value.name = "Abundance", variable.name = "Sample") %>% # to long
      #    .[Protein %in% toptb$Protein[1:topN]] %>%  # only top 20 with smallest p.values
      #    merge(., pdt)
      #
      #
      # pl <- ggplot(
      #    dt.plot,
      #    aes(x = Day, y = Abundance, col = Condition)) +
      #    geom_point() +
      #    # geom_line(alpha = 0.9, aes(col = Condition, group = Condition)) + #
      #    # facet_grid(rows = vars(cluster), cols = vars(Condition)) +
      #    facet_wrap(facets = "Protein") +
      #    theme_bw() +
      #
      #    scale_colour_manual(
      #       values = c("#999999", "orange2"),
      #       breaks = c("ctr", "mg")
      #    ) +
      #    geom_spline(aes(x = Day, y = Abundance, group = Condition, col = Condition),
      #                df = df, size = 1.
      #    ) +
      #    guides(alpha = "none") +
      #    theme(axis.text.x = element_text(angle = 90, size = 8))
      #
      # ggsave(pl, filename = file.path(ProfilePlots.path, paste0("SplineFit_df", df,"_top",topN,"_", im, ".pdf")),
      #        width = 9, height = 8)

   }
}

# save masterList
write.csv(MasterList, file = file.path(topdir, "MasterList.csv"))


