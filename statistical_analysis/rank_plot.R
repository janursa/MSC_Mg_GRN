

library(data.table)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(magrittr)
library(GGally)
library(viridis)
library(ggpubr)


# all, rank product and rank plot
#####################

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# read masterList
mlist <- fread(file = "MasterList_imps_days.csv", drop = 1)

mthd <- c("KNN", "MinProb");

# remove junk proteins (all starting with 'p_')
mstl <- mlist[method %in% mthd] %>%
   .[!Protein %like% "p_"] %>%
   unique()

# order by p, and rank
mstl %<>% .[order(days, method, P.Value)] %>%
   .[, rank := seq_along(P.Value), by = .(method, days)]


# produce rank matrixs
f.rnk <- function(mt) {
   mstl[method == mt] %>%
      dcast(., Protein ~ days, value.var = "rank") %>%
      .[order(day1_04)] %>% # reference: day1-4
      .[, short.term.prod := day1_04 * day1_07 * day1_09] %>% # short term effect
      .[, long.term.prod := day1_11 * day1_14 * day1_21] %>% # long term effect
      na.omit %>% # omit rows containg NA (losing almost 30% of proteins, no work-around)
      .[, short.term.prod.rank := rank(short.term.prod)] %>% # rank base on day-rank products
      .[, long.term.prod.rank := rank(long.term.prod)] %>%
      .[, method := mt]
}

rnk.dt.knn <- f.rnk('KNN')
rnk.dt.minprob <- f.rnk('MinProb')


# combining both methods. Order does not mater, as it's a product
rnk.dt.both <- rbindlist(
   list(rnk.dt.minprob[, .(Protein, short.term.prod.rank, long.term.prod.rank, method)],
        rnk.dt.knn[, .(Protein, short.term.prod.rank, long.term.prod.rank, method)]
        )) %>%
   melt(id.vars = c("Protein", "method")) %>%
   .[, eff := paste(method, variable, sep = "_")] %>%
   .[, c("method", "variable") := NULL] %>%
   dcast(., Protein ~ eff) %>%
   .[, short.term.prod := KNN_short.term.prod.rank * MinProb_short.term.prod.rank] %>%
   .[, long.term.prod := KNN_long.term.prod.rank * MinProb_long.term.prod.rank] %>%
   .[, short.term.prod.rank := rank(short.term.prod)] %>%
   .[, long.term.prod.rank := rank(long.term.prod)]



## save ranked lists
#write.csv(rnk.dt.knn, file = file.path(topdir, paste0("knn_rankedList.csv")))
#write.csv(rnk.dt.minprob, file = file.path(topdir, paste0("MinProb_rankedList.csv")))
#write.csv(rnk.dt.both, file = file.path(topdir, paste0("BothMethods_rankedList.csv")))



#function to plot top n ranked genes
f.rnk.plot <- function(dat, n, cls, gcl, ttl, xlab, ylab) {
   dat %>%
      .[order(short.term.prod.rank)] %>%
      .[1:n, ] %>%
   ggparcoord(.,
              boxplot = F,
              missing = "exclude",
              columns = cls,
              groupColumn = gcl,
              showPoints = TRUE,
              alphaLines = 0.3,
              scale = "std" # look at the help desc
   ) +
      scale_color_viridis(discrete = F) +
      labs(title = ttl, x = xlab, y = ylab) +
      theme_light()
}

ylabel <- "Rank (standardised)" 
xlabel <- ""
title <- ""
# method knn
rpl.knn <-  f.rnk.plot(dat = rnk.dt.knn, n = 100, cls = 10:11, gcl = 10,
                        ttl = title,
                        xlab =  xlabel,
                        ylab = ylabel)
# method MinProb
rpl.minprob <- f.rnk.plot(dat = rnk.dt.minprob, n = 100, cls = 10:11, gcl = 10,
                          ttl = title,
                          xlab =  xlabel,
                          ylab = ylabel)

# both method
rpl.both <- f.rnk.plot(dat = rnk.dt.both, n = 100, cls = 8:9, gcl = 8,
                       ttl = title,
                       xlab =  xlabel,
                       ylab = ylabel)



# save knn and minprob together
twoplots <- ggpubr::ggarrange(
   rpl.knn, rpl.minprob, labels = c("A", "B"),
   common.legend = TRUE, legend = "bottom"
)
ggexport(twoplots, filename = "rankplot_knn_minprob.pdf", width = 7, height = 6)


# save combined methods
ggsave(plot = rpl.both, filename = "rankplot_bothMethods.pdf",
       device = "pdf", width = 6, height = 6)
