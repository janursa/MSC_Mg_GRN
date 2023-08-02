library(edgeR)
counts <- read.delim("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2018-June-RNA-Seq-Workshop/master/thursday/all_counts.txt")
d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
snames <- colnames(counts) # Sample names
cultivar <- substr(snames, 1, nchar(snames) - 2) 
time <- substr(snames, nchar(snames) - 1, nchar(snames) - 1)
group <- interaction(cultivar, time)
plotMDS(d, col = as.numeric(group))
mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T)
tmp <- voom(d0, mm, plot = T)
fit <- lmFit(y, mm)
contr <- makeContrasts(groupI5.9 - groupI5.6, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr) #Estimate contrast for each gene
tmp <- eBayes(tmp) 
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05)) #How many DE genes are there?
top.table$Gene <- rownames(top.table) #Write top.table to a file
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
write.table(top.table, file = "time9_v_time6_I5.txt", row.names = F, sep = "\t", quote = F)
