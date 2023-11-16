# meta.dt: meta data, simmilar to the excel file I send before, including time, group and techrep (like donor)
# mat: input matrix

# important: matrix column order MUST match the meta.dt rows
# check:
identical(rownames(meta.dt), colnames(mat)) # must return TRUE



time <- as.numeric(factor(meta.dt$time, levels = c("t1", "t2", "t3", "t4", "t5"))) # time points in the right order
group <- factor(meta.dt$group, levels = c("mg", "ctrl"))
donor <- factor(meta.dt$donor)
X <- splines::ns(time, df = 3) # spline fit: degree of freedom. usually 3 to 5, depending on the number of time points


## Donor as FIXED EFFECT (if design is PAIRED)
#####################
design <- model.matrix(~ group * X + Donor)
colnames(design) <- make.names(colnames(design))
fit <- lmFit(mat, design)
fit <- eBayes(fit)
#####################


## Donoras as MIXED EFFECT (if design is mot PAIRED)
#####################
design <- model.matrix(~ group * X)
dupcor <- duplicateCorrelation(mat, design, block = donor)
colnames(design) %<>% make.names(.)
fit0 <- lmFit(mat, 
              design,
              block = donor,
              correlation = dupcor$consensus.correlation)
fit <- eBayes(fit)
