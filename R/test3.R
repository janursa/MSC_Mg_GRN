library(limma)
library(Matrix)

## Create the experimental variables and put them in a data frame
animal <- factor(rep(sprintf("Animal %s", 1:6), each=2))
timep <- factor(c("Pre", "Post"), levels=c("Pre", "Post"))
treat <- factor(c("control", "control", "reagent", "reagent"))
df <- data.frame(animal, timep, treat)
df$replicate <- as.numeric(droplevels(df$animal:df$timep:df$treat))
## Create technical replicates by duplicating each row
df <- df[rep(1:nrow(df), each=2),]
rownames(df) <- as.character(1:nrow(df))
# 
## Manipulate the contrasts of the animal variable. You can view the
## result using print(contrasts(df$animal))
contrasts(df$animal) <- as.matrix(bdiag(contr.sum(3), matrix(nrow=0, ncol=1), contr.sum(3)))
# 
# design <- model.matrix(terms(~0+treat:timep + animal, keep.order=TRUE), df)
# ## Drop the column of all zeros from the design matrix
# design <- design[,colSums(abs(design)) > 0]
# stopifnot(is.fullrank(design))
# ## Make the colnames usable as variable names
# colnames(design) <- make.names(colnames(design))
# 
# ## Generate some random data for demonstration purposes
# x <- matrix(rnorm(nrow(df) * 100), nrow=100, ncol=nrow(df),
#             dimnames=list(gene=sprintf("gene%s", 1:100), sample=rownames(df)))
# 
# ## Fit model
# corfit <- duplicateCorrelation(x, design=design, block=df$replicate)
# fit <- lmFit(x, design, correlation = corfit$consensus.correlation)
# 
# ## Prepare contrast matrices for differential expression tests
# mytests <- lapply(c(
#   control.PreVsPost="treatcontrol.timepPost - treatcontrol.timepPre",
#   reagent.PreVsPost="treatreagent.timepPost - treatreagent.timepPre",
#   Pre.controlVsreagent="treatreagent.timepPre - treatcontrol.timepPre",
#   Post.controlVsreagent="treatreagent.timepPost - treatcontrol.timepPost",
#   Interaction="(treatreagent.timepPost - treatreagent.timepPre) - (treatcontrol.timepPost - treatcontrol.timepPre)"),
#   FUN=function(x) makeContrasts(contrasts=x, levels=design))
# 
# ## Perform each test and collect the results
# results <- lapply(mytests, function(test) topTable(eBayes(contrasts.fit(fit, contrasts=test))))
# 
# ## One data frame for each test
# print(results)