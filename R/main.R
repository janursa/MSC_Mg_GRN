library(edgeR)
df <- read.csv("C:/Users/nourisa/Downloads/testProjs/MSC_Mg_omics/data/omics.csv") # remember to import with no space for the first row
# data <- DGEList(data)
# data <- calcNormFactors(data)

df['ctr_11'] <- df['ctr_1']
df['ctr_12'] <- df['ctr_1']+.2
df['Mg_11'] <- df['Mg_1']
df['Mg_12'] <- df['Mg_1']+0.2
df$ctr_1 <- NULL
df$Mg_1 <- NULL
snames <- colnames(df) 
cultivar <- substr(snames, 1, nchar(snames) - 2)
time <- substr(snames, nchar(snames)-1 , nchar(snames)-1 )
group <- interaction(cultivar, time)
design <- model.matrix(~0 + group)
y <- voom(df, design, plot = T)
fit <- lmFit(df, design)
fit <- eBayes(fit)
# fit <- treat(fit,lfc=0.1)
topTable(fit,coef=2)
plotMDS(df, col = as.numeric(group))


