#edgeR tutorial
#from edgeR users guide
#RNA-Seq of pathogen inoculated Arabidopsis with batch effects
setwd("/Users/Cody_2/git.repos/edgeR")

source("http://www.bioconductor.org/biocLite.R")
biocLite("edgeR")
library(edgeR)


#I had trouble compiling the NBPSeq package on my OS so...
#I took arab.rda from the NBPSeq package and placed a copy in the above
#working directory.
load("arab.rda")
head(arab)

#prep data
Treat <- factor(substring(colnames(arab), 1, 4))
head(Treat)
?substring
Treat <- relevel(Treat, ref = "mock")
head(Treat)

Time  <- factor(substring(colnames(arab), 5, 5))
head(Time)

#filter genes that are expressed in both conditions
keep <- rowSums(cpm(arab) > 2) >= 3
arab <- arab[keep, ]
table(keep)

#Create a DGEList and Apply TMM normalization
y <- DGEList(counts = arab, group = Treat)
y <- calcNormFactors(y)
y$samples

#Examine the similarity between samples of the same treatment
plotMDS(y)

design <- model.matrix(~Time+Time:Treat)
head(design)
logFC <- predFC(y, design, prior.count = 1, dispersion = 0.05)
head(logFC)
#are the 3 timepoints correlated as we would expect?

cor(logFC[,4:6])
#the correlation is highest between time 2 and time 3

#to create an additive linear model
design <- model.matrix(~Time+Treat)
rownames(design) <- colnames(y)
design
#prints design matrix and other parts of the design object

#estimating the dispersion
y <- estimateGLMCommonDisp(y, design, verbose = TRUE)
# returns Disp = 0.07051 , BCV = 0.2655 
# where the sqroot of dispersion is the coefficient of Biological variation
# this value is claimed to be on the "high side", given gentically identical 
# individuals used in the experiment

# now estimate the genewise dispersion estimates
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
plotBCV(y)
#at very low logCPM values, the dispersions are very large
# what does one do about this? Just know that it is there?


#Differential Expression
fit <- glmFit(y, design)

# first check if there is a need to adjust for experiment times
lrt <- glmLRT(fit, coef = 2:3)
topTags(lrt)
FDR <- p.adjust(lrt$table$PValue, method = "BH")
sum(FDR < 0.05)
#3143

# now preform test on the last coefficient in the design matrix, aka treatment
lrt <- glmLRT(fit)
topTags(lrt)

#take a look at the top genes across the three replicates
top <- rownames(topTags(lrt))
cpm(y)[top, ]

#the total number of genes sig up or down regulated at 5% FDR:
summary(dt <- decideTestsDGE(lrt))
#prints
# -1  1199
# 0  14088
# 1   1239

#get the DE genes
isDE <- as.logical(dt)
DEnames <- rownames(y)[isDE] 
head(DEnames)

#plot all of the logFC against average count size
plotSmear(lrt, de.tags = DEnames)
abline(h=c(-1, 1), col = "blue")
#blue lines indicate 2 fold change up or down


























