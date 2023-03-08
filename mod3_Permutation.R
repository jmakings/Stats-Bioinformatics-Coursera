##### Permutation #####

library(devtools)
library(Biobase)
library(limma)
library(edge)
library(genefilter)

# load data
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
pdata=pData(bot)
edata=as.matrix(exprs(bot))
fdata = fData(bot)

# transform
edata = log2(as.matrix(edata) + 1)
edata = edata[rowMeans(edata) > 10, ]

# calculate t statistic for rows with strain variable
tstats_obj = rowttests(edata,pdata$strain)
hist(tstats_obj$statistic,col=2,xlim=c(-5,2))

# permute the sample labels 
set.seed(135)
strain = pdata$strain
strain0 = sample(strain)
tstats_obj0 = rowttests(edata,strain0)
par(mfrow = c(1,2))
hist(tstats_obj0$statistic,col=2,xlim=c(-5,2))

# view quantiles of these different permutations
quantile(tstats_obj0$statistic)

quantile(tstats_obj$statistic)




