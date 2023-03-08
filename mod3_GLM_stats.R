
##### Generalized linear models #####

#install.packages('MASS')

library(devtools)
library(Biobase)
library(snpStats)
library(broom)
library(MASS)
library(DESeq2)

# get data 
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]

# calculate eigenvectors and pcs 
xxmat <- xxt(sub.10, correct.for.missing=FALSE)
evv <- eigen(xxmat, symmetric=TRUE)
pcs <- evv$vectors[,1:5]

# single logistic regression(additive model)
snpdata = sub.10@.Data
status = subject.support$cc
snp1 = as.numeric(snpdata[,1])

# this is so R knows which values are NA
snp1[snp1==0] = NA
glm1 = glm(status ~ snp1,family="binomial")
tidy(glm1)

# If we want the dominant model
snp1_dom = (snp1 == 1)
glm1_dom = glm(status ~ snp1_dom,family="binomial")

# compare the two
tidy(glm1_dom)
tidy(glm1)

# Adjust for first 5 PCs 
glm2 = glm(status ~ snp1 + pcs[,1:5],family="binomial")
tidy(glm2)

# fit many GLMs at once
glm_all = snp.rhs.tests(status ~ 1,snp.data=sub.10)
slotNames(glm_all)

# qqplot of these chi squared values
qq.chisq(chi.squared(glm_all),df=1)

# new glm but adjusted for PCs
glm_all_adj = snp.rhs.tests(status ~ pcs,snp.data=sub.10)

# this qqplot shows chi squared that matches with expected much better
qq.chisq(chi.squared(glm_all_adj),df=1)

# Poisson/negative binomial regression 

# data
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
pdata=pData(bot)
edata=as.matrix(exprs(bot))
fdata = fData(bot)

# transfrom the data
edata = edata[rowMeans(edata) > 10, ]

# single poisson regression
glm3 = glm(edata[1, ] ~ pdata$strain,family="poisson")
tidy(glm3)

# single negative binomial model
glm.nb1 = glm.nb(edata[1, ] ~ pdata$strain)
tidy(glm.nb1)
# similar estimates but bigger standard error

# use DESeq2 to perform many negative binomial regressions at once
de = DESeqDataSetFromMatrix(edata, pdata, ~strain)
glm_all_nb = DESeq(de)
result_nb = results(glm_all_nb)
hist(result_nb$stat)

##### Calculating statistics in R #####

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

# log transform remove lowly expressed genes
edata = log2(as.matrix(edata) + 1)
edata = edata[rowMeans(edata) > 10, ]

# calculate t statistic rapidly
tstats_obj = rowttests(edata,pdata$strain)
names(tstats_obj)

# histogram of t statistic per gene
hist(tstats_obj$statistic,col=2)

# F statistic calc for lane number
fstats_obj = rowFtests(edata,as.factor(pdata$lane.number))
names(fstats_obj)

# hist of this
hist(fstats_obj$statistic,col=2, breaks = 200)

# fit many stats with limma
mod = model.matrix(~ pdata$strain)
fit_limma = lmFit(edata,mod)
ebayes_limma = eBayes(fit_limma)
head(ebayes_limma$t)

# plot 
plot(ebayes_limma$t[,2],-tstats_obj$statistic,col=4,
     xlab="Moderated T-stat",ylab="T-stat")
abline(c(0,1),col="darkgrey",lwd=3)

# fit statistics adjusted for lane number with limma 
mod_adj = model.matrix(~ pdata$strain + as.factor(pdata$lane.number))
fit_limma_adj = lmFit(edata,mod_adj)
ebayes_limma_adj = eBayes(fit_limma_adj)
head(ebayes_limma_adj$t)

# plot t stat of this
plot(ebayes_limma_adj$t[,2],-tstats_obj$statistic,col=3,
     xlab="Moderated T-stat",ylab="T-stat")
abline(c(0,1),lwd=3,col="darkgrey")

# compare null to alternative model (including lane effects)
mod_lane = model.matrix(~ as.factor(pdata$lane.number))
fit_limma_lane = lmFit(edata,mod_lane)
ebayes_limma_lane = eBayes(fit_limma_lane) 
head(ebayes_limma_lane$t)

# find differences in gene expression in lane 
top_lane = topTable(ebayes_limma_lane, coef=2:7,
                    number=dim(edata)[1],sort.by="none")
head(top_lane)

# Plot F-statistics
plot(top_lane$F,fstats_obj$statistic,
     xlab="Moderated F-statistic",ylab="F-statistic",col=3)

# Nested comparison in Edge
edge_study = build_study(edata, grp = as.factor(pdata$lane.number))
de_obj = lrt(edge_study)
qval = qvalueObj(de_obj)
plot(qval$stat,fstats_obj$statistic,col=4,
     xlab="F-stat from edge",ylab="F-stat from genefilter")

# adjust for vairables by adding adj.var 
edge_study2 = build_study(edata, grp = as.factor(pdata$lane.number),
                          adj.var=pdata$strain)
de_obj2 = lrt(edge_study2)
qval2 = qvalueObj(de_obj2)
plot(qval2$stat,fstats_obj$statistic,col=4,
     xlab="F-stat from edge",ylab="F-stat from genefilter")


