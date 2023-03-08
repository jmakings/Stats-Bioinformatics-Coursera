
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





