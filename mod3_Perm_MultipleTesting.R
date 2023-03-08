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

##### Multiple Testing 1 ######
library(qvalue)

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
pdata=pData(bot)
edata=as.matrix(exprs(bot))
fdata = fData(bot)

edata = log2(as.matrix(edata) + 1)
edata = edata[rowMeans(edata) > 10, ]

fstats_obj = rowFtests(edata,as.factor(pdata$strain))
hist(fstats_obj$p.value,col=2)

edge_study = build_study(edata, grp = pdata$strain, 
                         adj.var = as.factor(pdata$lane.number))
de_obj = lrt(edge_study)
qval = qvalueObj(de_obj)
hist(qval$pvalues,col=3)

mod = model.matrix(~ pdata$strain + pdata$lane.number)
fit_limma = lmFit(edata,mod)
ebayes_limma = eBayes(fit_limma)
limma_pvals = topTable(ebayes_limma,number=dim(edata)[1])$P.Value
hist(limma_pvals,col=4)

# calculating empirical permutation p-values with edge
set.seed(3333)
B = 1000
tstats_obj = rowttests(edata,pdata$strain)
tstat0 = matrix(NA,nrow=dim(edata)[1],ncol=B)
tstat = tstats_obj$statistic
strain = pdata$strain
for(i in 1:B){
  strain0 = sample(strain)
  tstat0[,i] = rowttests(edata,strain0)$statistic
}

emp_pvals = empPvals(tstat,tstat0)
hist(emp_pvals,col=2)

##### Multiple Testing 2 #####

# multiple testing corrected p-value 
fp_bonf = p.adjust(fstats_obj$p.value,method="bonferroni")
hist(fp_bonf,col=3)

# bonferroni corrected quantiles
quantile(fp_bonf)

# Benjamini-Hochberg
fp_bh = p.adjust(fstats_obj$p.value,method="BH")
hist(fp_bh,col=3)

# nothing significant lol
sum(fp_bh < 0.05)

# adjusted p-values from limma
limma_pvals_adj = topTable(ebayes_limma,number=dim(edata)[1])$adj.P.Val
hist(limma_pvals_adj,col=2)

# only 2 found
sum(limma_pvals_adj < 0.05)

quantile(limma_pvals_adj)

# direct q-values 
qval_limma = qvalue(limma_pvals)
summary(qval_limma)

# q value using edge
qval = qvalueObj(de_obj)
summary(qval)




