##### Dimension Reduction for Genomics ######

library(devtools)
library(Biobase)

# load data
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)
ls()

# reduce size of data set, log transform, and center data
edata = edata[rowMeans(edata) > 100, ]
edata = log2(edata + 1)
edata_centered = edata - rowMeans(edata)

# apply svd (3 matrices)
svd1 = svd(edata_centered)
names(svd1)

# look at d matrix to see variance explained by each singular value
plot(svd1$d,ylab="Singular value",col=2)

# use formula to convert this to percentage of variance explained
plot(svd1$d^2/sum(svd1$d^2),ylab="Percent Variance Explained",col=2)

# (first value explains over 50% of variance)

# plot top two SVs 
par(mfrow=c(1,2))
plot(svd1$v[,1],col=2,ylab="1st SV")
plot(svd1$v[,2],col=2,ylab="2nd SV")

# Plot SV1 vs SV2
plot(svd1$v[,1],svd1$v[,2],col=2,ylab="2nd SV",xlab="1st SV")

# add color based on which study data is from
plot(svd1$v[,1],svd1$v[,2],ylab="2nd SV",
     xlab="1st SV",col=as.numeric(pdata$study))

# apply boxplots 
boxplot(svd1$v[,1] ~ pdata$study,border=c(1,2))
points(svd1$v[,1] ~ jitter(as.numeric(pdata$study)),col=as.numeric(pdata$study))

# PCs vs SVs
pc1 = prcomp(edata)
plot(pc1$rotation[,1],svd1$v[,1])

# To get the actual PCs you have to subtract the column means rather 
# than the row means when normalizing.
edata_centered2 = t(t(edata) - colMeans(edata))
svd2 = svd(edata_centered2)

# here to show that svd2 is now equivalent to the PCs
plot(pc1$rotation[,1],svd2$v[,1],col=2)

# lets introduce a single outlier to see what happens
edata_outlier = edata_centered
edata_outlier[1,] = edata_centered[1,] * 10000
svd3 = svd(edata_outlier)
par(mfrow=c(1,2))
# this single outlier has changed our singular values so they are no longer equivalent
plot(svd1$v[,1],col=1,main="Without outlier")
plot(svd3$v[,1],col=2,main="With outlier")
plot(svd1$v[,1], svd3$v[,1], xlab="Without outlier", ylab="With outlier")

# The top singular vector is perfectly correlated with the outlying gene
plot(svd3$v[,1],edata_outlier[1,],col=4)

##### Quantile Normalization ######

#BiocManager::install("preprocessCore")
library(devtools)
library(Biobase)
library(preprocessCore)

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)
ls()

# log transform data and remove lower values
edata = log2(edata + 1)
edata = edata[rowMeans(edata) > 3, ]

# density plot for first 20 samples
colramp = colorRampPalette(c(3,"white",2))(20)
plot(density(edata[,1]),col=colramp[1],lwd=3,ylim=c(0,.30))
for(i in 2:20){lines(density(edata[,i]),lwd=3,col=colramp[i])}

# quantile normalization 
norm_edata = normalize.quantiles(as.matrix(edata))
plot(density(norm_edata[,1]),col=colramp[1],lwd=3,ylim=c(0,.20))
for(i in 2:20){lines(density(norm_edata[,i]),lwd=3,col=colramp[i])}

# check distributions
plot(norm_edata[1,],col=as.numeric(pdata$study))

# plot principal components (still differences between studies)
svd1 = svd(norm_edata - rowMeans(norm_edata))
plot(svd1$v[,1],svd1$v[,2],xlab="PC1",ylab="PC2",
     col=as.numeric(pdata$study))

# takeaway from this is even after quantile normalization,
# there can still be artifacts (batch effects, experiemnte differences)





