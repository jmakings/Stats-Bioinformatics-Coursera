# Data Transforms

library(devtools)
library(Biobase)

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
pdata=pData(bm)
edata=as.data.frame(exprs(bm))
fdata = fData(bm)
ls()

hist(rnorm(1000),col=2)

# normal data
hist(edata[,1],col=2,breaks=100)

# log scaled data
hist(log(edata[,1]),col=2,breaks=100)

# however 
min(log(edata))
quantile(log(edata[,1]))

# therefore we must add 1 to not have negative infinity
min(log(edata[,1] + 1)) 
quantile(log(edata[,1]+1))

# log transformed
hist(log(edata[,1] + 1),breaks=100,col=2)

# log2 transformed
hist(log2(edata[,1] + 1),breaks=100,col=2)

# look at values greater than 0
hist(log2(edata[,1] + 1),breaks=100,col=2,xlim=c(1,15),ylim=c(0,400))

# how many rows sums have expression data equal to 0
hist(rowSums(edata==0),col=2)

# filter out genes with low expression across all samples
low_genes = rowMeans(edata) < 5
table(low_genes)

filt_edata = filter(as.data.frame(edata), !low_genes)
dim(filt_edata)
summary(edata)

# compares means and medians
low_genes2 = rowMedians(as.matrix(edata)) < 5
table(low_genes2,low_genes)

filt_edata2 = filter(edata,!low_genes2)
dim(filt_edata2)

hist((filt_edata[,1] + 1),col=2,xlim=c(0,20),ylim=c(0,1500))

     