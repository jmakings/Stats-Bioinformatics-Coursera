if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("Biobase","sva","bladderbatch","snpStats"))

library(devtools)
library(Biobase)
library(sva)
library(bladderbatch)
library(snpStats)

# get data
data(bladderdata)

pheno = pData(bladderEset)
edata = exprs(bladderEset)

# use linear model to adjust for batch effects
mod = model.matrix(~as.factor(cancer) + as.factor(batch),data=pheno)
fit = lm.fit(mod,t(edata))
hist(fit$coefficients[2,],col=2,breaks=100)

table(pheno$cancer,pheno$batch)

# clean batch effects using Combat
batch = pheno$batch
modcombat = model.matrix(~1, data=pheno)
modcancer = model.matrix(~cancer, data=pheno)
combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=TRUE)

# do linear model using cleaned data 
combat_fit = lm.fit(modcancer,t(combat_edata))
hist(combat_fit$coefficients[2,],col=2,breaks=100)

# compare combat and linear adjustment
plot(fit$coefficients[2,],combat_fit$coefficients[2,],col=2,
     xlab="Linear Model",ylab="Combat",xlim=c(-5,5),ylim=c(-5,5))
abline(c(0,1),col=1,lwd=3)

# adjusting for batch effects using SVA
mod = model.matrix(~cancer,data=pheno)
mod0 = model.matrix(~1, data=pheno)
sva1 = sva(edata,mod,mod0,n.sv=2)

summary(lm(sva1$sv ~ pheno$batch))

# box plot of surrogate variables
boxplot(sva1$sv[,2] ~ pheno$batch)
points(sva1$sv[,2] ~ jitter(as.numeric(pheno$batch)),col=as.numeric(pheno$batch))

# add surrogate variables to the model 
modsv = cbind(mod,sva1$sv)
fitsv = lm.fit(modsv,t(edata))

# compare combat vs linear model 
par(mfrow=c(1,2))
plot(fitsv$coefficients[2,],combat_fit$coefficients[2,],col=2,
     xlab="SVA",ylab="Combat",xlim=c(-5,5),ylim=c(-5,5))
abline(c(0,1),col=1,lwd=3)
plot(fitsv$coefficients[2,], fit$coefficients[2,],col=2,
     xlab="SVA",ylab="linear model",xlim=c(-5,5),ylim=c(-5,5))
abline(c(0,1),col=1,lwd=3)

# principal components for population structure (Part 2)

# load data
data(for.exercise)
controls <- rownames(subject.support)[subject.support$cc==0]
use <- seq(1, ncol(snps.10), 10)
ctl.10 <- snps.10[controls,use]

# calculate the PCs 
xxmat <- xxt(ctl.10, correct.for.missing=FALSE)
evv <- eigen(xxmat, symmetric=TRUE)
pcs <- evv$vectors[,1:5]

# compare PCs to population labels
pop <- subject.support[controls,"stratum"]
plot(pcs[,1],pcs[,2],col=as.numeric(pop),
     xlab="PC1",ylab="PC2")
legend(0,0.15,legend=levels(pop),pch=19,col=1:2)


















