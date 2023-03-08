install.packages(c("devtools","broom"))

library(devtools)
library(Biobase)
library(broom)

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
pdata=pData(bm)
edata=as.data.frame(exprs(bm))
fdata = fData(bm)
ls()

# simple linear regression 
edata = as.matrix(edata)
lm1 = lm(edata[1,] ~ pdata$age)
tidy(lm1)

# plot linear regression 
plot(pdata$age,edata[1,], col=1)
abline(lm1$coeff[1],lm1$coeff[2], col=2,lwd=3)

# doing this for categorical variables
pdata$gender
table(pdata$gender)

# visualize difference between genders 
boxplot(edata[1,] ~ pdata$gender)
points(edata[1,] ~ jitter(as.numeric(pdata$gender)),
       col=as.numeric(pdata$gender))

# dummy variables to quantify differences between gender
dummy_m = pdata$gender=="M"
dummy_m

dummy_f = pdata$gender=="F"
dummy_f

# but R can do this for you
lm2 = lm(edata[1,] ~ pdata$gender)
tidy(lm2)

# multiple coefficients for one variable
tidy(lm(edata[1,] ~ pdata$tissue.type ))

# expanding the model to include age and gender
lm3 = lm(edata[1,] ~ pdata$age + pdata$gender)
tidy(lm3)

# fitting interaction term (age interacting with gender)
lm4 = lm(edata[1,] ~ pdata$age*pdata$gender)
tidy(lm4)

# viewing outliers

# in this plot, there is a significant outlier, but it has small impact
lm4 = lm(edata[6,] ~ pdata$age)
plot(pdata$age,edata[6,],col=2)
abline(lm4,col=1,lwd=3)

# but in this plot, it has a huge impact
index = 1:19
lm5 = lm(edata[6,] ~ index)
plot(index,edata[6,],col=2)
abline(lm5,col=1,lwd=3)

lm6 = lm(edata[6,-19] ~ index[-19])
abline(lm6,col=3,lwd=3)

legend(5,1000,c("With outlier","Without outlier"),col=c(1,3),lwd=3)

# In general, we want residuals to be normally distributed, but they arent

# this to change layout of plots
par(mfrow=c(2,2))

# not good looking residuals
hist(lm6$residuals,col=2)
hist(lm5$residuals,col=3)

# transform data before regression

gene1 = log2(edata[1,]+1)
lm7 = lm(gene1 ~ index)
# now our residuals look closer to a normal distribution
hist(lm7$residuals,col=4)

# too many covariates in this model so R will fail at this. 
# can also happen when two variables are highly correlated
lm8 = lm(gene1 ~ pdata$tissue.type + pdata$age)
tidy(lm8)

# fitting 16 data points to 18 parameters
dim(model.matrix( ~ pdata$tissue.type + pdata$age))

# plot the residuals by 
par(mfrow=c(1,1))
colramp = colorRampPalette(1:4)(17)
lm9 = lm(edata[2,] ~ pdata$age)
plot(lm9$residuals,col=colramp[as.numeric(pdata$tissue.type)])

# many regression models

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edge")

library(devtools)
library(Biobase)
library(limma)
library(edge)

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
pdata=pData(bot)
edata=as.matrix(exprs(bot))
fdata = fData(bot)
ls()

# transform the data
edata = log2(as.matrix(edata) + 1)
edata = edata[rowMeans(edata) > 10, ]
dim(t(edata))

# build many regression models at once
mod = model.matrix(~ pdata$strain)
fit = lm.fit(mod,t(edata))
names(fit)

# view coefficients from multi vs single regression
fit$coefficients[,1]

tidy(lm(as.numeric(edata[1, ]) ~ pdata$strain))
# as you can see, the intercept and coefficient is the same for both

# view coefficient histograms
par(mfrow=c(1,2))
hist(fit$coefficients[1,],breaks=100,col=2,xlab="Intercept")
hist(fit$coefficients[2,],breaks=100,col=2,xlab="Strain")
abline(v=0,lwd=3,col=3)

# view residuals
par(mfrow=c(1,2))
plot(fit$residuals[,1],col=2)
plot(fit$residuals[,2],col=2)

# fit with adjustment for lane number
mod_adj = model.matrix(~ pdata$strain + as.factor(pdata$lane.number))
fit_adj = lm.fit(mod_adj,t(edata))
fit_adj$coefficients[,1]

# using limma package
fit_limma = lmFit(edata,mod_adj)
names(fit_limma)

fit_limma$coefficients[1,]

fit_adj$coefficients[,1]

# edge package to fit many regressions
edge_study = build_study(data=edata,grp=pdata$strain,adj.var=as.factor(pdata$lane.number))
fit_edge = fit_models(edge_study)
summary(fit_edge)

fit_edge@beta.coef[1,]

fit_limma$coefficients[1,]

