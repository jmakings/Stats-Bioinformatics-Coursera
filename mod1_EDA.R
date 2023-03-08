# Exploratory Data Analysis of genomic data in R

##### Lecture 1 #####
# https://www.coursera.org/learn/statistical-genomics/lecture/pOrnd/exploratory-analysis-in-r-part-i-7-22

# load pallete and update size of dots
tropical = c('darkorange','dodgerblue','hotpink','limegreen','yellow')
palette(tropical)
par(pch=19)

# Install and load required packages

# devtools::install_github('alyssafrazee/RSkittleBrewer')
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("org.Hs.eg.db")
# source("https://bioconductor.org/biocLite.R")
# biocLite("org.Hs.eg.db")
install.packages("tidyverse")
library(gplots)
library(devtools)
library(Biobase)
library(RSkittleBrewer)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(tidyverse)

# load data
con = url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)

# extract out phenotype, expression, and feature data with commands
bm = bodymap.eset
pdata=pData(bm)
edata=exprs(bm)
fdata = fData(bm)
ls()

# make tables
table(pdata$gender)
table(pdata$gender, pdata$race)
summary(edata)

# check NAs 
table(pdata$age,useNA="ifany")
table(is.na(pdata$age))
sum(pdata$age==" ")
sum(is.na(edata))

gene_na = rowSums(is.na(edata))
table(gene_na)

sample_na = rowSums(is.na(edata))
table(sample_na)

# make sure dimensions match up
dim(fdata)
dim(pdata)
dim(edata)

##### Lecture 2 #####
# https://www.coursera.org/learn/statistical-genomics/lecture/og2ns/exploratory-analysis-in-r-part-ii-10-07

# boxplot of log transformed expression data
# shows data is very skewed
boxplot(log2(edata+1),col=2,range=0)

# histogram for 2 samples
par(mfrow=c(1,2))
hist(log2(edata[,1]+1),col=2)
hist(log2(edata[,2]+1),col=2)

# density plot for 2 samples
plot(density(log2(edata[,1]+1)),col=2)
# this overlays a second plot on top of the first
lines(density(log2(edata[,2]+1)),col=3)

# qqplot of distributions of measurements between 2 samples
qqplot(log2(edata[,1]+1), log2(edata[,2]+1),col=3)
abline(c(0,1))

# M-A plot between two samples
mm = log2(edata[,1]+1) - log2(edata[,2]+1)
aa = log2(edata[,1]+1) + log2(edata[,2]+1)
plot(aa,mm,col=2)

# remove rows that are almost zero
# this doesnt work bc there is NA values that were not removed previosuly
edata = as.data.frame(edata)
filt_edata = filter(edata,rowMeans(edata) > 1)
dim(filt_edata)
dim(edata)
boxplot(as.matrix(log2(filt_edata+1)),col=2)

##### Lecture 3 #####
# https://www.coursera.org/learn/statistical-genomics/lecture/ccoiM/exploratory-analysis-in-r-part-iii-7-26

aeid = as.character(fdata[,1])
chr = AnnotationDbi::select(org.Hs.eg.db,keys=aeid,keytype="ENSEMBL",columns="CHR")
head(chr)

dim(chr)

dim(edata)

chr = chr[!duplicated(chr[,1]),]
all(chr[,1] == rownames(edata))

# Select the chromosome Y samples
edatay = dplyr::filter(edata,chr$CHR=="Y")
