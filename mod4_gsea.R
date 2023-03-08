# Gene Set Enrichment 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("goseq")

library(devtools)
library(Biobase)
library(goseq)
library(DESeq2)

# different databases on goseq 
head(supportedGenomes())

# different supported gene IDs 
head(supportedGeneIDs())

# goseq analysis

# load example from goseq package
temp_data =read.table(system.file("extdata","Li_sum.txt",
                                  package="goseq"),sep="\t",
                      header=TRUE,
                      stringsAsFactors=FALSE)
expr= temp_data[,-1]
rownames(expr) = temp_data[,1]
expr = expr[rowMeans(expr) > 5,]
grp=factor(rep(c("Control","Treated"),times=c(4,3)))
pdata  = data.frame(grp)

# perform differential expression with DESeq2
de = DESeqDataSetFromMatrix(expr, pdata, ~grp)
de_fit = DESeq(de)
de_results = results(de_fit)

# differentially expressed results
genes = as.integer(de_results$padj < 0.05)
not_na = !is.na(genes)
names(genes) = rownames(expr)
genes = genes[not_na]

# Calculate probability weight function
pwf=nullp(genes,"hg19","ensGene")

head(pwf)

# parametric test to look for enrichment differences
# with respect to different categories
GO.wall=goseq(pwf,"hg19","ensGene")

# shows over and under-represented p-values, as well 
# as differentially expressed genes for each gene set
head(GO.wall)

# limiting to a singular category
GO.MF=goseq(pwf,"hg19","ensGene",test.cats=c("GO:MF"))
head(GO.MF)

