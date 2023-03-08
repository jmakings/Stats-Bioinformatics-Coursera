
BiocManager::install('MatrixEQTL')

library(devtools)
library(Biobase)
library(MatrixEQTL)

# load data
base.dir = find.package("MatrixEQTL")
SNP_file_name = paste(base.dir, "/data/SNP.txt", sep="");
expression_file_name = paste(base.dir, "/data/GE.txt", sep="")
covariates_file_name = paste(base.dir, "/data/Covariates.txt", sep="")
output_file_name = tempfile()

# create expr variable
expr = read.table(expression_file_name,sep="\t",
                  header=T,row.names=1)
expr[1,]

# snps variable read in
snps = read.table(SNP_file_name,sep="\t",
                  header=T,row.names=1)
snps[1,]

# covariates variable read in 
cvrt = read.table(covariates_file_name,sep="\t",
                  header=T,row.names=1)

# eQTL is a simple linear regression for SNP/gene pairs
e1 = as.numeric(expr[1,])
s1 = as.numeric(snps[1,])
lm1 = lm(e1 ~ s1)
tidy(lm1)

# plot expression level vs genotype
par(mfrow=c(1,1))
plot(e1 ~ jitter(s1),
     col=(s1+1),xaxt="n",xlab="Genotype",ylab="Expression")
axis(1,at=c(0:2),labels=c("AA","Aa","aa"))
lines(lm1$fitted ~ s1,type="b",pch=15,col="darkgrey")

# fitting many eQTL models

# set up variables beforehand
pvOutputThreshold = 1e-2
errorCovariance = numeric()
useModel = modelLINEAR

# set up data in special format for snps
snps = SlicedData$new()
snps$fileDelimiter = "\t"     # the TAB character
snps$fileOmitCharacters = "NA" # denote missing values;
snps$fileSkipRows = 1          # one row of column labels
snps$fileSkipColumns = 1       # one column of row labels
snps$fileSliceSize = 2000     # read file in pieces of 2,000 rows
snps$LoadFile( SNP_file_name )

# set up data in special format for genes
gene = SlicedData$new()
gene$fileDelimiter = "\t"      # the TAB character
gene$fileOmitCharacters = "NA" # denote missing values;
gene$fileSkipRows = 1          # one row of column labels
gene$fileSkipColumns = 1      # one column of row labels
gene$fileSliceSize = 2000      # read file in pieces of 2,000 rows
gene$LoadFile(expression_file_name)

cvrt = SlicedData$new()

# main command for running matrix_eQTL
me = Matrix_eQTL_engine(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = NULL,
  pvOutputThreshold = pvOutputThreshold,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)

# plot it 
plot(me)

# extra info
names(me$all)
me$all$neqtls
me$all$eqtls




