# Purpose: Split columns from GeneCount.ReCalcRPKM.tsv into individual samples
# Created: March 5, 2019
# Author: Adam Faranda

# Setup Environment
options(echo=F)
library(dplyr)
library(openxlsx)

wd<-getwd()
outDir<-paste(wd,'GEO_Split_Files', sep="/")
dataFile<-"GeneCount.ReCalcRPKM.tsv"

# Import data tables
setwd(wd)
countData<-read.table(dataFile, quote="", sep="\t", header=T)

commonCols<-c("Chr", "GeneID", "Start", "Stop", "CodingLength")
setwd(outDir)
for ( n in names(countData)[grep("_GeneCount", names(countData))]){
	m<-(gsub("_GeneCount","_RPKM",n))
	print(paste(gsub("_RPKM", ".tab.txt",m), sep=""))
	print(head(countData[c(commonCols, n, m)]))
	write.table(
		countData[c(commonCols, n, m)],
		file=paste(gsub("_RPKM", ".tab.txt",m), sep=""),
		quote=F, sep='\t', row.names=F
	)
}