################################################################################
# File: TimeSeriesAnalysis.R                                                   #
# Purpose: Identify gene clusters that correspond to temporal patterns.        #
# Created: May 1, 2019                                                         #
# Author: Adam Faranda                                                         #
################################################################################

############################ Setup Environment #################################

# CHANGE THIS DIRECTORY
setwd('/home/adam/Documents/LEC_Time_Series')
#load("GeneLengthTable.Rdata")
library(dplyr)
library(cluster)
library(reshape2)
library(pheatmap)

wd<-getwd()
source('transcriptomic_analysis_scripts/BuildDataMatrix.R')
source('transcriptomic_analysis_scripts/PreprocessingFunctions.R')
source('transcriptomic_analysis_scripts/PrincipalComponents.R')
source('transcriptomic_analysis_scripts/ClusteringFunctions.R')

############################ Load in Data Files ################################

# CHANGE THIS DIRECTORY
dl<-"~/Documents/LEC_Time_Series_HTSeq_Counts"

# ft<-hc_getFileTable(
#   dirList=dl, filename = "HTSeq_GeneCounts_All.csv"
# )

ft<-hc_getFileTable(
  dirList=dl, filename = "HTSeq_GeneCounts_Wildtype.csv"
)

# ft<-hc_getFileTable(
#   dirList=dl, filename = "H TSeq_GeneCounts_Wildtype_0-48.csv"
# )

ds<-hc_loadFiles(ft)
ft<-hc_identifierConsistency(ds, ft)

library('AnnotationHub')
lt<-read.table(
  paste(dl,"gene_coding_lengths.txt", sep='/'),
  header = T, stringsAsFactors = F
)
ah<-AnnotationHub()
# Run Query to find proper Annotation Set:
#AnnotationHub::query(ah, pattern=c("EnsDb", "Mus musculus", "98"))
edb<-ah[['AH75036']]
lt<-merge(
  lt, AnnotationDbi::select(
    edb, keys=lt$gene_id, 
    columns = c("SYMBOL", "DESCRIPTION"), 
    keytype = "GENEID"
  ),
  by.x = 'gene_id', by.y='GENEID'
)
row.names(lt)<-lt$gene_id
rm(ah, edb)
detach(package:AnnotationHub, unload=T)
detach(package:ensembldb, unload=T)
detach(package:AnnotationFilter, unload=T)

# Setup an edgeR DGE List to keep everything organized
htseq_count<-hc_buildDataFrame(ds, ft)

dge<-DGEList(
  htseq_count, 
  samples=ft, 
  genes=lt[row.names(htseq_count),]
)

# Reorder grouping factor, drop unused levels
dge$samples$group<-droplevels(
  factor(
    paste(dge$samples$hours.pcs, 'H', sep=''),
    levels = c('0H', '6H', '24H','48H', '120H')
  )
)

# Build design matrix with batch coveriate
odg<-dge
dge<-dge[filterByExpr(dge), , keep.lib.sizes=FALSE]
dge$genes$Var<-apply(dge$counts, 1, var)
design<-model.matrix(~group + batch, dge$samples)

# Get statistical significance
dge<-calcNormFactors(dge)
dge<-estimateCommonDisp(dge)
dge<-estimateTagwiseDisp(dge)

# Calcuate Statistical Significance (Gene DE at ANY Timepoint)
fit<-glmFit(dge, design)
qlf<-glmLRT(fit, coef=2:5)
deg<-as.data.frame(topTags(qlf, n=Inf))


############# Generate Cpm Matrix for clustering ############################
 
mat<-as.data.frame(dge$counts)
ecpm<-edgeRcpm(mat)                     # Normalize using edgeR's TMM method



############## Apply LogFC Filters; Plot Principal Components ###############
varRanks <-c(10, 50, 100, 250, 500)           # Try different variance filters
for( v in varRanks){
	print(v)
	ecpm.filter <-varianceFilter(ecpm, threshold=v)
	
	# Plot results for TMM Normalized Data -- selected by variance
	f1<-paste('ECPM_Samples_Top_', v,'_Variance_Class.png')
	f2<-paste('ECPM_Samples_Top_', v,'_Variance_Lab.png')
	f3<-paste('ECPM_Samples_Top_', v,'_Variance_Clusters.png')
	png(f1, width=480, height=300)
		print(plotPrinComp(ecpm.filter, dge$samples, groupCol=1, idCol=0))
	dev.off()
	
	png(f2, width=480, height=300)
		print(plotPrinComp(ecpm.filter, dge$samples, groupCol=7, idCol=0))
	dev.off()
	
	png(f3, width=1200, height = 800)
	pheatmap(
	  log(ecpm.filter), 
	  annotation_col = dge$samples[,c('group', 'batch')],
	  show_rownames = F, fontsize=20, cellwidth = 20
	)
	dev.off()
}

for( lfc in c(2, 5, 7)){

	# Get list of genes meeting key criteria
	like<-as.character(
	  (
	    deg %>% 
	      filter(logCPM > 2 & FDR < 0.05) %>%
	      filter(
	        abs(logFC.group6H) > lfc |
	          abs(logFC.group24H) > lfc |
	          abs(logFC.group48H) > lfc |
	          abs(logFC.group120H) > lfc
	      )
	  )$gene_id
	)
	# Plot results selected by fold change level
	f1<-paste('ECPM_Samples_logFC_', lfc,'_',length(like),'_Significant_genes_Class.png', sep='')
	f2<-paste('ECPM_Samples_logFC_', lfc,'_',length(like),'_Significant_genes_Lab.png', sep='')
	f3<-paste('ECPM_Samples_logFC_', lfc,'_',length(like),'_Significant_genes_Clusters.png', sep='')
	f4<-paste('ECPM_Samples_logFC_', lfc,'_',length(like),'_Significant_genes_Matrix.txt', sep='')
	
	png(f1, width=480, height=300)
		print(plotPrinComp(ecpm[like,], dge$samples, groupCol=1, idCol=0))
	dev.off()
	png(f2, width=480, height=300)
		print(plotPrinComp(ecpm[like,], dge$samples, groupCol=7, idCol=0))
	dev.off()
	
	png(f3, width=1200, height = 800)
	  pheatmap(
	    log(ecpm[like,]), 
	    annotation_col = dge$samples[,c('group', 'batch')],
	    show_rownames = F, fontsize=20, cellwidth = 20
	  )
	  
	  sym<-dge$genes$SYMBOL
	  names(sym)<-dge$genes$gene_id
	  
	  mat<-as.data.frame(log(ecpm[like,]))
	  mat$gene<-sym[rownames(mat)]
	  mat<-mat[,c('gene', colnames(ecpm))]
	  write.table(
	    mat, f4, sep='\t', quote=F, row.names = F
	  )
	dev.off()
}

# Write Result Tables
htseq_count<-as.data.frame(htseq_count)
htseq_count$gene_id <-row.names(htseq_count)
write.table(
  htseq_count, "Injury_Model_Raw_Counts_All_Genes.txt", row.names=F, 
  sep = '\t', quote=F
)

htseq_fpkm<-as.data.frame(rpkm(odg))
write.table(
  htseq_fpkm, "Injury_Model_FPKM_All_Genes.txt", row.names=F, 
  sep = '\t', quote=F
)

# Filtered For Abundance: 
# at least 10 reads in three samples, sum accross all samples at least 15
# and at least 10 reads in 70% of samples in the smallest group
ecpm<-as.data.frame(ecpm)
ecpm$gene_id <-row.names(ecpm)
write.table(
  ecpm, "Injury_Model_TMM_Normalized_Present_Genes.txt", row.names=F, 
  sep = '\t', quote=F
)

write.table(
  dge$samples, "Sample_Metadata.txt", row.names=F, 
  sep = '\t', quote=F
)

write.table(
  dge$genes, "Gene_Metadata.txt", row.names=F, 
  sep = '\t', quote=F
)





