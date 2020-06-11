################################################################################
# File: TimeSeriesAnalysis.R                                                   #
# Purpose: Identify gene clusters that correspond to temporal patterns.        #
# Created: May 1, 2019                                                         #
# Author: Adam Faranda                                                         #
################################################################################

############################ Setup Environment #################################
setwd('/home/adam/Documents/LEC_Time_Series')
#load("GeneLengthTable.Rdata")
library(dplyr)
library(cluster)
library(reshape2)
wd<-getwd()
source('transcriptomic_analysis_scripts/BuildDataMatrix.R')
source('transcriptomic_analysis_scripts/PreprocessingFunctions.R')
source('transcriptomic_analysis_scripts/PrincipalComponents.R')
source('transcriptomic_analysis_scripts/ClusteringFunctions.R')

############################ Load in Data Files ################################

dl<-"~/Documents/LEC_Time_Series_HTSeq_Counts"

# ft<-hc_getFileTable(
#   dirList=dl, filename = "HTSeq_GeneCounts_All.csv"
# )
# ft<-hc_getFileTable(
#   dirList=dl, filename = "HTSeq_GeneCounts_Wildtype.csv"
# )
ft<-hc_getFileTable(
  dirList=dl, filename = "HTSeq_GeneCounts_Wildtype.csv"
)

ds<-hc_loadFiles(ft)
ft<-hc_identifierConsistency(ds, ft)
lt<-read.table(
  paste(dl,"gene_coding_lengths.txt", sep='/'),
  header = T
)

# Setup an edgeR DGE List to keep everything organized
row.names(lt)<-lt$gene_id
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
design<-model.matrix(~group + batch, dge$samples)
dge<-dge[filterByExpr(dge, design=design),]
dge$genes$Var<-apply(dge$counts, 1, var)


# Get statistical significance

dge<-estimateCommonDisp(dge)
dge<-estimateTagwiseDisp(dge)
dge<-calcNormFactors(dge)

# Calcuate Statistical Significance (Gene DE at ANY Timepoint)
fit<-glmFit(dge, design)
qlf<-glmLRT(fit, coef=2:4)
deg<-as.data.frame(topTags(qlf, n=Inf))

# Get list of genes meeting key criteria
lfc<-7
like<-as.character(
  (
    deg %>% 
      filter(logCPM > 2 & FDR < 0.05) %>%
      filter(
        abs(logFC.group6H) > lfc |
          abs(logFC.group24H) > lfc |
          abs(logFC.group48H) > lfc 
      )
  )$gene_id
)

######################### Apply Preprocessing Steps ############################
 
mat<-as.data.frame(dge$counts)

ecpm<-edgeRcpm(mat)                     # Normalize using edgeR's TMM method
etpm<-normGeneLength(                   # Normalize to Transcripts per million
  ft=ft, mat=mat, gt=dge$genes, 
  gt_ln='coding_length'
)

ecmb<-ComBat(                                  # Apply Combat Batch Correction
  ecpm, mod=model.matrix(~group, dge$samples),  # to TMM normalized CPM values
  batch=as.factor(dge$samples$batch)
)

ecmb<-fixCombatNegatives(ecmb, idCol=0)             # replace negative values with min +ve
etpm[etpm==0]<-0.0001
etmb<-ComBat(
  dat=as.matrix(etpm), 
  batch=ft$batch,
  mod = model.matrix(~1, ft)
)
etmb<-fixCombatNegatives(etmb, idCol =0)

############## Apply Variance Filters; Plot Principal Components ###############
varRanks <-c(10, 50, 100, 10000)                # Try different variance filters
for( v in varRanks){
	print(v)
	ecpm.filter <-varianceFilter(ecpm, threshold=v)
	ecmb.filter <-varianceFilter(ecmb, threshold=v)
	
	# Plot results for TMM Normalized Data
	f1<-paste('ECPM_Samples_Top_', v,'_Ranked_Class.png')
	f2<-paste('ECPM_Samples_Top_', v,'_Ranked_Lab.png')
	png(f1, width=240, height=150)
		print(plotPrinComp(ecpm.filter, dge$samples, groupCol=1, idCol=0))
	dev.off()
	
	png(f2, width=240, height=150)
		print(plotPrinComp(ecpm.filter, dge$samples, groupCol=7, idCol=0))
	dev.off()

	# Plot results for batch adjusted TMM data
	f1<-paste('ECMB_Samples_Top_', v,'_Ranked_Class.png')
	f2<-paste('ECMB_Samples_Top_', v,'_Ranked_Lab.png')
	png(f1, width=240, height=150)
		print(plotPrinComp(ecmb.filter, dge$samples, groupCol=1, idCol=0))
	dev.off()
	png(f2, width=240, height=150)
		print(plotPrinComp(ecmb.filter, dge$samples, groupCol=7, idCol=0))
	dev.off()
}

############ Analyze Sample Clusters at desired Variance Threshold #############
distm <-c('euclidean', 'manhattan')			  # Try different distance methods
linkm <-c('complete', 'average', 'single')    #  Try different linkage methods
trees <-c(1,2,3,4,5,6)                        # Different levels k
v =200
ecpm.filter<-varianceFilter(ecpm, threshold=v)
ecmb.filter<-varianceFilter(ecmb, threshold=v)

clustStats<-rbind(
	summarizeSampleClusters(
		data=ecpm, distm=distm, linkm=linkm, v=v, label='ecpm'
		
	),
	summarizeSampleClusters(
		data=ecmb, distm=distm, linkm=linkm, v=v, label='ecmb'
	)
)

write.csv(clustStats, 'Sample_Cluster_Statistics.csv')

################# Analyze Gene Clusters at Variance Threshold ####################
v = 200
mat<-varianceFilter(etmb, threshold=v)
mat.l<-log(mat)
mat.s<-scaleCenterByRow(mat)

d.meth='manhattan'
h.method='complete'
h<-wrapHclust(mat.s, idCol=0, transpose = F)
ktable<-tabulate_H_Clusters(h)
ft<-dge$samples
ft$Sample<-row.names(dge$samples)
ft$Class<-ft$group
ft$Lab<-ft$batch
rt<-reshapeClusterTable(mat.s, ktable, 5, ft=ft)

rt<-inner_join(rt, ft[c('Sample','group', 'batch')], by='Sample' )
mt<- rt %>% group_by(ID, group) %>% summarize(nth(Cluster,1), mean(value))

summarizeGeneClusters(mat.l, label='Scaled_Centered_v200', nclust=5)


v = 2000
ecmb.filter<-varianceFilter(ecmb, threshold=v)
mat<-as.matrix(ecmb.filter[,2:ncol(ecmb.filter)])
mat.l<-log(mat)
mat.s<-scaleCenterByRow(mat)
summarizeGeneClusters(mat.l, label='Scaled_Centered_v2000')
  
design<-model.matrix(~Interval+Lab, ft)
rownames(design)<-colnames(y)

# Remove genes with fewer than 10 counts in 'n' samples based on design matrix
z<-y[filterByExpr(y, design=design), , keep.lib.sizes=FALSE]

# Calculate Parameter Estimates / Normalization Factors
#z<-estimateDisp(z, design=design)
dge<-estimateCommonDisp(dge)
dge<-estimateTagwiseDisp(dge)
dge<-calcNormFactors(dge)
  





