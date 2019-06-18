################################################################################
# File: TimeSeriesAnalysis.R												       #
# Purpose: Identify gene clusters that correspond to temporal patterns.		   #
# Created: May 1, 2019 														   #
# Author: Adam Faranda														   #
################################################################################

############################ Setup Environment #################################
setwd('/Users/afaranda/Documents/CISC683_Project_Report')
library(dplyr)
library(cluster)
library(reshape2)
# wd<-getwd()
source('BuildDataMatrix.R')
source('PreprocessingFunctions.R')
source('PrincipalComponents.R')
source('ClusteringFunctions.R')
######################### Apply Preprocessing Steps ############################

# Add Class attribute to feature definition table
ft<-raw[[1]]
ft$Class<-as.factor(paste("Hour",ft$Hours_PCS,sep=""))
ft$Class<-factor(
	ft$Class, 
	levels=unique(
		as.character(ft$Class[order(ft$Hours_PCS)])
	)
)


counts<-(raw[[2]][6:nrow(raw[[2]]),])      # Extract Raw Counts
ecpm<-edgeRcpm(counts)                     # Normalize using edgeR's TMM method
ecmb<-wrapCombat_intOnly(ecpm, ft)         # Correct for batch effects
ecmb<-fixCombatNegatives(ecmb)             # replace negative values with min +ve


############## Apply Variance Filters; Plot Principal Components ###############
varRanks <-c(10, 50, 100, 200)                # Try different variance filters
for( v in varRanks){
	print(v)
	ecpm.filter <-varianceFilter(ecpm, threshold=v)
	ecmb.filter <-varianceFilter(ecmb, threshold=v)
	
	# Plot results for TMM Normalized Data
	f1<-paste('ECPM_Samples_Top_', v,'_Ranked_Class.png')
	f2<-paste('ECPM_Samples_Top_', v,'_Ranked_Lab.png')
	png(f1, width=240, height=150)
		print(plotPrinComp(ecpm.filter, ft, groupCol=8, idCol=1))
	dev.off()
	png(f2, width=240, height=150)
		print(plotPrinComp(ecpm.filter, ft, groupCol=4, idCol=1))
	dev.off()

	# Plot results for batch adjusted TMM data
	f1<-paste('ECMB_Samples_Top_', v,'_Ranked_Class.png')
	f2<-paste('ECMB_Samples_Top_', v,'_Ranked_Lab.png')
	png(f1, width=240, height=150)
		print(plotPrinComp(ecmb.filter, ft, groupCol=8, idCol=1))
	dev.off()
	png(f2, width=240, height=150)
		print(plotPrinComp(ecmb.filter, ft, groupCol=4, idCol=1))
	dev.off()
}

############ Analyze Sample Clusters at desired Variance Threshold #############
distm <-c('euclidean', 'manhattan')			  # Try different distance methods
linkm <-c('complete', 'average', 'single')    # Try different linkage methods
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
ecmb.filter<-varianceFilter(ecmb, threshold=v)
mat<-as.matrix(ecmb.filter[,2:ncol(ecmb.filter)])
mat.l<-log(mat)
mat.s<-scaleCenterByRow(mat)

d.meth='manhattan'
h.method='complete'

summarizeGeneClusters(mat.s, label='Scaled_Centered_v200')


v = 2000
ecmb.filter<-varianceFilter(ecmb, threshold=v)
mat<-as.matrix(ecmb.filter[,2:ncol(ecmb.filter)])
mat.l<-log(mat)
mat.s<-scaleCenterByRow(mat)
summarizeGeneClusters(mat.s, label='Scaled_Centered_v200')




