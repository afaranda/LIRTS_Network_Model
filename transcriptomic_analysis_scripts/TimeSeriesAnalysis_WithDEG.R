################################################################################
# File: TimeSeriesAnalysis.R												       #
# Purpose: Identify gene clusters that correspond to temporal patterns.		   #
# Created: May 1, 2019 														   #
# Author: Adam Faranda														   #
################################################################################

############################ Setup Environment #################################
setwd('/Users/afaranda/Documents/CISC683_Project')
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
ecpm.log<-log(ecmb[,2:19])				   # Apply log transformation Combat

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
v = 500
ecpm.filter<-varianceFilter(ecpm, threshold=v)
ecmb.filter<-varianceFilter(ecmb, threshold=v)

clustStats<-rbind(
	summarizeSampleClusters(
		data=ecpm, distm=distm, linkm=linkm, v=50, label='ecpm'
		
	),
	summarizeSampleClusters(
		data=ecmb, distm=distm, linkm=linkm, v=50, label='ecmb'
	)
)

write.csv(clustStats, 'Sample_Cluster_Statistics.csv')

################# Analyze Gene Clusters Based on Variance ######################
v = 500
ecmb.filter<-varianceFilter(ecmb, threshold=v)
row.names(ecmb.filter)<-ecmb.filter$ID
mat<-ecmb.filter[,2:ncol(ecpm.filter)]

h<-wrapHclust(log(mat), idCol=0, transpose=F, d.meth='euclidean', h.method='complete')
ktable<-tabulate_H_Clusters(h, ks=1:20)
ct<-reshapeClusterTable(log(mat), ktable, ft, k=20)

cl<-ct %>% 
	group_by(Cluster) %>%
	summarize(Count = n()/18, Mean_logCPM=mean(value), STDev_logCPM = sd(value)) %>%
	filter(Count > 5, STDev_logCPM < 1)
	
cl<-cl$Cluster

bp<-ggplot(data=ct[ct$Cluster %in% cl,], mapping=aes(x=Class, y=value, color=Class)) + 
	geom_boxplot() + facet_grid( . ~Cluster)
bp

########### Analyze Gene Clusters Based Differential Expression ################
# Get differentially expressed genes and use that to guide analysis
dbi.contrasts <-list(c("Hour0", "Hour24"), c("Hour0", "Hour48"), c("Hour24", "Hour48"))
dna.contrasts <-list(c("Hour0", "Hour6"), c("Hour0", "Hour24"), c("Hour6", "Hour24"))
for( i in unique(ft$Seq_Lab)){
	if(i == 'DNA'){
		for(j in dna.contrasts){
			g1<-ft[ft$Seq_Lab == i & ft$Class == j[1],'Sample_Number']
			g2<-ft[ft$Seq_Lab == i & ft$Class == j[2],'Sample_Number']
			deg<-edgeRPairwise(raw[[2]][6:nrow(raw[[2]]),], ft, group.dn=g1, group.up=g2)[[2]]
			deg$Lab <- i
			deg$DownGroup = j[1]
			deg$UpGroup = j[2]
			if(!exists('deg_master')){ deg_master<-deg}
			else{deg_master<-rbind(deg_master, deg)}
		}
	}
	else{
		for(j in dbi.contrasts){
			g1<-ft[ft$Seq_Lab == i & ft$Class == j[1],'Sample_Number']
			g2<-ft[ft$Seq_Lab == i & ft$Class == j[2],'Sample_Number']
			deg<-edgeRPairwise(raw[[2]][6:nrow(raw[[2]]),], ft, group.dn=g1, group.up=g2)[[2]]
			deg$Lab <- i
			deg$DownGroup = j[1]
			deg$UpGroup = j[2]
			if(!exists('deg_master')){deg_master<-deg}
			else{deg_master<-rbind(deg_master, deg)}
		
		}
	}
}


# Get DBI Degs -- Try Clustering Biclustering -- not much better
deg_master %>%
	filter(FDR < 0.05, abs(logFC) > 1) %>%
	group_by(Lab, DownGroup, UpGroup) %>%
	summarize(n(), Upregulated=sum(logFC > 1), Downregulated=sum(logFC < 1))

# Get Genes that are significant in DBI based on pairwise contrasts
dbi.gl<-inner_join(
	deg_master %>% 
	filter(FDR < 0.05, abs(logFC) > 1) %>%
	filter(DownGroup == 'Hour0', UpGroup == 'Hour24', Lab=='DBI') %>%
	dplyr::select(genes, logFC) %>%
	rename( H0vsH24_logFC= 'logFC'),
	
	deg_master %>% 
	filter(FDR < 0.05, abs(logFC) > 1) %>%
	filter(DownGroup == 'Hour24', UpGroup == 'Hour48', Lab=='DBI') %>%
	dplyr::select(genes, logFC) %>%
	rename(H24vsH48_logFC = 'logFC'),
	by='genes'
)

dna.gl<-inner_join(
	deg_master %>% 
	filter(FDR < 0.05, abs(logFC) > 1) %>%
	filter(DownGroup == 'Hour0', UpGroup == 'Hour6', Lab=='DNA') %>%
	select(genes, logFC) %>%
	rename( H0vsH24_logFC= 'logFC'),
	
	deg_master %>% 
	filter(FDR < 0.05, abs(logFC) > 1) %>%
	filter(DownGroup == 'Hour6', UpGroup == 'Hour24', Lab=='DNA') %>%
	select(genes, logFC) %>%
	rename(H24vsH48_logFC = 'logFC'),
	by='genes'
)

##############################################################
# Note -- Log transformation DRAMATICALLY improves clustering#
##############################################################
ecpm.filter<-ecpm[ecpm$ID %in% dbi.gl$genes, 1:10]
row.names(ecpm.filter)<-ecpm.filter$ID
ecpm.filter<-ecpm.filter[,2:10]
mat<-as.matrix(ecpm.filter)

h<-wrapHclust(log(mat), d.meth = 'maximum', h.method='single', idCol=0, transpose=F)
ktable<-tabulate_H_Clusters(h, ks=1:50)
plotGeneCluster(mat, ktable, k=20, c=1, ft, groupCol=8)
agnes()

# Get Clusters for DNALink Genes - - - - - 
ecpm.filter<-ecpm[ecpm$ID %in% dna.gl$genes, c(1,11:19)]
row.names(ecpm.filter)<-ecpm.filter$ID
ecpm.filter<-ecpm.filter[,2:10]
mat<-as.matrix(ecpm.filter)


ktable<-tabulate_k_means(log(mat), idCol=0, transpose=F, ks=1:20)
h<-wrapHclust(mat, d.meth = 'manhattan', h.method='complete', idCol=0, transpose=F)
ktable<-tabulate_H_Clusters(h, ks=1:50)
plotGeneCluster(mat, ktable, k=20, c=1, ft, groupCol=8)

# Test some Plotting Functions
l<-plotGeneCluster(mat, ktable, k=20, c=1, ft, groupCol=8)

x<-reshapeClusterTable(mat, ktable, ft, k=20)
#x$Class<-sapply(x$variable, function(s) ft[ft$Sample_Number == s, 'Class'])
bp<-ggplot(data=x[x$Cluster %in% c(10),], mapping=aes(x=Class, y=value, color=Class)) + 
	geom_boxplot() + 
	facet_grid( . ~Cluster)


