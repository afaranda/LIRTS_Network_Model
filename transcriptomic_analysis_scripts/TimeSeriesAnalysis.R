################################################################################
# File: TimeSeriesAnalysis.R                                                   #
# Purpose: Identify gene clusters that correspond to temporal patterns.        #
# Created: May 1, 2019                                                         #
# Author: Adam Faranda                                                         #
################################################################################

############################ Setup Environment #################################
setwd('/home/adam/Documents/LEC_Time_Series')
library(dplyr)
library(cluster)
library(reshape2)
wd<-getwd()
source('transcriptomic_analysis_scripts/BuildDataMatrix.R')
source('transcriptomic_analysis_scripts/PreprocessingFunctions.R')
source('transcriptomic_analysis_scripts/PrincipalComponents.R')
source('transcriptomic_analysis_scripts/ClusteringFunctions.R')

############################ Load in Data Files ################################

gtfpath<-"/home/adam/Documents/LTS_Data/Mus_musculus.GRCm38.96.chr.gtf"
orig_gtf<-readGFF(gtfpath)

# Manually Annotate the File Table tf wisth the followieng data:
#               Group Information: Genotype, Time Point,  Sequencing Lab
#               Read Information: Average Length, library type (paired / single)

dl<-c(
  "/home/adam/Documents/LTS_Data/DBI_NoTrim_HTSeq_Count_Gene", 
  "/home/adam/Documents/LTS_Data/DNA_Link_NoTrim_HTSeq_Count_Gene"
)

if(length(list.files(pattern = "GeneLengthTable.Rdata"))>0){
  load("GeneLengthTable.Rdata")
} else {
  lt<-lengthTable(gtfpath)
}
ft<-hc_getFileTable(dirList=dl)
ds<-hc_loadFiles(ft)
ft<-hc_identifierConsistency(ds, ft)
htseq_count<-hc_buildDataFrame(ds, ft)
htseq_count<-hc_dropSamples(
  ft = htseq_count[[1]],
  mat = htseq_count[[2]],
  samples = ft[ft$Genotype !="WT", 1]
)
htseq_dge<-buildDGE(htseq_count[[2]], ft=htseq_count[[1]], 
                    gt=data.frame(lt %>% 
                                    group_by(gene_id) %>% 
                                    summarise(Union_Length = dplyr::first(Union_Length)
                                    )
                    )
)

# Import Stringtie Data
dll<-c(
  "/home/adam/Documents/LTS_Data/DBI_NoTrim_StringTie", 
  "/home/adam/Documents/LTS_Data/DNA_Link_NoTrim_StringTie"
)

ftt<-st_getFileTable(dll, wd=".")
stg<-st_loadFiles(ftt, fnCol=4)
bg<-st_sewBallGown(ftt)
stringtie_tran_tpm<-st_buildTranscriptMatrix(ds=stg, idCol = 10, measCol = 14)
stringtie_gene_tpm<-st_buildGeneMatrix(ds=stg, g_idCol=9, measCol = 14)
ls()


######################### Apply Preprocessing Steps ############################

# Add Class attribute to feature definition table
ft<-htseq_count[[1]]
ft$Class<-as.factor(paste("Hour",ft$Hours_PCS,sep=""))
ft$Class<-factor(
	ft$Class, 
	levels=unique(
		as.character(ft$Class[order(ft$Hours_PCS)])
	)
)
mat<-htseq_count[[2]]
genes<-lt %>% group_by(gene_id) %>% filter(row_number()==1)

ecpm<-edgeRcpm(mat)                     # Normalize using edgeR's TMM method
etpm<-normGeneLength(ft=ft, mat=mat, gt=genes)
ecmb<-wrapCombat(ecpm, ft, groupCol = 9, batchCol = 6, idCol = 0)# Correct for batch effects
ecmb<-fixCombatNegatives(ecmb, idCol=0)             # replace negative values with min +ve
etpm[etpm==0]<-0.0001
etmb<-ComBat(
  dat=as.matrix(etpm), 
  batch=ft$Lab,
  mod = model.matrix(~1, ft)
)
etmb<-fixCombatNegatives(etmb, idCol =0)

############## Apply Variance Filters; Plot Principal Components ###############
varRanks <-c(10, 50, 100, 10000)                # Try different variance filters
for( v in varRanks){
	print(v)
	ecpm.filter <-varianceFilter(etpm, threshold=v)
	ecmb.filter <-varianceFilter(etmb, threshold=v)
	
	# Plot results for TMM Normalized Data
	f1<-paste('ECPM_Samples_Top_', v,'_Ranked_Class.png')
	f2<-paste('ECPM_Samples_Top_', v,'_Ranked_Lab.png')
	png(f1, width=240, height=150)
		print(plotPrinComp(ecpm.filter, ft, groupCol=9, idCol=0))
	dev.off()
	
	png(f2, width=240, height=150)
		print(plotPrinComp(ecpm.filter, ft, groupCol=6, idCol=0))
	dev.off()

	# Plot results for batch adjusted TMM data
	f1<-paste('ECMB_Samples_Top_', v,'_Ranked_Class.png')
	f2<-paste('ECMB_Samples_Top_', v,'_Ranked_Lab.png')
	png(f1, width=240, height=150)
		print(plotPrinComp(ecmb.filter, ft, groupCol=9, idCol=0))
	dev.off()
	png(f2, width=240, height=150)
		print(plotPrinComp(ecmb.filter, ft, groupCol=6, idCol=0))
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
mat<-varianceFilter(etmb, threshold=v)
mat.l<-log(mat)
mat.s<-scaleCenterByRow(mat)

d.meth='manhattan'
h.method='complete'
h<-wrapHclust(mat.s, idCol=0, transpose = F)
ktable<-tabulate_H_Clusters(h)
rt<-reshapeClusterTable(mat.s, ktable, 5, ft=ft)

summarizeGeneClusters(mat.s, label='Scaled_Centered_v200', nclust=5)


v = 2000
ecmb.filter<-varianceFilter(ecmb, threshold=v)
mat<-as.matrix(ecmb.filter[,2:ncol(ecmb.filter)])
mat.l<-log(mat)
mat.s<-scaleCenterByRow(mat)
summarizeGeneClusters(mat.s, label='Scaled_Centered_v2000')
dev.



