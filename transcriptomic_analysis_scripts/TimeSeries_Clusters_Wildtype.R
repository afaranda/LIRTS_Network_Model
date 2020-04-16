##############################################################################
# File: TimeSeries_Clusters_Wildtype.R                                       #
# Created: April 14, 2020                                                    #
# Author: Adam Faranda                                                       #
# Purpose:                                                                   #   
#         Apply various clustering methods to try and identify gene          #
#         Gene clusters with distinct temporal profiles                      #   
#                                                                            #
##############################################################################

############################ Setup Environment ###############################

library(edgeR)
library(pheatmap)
library(TMixClust)
setwd('/home/adam/Documents/LEC_Time_Series')
source('transcriptomic_analysis_scripts/ClusteringFunctions.R')
min_lfc=4                                          # log fold change threshold
min_cpm=0.50                                       # Minimum overall abundance
max_cpm=20                                        # Maximum overall abundance
  output_prefix="lfc4"                        # Prefix for output files
dgeFile="LTS_DGEList.Rdata"                   # Rdata file with dgelist object

## Helper Function filters a count matrix by given criteria
filterCPMmat<-function(
  object = dge, samples = wt_samples, func=cpm,
  lfc = 3, min_log_cpm = 0, max_log_cpm = 100, include_genes = biosig$gene_id, 
  lfc_columns = paste('logFC', LETTERS[1:4], sep="."), deg_table = deg
){
  ecpm<-func(object, log = T)        # Generate a gene x sample log CPM matrix
  
  # Filter differential expression table by fold change and abundance
  dt<-deg_table[
    apply(deg_table[,lfc_columns], 
          1, function(x, l=lfc) any(abs(x) > l)),
    ]
  dt<-dt[dt$logCPM > min_log_cpm & dt$logCPM < max_log_cpm, ]
  dt<-dt[dt$FDR < 0.05, ]
  
  # If a subset of genes is provided, take the intersection
  if (!is.null(include_genes)){
    genes<-intersect(dt$gene_id, include_genes)
  } else{
    genes<-dt$gene_id
  }
  return(ecpm[genes, samples])
}


######################### Load and Prepare Data ##############################
load(dgeFile)
s<-master$samples[master$samples$genotype == "WT",'sample']
dge<-master[,s]
design<-model.matrix(~interval + batch, dge$samples)
colnames(design)<-gsub("interval", '', colnames(design))
keep<-filterByExpr(dge, design)
dge<-dge[keep,,keep.lib.sizes=F]

# Normalize and Estimate Dispersions
dge<-calcNormFactors(dge)
dge<-estimateDisp(dge, design, robust = T)

# Calcuate Statistical Significance (Gene DE at ANY Timepoint)
fit<-glmQLFit(dge, design)
qlf<-glmQLFTest(fit, coef=2:5)
deg<-as.data.frame(topTags(qlf, n=Inf))

# construct datExpr and datTraits using filteredGenes
datExpr<-t(
  filterCPMmat(
  object = dge, samples=dge$samples$sample, lfc=min_lfc,
  min_log_cpm = min_cpm, max_log_cpm = max_cpm, func=cpm,
  lfc_columns = paste('logFC', c('6H', '24H', '48H', '120H'), sep="."),
  deg_table = deg, include_genes = NULL
  )
)

# remove columns that hold information we do not need.
allTraits = dge$samples[, -c(1,2,3,5,6,8,9,10,11)];
row.names(allTraits)<-allTraits$sample
modbatch<-model.matrix(~ 0 + batch, allTraits)
modhours<-model.matrix(~ 0 + interval, allTraits)

allTraits<-merge(allTraits, modhours, by='row.names')
row.names(allTraits)<-allTraits$Row.names
allTraits<-allTraits[-1]

allTraits<-merge(allTraits, modbatch,by='row.names')
rownames(allTraits)<-allTraits$sample
allTraits<-allTraits[-1]

allTraits<-allTraits[,-c(1,2,3)]
names(allTraits)<-gsub("batch", "", names(allTraits))
names(allTraits)<-gsub("interval", "", names(allTraits))


####################### Generate Heirarchical Clusters #######################

h<-hclust(dist(t(datExpr), method = "manhattan"), method = 'complete' )
tab<-tabulate_H_Clusters(h,ks=c(1:10))

#################### Generate Line Plots For Each Module #####################
library(data.table)
library(ggplot2)
dt<-data.table(merge(datExpr, dge$samples, by='row.names'))
for(c in 1:10){
  g<-row.names(tab[tab$k_eq_10 == c,])
  x<-melt(
    dt[,c(g, 'interval'),with=F][,lapply(.SD, mean), by=interval],
    id.vars = 'interval'
  )
  fn<-paste(
    output_prefix,"Cluster_Number",
    c, "HClust_Lineplot.png",
    sep="_"
  )
  ttl<-paste("cluster", c, "size", length(g))
  ggplot(x, aes(x=interval, y=value, color=variable, group=variable)) + 
    geom_line() + theme(legend.position = 'none') + ggtitle(ttl)
  ggsave(fn)
}

######################### Run A TMixClust Analysis ###########################
ts_dat <- cpmByGroup(dge[colnames(datExpr),], group=dge$samples$interval)

ts_clust_3 <-TMixClust(ts_dat, nb_clusters = 3)
ts_clust_4 <-TMixClust(ts_dat, nb_clusters = 4)
ts_clust_5 <-TMixClust(ts_dat, nb_clusters = 5)
ts_clust_6 <-TMixClust(ts_dat, nb_clusters = 6)

save(
  ts_clust_3, ts_clust_4, ts_clust_5, ts_clust_6,
  file = "TMix_LTS_Clusters.Rdata"
)






