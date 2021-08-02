##############################################################################
# File: TimeSeries_Clusters_Wildtype.R                                       #
# Created: April 14, 2020                                                    #
# Author: Adam Faranda                                                       #
# Purpose: Use several clustering algorithms (heirarchical, SOTA) 
#                                                                            #
##############################################################################

############################ Setup Environment ###############################
# CHANGE THIS DIRECTORY
setwd('~/Documents/LEC_Time_Series')
#setwd("~/Documents/Adam_LEC_Time_Series_DEG_Analysis_30_Apr_2021")
#load("GeneLengthTable.Rdata")
library(dplyr)
library(cluster)
library(reshape2)
library(pheatmap)
library(sva)
library(synapser)
wd<-getwd()

synapser::synLogin()
syn_project <- "syn25579691"                 ## Synapse ID for this project
syn_count_dir <- "syn25979754"               ## Synapse folder with counts
syn_code_dir <- "syn25976329"                ## Synapse folder with code
syn_sample_table <- "syn25582010"            ## Synapse sample table
syn_gene_meta <- "syn25976328"               ## Gene Annotations
syn_expr_mat <- "syn25992507"                ## Foldeer with FPKM Matrices

local_data_dir <- paste0(wd,"/LIRTS_Raw_Data")        ## Local data directory

source('transcriptomic_analysis_scripts/Overlap_Comparison_Functions.R')
source('transcriptomic_analysis_scripts/LIRTS_Wrap_DEG_Functions.R')

# Create directory to store results (if it does not yet exist)
if(!dir.exists('LIRTS_Cluster_Analysis_results'))
{
  dir.create('LIRTS_Cluster_Analysis_results')
}

## Analysis parameters for feature selection on global wildtype samples
cpm_weight = 0.9           ## Weight applied to minimum cpm threshold
trend_disp_weight = 1.7    ## Weight applied to trended dispersion threshold


######################## Fetch Expression Matrices ###########################
fn <- "GWT_AllPresent_TMM-FPKM_Matrix.csv"
fpkm <- synFindEntityId(name=fn, parent=syn_expr_mat)
if(file.exists(paste0("LIRTS_DEG_Analysis_results/",fn))){
  print("Loading local fpkm matrix")
  fpkm <- read.csv(
    paste0("LIRTS_DEG_Analysis_results/",fn)
  )
} else if(!is.null(syn_dge)){
  print("Fetching fpkm matrix from Synapse")
  synGet(
    fpkm,
    downloadLocation = "LIRTS_DEG_Analysis_results/"
  )
  fpkm <- read.csv(
    paste0("LIRTS_DEG_Analysis_results/",fn)
  )

} else {
  print("Generating fpkm matrix")
  source("transcriptomic_analysis_scripts/LIRTS_DEG_Analysis.R")
  fpkm <- read.csv(
    paste0("LIRTS_DEG_Analysis_results/",fn)
  )
}

fn <- "GWT_APComBatSeq_TMM-FPKM_Matrix.csv"
fpkm_bat <- synFindEntityId(name=fn, parent=syn_expr_mat)
if(file.exists(paste0("LIRTS_DEG_Analysis_results/",fn))){
  print("Loading local fpkm matrix")
  fpkm_bat <- read.csv(
    paste0("LIRTS_DEG_Analysis_results/",fn)
  )
} else if(!is.null(syn_dge)){
  print("Fetching fpkm matrix from Synapse")
  synGet(
    fpkm_bat,
    downloadLocation = "LIRTS_DEG_Analysis_results/"
  )
  fpkm_bat <- read.csv(
    paste0("LIRTS_DEG_Analysis_results/",fn)
  )
  
} else {
  print("Generating fpkm matrix")
  source("transcriptomic_analysis_scripts/LIRTS_DEG_Analysis.R")
  fpkm_bat <- read.csv(
    paste0("LIRTS_DEG_Analysis_results/",fn)
  )
}


##################### Fetch Feature Selection Tables #########################

fn <- "GWT_AllPresent_Feature_Selection_Table.csv"
feature_selection <- synFindEntityId(name=fn, parent=syn_gene_meta)
if(file.exists(paste0("LIRTS_DEG_Analysis_results/",fn))){
  print("Loading local feature selection table")
  feature_selection <- read.csv(
    paste0("LIRTS_DEG_Analysis_results/",fn)
  )
} else if(!is.null(syn_dge)){
  print("Fetching feature selection table from Synapse")
  synGet(
    syn_dge,
    downloadLocation = "LIRTS_DEG_Analysis_results/"
  )
  feature_selection <- read.csv(
    paste0("LIRTS_DEG_Analysis_results/",fn)
  )
  
} else {
  print("Generating feature selection table")
  source("transcriptomic_analysis_scripts/LIRTS_DEG_Analysis.R")
  feature_selection <- read.csv(
    paste0("LIRTS_DEG_Analysis_results/",fn)
  )
}



fn <- "GWT_APComBatSeq_Feature_Selection_Table.csv"
feature_selection <- synFindEntityId(name=fn, parent=syn_gene_meta)
if(file.exists(paste0("LIRTS_DEG_Analysis_results/",fn))){
  print("Loading local feature selection table")
  feature_selection <- read.csv(
    paste0("LIRTS_DEG_Analysis_results/",fn)
  )
} else if(!is.null(syn_dge)){
  print("Fetching feature selection table from Synapse")
  synGet(
    syn_dge,
    downloadLocation = "LIRTS_DEG_Analysis_results/"
  )
  feature_selection <- read.csv(
    paste0("LIRTS_DEG_Analysis_results/",fn)
  )
  
} else {
  print("Generating feature selection table")
  source("transcriptomic_analysis_scripts/LIRTS_DEG_Analysis.R")
  feature_selection <- read.csv(
    paste0("LIRTS_DEG_Analysis_results/",fn)
  )
}


######################### Load and Prepare Data ##############################

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






