##############################################################################
# File: LIRTS_DEG_Analysis.R                                                 #
# Purpose:                                                                   #
#       Analyze the LEC Injury time series for differential gene expression. #
#       Use the edgeR exactTest and Quasi Liklihood testing frameowrks       #
#       Generate results for several different cuts of the data.             #
#           - DBI Experiment: 3 contrasts, 24v0, 48v0, 48v24                 #
#           - DNA Link Experiment 1: 3 contrasts, 6v0, 24v0, 24v6            #
#           - DNA Link Experiment 2: 1 contrast, 120v0                       #
#           - DNA Link Experiment 3: 1 contrast, 72vs0                       #
#           - 0 and 24 Hours: 2way Factor Interaction model DBI & DNA1       #
#           - Global Wildtype: All subseqeuent intervals vs 0                #
#           - Fibronectin: DBI samples, WT vs FNcKO, 0 vs 48 Hours           #
#           - Beta 8 Integrin: DNA1 samples, WT vs B8cKO, 0 vs 24 Hours      #
#                                                                            #
#       For each of the above experiments generate a DEG Table, a summary    #
#       table, BCV and PCA plots.                                            #
#                                                                            #
#       This script also generates several TMM scaled FPKM Matrices for      #
#       several subsets of the global wildtype analysis including:           #
#           - All genes meeting minimum detection criteria under both        #
#             library chemistries                                            #
#           - Genes in the Upper 50th percentile of tagwise dispersion       #
#           - Excluding genes with a significant batch dependent response    #
#             greater than 4 fold between DBI samples and DNA Link samples   #
#                                                                            #
#                                                                            #
#                                                                            #
#                                                                            #
# Created: June 11, 2019                                                     #
# Updates:                                                                   #
#     April 30, 2021 Major Changes to All Sections                           #
#     Switched to counts generated using Rn45s ribofiltered alignments.      #
#     Updated the Ensmbl annotation hub target from Mm 98 to Mm 101. Updated #
#     the Annotation section to reflect use the Mm101 length table.          #
#     Previously, FPKM were calculated prior to estimation of size factors.  #
#     In this analysis, FPKM are calculated AFTER size factors have been     #
#     estimated.                                                             #
#                                                                            #
#     April 30, 2021 Section: Global Analysis -- Wildype Time Series         #
#     Use filter by expression to exclude genes that are not adequately      #
#     measured in all contrasts.  Improve visualzation of PCA plots by       #
#     assigning batch to shape and interval to color. Generate PCA plots     #
#     for two different filtering strategies: greatest time dependent        #
#     differential expression, least batch dependent differential            #
#     expression. Also plot batch corrected PCA.                             #
#                                                                            #
#     May 3, 2021 Section: Global Analysis -- Wiltype Time Series            #
#     Added heatmap plots.                                                   #
#                                                                            #
#     May 10, 2021 Section: Global Analysis -- Wildtype Time Series          #
#     Revised filter by expression step to consider DNA1 and DNA2 at the     #
#     same time, in order to avoid excluding genes that downregulate by      #
#     120 hours, but upregulate between 6 and 24 hours. Added "library"      #
#     attribute to distinguish between single-end polyA libraries generated  #
#     by DBI and paired end ribodepletion libraries from DNA Link            #
#                                                                            #
#     July 10, 2021 Major Change -- Synapse Integration                      #
#     Changed approach to DGEList construction to allow for Synapse          #
#     integration. This script will check for the existence of a master      #
#     DGEList, and call LIRTS_Push_Synapse.R if none exists                  #
#                                                                            #
#     July 10, 2021 Major Change -- Updated Genome to Ensembl v104 / GRCm39  #
#
#     July 11                                                                #
#                                                                            #
#     TODO:                                                                  #
#     Define Synapse organization, and setup integration                     #
#                                                                            #
#                                                                            #                                      
# Author: Adam Faranda                                                       #
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
local_data_dir <- paste0(wd,"/LIRTS_Raw_Data")        ## Local data directory

source('transcriptomic_analysis_scripts/Overlap_Comparison_Functions.R')
source('transcriptomic_analysis_scripts/LIRTS_Wrap_DEG_Functions.R')


# Create directory to store results (if it does not yet exist)
if(!dir.exists('LIRTS_DEG_Analysis_results'))
{
  dir.create('LIRTS_DEG_Analysis_results')
}

########################### Fetch Master DGE List ############################

fn <- "LIRTS_master_dgelist.Rdata"
syn_dge <- synFindEntityId(name=fn, parent=syn_count_dir)
if(file.exists(paste0("LIRTS_Raw_Data/",fn))){
  print("Loading Local DGEList")
  load(paste0("LIRTS_Raw_Data/",fn))
} else if(!is.null(syn_dge)){
  print("Fetching DGEList from Synapse")
  synGet(
    syn_dge,
    downloadLocation = "LIRTS_Raw_Data"
  )
  load(paste0("LIRTS_Raw_Data/",fn))
} else {
  print("Assembling DGEList from Counts")
  source("transcriptomic_analysis_scripts/LIRTS_Fetch_Count_Data.R")
}

######################### Analyze Wildtype Samples ###########################
# DBI Analysis ####
dge<-master[,
            master$samples$genotype == 'WT' &    # Select Samples
              master$samples$batch == 'DBI' 
            ]
dge$samples$hours_pcs<-droplevels(dge$samples$hours_pcs)
dge$samples$batch<-droplevels(factor(dge$samples$batch))
dge$samples$genotype<-droplevels(factor(dge$samples$genotype))


# Define Experimental Design
design<-model.matrix(~0+group, dge$samples)     
colnames(design)<-gsub(
  'group', '', 
  colnames(design)
)

# Define Contrasts
cntmat<-makeContrasts(                         
  WT24vs0H = WT_24H_DBI - WT_0H_DBI,
  WT48vs0H = WT_48H_DBI - WT_0H_DBI,
  WT48vs24H = WT_48H_DBI - WT_24H_DBI,
  levels = design
)

# Normalize DGEList and fit model
obj <- process_edgeR_ByDesign(
  dge,
  design=design
)

# Create Diagnostic Figures
diagnostic_plots(
  dge=obj$dge, prefix = "DBI_Wildtype"
)

## Iterate over contrasts and prepare summaries / DEG tables
df <- data.frame()
deg <- data.frame()
res <- iterate_edgeR_pairwise_contrasts(
  obj[[1]], obj[[2]], cntmat, df=df, design=design,
  deg=deg, prefix="DBI_Wildtype"
)

# DNA Link Analysis -- First Experiment ####
dge<-master[,
            master$samples$genotype == 'WT' &    # Select Samples
              master$samples$batch == 'DNA1' 
            ]
dge$samples$hours_pcs<-droplevels(dge$samples$hours_pcs)
dge$samples$batch<-droplevels(factor(dge$samples$batch))
dge$samples$genotype<-droplevels(factor(dge$samples$genotype))

design<-model.matrix(~0+group, dge$samples)     # Define Experimental Design
colnames(design)<-gsub(
  'group', '', 
  colnames(design)
)
cntmat<-makeContrasts(
  WT6vs0H = WT_6H_DNA1 - WT_0H_DNA1,
  WT24vs0H = WT_24H_DNA1 - WT_0H_DNA1,
  WT24vs6H = WT_24H_DNA1 - WT_6H_DNA1,
  levels = design
)

# Normalize DGEList and fit model
obj <- process_edgeR_ByDesign(
  dge,
  design=design
)

# Create Diagnostic Figures
diagnostic_plots(
  dge=obj$dge, prefix = "DNA1_Wildtype"
)

## Iterate over contrasts and prepare DEG / summary tables
df<-res[[1]]
deg <- res[[2]]
res <- iterate_edgeR_pairwise_contrasts(
  obj[[1]], obj[[2]], cntmat, df=df, design=design,
  deg=deg, prefix="DNA1_Wildtype"
)

# DNA Link Analysis --  Second Experiment ####
dge<-master[,
            master$samples$genotype == 'WT' &    # Select Samples
              master$samples$batch == 'DNA2' 
            ]
dge$samples$hours_pcs<-droplevels(dge$samples$hours_pcs)
dge$samples$batch<-droplevels(factor(dge$samples$batch))
dge$samples$genotype<-droplevels(factor(dge$samples$genotype))


design<-model.matrix(~0+group, dge$samples)     # Define Experimental Design
colnames(design)<-gsub(
  'group', '', 
  colnames(design)
)
cntmat<-makeContrasts(
  WT120vs0H = WT_120H_DNA2 - WT_0H_DNA2,
  levels = design
)

# Normalize DGEList and fit model
obj <- process_edgeR_ByDesign(
  dge,
  design=design
)

# Create Diagnostic Figures
diagnostic_plots(
  dge=obj$dge, prefix = "DNA2_Wildtype"
)

## Iterate over contrasts and prepare summary figures / tables
df <- res[[1]]
deg <- res[[2]]
res <- iterate_edgeR_pairwise_contrasts(
  obj[[1]], obj[[2]], cntmat, df=df, design=design,
  deg=deg, prefix="DNA2_Wildtype"
)

# 0 vs 24 Hours at 2 labs -- Wildtype Time Series ####
dge<-master[,
            master$samples$genotype == 'WT' &                # Select Samples
              master$samples$hours_pcs %in% c('0H', '24H') &
              master$samples$batch %in% c('DBI', 'DNA1')
            ]
dge$samples$hours_pcs<-droplevels(dge$samples$hours_pcs)
dge$samples$batch<-droplevels(factor(dge$samples$batch))
dge$samples$genotype<-droplevels(factor(dge$samples$genotype))

# Define Experimental Design: Interval, Batch, and Batch Interaction.
design<-model.matrix(~hours_pcs*batch, dge$samples)     

# Normalize DGEList and fit model
obj <- process_edgeR_ByDesign(
  dge, design=design
)

# Create Diagnostic Figures
diagnostic_plots(
  dge=obj$dge, prefix = "24_Hour_LabComp_TimeBatchInx",
  color_attrib = "hours_pcs", shape_attrib = "batch"
)

## Iterate over contrasts and prepare summary figures / tables
df<-res[[1]]
deg <- res[[2]]
res <- iterate_edgeR_design_coefficients(
  obj[[1]], obj[[2]], coefs=c(2,3,4), df=df, design=design,
  deg=deg, prefix="24_Hour_LabComp_TimeBatchInx", group_label_list = list(
    c("0H", "24H"), c("DBI", "DNA1"), c("Intercept", "InxDNA1_24H")
  )
)

# Global Analysis -- Wildtype Time Series ####

## Setup library variable (as that seems to be the main batch effect driver)
master$samples$library <- factor(
  ifelse(master$samples$batch == "DBI", "single", "paired"),
  levels=c("single", "paired")
)

## Extract a set of genes that are minimally detected under both library
## preparation schemes.
genes <- master$genes$gene_id
for(b in unique(master$samples$library)){
  dge<-master[,
              master$samples$library == b & 
                master$samples$genotype == 'WT'
              ]    # Select Samples
  dge$samples$group<-droplevels(  
    dge$samples$hours_pcs
  )
  print(unique(dge$samples$library))
  design<-model.matrix(~hours_pcs, dge$samples)    
  colnames(design)<-gsub("hours_pcs", '', colnames(design))
  
  genes<-intersect(
    genes,
    (dge[filterByExpr(dge, design), ,keep.lib.sizes=F])$genes$gene_id
  )
  print(length(genes))
}
dge <- master[,master$samples$genotype=="WT"]
dge$samples$hours_pcs<-droplevels(dge$samples$hours_pcs)
dge$samples$batch<-droplevels(factor(dge$samples$batch))
dge$samples$genotype<-droplevels(factor(dge$samples$genotype))

# Define Experimental Design (full model)
design<-model.matrix(~library + hours_pcs, dge$samples)   
#colnames(design)<-gsub("hours_pcs", '', colnames(design))

obj <- process_selected_features(
  dge=dge, design=design, genes=genes, color_attrib = "hours_pcs",
  shape_attrib = "batch", prefix="Global_Wildtype"
)

# Iterate over model coefficients and generate DEG tables
df<-data.frame()
deg <- data.frame()
res <- iterate_edgeR_design_coefficients(
  obj[[1]], obj[[2]], coefs=2:7, df=df, design=design,
  deg=deg, prefix="Global_Wildtype", group_label_list = list(
    c("single", "paired"),c("0H", "6H"),
    c("0H", "24H"),c("0H", "48H"), c("0H", "72H"),c("0H", "120H")
  )
)

## Construct Feature Selection Table
feature_selection <- obj$dge$genes %>% 
  dplyr::select(gene_id, matches("_FPKM"))
feature_selection$tagwise.dispersion <- obj$dge$tagwise.dispersion
feature_selection <- feature_selection%>%
  inner_join(
    res[[2]] %>% filter(Samples=="Global_Wildtype") %>%
      mutate(Contrast=paste0(Group_2,"v", Group_1)) %>%
      dplyr::select(gene_id, Contrast, logFC, FDR) %>%
      tidyr::pivot_wider(
        id_cols="gene_id",
        names_from="Contrast",
        values_from=c("logFC", "FDR")
      ) %>% as.data.frame(),
    by="gene_id"
  )

## Extract Genes in the top 50th percentile of tag-wise dispersion
genes.disp <- feature_selection %>%
  filter(tagwise.dispersion > quantile(tagwise.dispersion, 0.5)) %>% 
  pull("gene_id")

process_selected_features(
  dge=dge, genes = genes.disp, design=design,
  prefix="Global_Wildtype_Upper50_Disp"
)

## Exclude genes with a significant difference between the DBI batch
## and DNA link batches.
# genes.batch <- union(
#   union(
#     feature_selection %>%
#       #filter(tagwise.dispersion > quantile(tagwise.dispersion, 0.9)) %>% 
#       filter(abs(logFC_DNA1vDBI) > 1 & FDR_DNA1vDBI < 0.1) %>% 
#       pull("gene_id"),
#     feature_selection %>%
#       #filter(tagwise.dispersion > quantile(tagwise.dispersion, 0.9)) %>% 
#       filter(abs(logFC_DNA2vDBI) > 1 & FDR_DNA2vDBI < 0.1) %>% 
#       pull("gene_id")
#   ),
#   feature_selection %>%
#     #filter(tagwise.dispersion > quantile(tagwise.dispersion, 0.9)) %>% 
#     filter(abs(logFC_DNA3vDBI) > 1 & FDR_DNA3vDBI < 0.1) %>% 
#     pull("gene_id")
# )

genes.batch <-  feature_selection %>%
  filter(abs(logFC_pairedvsingle) > 0.5 & FDR_pairedvsingle < 0.1) %>%
  pull("gene_id")

process_selected_features(
  dge=dge, genes = setdiff(genes, genes.batch), design=design,
  prefix="Global_Wildtype_NoBatchSig"
)

## Select genes with a significant injury response at any interval
genes.injury <- feature_selection %>%
  filter(
    (abs(logFC_6Hv0H) > 1 & FDR_6Hv0H < 0.05) |
      (abs(logFC_24Hv0H) > 1 & FDR_24Hv0H < 0.05) |
      (abs(logFC_48Hv0H) > 1 & FDR_48Hv0H < 0.05) |
      (abs(logFC_72Hv0H) > 1 & FDR_72Hv0H < 0.05) |
      (abs(logFC_120Hv0H) > 1 & FDR_120Hv0H < 0.05)
  )%>%
  pull("gene_id")

process_selected_features(
  dge=dge, genes = genes.injury, design=design,
  prefix="Global_Wildtype_InjurySig"
)

## Use Combat-Seq to generate batch corrected counts
design.exp <- setdiff(
  colnames(design), c('librarypaired')
)

bat.counts <- ComBat_seq(
  dge$counts,
  batch=dge$samples$batch,
  covar_mod = design[,design.exp],
  full_mod = T
)

obj.bat <- process_selected_features(
  dge=dge, genes = genes, design=design,
  prefix="Global_Wildtype_CombatSeq",
  counts=bat.counts
)

# Store Tagwise dispersions for batch corrected counts in 'twd'
twd <- obj.bat$dge$tagwise.dispersion
genes.bat.disp <- obj.bat$dge$genes$gene_id[twd > quantile(twd, 0.5)]
process_selected_features(
  dge=dge, genes = genes.bat.disp, design=design,
  prefix="Global_Wildtype_CombatSeq_Upper50_Disp",
  counts=bat.counts
)


######################## Analyze Mutation effects ############################
# DBI -- Fibronectin vs Wildtype 0 and 48 hours ####
dge<-master[,
              master$samples$hours_pcs %in% c('0H', '48H')&    # Select Samples
              master$samples$batch == 'DBI',
            ]
dge$samples$genotype<-droplevels(dge$samples$genotype)
dge$samples$hours_pcs<-droplevels(dge$samples$hours_pcs)
dge$samples$batch<-droplevels(factor(dge$samples$batch))

# Define Experimental Design for pairwise contrasts by group
design<-model.matrix(~0+group, dge$samples) 
colnames(design)<-gsub(
  'group', '', 
  colnames(design)
)

cntmat<-makeContrasts(
  WT48vs0H = WT_48H_DBI - WT_0H_DBI,
  FN48vs0H = FN_48H_DBI - FN_0H_DBI,
  FNvsWT_0H = FN_0H_DBI - WT_0H_DBI,
  FNvsWT_48H = FN_48H_DBI - WT_48H_DBI,
  levels = design
)

# Normalize DGEList and fit model
obj <- process_edgeR_ByDesign(
  dge,
  design=design
)

# Create Diagnostic Figures
diagnostic_plots(
  dge=obj$dge, prefix = "Fn1_Groupwise"
)

## Iterate over contrasts and prepare DEG / summary tables
df<-res[[1]]
deg <- res[[2]]
res <- iterate_edgeR_pairwise_contrasts(
  obj[[1]], obj[[2]], cntmat, df=df, design=design,
  deg=deg, prefix="Fn1_Groupwise"
)

# Define Experimental Design for factorial design hours_pcs by genotype
# With an interaction term
design<-model.matrix(~hours_pcs * genotype, dge$samples)   

obj <- process_edgeR_ByDesign(
  dge,
  design = design
)

# Create Diagnostic Figures
diagnostic_plots(
  dge=obj$dge, prefix = "Fn1_FactorInx",
  color_attrib = "hours_pcs", shape_attrib = "genotype"
)

# Iterate over model coefficients and generate DEG tables
df<-res[[1]]
deg <- res[[2]]
res <- iterate_edgeR_design_coefficients(
  obj[[1]], obj[[2]], coefs=2:4, df=df, design=design,
  deg=deg, prefix="Fn1_FactorInx", group_label_list = list(
    c("0H", "48H"), c("WT", "FN"),  c("Intercept", "Inx48H_FN")
  )
)

# DNA Link -- Beta 8 Integrin vs Wildtype 0 and 48 hours ####
dge<-master[,
            master$samples$hours_pcs %in% c('0H', '24H')&    # Select Samples
              master$samples$batch == 'DNA1',
            ]
dge$samples$genotype<-droplevels(dge$samples$genotype)
dge$samples$hours_pcs<-droplevels(dge$samples$hours_pcs)
dge$samples$batch<-droplevels(factor(dge$samples$batch))

design<-model.matrix(~0+group, dge$samples)     # Define Experimental Design
colnames(design)<-gsub(
  'group', '', 
  colnames(design)
)

cntmat<-makeContrasts(
  WT24vs0H = WT_24H_DNA1 - WT_0H_DNA1,
  B824vs0H = B8_24H_DNA1 - B8_0H_DNA1,
  B8vsWT_0H = B8_0H_DNA1 - WT_0H_DNA1,
  B8vsWT_24H = B8_24H_DNA1 - WT_24H_DNA1,
  levels = design
)

# Normalize DGEList and fit model
obj <- process_edgeR_ByDesign(
  dge,
  design=design
)

# Create Diagnostic Figures
diagnostic_plots(
  dge=obj$dge, prefix = "Itgb8_Groupwise"
)

## Iterate over contrasts and prepare DEG / summary tables
df<-res[[1]]
deg <- res[[2]]
res <- iterate_edgeR_pairwise_contrasts(
  obj[[1]], obj[[2]], cntmat, df=df, design=design,
  deg=deg, prefix="Itgb8_Groupwise"
)

# Define Experimental Design for factorial design hours_pcs by genotype
# With an interaction term
design<-model.matrix(~hours_pcs * genotype, dge$samples)   

obj <- process_edgeR_ByDesign(
  dge,
  design = design
)

# Create Diagnostic Figures
diagnostic_plots(
  dge=obj$dge, prefix = "Itgb8_FactorInx",
  color_attrib = "hours_pcs", shape_attrib = "batch"
)

# Iterate over model coefficients and generate DEG tables
df<-res[[1]]
deg <- res[[2]]
res <- iterate_edgeR_design_coefficients(
  obj[[1]], obj[[2]], coefs=2:4, df=df, design=design,
  deg=deg, prefix="Itgb8_FactorInx", group_label_list = list(
    c("0H", "24H"), c("WT", "B8"),  c("Intercept", "Inx24H_B8")
  )
)

save(res, file="LIRTS_DEG_Analysis_Results/LIRTS_Master_DEG_Table.Rdata")

