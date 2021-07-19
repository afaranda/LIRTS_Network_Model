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
#                                                                            #
#     July 11 - 14, 2021 Major Refactoring                                   #
#     Wrapped most analysis steps in helper functions in the script          #
#     LIRTS_Wrap_DEG_Functions.R that are called from this script.           #
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

## Analysis parameters for feature selection on global wildtype samples
cpm_weight = 0.9           ## Weight applied to minimum cpm threshold
trend_disp_weight = 1.7    ## Weight applied to trended dispersion threshold



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
df <- res[[1]]
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
  shape_attrib = "batch", prefix="GWT_AllPresent"
)

# Iterate over model coefficients and generate DEG tables
df <- res[[1]]
deg <- res[[2]]
res <- iterate_edgeR_design_coefficients(
  obj[[1]], obj[[2]], coefs=2:7, df=df, design=design,
  deg=deg, prefix="GWT_AllPresent", group_label_list = list(
    c("single", "paired"),c("0H", "6H"),
    c("0H", "24H"),c("0H", "48H"), c("0H", "72H"),c("0H", "120H")
  )
)

## Construct Feature Selection Table

# Fetch list of genes that we have previously considered important
fn <- "Key_Injury_Response_Genelist.txt"
key_lens_injury <- synFindEntityId(fn, parent = syn_gene_meta)
if(file.exists(paste0("LIRTS_Raw_Data/",fn))){
  key_lens_injury <- read.table(paste0("LIRTS_Raw_Data/", fn))[,1]
} else{
  synGet(
    key_lens_injury, 
    downloadLocation = "LIRTS_Raw_Data"
  )
  key_lens_injury <- read.table(paste0("LIRTS_Raw_Data/", fn))[,1]
}

# Fetch list of genes that are assigned to the GO term
# Lens Development in Camera Type Eye
fn <- "LensCameraTypeEye_GO_0002088_Genelist.txt"
lens_dev <- synFindEntityId(fn, parent = syn_gene_meta)
if(file.exists(paste0("LIRTS_Raw_Data/",fn))){
  lens_dev <- read.table(paste0("LIRTS_Raw_Data/", fn))[,1]
} else{
  synGet(
    lens_dev, 
    downloadLocation = "LIRTS_Raw_Data"
  )
  lens_dev <- read.table(paste0("LIRTS_Raw_Data/", fn))[,1]
}

## Define a table with potential parameters for feature selection
feature_selection <- obj$dge$genes %>% 
  dplyr::select(gene_id, SYMBOL, matches("_FPKM")) %>%
  mutate(
    Tag_disp = obj$dge$tagwise.dispersion,
    AveLogCPM = obj$dge$AveLogCPM,
    lens_dev = ifelse(SYMBOL %in% lens_dev, TRUE, FALSE),
    key_lens_injury = ifelse(SYMBOL %in% key_lens_injury, TRUE, FALSE),
    Trend_disp=obj$dge$trended.dispersion
  ) %>%
  inner_join(
    res[[2]] %>% filter(Samples=="GWT_AllPresent") %>%
      mutate(Contrast=paste0(Group_2,"v", Group_1)) %>%
      dplyr::select(gene_id, Contrast, logFC, FDR) %>%
      tidyr::pivot_wider(
        id_cols="gene_id",
        names_from="Contrast",
        values_from=c("logFC", "FDR")
      ) %>% as.data.frame(),
    by="gene_id"
  ) %>%
  mutate(
    Sig_Any_Interval=(
      abs(logFC_6Hv0H) > 1 & FDR_6Hv0H < 0.05)|
      (abs(logFC_24Hv0H) > 1 & FDR_24Hv0H < 0.05)|
      (abs(logFC_48Hv0H) > 1 & FDR_48Hv0H < 0.05)|
      (abs(logFC_72Hv0H) > 1 & FDR_72Hv0H < 0.05)|
      (abs(logFC_120Hv0H) > 1 & FDR_120Hv0H < 0.05)
  )

## Check that all features where "Sig_Any_Interval == TRUE have
## at least one observed absolute log2 fold change > 1
print(
  feature_selection %>%
    filter(Sig_Any_Interval) %>%
    select(gene_id, matches("logFC_[12467]")) %>%
    rowwise() %>%
    mutate(
      MaxAbs=max(
        abs(logFC_6Hv0H), 
        abs(logFC_24Hv0H),
        abs(logFC_48Hv0H),
        abs(logFC_72Hv0H),
        abs(logFC_120Hv0H)
      )
    ) %>% pull(MaxAbs) %>% min()
)

## Determine the minimum value of AveLogCPM for genes that we have 
## previously queried, or lens development genes. 
key_genes <- feature_selection %>%
  filter(key_lens_injury  | lens_dev ) %>%
  pull('gene_id')

cpm_cut <- feature_selection%>%
  filter(gene_id %in% key_genes & Sig_Any_Interval) %>%
  pull("AveLogCPM") %>% min()


## Add selection attribute to feature selection table (Values = TRUE, FALSE)
feature_selection <- feature_selection %>%
  mutate(
    Selected = ifelse(
      (AveLogCPM > cpm_cut *cpm_weight) & 
        (Tag_disp > trend_disp_weight * Trend_disp),
      TRUE, FALSE
    )
  )

write.csv(
  feature_selection,
  file="LIRTS_DEG_Analysis_results/GWT_AllPresent_Feature_Selection_Table.csv"
)

## Highlight key genes with significant differential expression (lfc > 1)
pdf(
  file=paste0(
    "LIRTS_DEG_Analysis_results/",
    "GWT_AllPresent_LFC1_Significant_Key_Genes_BCV_Plot.pdf"
  ), 
  width = 6, height = 6
)
plot_BCV_flag(
  obj$dge,
  cex=c(1, 0.1),
  flag=(
    feature_selection$Sig_Any_Interval & 
      feature_selection$gene_id %in% key_genes
  ),
  flag_labels = c("Tagwise - Key/Sig.", "Tagwise - Other")
)
dev.off()

## Show genes with a strong response at any post-surgical interval
pdf(
  file="LIRTS_DEG_Analysis_results/GWT_LFC3_Significant_Genes_BCV_Plot.pdf", 
  width = 6, height = 6
)
plot_BCV_flag(
  obj$dge,
  cex=c(1, 0.1),
  flag=(
    (abs(feature_selection$logFC_6Hv0H) > 3 & feature_selection$FDR_6Hv0H < 0.05)|
      (abs(feature_selection$logFC_24Hv0H) > 3 & feature_selection$FDR_24Hv0H < 0.05)|
      (abs(feature_selection$logFC_48Hv0H) > 3 & feature_selection$FDR_48Hv0H < 0.05)|
      (abs(feature_selection$logFC_72Hv0H) > 3 & feature_selection$FDR_72Hv0H < 0.05)|
      (abs(feature_selection$logFC_120Hv0H) > 3 & feature_selection$FDR_120Hv0H < 0.05)
  ),
  flag_labels = c("Tagwise - DEG", "Tagwise Other Gene"),
  main="Genes with 8-fold or greater DE at any interval after 0H"
)
dev.off()


## Query on "feature_selection Summary table for DEG counts at different
## logFC thresholds.
write.csv(
  feature_selection %>%
    select(
      gene_id, matches("logFC_[12467]"),
      matches("FDR_[12467]")
    ) %>%
    tidyr::pivot_longer(
      cols=matches('(logFC|FDR)_[12467]'),
      names_to=c(".value", "Interval"),
      names_sep="_"
    ) %>%
    group_by(Interval) %>%
    mutate(
      Interval=factor(
        gsub("logFC_","",Interval),
        levels=c(
          "6Hv0H",
          "24Hv0H",
          "48Hv0H",
          "72Hv0H",
          "120Hv0H",
          "Any"
        )
      )
    )%>%
    summarise(
      LogFC_1=sum(abs(logFC)>1 & FDR < 0.05),
      LogFC_3=sum(abs(logFC)>3 & FDR < 0.05),
      LogFC_5=sum(abs(logFC)>5 & FDR < 0.05)
      
    )%>%
    bind_rows(
      feature_selection %>%
        select(
          gene_id, matches("logFC_[12467]"),
          matches("FDR_[12467]")
        ) %>%
        tidyr::pivot_longer(
          cols=matches('(logFC|FDR)_[12467]'),
          names_to=c(".value", "Interval"),
          names_sep="_"
        ) %>%
        filter(FDR < 0.05) %>%
        group_by(gene_id) %>%
        filter(abs(logFC)==max(abs(logFC))) %>%
        filter(row_number() ==1) %>%
        group_by()%>%
        summarise(
          Interval= "Any",
          LogFC_1 = sum(abs(logFC)>1),
          LogFC_3 = sum(abs(logFC)>3),
          LogFC_5 = sum(abs(logFC)>5)
        )
    ), 
  file="LIRTS_DEG_Analysis_results/GWT_DEG_Count_Totals_by_lfc_threshold.csv"
)

process_selected_features(
  dge, genes = feature_selection %>%filter(Selected) %>%pull("gene_id"),
  design=design, prefix = "GWT_Selected"
)

pdf(
  file=paste0(
    "LIRTS_DEG_Analysis_results/",
    "GWT_AllPresent_Selected_Genes_BCV_Plot.pdf"
  ), 
  width = 6, height = 6
)
plot_BCV_flag(
  obj$dge,
  cex=c(1, 0.1),
  flag=(
    feature_selection$Sig_Any_Interval & 
      feature_selection$gene_id %in% key_genes
  ),
  flag_labels = c("Tagwise - Selected", "Tagwise - Other")
)
dev.off()

## Extract Genes in the top 50th percentile of tag-wise dispersion
genes.disp <- feature_selection %>%
  filter(Tag_disp > quantile(Tag_disp, 0.5)) %>% 
  pull("gene_id")

process_selected_features(
  dge=dge, genes = genes.disp, design=design,
  prefix="GWT_Upper50_Disp"
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
  prefix="GWT_NoBatchSig"
)

## Select genes with an 8 fold injury response at any interval
genes.injury <- feature_selection %>%
  filter(
    (abs(logFC_6Hv0H) > 3 & FDR_6Hv0H < 0.05) |
      (abs(logFC_24Hv0H) > 3 & FDR_24Hv0H < 0.05) |
      (abs(logFC_48Hv0H) > 3 & FDR_48Hv0H < 0.05) |
      (abs(logFC_72Hv0H) > 3 & FDR_72Hv0H < 0.05) |
      (abs(logFC_120Hv0H) > 3 & FDR_120Hv0H < 0.05)
  )%>%
  pull("gene_id")

process_selected_features(
  dge=dge, genes = genes.injury, design=design,
  prefix="GWT_InjurySig"
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
  dge=dge, genes = genes, design=design[,design.exp],
  prefix="GWT_APCombatSeq",
  counts=bat.counts
)

df<-res[[1]]
deg <-res[[2]]
res <- iterate_edgeR_design_coefficients(
  obj.bat[[1]], obj.bat[[2]], coefs=2:6, df=df, design=design[,design.exp],
  deg=deg, prefix="GWT_APCombatSeq", group_label_list = list(
    c("0H", "6H"),c("0H", "24H"),c("0H", "48H"),
    c("0H", "72H"),c("0H", "120H")
  )
)

## Define a table with potential parameters for feature selection
feature_selection_bat <- obj.bat$dge$genes %>% 
  dplyr::select(gene_id, SYMBOL, matches("_FPKM")) %>%
  mutate(
    Tag_disp = obj.bat$dge$tagwise.dispersion,
    AveLogCPM = obj.bat$dge$AveLogCPM,
    lens_dev = ifelse(SYMBOL %in% lens_dev, TRUE, FALSE),
    key_lens_injury = ifelse(SYMBOL %in% key_lens_injury, TRUE, FALSE),
    Trend_disp=obj.bat$dge$trended.dispersion
  ) %>%
  inner_join(
    res[[2]] %>% filter(Samples=="GWT_APCombatSeq") %>%
      mutate(Contrast=paste0(Group_2,"v", Group_1)) %>%
      dplyr::select(gene_id, Contrast, logFC, FDR) %>%
      tidyr::pivot_wider(
        id_cols="gene_id",
        names_from="Contrast",
        values_from=c("logFC", "FDR")
      ) %>% as.data.frame(),
    by="gene_id"
  ) %>%
  mutate(
    Sig_Any_Interval=(
      abs(logFC_6Hv0H) > 1 & FDR_6Hv0H < 0.05)|
      (abs(logFC_24Hv0H) > 1 & FDR_24Hv0H < 0.05)|
      (abs(logFC_48Hv0H) > 1 & FDR_48Hv0H < 0.05)|
      (abs(logFC_72Hv0H) > 1 & FDR_72Hv0H < 0.05)|
      (abs(logFC_120Hv0H) > 1 & FDR_120Hv0H < 0.05)
  )

key_genes <- feature_selection_bat %>%
  filter(key_lens_injury  | lens_dev ) %>%
  pull('gene_id')

cpm_cut <- feature_selection_bat %>%
  filter(gene_id %in% key_genes & Sig_Any_Interval) %>%
  pull("AveLogCPM") %>% min()


## Add selection attribute to feature selection table (Values = TRUE, FALSE)
feature_selection_bat <- feature_selection_bat %>%
  mutate(
    Selected = ifelse(
      (AveLogCPM > cpm_cut *cpm_weight) & 
        (Tag_disp > trend_disp_weight * Trend_disp),
      TRUE, FALSE
    )
  )

write.csv(
  feature_selection,
  file=paste0(
    "LIRTS_DEG_Analysis_results/",
    "GWT_APComBatSeq_Feature_Selection_Table.csv"
  )
)


## Tabulate total number of DEG selected by thresholding on a multiple of the
## trended dispersion (Tagwise > 1.7 X Trended)
write.csv(
  bind_rows(
    feature_selection %>%
      summarise(
        Total_Selected = sum(Selected),
        Total_Significant_DE = sum(
          Selected & Sig_Any_Interval
        )
      ) %>%
      mutate(
        Correction="No Batch Correction"
      ),
    feature_selection_bat %>%
      summarise(
        Total_Selected = sum(Selected),
        Total_Significant_DE = sum(
          Selected & Sig_Any_Interval
        )
      ) %>%
      mutate(
        Correction="Batch Corrected"
      ),
    data.frame(
      Total_Selected=length(
        intersect(
          feature_selection %>% filter(Selected) %>% pull("gene_id"),
          feature_selection_bat %>% filter(Selected) %>% pull("gene_id")
        )
      ),
      Total_Significant_DE=length(
        intersect(
          feature_selection %>% filter(Selected & Sig_Any_Interval)%>% 
            pull("gene_id"),
          feature_selection_bat %>% filter(Selected & Sig_Any_Interval) %>% 
            pull("gene_id")
        )
      ),
      Correction="In Both"
    )
  ),
  file="LIRTS_DEG_Analysis_results/GWT_Selected_Features_Summary.csv"
)

process_selected_features(
  dge, genes = feature_selection %>%filter(Selected) %>%pull("gene_id"),
  design=design, prefix = "GWT_APComBatSeq_Selected",
  counts=bat.counts
)

## Highlight key genes with significant differential expression (lfc > 1)
pdf(
  file=paste0(
    "LIRTS_DEG_Analysis_results/",
    "GWT_APComBatSeq_LFC1_Significant_Key_Genes_BCV_Plot.pdf"
  ), 
  width = 6, height = 6
)
plot_BCV_flag(
  obj.bat$dge,
  cex=c(1, 0.1),
  flag=(
    feature_selection_bat$Sig_Any_Interval & 
      feature_selection_bat$gene_id %in% key_genes
  ),
  flag_labels = c("Tagwise - Selected", "Tagwise - Other")
)
dev.off()

# Store Tagwise dispersions for batch corrected counts in 'twd'
twd <- obj.bat$dge$tagwise.dispersion
genes.bat.disp <- obj.bat$dge$genes$gene_id[twd > quantile(twd, 0.5)]
process_selected_features(
  dge=dge, genes = genes.bat.disp, design=design,
  prefix="GWT_CombatSeq_Upper50_Disp",
  counts=bat.counts
)

### Evaluate DNA-Link Only samples
dge<-master[,
            master$samples$genotype == 'WT' &                # Select Samples
              master$samples$batch %in% c('DNA1', 'DNA2', 'DNA3')
            ]
dge$samples$hours_pcs<-droplevels(dge$samples$hours_pcs)
dge$samples$batch<-droplevels(factor(dge$samples$batch))
dge$samples$genotype<-droplevels(factor(dge$samples$genotype))

# Define Experimental Design: Interval, Batch, and Batch Interaction.
design<-model.matrix(~hours_pcs+batch, dge$samples)     

# Normalize DGEList and fit model
obj <- process_edgeR_ByDesign(
  dge, design=design
)

# Create Diagnostic Figures
diagnostic_plots(
  dge=obj$dge, prefix = "DNA_LinkSamples",
  color_attrib = "hours_pcs", shape_attrib = "batch"
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
  dge=obj$dge, prefix = "Fn1_Groupwise",
  color_attrib = "hours_pcs", shape_attrib = "genotype"
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
  dge=obj$dge, prefix = "Itgb8_Groupwise",
  color_attrib = "hours_pcs", shape_attrib = "genotype"
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

###################### Push Analysis Results to Synapse ######################
# Main Analysis Script
syn_main_script <- File(
  path="transcriptomic_analysis_scripts/LIRTS_DEG_Analysis.R",
  parent=syn_code_dir
)
syn_main_script <- synStore(
  syn_main_script
)

# Wrappers for DEG Analysis Functions
syn_wrap_script <- File(
  path="transcriptomic_analysis_scripts/LIRTS_Wrap_DEG_Functions.R",
  parent=syn_code_dir
)
syn_wrap_script <- synStore(
  syn_wrap_script
)

# Helper script for deg summaries
syn_ovl_script <- File(
  path="transcriptomic_analysis_scripts/Overlap_Comparison_Functions.R",
  parent=syn_code_dir
)
syn_ovl_script <- synStore(
  syn_ovl_script
)

## Construct activity to set provenance on
syn_act <- Activity(
  name="DEG_Analysis",
  description=paste(
    "Analyze various contrasts in the LIRTS",
    "data set for differential expression"
  ),
  executed=c(
    syn_main_script,
    syn_wrap_script,
    syn_ovl_script
  ),
  used=syn_dge
)

# Add a LIRTS_DEG_Tables Directory to store individual DEG Tables
syn_fpkm_dir <- Folder(name="LIRTS_Expression_Matrices", parent=syn_project)
syn_fpkm_dir <- synStore(
  syn_fpkm_dir
)

for(f in list.files("LIRTS_DEG_Analysis_results/", pattern="FPKM_Matrix")){
  syn_fpkm <- File(
    path=paste0(
      "LIRTS_DEG_Analysis_results/", f
    ),
    parent=syn_fpkm_dir
  )
  synStore(syn_fpkm, activity = syn_act)
}

# Add Feature Selection tables to the Genes folder
for(
  f in list.files(
    "LIRTS_DEG_Analysis_results/", 
    pattern="Feature_Selection_Table")
){
  syn_fst <- File(
    path=paste0(
      "LIRTS_DEG_Analysis_results/", f
    ),
    parent=syn_gene_meta
  )
  synStore(syn_fst, activity = syn_act)
}

# Add a LIRTS_DEG_Tables Directory to store individual DEG Tables
syn_deg_dir <- Folder(name="LIRTS_DEG_Tables", parent=syn_project)
syn_deg_dir <- synStore(
  syn_deg_dir
)

for(f in list.files("LIRTS_DEG_Analysis_results/", pattern="Test_DEG")){
  syn_deg <- File(
    path=paste0(
      "LIRTS_DEG_Analysis_results/", f
    ),
    parent=syn_deg_dir
  )
  synStore(syn_deg, activity = syn_act)
}

# Add a LIRTS_DEG_Tables Directory to store Princapal Components and
# BCV Plots (Diagnostic Plots)
syn_fig_dir <- Folder(name="LIRTS_DEG_PCA_BCV_Plots", parent=syn_project)
syn_fig_dir <- synStore(
  syn_fig_dir
)

for(f in list.files("LIRTS_DEG_Analysis_results/", pattern="PCA_Plot")){
  syn_pca <- File(
    path=paste0(
      "LIRTS_DEG_Analysis_results/", f
    ),
    parent=syn_fig_dir
  )
  synStore(syn_pca, activity = syn_act)
}

for(f in list.files("LIRTS_DEG_Analysis_results/", pattern="BCV_Plot")){
  syn_pca <- File(
    path=paste0(
      "LIRTS_DEG_Analysis_results/", f
    ),
    parent=syn_fig_dir
  )
  synStore(syn_pca, activity = syn_act)
}

## Create master DEG table in Synapse (if none exists)
syn_deg_table <- synFindEntityId(
  "LIRTS_DEG_Master_Table", 
  parent=syn_project
)

if(is.null(syn_deg_table)){
  syn_cols <- list(                                 
    Column(name="gene_id", columnType='STRING', maximumSize =100), 
    Column(name='logFC', columnType='DOUBLE'),
    Column(name='logCPM', columnType='DOUBLE'),
    Column(name='PValue', columnType='DOUBLE'),
    Column(name='FDR', columnType='DOUBLE'),
    Column(name='Avg1', columnType='DOUBLE'),
    Column(name='Avg2', columnType='DOUBLE'),
    Column(name='Group_1', columnType='STRING', maximumSize=50),
    Column(name='Group_2', columnType='STRING', maximumSize=50),
    Column(name='Test', columnType='STRING', maximumSize=100),
    Column(name='Samples', columnType='STRING', maximumSize=100)
  )
  syn_schema <- Schema(
    name="LIRTS_DEG_Master_Table",
    columns=syn_cols,
    parent=syn_project
  )
  syn_table <- Table(
    schema = syn_schema,
    values=res[[2]]
  )
  
  syn_table <- synStore(syn_table, activity = syn_act)
}