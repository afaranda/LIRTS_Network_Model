##############################################################################
# File: LTS_DEG_Analysis.R                                                   #
# Purpose:                                                                   #
#       Analyze the LEC Injury time series for differential gene expression. #
#       Use the edgeR exactTest and Quasi Liklihood testing frameowrks       #
#       Generate results for several different cuts of the data.             #
#           - DBI Experiment: 3 contrasts, 24v0, 48v0, 48v24                 #
#           - DNA Link Experiment 1: 3 contrasts, 6v0, 24v0, 24v6            #
#           - DNA Link Experiment 2: 1 contrast, 120v0                       #
#           - DNA Link Experiment 3: 1 contrast, 72vs0                       #
#           - 0 and 24 Hours: 6 contrasts permuted from 24v0, DNALink v DBI  #
#           - Global Wildtype: All subseqeuent intervals vs 0                #
#           - Fibronectin: DBI samples, WT vs FNcKO, 0 vs 48 Hours           #
#           - Beta 8 Integrin: DNA1 samples, WT vs B8cKO, 0 vs 24 Hours      #
#                                                                            #
#       For each of the above experiments generate a DEG Table, a summary    #
#       table, BCV, PCA and MDS plots.                                       #
#                                                                            #
#       This script also generates several additional data sets including:   #
#           - matrices (sample X gene) of raw counts (All genes)             #
#           - matrices of FPKM values (All genes)                            #
#           - matrices of log 2 transformed, TMM normalized CPM values       #
#           - Sample and Gene metadata tables                                #
#           - The 'master' DGEList object                                    #
#                                                                            #
#       TMM normalized CPM matrices were generated for the Wildtype subset   #
#       as well as the full data set. Only genes that passed through         #
#       'filterByExpression' were included in these matrices                 #
#                                                                            #
#                                                                            #
#                                                                            #
# Created: June 11, 2019                                                     #
# Updates:                                                                   #
#     Aug 28, 2020 Section: Global Analysis -- Wildtype Time Series          #
#     Added generation of time series plots and summary statistics for       #
#     selected genes.  these files are stored in the folder                  #
#     "wildtype_time_series_plots"                                           #
#                                                                            #
#     Aug 31, 2020 All Sections                                              #
#     Major change in FPKM Calculations. Previously, FPKM were calculated    #
#     using total library size based on all GTF Features. In this analysis   #
#     FPKM are calculated AFTER filterByExpression has been used to remove   #
#     features below the limit of detection the majority of samples          #
#                                                                            #
#     October 28, 2020 Section: Global Analysis -- Wildype Time Series       #
#     Added several genes to the selection list for time series plots        #
#                                                                            #
#     October 28, 2020 All Sections                                          #
#     Added code to generate length bias plots for all pairwise contrasts.   #
#     The exact test fold changes were used for all contrasts except those   # 
#     estimated in the section "Global Analysis Wildtype time series", where #
#     the quasi likelihood test was used instead.                            #
#                                                                            #
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
#     TODO:                                                                  #
#     Create separate script that wraps edgeR calling sequences and          #
#     plotting functions into functions that take a dgelist and key sample   #
#     information to generate an analysis.                                   #
#                                                                            #
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
library(org.Mm.eg.db)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)
library(synapser)
wd<-getwd()

synapser::synLogin()
syn_project <- "syn25579691"                 ## Synapse ID for this project
syn_count_dir <- "syn25979754"               ## Synapse folder with counts
syn_code_dir <- "syn25976329"                ## Synapse folder with code
syn_sample_table <- "syn25582010"            ## Synapse sample table
syn_gene_meta <- "syn25976328"               ## Gene Annotations
local_data_dir <- paste0(wd,"/LIRTS_Raw_Data")        ## Local data directory

source('transcriptomic_analysis_scripts/PreprocessingFunctions.R')
source('transcriptomic_analysis_scripts/PrincipalComponents.R')
source('transcriptomic_analysis_scripts/ClusteringFunctions.R')
source('transcriptomic_analysis_scripts/Overlap_Comparison_Functions.R')
source('transcriptomic_analysis_scripts/LIRTS_Wrap_DEG_Functions.R')


# Create directory to store results (if it does not yet exist)
if(!dir.exists('LIRTS_DEG_Analysis_results'))
{
  dir.create('LIRTS_DEG_Analysis_results')
}

########################## Define Helper Functions ###########################

# Normalize and fit a dgelist based on a design & contrast matrix

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
dge$samples$genotype<-droplevels(dge$samples$genotype)

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

write.csv(
  res[[1]],"LIRTS_DEG_Analysis_results/DBI_Experiment_DEG_Summary_Tables.csv"
)


# DNA Link Analysis -- First Experiment ####
dge<-master[,
            master$samples$genotype == 'WT' &    # Select Samples
              master$samples$batch == 'DNA1' 
            ]
dge$samples$hours_pcs<-droplevels(dge$samples$hours_pcs)

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

## Iterate over contrasts and prepare summary figures / tables
df<-data.frame()
deg <- res[[2]]
res <- iterate_edgeR_pairwise_contrasts(
  obj[[1]], obj[[2]], cntmat, df=df, design=design,
  deg=deg, prefix="DNA1_Wildtype"
)

write.csv(
  res[[1]],
  "LIRTS_DEG_Analysis_results/DNA_Link_Experiment_1_DEG_Summary_Tables.csv"
)


# DNA Link Analysis --  Second Experiment ####
dge<-master[,
            master$samples$genotype == 'WT' &    # Select Samples
              master$samples$batch == 'DNA2' 
            ]
dge$samples$hours_pcs<-droplevels(dge$samples$hours_pcs)


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
df<-data.frame()
deg <- res[[2]]
res <- iterate_edgeR_pairwise_contrasts(
  obj[[1]], obj[[2]], cntmat, df=df, design=design,
  deg=deg, prefix="DNA2_Wildtype"
)

write.csv(
  res[[1]],
  "LIRTS_DEG_Analysis_results/DNA_Link_Experiment_2_DEG_Summary_Tables.csv"
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
  dge=obj$dge, prefix = "24_Hour_LabComp",
  color_attrib = "hours_pcs", shape_attrib = "batch"
)

## Iterate over contrasts and prepare summary figures / tables
df<-data.frame()
deg <- res[[2]]
res <- iterate_edgeR_design_coefficients(
  obj[[1]], obj[[2]], coefs=c(2,3,4), df=df, design=design,
  deg=deg, prefix="24_Hour_LabComp", group_label_list = list(
    c("0H", "24H"), c("DBI", "DNA1"), c("Intercept", "InxDNA1_24H")
  )
)

write.csv(
  res[[1]],
  paste0(
    "LIRTS_DEG_Analysis_results/",
    "24vs0_Hour_Lab_Comparison_DEG_Summary_Tables.csv"
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

# Define Experimental Design (full model)
design<-model.matrix(~batch + hours_pcs, dge$samples)   
#colnames(design)<-gsub("hours_pcs", '', colnames(design))

obj <- process_edgeR_ByDesign(
  dge, genes = genes,
  design = design
)

# Iterate over model coefficients and generate DEG tables
df<-data.frame()
deg <- data.frame()
res <- iterate_edgeR_design_coefficients(
  obj[[1]], obj[[2]], coefs=5:9, df=df, design=design,
  deg=deg, prefix="Global_Wildtype", group_label_list = list(
    c("0H", "6H"), c("0H", "24H"),c("0H", "48H"), c("0H", "72H"),
    c("0H", "120H")
  )
)

# Identify genes with a significant batch dependent fold change
# between DBI samples and the three DNA Link Batches
bat<-glmQLFTest(obj$fit, coef=2:4)          # Estimate Batch DE
bdg<-as.data.frame(topTags(bat, n=Inf)) %>%
  dplyr::select(gene_id, logCPM, PValue, FDR, matches("logFC"))



# Get FPKM Matrix for individual samples
fpkm<-as.data.frame(rpkm(dge, gene.length="eu_length"))
colnames(fpkm) <- obj$dge$samples$label
fpkm$gene_id<-row.names(fpkm)

makeContrasts(
  DBIvsDNA1=Intercept - batchDNA1,
  levels=design
)

## Generate diagnostic plots. 

## Iterate over contrasts and prepare summary figures / tables
df<-data.frame()
coefs=c(WT6='6H', WT24='24H', WT48='48H', WT72='72H', WT120='120H')

## Generate Time Series Plots for Selected Genes
#List of Genes over which to tabulate time series
# Un-comment genes that you would like to generate a table for
genes<-c(
  Egr1="ENSMUSG00000038418", Egr2="ENSMUSG00000037868", Egr3="ENSMUSG00000033730", 
  Grem1="ENSMUSG00000074934", Fn1="ENSMUSG00000026193",
  Actn1="ENSMUSG00000015143", Acta2="ENSMUSG00000035783", Col1a1="ENSMUSG00000001506",
  Atf3="ENSMUSG00000026628", Pitx3="ENSMUSG00000025229", Klf2="ENSMUSG00000055148",
  Emp1="ENSMUSG00000030208",  
  #Ppbp="ENSMUSG00000029372", #Cxcl9="ENSMUSG00000029417",
  Cxcr2="ENSMUSG00000026180", #Cxcr4="ENSG00000121966", Runx1="ENSMUSG00000022952",
  Tnc="ENSMUSG00000028364", Fos="ENSMUSG00000021250", Ier2="ENSMUSG00000053560",
  Fosb="ENSMUSG00000003545", Jun="ENSMUSG00000052684", Junb="ENSMUSG00000052837",
  Timp1 = 'ENSMUSG00000001131', Serpine1 = 'ENSMUSG00000037411', Hmga1 = 'ENSMUSG00000046711',
  Thbs1 ='ENSMUSG00000040152', Ier3 = 'ENSMUSG00000003541', Esm1 = 'ENSMUSG00000042379',
  Zc3h12a = 'ENSMUSG00000042677', Dusp6 = 'ENSMUSG00000019960', Tcim  = 'ENSMUSG00000056313',
  Mrtfa = 'ENSMUSG00000042292', Ch25h = 'ENSMUSG00000050370', Tagln = 'ENSMUSG00000032085',
  Tagln2 = 'ENSMUSG00000026547', Hsbp1 = 'ENSMUSG00000031839',
  Dkk3 = 'ENSMUSG00000030772', Ccl2 = "ENSMUSG00000035385", 
  Jund = 'ENSMUSG00000071076', Cdk5r1='ENSMUSG00000048895', F3='ENSMUSG00000028128',
  Fas = 'ENSMUSG00000024778', Gadd45a='ENSMUSG00000036390', Gadd45b='ENSMUSG00000015312',
  Icam1 = 'ENSMUSG00000037405', Nr4a1='ENSMUSG00000023034', Rcan1='ENSMUSG00000022951', 
  Runx1 = 'ENSMUSG00000022952', Itgb8="ENSMUSG00000025321", Cav1 ="ENSMUSG00000007655",
  Col1a1 = "ENSMUSG00000001506", Ecm1="ENSMUSG00000028108", Ltbp2="ENSMUSG00000002020",
  Vcan = "ENSMUSG00000021614", Nfkb1="ENSMUSG00000028163",Bcl3="ENSMUSG00000053175", 
  Ldb1="ENSMUSG00000025223", Spi1="ENSMUSG00000002111", Smad6="ENSMUSG00000036867",
  Igfbp3="ENSMUSG00000020427", Cebpa="ENSMUSG00000034957", Fbln2="ENSMUSG00000064080",
  Col5a3="ENSMUSG00000004098", Col5a2="ENSMUSG00000026042", Col7a1="ENSMUSG00000025650",
  Col16a1="ENSMUSG00000040690", Atf4="ENSMUSG00000042406", Serpina3n="ENSMUSG00000021091",
  Serpina3h="ENSMUSG00000041449", Tdrd7="ENSMUSG00000035517", Celf1="ENSMUSG00000005506", 
  Rbm3="ENSMUSG00000031167", Ptx3="ENSMUSG00000027832", Cxcl1="ENSMUSG00000029380",
  Cxcl2="ENSMUSG00000058427", Cxcl5="ENSMUSG00000029371", Ppbp="ENSMUSG00000029372",
  Tagln3 = 'ENSMUSG00000022658', Wnt16 = 'ENSMUSG00000029671', Serpina3i="ENSMUSG00000079014"
)

#  This gene was observed with counts from all alignable reads, it's
#   counts plummet after ribofiltering Serpina3g="ENSMUSG00000041481"

#  These genes were lost after excluding features that were not measured 
#  in all three batches: 
# lost_genes <- c(Ptx3="ENSMUSG00000027832", Cxcl1="ENSMUSG00000029380", Cxcl2="ENSMUSG00000058427",
#     Cxcl5="ENSMUSG00000029371", Ppbp="ENSMUSG00000029372", Tagln3 = 'ENSMUSG00000022658',
#     Wnt16 = 'ENSMUSG00000029671', Serpina3c="ENSMUSG00000066361", Serpina3f="ENSMUSG00000066363",
#     Serpina3m="ENSMUSG00000079012", Serpina3i="ENSMUSG00000079014"
# )
## After switching from batch to library type, most of these came back except for:
## Serpina3c, Serpina3f, and Serpina3m


# Check entrez ID's against symbols
for(g in names(genes)){
  eid<- master$genes %>% filter(SYMBOL == g) %>% pull(gene_id)
  if(eid != genes[g]){
    print(paste(g, eid, genes[g]))
  } else {
    print("valid")
  }
}

# Plot RPKM for EGR1, FosB, IER2, JunB, and Fos/c-Fos
for (g in names(genes)){
  plt<-inner_join(
    dge$samples, 
    reshape2::melt(fpkm %>% filter(gene_id == genes[g])) %>%
      filter(variable != "eu_length") %>%
      dplyr::select(label=variable, FPKM=value),
    by="label"
  )
  plt$gene_symbol<-g
  plt$ensembl_id<-genes[g]
  ggplot(plt, aes(x=hours_pcs, y=FPKM, group=hours_pcs)) + 
    geom_boxplot() + ggtitle(g)
  ggsave(
    filename=paste(
      "LIRTS_DEG_Analysis_results/",g,
      "_fpkm_box.png", sep=""
    ), width=5, height=3
  )
  
  ggplot(plt, aes(x=hours_pcs, y=FPKM, group=hours_pcs)) + 
    geom_point(aes(color=batch)) + ggtitle(g)
  ggsave(
    filename=paste(
      "LIRTS_DEG_Analysis_results/",g,
      "_fpkm_pnt.png", sep=""
    ), width=5, height=3
  )
  
  fn<- paste("LIRTS_DEG_Analysis_results/",g,"_FPKM_Time.tsv", sep="")
  write.table(
    plt %>%
      arrange(interval, batch) %>%
      dplyr::select(
        gene_symbol, ensembl_id, 
        sample, batch, hours.pcs, 
        interval, FPKM
      ), file = fn, row.names = F,
    sep="\t", quote=F
  )
}


# Generate Tables For Further Plotting, with statistics
write.csv(df,"LIRTS_DEG_Analysis_results/Global_Wildtype_DEG_Summary_Tables.csv")
write.csv(deg,"LIRTS_DEG_Analysis_results/Global_Wildtype_QLFTest_DEG.csv")
write.csv(
  fpkm, 
  "LIRTS_DEG_Analysis_results/Global_Wildtype_QLFTest_FPKM_Matrix.csv"
)

# edgeR TMM/cpm (log2 transformed) used for clustering (whole data set)
c<-cpm(dge, log=T)
write.csv(
  c, 
  "LIRTS_DEG_Analysis_results/Global_Wildtype_TMM_Normalized_Present_Genes.txt", 
  row.names=F
)

### Plot PCA for filtered data

## Filtered by top dispersion:
dge$genes$twd <-(
  dge$tagwise.dispersion> quantile(dge$tagwise.dispersion, 0.5)
)
pca <- prcomp(
  t(cpm(dge[dge$genes$twd,], log=T)),scale=T
)
ggsave(
  paste0(
    "LIRTS_DEG_Analysis_results/",
    "Global_Wildtype_Above_Median_Dispersion_PCA_Plot.png"
  ),
  autoplot(
    pca, data=dge$samples,
    colour="interval", shape="batch", size=3
  ), height =4, width=6
)

## Filtered by top Variance:
dge$genes$variance <- apply(cpm(dge), 1, var)
byvar <- dge$genes$variance > quantile(dge$genes$variance, 0.5)
pca <- prcomp(
  t(cpm(dge[byvar,], log=T)),scale=T
)
ggsave(
  paste0(
    "LIRTS_DEG_Analysis_results/",
    "Global_Wildtype_Above_Median_Variance_PCA_Plot.png"
  ),
  autoplot(
    pca, data=dge$samples,
    colour="interval", shape="batch", size=3
  ), height =4, width=6
)


## Exclude features with strong batch dependence
batch_genes <- bdg %>%
  filter( FDR < 0.05) %>% 
  filter(abs(logFC.batchDNA1)>1 | abs(logFC.batchDNA2)> 1) %>%
  pull(gene_id)
dgf <- dge[dge$genes$twd,]
pca <- prcomp(
  t(cpm(dgf[setdiff(dgf$genes$gene_id, batch_genes),], log=T)),scale=T
)
ggsave(
  paste0(
    "LIRTS_DEG_Analysis_results/",
    "Global_Wildtype_Exclude_Batch_Sensetive_PCA_Plot.png"
  ),
  autoplot(
    pca, data=dge$samples,
    colour="interval", shape="batch", size=3
  ), height =4, width=6
)



## Exclude features with strong batch dependence and only consider genes with
## two fold differential expression between 0 hours and any one interval
injury_genes <- deg %>%
  filter( FDR < 0.05) %>% 
  filter(
    abs(logFC.6H )>1 |
      abs(logFC.24H)> 1 |
      abs(logFC.48H)> 1 |
      abs(logFC.120H)>1
  ) %>%
  pull(gene_id)

pca <- prcomp(
  t(cpm(dge[setdiff(injury_genes, batch_genes),], log=T)),scale=T
)
ggsave(
  paste0(
    "LIRTS_DEG_Analysis_results/",
    "Global_Wildtype_Include_Injury_Responsive_PCA_Plot.png"
  ),
  autoplot(
    pca, data=dge$samples,
    colour="interval", shape="batch", size=3
  ), height =4, width=6
)

## Show Correlation Structure for all detectable genes
cm <- cor(cpm(dge[,], log=T), method = "spearman")
pheatmap(
  cm, annotation_col = dge$samples[,c("interval", "batch")],
  filename = paste0(
    "LIRTS_DEG_Analysis_results/","
    Global_Wildtype_cormat.png"
  ),height = 6, width=8
)

## Show Correlation Structure for genes above median dispersion
cm <- cor(cpm(dge[dge$genes$twd,], log=T), method = "spearman")
pheatmap(
  cm, annotation_col = dge$samples[,c("interval", "batch")],
  filename = paste0(
    "LIRTS_DEG_Analysis_results/","
    Global_Wildtype_Above_Median_Dispersion_cormat.png"
  ),height = 6, width=8
)

## Show Correlation Structure for genes above median variance
cm <- cor(cpm(dge[byvar,], log=T), method = "spearman")
pheatmap(
  cm, annotation_col = dge$samples[,c("interval", "batch")],
  filename = paste0(
    "LIRTS_DEG_Analysis_results/","
    Global_Wildtype_Above_Median_Variance_cormat.png"
  ),height = 6, width=8
)


## Show Correlation Structure after excluding batch sensetive genes
cm <- cor(
  cpm(dge[setdiff(dge$genes$gene_id, batch_genes),], log=T), 
  method = "spearman"
)
pheatmap(
  cm, annotation_col = dge$samples[,c("interval", "batch")],
  filename = paste0(
    "LIRTS_DEG_Analysis_results/",
    "Global_Wildtype_Exclude_Batch_Sensetive_cormat.png"
  ),height = 6, width=8
)

## Show Correlation Structure after excluding batch sensetive genes and
## restricting analysis to top injury response genes
cm <- cor(
  cpm(dge[setdiff(injury_genes, batch_genes),], log=T), 
  method = "spearman"
)
pheatmap(
  cm, annotation_col = dge$samples[,c("interval", "batch")],
  filename = paste0(
    "LIRTS_DEG_Analysis_results/",
    "Global_Wildtype_Include_Injury_Responsive_cormat.png"
  ),height = 6, width=8
)

# edgeR TMM/cpm (batch corrected) used for clustering (whole data set)
design.exp<-design[,setdiff(colnames(design), c('batchDNA1', 'batchDNA2'))]
c<-cpm(dge[setdiff(dge$genes$gene_id, batch_genes),], log=T)
c.bat<-removeBatchEffect(
  c, design = design.exp, batch=dge$samples$batch
)
pca <- prcomp(
  t(c.bat),scale=T
)
ggsave(
  paste0(
    "LIRTS_DEG_Analysis_results/",
    "Global_Wildtype_Batch_Corrected_PCA_Plot.png"
  ),
  autoplot(
    pca, data=dge$samples,
    colour="interval", shape="batch", size=3
  ), height =4, width=6
)

write.csv(
  c.bat, 
  "LIRTS_DEG_Analysis_results/Global_Wildtype_TMM_Batch_Corrected_Present_Genes.txt",
  row.names=T
)

cm <- cor(c.bat, method = "spearman")
pheatmap(
  cm, annotation_col = dge$samples[,c("interval", "batch")],
  filename = paste0(
    "LIRTS_DEG_Analysis_results/",
    "Global_Wildtype_Batch_Corrected_cormat.png"
  ), height = 6, width=8
)

######################## Analyze Mutation effects ############################
# DBI -- Fibronectin vs Wildtype 0 and 48 hours ####
dge<-master[,
              master$samples$hours_pcs %in% c('0H', '48H')&    # Select Samples
              master$samples$batch == 'DBI',
            ]
dge$samples$hours_pcs<-droplevels(dge$samples$hours_pcs)
dge$samples$batch<-droplevels(factor(dge$samples$batch))

design<-model.matrix(~0+group, dge$samples)     # Define Experimental Design
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

dge<-dge[filterByExpr(dge, design), ,keep.lib.sizes=F] # drop low features
dge<-calcNormFactors(dge)                    # Factors & Dispersion
dge<-estimateDisp(dge, design, robust = T)

rbg<-as.data.frame(                            # Add Group Mean FPKMs
  rpkmByGroup(dge, gene.length = "eu_length")
) 
rbg$gene_id<-row.names(rbg)
dge$genes<-merge(
  dge$genes, rbg,
  by='gene_id'
)
row.names(dge$genes)<-dge$genes$gene_id

fit<-glmQLFit(dge, design, robust = T)  # Run Model and estimate DE
qlf<-glmQLFTest(fit, contrast = cntmat)
deg<-as.data.frame(topTags(qlf, n=Inf))

## Generate diagnostic plots. 
png("LIRTS_DEG_Analysis_results/Fibronectin_Experiment_BCV_Plot.png")
plotBCV(dge)                                                 # BCV Plot
dev.off()

png(
  "LIRTS_DEG_Analysis_results/Fibronectin_Experiment_PCA_Plot.png",  # PCA Plot
  width=600, height = 400
)
plotPrinComp(
  cpm(dge, log=T), ft=dge$samples,                  
  idCol = 0, groupCol = 'group'
)
dev.off()

png(
  "LIRTS_DEG_Analysis_results/Fibronectin_Experiment_MDS_Plot.png",  # MDS Plot
  width=600, height=500
)
par(mar=c(8.5,5,4.1,2.1))
plotMDS(
  dge, cex=1.5, cex.lab = 1.5,
  cex.axis = 1.5, cex.main = 2,
  pch = rep(15:18, each=3),
  main = "Fibronecitn Experiment",
  col = rep(c('red', 'blue','green4','black'), each=3)
)
legend(
  -3, -2.4, 
  legend = unique(dge$samples$group), 
  pch = 15:18, cex = 1.5,
  col=c('red', 'blue','green4', 'black'), xpd=T
)
dev.off()

## Iterate over contrasts and prepare summary figures / tables
df<-data.frame()
for (c in colnames(cntmat)) {             
  # Run Exact Tests
  pair<-c(
    names(cntmat[,c])[cntmat[,c] == -1], 
    names(cntmat[,c])[cntmat[,c] == 1]
  )
  deg.et<-as.data.frame(topTags(exactTest(dge, pair = pair), n=Inf))
  deg.qt<-as.data.frame(topTags(glmQLFTest(fit,contrast=cntmat[,c]), n=Inf))
  
  fn<-paste(
    "LIRTS_DEG_Analysis_results/Fibronectin_Experiment_",
    c,"_Exact_Test_DEG.csv", sep="")
  write.csv(deg.et, fn)
  
  fn<-paste(
    "LIRTS_DEG_Analysis_results/Fibronectin_Experiment_",
    c,"_QLFTest_DEG.csv", sep="")
  write.csv(deg.qt, fn)
  
  # Generate volcano plots 
  fn<-paste(
    "LIRTS_DEG_Analysis_results/Fibronectin_Experiment_",
    c,"_QLF_Volcano.png", sep="")
  ggsave(
    fn,
    EnhancedVolcano(
      toptable = deg.qt, 
      lab = deg.qt$SYMBOL,
      x = "logFC",
      y="FDR",
      ylim =  c(0, max(-log10(deg.qt[,'FDR']), na.rm=TRUE) + 2),
      title = "Fibronectin Experiment",
      subtitle = c,
      pCutoff = 0.05
    )
  )
  
  fn<-paste(
    "LIRTS_DEG_Analysis_results/Fibronectin_Experiment_",
    c,"_Exact_Volcano.png", sep="")
  ggsave(
    fn,
    EnhancedVolcano(
      toptable = deg.et, 
      lab = deg.et$SYMBOL,
      x ="logFC",
      y="FDR",
      ylim =  c(0, max(-log10(deg.et[,'FDR']), na.rm=TRUE) + 2),
      title = "Fibronectin Experiment",
      subtitle = c,
      pCutoff = 0.05
    )
  )
  
  dg<-degSummary(                         # Generate QLF Test Summary tables 
    deg.qt,
    lfc = "logFC",
    fdr = 'FDR', 
    Avg1 = rownames(cntmat)[cntmat[,c] !=0][1],
    Avg2 = rownames(cntmat)[cntmat[,c] !=0][2]
  )
  dg$contrast<-c
  dg$test<-"QLFTest"
  df<-bind_rows(df, dg)
  print(names(deg.et))
  
  dg<-degSummary(                         # Generate Exact Test Summary tables
    deg.et,
    lfc = 'logFC',
    fdr = 'FDR', 
    Avg1 = rownames(cntmat)[cntmat[,c] !=0][1],
    Avg2 = rownames(cntmat)[cntmat[,c] !=0][2]
  )
  dg$contrast<-c
  dg$test<-"Exact Test"
  df<-bind_rows(df, dg)
  
  fn<-paste(
    "LIRTS_DEG_Analysis_results/Fibronectin_Experiment_",
    c,"_Length_Bias.png", sep=""
  )
  mn<-paste("Length Bias in ",c,sep='')
  yl <-paste(                                    # Construct y axis label
    "Log2 Fold Change in", pair[2],          
    "vs",pair[1]
  )
  png(fn, width=6, height = 5, units="in", res=1200)
  plot(
    x=log(deg.et$eu_length,2), y=deg.et$logFC, main=mn,
    xlab = "Log2 Gene Length (exon-union)",
    ylab = yl,
    col = ifelse(abs(deg.et$logFC) > 1 & deg.et$FDR < 0.05, "red", "black")
  )
  
  # Add Correlation Coefficients to tests
  ct=cor.test(deg.et$logFC, log(deg.et$eu_length,2), method="spearman")
  rho=paste("rho:", round(ct$estimate,3))
  sig=paste("p value:", signif(ct$p.value,3), sep="")
  abline(lm(deg.et$logFC~log(deg.et$eu_length,2)), col="red", lwd=2)
  text(6.5,max(deg.et$logFC)-1,rho)
  text(7,max(deg.et$logFC)-3,sig)
  dev.off()
}
write.csv(
  df,"LIRTS_DEG_Analysis_results/Fibronectin_Experiment_DEG_Summary_Tables.csv")
write.csv(
  deg,
  "LIRTS_DEG_Analysis_results/Fibronectin_Experiment_QLFTest_All_Contrasts_DEG.csv"
)


# DNA Link -- Beta 8 Integrin vs Wildtype 0 and 48 hours ####
dge<-master[,
            master$samples$hours_pcs %in% c('0H', '24H')&    # Select Samples
              master$samples$batch == 'DNA1',
            ]

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

dge<-dge[filterByExpr(dge, design), ,keep.lib.sizes=F] # drop low features
dge<-calcNormFactors(dge)                      # Factors & Dispersion
dge<-estimateDisp(dge, design, robust = T)

rbg<-as.data.frame(                            # Add Group Mean FPKMs
  rpkmByGroup(dge, gene.length = "eu_length")
) 
rbg$gene_id<-row.names(rbg)
dge$genes<-merge(
  dge$genes, rbg,
  by='gene_id'
)
row.names(dge$genes)<-dge$genes$gene_id

fit<-glmQLFit(dge, design, robust = T)  # Run Model and estimate DE
qlf<-glmQLFTest(fit, contrast = cntmat)
deg<-as.data.frame(topTags(qlf, n=Inf))

## Generate diagnostic plots. 
png("LIRTS_DEG_Analysis_results/Beta8_Integrin_Experiment_BCV_Plot.png")
plotBCV(dge)                                                 # BCV Plot
dev.off()

png(
  "LIRTS_DEG_Analysis_results/Beta8_Integrin_Experiment_PCA_Plot.png",# PCA Plot
  width=600, height = 400
)
plotPrinComp(
  cpm(dge, log=T), ft=dge$samples,                  
  idCol = 0, groupCol = 'group'
)
dev.off()

png(
  "LIRTS_DEG_Analysis_results/Beta8_Integrin_Experiment_MDS_Plot.png",# MDS Plot
  width=600, height=500
)
par(mar=c(8.5,5,4.1,2.1))
plotMDS(
  dge, cex=1.5, cex.lab = 1.5,
  cex.axis = 1.5, cex.main = 2,
  pch = rep(15:18, each=3),
  main = "Beta8 Integrin Experiment",
  col = rep(c('red', 'blue', 'green4', 'black'), each=3)
)
legend(
  -3, -2.4, 
  legend = unique(dge$samples$group), 
  pch = 15:17, cex = 1.5,
  col=c('red', 'blue', 'green4', 'black'), xpd=T
)
dev.off()

## Iterate over contrasts and prepare summary figures / tables
df<-data.frame()
for (c in colnames(cntmat)) {             
  # Run Exact Tests
  pair<-c(
    names(cntmat[,c])[cntmat[,c] == -1], 
    names(cntmat[,c])[cntmat[,c] == 1]
  )
  deg.et<-as.data.frame(topTags(exactTest(dge, pair = pair), n=Inf))
  deg.qt<-as.data.frame(topTags(glmQLFTest(fit,contrast=cntmat[,c]), n=Inf))
  
  fn<-paste(
    "LIRTS_DEG_Analysis_results/Beta8_Integrin_Experiment_",
    c,"_Exact_Test_DEG.csv", sep="")
  write.csv(deg.et, fn)
  
  fn<-paste(
    "LIRTS_DEG_Analysis_results/Beta8_Integrin_Experiment_",
    c,"_QLFTest_DEG.csv", sep="")
  write.csv(deg.qt, fn)
  
  # Generate volcano plots 
  fn<-paste(
    "LIRTS_DEG_Analysis_results/Beta8_Integrin_Experiment_",
    c,"_QLF_Volcano.png", sep="")
  ggsave(
    fn,
    EnhancedVolcano(
      toptable = deg.qt, 
      lab = deg.qt$SYMBOL,
      x = "logFC",
      y="FDR",
      ylim =  c(0, max(-log10(deg.qt[,'FDR']), na.rm=TRUE) + 2),
      title = "Beta8 Integrin Experiment",
      subtitle = c,
      pCutoff = 0.05
    )
  )
  
  fn<-paste(
    "LIRTS_DEG_Analysis_results/Beta8_Integrin_Experiment_",
    c,"_Exact_Volcano.png", sep="")
  ggsave(
    fn,
    EnhancedVolcano(
      toptable = deg.et, 
      lab = deg.et$SYMBOL,
      x ="logFC",
      y="FDR",
      ylim =  c(0, max(-log10(deg.et[,'FDR']), na.rm=TRUE) + 2),
      title = "Beta8 Integrin Experiment",
      subtitle = c,
      pCutoff = 0.05
    )
  )
  
  dg<-degSummary(                         # Generate QLF Test Summary tables 
    deg.qt,
    lfc = "logFC",
    fdr = 'FDR', 
    Avg1 = rownames(cntmat)[cntmat[,c] !=0][1],
    Avg2 = rownames(cntmat)[cntmat[,c] !=0][2]
  )
  dg$contrast<-c
  dg$test<-"QLFTest"
  df<-bind_rows(df, dg)
  print(names(deg.et))
  
  dg<-degSummary(                         # Generate Exact Test Summary tables
    deg.et,
    lfc = 'logFC',
    fdr = 'FDR', 
    Avg1 = rownames(cntmat)[cntmat[,c] !=0][1],
    Avg2 = rownames(cntmat)[cntmat[,c] !=0][2]
  )
  dg$contrast<-c
  dg$test<-"Exact Test"
  df<-bind_rows(df, dg)
  
  fn<-paste(
    "LIRTS_DEG_Analysis_results/Beta8_Integrin_Experiment_",
    c,"_Length_Bias.png", sep=""
  )
  mn<-paste("Length Bias in ",c,sep='')
  yl <-paste(                                    # Construct y axis label
    "Log2 Fold Change in", pair[2],          
    "vs",pair[1]
  )
  
  png(fn, width=6, height = 5, units="in", res=1200)
  plot(
    x=log(deg.et$eu_length,2), y=deg.et$logFC, main=mn,
    xlab = "Log2 Gene Length (exon-union)",
    ylab = yl,
    col = ifelse(abs(deg.et$logFC) > 1 & deg.et$FDR < 0.05, "red", "black")
  )
  
  # Add Correlation Coefficients to tests
  ct=cor.test(deg.et$logFC, log(deg.et$eu_length,2), method="spearman")
  rho=paste("rho:", round(ct$estimate,3))
  sig=paste("p value:", signif(ct$p.value,3), sep="")
  abline(lm(deg.et$logFC~log(deg.et$eu_length,2)), col="red", lwd=2)
  text(6.5,max(deg.et$logFC)-1,rho)
  text(7,max(deg.et$logFC)-3,sig)
  dev.off()
}
write.csv(
  df,
  "LIRTS_DEG_Analysis_results/Beta8_Integrin_Experiment_DEG_Summary_Tables.csv")
write.csv(
  deg,
  "LIRTS_DEG_Analysis_results/Beta8_Integrin_Experiment_QLFTest_All_Contrasts_DEG.csv"
)

####### Write Out Full Count and RPKM Tables, and Filtered CPM Tables ########
dge<-master
design<-model.matrix(~0 + genotype:interval + batch, dge$samples)
design<-design[,apply(design, 2, sum) > 0]
design<-design[, -3]                           # Drop redundant DNA2 column
colnames(design)<-c(
  'DBI', 'DNA1',
  'WT_0H', 'FN_0H', 'B8_0H',
  'WT_6H', 'WT_24H', 'B8_24H',
  'WT_48H', 'FN_48H', 'WT_120H'
)
dge<-dge[filterByExpr(dge, design), , keep.lib.size=FALSE]

# Calculate size factors and Dispersion
dge<-calcNormFactors(dge)

# Raw Counts for Complete Expermient
m<-as.data.frame(master$counts)
s<-colnames(m)
m$gene_id<-row.names(m)
m<-m[,c('gene_id', s)]
write.csv(
  m, 
  "LIRTS_DEG_Analysis_results/Injury_Model_Raw_Counts_All_Genes.txt", 
  row.names=F
)

# Length normalized, log2 transformed FPKM values (no TMM scaling)
r<-as.data.frame(rpkm(master, gene.length = "eu_length", log=T))
s<-colnames(r)
r$gene_id<-row.names(r)
r<-r[,c('gene_id', s)]

write.csv(
  r, 
  "LIRTS_DEG_Analysis_results/Full_Injury_Model_FPKM_All_Genes.csv", 
  row.names = F
)

# Sample ID's and Covariates
write.csv(
  master$samples, 
  "LIRTS_DEG_Analysis_results/Full_Injury_Model_Sample_Metadata.txt", 
  row.names=F
)

# Gene ID's with Symbol, name and Length Used for RPKM normalization
write.csv(
  master$genes, 
  "LIRTS_DEG_Analysis_results/Full_Injury_Model_Gene_Metadata.txt", 
  row.names=F
)

# edgeR TMM/cpm (log2 transformed) used for clustering (whole data set)
c<-cpm(dge, log=T)
write.csv(
  c, 
  "LIRTS_DEG_Analysis_results/Full_Injury_Model_TMM_Normalized_Present_Genes.txt", 
  row.names=F
)

png(
  "LIRTS_DEG_Analysis_results/Full_Injury_Model_PCA_Plot_ByInterval.png",     # PCA Plot
  width=600, height = 400
)
plotPrinComp(
  c, ft=dge$samples,                  
  idCol = 0, groupCol = 'group'
)
dev.off()

png(
  "LIRTS_DEG_Analysis_results/Full_Injury_Model_PCA_Plot_ByBatch.png",        # PCA Plot
  width=600, height = 400
)
plotPrinComp(
  c, ft=dge$samples,                  
  idCol = 0, groupCol = 'batch'
)
dev.off()


# edgeR TMM/cpm (batch corrected) used for clustering (whole data set)
design.exp<-design[,setdiff(colnames(design), c('DBI', 'DNA1'))]
c.bat<-removeBatchEffect(
  c, design = design.exp, batch=dge$samples$batch
)

write.csv(
  c.bat, 
  "LIRTS_DEG_Analysis_results/Full_Injury_Model_TMM_Batch_Corrected_Present_Genes.txt",
  row.names=F
)

png(
  "LIRTS_DEG_Analysis_results/Full_Injury_Model_BatchCorrected_PCA_Plot_ByInterval.png",
  width=600, height = 400
)
plotPrinComp(
  c.bat, ft=dge$samples,                  
  idCol = 0, groupCol = 'group'
)
dev.off()

png(
  "LIRTS_DEG_Analysis_results/Full_Injury_Model_BatchCorrected_PCA_Plot_ByBatch.png",
  width=600, height = 400
)
plotPrinComp(
  c.bat, ft=dge$samples,                  
  idCol = 0, groupCol = 'batch'
)
dev.off()



