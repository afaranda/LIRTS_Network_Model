##############################################################################
# File: LTS_DEG_Analysis.R                                                   #
# Purpose:                                                                   #
#       Analyze the LEC Injury time series for differential gene expression. #
#       Use the edgeR exactTest and Quasi Liklihood testing frameowrks       #
#       Generate results for several different cuts of the data.             #
#           - DBI Experiment: 3 contrasts, 24v0, 48v0, 48v24                 #
#           - DNA Link Experiment 1: 3 contrasts, 6v0, 24v0, 24v6            #
#           - DNA Link Experiment 2: 1 contrast, 120v0                       #
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
#                                                                            #
#     TODO                                                                   #
#        -  When filtering by expression, take DNA1 and DNA2 together.       #
#           Instead of batch, use sequeuncing technology as the determining  #
#           factor                                                           #
#                                                                            #                                      
# Author: Adam Faranda                                                       #
##############################################################################

############################ Setup Environment ###############################

# CHANGE THIS DIRECTORY
#setwd('~/Documents/LEC_Time_Series')
setwd("~/Documents/Adam_LEC_Time_Series_DEG_Analysis_30_Apr_2021")
#load("GeneLengthTable.Rdata")
library(dplyr)
library(cluster)
library(reshape2)
library(org.Mm.eg.db)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)
wd<-getwd()

source('transcriptomic_analysis_scripts/BuildDataMatrix.R')
source('transcriptomic_analysis_scripts/PreprocessingFunctions.R')
source('transcriptomic_analysis_scripts/PrincipalComponents.R')
source('transcriptomic_analysis_scripts/ClusteringFunctions.R')
source('transcriptomic_analysis_scripts/Overlap_Comparison_Functions.R')
#source('transcriptomic_analysis_scripts/Wrap_edgeR_Functions.R')


# Create directory to store results (if it does not yet exist)
if(!dir.exists('LTS_DEG_Analysis_results'))
{
  dir.create('LTS_DEG_Analysis_results')
}
############################ Load in Data Files ##############################

# CHANGE THIS DIRECTORY
dl<-paste(
  wd, "/LEC_Time_Series_HTSeq_Counts",
  sep = ""
)
ft<-hc_getFileTable(
  dirList=dl, filename = "HTSeq_GeneCounts_All.csv"
)

ds<-hc_loadFiles(ft)
ft<-hc_identifierConsistency(ds, ft)

###################### Annotate Genes in table 'lt' ##########################

# Import Annotations, append to "Gene Length" Table
fn<-paste(dl,'Gene_Annotations.csv', sep='/')
if(file.exists(fn)){
  lt<-read.csv(fn, row.names = 1)
  print(TRUE)
} else {
  library('AnnotationHub')
  lt<-read.table(
    paste(dl,'Ens_Mm101_Length_GC.txt', sep="/"),
    header=T, quote="", sep="\t", 
    stringsAsFactors = F
  )
  ah<-AnnotationHub()
  # Run Query to find proper Annotation Set:
  #AnnotationHub::query(ah, pattern=c("EnsDb", "Mus musculus", "101"))
  edb<-ah[['AH83247']]
  lt<-merge(
    lt, AnnotationDbi::select(
      edb, keys=lt$gene_id, 
      columns = c("SYMBOL", "DESCRIPTION", "GENEBIOTYPE", "SEQNAME"), 
      keytype = "GENEID"
    ),
    by.x = 'gene_id', by.y='GENEID'
  )
  
  library(org.Mm.eg.db)
  lt<-merge(
    lt, AnnotationDbi::select(
      org.Mm.eg.db, keys=unique(lt$SYMBOL), 
      columns = c("ENTREZID"), 
      keytype = "SYMBOL"
    ), by='SYMBOL'
  )
  row.names(lt)<-lt$gene_id
  rm(ah, edb)
  detach(package:org.Mm.eg.db, unload = T)
  detach(package:AnnotationHub, unload=T)
  detach(package:ensembldb, unload=T)
  detach(package:AnnotationFilter, unload=T)
  write.csv(lt, fn)
} 

########################### Setup Master DGE List ############################
htseq_count<-hc_buildDataFrame(ds, ft)
rownames(ft)<-ft$sample
ft$sample<-paste(
  ft$genotype, 
  paste(ft$hours.pcs,"H",sep=''), 
  ft$batch,1:36, sep="_"
)
colnames(htseq_count)<-ft[colnames(htseq_count), 'sample']
row.names(ft)<-ft$sample

master<-DGEList(
  htseq_count, 
  samples=ft, 
  genes=lt[row.names(htseq_count),]
)

# Reorder grouping factor, drop unused levels
master$samples$interval<-droplevels(
  factor(
    paste(master$samples$hours.pcs, 'H', sep=''),
    levels = c('0H', '6H', '24H','48H', '120H')
  )
)
master$samples$genotype<-droplevels(
  factor(
    master$samples$genotype,
    levels = c('WT', 'FN', 'B8')
  )
)
master$samples$library <- droplevels(
  factor(
    ifelse(master$samples$batch =="DBI", "single", "paired"),
    levels=c("single", "paired")
  )
)

master$samples$group<-factor(gsub("_[0-9]+$",'',master$samples$sample))
save(master, file = "LTS_DEG_Analysis_results/LTS_DGEList.Rdata")

######################### Analyze Wildtype Samples ###########################
# Global Analysis -- Wildtype Time Series ####

# Get batch subsets to find equally represented genes
### Todo -- CHANGE THIS TO DNA vs DBI -- the differences are
### Likely more related to sequencing methodology!!!


genes <- master$genes$gene_id
for(b in unique(master$samples$library)){
  dge<-master[,
              master$samples$library == b & 
                master$samples$genotype == 'WT'
              ]    # Select Samples
  dge$samples$group<-droplevels(  
    dge$samples$interval
  )
  print(unique(dge$samples$library))
  design<-model.matrix(~interval, dge$samples)    
  colnames(design)<-gsub("interval", '', colnames(design))
  
  genes<-intersect(
    genes,
    (dge[filterByExpr(dge, design), ,keep.lib.sizes=F])$genes$gene_id
  )
  print(length(genes))
}

dge<-master[, master$samples$genotype == 'WT']    # Select Samples
dge$samples$group<-droplevels(  
  dge$samples$interval
)
head(dge)

design<-model.matrix(~interval + batch, dge$samples)    
colnames(design)<-gsub("interval", '', colnames(design))

dge <- dge[genes,,keep.lib.sizes=F]
dge <- calcNormFactors(dge)

# NOTE -- FPKM Calculation Fails to account for VOOM precision Weights !!!
# THIS SHOULD BE FIXED if VOOM is intended for use !!!
# fpkm<-rpkm(dge, normalized.lib.sizes = T, gene.length = "eu_length")
# for(g in unique(dge$samples$group)){
#   s<-row.names(dge$samples[dge$samples$group == g,])
#   dge$genes[paste0(g,"_Avg_FPKM")]<-apply(fpkm[,s],1,mean)
# }

design<-model.matrix(~interval + batch, dge$samples)  
v <- voom(dge, design=design)
fit <- lmFit(v, design=design)
cft <- contrasts.fit(
  fit, contrasts=makeContrasts(
    contrasts = paste("6H", group2, sep="-"),
    levels=colnames(design)
  )
)
ebs <- eBayes(fit)

design.nb<-model.matrix(~interval, dge$samples)  
vnb <- voom(dge, design=design.nb)
fit <- lmFit(vnb, design=design.nb)
cft <- contrasts.fit(
  fit, contrasts=makeContrasts(
    contrasts = paste("6H", group2, sep="-"),
    levels=colnames(design)
  )
)
ebs <- eBayes(fit)


## Plot Log Fold Change
fold_change_table <- bind_rows(
  topTable(ebs, coef=2, confint = T, n=Inf) %>% 
    #filter(abs(logFC) > 1 & adj.P.Val < 0.05 ) %>%
    filter(
      SYMBOL %in% c(
        "Runx1", "Serpine1", "Timp1", 
        "Thbs1", "Acta2", "Fn1"
      )
    ) %>%
    rowwise() %>%
    mutate(
      interval="0H",
      logFC=0,
      CI.L=0,
      CI.R=0
    ) %>%
    mutate(Fold_Change=ifelse(logFC > 0, 2^logFC, 0-(2^abs(logFC)))),
  topTable(ebs, coef=2, confint = T, n=Inf) %>% 
    #filter(abs(logFC) > 1 & adj.P.Val < 0.05 ) %>%
    filter(
      SYMBOL %in% c(
        "Runx1", "Serpine1", "Timp1", 
        "Thbs1", "Acta2", "Fn1"
      )
    ) %>%
    rowwise() %>%
    mutate(
      interval="6H",
      logFC=ifelse(SYMBOL %in% c("Fn1", "Acta2"), NA, logFC),
      CI.L=ifelse(SYMBOL %in% c("Fn1", "Acta2"), NA, CI.L),
      CI.R=ifelse(SYMBOL %in% c("Fn1", "Acta2"), NA, CI.R)
    ) %>%
    mutate(Fold_Change=ifelse(logFC > 0, 2^logFC, 0-(2^abs(logFC)))),
  topTable(ebs, coef=3, confint = T, n=Inf) %>% 
    #filter(abs(logFC) > 1 & adj.P.Val < 0.05 ) %>%
    filter(
      SYMBOL %in% c(
        "Runx1", "Serpine1", "Timp1", 
        "Thbs1", "Acta2", "Fn1"
      )
    ) %>%
    mutate(interval="24H")%>%
    mutate(Fold_Change=ifelse(logFC > 0, 2^logFC, 0-(2^abs(logFC)))),
  topTable(ebs, coef=4, confint = T, n=Inf) %>% 
    #filter(abs(logFC) > 1 & adj.P.Val < 0.05 ) %>%
    filter(
      SYMBOL %in% c(
        "Runx1", "Serpine1", "Timp1", 
        "Thbs1", "Acta2", "Fn1"
      )
    ) %>% mutate(interval="48H")%>%
    mutate(Fold_Change=ifelse(logFC > 0, 2^logFC, 0-(2^abs(logFC)))),
  topTable(ebs, coef=5, confint = T, n=Inf) %>% 
    #filter(abs(logFC) > 1 & adj.P.Val < 0.05 ) %>%
    filter(
      SYMBOL %in% c(
        "Runx1", "Serpine1", "Timp1", 
        "Thbs1", "Acta2", "Fn1"
      )
    ) %>% mutate(interval="120H")%>%
    mutate(Fold_Change=ifelse(logFC > 0, 2^logFC, 0-(2^abs(logFC))))
) %>% mutate(
  interval=factor(interval, levels=c("0H","6H", "24H", "48H", "120H")),
  SYMBOL=factor(
    SYMBOL, levels=c(
      "Runx1", "Serpine1", "Timp1", 
      "Thbs1", "Acta2", "Fn1"
    )
  ),
  FCI.L=ifelse(CI.L > 0, 2^CI.L, 0-(2^abs(CI.L))),
  FCI.R=ifelse(CI.R > 0, 2^CI.R, 0-(2^abs(CI.R)))
) %>% #dplyr::select(SYMBOL, interval, Fold_Change) %>%View()
  arrange(interval, SYMBOL)

fold_change_table%>%
  ggplot(aes(y=logFC, x=interval, fill=SYMBOL)) +
  geom_errorbar(aes(ymin=CI.L, ymax=CI.R), position = position_dodge(width=0.9), width=0.4)+
  geom_bar( stat="identity", position=position_dodge(width=0.9)) + 
  ylab("Fold change in expression after surgery \n (log base 2)") + 
  xlab("Post surgical interval")


fold_change_table%>%
  #filter(SYMBOL %in% c("Runx1", "Acta2", "Fn1")) %>%
  ggplot(aes(y=Fold_Change, x=interval, fill=SYMBOL)) +
  geom_errorbar(aes(ymin=FCI.L, ymax=FCI.R), position = position_dodge(width=0.9), width=0.4)+
  geom_bar( stat="identity", position=position_dodge(width=0.9)) + 
  ylab("Fold change in expression after surgery \n (Scaled to log base 10)") +
  scale_y_log10() +
  xlab("Post surgical interval") +
  labs(fill="Gene")

ggsave("~/Desktop/Injury_response_Runx1_Acta2_Fn1.tiff", width=8, height=4)
fold_change_table%>%
  filter(SYMBOL %in% c("Timp1", "Serpine1", "Thbs1")) %>%
  ggplot(aes(y=Fold_Change, x=interval, fill=SYMBOL)) +
  geom_errorbar(aes(ymin=FCI.L, ymax=FCI.R), position = position_dodge(width=0.9), width=0.4)+
  geom_bar( stat="identity", position=position_dodge(width=0.9)) + 
  ylab("Fold change in expression after surgery \n (Scaled to log base 10)") +
  scale_y_log10() +
  xlab("Post surgical interval") +
  labs(fill="Gene")

ggsave("~/Desktop/Injury_response_Timp1_Serpine1_Thbs1.tiff", width=8, height=4)

targets <- c("Runx1", "Serpine1", "Timp1", "Thbs1", "Acta2", "Fn1")
measurements.d <- rpkm(dge, log=F, gene.length="eu_length")[v$genes$SYMBOL %in% targets,] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene_id") %>%
  inner_join(v$genes %>% dplyr::select(gene_id, SYMBOL), by="gene_id") %>%
  dplyr::select(-gene_id) %>%
  tidyr::pivot_longer(cols=!matches("SYMBOL"), names_to="sample", values_to="CPM") %>%
  inner_join(
    v$targets %>% dplyr::select(sample, group),
    by="sample"
  ) %>%
  group_by(SYMBOL, group) %>%
  summarize(Mean_CPM=mean(CPM), StDEV=sd(CPM), num_obs=n())

measurements.v <- v$E[v$genes$SYMBOL %in% targets,] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene_id") %>%
  inner_join(v$genes %>% dplyr::select(gene_id, SYMBOL), by="gene_id") %>%
  dplyr::select(-gene_id) %>%
  tidyr::pivot_longer(cols=!matches("SYMBOL"), names_to="sample", values_to="CPM") %>%
  inner_join(
    v$targets %>% dplyr::select(sample, group),
    by="sample"
  ) %>%
  group_by(SYMBOL, group) %>%
  summarize(Mean_CPM=mean(2^CPM), StDEV=sd(2^CPM), num_obs=n())


limits <- aes(
  ymax=Mean_CPM + (StDEV/1), 
  ymin=max((Mean_CPM - (StDEV/1))/2, 10)
)

f <- function(x, y){
  return(max(x - (y/1)/2,0.8))
}

measurements.d%>%
ggplot(aes(y=Mean_CPM, fill=group, x=group)) +
  geom_errorbar( aes(
    ymax=Mean_CPM + (StDEV/1), 
    ymin=f( Mean_CPM - (StDEV/1))/2)
  ),  width=0.4) +
  geom_bar( stat="identity") +
  facet_wrap(~SYMBOL, scales="free")
  
measurements.v%>%
  ggplot(aes(y=Mean_CPM, fill=group, x=group)) +
  geom_errorbar(limits,  width=0.4) +
  geom_bar( stat="identity") +
  facet_wrap(~SYMBOL, scales="free")




                