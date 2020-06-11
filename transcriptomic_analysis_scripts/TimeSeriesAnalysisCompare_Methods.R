############################### H E A D E R ####################################
# File: TimeSeriesAnalysis.R                                                   #
# Purpose: Compare Results from differerent edgeR Calls on the same data.      #
# Created: May 1, 2019                                                         #
# Author: Adam Faranda                                                         #
################################################################################

############################ Setup Environment #################################
  
# CHANGE THIS DIRECTORY
setwd('/home/adam/Documents/LEC_Time_Series')
#load("GeneLengthTable.Rdata")
library(dplyr)
library(cluster)
library(reshape2)
library(pheatmap)

wd<-getwd()
source('transcriptomic_analysis_scripts/BuildDataMatrix.R')
source('transcriptomic_analysis_scripts/PreprocessingFunctions.R')
source('transcriptomic_analysis_scripts/PrincipalComponents.R')
source('transcriptomic_analysis_scripts/ClusteringFunctions.R')

############################ Load in Data Files ################################

# CHANGE THIS DIRECTORY #
dl<-"~/Documents/LEC_Time_Series_HTSeq_Counts"

ft<-hc_getFileTable(
  dirList=dl, filename = "HTSeq_GeneCounts_All.csv"
)

# ft<-hc_getFileTable(
#   dirList=dl, filename = "HTSeq_GeneCounts_Wildtype.csv"
# )

# ft<-hc_getFileTable(
#   dirList=dl, filename = "H TSeq_GeneCounts_Wildtype_0-48.csv"
# )

ds<-hc_loadFiles(ft)
ft<-hc_identifierConsistency(ds, ft)
######################### Annotate Gene Length File ##########################     
library('AnnotationHub')
lt<-read.table(
  paste(dl,"gene_coding_lengths.txt", sep='/'),
  header = T, stringsAsFactors = F
)
ah<-AnnotationHub()
# Run Query to find proper Annotation Set:
#AnnotationHub::query(ah, pattern=c("EnsDb", "Mus musculus", "98"))
edb<-ah[['AH75036']]
lt<-merge(
  lt, AnnotationDbi::select(
    edb, keys=lt$gene_id, 
    columns = c("SYMBOL", "DESCRIPTION"), 
    keytype = "GENEID"
  ),
  by.x = 'gene_id', by.y='GENEID'
)

row.names(lt)<-lt$gene_id
rm(ah, edb)
detach(package:AnnotationHub, unload=T)
detach(package:ensembldb, unload=T)
detach(package:AnnotationFilter, unload=T)

####################### Setup master edgeR DGEList ###########################
htseq_count<-hc_buildDataFrame(ds, ft)
master<-DGEList(
  htseq_count, 
  samples=ft, 
  genes=lt[row.names(htseq_count),]
)
master$genes$Var<-apply(master$counts, 1, var)

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


############ Wildtype andBeta 8 Integrin 0vs24 Hours Combined ################
dge<-master[,
               master$samples$batch == "DNA1" & 
                 master$samples$hours.pcs %in% c(0,24)
               ]
dge$samples$genotype<-dropEmptyLevels(dge$samples$genotype)
dge$samples$interval<-dropEmptyLevels(dge$samples$interval)
dge$samples$group<-factor(
  paste(
    dge$samples$genotype, 
    dge$samples$interval, sep="_" 
  ), levels=c("WT_0H","WT_24H","B8_0H", "B8_24H")
)
design<-model.matrix(~0+group, dge$samples)
colnames(design)<-gsub('group', '', colnames(design))
dge<-dge[filterByExpr(dge, design), , keep.lib.size=FALSE]

cmat<-makeContrasts(
  B8_0H-WT_0H, B8_24H-WT_24H,
  WT_24H-WT_0H, B8_24H-B8_0H,
  levels=design)
colnames(cmat)<-c(
  "Hour_0", "Hour_24",
  "WT", "B8"
)
dge<-calcNormFactors(dge)
dge<-estimateDisp(dge, design, robust=T)
qlf<-glmQLFit(dge, design)
qlt<-glmQLFTest(qlf, contrast = cmat)
deg<-as.data.frame(topTags(qlt, n=Inf))
degQ<-data.table(deg)

et<-exactTest(dge, pair=c('WT_0H', 'WT_24H'))
deg<-as.data.frame(topTags(et, n=Inf))
degE<-data.table(deg)

qt_genes<-degQ[FDR < 0.05 & abs(logFC.WT) > 1, ]$gene_id
et_genes<-degE[FDR < 0.05 & abs(logFC) > 1, ]$gene_id

q_not_e<-setdiff(qt_genes, et_genes)
plot(
  degQ[gene_id %in% q_not_e,][order(gene_id)]$logFC.WT,
  degE[gene_id %in% q_not_e,][order(gene_id)]$logFC
)

e_not_q<-setdiff(et_genes, qt_genes)
plot(
  degQ[gene_id %in% e_not_q,][order(gene_id)]$logFC.WT,
  degE[gene_id %in% e_not_q,][order(gene_id)]$logFC
)

## Pull out the 0 vs 24 Hour contrast into 
dge<-master[,
            master$samples$batch == "DNA1" &
              master$samples$genotype == "WT" &
              master$samples$hours.pcs %in% c(0,24)
            ]
dge$samples$genotype<-dropEmptyLevels(dge$samples$genotype)
dge$samples$interval<-dropEmptyLevels(dge$samples$interval)
dge$samples$group<-factor(
  paste(
    dge$samples$genotype, 
    dge$samples$interval, sep="_" 
  ), levels=c("WT_0H","WT_24H")
)
design<-model.matrix(~0+group, dge$samples)
colnames(design)<-gsub('group', '', colnames(design))
dge<-dge[filterByExpr(dge, design), , keep.lib.size=FALSE]

cmat<-makeContrasts(WT_24H-WT_0H,levels=design)
colnames(cmat)<-c("WT")

dge<-calcNormFactors(dge)
dge<-estimateDisp(dge, design, robust=T)
qlf<-glmQLFit(dge, design)
qlt<-glmQLFTest(qlf, contrast = cmat)
deg<-as.data.frame(topTags(qlt, n=Inf))
degQb<-data.table(deg)

et<-exactTest(dge, pair=c('WT_0H', 'WT_24H'))
deg<-as.data.frame(topTags(et, n=Inf))
degEb<-data.table(deg)

qt_genes_b<-degQb[FDR < 0.05 & abs(logFC) > 1, ]$gene_id
et_genes_b<-degEb[FDR < 0.05 & abs(logFC) > 1, ]$gene_id

q_not_e<-setdiff(qt_genes, et_genes)
plot(
  degQ[gene_id %in% q_not_e,][order(gene_id)]$logFC.WT,
  degE[gene_id %in% q_not_e,][order(gene_id)]$logFC
)

e_not_q<-setdiff(et_genes, qt_genes)
plot(
  degQ[gene_id %in% e_not_q,][order(gene_id)]$logFC.WT,
  degE[gene_id %in% e_not_q,][order(gene_id)]$logFC
)

qbg<-setdiff(qt_genes_b, qt_genes)
qag<-setdiff(qt_genes, qt_genes_b)




### How many different ways to get pairwise DEG
#
# Robust: True or False
#  - True to reduce outlier shrinkage
#
# Test: Exact or Quasi Likelihood
#  - Use Quasi Likelihood 
#
# Data: Separate or Full data set. 
#
#











