##############################################################################
# File: PCA_24vs0_Hours.R                                                    #
# Created: April 14, 2020                                                    #
# Author: Adam Faranda                                                       #
# Purpose:                                                                   #   
#         Apply various clustering methods to try and identify gene          #
#         Gene clusters with distinct temporal profiles                      #   
#                                                                            #
##############################################################################

############################ Setup Environment ###############################
setwd("~/Documents/LEC_Time_Series")
#load("GeneLengthTable.Rdata")

library(dplyr)
library(cluster)
library(reshape2)
library(org.Mm.eg.db)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)
library(RUVSeq)
wd<-getwd()

source('transcriptomic_analysis_scripts/BuildDataMatrix.R')
source('transcriptomic_analysis_scripts/PreprocessingFunctions.R')
source('transcriptomic_analysis_scripts/PrincipalComponents.R')
source('transcriptomic_analysis_scripts/ClusteringFunctions.R')
source('transcriptomic_analysis_scripts/Overlap_Comparison_Functions.R')
dgeFile="LTS_DEG_Analysis_results/LTS_DGEList.Rdata"      # Rdata file with dgelist object
load(dgeFile)
# 0 vs 24 Hours at 2 labs -- Wildtype Time Series ####
dge<-master[,
            master$samples$genotype == 'WT' &                # Select Samples
              master$samples$interval %in% c('0H', '24H') &
              master$samples$batch %in% c('DBI', 'DNA1')
            ]

dge$samples$interval<-droplevels(dge$samples$interval)
dge$samples$batch<-droplevels(factor(dge$samples$batch))


design<-model.matrix(~0+group, dge$samples)     # Define Experimental Design
colnames(design)<-gsub(
  'group', '', 
  colnames(design)
)

# cntmat<-makeContrasts(
#   WT24vs0H = (WT_24H_DBI + WT_24H_DNA1) - (WT_0H_DBI + WT_0H_DNA1),
#   DNAvsDBI = (WT_0H_DNA1 + WT_24H_DNA1) - (WT_24H_DBI + WT_0H_DBI ),
#   WT0H_DNAvsDBI = WT_0H_DNA1 - WT_0H_DBI,
#   WT24H_DNAvsDBI = WT_24H_DNA1 - WT_24H_DBI,
#   DBI_24vs0 = WT_24H_DBI - WT_0H_DBI,
#   DNA_24vs0 = WT_24H_DNA1 - WT_0H_DNA1,
#   levels = design
# )

cntmat<-makeContrasts(
  #WT24vs0H = (WT_24H_DBI + WT_24H_DNA1) - (WT_0H_DBI + WT_0H_DNA1),
  DNAvsDBI = (WT_0H_DNA1 + WT_24H_DNA1) - (WT_24H_DBI + WT_0H_DBI ),
  levels = design
)


dge<-dge[filterByExpr(dge, design), ,keep.lib.sizes=F] # drop low features
rbg<-as.data.frame(rpkmByGroup(dge))             # Add Group Mean RPKMs 
rbg$gene_id<-row.names(rbg)
dge$genes<-merge(
  dge$genes, rbg,
  by='gene_id'
)
row.names(dge$genes)<-dge$genes$gene_id

dge<-calcNormFactors(dge)                    # Factors & Dispersion
dge<-estimateDisp(dge, design, robust = T)

fit<-glmQLFit(dge, design, robust = T)  # Run Model and estimate DE
qlf<-glmQLFTest(fit, contrast = cntmat)
deg<-as.data.frame(topTags(qlf, n=Inf))
g <- deg %>% filter(abs(logFC) < 2 & FDR > 0.01) %>% pull(gene_id)

# pdf(
#   "~/Desktop/24vs0_Hour_Lab_Comparison_PCA_Plot.pdf", 
#   width=600, height = 400
# )
plotPrinComp(
  cpm(y[,], log=T), ft=dge$samples,                               # PCA Plot
  idCol = 0, groupCol = 'group'
)
ggsave(
  "~/Desktop/24vs0_Hour_Lab_Comparison_PCA_Plot.pdf", 
  width=8, height = 5
)


z <- master[,
            master$samples$genotype == 'WT' &                # Select Samples
              master$samples$interval %in% c('0H', '24H') &
              master$samples$batch %in% c('DBI', 'DNA1')
            ]
z$samples$interval <- droplevels(z$samples$interval)
y <- z

design<-model.matrix(~interval, z$samples)
colnames(design) <-gsub("group","",colnames(design))
z <- z[filterByExpr(z, design),,keep.lib.sizes=F]
z <- calcNormFactors(z)
z <- estimateDisp(z, design, robust=T)
z <- glmQLFit(z, design, robust = T)
print(head(residuals(z, typ="deviance")))

# Estimate nuisance parameter W_1 using the RUVr method
seq <- newSeqExpressionSet(counts = z$counts)
seqUQ <- betweenLaneNormalization(seq, which="upper")
ruvr <- RUVr(
  seqUQ, cIdx = row.names(z),
  k = 1, residuals = residuals(z, type="deviance")
)
y$samples$W_1 <- pData(ruvr)[row.names(y$samples), "W_1"]
design<-model.matrix(~interval + W_1, y$samples)
#colnames(design) <-gsub("inter","",colnames(design))

y$counts <- normCounts(ruvr)
y <- y[filterByExpr(y, design),,keep.lib.sizes=F]
y <- calcNormFactors(y)
y <- estimateDisp(y, design, robust=T)








dge<-master[, master$samples$genotype == 'WT']    # Select Samples
dge$samples$group<-droplevels(  
  dge$samples$interval
)

# Define Experimental Design
design<-model.matrix(~interval + batch, dge$samples)    
colnames(design)<-gsub("interval", '', colnames(design))

dge<-dge[filterByExpr(dge, design), ,keep.lib.sizes=F] 
rbg<-as.data.frame(rpkmByGroup(dge))      
rbg$gene_id<-row.names(rbg)
dge$genes<-merge(
  dge$genes, rbg,
  by='gene_id'
)
row.names(dge$genes)<-dge$genes$gene_id
fpkm<-as.data.frame(rpkm(dge, gene.length="coding_length"))
fpkm$gene_id<-row.names(fpkm)

dge<-calcNormFactors(dge)                    # Factors & Dispersion
dge<-estimateDisp(dge, design, robust = T)



fit<-glmQLFit(dge, design, robust = T)  # Run Model and estimate DE
qlf<-glmQLFTest(fit, contras=c(0,0,0,0,0,1,1))
deg<-as.data.frame(topTags(qlf, n=Inf))
g <- deg %>% filter(abs(logFC) < 10 & FDR > 0.05) %>% pull(gene_id)

plotPrinComp(
  cpm(dge[g,], log=T), ft=dge$samples,                               # PCA Plot
  idCol = 0, groupCol = 'group'
)
ggsave(
  "~/Desktop/24vs0_Hour_Lab_Comparison_PCA_Plot.pdf", 
  width=8, height = 5
)

pca <- prcomp((cpm(dge[g,], log=T)))
length(pca$rotation[,'PC1'])
nrow((cpm(dge[g,], log=T)))

rot <- cpm(dge[g,], log=T) * pca$rotation[,'PC1']
r1 <- colSums(rot)


rot <- cpm(dge[g,], log=T)r1

r2 <- colSums(rot)

x <- 1:5

y <- data.frame(
  V1 = c(
    1,1,1,2,3
  ), 
  V2 = c(
    2,1,1,2,3
  ), 
  V3 = c(
    3,2,1,2,1
  ), 
  V4 = c(
    4, 3,2,2,3
  ),
  V5 = c(
    5,4,3,3,3
  )
)
