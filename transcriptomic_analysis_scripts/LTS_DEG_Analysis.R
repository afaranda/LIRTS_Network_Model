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
#                                                                            #
#                                                                            #                                      
# Author: Adam Faranda                                                       #
##############################################################################

############################ Setup Environment ###############################

# CHANGE THIS DIRECTORY
setwd('~/Documents/LEC_Time_Series')
#setwd("~/Documents/Adam_LEC_Time_Series_DEG_Analysis_3_Oct_2020")
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

library('AnnotationHub')
lt<-read.table(
  paste(dl,"gene_coding_lengths.txt", sep='/'),
  header = T, stringsAsFactors = F
)

###################### Annotate Genes in table 'lt' ##########################
ah<-AnnotationHub()
# Run Query to find proper Annotation Set:
#AnnotationHub::query(ah, pattern=c("EnsDb", "Mus musculus", "98"))
edb<-ah[['AH75036']]
lt<-merge(
  lt, AnnotationDbi::select(
    edb, keys=lt$gene_id, 
    columns = c("SYMBOL", "DESCRIPTION", "GENEBIOTYPE", "SEQCOORDSYSTEM"), 
    keytype = "GENEID"
  ),
  by.x = 'gene_id', by.y='GENEID'
)

lt<-merge(
  lt, AnnotationDbi::select(
    org.Mm.eg.db, keys=unique(lt$SYMBOL), 
    columns = c("ENTREZID"), 
    keytype = "SYMBOL"
  ), by='SYMBOL'
)
row.names(lt)<-lt$gene_id
rm(ah, edb)
detach(package:AnnotationHub, unload=T)
detach(package:ensembldb, unload=T)
detach(package:AnnotationFilter, unload=T)

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
master$samples$group<-factor(gsub("_[0-9]+$",'',master$samples$sample))
save(master, file = "LTS_DEG_Analysis_results/LTS_DGEList.Rdata")

######################### Analyze Wildtype Samples ###########################
# DBI Analysis ####
dge<-master[,
            master$samples$genotype == 'WT' &    # Select Samples
              master$samples$batch == 'DBI' 
            ]
dge$samples$interval<-droplevels(dge$samples$interval)

design<-model.matrix(~0+group, dge$samples)     # Define Experimental Design
colnames(design)<-gsub(
  'group', '', 
  colnames(design)
)
cntmat<-makeContrasts(
  WT24vs0H = WT_24H_DBI - WT_0H_DBI,
  WT48vs0H = WT_48H_DBI - WT_0H_DBI,
  WT48vs24H = WT_48H_DBI - WT_24H_DBI,
  levels = design
)

dge<-dge[filterByExpr(dge, design), ,keep.lib.sizes=F] # Drop low features
rbg<-as.data.frame(rpkmByGroup(dge))             # Add Group Mean RPKMs 
rbg$gene_id<-row.names(rbg)
dge$genes<-merge(
  dge$genes, rbg,
  by='gene_id'
)
row.names(dge$genes)<-dge$genes$gene_id

dge<-calcNormFactors(dge)                  # Factors and Dispersion
dge<-estimateDisp(dge, design, robust = T)

fit<-glmQLFit(dge, design, robust = T)  # Run Model and estimate DE
qlf<-glmQLFTest(fit, contrast = cntmat)
deg<-as.data.frame(topTags(qlf, n=Inf))

## Generate diagnostic plots. 
png("LTS_DEG_Analysis_results/DBI_Experiment_BCV_Plot.png")
plotBCV(dge)                                                 # BCV Plot
dev.off()

png(
  "LTS_DEG_Analysis_results/DBI_Experiment_PCA_Plot.png",   # PCA Plot
  width=600, height = 400
)
plotPrinComp(
  cpm(dge, log=T), ft=dge$samples,                  
  idCol = 0, groupCol = 'group'
)
dev.off()

png(
  "LTS_DEG_Analysis_results/DBI_Experiment_MDS_Plot.png",  # MDS Plot
  width=600, height=500
)
par(mar=c(8.5,5,4.1,2.1))
plotMDS(
  dge, cex=1.5, cex.lab = 1.5,
  cex.axis = 1.5, cex.main = 2,
  pch = rep(15:18, each=3),
  main = "DBI_Experiment",
  col = rep(c('red', 'blue','black'), each=3)
)
legend(
  -3, -2.4, 
  legend = unique(dge$samples$group), 
  pch = 15:17, cex = 1.5,
  col=c('red', 'blue','black'), xpd=T
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
    "LTS_DEG_Analysis_results/DBI_",c,"_Exact_Test_DEG.csv", sep="")
  write.csv(deg.et, fn)
  
  fn<-paste(
    "LTS_DEG_Analysis_results/DBI_",c,"_QLFTest_DEG.csv", sep="")
  write.csv(deg.qt, fn)
  
  # Generate volcano plots 
  fn<-paste(
    "LTS_DEG_Analysis_results/DBI_",c,"_QLF_Volcano.png", sep="")
  ggsave(
    fn,
    EnhancedVolcano(
      toptable = deg.qt, 
      lab = deg.qt$SYMBOL,
      x = "logFC",
      y="FDR",
      ylim =  c(0, max(-log10(deg.qt[,'FDR']), na.rm=TRUE) + 2),
      title = "DBI Experiment",
      subtitle = c,
      pCutoff = 0.05
    )
  )
  
  fn<-paste("LTS_DEG_Analysis_results/DBI_",c,"_Exact_Volcano.png", sep="")
  ggsave(
    fn,
    EnhancedVolcano(
      toptable = deg.et, 
      lab = deg.et$SYMBOL,
      x ="logFC",
      y="FDR",
      ylim =  c(0, max(-log10(deg.et[,'FDR']), na.rm=TRUE) + 2),
      title = "DBI Experiment",
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
  
  fn<-paste("LTS_DEG_Analysis_results/DBI_",c,"_Length_Bias.png", sep="")
  mn<-paste("Length Bias in ",c,sep='')
  yl <-paste(                                    # Construct y axis label
    "Log2 Fold Change in", pair[2],          
    "vs",pair[1]
  )
  
  png(fn, width=6, height = 5, units="in", res=1200)
  plot(
    x=log(deg.et$coding_length,2), y=deg.et$logFC, main=mn,
    xlab = "Log2 Gene Length (exon-union)",
    ylab = yl,
    col = ifelse(abs(deg.et$logFC) > 1 & deg.et$FDR < 0.05, "red", "black")
  )
  
  # Add Correlation Coefficients to tests
  ct=cor.test(deg.et$logFC, log(deg.et$coding_length,2), method="spearman")
  rho=paste("rho:", round(ct$estimate,3))
  sig=paste("p value:", signif(ct$p.value,3), sep="")
  abline(lm(deg.et$logFC~log(deg.et$coding_length,2)), col="red", lwd=2)
  text(6.5,max(deg.et$logFC)-1,rho)
  text(7,max(deg.et$logFC)-3,sig)
  dev.off()
  
  
  
  dg$contrast<-c
  dg$test<-"Exact Test"
  df<-bind_rows(df, dg)
}
write.csv(df,"LTS_DEG_Analysis_results/DBI_Experiment_DEG_Summary_Tables.csv")
write.csv(
  deg,
  "LTS_DEG_Analysis_results/DBI_Experiment_QLFTest_Wildtype_Contrasts_DEG.csv"
)

# DNA Link Analysis -- First Experiment ####
dge<-master[,
            master$samples$genotype == 'WT' &    # Select Samples
              master$samples$batch == 'DNA1' 
            ]
dge$samples$interval<-droplevels(dge$samples$interval)

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

dge<-dge[filterByExpr(dge, design), ,keep.lib.sizes=F] # Drop low features
rbg<-as.data.frame(rpkmByGroup(dge))             # Add Group Mean RPKMs 
rbg$gene_id<-row.names(rbg)
dge$genes<-merge(
  dge$genes, rbg,
  by='gene_id'
)
row.names(dge$genes)<-dge$genes$gene_id

dge<-calcNormFactors(dge)                       # Factors & Dispersion
dge<-estimateDisp(dge, design, robust = T)

fit<-glmQLFit(dge, design, robust = T)  # Run Model and estimate DE
qlf<-glmQLFTest(fit, contrast = cntmat)
deg<-as.data.frame(topTags(qlf, n=Inf))


## Generate diagnostic plots. 
png("LTS_DEG_Analysis_results/DNA_Link_Experiment_1_BCV_Plot.png")
plotBCV(dge)                                              # BCV Plot
dev.off()

png(
  "LTS_DEG_Analysis_results/DNA_Link_Experiment_1_PCA_Plot.png", 
  width=600, height = 400
)
plotPrinComp(
  cpm(dge, log=T), ft=dge$samples,                               # PCA Plot
  idCol = 0, groupCol = 'group'
)
dev.off()

png(
  "LTS_DEG_Analysis_results/DNA_Link_Experiment_1_MDS_Plot.png",           # MDS Plot
  height = 500, width=600
)
par(mar=c(8.5,5,4.1,2.1))
plotMDS(
  dge, cex=1.5, cex.lab = 1.5,
  cex.axis = 1.5, cex.main = 2,
  pch = rep(15:18, each=3),
  main = "DNA_Link_Experiment_1",
  col = rep(c('red', 'blue','black'), each=3)
)
legend(
  -2, -2.65, 
  legend = unique(dge$samples$group), 
  pch = 15:17, cex =1.5,
  col=c('red', 'blue','black'), xpd=T
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
  fn<-paste("LTS_DEG_Analysis_results/DNA_Link_Experiment_1_",c,"_Exact_Test_DEG.csv", sep="")
  write.csv(deg.et, fn)
  
  fn<-paste("LTS_DEG_Analysis_results/DNA_Link_Experiment_1_",c,"_QLFTest_DEG.csv", sep="")
  write.csv(deg.qt, fn)
  
  # Generate volcano plots 
  fn<-paste("LTS_DEG_Analysis_results/DNA_Link_Experiment_1_",c,"_QLF_Volcano.png", sep="")
  ggsave(
    fn,
    EnhancedVolcano(
      toptable = deg.qt, 
      lab = deg.qt$SYMBOL,
      x = "logFC",
      y="FDR",
      ylim =  c(0, max(-log10(deg.qt[,'FDR']), na.rm=TRUE) + 2),
      title = "DNA Link Experiment 1",
      subtitle = c,
      pCutoff = 0.05
    )
  )
  
  fn<-paste("LTS_DEG_Analysis_results/DNA_Link_Experiment_1_",c,"_Exact_Volcano.png", sep="")
  ggsave(
    fn,
    EnhancedVolcano(
      toptable = deg.et, 
      lab = deg.et$SYMBOL,
      x ="logFC",
      y="FDR",
      ylim =  c(0, max(-log10(deg.et[,'FDR']), na.rm=TRUE) + 2),
      title = "DNA Link Experiment 1",
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
    "LTS_DEG_Analysis_results/DNA_Link_Experiment_1_",
    c,"_Length_Bias.png", sep=""
  )
  mn<-paste("Length Bias in ",c,sep='')
  yl <-paste(                                    # Construct y axis label
    "Log2 Fold Change in", pair[2],          
    "vs",pair[1]
  )
  
  png(fn, width=6, height = 5, units="in", res=1200)
  plot(
    x=log(deg.et$coding_length,2), y=deg.et$logFC, main=mn,
    xlab = "Log2 Gene Length (exon-union)",
    ylab = yl,
    col = ifelse(abs(deg.et$logFC) > 1 & deg.et$FDR < 0.05, "red", "black")
  )
  
  # Add Correlation Coefficients to tests
  ct=cor.test(deg.et$logFC, log(deg.et$coding_length,2), method="spearman")
  rho=paste("rho:", round(ct$estimate,3))
  sig=paste("p value:", signif(ct$p.value,3), sep="")
  abline(lm(deg.et$logFC~log(deg.et$coding_length,2)), col="red", lwd=2)
  text(6.5,max(deg.et$logFC)-1,rho)
  text(7,max(deg.et$logFC)-3,sig)
  dev.off()
}

write.csv(df,"LTS_DEG_Analysis_results/DNA_Link_Experiment_1_DEG_Summary_Tables.csv")
write.csv(deg,"LTS_DEG_Analysis_results/DNA_Link_Experiment_1_QLFTest_Wildtype_Contrasts_DEG.csv")

# DNA Link Analysis --  Second Experiment ####
dge<-master[,
            master$samples$genotype == 'WT' &    # Select Samples
              master$samples$batch == 'DNA2' 
            ]
dge$samples$interval<-droplevels(dge$samples$interval)


design<-model.matrix(~0+group, dge$samples)     # Define Experimental Design
colnames(design)<-gsub(
  'group', '', 
  colnames(design)
)
cntmat<-makeContrasts(
  WT120vs0H = WT_120H_DNA2 - WT_0H_DNA2,
  levels = design
)

dge<-dge[filterByExpr(dge, design), ,keep.lib.sizes=F] 
rbg<-as.data.frame(rpkmByGroup(dge))             # Add Group Mean RPKMs 
rbg$gene_id<-row.names(rbg)
dge$genes<-merge(
  dge$genes, rbg,
  by='gene_id'
)
row.names(dge$genes)<-dge$genes$gene_id

dge<-calcNormFactors(dge)                   # Factors & Dispersion
dge<-estimateDisp(dge, design, robust = T)

fit<-glmQLFit(dge, design, robust = T)  # Run Model and estimate DE
qlf<-glmQLFTest(fit, contrast = cntmat)
deg<-as.data.frame(topTags(qlf, n=Inf))

## Generate diagnostic plots. 
png("LTS_DEG_Analysis_results/DNA_Link_Experiment_2_BCV_Plot.png")
plotBCV(dge)                       # BCV Plot
dev.off()

png(
  "LTS_DEG_Analysis_results/DNA_Link_Experiment_2_PCA_Plot.png", 
  width=600, height = 400
)
plotPrinComp(
  cpm(dge, log=T), ft=dge$samples,                               # PCA Plot
  idCol = 0, groupCol = 'group'
)
dev.off()

png(
  "LTS_DEG_Analysis_results/DNA_Link_Experiment_2_MDS_Plot.png",           # MDS Plot
  width=600, height=500
)
par(mar=c(8,5,4.1,2.1))
plotMDS(
  dge, cex=1.5, cex.lab = 1.5,
  cex.axis = 1.5, cex.main = 2,
  pch = rep(15:16, each=3),
  main = "DNA Link Experiment 2",
  col = rep(c('red', 'blue'), each=3)
)
legend(
  -3, -1.4, 
  legend = unique(dge$samples$group), 
  pch = 15:16, cex =1.5,
  col=c('red', 'blue'), xpd=T
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
  fn<-paste("LTS_DEG_Analysis_results/DNA_Link_Experiment_2_",c,"_Exact_Test_DEG.csv", sep="")
  write.csv(deg.et, fn)
  
  fn<-paste("LTS_DEG_Analysis_results/DNA_Link_Experiment_2_",c,"_QLFTest_DEG.csv", sep="")
  write.csv(deg.qt, fn)
  
  # Generate volcano plots 
  fn<-paste("LTS_DEG_Analysis_results/DNA_Link_Experiment_2_",c,"_QLF_Volcano.png", sep="")
  ggsave(
    fn,
    EnhancedVolcano(
      toptable = deg.qt, 
      lab = deg.qt$SYMBOL,
      x = "logFC",
      y="FDR",
      ylim =  c(0, max(-log10(deg.qt[,'FDR']), na.rm=TRUE) + 2),
      title = "DNA Link Experiment 2",
      subtitle = c,
      pCutoff = 0.05
    )
  )
  
  fn<-paste("LTS_DEG_Analysis_results/DNA_Link_Experiment_2_",c,"_Exact_Volcano.png", sep="")
  ggsave(
    fn,
    EnhancedVolcano(
      toptable = deg.et, 
      lab = deg.et$SYMBOL,
      x ="logFC",
      y="FDR",
      ylim =  c(0, max(-log10(deg.et[,'FDR']), na.rm=TRUE) + 2),
      title = "DNA Link Experiment 2",
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
  
  # Generate Length Bias Plots
  fn<-paste(
    "LTS_DEG_Analysis_results/DNA_Link_Experiment_2_",
    c,"_Length_Bias.png", sep=""
  )
  mn<-paste("Length Bias in ",c,sep='')
  yl <-paste(                                    # Construct y axis label
    "Log2 Fold Change in", pair[2],          
    "vs",pair[1]
  )
  png(fn, width=6, height = 5, units="in", res=1200)
  plot(
    x=log(deg.et$coding_length,2), y=deg.et$logFC, main=mn,
    xlab = "Log2 Gene Length (exon-union)",
    ylab = yl,
    col = ifelse(abs(deg.et$logFC) > 1 & deg.et$FDR < 0.05, "red", "black")
  )
  
  # Add Correlation Coefficients to tests
  ct=cor.test(deg.et$logFC, log(deg.et$coding_length,2), method="spearman")
  rho=paste("rho:", round(ct$estimate,3))
  sig=paste("p value:", signif(ct$p.value,3), sep="")
  abline(lm(deg.et$logFC~log(deg.et$coding_length,2)), col="red", lwd=2)
  text(6.5,max(deg.et$logFC)-1,rho)
  text(7,max(deg.et$logFC)-3,sig)
  dev.off()
  
}

write.csv(df,"LTS_DEG_Analysis_results/DNA_Link_Experiment_2_DEG_Summary_Tables.csv")
write.csv(deg,"LTS_DEG_Analysis_results/DNA_Link_Experiment_2_QLFTest_DEG.csv")

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

cntmat<-makeContrasts(
  WT24vs0H = (WT_24H_DBI + WT_24H_DNA1) - (WT_0H_DBI + WT_0H_DNA1),
  DNAvsDBI = (WT_0H_DNA1 + WT_24H_DNA1) - (WT_24H_DBI + WT_0H_DBI ),
  WT0H_DNAvsDBI = WT_0H_DNA1 - WT_0H_DBI,
  WT24H_DNAvsDBI = WT_24H_DNA1 - WT_24H_DBI,
  DBI_24vs0 = WT_24H_DBI - WT_0H_DBI,
  DNA_24vs0 = WT_24H_DNA1 - WT_0H_DNA1,
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

## Generate diagnostic plots. 
png(
  "LTS_DEG_Analysis_results/24vs0_Hour_Lab_Comparison_BCV_Plot.png"
)
plotBCV(dge)                                              # BCV Plot
dev.off()

png(
  "LTS_DEG_Analysis_results/24vs0_Hour_Lab_Comparison_PCA_Plot.png", 
  width=600, height = 400
)
plotPrinComp(
  cpm(dge, log=T), ft=dge$samples,                               # PCA Plot
  idCol = 0, groupCol = 'group'
)
dev.off()

png(
  "LTS_DEG_Analysis_results/24vs0_Hour_Lab_Comparison_MDS_Plot.png",       # MDS Plot
  width = 600, height = 500
)
par(mar=c(10,5,4.1,2.1))
plotMDS(
  dge, cex=1.5, cex.lab = 1.5,
  cex.axis = 1.5, cex.main = 2,
  pch = rep(15:18, each=3),
  main = "24 vs 0 Hours at DBI and DNA Link",
  col = rep(c('red', 'blue','green4','black'), each=3)
)
legend(
  -3.4, -3.3, 
  legend = unique(dge$samples$group), 
  pch = 15:18, cex =1.5,
  col=c('red', 'blue','green4','black'), xpd=T
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
  g<-dge$samples$group
  if(c == 'WT24vs0H'){
    dge$samples$group <- factor(
      dge$samples$interval,
      levels = c("0H", "24H")
    )
    deg.et<-as.data.frame(
      topTags(exactTest(dge, pair = c("0H", "24H")), n=Inf)
    )
    dge$samples$group<-g
  } else if(c == 'DNAvsDBI'){
    dge$samples$group <- factor(
      dge$samples$batch,
      levels = c("DBI", "DNA1")
    )
    deg.et<-as.data.frame(
      topTags(exactTest(dge, pair = c("DBI", "DNA1")), n=Inf)
    )
    dge$samples$group<-g
  } else {
    pair<-c(
      names(cntmat[,c])[cntmat[,c] == -1], 
      names(cntmat[,c])[cntmat[,c] == 1]
    )
    deg.et<-as.data.frame(
      topTags(exactTest(dge, pair = pair), n=Inf)
    )
  }
  deg.qt<-as.data.frame(topTags(glmQLFTest(fit,contrast=cntmat[,c]), n=Inf))
  fn<-paste("LTS_DEG_Analysis_results/24vs0_Hour_Lab_Comparison_",c,"_Exact_Test_DEG.csv", sep="")
  write.csv(deg.et, fn)
  
  fn<-paste("LTS_DEG_Analysis_results/24vs0_Hour_Lab_Comparison_",c,"_QLFTest_DEG.csv", sep="")
  write.csv(deg.qt, fn)
  
  # Generate volcano plots 
  fn<-paste("LTS_DEG_Analysis_results/24vs0_Hour_Lab_Comparison_",c,"_QLF_Volcano.png", sep="")
  ggsave(
    fn,
    EnhancedVolcano(
      toptable = deg.qt, 
      lab = deg.qt$SYMBOL,
      x = "logFC",
      y="FDR",
      ylim =  c(0, max(-log10(deg.qt[,'FDR']), na.rm=TRUE) + 2),
      title = "24 vs 0 Hour Lab Comparison",
      subtitle = c,
      pCutoff = 0.05
    )
  )
  
  fn<-paste("LTS_DEG_Analysis_results/24vs0_Hour_Lab_Comparison_",c,"_Exact_Volcano.png", sep="")
  ggsave(
    fn,
    EnhancedVolcano(
      toptable = deg.et, 
      lab = deg.et$SYMBOL,
      x ="logFC",
      y="FDR",
      ylim =  c(0, max(-log10(deg.et[,'FDR']), na.rm=TRUE) + 2),
      title = "24 vs 0 Hour Lab Comparison",
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
  
  # Generate Length Bias Plots
  fn<-paste(
    "LTS_DEG_Analysis_results/24vs0_Hour_Lab_Comparison_",
    c,"_Length_Bias.png", sep=""
  )
  mn<-paste("Length Bias in ",c,sep='')
  yl <-paste(                                    # Construct y axis label
    "Log2 Fold Change in", pair[2],          
    "vs",pair[1]
  )
  
  png(fn, width=6, height = 5, units="in", res=1200)
  plot(
    x=log(deg.et$coding_length,2), y=deg.et$logFC, main=mn,
    xlab = "Log2 Gene Length (exon-union)",
    ylab = yl,
    col = ifelse(abs(deg.et$logFC) > 1 & deg.et$FDR < 0.05, "red", "black")
  )
  
  # Add Correlation Coefficients to tests
  ct=cor.test(deg.et$logFC, log(deg.et$coding_length,2), method="spearman")
  rho=paste("rho:", round(ct$estimate,3))
  sig=paste("p value:", signif(ct$p.value,3), sep="")
  abline(lm(deg.et$logFC~log(deg.et$coding_length,2)), col="red", lwd=2)
  text(6.5,max(deg.et$logFC)-1,rho)
  text(7,max(deg.et$logFC)-3,sig)
  dev.off()
}

write.csv(df,"LTS_DEG_Analysis_results/24vs0_Hour_Lab_Comparison_DEG_Summary_Tables.csv")
write.csv(deg,"LTS_DEG_Analysis_results/24vs0_Hour_Lab_Comparison_QLFTest_DEG.csv")


# Global Analysis -- Wildtype Time Series ####
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
qlf<-glmQLFTest(fit, coef=2:5)
deg<-as.data.frame(topTags(qlf, n=Inf))


## Generate diagnostic plots. 
png(
  "LTS_DEG_Analysis_results/Global_Wildtype_BCV_Plot.png"
)
plotBCV(dge)                                              # BCV Plot
dev.off()

png(
  "LTS_DEG_Analysis_results/Global_Wildtype_PCA_Plot.png", 
  width=600, height = 400
)
plotPrinComp(
  cpm(dge, log=T), ft=dge$samples,                               # PCA Plot
  idCol = 0, groupCol = 'group'
)
dev.off()

png(
  "LTS_DEG_Analysis_results/Global_Wildtype_MDS_Plot.png",                      # MDS Plot
  width = 600, height = 500
)
par(mar=c(10,5,4.1,2.1))
plotMDS(
  dge, cex=1.5, cex.lab = 1.5,
  cex.axis = 1.5, cex.main = 2,
  pch = rep(15:19, each=3),
  main = "All intervals, DBI and DNA Link",
  col = rep(c('red', 'blue','green4','black', 'grey'), each=3)
)
legend(
  -3.4, -3.3, 
  legend = unique(dge$samples$group), 
  pch = 15:19, cex =1.5,
  col=c('red', 'blue','green4','black', 'grey'), xpd=T
)
dev.off()

## Iterate over contrasts and prepare summary figures / tables
df<-data.frame()
coefs=c(WT6='6H', WT24='24H', WT48='48H', WT120='120H')
for (c in c('WT6', 'WT24', 'WT48', 'WT120')) {             
  deg.qt<-as.data.frame(
    topTags(glmQLFTest(fit, coef=coefs[c]), n=Inf)
  )

  fn<-paste("LTS_DEG_Analysis_results/Global_Wildtype_",c,"_QLFTest_DEG.csv", sep="")
  write.csv(deg.qt, fn)

  fn<-paste(                                          # Generate volcano plots
    "LTS_DEG_Analysis_results/Global_Wildtype_",c,
    "vs0H_QLF_Volcano.png", sep=""
  )
  ggsave(
    fn,
    EnhancedVolcano(
      toptable = deg.qt, 
      lab = deg.qt$SYMBOL,
      x = "logFC",
      y="FDR",
      ylim =  c(0, max(-log10(deg[,'FDR']), na.rm=TRUE) + 2),
      title = "Global Wildtype Analysis",
      subtitle = paste(c,"0H",sep="vs" ),
      pCutoff = 0.05
    )
  )
  names(deg.qt)[names(deg.qt) %in% c('0H', coefs[c])] <-c('WT0', c)
  dg<-degSummary(                         # Generate QLF Test Summary tables
    deg.qt,
    lfc = "logFC",
    fdr = 'FDR',
    Avg1 = 'WT0',
    Avg2 = c
  )
  dg$contrast<-paste(c,"0H", sep="vs")
  dg$test<-"QLFTest"
  df<-bind_rows(df, dg)
  
  # Generate Length Bias Plots
  fn<-paste(
    "LTS_DEG_Analysis_results/Global_Wildtype_",
    c,"vs0H_Length_Bias.png", sep=""
  )
  mn<-paste("Length Bias in ",c,sep='')
  yl <-paste(                                    # Construct y axis label
    "Log2 Fold Change in", c,          
    "vs 0H"
  )
  
  png(fn, width=6, height = 5, units="in", res=1200)
  plot(
    x=log(deg.qt$coding_length,2), y=deg.qt$logFC, main=mn,
    xlab = "Log2 Gene Length (exon-union)",
    ylab = yl,
    col = ifelse(abs(deg.qt$logFC) > 1 & deg.qt$FDR < 0.05, "red", "black")
  )
  
  # Add Correlation Coefficients to tests
  ct=cor.test(deg.qt$logFC, log(deg.qt$coding_length,2), method="spearman")
  rho=paste("rho:", round(ct$estimate,3))
  sig=paste("p value:", signif(ct$p.value,3), sep="")
  abline(lm(deg.qt$logFC~log(deg.qt$coding_length,2)), col="red", lwd=2)
  text(6.5,max(deg.qt$logFC)-1,rho)
  text(7,max(deg.qt$logFC)-3,sig)
  dev.off()
  
}
## Generate Time Series Plots for Selected Genes
#List of Genes over which to tabulate time series
# Un-comment genes that you would like to generate a table for
genes<-c(
  Egr1="ENSMUSG00000038418", Egr2="ENSMUSG00000037868", Egr3="ENSMUSG00000033730", 
  Grem1="ENSMUSG00000074934", Ptx3="ENSMUSG00000027832", Fn1="ENSMUSG00000026193",
  Actn1="ENSMUSG00000015143", Acta2="ENSMUSG00000035783", Col1a1="ENSMUSG00000001506",
  Atf3="ENSMUSG00000026628", Pitx3="ENSMUSG00000025229", Klf2="ENSMUSG00000055148",
  Emp1="ENSMUSG00000030208", Cxcl1="ENSMUSG00000029380", Cxcl2="ENSMUSG00000058427",
  Cxcl5="ENSMUSG00000029371", Ppbp="ENSMUSG00000029372", #Cxcl9="ENSMUSG00000029417",
  Cxcr2="ENSMUSG00000026180", #Cxcr4="ENSG00000121966", Runx1="ENSMUSG00000022952",
  Tnc="ENSMUSG00000028364", Fos="ENSMUSG00000021250", Ier2="ENSMUSG00000053560",
  Fosb="ENSMUSG00000003545", Jun="ENSMUSG00000052684", Junb="ENSMUSG00000052837",
  Timp1 = 'ENSMUSG00000001131', Serpine1 = 'ENSMUSG00000037411', Hmga1 = 'ENSMUSG00000046711',
  Thbs1 ='ENSMUSG00000040152', Ier3 = 'ENSMUSG00000003541', Esm1 = 'ENSMUSG00000042379',
  Zc3h12a = 'ENSMUSG00000042677', Dusp6 = 'ENSMUSG00000019960', Tcim  = 'ENSMUSG00000056313',
  Mrtfa = 'ENSMUSG00000042292', Ch25h = 'ENSMUSG00000050370', Tagln = 'ENSMUSG00000032085',
  Tagln2 = 'ENSMUSG00000026547', Tagln3 = 'ENSMUSG00000022658', Hsbp1 = 'ENSMUSG00000031839',
  Dkk3 = 'ENSMUSG00000030772', Wnt16 = 'ENSMUSG00000029671', Ccl2 = "ENSMUSG00000035385", 
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
  Serpina3c="ENSMUSG00000066361", Serpina3h="ENSMUSG00000041449", Serpina3f="ENSMUSG00000066363",
  Serpina3m="ENSMUSG00000079012", Serpina3i="ENSMUSG00000079014", Serpina3g="ENSMUSG00000041481"
)

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
      filter(variable != "coding_length") %>%
      dplyr::select(sample=variable, FPKM=value),
    by="sample"
  )
  plt$gene_symbol<-g
  plt$ensembl_id<-genes[g]
  ggplot(plt, aes(x=hours.pcs, y=FPKM, group=interval)) + 
    geom_boxplot() + ggtitle(g)
  ggsave(
    filename=paste(
      "LTS_DEG_Analysis_results/",g,
      "_fpkm_box.png", sep=""
    ), width=6, height=2
  )
  
  ggplot(plt, aes(x=hours.pcs, y=FPKM, group=interval)) + 
    geom_point(aes(color=batch)) + ggtitle(g)
  ggsave(
    filename=paste(
      "LTS_DEG_Analysis_results/",g,
      "_fpkm_pnt.png", sep=""
    ), width=6, height=2
  )
  
  fn<- paste("LTS_DEG_Analysis_results/",g,"_FPKM_Time.tsv", sep="")
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
write.csv(df,"LTS_DEG_Analysis_results/Global_Wildtype_DEG_Summary_Tables.csv")
write.csv(deg,"LTS_DEG_Analysis_results/Global_Wildtype_QLFTest_DEG.csv")
write.csv(
  fpkm, 
  "LTS_DEG_Analysis_results/Global_Wildtype_QLFTest_FPKM_Matrix.csv"
)

# edgeR TMM/cpm (log2 transformed) used for clustering (whole data set)
c<-cpm(dge, log=T)
write.csv(
  c, 
  "LTS_DEG_Analysis_results/Global_Wildtype_TMM_Normalized_Present_Genes.txt", 
  row.names=F
)

png(
  "LTS_DEG_Analysis_results/Global_Wildtype_PCA_Plot_ByBatch.png",                     # PCA Plot
  width=600, height = 400
)
plotPrinComp(
  c, ft=dge$samples,                  
  idCol = 0, groupCol = 'batch'
)
dev.off()


# edgeR TMM/cpm (batch corrected) used for clustering (whole data set)
design.exp<-design[,setdiff(colnames(design), c('batchDNA1', 'batchDNA2'))]
c.bat<-removeBatchEffect(
  c, design = design.exp, batch=dge$samples$batch
)

write.csv(
  c.bat, 
  "LTS_DEG_Analysis_results/Global_Wildtype_TMM_Batch_Corrected_Present_Genes.txt",
  row.names=F
)

png(
  "LTS_DEG_Analysis_results/Global_Wildtype_Model_BatchCorrected_PCA_Plot_ByInterval.png",
  width=600, height = 400
)
plotPrinComp(
  c.bat, ft=dge$samples,                  
  idCol = 0, groupCol = 'group'
)
dev.off()

png(
  "LTS_DEG_Analysis_results/Global_Wildtype_BatchCorrected_PCA_Plot_ByBatch.png",
  width=600, height = 400
)
plotPrinComp(
  c.bat, ft=dge$samples,                  
  idCol = 0, groupCol = 'batch'
)
dev.off()

######################## Analyze Mutation effects ############################
# DBI -- Fibronectin vs Wildtype 0 and 48 hours ####
dge<-master[,
              master$samples$interval %in% c('0H', '48H')&    # Select Samples
              master$samples$batch == 'DBI',
            ]
dge$samples$interval<-droplevels(dge$samples$interval)
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

## Generate diagnostic plots. 
png("LTS_DEG_Analysis_results/Fibronectin_Experiment_BCV_Plot.png")
plotBCV(dge)                                                 # BCV Plot
dev.off()

png(
  "LTS_DEG_Analysis_results/Fibronectin_Experiment_PCA_Plot.png",  # PCA Plot
  width=600, height = 400
)
plotPrinComp(
  cpm(dge, log=T), ft=dge$samples,                  
  idCol = 0, groupCol = 'group'
)
dev.off()

png(
  "LTS_DEG_Analysis_results/Fibronectin_Experiment_MDS_Plot.png",  # MDS Plot
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
    "LTS_DEG_Analysis_results/Fibronectin_Experiment_",
    c,"_Exact_Test_DEG.csv", sep="")
  write.csv(deg.et, fn)
  
  fn<-paste(
    "LTS_DEG_Analysis_results/Fibronectin_Experiment_",
    c,"_QLFTest_DEG.csv", sep="")
  write.csv(deg.qt, fn)
  
  # Generate volcano plots 
  fn<-paste(
    "LTS_DEG_Analysis_results/Fibronectin_Experiment_",
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
    "LTS_DEG_Analysis_results/Fibronectin_Experiment_",
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
    "LTS_DEG_Analysis_results/Fibronectin_Experiment_",
    c,"_Length_Bias.png", sep=""
  )
  mn<-paste("Length Bias in ",c,sep='')
  yl <-paste(                                    # Construct y axis label
    "Log2 Fold Change in", pair[2],          
    "vs",pair[1]
  )
  png(fn, width=6, height = 5, units="in", res=1200)
  plot(
    x=log(deg.et$coding_length,2), y=deg.et$logFC, main=mn,
    xlab = "Log2 Gene Length (exon-union)",
    ylab = yl,
    col = ifelse(abs(deg.et$logFC) > 1 & deg.et$FDR < 0.05, "red", "black")
  )
  
  # Add Correlation Coefficients to tests
  ct=cor.test(deg.et$logFC, log(deg.et$coding_length,2), method="spearman")
  rho=paste("rho:", round(ct$estimate,3))
  sig=paste("p value:", signif(ct$p.value,3), sep="")
  abline(lm(deg.et$logFC~log(deg.et$coding_length,2)), col="red", lwd=2)
  text(6.5,max(deg.et$logFC)-1,rho)
  text(7,max(deg.et$logFC)-3,sig)
  dev.off()
}
write.csv(
  df,"LTS_DEG_Analysis_results/Fibronectin_Experiment_DEG_Summary_Tables.csv")
write.csv(
  deg,
  "LTS_DEG_Analysis_results/Fibronectin_Experiment_QLFTest_All_Contrasts_DEG.csv"
)


# DNA Link -- Beta 8 Integrin vs Wildtype 0 and 48 hours ####
dge<-master[,
            master$samples$interval %in% c('0H', '24H')&    # Select Samples
              master$samples$batch == 'DNA1',
            ]

dge$samples$interval<-droplevels(dge$samples$interval)
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
rbg<-as.data.frame(rpkmByGroup(dge))             # Add Group Mean RPKMs 
rbg$gene_id<-row.names(rbg)
dge$genes<-merge(
  dge$genes, rbg,
  by='gene_id'
)
row.names(dge$genes)<-dge$genes$gene_id

dge<-calcNormFactors(dge)                      # Factors & Dispersion
dge<-estimateDisp(dge, design, robust = T)

fit<-glmQLFit(dge, design, robust = T)  # Run Model and estimate DE
qlf<-glmQLFTest(fit, contrast = cntmat)
deg<-as.data.frame(topTags(qlf, n=Inf))

## Generate diagnostic plots. 
png("LTS_DEG_Analysis_results/Beta8_Integrin_Experiment_BCV_Plot.png")
plotBCV(dge)                                                 # BCV Plot
dev.off()

png(
  "LTS_DEG_Analysis_results/Beta8_Integrin_Experiment_PCA_Plot.png",# PCA Plot
  width=600, height = 400
)
plotPrinComp(
  cpm(dge, log=T), ft=dge$samples,                  
  idCol = 0, groupCol = 'group'
)
dev.off()

png(
  "LTS_DEG_Analysis_results/Beta8_Integrin_Experiment_MDS_Plot.png",# MDS Plot
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
    "LTS_DEG_Analysis_results/Beta8_Integrin_Experiment_",
    c,"_Exact_Test_DEG.csv", sep="")
  write.csv(deg.et, fn)
  
  fn<-paste(
    "LTS_DEG_Analysis_results/Beta8_Integrin_Experiment_",
    c,"_QLFTest_DEG.csv", sep="")
  write.csv(deg.qt, fn)
  
  # Generate volcano plots 
  fn<-paste(
    "LTS_DEG_Analysis_results/Beta8_Integrin_Experiment_",
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
    "LTS_DEG_Analysis_results/Beta8_Integrin_Experiment_",
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
    "LTS_DEG_Analysis_results/Beta8_Integrin_Experiment_",
    c,"_Length_Bias.png", sep=""
  )
  mn<-paste("Length Bias in ",c,sep='')
  yl <-paste(                                    # Construct y axis label
    "Log2 Fold Change in", pair[2],          
    "vs",pair[1]
  )
  
  png(fn, width=6, height = 5, units="in", res=1200)
  plot(
    x=log(deg.et$coding_length,2), y=deg.et$logFC, main=mn,
    xlab = "Log2 Gene Length (exon-union)",
    ylab = yl,
    col = ifelse(abs(deg.et$logFC) > 1 & deg.et$FDR < 0.05, "red", "black")
  )
  
  # Add Correlation Coefficients to tests
  ct=cor.test(deg.et$logFC, log(deg.et$coding_length,2), method="spearman")
  rho=paste("rho:", round(ct$estimate,3))
  sig=paste("p value:", signif(ct$p.value,3), sep="")
  abline(lm(deg.et$logFC~log(deg.et$coding_length,2)), col="red", lwd=2)
  text(6.5,max(deg.et$logFC)-1,rho)
  text(7,max(deg.et$logFC)-3,sig)
  dev.off()
}
write.csv(
  df,
  "LTS_DEG_Analysis_results/Beta8_Integrin_Experiment_DEG_Summary_Tables.csv")
write.csv(
  deg,
  "LTS_DEG_Analysis_results/Beta8_Integrin_Experiment_QLFTest_All_Contrasts_DEG.csv"
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
  "LTS_DEG_Analysis_results/Injury_Model_Raw_Counts_All_Genes.txt", 
  row.names=F
)

# Length normalized, log2 transformed FPKM values (no TMM scaling)
r<-as.data.frame(rpkm(master, log=T))
s<-colnames(r)
r$gene_id<-row.names(r)
r<-r[,c('gene_id', s)]

write.csv(
  r, 
  "LTS_DEG_Analysis_results/Full_Injury_Model_FPKM_All_Genes.csv", 
  row.names = F
)

# Sample ID's and Covariates
write.csv(
  master$samples, 
  "LTS_DEG_Analysis_results/Full_Injury_Model_Sample_Metadata.txt", 
  row.names=F
)

# Gene ID's with Symbol, name and Length Used for RPKM normalization
write.csv(
  master$genes, 
  "LTS_DEG_Analysis_results/Full_Injury_Model_Gene_Metadata.txt", 
  row.names=F
)

# edgeR TMM/cpm (log2 transformed) used for clustering (whole data set)
c<-cpm(dge, log=T)
write.csv(
  c, 
  "LTS_DEG_Analysis_results/Full_Injury_Model_TMM_Normalized_Present_Genes.txt", 
  row.names=F
)

png(
  "LTS_DEG_Analysis_results/Full_Injury_Model_PCA_Plot_ByInterval.png",     # PCA Plot
  width=600, height = 400
)
plotPrinComp(
  c, ft=dge$samples,                  
  idCol = 0, groupCol = 'group'
)
dev.off()

png(
  "LTS_DEG_Analysis_results/Full_Injury_Model_PCA_Plot_ByBatch.png",        # PCA Plot
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
  "LTS_DEG_Analysis_results/Full_Injury_Model_TMM_Batch_Corrected_Present_Genes.txt",
  row.names=F
)

png(
  "LTS_DEG_Analysis_results/Full_Injury_Model_BatchCorrected_PCA_Plot_ByInterval.png",
  width=600, height = 400
)
plotPrinComp(
  c.bat, ft=dge$samples,                  
  idCol = 0, groupCol = 'group'
)
dev.off()

png(
  "LTS_DEG_Analysis_results/Full_Injury_Model_BatchCorrected_PCA_Plot_ByBatch.png",
  width=600, height = 400
)
plotPrinComp(
  c.bat, ft=dge$samples,                  
  idCol = 0, groupCol = 'batch'
)
dev.off()



