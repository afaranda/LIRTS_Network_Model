################################################################################
# File: TimeSeriesAnalysis.R                                                   #
# Purpose: Identify gene clusters that correspond to temporal patterns.        #
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
library(RColorBrewer)
library(TMixClust)
wd<-getwd()
source('transcriptomic_analysis_scripts/BuildDataMatrix.R')
source('transcriptomic_analysis_scripts/PreprocessingFunctions.R')
source('transcriptomic_analysis_scripts/PrincipalComponents.R')
source('transcriptomic_analysis_scripts/ClusteringFunctions.R')

############################ Load in Data Files ################################

# CHANGE THIS DIRECTORY
dl<-"~/Documents/LEC_Time_Series_HTSeq_Counts"

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
names(lt)[2]<-"length"

###################### Annotate Genes in table 'lt' ##########################
ah<-AnnotationHub(localHub = T)
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

########################### Setup Master DGE List ############################
htseq_count<-hc_buildDataFrame(ds, ft)
rownames(ft)<-ft$sample
ft$sample<-paste(ft$genotype, paste(ft$hours.pcs,"H",sep=''), ft$batch,1:36, sep="_")
colnames(htseq_count)<-ft[colnames(htseq_count), 'sample']
row.names(ft)<-ft$sample

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
save(master, file = "LTS_DGEList.Rdata")


######################## Analyze Complete Data Set ###########################
# Build design matrix with batch coveriates
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

# Define Full Contrast Matrix.
# Note that each additional contrast increases overall sensetivity
cmat<-makeContrasts( 
  A=WT_6H-WT_0H, 
  B=WT_24H-WT_0H, 
  C=WT_48H-WT_0H, 
  D=WT_120H-WT_0H,
  E=FN_48H-FN_0H,
  `F`=FN_0H-WT_0H,
  G=FN_48H-WT_48H,
  H=B8_24H-B8_0H,
  I=B8_0H-WT_0H,
  J=B8_24H-WT_24H,
  levels = design
)

dge<-dge[filterByExpr(dge, design), , keep.lib.size=FALSE]

# Calculate size factors and Dispersion
dge<-calcNormFactors(dge)
dge<-estimateDisp(dge, design, robust = T)

# Fit Quasi Likelihood Model 
fit<-glmQLFit(dge, design)
qlf<-glmQLFTest(fit, contrast = cmat)     # Test each timepoint vs time zero
deg <-as.data.frame(topTags(qlf, n=Inf))  # Tabulate DEG from each time point

########## Tabulate WT Group Mean RPKM for Biological significance ###########
biosig<-data.frame(
  gene_id=row.names(rpkm(dge)),
  `0H`=apply(rpkm(
    dge[,dge$samples$genotype == "WT" & 
          dge$samples$interval == "0H"
        ]), 1, mean
  ),
  `6H`=apply(rpkm(
    dge[,dge$samples$genotype == "WT" & 
          dge$samples$interval == "6H"
        ]), 1, mean
  ),
  `24H`=apply(rpkm(
    dge[,dge$samples$genotype == "WT" & 
          dge$samples$interval == "24H"
        ]), 1, mean
  ),
  `48H`=apply(rpkm(
    dge[,dge$samples$genotype == "WT" & 
          dge$samples$interval == "48H"
        ]), 1, mean
  ),
  `120H`=apply(rpkm(
    dge[,dge$samples$genotype == "WT" & 
          dge$samples$interval == "120H"
        ]), 1, mean
  ), stringsAsFactors=F
)
colnames(biosig)<-gsub('X', '', colnames(biosig))
biosig<-biosig %>% filter(`0H` > 2 | `6H` > 2 | `24H` > 2 | `48H` > 2 | `120H` > 2) %>%
  filter(
    abs(`0H` - `6H`) > 2 |
      abs(`0H` - `24H`) > 2 |
      abs(`0H` - `48H`) > 2 |
      abs(`0H` - `120H`) > 2
  )
######################### Get Batch Dependent Genes ##########################
bta<-master
des<-model.matrix(~batch, bta$samples)
colnames(des)<-c(
  '(Intercept)', 'DNA1','DNA2'
)

# Define Full Contrast Matrix.
# Note that each additional contrast increases overall sensetivity

bta<-bta[filterByExpr(bta, des), , keep.lib.size=FALSE]

# Calculate size factors and Dispersion
bta<-calcNormFactors(bta)
bta<-estimateDisp(bta, design, robust = T)

# Fit Quasi Likelihood Model 
fit<-glmQLFit(bta, des)
qlf<-glmQLFTest(fit, coef=2:3)     # Test each timepoint vs time zero
bat <-as.data.frame(topTags(qlf, n=Inf))  # Tabulate DEG from each time point

bg<-(
  bat %>% filter(
    abs(logFC.DNA1) > 1 & 
      abs(logFC.DNA2) > 1 &
      FDR < 0.05)
)$gene_id

########################## Functionalize Filtering ###########################
#
#  The filterCPMmat function returns a matrix that has been filtered
#  to include only genes with an FDR < 0.05, a log fold change > than
#  a given threshold in any of the specified 'lfc_columns', that has an
#  average value (accross all samples) between 'min_log_cpm' and 'max_log_cpm'
#  and is included in a list of 'include_genes' if this list is provided. 
#
wt_samples<-dge$samples[dge$samples$genotype %in% c('WT'),'sample']
filterCPMmat<-function(
  object = dge, samples = wt_samples,
  lfc = 3, min_log_cpm = 0, max_log_cpm = 100, include_genes = biosig$gene_id, 
  lfc_columns = paste('logFC', LETTERS[1:4], sep="."), deg_table = deg
){
  ecpm<-cpm(object, log = T)        # Generate a gene x sample log CPM matrix
  
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

########### Filter by Variance; Plot Clusters and PCA for Wildtype ###########
dge$genes$Var<-apply(dge$counts[,wt_samples], 1, var)
varRanks <-c(10, 50, 100, 250, 500)           # Try different variance filters
for( v in varRanks){
  print(v)
  
  ecpm.filter <- cpm(
    dge[order(dge$genes$Var, decreasing=T)[1:v],wt_samples], log=T
  )
  
  # Plot results for TMM Normalized Data -- selected by 
  f1<-paste('plots/ECPM_Samples_Top_', v,'_Variance_Class.png', sep='')
  f2<-paste('plots/ECPM_Samples_Top_', v,'_Variance_Lab.png', sep='')
  f3<-paste('plots/ECPM_Samples_Top_', v,'_Variance_Clusters.png', sep='')
  png(f1, width=480, height=300)
  print(
    plotPrinComp(ecpm.filter, dge$samples[wt_samples,], groupCol=1, idCol=0))
  dev.off()
  
  png(f2, width=480, height=300)
  print(
    plotPrinComp(ecpm.filter, dge$samples[wt_samples,], groupCol=7, idCol=0))
  dev.off()
  
  png(f3, width=1200, height = 800)
  pheatmap(
    ecpm.filter, 
    annotation_col = dge$samples[wt_samples,c('interval', 'batch')],
    show_rownames = F, fontsize=20, cellwidth = 20
  )
  dev.off()
}

########## Filter by Fold Change; Plot Clusters and PCA for Wildtype #########
icol<-c(RColorBrewer::brewer.pal(5, 'YlGnBu'))
bcol<-c(RColorBrewer::brewer.pal(3, 'YlGnBu'))
mcol<-c(RColorBrewer::brewer.pal(3, 'YlGnBu'))
ccol<-colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100)
names(icol)<-levels(dge$samples$interval)
names(bcol)<-levels(as.factor(dge$samples$batch))
names(mcol)<-levels(as.factor(dge$samples$genotype))

key_genes<-c(
  Egr1="ENSMUSG00000038418",
  Grem1="ENSMUSG00000074934",
  Ptx3="ENSMUSG00000027832",
  Fn1="ENSMUSG00000026193",
  Actn1="ENSMUSG00000015143",
  Acta2="ENSMUSG00000035783",
  Col1a1="ENSMUSG00000001506",
  Atf3="ENSMUSG00000026628",
  Pitx3="ENSMUSG00000025229",
  Klf2="ENSMUSG00000055148",
  Emp1="ENSMUSG00000030208",
  Cxcl1="ENSMUSG00000029380",
  Cxcl2="ENSMUSG00000058427",
  Cxcl5="ENSMUSG00000029371",
  `Ppbp (Cxcl7)`="ENSMUSG00000029372",
  #Cxcl9="ENSMUSG00000029417",
  Cxcr2="ENSMUSG00000026180",
  #Cxcr4="ENSG00000121966",
  Runx1="ENSMUSG00000022952",
  Tnc="ENSMUSG00000028364",
  Fos="ENSMUSG00000021250",
  Ier2="ENSMUSG00000053560",
  FosB="ENSMUSG00000003545",
  Jun="ENSMUSG00000052684",
  JunB="ENSMUSG00000052837",
  Tgfb2="ENSMUSG00000039239",
  Tgfb1="ENSMUSG00000002603",
  Pdfgb="ENSMUSG00000000489",
  Pdgfrb="ENSMUSG00000024620",
  Asah1="ENSMUSG00000031591"
)

## Original Parameter List
# parameters<-list(
#   P1=list(lfc=3, minlc=0.1, maxlc=1, inc=NULL),
#   P2=list(lfc=3, minlc=0.1, maxlc=5, inc=NULL),
#   P3=list(lfc=3, minlc=0.1, maxlc=1, inc=biosig$gene_id),
#   P4=list(lfc=3, minlc=0.1, maxlc=5, inc=biosig$gene_id),
#   P5=list(lfc=7, minlc=0.5, maxlc=2, inc=NULL),
#   P6=list(lfc=7, minlc=0.5, maxlc=5, inc=NULL),
#   P7=list(lfc=7, minlc=0.5, maxlc=2, inc=biosig$gene_id),
#   P8=list(lfc=7, minlc=0.5, maxlc=5, inc=biosig$gene_id),
#   P9=list(lfc=1, minlc=0.5, maxlc=20, inc=NULL),
#   P10=list(lfc=5, minlc=0.5, maxlc=20, inc=NULL),
#   P11=list(lfc=1, minlc=0.5, maxlc=20, inc=biosig$gene_id),
#   P12=list(lfc=5, minlc=0.5, maxlc=20, inc=biosig$gene_id),
#   P13=list(lfc=7, minlc=0.5, maxlc=20, inc=NULL),
#   P14=list(lfc=7, minlc=0.5, maxlc=20, inc=biosig$gene_id),
#   P15=list(lfc=9, minlc=0.5, maxlc=20, inc=NULL),
#   P16=list(lfc=9, minlc=0.5, maxlc=20, inc=biosig$gene_id),
#   P17=list(lfc=2, minlc=0.5, maxlc=10, inc=NULL),
#   P18=list(lfc=5, minlc=0.5, maxlc=10, inc=NULL),
#   P19=list(lfc=2, minlc=0.5, maxlc=10, inc=biosig$gene_id),
#   P20=list(lfc=5, minlc=0.5, maxlc=10, inc=biosig$gene_id),
#   P21=list(lfc=2, minlc=1, maxlc=12, inc=NULL),
#   P22=list(lfc=2, minlc=2, maxlc=12, inc=NULL),
#   P23=list(lfc=2, minlc=1, maxlc=12, inc=biosig$gene_id),
#   P24=list(lfc=2, minlc=2, maxlc=12, inc=biosig$gene_id),
#   P25=list(lfc=3, minlc=1, maxlc=12, inc=NULL),
#   P26=list(lfc=3, minlc=2, maxlc=12, inc=NULL),
#   P27=list(lfc=3, minlc=1, maxlc=12, inc=biosig$gene_id),
#   P28=list(lfc=3, minlc=2, maxlc=12, inc=biosig$gene_id),
#   P29=list(lfc=4, minlc=0.1, maxlc=12, inc=NULL),
#   P30=list(lfc=4, minlc=0.2, maxlc=12, inc=NULL),
#   P31=list(lfc=4, minlc=0.1, maxlc=12, inc=biosig$gene_id),
#   P32=list(lfc=4, minlc=0.2, maxlc=12, inc=biosig$gene_id),
#   P33=list(lfc=3, minlc=1, maxlc=6, inc=NULL),
#   P34=list(lfc=3, minlc=2, maxlc=6, inc=NULL),
#   P35=list(lfc=3, minlc=1, maxlc=6, inc=biosig$gene_id),
#   P36=list(lfc=3, minlc=2, maxlc=6, inc=biosig$gene_id),
#   P37=list(lfc=2.9, minlc=0.5, maxlc=9, inc=biosig$gene_id),
#   P38=list(lfc=2.9, minlc=0.1, maxlc=9, inc=biosig$gene_id),
#   P39=list(lfc=2.9, minlc=0.5, maxlc=12, inc=biosig$gene_id),
#   P40=list(lfc=2.9, minlc=0.1, maxlc=12, inc=biosig$gene_id),
#   P41=list(lfc=2.5, minlc=0.5, maxlc=9, inc=biosig$gene_id),
#   P42=list(lfc=2.5, minlc=0.1, maxlc=9, inc=biosig$gene_id),
#   P43=list(lfc=2.5, minlc=0.5, maxlc=12, inc=biosig$gene_id),
#   P44=list(lfc=2.5, minlc=0.1, maxlc=12, inc=biosig$gene_id),
#   P45=list(
#     lfc=3, minlc=0.1, maxlc=5, inc=
#       setdiff(deg$gene_id, biosig$gene_id)
#   ),
#   P46=list(lfc=2.5, minlc=0.1, maxlc=12, inc=NULL),
#   P47=list(lfc=2.5, minlc=0.1, maxlc=9, inc=NULL),
#   P48=list(lfc=2.5, minlc=0.1, maxlc=8, inc=NULL),
#   P49=list(lfc=2.5, minlc=0.1, maxlc=7, inc=NULL),
#   P50=list(lfc=2.5, minlc=2, maxlc=20, inc=NULL),
#   P51=list(lfc=2.5, minlc=4, maxlc=20, inc=NULL),
#   P52=list(lfc=2.5, minlc=6, maxlc=20, inc=NULL)
# )

# evalClusters<-data.frame(
#   Params=character(),
#   lfc=numeric(),
#   minlc=numeric(),
#   maxlc=numeric(),
#   inclen=numeric(),
#   timeRand=numeric(),
#   batchRand=numeric(),
#   numGenes=numeric(),
#   keyGens=numeric(),
#   stringsAsFactors=F
# )
# 
# 


wt_params<-list(
  P1=list(lfc=3, minlc=0.1, maxlc=1, inc=NULL),
  P2=list(lfc=3, minlc=0.1, maxlc=5, inc=NULL),
  P44=list(lfc=2.5, minlc=0.1, maxlc=12, inc=biosig$gene_id),
  P45=list(
    lfc=3, minlc=0.1, maxlc=5, inc=
      setdiff(deg$gene_id, biosig$gene_id)
  ),
  P46=list(lfc=2.5, minlc=0.1, maxlc=12, inc=NULL),
  P49=list(lfc=2.5, minlc=0.1, maxlc=7, inc=NULL),
  P51=list(lfc=2.5, minlc=0, maxlc=20, inc=NULL),
  P52=list(lfc=2.5, minlc=0, maxlc=20, inc=biosig$gene_id)
)

for( pid in names(wt_params)){
  p<-wt_params[[pid]]
  # Get list of genes meeting key criteria
  m <-filterCPMmat(
    object = dge, samples=wt_samples, lfc=p[['lfc']],
    min_log_cpm = p[['minlc']], max_log_cpm = p[['maxlc']],
    include_genes = p[['inc']]
  )
  print(
    paste("pid:",pid,"size:",nrow(m),"range:", min(m), '-', max(m)))
  print(pid)
  print(nrow(m))
  print(paste("Number of Key genes:", length(intersect(key_genes, row.names(m)))))
  print(names(key_genes)[key_genes %in% intersect(key_genes, row.names(m))])
  # Plot results selected by fold change level
  f1<-paste(
    'plots/ECPM_Samples_Params_', pid,'_',nrow(m),
    '_Significant_genes_Class.png', sep=''
  )
  f2<-paste(
    'plots/ECPM_Samples_Params_', pid,'_',nrow(m),
    '_Significant_genes_Lab.png', sep=''
  )
  f3<-paste(
    'plots/WT_Samples_Params_', pid,'_genes_',
    nrow(m),'_Clustered.png', sep=''
  )
  f4<-gsub('Clustered.png', 'Matrix.txt', f3)
  
  # Plot Principal Components for Wildtype samples
  png(f1, width=480, height=300)
  print(
    plotPrinComp(m, dge$samples[wt_samples,], groupCol='interval', idCol=0))
  dev.off()
  
  png(f2, width=480, height=300)
  print(
    plotPrinComp(m, dge$samples[wt_samples,], groupCol=7, idCol=0))
  dev.off()
  
  # Plot clusters for wildtype animals
  png(f3, width=600, height = 700)
  pheatmap(
    m, 
    annotation = dge$samples[wt_samples,c('interval','batch')],
    show_rownames = F, show_colnames = T, color = ccol,
    annotation_colors = list(interval = icol,batch = bcol),
    treeheight_row = 0, fontsize = 20
  )
  dev.off()
  
  m<-as.data.frame(m)
  m$gene_id<-row.names(m)
  m<-m[,c('gene_id', wt_samples)]
  
  write.table(m, f4, sep="\t", row.names = F, quote=F)
  f5<-gsub('WT_', 'All_', f4)
  
  n<-as.data.frame(cpm(dge[row.names(m),],  log=T))
  s<-colnames(n)
  n$gene_id<-row.names(n)
  write.table(
    n[,c('gene_id', s)],
    f5, sep="\t", row.names = F, quote=F
  )
}

##################### Try and Identify Mutant Clusters #######################

b8_samples<-dge$samples[
  dge$samples$genotype %in% c('WT', 'B8') &
    dge$samples$batch == 'DNA1' &
    dge$samples$hours.pcs %in% c(0, 24),'sample'
  ]

fn_samples<-dge$samples[
  dge$samples$genotype %in% c('WT', 'FN') &
    dge$samples$batch == 'DBI' &
    dge$samples$hours.pcs %in% c(0, 48),'sample'
  ]

p<-wt_params[['P51']]
m <-filterCPMmat(
  object = dge, samples=wt_samples, lfc=p[['lfc']],
  min_log_cpm = p[['minlc']], max_log_cpm = p[['maxlc']],
  include_genes = p[['inc']]
)

n <-filterCPMmat(
  object = dge, samples=fn_samples, lfc=2,
  min_log_cpm = 1, max_log_cpm = 100,
  include_genes = NULL,
  lfc_columns = paste('logFC',c('F','G'), sep='.')
)

ix<-intersect(row.names(m), row.names(n))
f3<-paste(
  'plots/fibronectin_',nrow(n),"_Clustered.png", sep=""
)
png(f3, width=600, height = 800)
pheatmap(
  n, 
  annotation = dge$samples[b8_samples, c('genotype','interval')],
  show_rownames = F, show_colnames = T, color = ccol,
  annotation_colors = list(genotype = mcol, interval = icol),
  treeheight_row = 0, fontsize = 20
)
dev.off()

f3<-paste(
  'plots/fibronectin_',length(ix),'_of_',nrow(n),"_Clustered.png", sep=""
)
png(f3, width=600, height = 800)
pheatmap(
  n[ix,], 
  annotation = dge$samples[b8_samples, c('genotype','interval')],
  show_rownames = F, show_colnames = T, color = ccol,
  annotation_colors = list(genotype = mcol, interval = icol),
  treeheight_row = 0, fontsize = 20
)
dev.off()

n <-filterCPMmat(
  object = dge, samples=b8_samples, lfc=2,
  min_log_cpm = 1, max_log_cpm = 100,
  include_genes = NULL,
  lfc_columns = paste('logFC',c('I','J'), sep='.')
)
  
ix<-intersect(row.names(m), row.names(n))
f3<-paste(
  'plots/beta_8_',nrow(n),"_Clustered.png", sep=""
)
png(f3, width=600, height = 800)
pheatmap(
  n, 
  annotation = dge$samples[b8_samples, c('genotype','interval')],
  show_rownames = F, show_colnames = T, color = ccol,
  annotation_colors = list(genotype = mcol, interval = icol),
  treeheight_row = 0, fontsize = 20
)
dev.off()

f3<-paste(
  'plots/beta_8_',length(ix),'_of_',nrow(n),"_Clustered.png", sep=""
)
png(f3, width=600, height = 800)
pheatmap(
  n[ix,], 
  annotation = dge$samples[b8_samples, c('genotype','interval')],
  show_rownames = F, show_colnames = T, color = ccol,
  annotation_colors = list(genotype = mcol, interval = icol),
  treeheight_row = 0, fontsize = 20
)
dev.off()
############################ Plot Full Matrix ################################
m <-filterCPMmat(
  object = dge, samples=dge$samples$sample, lfc=3.5,
  min_log_cpm = 0, max_log_cpm = 20,
  include_genes = NULL, lfc_columns = grep("logFC", names(deg))
)

f3<-paste(
  'plots/All_Samples_lfc_3.5_any_contrast_',nrow(m),'_genes.png',
  sep=''
)
png(f3, width=700, height = 700)
pheatmap(
  m,
  annotation = dge$samples[, c('genotype','interval', 'batch')],
  show_rownames = F, show_colnames = T, color = ccol,
  annotation_colors = list(batch=bcol, genotype = mcol, interval = icol),
  treeheight_row = 0, fontsize = 10
)
dev.off()

f3<-gsub('_genes.png','_Matrix.txt', f3)
m<-as.data.frame(m)
m$gene_id<-row.names(m)
m<-m[,c('gene_id', setdiff(colnames(m), 'gene_id'))]
write.table(
  m, f3, sep="\t", quote=F, row.names = F
)

####### Write Out Full Count and RPKM Tables, and Filtered CPM Tables ########

# Raw Counts for Complete Expermient
m<-as.data.frame(master$counts)
s<-colnames(m)
m$gene_id<-row.names(m)
m<-m[,c('gene_id', s)]
write.table(
  m, "plots/Injury_Model_Raw_Counts_All_Genes.txt", row.names=F, 
  sep = '\t', quote=F
)

# Length normalized FPKM values (no TMM scaling)
r<-as.data.frame(rpkm(master))
s<-colnames(r)
r$gene_id<-row.names(r)
r<-r[,c('gene_id', s)]

write.table(
  r, "plots/Injury_Model_FPKM_All_Genes.txt", row.names=F, 
  sep = '\t', quote=F
)

# Sample ID's and Covariates
write.table(
  dge$samples, "plots/Sample_Metadata.txt", row.names=F, 
  sep = '\t', quote=F
)

# Gene ID's with Symbol, name and Length Used for RPKM normalization
write.table(
  dge$genes, "plots/Gene_Metadata.txt", row.names=F, 
  sep = '\t', quote=F
)

# edgeR TMM/cpm (log2 transformed) used for clustering (whole data set)
c<-cpm(dge, log=T)
write.table(
  c, "plots/Injury_Model_TMM_Normalized_Present_Genes.txt", row.names=F, 
  sep = '\t', quote=F
)




