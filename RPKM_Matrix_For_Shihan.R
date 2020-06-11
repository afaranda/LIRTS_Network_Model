setwd('~/Documents/LEC_Time_Series')
load('LTS_DGEList.Rdata')
source('transcriptomic_analysis_scripts/BuildDataMatrix.R')
source('transcriptomic_analysis_scripts/PreprocessingFunctions.R')
source('transcriptomic_analysis_scripts/PrincipalComponents.R')
source('transcriptomic_analysis_scripts/ClusteringFunctions.R')
source('transcriptomic_analysis_scripts/Overlap_Comparison_Functions.R')

dge<-master[,
            master$samples$batch == 'DNA1' &
              master$samples$genotype %in% c("WT", "B8") &
              master$samples$interval %in% c('0H', '24H')
            ]
# Reorder grouping factor, drop unused levels
dge$samples$interval<-droplevels(
  factor(
    paste(dge$samples$hours.pcs, 'H', sep=''),
    levels = c('0H', '6H', '24H','48H', '120H')
  )
)

dge$samples$genotype<-droplevels(
  factor(
    dge$samples$genotype,
    levels = c('WT', 'FN', 'B8')
  )
)

dge$samples$group<-gsub("_DNA1_[0-9].","",dge$samples$sample)

# Add Average FPKMs to gene table in master DGEList (based on full library)
for(g in unique(dge$samples$group)){
  cn<-paste(g, "Avg", sep="_")
  s<-dge$samples[dge$samples$group == g, 'sample']
  dge$genes[cn] <- apply(rpkm(dge)[,s], 1, mean)   
}

# Drop genes with low detection
dge<-dge[filterByExpr(dge, group=dge$samples),,keep.lib.sizes=F]

## Drop remaining duplicates due to lack of biological significance
d<-duplicated(dge$genes$SYMBOL)
u<-unique(dge$genes[d,'gene_id'])
dge<-dge[!dge$genes$gene_id %in% u,]

dge<-calcNormFactors(dge)
dge<-estimateDisp(dge, robust=T)

dg1<-uniqueMaxLfc(as.data.frame(
  topTags(exactTest(dge, pair = c('WT_0H', 'WT_24H')), n=Inf)
), idc='SYMBOL')
row.names(dg1)<-dg1$gene_id

dg2<-uniqueMaxLfc(as.data.frame(
  topTags(exactTest(dge, pair = c('WT_0H', 'B8_0H')), n=Inf)
), idc='SYMBOL')
row.names(dg2)<-dg2$gene_id

dg3<-uniqueMaxLfc(as.data.frame(
  topTags(exactTest(dge, pair = c('B8_0H', 'B8_24H')), n=Inf)
), idc='SYMBOL')
row.names(dg3)<-dg3$gene_id

dg4<-uniqueMaxLfc(as.data.frame(
  topTags(exactTest(dge, pair = c('WT_24H', 'B8_24H')), n=Inf)
), idc='SYMBOL')
row.names(dg4)<-dg4$gene_id


dge$genes$WT0HvsWT24H_logFC<-dg1[dge$genes$gene_id, 'logFC']
dge$genes$WT0HvsWT24H_FDR<-dg1[dge$genes$gene_id, 'FDR']

dge$genes$WT0HvsB80H_logFC<-dg2[dge$genes$gene_id, 'logFC']
dge$genes$WT0HvsB80H_FDR<-dg2[dge$genes$gene_id, 'FDR']

dge$genes$B80HvsB824H_logFC<-dg3[dge$genes$gene_id, 'logFC']
dge$genes$B80HvsB824H_FDR<-dg3[dge$genes$gene_id, 'FDR']

dge$genes$WT24HvsB824H_logFC<-dg4[dge$genes$gene_id, 'logFC']
dge$genes$WT24HvsB824H_FDR<-dg4[dge$genes$gene_id, 'FDR']

r<-as.data.frame(rpkm(dge, log=T))
r$gene_id<-row.names(r)

r<-inner_join(
  dge$genes,
  r, by='gene_id'
)
write.csv(r, "RPKM_Matrix_For_Shihan.txt")
