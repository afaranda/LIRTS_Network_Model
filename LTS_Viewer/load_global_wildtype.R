
# Set Up Environment
library(dplyr)
library(cluster)
library(reshape2)
library(tidyr)
library(ggplot2)

# Import master DGE_List object
load("data/LTS_DGEList.Rdata")

# Run Global Wildtype Analysis
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

plt<-fpkm %>% 
  pivot_longer(!gene_id, names_to="sample", values_to="fpkm") %>%
  inner_join(dge$samples, by="sample")

# Load Transcription Factor tables
hs.tf<-read.table(
  "data/trrust_rawdata.human.tsv", sep="\t", quote="", header=F
)
# Quick and Dirty Human to Mouse Conversion
hs.tf$V1<-paste(
  substr(hs.tf$V1,1,1),
  tolower(substr(hs.tf$V1,2,nchar(hs.tf$V1))),
  sep=""
)

hs.tf$V2<-paste(
  substr(hs.tf$V2,1,1),
  tolower(substr(hs.tf$V2,2,nchar(hs.tf$V2))),
  sep=""
)

mm.tf<-read.table(
  "data/trrust_rawdata.mouse.tsv", sep="\t", quote="", header=F
)

tf_tgt<-bind_rows(
  hs.tf %>% select(Factor=V1,Target=V2),
  mm.tf %>% select(Factor=V1,Target=V2)
) %>% distinct() %>% inner_join(
  deg %>% select(SYMBOL, FDR ,tgt_gid=gene_id), by=c(Target="SYMBOL")
) %>% filter(FDR < 0.05)  %>%
  group_by(Factor) %>%
  filter(n() > 10)

# Setp Annotations
annot<-plt %>% 
  select( sample, interval) %>% distinct() %>% as.data.frame()
row.names(annot)<-annot$sample




