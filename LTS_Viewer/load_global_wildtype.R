## Need Major Change, make this script pull data from synapse

# Set Up Environment
library(dplyr)
library(cluster)
library(reshape2)
library(tidyr)
library(ggplot2)


# Import DEG Table and FPKM matrix
load("data/LIRTS_master_dgelist.Rdata")
load("data/LIRTS_Master_DEG_Table.Rdata")
gene_meta <- read.csv(
  file="data/Mouse_Gene_Annotations.csv",
  header=T, row.names = 1
)
fpkm <- read.csv(
  file="data/Global_Wildtype_TMM-FPKM_Matrix.csv",
  header=T
)
plt<-fpkm %>% 
  pivot_longer(!gene_id, names_to="sample", values_to="fpkm") %>%
  mutate(fpkm=2^fpkm) %>%
  inner_join(
    master$samples, by=c("sample"="label")
  )

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
  res[[2]] %>%
    filter(Samples =="Global_Wildtype" & Group_1 !="single") %>%
    inner_join(gene_meta %>% select(gene_id, SYMBOL), by="gene_id") %>%
    select(SYMBOL, FDR ,tgt_gid=gene_id), 
  by=c(Target="SYMBOL")
) %>% filter(FDR < 0.05)  %>%
  group_by(Factor) %>%
  filter(n() > 10) %>%
  select(-FDR) %>%
  distinct(Factor, Target, tgt_gid) %>%
  group_by(Factor,Target, tgt_gid) %>%
  summarise(Rows=n()) %>%
  as.data.frame()

# Setp Annotations
annot<-plt %>% 
  select( sample, hours_pcs) %>% distinct() %>% as.data.frame()
row.names(annot)<-annot$sample




