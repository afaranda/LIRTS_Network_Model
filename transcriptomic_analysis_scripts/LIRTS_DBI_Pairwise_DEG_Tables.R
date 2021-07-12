##############################################################################
# File: LIRTS_DBI_Pairwise_DEG_Tables.R                                      #
# Purpose:                                                                   #
#      Construct differential gene expression (DEG) tables for pairwise      #
#      contrasts from the Runx1 Study transcriptome profiles.  Counts were   #
#      normalized (TMM) according to the specific contrast being analyzed,   #
#      such that only the subset of samples  participating in a given        #
#      were normalized together.  Abundance estimates averaged by group for  #
#      each gene, are reported as FPKM. FPKM values were calculated using    #
#      the exon-union length, after removing low-count genes (filterByExpr), #
#      using library sizes scaled by the TMM normalization factor.           #
#      DEG Tables capturing  differential expression estimates were          #
#      generated using an exact test and a Quasi-Likelihood F-Test           #
#      (in separate tables).                                                 # 
#                                                                            #
#      Labels for Sample Groups:                                             #
#        WT0H -- Zero Hours PCS, Cre Negative LEC (Intact Runx1)             #
#        WT24H -- 24 Hours PCS, Cre Negative LEC  (Intact Runx1)             #
#        R10H -- Zero Hours PCS, Cre Positive LEC (Excised Runx1)            #
#        R124H -- 24 Hours PCS, Cre Positive LEC (Excised Runx1)             #
#                                                                            #
#                                                                            #
#                                                                            #
#                                                                            #                                      
# Author: Adam Faranda                                                       #
# Created: June 2, 2021                                                      #
##############################################################################

############################ Setup Environment ###############################

# CHANGE THIS DIRECTORY
setwd("~/Documents/Runx1_Study")
#load("GeneLengthTable.Rdata")
library(dplyr)
library(cluster)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)
library(ggfortify)
library(synapser)
wd<-getwd()

synapser::synLogin()
syn_project <- "syn25831542"                 ## Synapse Project ID
syn_code_dir <- "syn25831545"     ## Synapse folder with code
syn_figures_dir <- "syn25831833"             ## Synapse folder with figures
local_data_dir <- paste0(wd,"/data")         ## Local data directory
deg_master <- data.frame()                   ## Table with all DE contrasts
deg_summary <- data.frame()                  ## Tabulate all DE summary stats

## Create Folder to store DEG tables from pairwise contrasts evaluated
## using edgeR's exactTest
syn_deg_dir <- Folder(
  name="DEG_Tables_Pairwise_Analysis",
  parent = syn_project
)
syn_deg_dir <- synStore(
  syn_deg_dir
)

# Add this script to the code dir
synapse_push <- File(
  path="scripts/Runx1_Pairwise_DEG_Tables.R",
  parent=syn_code_dir
)
synapse_push <- synStore(
  synapse_push
)

syn_act <- Activity(
  name="upload_analysis_results",
  description="upload analysis results"
)
syn_act$executed(synapse_push)


########################## Load in Master DGEList ############################

if(!file.exists('data/Runx1_master_dgelist.Rdata')){
  source('scripts/Runx1_Fetch_Count_Data.R')
} else {
  load('data/Runx1_master_dgelist.Rdata')
}

#####################  Setup group and contrast iterators ####################
# Define Groups As a list of the form: "Group1=c('S1', 'S2', 'S3')
# Each list element is the name of the Group and points to a
# Vector of sample labels, drop samples without a group. 
# A sample can only be assigned to one group.
groups=list()
for(g in unique(master$samples$group)){
  groups[[g]] <- which(master$samples$group == g)
}

# Define A list of named contrasts; each element points to a vector with
# a pair of group labels. Positive fold changes will be associated
# with the second group listed. 
contrasts=list(
  Wildtype_72vs0H=c('WT0H', 'WT72H'),
  Runx1_72vs0H=c('R10H', 'R172H'),
  Runx1vsWT_0H=c('WT0H', 'R10H'),
  Runx1vsWT_72H=c('WT72H', 'R172H')
)

######################## Evaluate Pairwise Contrasts #########################
for( c in names(contrasts)){

  ## Subset DGEList and drop unused covariate levels
  gr1<-groups[[contrasts[[c]][1]]]
  gr2<-groups[[contrasts[[c]][2]]]
  dge <- master[,c(gr1,gr2), keep.lib.sizes=F]
  dge$samples$group <- factor(
    dge$samples$group,
    levels=c(
      contrasts[[c]][1],
      contrasts[[c]][2]
    )
  )
  for(covariate in names(dge$samples)){
    if(is.factor(dge$samples[,covariate])){
      dge$samples[covariate] <- droplevels(dge$samples[,covariate])
    }
  }
  

  ## Drop undetectable genes, Normalize and Fit model
  design <- model.matrix(~ group, dge$samples)
  colnames(design) <- gsub("group","", colnames(design))
  dge <- dge[filterByExpr(dge, design),,keep.lib.sizes=F]
  dge <- calcNormFactors(dge)
  dge <- estimateDisp(dge, design, robust = T)
  
  # Calculate Group Average FPKM and join to genes table
  fpkm <- as.data.frame(
    rpkmByGroup(
      dge, gene.length = "eu_length"
    )[,contrasts[[c]]]
  )
  dge$genes <- cbind(dge$genes, fpkm[dge$genes$gene_id,])
  
  # Fit Robust Quasi-likelihood model
  fit <- glmQLFit(dge, design, robust=T)
  
  ## Validate Groups / Contrast assignments
  print(paste("Contrast:",c, " Groups: ", contrasts[[c]]))
  print(paste("Group 1:",paste(gr1, collapse=", ")))
  print(paste("Group 2:",paste(gr2, collapse=", ")))
  print("Design Matrix:")
  print(design)
  print(
    dge$samples[,c("sample", "label", "genotype", "interval")]
  )
  
  # Calculate Differential Expression Statistics
  deg.et<-as.data.frame(topTags(exactTest(dge, pair = contrasts[[c]]), n=Inf))
  deg.et$FDR <- qvalue::qvalue(deg.et$PValue)$qvalues
  names(deg.et) <- gsub(
    paste0("^",contrasts[[c]][1],"$"),
    "Avg1", names(deg.et)
  )
  names(deg.et) <- gsub(
    paste0("^",contrasts[[c]][2],"$"),
    "Avg2", names(deg.et)
  )
  deg.et$Group_1 <- contrasts[[c]][1]
  deg.et$Group_2 <- contrasts[[c]][2]
  
  
  deg.qt<-as.data.frame(topTags(glmQLFTest(fit,coef=2), n=Inf))
  deg.qt$FDR <- qvalue::qvalue(deg.qt$PValue)$qvalues
  names(deg.qt) <- gsub(
    paste0("^",contrasts[[c]][1],"$"),
    "Avg1", names(deg.qt)
  )
  names(deg.qt) <- gsub(
    paste0("^",contrasts[[c]][2],"$"),
    "Avg2", names(deg.qt)
  )
  deg.qt$Group_1 <- contrasts[[c]][1]
  deg.qt$Group_2 <- contrasts[[c]][2]
  
  
  # Push ExactTest Results to Synapse
  fn<-paste(
    "results/R1S_",c,"_Pairwise_Exact_Test_DEG.tsv", sep="")
  write.table(
    deg.et, fn, col.names =T, quote = F, sep="\t", row.names = F
  )
  syn_deg_file <- File(
    path=fn,
    parent=syn_deg_dir
  )
  syn_deg_file <- synStore(
    syn_deg_file
  )
  synSetProvenance(
    syn_deg_file,
    syn_act
  )
  
  # Push Quasi Likelihood Results to Synapse
  fn<-paste(
    "results/R1S_",c,"_Pairwise_QLFTest_DEG.tsv", sep="")
  write.table(
    deg.qt, fn, col.names=T, quote = F, sep="\t", row.names = F
  )
  syn_deg_file <- File(
    path=fn,
    parent=syn_deg_dir
  )
  syn_deg_file <- synStore(
    syn_deg_file
  )
  synSetProvenance(
    syn_deg_file,
    syn_act
  )
  
  ## Generate diagnostic plots. 
  png(
    paste0("results/R1S_",c,"_Pairwise_BCV_Plot.png")
  )
  plotBCV(dge)                                              # BCV Plot
  dev.off()
  synapse_push <- File(
    path=paste0("results/R1S_",c,"_Pairwise_BCV_Plot.png"),
    parent=syn_figures_dir
  )
  synapse_push <- synStore(
    synapse_push
  )
  synSetProvenance(
    synapse_push,
    syn_act
  )
  
  ## Print summary tables for DEG analysis
  deg_summary<-bind_rows(
    deg_summary,
    deg.et %>%
      # filter(abs(logFC) > 1 & FDR < 0.05) %>%
      summarize(
        Samples = "Pairwise",
        Contrast = c,
        Method="Exact", 
        Total = sum(abs(logFC) > 1 & FDR < 0.05), 
        Up = sum(logFC > 1 & FDR < 0.05), 
        Down = sum(logFC < -1 & FDR < 0.05),
        EU_Length_Bias = cor(
          logFC, log(eu_length,2), method = "spearman"
        ),
        EU_GC_Bias = cor(
          logFC, eu_gc, method = "spearman"
        ),
        PL_Length_Bias = cor(
          logFC, log(pl_length,2), method = "spearman"
        ),
        PL_GC_Bias = cor(
          logFC, pl_gc, method = "spearman"
        )
      ), 
    deg.qt %>%
      # filter(abs(logFC) > 1 & FDR < 0.05) %>%
      summarize(
        Samples = "Pairwise",
        Contrast = c,
        Method="QLF",
        Total = sum(abs(logFC) > 1 & FDR < 0.05), 
        Up = sum(logFC > 1 & FDR < 0.05), 
        Down = sum(logFC < -1 & FDR < 0.05),
        EU_Length_Bias = cor(
          logFC, log(eu_length,2), method = "spearman"
        ),
        EU_GC_Bias = cor(
          logFC, eu_gc, method = "spearman"
        ),
        PL_Length_Bias = cor(
          logFC, log(pl_length,2), method = "spearman"
        ),
        PL_GC_Bias = cor(
          logFC, pl_gc, method = "spearman"
        )
      )
  )
  print(deg_summary)
  
  ## Add DEG results to master DEG Table. 
  deg_master<-bind_rows(
    deg_master,
    deg.et %>% 
      tibble::remove_rownames() %>%
      select(
        gene_id, logFC, PValue,
        FDR, Avg1, Avg2
      ) %>% mutate(
        Group_1 = contrasts[[c]][1], 
        Group_2 = contrasts[[c]][2],
        Test="Exact", 
        Samples="Pairwise"
      ),
    deg.qt  %>% select(
        gene_id, logFC, PValue,
        FDR, Avg1, Avg2
      ) %>% 
      tibble::remove_rownames() %>% 
      mutate(
        Group_1 = contrasts[[c]][1], 
        Group_2 = contrasts[[c]][2],
        Test="QLF", 
        Samples="Pairwise"
      )
  )
  
  
  ## Generate Bias Plots (Using Exact Test)
  
  yl <-paste(                                    # Construct y axis label
    "Log2 Fold Change in", contrasts[[c]][2],          
    "vs", contrasts[[c]][1]
  )
  
  fn<-paste("results/EU_Length_Bias_in_Pairwise_",c,".png",sep='')
  mn<-paste("Length Bias in ",c,sep='')
  png(fn, width=6, height = 5, units="in", res=1200)
  plot(
    x=log(deg.et$eu_length,2), y=deg.et$logFC, main=mn,
    xlab = "Log2 Gene Length (exon-union)",
    ylab = yl,
    col = ifelse(
      deg.et[,"Avg1"] == 0,
      "red",
      ifelse(deg.et[,"Avg2"] == 0, "red", "black")
    )
  )
  
  
  # Add Correlation Coefficients to tests
  ct=cor.test(deg.et$logFC, log(deg.et$eu_length,2), method="spearman")
  rho=paste("rho:", round(ct$estimate,3))
  sig=paste("p value:", signif(ct$p.value,3), sep="")
  abline(lm(deg.et$logFC~log(deg.et$eu_length,2)), col="red", lwd=2)
  text(7,max(deg.et$logFC)-1,rho)
  text(7,max(deg.et$logFC)-3,sig)
  dev.off()
  
  synapse_push <- File(
    path=fn,
    parent=syn_figures_dir
  )
  synapse_push <- synStore(
    synapse_push
  )
  synSetProvenance(
    synapse_push,
    syn_act
  )
  
  fn<-paste("results/EU_GC_Bias_in_Pairwise_",c,".png",sep='')
  mn<-paste("GC Bias in ",c,sep='')
  png(fn, width=6, height = 5, units="in", res=1200)
  plot(
    x=deg.et$eu_gc, y=deg.et$logFC, main=mn,
    xlab = "Fractional GC content (exon-union)",
    ylab = yl,
    col = ifelse(
      deg.et[,"Avg1"] == 0,
      "red",
      ifelse(deg.et[,"Avg2"] == 0, "red", "black")
    )
  )
  
  ct=cor.test(deg.et$logFC, log(deg.et$eu_gc,2), method="spearman")
  rho=paste("rho:", round(ct$estimate,3))
  sig=paste("p value:", signif(ct$p.value,3), sep="")
  abline(lm(deg.et$logFC~deg.et$eu_gc), col="red", lwd=2)
  text(0.35,max(deg.et$logFC)-1,rho)
  text(0.35,max(deg.et$logFC)-3,sig)
  dev.off()
  
  synapse_push <- File(
    path=fn,
    parent=syn_figures_dir
  )
  synapse_push <- synStore(
    synapse_push
  )
  synSetProvenance(
    synapse_push,
    syn_act
  )
}

## Push summary and master DEG table to synapse
syn_cols <- list(                                 
  Column(name="Samples", columnType='STRING', maximumSize =100), 
  Column(name="Contrast", columnType='STRING', maximumSize =100),
  Column(name="Method", columnType='STRING', maximumSize =100),
  Column(name='Total', columnType='INTEGER'),
  Column(name="Up", columnType='INTEGER'),
  Column(name='Down', columnType='INTEGER'),
  Column(name='EU_Length_Bias', columnType='DOUBLE'),
  Column(name="EU_GC_Bias", columnType='DOUBLE'),
  Column(name="PL_Length_Bias", columnType='DOUBLE'),
  Column(name="PL_GC_Bias", columnType='DOUBLE')
)
syn_schema <- Schema(
  name="Pairwise_DEG_Summary_Stats",
  columns=syn_cols,
  parent=syn_project
)
syn_table <- Table(
  schema = syn_schema,
  values=deg_summary
)

syn_table <- synStore(syn_table)

syn_act <- Activity(
  name="upload_deg_results",
  description="upload results from differential expression analysis"
)

syn_act$executed(synapse_push)

synSetProvenance(
  syn_table$tableId,
  syn_act
)

syn_cols <- list(                                 
  Column(name="gene_id", columnType='STRING', maximumSize =100), 
  Column(name='logFC', columnType='DOUBLE'),
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
  name="Pairwise_DEG_Master_Table",
  columns=syn_cols,
  parent=syn_project
)
syn_table <- Table(
  schema = syn_schema,
  values=deg_master
)

syn_table <- synStore(syn_table)

syn_act <- Activity(
  name="upload_deg_results",
  description="upload results from differential expression analysis"
)

syn_act$executed(synapse_push)

synSetProvenance(
  syn_table$tableId,
  syn_act
)



