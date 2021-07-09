########################         HEADER BLOCK      ###########################
#  File:    Wrap_edgeR_Functions.R                                           #
#  Purpose: Define wrapers to streamline automated edgeR analysis            #
#  Created: Feb 2, 2021                                                      #
#  Author:  Adam Faranda                                                     #
#                                                                            #
##############################################################################

#########      Load Libraries and import expression data          ############
library(edgeR)
library(dplyr)
library(ggplot2)
setwd("~/Desktop/FCM_Methods_Paper_Analysis")
wd<-getwd()
########  Wrapper function to fit a QLF model to a design matrix  ############
processByDesign <- function(y, design, rob=T){
  y <- y[filterByExpr(y, design),,keep.lib.sizes=F]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design, robust = rob)
  
  fpkm<-rpkm(y, normalized.lib.sizes = T, gene.length = "eu_length")
  for(g in unique(y$samples$group)){
    s<-row.names(y$samples[y$samples$group == g,])
    y$genes[paste0(g,"_Avg_FPKM")]<-apply(fpkm[,s],1,mean)
  }
  
  fit <- glmQLFit(y, design, robust = rob)
  return(list(dge=y, fit=fit))
}

######## Split and Process DGE List based on a set of sample groups ##########
subsetDGEListByGroups<-function(y, groups=c("GR1", "GR2"), norm="TMM"){
  
  # Get samples associated with two groups
  group_subset<-list()
  for(gr in groups){
    group_subset[[gr]] <- row.names(
      y$samples %>% dplyr::filter(group == gr)
    )
  }
  print(group_subset)
  # Extract Subset DGEList 
  y<-y[, unlist(group_subset)]
  
  # Relevel Factors as needed
  for(covariate in colnames(y$samples)){
    if(is.factor(y$samples[covariate])){
      y$samples[covariate] <-droplevels(y$samples[covariate])
    }
  }
  
  # Fix the first group as the reference level / Intercept
  y$samples$group <- relevel(y$samples$group, groups[1])
  
  if(norm == "TMM"){
    design<-model.matrix(~group, y$samples)
    colnames(design) <-gsub("group","",colnames(design))
    y <- y[filterByExpr(y, design),,keep.lib.sizes=F]
    y <- calcNormFactors(y)
    y <- estimateDisp(y, design, robust=T)
    
    fpkm<-rpkm(y, normalized.lib.sizes = T, gene.length = "eu_length")
    for(g in unique(y$samples$group)){
      s<-row.names(y$samples[y$samples$group == g,])
      y$genes[paste0(g,"_Avg_FPKM")]<-apply(fpkm[,s],1,mean)
    }
    
    fit <- glmQLFit(y, design, robust = T)
    return(
      list(
        dge=y, fit=fit, design=design
      )
    )
  } else if (norm == "VOOM") {
    print(levels(y$samples$group))
    design<-model.matrix(~group, y$samples)
    colnames(design) <- gsub("(\\(|\\))", "", colnames(design))
    colnames(design) <-gsub("group","",colnames(design))
    y <- y[filterByExpr(y, design),,keep.lib.sizes=F]
    y <- calcNormFactors(y)
    
    # NOTE -- FPKM Calculation Fails to account for VOOM precision Weights !!!
    # THIS SHOULD BE FIXED if VOOM is intended for use !!!
    fpkm<-rpkm(y, normalized.lib.sizes = T, gene.length = "eu_length")
    for(g in unique(y$samples$group)){
      s<-row.names(y$samples[y$samples$group == g,])
      y$genes[paste0(g,"_Avg_FPKM")]<-apply(fpkm[,s],1,mean)
    }

    v <- voom(y, design=design)
    return(list(v=v, design=design))
  } else if (norm == "RUVs"){
    
    design<-model.matrix(~group, y$samples)
    colnames(design) <-gsub("group","",colnames(design))
    y <- y[filterByExpr(y, design),,keep.lib.sizes=F]
    controls <- row.names(y)
    
    # Assemble "Samples" matrix required by RUVs
    differences <- matrix(-1,
        nrow = length(levels(y$samples$group)), 
        ncol = max(table(y$samples$group))
    )
    for(gri in 1:length(levels(y$samples$group))){
      level <- levels(y$samples$group)[gri]
      
      s <- which(y$samples$group == level)
      for(j in 1:length(s)){
        differences[gri, j] <- s[j]
      }
    }
    print("Difference Sets for RUVs")
    print(differences)
    
    # Estimate nuisance parameter W_1
    seq <- newSeqExpressionSet(counts = y$counts)
    ruvs <- RUVs(
      seq, cIdx = row.names(y),
      k =1, scIdx=differences
    )
    y$samples$W_1 <- pData(ruvs)[row.names(y$samples), "W_1"]
    design<-model.matrix(~group + W_1, y$samples)
    colnames(design) <-gsub("group","",colnames(design))
    
    y <- y[filterByExpr(y, design),,keep.lib.sizes=F]
    y <- calcNormFactors(y)
    y <- estimateDisp(y, design, robust=T)
    
    # NOTE -- FPKM Calculation Fails to account for RUV adjustment !!!
    # THIS may need to be FIXED if RUV is intended for use !!!
    fpkm<-rpkm(y, normalized.lib.sizes = T, gene.length = "eu_length")
    for(g in unique(y$samples$group)){
      s<-row.names(y$samples[y$samples$group == g,])
      y$genes[paste0(g,"_Avg_FPKM")]<-apply(fpkm[,s],1,mean)
    }
    
    fit <- glmQLFit(y, design, robust = T)
    return(
      list(
        dge=y, fit=fit, ruvs=ruvs, design=design
      )
    )
    
  } else if (norm == "RUVr"){
    
    # Fit a typical edgeR model
    z <- y
    design<-model.matrix(~group, z$samples)
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
    design<-model.matrix(~group + W_1, y$samples)
    colnames(design) <-gsub("group","",colnames(design))
    
    y <- y[filterByExpr(y, design),,keep.lib.sizes=F]
    y <- calcNormFactors(y)
    y <- estimateDisp(y, design, robust=T)
    
    # NOTE -- FPKM Calculation Fails to account for RUV adjustment !!!
    # THIS may need to be FIXED if RUV is intended for use !!!
    fpkm<-rpkm(y, normalized.lib.sizes = T, gene.length = "eu_length")
    for(g in unique(y$samples$group)){
      s<-row.names(y$samples[y$samples$group == g,])
      y$genes[paste0(g,"_Avg_FPKM")]<-apply(fpkm[,s],1,mean)
    }
    
    fit <- glmQLFit(y, design, robust = T)
    return(
      list(
        dge=y, fit=fit, ruvr=ruvr, design=design
      )
    )
    
  } else {
    return(y)
  }
}

#### Define function to generate DEG Tables from a fit or DGEList ############
genDegTable<-function(y, group1, group2, design){
  print(
    paste(
      group1, "Samples:",
      paste(
        y$samples %>% filter(group == group1) %>%row.names(),
        collapse = ", "
      )
    )
  )
  print(
    paste(
      group2, "Samples:",
      paste(
        y$samples %>% filter(group == group2) %>%row.names(),
        collapse = ", "
      )
    )
  )
  
  if(class(y) == "DGEGLM"){
    cn<-sapply(
      colnames(design), 
      function(n){
        ifelse(
          n == group1, -1, 
          ifelse(n == group2, 1, 0
          )
        )
      }
    )
    print(cn)
    return(
      as.data.frame(
        topTags(
          glmQLFTest(y, contrast=cn), n=Inf
        )
      ) %>% 
        dplyr::mutate(
          Test = "QLFTest",
          Group_1 = group1, 
          Group_2 = group2
        ) %>%
        dplyr::select(
          Test, Group_1, Group_2,
          gene_id, logFC, logCPM, PValue, FDR, 
          Avg1 = as.name(paste0(group1, "_Avg_FPKM")),
          Avg2 = as.name(paste0(group2, "_Avg_FPKM"))
        )
    )
  } else if(class(y) == "DGEList") {
    return(
      as.data.frame(
        topTags(
          exactTest(y, pair=c(group1, group2)), n=Inf
        )
      ) %>% 
        dplyr::mutate(
          Test = "ExactTest",
          Group_1 = group1, 
          Group_2 = group2
        ) %>%
        dplyr::select(
          Test, Group_1, Group_2,
          gene_id, logFC, logCPM, PValue, FDR, 
          Avg1 = as.name(paste0(group1, "_Avg_FPKM")),
          Avg2 = as.name(paste0(group2, "_Avg_FPKM"))
        )
    )
  } else if(class(y) == "EList"){
    fit <- lmFit(y, design=design)
    cft <- contrasts.fit(
      fit, contrasts=makeContrasts(
        contrasts = paste(group2, group2, sep="-"),
        levels=colnames(design)
      )
    )
    ebs <- eBayes(cft)
    return(
      as.data.frame(
        topTable(cft, n=Inf)
      ) %>% 
        dplyr::mutate(
          Test = "LimmaVoom",
          Group_1 = group1, 
          Group_2 = group2
        ) %>%
        dplyr::select(
          Test, Group_1, Group_2,
          gene_id, logFC, logCPM=AveExper, PValue=P.Value, FDR=adj.P.Val, 
          Avg1 = as.name(paste0(group1, "_Avg_FPKM")),
          Avg2 = as.name(paste0(group2, "_Avg_FPKM"))
        )
    )
  }
}