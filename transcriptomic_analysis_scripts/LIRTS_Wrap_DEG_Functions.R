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
library(ggfortify)

############  Wrapper function to fit models to a design matrix  #############
process_edgeR_ByDesign <- function(y, genes=NULL, design, rob=T, norm="TMM"){
  if(is.null(genes)){
    y<-y[filterByExpr(y, design), ,keep.lib.sizes=F] # Drop low features
  } else {
    y<-y[genes, ,keep.lib.sizes=F] # Or keep selected features
  }
  y <- calcNormFactors(y, method = norm)
  y <- estimateDisp(y, design, robust = rob)
  
  fpkm<-rpkm(y, normalized.lib.sizes = T, gene.length = "eu_length")
  for(g in unique(y$samples$group)){
    s<-row.names(y$samples[y$samples$group == g,])
    y$genes[paste0(g,"_Avg_FPKM")]<-apply(fpkm[,s],1,mean)
  }
  
  fit <- glmQLFit(y, design, robust = rob)
  return(list(dge=y, fit=fit))
}

process_voom_ByDesign <- function(y, genes=NULL, design, rob=T, norm="TMM"){
  if(is.null(genes)){
    y<-y[filterByExpr(y, design), ,keep.lib.sizes=F] # Drop low features
  } else {
    y<-y[genes, ,keep.lib.sizes=F] # Or keep selected features
  }
  y <- calcNormFactors(y, method = norm)
  y <- estimateDisp(y, design, robust = rob)
  
  # Estimate precision weights and 
  v <- voom(y, design=design)
  
  fpkm_matrix <- apply(2^v$E, 2, function(x) x/(v$genes$eu_length/1000))
  
  # Calculate average values for each group
  for(g in unique(v$targets$group)){
    s<-row.names(v$targets[v$targets$group == g,])
    v$genes[paste0(g,"_Avg_FPKM")]<-apply(fpkm_matrix[,s],1,mean)
  }
  
  # Fit linear models to genes
  fit <- lmFit(v, design=design)
  return(list(dge=v, fit=fit))
}
############################ Plot Diagnostics ################################
## Generate diagnostic plots. 
diagnostic_plots <- function(
  dge=dge, color_attrib="group", shape_attrib=NULL, 
  respath="LIRTS_DEG_Analysis_results", prefix="DBI_Wildtype"
){
  # Colorblind friendly pallatte from 
  # https://bconnelly.net/
  
  colors=c(
    "#000000", "#E69F00", "#56B4E9", "#009E73",
    "#0072B2", "#D55E00", "#CC79A7"
  )

  png(
    paste(
      respath,"/",prefix,"_BCV_Plot.png", sep=""
    )
  )  
  plotBCV(dge)                                              # BCV Plot
  dev.off()
  
  # Plot sample projections in first two Principal Components
  pca <-prcomp(t(cpm(dge, log=T)), scale=T)
  fn <- paste(
    respath,"/",prefix,"_PCA_Plot.png", sep=""
  )
  
  if(is.null(color_attrib)){
    ggsave(
      fn,
      autoplot(
        pca, data=dge$samples,size=3
      ), height =4, width=6
    )
  } else if(is.null(shape_attrib)){
    ggsave(
      fn,
      autoplot(
        pca, data=dge$samples,
        colour=color_attrib, size=3
      ) + scale_color_manual(values=colors), 
      height =4, width=6
    )
  } else {
    ggsave(
      fn,
      autoplot(
        pca, data=dge$samples,
        colour=color_attrib, shape=shape_attrib, size=3
      ) + scale_color_manual(values=colors), 
      height =4, width=6
    )
  }
}


### Generate a BCV plot, highlighting Flagged genes in red
plot_BCV_flag<- function (
  y, xlab = "Average log CPM", ylab = "Biological coefficient of variation", 
  pch = 16, cex = c(0.2, 0.5), col.common = "grey", col.trend = "blue", 
  cols.tagwise = c("red","black"), flag=FALSE, flag_labels=c("Tagwise - Key Gene","Tagwise Other Gene"),...) 
{
  if (!is(y, "DGEList")) 
    stop("y must be a DGEList.")
  A <- y$AveLogCPM
  if (is.null(A)) 
    A <- aveLogCPM(y$counts, offset = getOffset(y))
  disp <- getDispersion(y)
  if (is.null(disp)) 
    stop("No dispersions to plot")
  if (attr(disp, "type") == "common") 
    disp <- rep_len(disp, length(A))
  plot(A, sqrt(disp), xlab = xlab, ylab = ylab, type = "n", 
       ...)
  labels <- cols <- lty <- pt <- NULL
  if (!is.null(y$tagwise.dispersion)) {
    points(A, sqrt(y$tagwise.dispersion), pch = pch, #cex = cex, 
           col = ifelse(
             flag,
             "red", "black"
           ),
           cex=ifelse(
             flag,
             cex[1], cex[2]
           )
    )
    labels <- c(labels, flag_labels)
    cols <- c(cols, cols.tagwise)
    lty <- c(lty, -1, -1)
    pt <- c(pt, pch, pch)
  }
  if (!is.null(y$common.dispersion)) {
    abline(h = sqrt(y$common.dispersion), col = col.common, 
           lwd = 2)
    labels <- c(labels, "Common")
    cols <- c(cols, col.common)
    lty <- c(lty, 1)
    pt <- c(pt, -1)
  }
  if (!is.null(y$trended.dispersion)) {
    o <- order(A)
    lines(A[o], sqrt(y$trended.dispersion)[o], col = col.trend, 
          lwd = 2)
    labels <- c(labels, "Trend")
    cols <- c(cols, col.trend)
    lty <- c(lty, 1)
    pt <- c(pt, -1)
  }
  legend("topright", legend = labels, lty = lty, pch = pt, 
         pt.cex = cex + c(0,0.5), lwd = 2, col = cols)
  invisible()
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

process_selected_features <- function(
  dge, design, genes=NULL, counts=NULL, 
  prefix="Global_Wildtype_Top_Genes",
  color_attrib="hours_pcs", 
  shape_attrib = "batch"
){
  if(!is.null(counts)){
    dge$counts <- counts
  }
  obj <- process_edgeR_ByDesign(y=dge, genes=genes, design=design)
  diagnostic_plots(
    obj$dge, prefix = prefix, 
    color_attrib = color_attrib,
    shape_attrib = shape_attrib
  )
  fpkm<-as.data.frame(
    rpkm(
      obj$dge, gene.length="eu_length", log=T,
      normalized.lib.sizes = T
    )
  )
  colnames(fpkm) <- obj$dge$samples$label
  
  ## Plot samplewise correlation matrix after excluding genes with 
  ## a strong batch response (DBI vs DNA1, DNA2 or DNA3)
  cm <- cor(fpkm, method = "spearman")
  annot <- obj$dge$samples[,c("hours_pcs", "batch", "label")]
  annot <- annot %>%
    tibble::remove_rownames()%>%
    tibble::column_to_rownames("label")
  
  pheatmap(
    cm, annotation_col =annot,
    filename = paste0(
      "LIRTS_DEG_Analysis_results/",
      prefix, "_cormat.png"
    ),height = 6, width=8
  )
  
  ## Save FPKM matrix for genes with a significant batch dependent logFC
  ## less than 2 (four fold difference) (DBI vs DNA1, DNA2 or DNA3)
  fpkm$gene_id<-row.names(fpkm)
  write.csv(
    fpkm, row.names = F,
    paste0(
      "LIRTS_DEG_Analysis_results/",
      prefix,"_TMM-FPKM_Matrix.csv"
    )
  )
  return(obj)
}

#### Define function to generate DEG Tables from a fit or DGEList ############
genPairwiseDegTable<-function(y, group1, group2, design){
  print(
    paste(
      group1, "Samples:",
      paste(
        y$samples %>% filter(group == group1) %>% pull("label"),
        collapse = ", "
      )
    )
  )
  print(
    paste(
      group2, "Samples:",
      paste(
        y$samples %>% filter(group == group2) %>% pull("label"),
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
  } else if(class(y) == "MArrayLM"){
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

## Generate a DEG table using coefficients to a design matrix
genDesignCoefDegTable<-function(y, design, coef, group_labels){
  
  # Print Experimental Design and samples associated with each coeffient
  print(design)
  for(i in c(1,coef)){
    print(colnames(design)[i])
    print(y$samples[which(design[,i] == 1),"label"])
  }
  
  if(class(y) == "DGEGLM"){
    print(coef)
    return(
      as.data.frame(
        topTags(
          glmQLFTest(y, coef=coef), n=Inf
        )
      ) %>% 
        dplyr::mutate(
          Test = "QLFTest",
          Group_1 = group_labels[1], 
          Group_2 = group_labels[2],
          Avg1 = 2^((logCPM - (logFC/2)) - log2(eu_length/1000)),
          Avg2 = 2^((logCPM + (logFC/2)) - log2(eu_length/1000))
        ) %>%
        dplyr::select(
          Test, Group_1, Group_2,
          gene_id, logFC, logCPM, PValue, FDR, 
          Avg1, Avg2
        )
    )
  } else if(class(y) == "DGEList") {
    print("Analysis Reqires a fitted model")
    return(NULL)

  } else if(class(y) == "MArrayLM"){
    ebs <- eBayes(y)
    return(
      x <- as.data.frame(
        topTable(ebs, n=Inf, coef=coef)
      ) %>% 
        dplyr::mutate(
          Test = "LimmaVoom",
          Group_1 = group_labels[1], 
          Group_2 = group_labels[2],
          Avg1 = 2^((AveExpr - (logFC/2)) - log2(eu_length/1000)),
          Avg2 = 2^((AveExpr + (logFC/2)) - log2(eu_length/1000))
        ) %>%
        dplyr::select(
          Test, Group_1, Group_2,
          gene_id, logFC, logCPM=AveExpr, PValue=P.Value, FDR=adj.P.Val, 
          Avg1, Avg2
        )
    )
  }
}



# Iterate over a set of contrasts, generate DEG Tables and
#  DEG Summary tables for each contrast. 
iterate_edgeR_pairwise_contrasts <- function(
  dge, fit, cntmat=cntmat, df=df, design=design,
  deg=deg,respath="LIRTS_DEG_Analysis_results",
  prefix="DBI"
){
  print(cntmat)
  for (c in colnames(cntmat)) {             
    # Run Exact Tests
    pair<-c(
      names(cntmat[,c])[cntmat[,c] == -1], 
      names(cntmat[,c])[cntmat[,c] == 1]
    )
    print(pair)

    deg.et<-genPairwiseDegTable(dge, pair[1], pair[2], design)
    deg.qt<-genPairwiseDegTable(fit, pair[1], pair[2], design)
    
    deg <- bind_rows(
      deg,
      deg.et %>%
        mutate(Samples=prefix) %>%
        tibble::remove_rownames(),
      deg.qt %>%
        mutate(Samples=prefix) %>%
        tibble::remove_rownames()
    )
    
    # Save DEG Tables
    fn<-paste(
      respath,"/",prefix,"_",c,"_Exact_Test_DEG.tsv", sep=""
    )
    write.table(
      deg.et, fn, col.names =T, quote = F, sep="\t", row.names = F
    )
    
    fn<-paste(
      respath,"/",prefix,"_",c,"_QLFTest_DEG.tsv", sep=""
    )
    write.table(
      deg.qt, fn, col.names=T, quote = F, sep="\t", row.names = F
    )
    
    dg<-degSummary(                         # Generate Exact Test Summary tables
      deg.et,
      lfc = 'logFC',
      fdr = 'FDR', 
      Avg1 = pair[1],
      Avg2 = pair[2]
    )
    
    dg$contrast<-c
    dg$test<-"Exact Test"
    df<-bind_rows(df, dg)
    
    dg<-degSummary(                         # Generate QLF Test Summary tables 
      deg.qt,
      lfc = "logFC",
      fdr = 'FDR', 
      Avg1 = pair[1],
      Avg2 = pair[2]
    )
    dg$contrast<-c
    dg$test<-"QLFTest"
    df<-bind_rows(df, dg)
    
  }
  return(list(df, deg))
}


# Iterate over a set of design matrix coefficients, 
# generate DEG Tables, DEG Summary tables
iterate_edgeR_design_coefficients <- function(
  dge, fit, coefs=c(2,3), df=df, design=design, deg=deg,
  respath="LIRTS_DEG_Analysis_results", prefix="DBI", 
  group_label_list=list(c("ctrl", "trt"), c("mut", "wt"))
){
  print(group_label_list)
  for (c in 1:length(coefs)) {             
    deg.qt<-genDesignCoefDegTable(
      fit, design, coef=coefs[c], 
      group_labels = group_label_list[[c]]
    )
    
    deg <- bind_rows(
      deg, deg.qt %>%
        mutate(Samples=prefix) %>%
        tibble::remove_rownames()
    )
    
    # Save DEG Tables
    fn<-paste(
      respath,"/",prefix,"_",c,"_QLFTest_DEG.tsv", sep=""
    )
    write.table(
      deg.qt, fn, col.names=T, quote = F, sep="\t", row.names = F
    )
    
    dg<-degSummary(                         # Generate QLF Test Summary tables 
      deg.qt,
      lfc = "logFC",
      fdr = 'FDR', 
      Avg1 = "Avg1",
      Avg2 = "Avg2"
    )
    dg$contrast<-colnames(design)[coefs[c]]
    dg$test<-"QLFTest"
    df<-bind_rows(df, dg)
    
  }
  return(list(df, deg))
}


## Generate Length Bias Plot
# fn<-paste(
#   respath,"/",prefix,"_",c,"_Length_Bias.png", sep=""
# )
# 
# mn<-paste("Length Bias in ",c,sep='')
# png(fn, width=6, height = 5, units="in", res=1200)
# plot(
#   x=log(deg.et$eu_length,2), y=deg.et$logFC, main=mn,
#   xlab = "Log2 Gene Length (exon-union)",
#   ylab = yl,
#   pch = ifelse(
#     deg.et[,"Avg1"] == 0,
#     1,
#     ifelse(deg.et[,"Avg2"] == 0, 1, 16)
#   ),
#   col=ifelse(abs(deg.et$logFC)> 1 & deg.et$FDR < 0.05, "red","black")
# )
# 
# # Add Correlation Coefficients to tests
# ct=cor.test(deg.et$logFC, log(deg.et$eu_length,2), method="spearman")
# rho=paste("rho:", round(ct$estimate,3))
# sig=paste("p value:", signif(ct$p.value,3), sep="")
# abline(lm(deg.et$logFC~log(deg.et$eu_length,2)), col="red", lwd=2)
# text(7,max(deg.et$logFC)-1,rho)
# text(7,max(deg.et$logFC)-3,sig)
# dev.off()
# }
# 
# ## Generate GC Bias Plot
# fn<-paste(
#   respath,"/",prefix,"_",c,"_GC_Bias.png", sep=""
# )
# 
# mn<-paste("GC Bias in ",c,sep='')
# png(fn, width=6, height = 5, units="in", res=1200)
# plot(
#   x=deg.et$eu_gc, y=deg.et$logFC, main=mn,
#   xlab = "Fractional GC content (exon-union)",
#   ylab = yl,
#   pch = ifelse(
#     deg.et[,"Avg1"] == 0,
#     1,
#     ifelse(deg.et[,"Avg2"] == 0, 1, 16)
#   ),
#   col=ifelse(abs(deg.et$logFC)> 1 & deg.et$FDR < 0.05, "red","black")
# )
# 
# ct=cor.test(deg.et$logFC, log(deg.et$eu_gc,2), method="spearman")
# rho=paste("rho:", round(ct$estimate,3))
# sig=paste("p value:", signif(ct$p.value,3), sep="")
# abline(lm(deg.et$logFC~deg.et$eu_gc), col="red", lwd=2)
# text(0.35,max(deg.et$logFC)-1,rho)
# text(0.35,max(deg.et$logFC)-3,sig)
# dev.off()
