library(edgeR)
library(dplyr)
gene_counts <- function(
  means=c(1), sizes=c(1), 
  samples_per_group=c(1)
){
  obs <- numeric()
  if(
    length(means)!=length(sizes) | length(sizes) !=length(samples_per_group)
  ){
    print("settings don't match")
  } else {
    for(i in 1:length(means)){
      obs <- c(
        obs, rnbinom(
          n=samples_per_group[i],
          mu=means[i], 
          size=sizes[i]
        )
      )
    }
    return(obs)
  }
}

balanced_param_list <- list(
  list(
    n=c(5,5),
    m=c(10,10),
    s=c(2,2)
  ),
  list(
    n=c(5,5),
    m=c(100,100),
    s=c(10,10)
  ),
  list(
    n=c(5,5),
    m=c(10,10),
    s=c(5,2)
  ),
  list(
    n=c(5,5),
    m=c(100,100),
    s=c(10,80)
  ),
  list(
    n=c(5,5),
    m=c(100,10),
    s=c(10,2)
  )
)
mat <- t(
  sapply(
    sample(balanced_param_list, replace=T,size=5000),
    function(x) gene_counts(
      means=x[['m']], 
      sizes=x[['s']], 
      samples_per_group = x[['n']]
    )
  )
)
dge <- DGEList(
  counts=mat,
  group=factor(rep(c("ctl", "trt"), each=5))
)

dge <- calcNormFactors(dge)
dge <- estimateDisp(dge)

et <- exactTest(dge)
tt <- topTags(et) %>% as.data.frame()

tt$Avg1 <- tt$logCPM - (tt$logFC / 2)
tt$Avg2 <- tt$logCPM + (tt$logFC / 2)

cpm(dge)[row.names(tt),]


