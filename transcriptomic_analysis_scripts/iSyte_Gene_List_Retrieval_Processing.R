# Purpose: Prepare Gene Lists for manual retrieval of iSyte Lens enrichment scores
# Created: June 11, 2019
# Author: Adam Faranda

# Setup Environment
library(dplyr)
library(reshape2)
library(openxlsx)




wd<-getwd()
outDir<-paste(wd,'DiffExp', sep="/")
iSyteDir<-paste(wd, 'iSyte_Data', sep="/")

# Import data tables
setwd(outDir)
degDataFiles<-list.files(pattern="expressedTags-all.txt$")


# Function that converts short names to descriptive filenames

descName<-function(x){
	x<-gsub('WT_0_Hour', 'Wildtype-0-Hour', x)
	x<-gsub('WT_48_Hour', 'Wildtype-48-Hour', x)
	x<-gsub('FN_0_Hour', 'FNcKO-0-Hour', x)
	x<-gsub('FN_48_Hour', 'FNcKO-48-Hour', x)
	
}

ag_master<-data.frame(stringsAsFactors = F)
ss_master<-data.frame(stringsAsFactors = F)
for (dataFile in degDataFiles){
	contrast<-gsub("_expressedTags-all.txt$","", dataFile)

	# Load in DEGs
	dg<-read.table(dataFile, sep="\t", quote="", header=T, stringsAsFactors = F)
	samples<-names(dg)[grep("_GeneCount.cpm", names(dg))]

	print(paste("############ FILTER DEG: ",contrast, "##########"))	

	# Add Columns for Filtering based on statistical Criteria
	print(paste("CPM Present Cols:", paste(samples, collapse=" "), sep= " "))
	dg$Cpm.Present<-apply(dg[, samples], 
		1, function(x) sum(x > 1) >=2
	)

	# Setup for RPKM Filtering Criteria
	samples<-gsub("_GeneCount.cpm", "_RPKM", samples)
	print(paste("RPKM Filtering Cols:", paste(samples, collapse=" "), sep=" "))
	dg$Avg1<-apply(dg[,samples[1:3]], 1, mean)
	dg$Avg2<-apply(dg[,samples[4:6]], 1, mean)
	
	dg$Group_1<-unlist(strsplit(contrast,"_vs_"))[1]
	dg$Group_2<-unlist(strsplit(contrast,"_vs_"))[2]
	
	#names(dg)[grep("^Avg", names(dg))]<-paste(unlist(strsplit(contrast,"_vs_")), "_Average_RPKM", sep="")
	dg$RPKM_gt2_Either_Cond<-apply(dg[,grep("^Avg", names(dg))], 1, function(x) (x[1] > 2) |(x[2] > 2) )
	dg$RPKM_Diff_gt2<-apply(dg[,grep("^Avg", names(dg))], 1, function(x) abs(x[1]-x[2]) > 2)
	
	# Get Statistically Significant Genes
	ds<-dg %>% filter(abs(logFC) > 1 & FDR < 0.05) # & Cpm.Present == T) -- CPM Filter already applied

	# Get Biologically significant genes based on the following criteria	
	db<-dg %>% filter(abs(logFC) > 1 & FDR < 0.05 & RPKM_gt2_Either_Cond == T & RPKM_Diff_gt2 == T)
	

	# Reorder Columns For readability, drop sample specific columns
	cols<-c("GeneID", 
	        "ensembl_gene_id", 
	        "description", 
	        "FC", "logFC", 
	        "FDR", "Group_1", 
	        "Group_2", "Avg1", "Avg2"
	       )
	ag_master<-bind_rows(
	  ag_master,
	  dg<-dg[,cols]
	)
	ss_master<-bind_rows(
	  ss_master,
	  ds<-ds[,cols]
	)

}

# Validate GeneID as primary key for each contrast (NOT FOR GLOBAL TABLE)

print("########## Confirm Unique Gene IDs in each Contrast (Stat)#########")
for( i in unique(ss_master$Group_1)){
  for(j in unique(ss_master$Group_2)){
    val<-ss_master %>% filter(Group_1 ==i, Group_2==j)
    if(nrow(val) > 0){
      print(
        paste(i,j,
              "Rows:", nrow(val),
              "Unique_Genes:", length(unique(val$GeneID))      
        )
      )
    }
  }
}

# Partition the complete list of unique Genes into blocks of 500
ugl=unique(ss_master$GeneID)
psize<-500
nparts<-1+(length(ugl) %/% psize)
partitions=list()
for( i in 1:nparts){
  g=(1+(i-1)*psize):(i*psize)
  g<-g[g<=length(ugl)]
  partitions[[i]]<-ugl[g]
}

# Write Partitions to columns in a spreadsheet, manually retrieve data for each partition
setwd(iSyteDir)
fn<-"iSyte_GeneLists.xlsx"
wb<-createWorkbook()
addWorksheet(wb,"GeneLists")
cells<-expand.grid(row=1:psize, col=1:nparts)
tx<-createStyle(numFmt="TEXT")
addStyle(wb, 1, rows=cells$row, cols=cells$col, style=tx)

for(i in 1:length(partitions)){
  writeData(wb,1,partitions[[i]], i,1)
}
saveWorkbook(wb, fn, overwrite =T )

# Import and process data from each partition
setwd(iSyteDir)
p<-read.table("6_iSyTE_Pvalue.txt", header=T, stringsAsFactors = F)

iSyteCleaner<-function(x, type="Pvalue"){
  q<-melt(x, id.vars = "Symbol")
  
  # Strip "Dev" tag out of melted column headers
  q$variable<-gsub("Dev","",as.character(q$variable))
  q$variable<-gsub(".PVal.", "", q$variable)
  # Strip the platform outy of melted headers to get the time interval
  q$Interval<-gsub(
    'affy430', '', gsub(
      'affy430a', '', gsub(
        'lumi1', '', gsub(
          'lumi2', '', 
          q$variable 
        )
      )
    )
  )
  
  q$Platform<-apply(
    q[,c('variable', 'Interval')], 1, 
    function(x)gsub(x[2], '', x[1])
  )
  q<-q[,c('Symbol', 'Platform', 'Interval', 'value')]
  names(q)[grep('value', names(q))] <-type
  
  q
}

fc<-data.frame()
for(i in list.files(pattern="FoldChange")){
  p<-read.table(i, header=T, stringsAsFactors = F)
  q<-iSyteCleaner(p, type = "FoldChange")
  fc<-bind_rows(fc, q)
}

pv<-data.frame()
for(i in list.files(pattern="Pvalue")){
  p<-read.table(i, header=T, stringsAsFactors = F)
  q<-iSyteCleaner(p, type = "Pvalue")
  pv<-bind_rows(pv, q)
}

# Helper function convert fold change to Log2 fold change
logify<-function(x, base=2){
  if(x >= 0){
    return (log(x, base))
  }
  else if( x < 0){
    return( log( 1/abs(x), base ))
  }
}

iSyte<-inner_join(pv, fc, by=c('Symbol', 'Platform', 'Interval'))
iSyte<-iSyte %>% 
  filter(Platform != "" & Interval != "Rank" & FoldChange !="-")
iSyte$Pvalue<-as.numeric(iSyte$Pvalue)
#iSyte$FoldChange<-as.numeric(iSyte$FoldChange)
iSyte$FoldChange<-sapply(as.numeric(iSyte$FoldChange), logify)

# Prepare Lens Enrichment Figures for various time points and contrasts
query<-function(
  isy=iSyte, deg=ss_master, 
  tp="P0", pl="affy430", 
  gr1="WT_0_Hour", gr2="WT_48_Hour",
  mle = 1
){
  isy<-isy %>% 
    filter(Platform == pl & Interval == tp &Pvalue < 0.05) %>%
    rename(
      LensEnrichment = FoldChange, 
      LensPvalue = Pvalue, 
      GeneID = Symbol
    ) %>%
    filter(abs(LensEnrichment) > mle)
  
  deg<-deg %>% filter(Group_1 == gr1 & Group_2 == gr2)
  print(nrow(deg))
  out<-inner_join(deg, isy, by ='GeneID')
  out
}


qplot(x=logFC, y=LensEnrichment, data=query(tp="P56", mle = 0)) +
  labs(y="Log2 Fold Change in the lens at P56", x = "Log2 Fold Change 0 vs 48 Hours Wildtype")


q1<-query(deg=ag_master %>% filter(FDR < 0.05), tp="P56", mle = 0)
q2<-query(deg=ag_master, tp="P56", mle = 0, gr1="FN_0_Hour", gr2 = "FN_48_Hour" )
  

x<-inner_join(
  q1, 
  q2 %>% select(GeneID, FN_0v48_logFC = logFC),
  by="GeneID"
)

qplot(x=logFC, y=LensEnrichment, color=FN_0v48_logFC, data=x) +
  labs(y="Log2 Fold Change in the lens at P56", 
       x = "Log2 Fold Change 0 vs 48 Hours Wildtype",
       color = "Log2FC 0 vs 48 Hours FNcKO"
  )







print("########## Confirm Unique Gene IDs in  each Contrast (Bio)#########")

# Get Overlapping Genes for each meaningful comparison between pairwise contrasts
comps<-list(
  `WT0vFN0_AND_WT0vWT48`=c("WT_0_Hour_vs_FN_0_Hour",  "WT_0_Hour_vs_WT_48_Hour"),
  `WT0vFN0_AND_FN0vFN48`=c("WT_0_Hour_vs_FN_0_Hour","FN_0_Hour_vs_FN_48_Hour"),
  `WT0vFN0_AND_WT48vFN48`=c("WT_0_Hour_vs_FN_0_Hour",  "WT_48_Hour_vs_FN_48_Hour"),
  `WT48vFN48_AND_WT0vWT48`=c("WT_48_Hour_vs_FN_48_Hour",  "WT_0_Hour_vs_WT_48_Hour"),
  `WT48vFN48_AND_FN0vFN48`=c("WT_48_Hour_vs_FN_48_Hour",  "FN_0_Hour_vs_FN_48_Hour"),
  `FN0vFN48_AND_WT0vWT48`=c("WT_0_Hour_vs_WT_48_Hour",  "FN_0_Hour_vs_FN_48_Hour")
)
# 
# for(comp in names(comps)){
#     c<-comps[[comp]]
#     
#     d1<-ss_master %>% filter(	
#           Group_1 == unlist(strsplit(c[1],"_vs_"))[1] & 
#           Group_2 == unlist(strsplit(c[1],"_vs_"))[2]
#     )
#     print(paste("First Contrast:", c[1], "Rows:", nrow(d1)))
#     
#     
#     d2<-ss_master %>% filter(	
#         Group_1 == unlist(strsplit(c[2],"_vs_"))[1] & 
#         Group_2 == unlist(strsplit(c[2],"_vs_"))[2]
#     )
#     
#     print(paste("Second Contrast:", c[2], "Rows:", nrow(d2)))
#    
#     
#     genes<-intersect(d1$GeneID, d2$GeneID)
#     print(paste("Genes In Common:", length(genes)))
#     
#     allResults<-inner_join(
#       d1[d1$GeneID %in% genes,],
#       d2[d2$GeneID %in% genes,],
#       by="GeneID",
#       suffix=c(paste("_",c[1], sep=''), paste("_",c[2], sep=''))
#     )
#     
#     degTable.G1<-data.frame(
#       criteria = c("Statistically Significant", "Biologically Significant"),
#       Total = c(
#         nrow(d1 %>% filter(abs(logFC) > 1, FDR < 0.05)), 
#         nrow(d1 %>% 
#                filter(Avg1 > 2  | Avg2 > 2) %>%
#                filter(abs(Avg1 - Avg2) > 2) %>%
#                filter(abs(logFC) > 1, FDR < 0.05)
#         )
#       ),
#       Up =  c(
#         nrow(d1 %>% filter(logFC > 0, FDR < 0.05)), 
#         nrow(d1 %>% 
#                filter(Avg1 > 2  | Avg2 > 2) %>%
#                filter(abs(Avg1 - Avg2) > 2) %>%
#                filter(logFC > 0, FDR < 0.05)
#         )
#       ),
#       Down =  c(
#         nrow(d1 %>% filter(logFC < 0, FDR < 0.05)), 
#         nrow(d1 %>% 
#                filter(Avg1 > 2  | Avg2 > 2) %>%
#                filter(abs(Avg1 - Avg2) > 2) %>%
#                filter(logFC < 0, FDR < 0.05)
#         )
#       )
#     )
#     
#     
#     degTable.G2<-data.frame(
#       criteria = c("Statistically Significant", "Biologically Significant"),
#       Total = c(
#         nrow(d2 %>% filter(abs(logFC) > 1, FDR < 0.05)), 
#         nrow(d2 %>% 
#                filter(Avg1 > 2  | Avg2 > 2) %>%
#                filter(abs(Avg1 - Avg2) > 2) %>%
#                filter(abs(logFC) > 1, FDR < 0.05)
#         )
#       ),
#       Up =  c(
#         nrow(d2 %>% filter(logFC > 0, FDR < 0.05)), 
#         nrow(d2 %>% 
#                filter(Avg1 > 2  | Avg2 > 2) %>%
#                filter(abs(Avg1 - Avg2) > 2) %>%
#                filter(logFC > 0, FDR < 0.05)
#         )
#       ),
#       Down =  c(
#         nrow(d2 %>% filter(logFC < 0, FDR < 0.05)), 
#         nrow(d2 %>% 
#                filter(Avg1 > 2  | Avg2 > 2) %>%
#                filter(abs(Avg1 - Avg2) > 2) %>%
#                filter(logFC < 0, FDR < 0.05)
#         )
#       )
#     )
#     
#     print(degTable.G1)
#     print(degTable.G2)
#   
# 
#     d1<-d1[,setdiff(colnames(d1), c("Group_1", "Group_2"))]
#     d2<-d2[,setdiff(colnames(d2), c("ensembl_gene_id", "description", "Group_1", "Group_2"))]
#     
#     
#     d1<-d1[,setdiff(colnames(d1), c("Group_1", "Group_2"))]
#     d2<-d2[,setdiff(colnames(d2), c("ensembl_gene_id", "description", "Group_1", "Group_2"))]
#     
#     allResults<-inner_join(
#       d1[d1$GeneID %in% genes,],
#       d2[d2$GeneID %in% genes,],
#       by="GeneID"
#     )
#     
#     bioResults<-allResults %>%
#       filter(Avg1.x > 2 | Avg2.x > 2) %>%
#       filter(abs(Avg1.x - Avg2.x) > 2) %>%
#       filter(Avg1.y > 2 | Avg2.y > 2) %>%
#       filter(abs(Avg1.y - Avg2.y) > 2)
#      
#     
#     names(allResults)<-gsub("Avg1.x",
#          paste(unlist(strsplit(c[1],"_vs_"))[1], "_Avg_RPKM.g1", sep=""), 
#          names(allResults)
#     )
#     names(allResults)<-gsub("Avg2.x",
#          paste(unlist(strsplit(c[1],"_vs_"))[2], "_Avg_RPKM.g1", sep=""), 
#          names(allResults)
#     )
#     names(allResults)<-gsub("Avg1.y",
#          paste(unlist(strsplit(c[2],"_vs_"))[1], "_Avg_RPKM.g2", sep=""), 
#          names(allResults)
#     )
#     names(allResults)<-gsub("Avg2.y",
#          paste(unlist(strsplit(c[2],"_vs_"))[2], "_Avg_RPKM.g2", sep=""), 
#          names(allResults)
#     )
# 
#     
#     names(bioResults)<-gsub("Avg1.x",
#                             paste(unlist(strsplit(c[1],"_vs_"))[1], "_Avg_RPKM.g1", sep=""), 
#                             names(bioResults)
#     )
#     names(bioResults)<-gsub("Avg2.x",
#                             paste(unlist(strsplit(c[1],"_vs_"))[2], "_Avg_RPKM.g1", sep=""), 
#                             names(bioResults)
#     )
#     names(bioResults)<-gsub("Avg1.y",
#                             paste(unlist(strsplit(c[2],"_vs_"))[1], "_Avg_RPKM.g2", sep=""), 
#                             names(bioResults)
#     )
#     names(bioResults)<-gsub("Avg2.y",
#                             paste(unlist(strsplit(c[2],"_vs_"))[2], "_Avg_RPKM.g2", sep=""), 
#                             names(bioResults)
#     )
#     
#     tables<-list(
#       `Stat Sig Intersection`=allResults,
#       `Bio Sig Intersection`=bioResults,
#       `SS Down Group_1 Up Group_2`=allResults %>% filter(logFC.x < 0,logFC.y > 0 ),
#       `SS Down Group_1 Down Group_2`=allResults %>% filter(logFC.x < 0,logFC.y < 0 ),
#       `SS Up Group_1 Up Group_2`=allResults %>% filter(logFC.x > 0,logFC.y > 0 ),
#       `SS Up Group_1 Down Group_2`=allResults %>% filter(logFC.x > 0,logFC.y < 0 ),
#       `BS Down Group_1 Up Group_2`=bioResults %>% filter(logFC.x < 0,logFC.y > 0 ),
#       `BS Down Group_1 Down Group_2`=bioResults %>% filter(logFC.x < 0,logFC.y < 0 ),
#       `BS Up Group_1 Up Group_2`=bioResults %>% filter(logFC.x > 0,logFC.y > 0 ),
#       `BS Up Group_1 Down Group_2`=bioResults %>% filter(logFC.x > 0,logFC.y < 0 )
#     )
#       
#     for(t in names(tables)){
#       names(tables[[t]])<-gsub(".x", paste("_", c[1], sep=""), names(tables[[t]]))
#       names(tables[[t]])<-gsub(".y", paste("_", c[2], sep=""), names(tables[[t]]))
#       #print(names(tables[[t]]))
#       print(nrow(tables[[t]]))
#     }
#     names(tables)<-gsub("Group_1", 
#                         unlist(strsplit(comp,"_AND_"))[1], 
#                         names(tables)
#     )
#     names(tables)<-gsub("Group_2", 
#                         unlist(strsplit(comp,"_AND_"))[2], 
#                         names(tables)
#     )
#     print(names(tables))
#     
#     wb<-loadWorkbook("Comparisons.xlsx")
#     writeData(wb, 1, degTable.G1, startCol=2, startRow=32,colNames=F)
#     writeData(wb, 1, degTable.G2, startCol=2, startRow=37,colNames=F)
#     
#     stat.Intersect<-data.frame(
#       U<-c(nrow(tables[[5]]), nrow(tables[[3]])),
#       D<-c(nrow(tables[[6]]), nrow(tables[[4]]))
#     )
#     writeData(wb, 1, stat.Intersect, startCol=3, startRow=42,colNames=F)
#     
#     biol.Intersect<-data.frame(
#       U<-c(nrow(tables[[9]]), nrow(tables[[7]])),
#       D<-c(nrow(tables[[10]]), nrow(tables[[8]]))
#     )
#     writeData(wb, 1, biol.Intersect, startCol=3, startRow=46,colNames=F)
#   
#     tx<-createStyle(numFmt="@")
#     for(i in names(tables)){
#       addWorksheet(wb,i)
#       print(nrow(tables[[i]]))
#       cells<-expand.grid(row=nrow(tables[[i]]), col=grep("GeneID", names(tables[[i]])))
#       addStyle(wb, i, rows=cells$row, cols=cells$col, style=tx)
#       writeData(wb, i, tables[[i]])
#     }
#     #saveWorkbook(wb, paste(comp, "Comparison.xlsx", sep="_"), overwrite=T)
#     
# }
print(sessionInfo())
setwd(wd)
#rm(list=ls())

