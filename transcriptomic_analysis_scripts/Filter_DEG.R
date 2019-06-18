# Purpose: Merge RPKM values into pairwise comparison tables
# Created: March 5, 2019
# Author: Adam Faranda

# Setup Environment
library(dplyr)
library(openxlsx)
#library(org.Mm.eg.db)

wd<-getwd()
outDir<-paste(wd,'DiffExp', sep="/")


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

for (dataFile in degDataFiles){
	contrast<-gsub("_expressedTags-all.txt$","", dataFile)

	# Load in DEGs
	dg<-read.table(dataFile, sep="\t", quote="", header=T)
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
	
	names(dg)[grep("^Avg", names(dg))]<-paste(unlist(strsplit(contrast,"_vs_")), "_Average_RPKM", sep="")
	dg$RPKM_gt2_Either_Cond<-apply(dg[,grep("_Average_RPKM", names(dg))], 1, function(x) (x[1] > 2) |(x[2] > 2) )
	dg$RPKM_Diff_gt2<-apply(dg[,grep("_Average_RPKM", names(dg))], 1, function(x) abs(x[1]-x[2]) > 2)
	
	# Reorder Columns For readability, drop merge column
	cols<-colnames(dg)
        reor<-c("GeneID", "ensembl_gene_id", "description")
	cols<-setdiff(c(reor, setdiff(cols, reor)),"uniqueId")
	dg<-dg[,cols]
	
	# Get Statistically Significant Genes
	ds<-dg %>% filter(abs(logFC) > 1 & FDR < 0.05) # & Cpm.Present == T) -- CPM Filter already applied
	
	# Get Biologically significant genes based on the following criteria	
	db<-dg %>% filter(abs(logFC) > 1 & FDR < 0.05 & RPKM_gt2_Either_Cond == T & RPKM_Diff_gt2 == T)


	print(
		paste(
			"Rowcount before join: ", nrow(df),
			"Rowcount after join: ", nrow(dg), 
			"Unique GeneID: ", length(unique(dg$GeneID))
		)
	)
	
	print(
		paste(
			"Rowcount Stat Sig.: ", nrow(ds), 
			"Stat Sig Unique GeneID: ", length(unique(ds$GeneID))
		)
	)
	
	print(
		paste(
			"Rowcount Biol Sig.: ", nrow(db), 
			"Biol Sig Unique GeneID: ", length(unique(db$GeneID))
		)
	)
	print("")
	degTable<-data.frame(
		Criteria=c("Statistically Significant", "Biologically Significant"),
		Total=c(nrow(ds), nrow(db)),
		Upregulated=c(
			nrow(ds[ds$logFC > 1,]),
			nrow(db[db$logFC > 1,])
		),

		Downregulated=c(
			nrow(ds[ds$logFC < -1,]),
			nrow(db[db$logFC < -1,])
		)
	)

	print(degTable)
	print(list.files("../", pattern=descName(contrast)))
	

	# Build Excel Workbook
	wb<-loadWorkbook(
		paste("..", 
			list.files("../", pattern=descName(contrast)), 
			sep="/"
		)
	)
	tx<-createStyle(numFmt="TEXT")
	writeData(wb, 1, degTable, startCol=2, startRow=24,colNames=F) # Add DEG counts to main page
	
	tables<-list(
		`All Present Genes`=dg,
		`Statistically Significant`=ds,
		`Biologically Significant`=db
	)

	for(i in names(tables)){
	      tables[[i]] <- tables[[i]][, setdiff(names(tables[[i]]), "Cpm.Present")] # Drop "Cpm.Present" Column
	      addWorksheet(wb,i)
	      print(nrow(tables[[i]]))
	      cells<-expand.grid(row=nrow(tables[[i]]), col=grep("GeneID", names(tables[[i]])))
	      addStyle(wb, i, rows=cells$row, cols=cells$col, style=tx)
	      writeData(wb, i, tables[[i]])
	}
	
	saveWorkbook(wb, paste(descName(contrast), "Differentially_Expressed_Genes.xlsx", sep="_"), overwrite=T)		
}

print(sessionInfo())
