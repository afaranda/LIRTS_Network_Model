library(edgeR)
library(biomaRt)
library(dplyr)
options(echo=F)

# Enter Working Directory and Load Raw Data
#setwd('/Users/afaranda/Desktop/LEC_Time_Series')
wd<-getwd()
outDir<-paste(wd,'DiffExp',sep='/')

df <- read.table("GeneCount.tsv", header=T, stringsAsFactors =F)

# Define Primary Key for Raw Data -- GeneID is not sufficient
# Use a composite of GeneID, Chromosome, Coding Length
df$PrimaryKey <- paste(df$GeneID, df$Chr, df$CodingLength)
print(paste("Validate Primary Key, Row count: ", nrow(df), ", Unique keys: ", length(unique(df$PrimaryKey)), sep=""))

# Helper Function to get version number from 'ensembl_gene_id_version'
getversion<-function(ens){
	v<-strsplit(ens, "\\.")[[1]][2]
	return(as.integer(v))
}

# Run Biomart Query to fetch descriptions and Ensembl ID's
query<-list.files(pattern="biomaRt_query")
if (length(query) > 0) {
   load(query[1])
} else {
	# Annotate Differentially Expressed Genes best possible ENSEMBL ID
		
	# 'attr': list of attributes to retrieve
	attr<-c("ensembl_gene_id", "external_gene_name", "description",
	        "ensembl_gene_id_version", "chromosome_name", "gene_biotype"
	        )
	# 'fc': attribute used for query by gene ID
	fc<-'external_gene_name'
	print(paste("Using attributes: ", attr, sep=""))
	print(paste("Searching on: ", fc, sep=""))
	
	# 'gq': list of unique 'GeneID' submitted as biomart query
	gq<-unique(df$GeneID)

	# Actual biomaRt query
	mart<-useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
	result<-getBM(mart=mart, attributes=attr, filters=fc, values=gq)
	
	# Get set of unique gene symbols for all duplicated rows
	dg<-unique(result[duplicated(result$external_gene_name),'external_gene_name'])
	
	# Filter biomaRt query to remove duplicate annotations
	result<-bind_rows(
		# All rows with a unique external gene name
		result %>% filter(!(external_gene_name %in% dg)),
		
		result %>% 
		# Restrict filter to rows with  duplicated 'external_gene_name' 
		filter(external_gene_name %in% dg) %>%
		
		# Group by external_gene_name and filter rows where chromosome_name is not
		# an actual chromosome; 1-18, X or Y
		group_by(external_gene_name) %>%
		filter(grepl("^[0-9]+$", chromosome_name) | grepl("^[X,Y]$", chromosome_name)) %>%
		
		# Group by 'external_gene_name' and choose the row with the highest version number
		group_by(external_gene_name) %>%
		filter(getversion(ensembl_gene_id_version) == max(getversion(ensembl_gene_id_version)) ) %>%
		
		# If there are any remaining duplicates, choose the one at the top of the group
		group_by(external_gene_name) %>%
		filter(row_number()==1)
	)
	save(result, file="biomaRt_query.Rdata")
}
print(paste("Total Annotation Rows: ", nrow(result), "Unique GeneID: ", length(unique(result$external_gene_name)), sep=''))

# Define Replicate Groups
sampleGroups<-list(
 	WT_0_Hour=c("S1", "S2", "S3"),
 	WT_48_Hour=c("S7", "S8", "S9"),
 	FN_0_Hour=c("S10", "S11", "S12"),
 	FN_48_Hour=c("S13", "S14", "S15")
)
 
pairwiseContrasts<-list(
	WT_0_Hour_vs_WT_48_Hour = c('WT_0_Hour', 'WT_48_Hour'),
	FN_0_Hour_vs_FN_48_Hour = c('FN_0_Hour', 'FN_48_Hour'),
	WT_0_Hour_vs_FN_0_Hour = c('WT_0_Hour', 'FN_0_Hour'),
	WT_48_Hour_vs_FN_48_Hour = c('WT_48_Hour', 'FN_48_Hour')
)

# Function Calculates fold change from logFC
calcFC<-function(lfc){
	if(!is.na(lfc) & lfc > 0){
		return(2^lfc)
	}
	else if(!is.na(lfc) & lfc < 0){
	     return(0 - 2^abs(lfc))
	}
	else { return(0)}
}

# Calculate RPKM Values based on counts in "GeneCount.tsv"
calcRPKM<-function(df, rSuffix="_RPKM", cSuffix="_GeneCount", lenCol="CodingLength"){
        countCols<-names(df)[grep(cSuffix, names(df))]
        for(i in countCols){
                rpkmCol<-gsub(cSuffix, rSuffix, i)
                print(paste(i, rpkmCol, sep=" "))
                s<-sum(df[,i])/1000000
                print(s)
                df[rpkmCol]<-(df[,i]/s)/(df[lenCol]/1000)
        }
        df
}

df<-calcRPKM(df)
write.table(df, "GeneCount.ReCalcRPKM.tsv", row.names=F, quote=F, sep="\t")

# Calculate Differential Expression based on 
for( g in pairwiseContrasts){
 	
  	# Get a list of Sample ID's in Group Order, corresponding to column names in "GeneCounts.tsv"
 	samples<-c(
		names(df)[sub("_GeneCount", "", names(df)) %in% sampleGroups[[g[1]]]],
		names(df)[sub("_GeneCount", "", names(df)) %in% sampleGroups[[g[2]]]]
	)
	
	# Define a factor for group assignment in the DGEList
	group<-factor(
		c(
			rep(1, length(sampleGroups[[g[1]]])), 
			rep(2, length(sampleGroups[[g[2]]]))
		), 
		labels=c(g[1], g[2])
	)
	
	# Build DGEList
	y<-DGEList(
		counts=df[,samples], 
		group=group, 
		genes=df[,c("PrimaryKey", "GeneID", "Chr", "CodingLength")]
	)
	row.names(y)<-y$genes$PrimaryKey

	# Filter DGEList: Remove genes where fewer than two samples have a cpm > 1
	keep <- rowSums(cpm(y) > 1) >= 2
	y.filter <-y[keep, ]
	y.filter$samples$lib.size <- colSums(y.filter$counts) # Fix library size after filtering
	
	# Calculate Normalization factors and dispersion estimates
	y.filter <-calcNormFactors(y.filter)
	y.filter <-estimateCommonDisp(y.filter, verbose=T)
	y.filter <-estimateTagwiseDisp(y.filter)
	
	# Calculate Differential Expression
	et<-exactTest(y.filter, pair=levels(group))
	degSet<-topTags(et, n=25000)@.Data[[1]]

	# Add Absolute Fold Change column "FC"
	degSet$FC<-sapply(degSet$logFC, calcFC)
	print(paste("Rows in degSet:", nrow(degSet)))
	
	# Join biomart query results to table of differentially expressed genes
	# on the gene symbol attribute. 
	degSet<-merge(degSet, result, by.x='GeneID', by.y='external_gene_name', all.x=T)
	file<-paste(outDir,"/",
		paste(g, collapse="_vs_"),
		"_expressedTags-all.txt", sep=""
	)
	print(paste("Rows in degSet after biomaRt:", nrow(degSet)))
	
	# Reorder Columns
	degSet<-degSet %>% dplyr::select(PrimaryKey,
				GeneID, ensembl_gene_id, description,
				FC, logFC, logCPM, PValue, FDR, CodingLength, Chr
	)
	
	# Join counts per million (calculated prior to filtering genes)
	cpmy<-as.data.frame(cpm(y))
	names(cpmy)<-gsub("Count","Count.cpm", names(cpmy))
	cpmy$PrimaryKey<-row.names(y)
	degSet<-merge(x=degSet, y=cpmy, by='PrimaryKey')
	print(paste("Rows after join cpm:", nrow(degSet)))

	# Join read counts
	counts<-as.data.frame(y$counts)
	names(counts)<-gsub("Count", "Count.counts", names(counts))
	counts$PrimaryKey<-row.names(counts)
	degSet<-merge(x=degSet, y=counts, by='PrimaryKey')
	print(paste("Rows after join counts:", nrow(degSet)))

	# Join rpkm
	degSet<-merge(x=degSet, y=df[,c("PrimaryKey", gsub("_GeneCount", "_RPKM", samples))], by="PrimaryKey")
	print(paste("Rows after join RPKM:", nrow(degSet)))
	
	write.table(degSet,
		file=file,
		sep="\t", row.names=F, quote=F
	)
}


print(sessionInfo())