#################################################################################
# File: BuildCountMatrix.R                                                      #
# Purpose: Implement functions to traverse a set of data directories, find      #
#          count data files for individual samples and join them into a data    #
#          matrix                                                               #
# Created: April 30, 2019                                                       #
# Changes:                                                                      #
#         18Jun2019                                                             #
#          - Add functionality for processing StringTie results                 #
#          - update function names to distinguish htseq methods from stringtie  #
#          - functions prefixed with "hc_" process htseq-count data             #
#          - functions prefixed with "st_" process stringtie data               #
#                                                                               #
# Author: Adam Faranda                                                          #
#################################################################################

# Iterate over a list of data directories containing htseq-count output files, 
# return a two column table of datafiles in each directory that match a target pattern
# dl<-c("DBI_NoTrim_HTSeq_Count", "DNA_Link_NoTrim_HTSeq_Count")

dl<-c(
  "/work/abf/LEC_Time_Series/DBI_NoTrim_HTSeq_Count_Gene", 
  "/work/abf/LEC_Time_Series/DNA_Link_NoTrim_HTSeq_Count_Gene"
)

hc_getFileTable <-function(
	dirList, pattern='GeneCount.txt', 
	filename="DataFileAnnotations.csv",
){
	if(filename %in% list.files()){
		ft<-read.csv(filename)
	} 
	else{
		ft<-data.frame(
			directory=character(),
			filename=character(),
			stringsAsFactors=F
		)
		
		for (d in dirList){
			 list.files(d, pattern)
			 x<-data.frame(
			 	directory = d,
			 	filename = list.files(d, pattern),
			 	stringsAsFactors=F
			 )
			 ft<-rbind(ft, x)
		}
		write.csv(ft, "DataFileAnnotations.csv", row.names=F)
	}
	return(ft)
}

# Given a file table Load datafiles into a list object
hc_loadFiles<-function(ft){
	dataSets<-list()
	for(i in 1:nrow(ft)){
		file<-paste(ft[i,1], ft[i,2], sep='/')
		print(file)
		dataSets[[file]]<-read.table(
			file, header=F, stringsAsFactors=F, sep='\t'
		)
	}
	dataSets
}

# Check Identifier Consistency by iterating over each datafile
hc_identifierConsistency<-function(ds, ft, idCol=1){
	ft$rowCount<-0
	ft$uniqueID<-0
	check<-ds[[1]][,idCol]
	for (i in 1:nrow(ft)){
		sample<-paste(ft[i,1], ft[i,2], sep='/')
		ft$rowCount<-nrow(ds[[sample]])
		ft$uniqueID<-length(
						unique(
							ds[[sample]][,idCol]
						)
		)
		if (sum(check != ds[[sample]][,idCol]) > 0){
			print(paste("ID mismatch between first sample and:", sample))
		}
	}
	ft
}

# Join count data columns into a matrix
hc_buildDataFrame<-function(ds, ft, idCol=1){
	ft$Sample_Number <-''
	df<-data.frame(Ensembl=ds[[1]][,idCol], stringsAsFactors=F)
	for( i in 1:nrow(ft)){
		sample<-paste(ft[i,1], ft[i,2], sep='/')
		if(!any(df$Ensembl != ds[[sample]][,idCol])){
			df<-merge(df, ds[[sample]], 
			by.x='Ensembl', by.y=1,
			sort = F)
			names(df)[grep('V2', names(df))]<-paste("Sample_", i, sep='')
			ft$Sample_Number[i]<-paste("Sample_", i, sep='')
		}
		else{
			print(sum(df$Ensembl != ds[[sample]][,idCol]))
			print(paste("Can't join sample:",sample, "ID Mismatch" ))
		}
	
	}
	list(ft, df[order(df[,idCol]),])
}

# Annotate File Table with Group Information -- manually add annotations
ft<-hc_getFileTable(dl)
ds<-hc_loadHtSeqFiles(ft)
ft<-hc_identifierConsistencyHtSeq(ds, ft)
raw<-hc_buildHtSeqDataFrame(ds, ft)








