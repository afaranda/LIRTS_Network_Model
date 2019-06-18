#################################################################################
# File: BuildCountMatrix.R                                                      #
# Purpose: Implement functions to traverse a set of data directories, find      #
#          count data files for individual samples and join them into a data    #
#          matrix                                                               #
# Created: April 30, 2019                                                       #
# Author: Adam Faranda                                                          #
#################################################################################

# Iterate over a list of data directories, return a two column table of 
# datafiles in each directory that match a target pattern
dl<-c("DBI_Trimmed_HTSeq_Count", "DNA_Link_HTSeq_Count")

getFileTable <-function(
	dirList, pattern='GeneCount.txt', 
	filename="DataFileAnnotations.csv"
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
loadFiles<-function(ft){
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
identifierConsistency<-function(ds, ft, idCol=1){
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
buildDataFrame<-function(ds, ft, idCol=1){
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
ft<-getFileTable(dl)
ds<-loadFiles(ft)
ft<-identifierConsistency(ds, ft)
raw<-buildDataFrame(ds, ft)








