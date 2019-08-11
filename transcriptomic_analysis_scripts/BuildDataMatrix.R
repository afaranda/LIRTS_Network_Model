#################################################################################
# File: BuildCountMatrix.R                                                      #
# Purpose: Implement functions to traverse a set of data directories, find      #
#          count data files for individual samples and join them into a data    #
#          matrix                                                               #
# Created: April 30, 2019                                                       #
# Changes:                                                                      #
#         June 2019                                                             #
#          - Add functionality for processing StringTie results                 #
#          - update function names to distinguish htseq methods from stringtie  #
#          - functions prefixed with "hc_" process htseq-count data             #
#          - functions prefixed with "st_" process stringtie data               #
#                                                                               #
# Author: Adam Faranda                                                          #
#################################################################################
library(dplyr)
library(ballgown)
library(rtracklayer)

#################################################################################
#   Function to tabulate gene and transcript lengths (exon union) from gtf      #
#################################################################################
# Given the path to a GTF File return a table of gene / transcript exon-union 
# lengths.  This takes a very long time to run. 
lengthTable<-function(gtfpath){
  orig_gtf<-readGFF(gtfpath)
  
  xs_len<-orig_gtf %>%
    mutate( length = abs(start - end) + 1) %>% 
    filter( type == "exon") %>% 
    group_by(transcript_id) %>%
    summarise(gene_id = unique(gene_id), Xscript_Exons = n(), Xscript_Length=sum(length))
  
  gn_len<-orig_gtf %>%
    mutate( length = abs(start - end) + 1) %>% 
    filter( type == "exon" & !duplicated(exon_id)) %>% 
    group_by(gene_id) %>%
    summarise( Union_Exons = n(), Union_Length=sum(length))
  
  as.data.frame(
    lenTable<-inner_join(
      xs_len, 
      gn_len,
      by ="gene_id"
    ),
    stringsAsFactors=F
  )
}

#################################################################################
#                 Functions to build HTSeq-Count Matrices                       #
#################################################################################

#Iterate over a list of data directories containing htseq-count output files, 
# return a two column table of datafiles in each directory that match a target 
# pattern

hc_getFileTable <-function(
  wd = ".",
	dirList, pattern='_GeneCount.txt', 
	filename="HtSeq_GeneCountFiles.csv"
){
	if(filename %in% list.files(wd)){
		ft<-read.csv(filename, stringsAsFactors = F)
	} 
	else{
		ft<-data.frame(
		  sample = character(),
			directory=character(),
			filename=character(),
			stringsAsFactors=F
		)
		
		for (d in dirList){
			 list.files(d, pattern)
			 x<-data.frame(
			  sample =  gsub(pattern, "", list.files(d, pattern)),
			 	directory = d,
			 	filename = list.files(d, pattern),
			 	stringsAsFactors=F
			 )
			 ft<-rbind(ft, x)
		}
		#write.csv(ft, filename, row.names=F)
	}
	return(ft)
}

# Given a file table Load datafiles into a list object
hc_loadFiles<-function(ft){
	dataSets<-list()
	for(i in 1:nrow(ft)){
		file<-paste(ft[i,2], ft[i,3], sep='/')
		print(ft[i,1])
		sample<-ft[i,1]
		print(file)
		dataSets[[sample]]<-read.table(
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
		ft$rowCount<-nrow(ds[[ft[i,1]]])
		ft$uniqueID<-length(
						unique(
							ds[[ft[i,1]]][,idCol]
						)
		)
		if (sum(check != ds[[ft[i,1]]][,idCol]) > 0){
			print(paste("ID mismatch between first sample and:", sample))
		}
	}
	ft
}

# Join count data columns into a matrix
hc_buildDataFrame<-function(ds, ft, idCol=1, measCol=2){
	df<-data.frame(Ensembl=ds[[1]][,idCol], stringsAsFactors=F)
	for( i in 1:nrow(ft)){
		if(!any(df$Ensembl != ds[[ft[i,1]]][,idCol])){
			df<-merge(df, ds[[ft[i,1]]][,c(idCol, measCol)], 
			by.x='Ensembl', by.y=1,
			sort = F)
			names(df)[grep('V2', names(df))]<-ft[i,1]
		}
		else{
			print(sum(df$Ensembl != ds[[ft[i,1]]][,idCol]))
			print(paste("Can't join sample:",sample, "ID Mismatch" ))
		}
	
	}
	df<- df %>% filter(!grepl("__",Ensembl))
	row.names(df)<-df$Ensembl
	df<-df[order(df[,"Ensembl"]),]
	dm<-as.matrix(df[,!grepl("Ensembl", names(df))])
	list(ft, dm)
}

# Compile expression data into a DGE list
buildDGE<-function(mat, ft=ft, ft_id="sample", gt=gt, gt_id="gene_id"){
  # Verify that all ft_id and gn_id are unique in their respective tables
  # and that both match the count matrix
  if(length(unique(ft[,ft_id])) == nrow(ft)){
    if(length(unique(gt[,gt_id])) == nrow(gt)){
      if(!any(colnames(mat) != ft[, ft_id])){
        if(!any(rownames(mat) != gt[,gt_id])){
          row.names(ft)<-ft[,ft_id]
          row.names(gt)<-gt[,gt_id]
          dgelist<-DGEList(
            counts = mat,
            samples = ft,
            genes = gt
          )
        } else if(length(setdiff(rownames(mat), genes[,gt_id])) == 0){
          print("Sorting rows by Gene ID ")
          row.names(ft)<-ft[,ft_id]
          row.names(genes)<-gt[,gt_id]
          gt<-genes[order(gt[,gt_id])]
          mat<-mat[order(rownames(mat)),]
          dgelist<-DGEList(
            counts = mat,
            samples = ft,
            genes = genes
          )
        } else {
          print("Genes given did not match data matrix")
          return(NULL)
        }
      } else if(length(setdiff(colnames(mat), ft[,ft_id])) == 0) {
        print("Sorting Columns by Sample ID and Genes by Gene ID")
        row.names(ft)<-ft[,id_col]
        row.names(genes)<-genes[,id_col]
        
        ft<-ft[order(ft[,ft_id]),]
        genes<-genes[order(genes[,id_col])]
        
        mat<-mat[,order(colnames(mat))]
        mat<-mat[order(rownames(mat)),]
        dgelist<-DGEList(
          counts = mat,
          samples = ft,
          genes = genes
        )
      } else {
        print("Columns in Sample data do not match columns in count matrix")
        return(NULL)
      }
    }
  }
}

# Drop HTSeq-Samples
hc_dropSamples<-function(mat, ft=ft, id_col=1,samples=c()){
  ft<-ft[!(ft[,id_col] %in% samples),]
  s<-ft[,id_col] 
  mat<-mat[,s]
  if(!any(colnames(mat)!= ft$sample)){
    return(list(ft, mat))
  } else {
    print("Column Matching Error")
  }
}
#################################################################################
#           Functions to build string-tie based matrices (TPM)                  #
#################################################################################

# Helper Function -- returns the last element in a delimited list
splitLast<-function(x, delim="/"){
  s<-strsplit(x, delim)
  s<-sapply(s, function(y) y[length(y)])
  s
}

st_getFileTable <-function(
  dirList,
  wd = ".",
  filename="StringTie_TPM_Files.csv"
){
  if(filename %in% list.files(wd)){
    ft<-read.csv(filename, stringsAsFactors=F)
  } 
  else{
    ft<-data.frame(
      sample=character(),
      directory=character(),
      filename=character(),
      gtffile=character(),
      stringsAsFactors=F
    )
    
    for (d in dirList){
      dr<-paste(d, list.files(d), sep = "/")
      
      x<-data.frame(
        sample = as.character(list.files(d)),
        directory = dr, 
        filename = character(length(dr)),
        gtffile =character(length(dr)),
        stringsAsFactors=F
      )
      for (r in 1:nrow(x)){
        d<-x[r, 'directory']
        fn<-paste(splitLast(d), ".txt", sep="")
        gt<-paste(splitLast(d), ".gtf", sep="")
        if(file.exists(paste(d, fn, sep="/"))){
          x[r, 'filename']<-fn
          x[r, 'gtffile'] <-gt
        }
      }
      ft<-bind_rows(ft, x)
    }
    ft<-ft%>% filter(filename != "")
    #write.csv(ft, filename, row.names=F)
  }

  return(ft)
}

# Returns a ballgown object based on data files listed in 
# the given file table (ft)
st_sewBallGown<-function(
  ft = ft,
  wd = ".",
  filename="bg.Rdata"
){
  if(filename %in% list.files(wd)){
    load(filename)
  }
  else{
    bg<-ballgown(ft$directory, pData=ft, meas="all")
    save(bg, file="bg.Rdata")
  }
  bg
}

# Imports GTF files for each sample into data frames, retruns 
# a list of data frames indexed by sample ID
st_loadFiles<-function(ft, fnCol='gtffile'){
  dataSets<-list()
  gtfNumCols<-c(
    "start", "end", "score",
    "cov", "FPKM", "TPM", "exon_number" 
  )
  for(i in 1:nrow(ft)){
    file<-paste(ft[i,2], ft[i,fnCol], sep='/')
    sample<-ft[i,1]
    print(file)
    if(grepl(".gtf",splitLast(file))){
      dataSets[[sample]]<-readGFF(file)
      dataSets[[sample]]$PK<-paste(
        dataSets[[sample]]$type,
        dataSets[[sample]]$start,
        dataSets[[sample]]$end,
        dataSets[[sample]]$transcript_id,
        sep=""
      )
      for(c in gtfNumCols){
        dataSets[[sample]][,c]<-as.numeric( dataSets[[sample]][,c])
      }
    }
    else if(grepl(".txt",splitLast(file))){
      dataSets[[sample]]<-read.table(file, header = T, stringsAsFactors = F, sep="\t")
    }
  }
  dataSets
}



# Calculate the sum of a measurement (FPKM, TPM) over genes in a data set
# returns a three column table: 
#           Gene ID, number of transcripts, summed measurement
sum_tr_over_gn<-function(df, g_idCol, measCol){
  df<-df[,c(g_idCol, measCol)]
  names(df)<-c("g_id", "meas")
  
  return(
    as.data.frame(
     df %>% 
        group_by(g_id) %>%
        summarise( Transcripts = n(), AvgMeasure = sum(as.numeric(meas)))
    )
  )
}

# Join Transcript level measurements (FPKM or TPM) into a data matrix
st_buildTranscriptMatrix<-function(ds, idCol=17, measCol=14){
  df<-data.frame(
    Ensembl=(ds[[1]] %>% filter(type == "transcript"))[,idCol],
    stringsAsFactors=F
  )
  for( i in 1:nrow(ft)){
    dg<-ds[[ft[i,1]]] %>% filter(type == 'transcript')
    if(
      length(intersect(df$Ensembl, dg[,idCol])) == length(unique(df$Ensembl))
    ){
      dg<-ds[[ft[i,1]]] %>% filter(type == 'transcript')
      df<-merge(df, dg[,c(idCol, measCol)], 
                by.x='Ensembl', by.y=1,
                sort = F)
      mcn<-names(ds[[ft[i,1]]][measCol])
      names(df)[grep(mcn, names(df))]<-ft[i,1]
      print(i)
    }
    else{
      print(paste("Can't join sample:",ft[i,1], "ID Mismatch" ))
    }
    
  }
  df<- df %>% filter(!grepl("__",Ensembl))
  row.names(df)<-df$Ensembl
  df<-df[order(df[,"Ensembl"]),]
  dm<-as.matrix(df[,!grepl("Ensembl", names(df))])
  list(ft, dm)
}

# Join Gene level measurements (FPKM or TPM) into a data matrix
# transcript level measurements are summed over their respective genes
st_buildGeneMatrix<-function(ds, g_idCol=9, measCol=14){
  df<-data.frame(
    Ensembl=unique(ds[[1]][,g_idCol]),
    stringsAsFactors=F
  )
  for( i in 1:nrow(ft)){
    dg<-ds[[ft[i,1]]] %>% filter(type == 'transcript')
    dg<-sum_tr_over_gn(dg, g_idCol=g_idCol, measCol=measCol)
    if(
      length(intersect(df$Ensembl, dg[,1])) == length(unique(df$Ensembl))
    ){
      
      df<-merge(df, dg[,c(1, 3)], 
                by.x='Ensembl', by.y=1,
                sort = F)
      mcn<-names(dg[3])
      names(df)[grep(mcn, names(df))]<-ft[i,1]
      print(i)
    }
    else{
      print(paste("Can't join sample:",ft[i,1], "ID Mismatch" ))
    }
    
  }
  row.names(df)<-df$Ensembl
  df<-df[order(df[,"Ensembl"]),]
  dm<-as.matrix(df[,!grepl("Ensembl", names(df))])
  list(ft, dm)
}


# Validate lengthTable Counts and lengths internally
#any(lt$Xscript_Exons > lt$Union_Exons)
#any(lt$Xscript_Length > lt$Union_Length)

# Iterate over genes and validate The number of unique exons in each 
# check<-data.frame(
#   gene_id=character(),
#   Ex = numeric(),
#   stringsAsFactors = F
# )
# for (gn in 5000:5500){
#   print(gn)
#   g<-gn_len$gene_id[gn]
#   check[gn,1]<-g
#   check[gn,2]<-length(
#     unique(
#       orig_gtf[orig_gtf$gene_id == g & orig_gtf$type == "exon", "exon_id"]
#       )
#     )
# }
# inner_join(
#   gn_len,
#   check,
#   by="gene_id"
# ) %>%   
#   mutate(ck=abs(Union_Exons - Ex)) %>%
#   summarize(sum(ck))


