################################################################################
# File: PrincipalComponents.R  												   #
# Purpose: Given a matrix of gene expression values, compute principal         #
#          components. Plot components and return corresponding data tables.   #
#          Also implement MDS methods from edgeR                               #
#                                                                              #
# Created: April 30, 2019 													   #
# Author: Adam Faranda														   #
################################################################################
library(ggfortify)

# Get the position(s) (column number) of (a) column(s) by name from a data frame
colNum<-function(df, n){
	positions<-c()
	for(i in n){
		positions<-c(
			positions,
			grep(i, names(df))
		)
	}
	positions
}

plotPrinComp<-function(df, ft, idCol=1, groupCol=3, transpose=T, legendTitle=""){
	idCol<-ifelse(is.character(idCol), colNum(df, idCol), idCol)
	groupCol<-ifelse(is.character(groupCol),
        groupCol,
        names(ft)[groupCol]
    )
    if(idCol !=0 ){
	    row.names(df)<-df[,idCol]
	    df<-df[, setdiff(1:ncol(df), idCol)]
	}
    if(transpose){
        df<-t(df)
    }
    pca<-prcomp(df, scale=T)
    p<-autoplot(pca, data = ft, colour=groupCol)
    p + labs(colour=legendTitle)
    p
    
}





?autoplot


