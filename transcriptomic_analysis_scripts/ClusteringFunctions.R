################################################################################
# File: ClusteringFunctions.R												   #
# Purpose: Wrappers for various clustering algorithms, and plotting functions. #
# Created: May 1, 2019														   #
# Author: Adam Faranda														   #
################################################################################
library(dplyr)
library(ggfortify)

wrapHclust<-function(df, idCol=1,
	transpose=T, d.meth="euclidean", d.p = 2,
	h.method = "average"
 ){
 	if(!is.null(idCol) & idCol !=0){
		idCol<-ifelse(is.character(idCol), colNum(x, idCol), idCol)
		row.names(df)<-df[,idCol]
		df<-df[,setdiff(1:ncol(df), idCol)]
	}
	if(transpose){df<-t(df)}
	d<-dist(df, method=d.meth, p=d.p)
	h<-hclust(d, method=h.method)
}

# Plot dendrograms for heirarchical clusters -- only use when clustering samples
plotHclust<-function(h, ft, sampleCol=7, labelCol=8, colorCol=4, main=''){
	labColor <-ifelse(ft[h$order,colorCol] == 'DNA', "red", "black")
	x<-1:nrow(ft)
	lab <- ft[h$order, labelCol]
	print(lab)
	plot(
		h, labels = F,
		main = main,
		xlab = '',
		sub = '', 
		hang = -1
	)
	#print(labColor)
	text(x=x, y=-1, labels = lab, col=labColor, srt = 90, xpd=NA, adj=c(1.2,0.5))
}


tabulate_H_Clusters<-function(h, ks=c(1:5) ){
	hc_table <- data.frame(
		Sample = h$labels
	)
	for (k in ks){
		print(k)
		if(k <= nrow(hc_table)){
			c=data.frame(c=cutree(h, k))
			c$Sample = row.names(c)
			hc_table<-merge(
				hc_table, c,
				by = 'Sample'
			)
			name<-paste("k_eq_", k, sep='')
			names(hc_table)[grep('c', names(hc_table))]<-name
		}
	}
	row.names(hc_table)<-hc_table$Sample
	hc_table<-hc_table[,setdiff(1:ncol(hc_table), grep('Sample', names(hc_table)))]
	hc_table
}


tabulate_k_means<-function(df, idCol=1, ks=c(1:5), transpose=T, method='Hartigan-Wong'){
	if(!is.null(idCol) & idCol != 0){
		idCol<-ifelse(is.character(idCol), colNum(x, idCol), idCol)
		row.names(df)<-df[,idCol]
		df<-df[,setdiff(1:ncol(df), idCol)]
	}
	if(transpose){df<-t(df)}
	hc_table <- data.frame(
		Sample = row.names(df)
	)
	for (k in ks){
		if(k <= nrow(hc_table)){
			c=data.frame(
				c=kmeans(df, nstart=nrow(hc_table), centers=k, 
				iter.max=30, algorithm=method )[[1]])
			c$Sample = row.names(c)
			hc_table<-merge(
				hc_table, c,
				by = 'Sample'
			)
			name<-paste("k_eq_", k, sep='')
			names(hc_table)[grep('c', names(hc_table))]<-name
		}
	}
	row.names(hc_table)<-hc_table$Sample
	hc_table<-hc_table[,setdiff(1:ncol(hc_table), grep('Sample', names(hc_table)))]
	hc_table
}

# Helper function to allow choice of statistics
apfunc<-function(x, applyFun='mean'){
	if(applyFun == 'mean'){
		return(mean(x, na.rm=T))
	} else if(applyFun == 'sd'){
		return(sd(x, na.rm=T))
	} else if(applyFun == 'max'){
		return(max(x, na.rm=T))
	} else if(applyFun == 'min'){
		return(min(x, na.rm=T))
	} else if(applyFun == 'length'){
		return(length(x))
	}
}

getClusterStat<-function(mat, ktable, k, applyFun='mean', transpose=F){
	if(transpose){mat<-t(mat)}
	mat<-as.data.frame(mat)
	ktable<-as.data.frame(ktable[k])
	mat$ID<-row.names(mat)
	kname<-as.character(names(ktable))
	ktable$ID<-row.names(ktable)
	
	x<-merge(mat, ktable, by='ID')
	row.names(x)<-x$ID
	x<-x[2:ncol(x)]
	
	for (i in unique(x[,kname])){
		y<-x[x[,kname] == i,setdiff(names(x), kname)]
		m<-as.list(apply(y, 2, apfunc, applyFun=applyFun))
		m[['count']]<-nrow(y)
		m[['cluster']]<-i
		#nc<-grep('kname', names(x))
		#x<-x[, c(nc, setdiff(1:ncol(x), nc))]
		if(exists('clust_table')){clust_table<-rbind(clust_table, m)}
		else{clust_table<-as.data.frame(m)}
	}
	clust_table<-as.data.frame(clust_table)
	row.names(clust_table)<-clust_table$cluster
	clust_table
}

reshapeClusterTable<-function(mat, ktable, k, ft, transpose=F){
	if(transpose){mat<-t(mat)}
	mat<-as.data.frame(mat)
	ktable<-as.data.frame(ktable[k])
	mat$ID<-row.names(mat)
	kname<-as.character(names(ktable))
	ktable$ID<-row.names(ktable)
	
	# Join Cluster table on Gene Expression Matrix
	x<-merge(mat, ktable, by='ID')
	names(x)[grep(paste('k_eq_',k,sep=''), names(x))]<-'Cluster'
	
	# Refactor gene expression matrix into 3NF compliant table
	x<-melt(x, id.vars=c('ID', 'Cluster'))
	names(x)[grep('variable', names(x))]<-'Sample'

	# Join Feature Table values on Sample
	x$Class <-droplevels(sapply(x$Sample, function(s) ft[ft$Sample_Number == s, 'Class']))
	x$Hours_PCS <-as.numeric(sapply(x$Sample, function(s) ft[ft$Sample_Number == s, 'Hours_PCS']))
	x$Lab <-droplevels(sapply(x$Sample, function(s) ft[ft$Sample_Number == s, 'Seq_Lab']))
	x
}

randIndex<-function(df, instanceCol=1, clusterCol=2, classCol=3){
	x.lev<-levels(as.factor(df[,classCol]))
	y.lev<-levels(as.factor(df[,clusterCol]))
	
	instance<-df[,instanceCol]
	pairs<-t(combn(instance, 2))
	r<-data.frame(
		Sample1 = pairs[,1],
		Sample2 = pairs[,2],
		matchClass = 0,
		matchCluster = 0
	)
	for(i in 1:nrow(r)){
		X1 <- df[df[,instanceCol]==r[i,'Sample1'],classCol]
		X2 <- df[df[,instanceCol]==r[i,'Sample2'],classCol]
		
		Y1 <- df[df[,instanceCol]==r[i,'Sample1'],clusterCol]
		Y2 <- df[df[,instanceCol]==r[i,'Sample2'],clusterCol]
		if(X1 == X2){
			r[i, 'matchClass']<-1
		}
		if(Y1 == Y2){
			r[i, 'matchCluster']<-1
		}
	}
	sum(r$matchClass == r$matchCluster)/choose(length(instance), 2)
}

summarizeSampleClusters<-function(data=ecpm, distm, linkm, v=50, label='ecpm'){
data<-varianceFilter(data, threshold=v)
row.names(data)<-data$ID
	for(d in distm){
		for(l in linkm){
			f1<-paste(label,'_Samples_Top_', v,'_',d,'_',l,'_cluster.png')
			mat<-data[,2:19]
			h1<-wrapHclust(mat, 
				idCol=0, transpose=T, d.meth=d, h.method = l
			)
			kt1<-tabulate_H_Clusters(h1, ks = trees)
			kt1$I<-row.names(kt1)
			kt1$TC<-sapply(kt1$I, function(i) ft[ft$Sample_Number ==i, 'Class'])
			kt1$TC<-sapply(kt1$I, function(i) ft[ft$Sample_Number ==i, 'Seq_Lab'])
			png(f1, width=240, height=200)
				plotHclust(h1, ft, sampleCol=7, labelCol=8, colorCol=4, main='')
			dev.off()
			for(t in trees){
				r<-randIndex(kt1, length(trees)+1, t, length(trees)+2)
				if(t > 1){
					sil<-as.data.frame(
						silhouette(
							cutree(h1, t), 
							dist(t(mat), method= d)
						)[,1:3]
					)
					sm<-mean(sil$sil_width)
					sil<-sil %>% 
						group_by(cluster) %>%
						summarize( mean(sil_width), n())
					slist<-paste(sil$`n()`, '(',round(sil$`mean(sil_width)`,3),')',sep='')
					print(slist)
				}
				else{ 
					sm<-0
					slist<-"None"
				}
				if(!exists('randTable')){
					randTable<-data.frame(
						Data = label,
						DistMethod = d,
						LinkMethod = l,
						NumClusters= t,
						RandIndex = r,
						MeanSilhouette=sm, 
						ClusterSilhouettes=paste(slist, collapse=', '),
						stringsAsFactors=F
					)
				}
				else{
					randTable<-rbind(
						randTable, 
						data.frame(
							Data = label,
							DistMethod = d,
							LinkMethod = l,
							NumClusters= t,
							RandIndex = r,
							MeanSilhouette=sm, 
							ClusterSilhouettes=paste(slist, collapse=', '),
							stringsAsFactors=F
						)
					)
				}
			}
		
		}
	}
 	randTable
}

summarizeGeneClusters<-function(
	m=mat, d.meth='manhattan', 
	h.method='complete', label='foo', 
	transpose=F, nclust=20, minGenes = 5,
	clustSD = 2000
){
	print('one')
	h<-wrapHclust(m, idCol=0, transpose=transpose, d.meth=d.meth, h.method=h.method)
	ktable<-tabulate_H_Clusters(h, ks=1:nclust)
	print('two')
	ct<-reshapeClusterTable(m, ktable, ft, k=nclust)
	print('there')
	s<-as.data.frame(
		silhouette(cutree(h, nclust), dist(m, method=d.meth))[,1:3]
	)
	s.means<-s %>% 
		group_by(cluster) %>%
		summarize(Mean_Silhouette=mean(sil_width))
	
	cstat<-inner_join(
		ct %>% 
			group_by(Cluster) %>%
			summarize(Count = n()/18, Mean_logCPM=mean(value), STDev_logCPM = sd(value)),		
		s %>%
			group_by(cluster) %>%
			dplyr::rename(Cluster = cluster) %>%
			summarize(Mean_Silhouette=mean(sil_width)),
			by='Cluster'
	)
	write.csv(cstat, paste(label, '_Gene_Cluster_Stats.csv', sep=''))
				
	cl<- cstat %>%
		filter(Count > minGenes, STDev_logCPM <clustSD)
	cl<-cl$Cluster
	
	bp<-ggplot(data=ct[ct$Cluster %in%cl,], mapping=aes(x=Class, y=value, color=Class)) + 
		geom_boxplot() + facet_grid( . ~Cluster)
	
	png(paste(label,'_Gene_Cluster_Profiles.png',sep=''), width=960, height=280)
		print(bp)
	dev.off()
	
	png(paste(label,'_Heatmap.png', sep=''))
		pheatmap(
			m, clustering_distance_cols=d.meth, 
			clustering_distance_rows=d.meth,
			clustering_method=h.method,
			scale='none'
		)
	dev.off()
		return(list(cstat, bp))
}






