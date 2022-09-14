wgcna <- getwd()
library(WGCNA, warn.conflicts=FALSE)
library(reshape2, warn.conflicts=FALSE)
main <- file.path(wgcna, '..', '..')
args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]

data_path <- paste(main, '/temp/', prefix, '_expression_data.csv', sep= "")
out_path <- paste(main, '/temp/', prefix,  '_edge_list.csv', sep = "")

geneData <- read.csv(data_path, header = TRUE, sep='\t', as.is=TRUE)
geneData <- geneData[-c(1)]

Corr <- adjacency(geneData, type = "signed", power = 1)
Adjacency <- signumAdjacencyFunction(Corr, threshold = 0.6)

#Obtaining the edge list from the provided matrix
gen.mat.to.edge.list<-function(mat,symmetric=TRUE,diagonal=FALSE,text=FALSE){
	
	mat<-as.matrix(mat)
	if(symmetric){mat[lower.tri(mat)]<-NA}
	if(!diagonal){diag(mat)<-NA}
	obj<-melt(mat)
	colnames(obj)<-c("source","target","value")
	obj<-obj[!is.na(obj$value),]
	obj<-obj[(obj$value != 0),]
	if(!text){obj$value<-as.numeric(as.character(obj$value))}
	return(obj)
}

Edges <- gen.mat.to.edge.list(Adjacency)
write.table(Edges, out_path, sep='\t')

