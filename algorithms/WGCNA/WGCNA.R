
wgcna <- getwd()
#if (!require("WGCNA", quietly = TRUE))
  #  install.packages("WGCNA", repos = "http://cran.us.r-project.org")
#if (!require("igraph", quietly = TRUE))
  # install.packages("igraph", repos = "http://cran.us.r-project.org")
print(wgcna)
#dir1 <- paste0(wgcna, "/WGCNA")
#dir2 <- paste0(wgcna, "/igraph")

#dir.create(dir1)
#dir.create(dir2)

#install.packages("WGCNA", lib = dir1, dependencies = TRUE, repos = "http://cran.us.r-project.org")
#install.packages("igraph", lib = dir2, dependencies = TRUE, repos = "http://cran.us.r-project.org")

#library.loc(WGCNA, loc = dir1)
#library.loc(igraph, loc = dir2)
#library(WGCNA, lib.loc = dir1)
#library(igraph, lib.loc = dir2)
library(WGCNA)
library(reshape2)
main <- file.path(wgcna, '..', '..')
args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]
#print(prefix)
data_path <- paste(main, '/temp/', prefix, '_expression_data.csv', sep= "")
out_path <- paste(main, '/temp/', prefix,  '_edge_list.csv', sep = "")
#adjacency_path <- paste(main, '/temp/', prefix, '_adjacency.csv', sep = "")

geneData <- read.csv(data_path, header = TRUE, sep='\t', as.is=TRUE)
geneData <- geneData[-c(1)]
#print(dim(geneData))
#print(geneData[1:10,1:10])

#Corr <- adjacency(geneData[c(1:dim(geneData)[1]),c(1:20)], type = "signed", power = 1)
Corr <- adjacency(geneData, type = "signed", power = 1)
#print(Corr[c(1:10),c(1:20)])

#Thresold <-pickHardThreshold.fromSimilarity(Corr)
Adjacency <- signumAdjacencyFunction(Corr, threshold = 0.6)

#Obtaining the edge list from the provided matrix
gen.mat.to.edge.list<-function(mat,symmetric=TRUE,diagonal=FALSE,text=FALSE){
	
	mat<-as.matrix(mat)
	if(symmetric){mat[lower.tri(mat)]<-NA}
	if(!diagonal){diag(mat)<-NA}
	obj<-melt(mat)
	colnames(obj)<-c("source","target","value")
	obj<-obj[!is.na(obj$value),]
	#obj$value[obj$value=="nna"]<-NA
	obj<-obj[(obj$value != 0),]
	if(!text){obj$value<-as.numeric(as.character(obj$value))}
	return(obj)
	# remove duplicates
}

#write.csv(Adjacency, adjacency_path, sep='\t')
Edges <- gen.mat.to.edge.list(Adjacency)
#(Edges)
#print(dim(Edges))

#print(dim(Adjacency))
#write.table(Adjacency, adjacency_path, sep='\t')
#Edges<<-graph_from_adjacency_matrix(Adjacency, mode = c("undirected"))
#Edges<<- as.data.frame(get.edgelist(Edges))
#print(dim(Edges))
#names(Edges)[1] <- 'node_lower'
#names(Edges)[2] <- 'node_upper'
#Edges$edges <<- paste(Edges[,c(1)], Edges[,c(2)], sep = "-")
write.table(Edges, out_path, sep='\t')

