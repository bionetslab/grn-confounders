
wgcna <- getwd()
if (!require("WGCNA", quietly = TRUE))
    install.packages("WGCNA", repos = "http://cran.us.r-project.org")
if (!require("igraph", quietly = TRUE))
   install.packages("igraph", repos = "http://cran.us.r-project.org")

library(WGCNA)
library(igraph)
enableWGCNAThreads()

#main <- file.path(wgcna, '..', '..')
args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]
print(prefix)
data_path <- paste(prefix, '/temp/', 'wgcna_expression_data.csv', sep= "")
out_path <- paste(prefix, '/temp/', 'wgcna_edge_list.csv', sep = "")
#adjacency_path <- paste(prefix, '_adjacency.csv', sep = "")

geneData <- read.csv(data_path, header = TRUE, sep='\t', as.is=TRUE)
geneData <- geneData[-c(1)]
#print(dim(geneData))
#print(geneData[1:10,1:10])

#Corr <- adjacency(geneData[c(1:dim(geneData)[1]),c(1:20)], type = "signed", power = 1)
Corr <- adjacency(geneData, type = "signed", power = 1)
#Corr[c(1:8),c(1:8)]
#print(dim(Corr))
Adjacency <- signumAdjacencyFunction(Corr, threshold = 0.6)
#print(Adjacency[1:10,1:10])
#print(dim(Adjacency))
#write.csv(Adjacency, adjacency_path, sep='\t')
Edges<<-graph_from_adjacency_matrix(Adjacency, mode = c("undirected"))
Edges<<- as.data.frame(get.edgelist(Edges))
#print(dim(Edges))
names(Edges)[1] <- 'node_lower'
names(Edges)[2] <- 'node_upper'
#Edges$edges <<- paste(Edges[,c(1)], Edges[,c(2)], sep = "-")
write.table(Edges, out_path, sep='\t')

