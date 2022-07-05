
wgcna <- getwd()
install.packages("WGCNA")
install.packages("igraph")
library(WGCNA)
library(igraph)
enableWGCNAThreads()

main <- file.path(wgcna, '..', '..')
args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]

data_path <- paste(main, '/temp/', prefix, '_expression_data.csv', sep="")
out_path <- paste(main, '/temp/', prefix, '_edge_list.csv', sep="")


geneData <- read.table(data_path, row.names=1, header = TRUE, sep='\t', as.is=TRUE)


Corr <- adjacency(geneData, type = "signed", power = 1)
#Corr[c(1:8),c(1:8)]
Adjacency <- signumAdjacencyFunction(Corr, threshold = 0.6)

Edges<<-graph_from_adjacency_matrix(Adjacency, mode = c("undirected"))
Edges<<- as.data.frame(get.edgelist(Edges))
names(Edges)[1] <- 'node_lower'
names(Edges)[2] <- 'node_upper'
#Edges$edges <<- paste(Edges[,c(1)], Edges[,c(2)], sep = "-")
write.table(Edges, out_path, sep='\t')

