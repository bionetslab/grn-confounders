
cemi <- getwd()
library("WGCNA")
library("CEMiTool")
library(igraph)

main <- file.path(cemi, '..', '..')
args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]

data_path <- paste(main, '/temp/', prefix, '_expression_data.csv', sep="")
out_path <- paste(main, '/temp/', prefix, '_edge_list.csv', sep="")

# gene data has to be of the form genes x samples
geneData <- read.table(data_path, row.names=1, header = TRUE, sep='\t', as.is=TRUE)

cem <- new_cem(geneData, filter=TRUE, apply_vst=FALSE)
cem <- get_adj(cem, beta=1)

adj <- adj_data(cem)
#print(adj)
Adjacency <- signumAdjacencyFunction(adj, threshold = 0.6)

Edges<<-graph_from_adjacency_matrix(Adjacency, mode = c("undirected"))
Edges<<- as.data.frame(get.edgelist(Edges))
names(Edges)[1] <- 'node_lower'
names(Edges)[2] <- 'node_upper'
#Edges$edges <<- paste(Edges[,c(1)], Edges[,c(2)], sep = "-")
write.table(Edges, out_path, sep='\t')

