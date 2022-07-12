
cemi <- getwd()
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!require("CEMiTool", quietly = TRUE))
    BiocManager::install("CEMiTool")
if (!require("WGCNA", quietly = TRUE))
    install.packages("WGCNA", repos = "http://cran.us.r-project.org")
if (!require("igraph", quietly = TRUE))
   install.packages("igraph", repos = "http://cran.us.r-project.org")

library(CEMiTool)
library(WGCNA)
library(igraph)

enableWGCNAThreads()

#main <- file.path(wgcna, '..', '..')
args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]

#data_path <- paste(main, '/temp/', prefix, '_expression_data.csv', sep="")
#out_path <- paste(main, '/temp/', prefix, '_edge_list.csv', sep="")

data_path <- paste(prefix, '/temp/', 'cemi_expression_data.csv', sep= "")
out_path <- paste(prefix, '/temp/', 'cemi_edge_list.csv', sep = "")
# gene data has to be of the form genes x samples
geneData <- read.csv(data_path, header = TRUE, sep='\t', row.names = 1, as.is=TRUE)
#geneData <- geneData[-c(1)]
#print(geneData[1:10,1:10])
data(geneData)
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

