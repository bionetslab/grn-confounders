get_netw <- function(prefix) {

library("WGCNA")
library("CEMiTool")
library(igraph)

args <- commandArgs(trailingOnly=TRUE)

data_path <- paste(main, prefix, '_expression_data.csv', sep="")
out_path <- paste(main, prefix,'_edge_list.csv', sep="")

# gene data has to be of the form genes x samples
geneData <- read.table(data_path, row.names=1, header = TRUE, sep='\t', as.is=TRUE)

#apply_vst recommended by authors for high-throughput mRNAseq data
cem <- new_cem(geneData, filter=TRUE, apply_vst=TRUE) # for the metabric BRCA dataset, set apply_vst=FALSE
beta <- get_beta_data(cem)
beta <- beta$powerEstimate
cem <- get_adj(cem, beta=beta)
adj <- adj_data(cem)

Edges<-graph_from_adjacency_matrix(adj, weighted=TRUE, diag=FALSE, mode = "undirected")
Edges <- cbind(get.edgelist(Edges), E(Edges)$weight)
colnames(Edges) <- c('source', 'target', 'score')
write.table(Edges[1:5005,], out_path, sep='\t')
}
