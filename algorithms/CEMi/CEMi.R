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

#apply_vst recommended by authors for high-throughput mRNAseq data
cem <- new_cem(geneData, filter=TRUE, apply_vst=TRUE)
beta <- get_beta_data(cem)
beta <- beta$powerEstimate
cem <- get_adj(cem, beta=beta)
#cem <- get_adj(cem, beta=8)
adj <- adj_data(cem)

Edges<<-graph_from_adjacency_matrix(adj, weighted=TRUE, diag=FALSE, mode = "undirected")
Edges <- cbind(get.edgelist(Edges), E(Edges)$weight)
colnames(Edges) <- c('source', 'target', 'score')
write.table(Edges, out_path, sep='\t')

