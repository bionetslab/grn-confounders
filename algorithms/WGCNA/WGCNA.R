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

# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Simulated-05-NetworkConstruction.pdf
adjacency=abs(cor(geneData,use="p"))^6

# Obtaining the edge list from the correlation matrix
adjacency<-as.matrix(adjacency)
adjacency[upper.tri(adjacency)]<-NA
diag(adjacency)<-NA
edges<-melt(adjacency)
colnames(edges)<-c("source","target","score")
edges<-edges[!is.na(edges$score),]
edges<-edges[(edges$score > 0.01),]
edges$score<-as.numeric(as.character(edges$score))

write.table(edges, out_path, sep='\t')

