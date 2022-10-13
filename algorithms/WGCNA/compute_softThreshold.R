wgcna <- getwd()
library(WGCNA, warn.conflicts=FALSE)
enableWGCNAThreads(nThreads = 20)
library(reshape2, warn.conflicts=FALSE)
main <- file.path(wgcna, '..', '..')
args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]

data_path <- paste(main, '/temp/', prefix, '_expression_data.csv', sep= "")
out_path <- paste(main, '/temp/', prefix,  '_edge_list.csv', sep = "")

datExpr <- read.csv(data_path, header = TRUE, sep='\t', as.is=TRUE)
datExpr <- datExpr[-c(1)]

# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Simulated-05-NetworkConstruction.pdf
#adjacency=abs(cor(geneData,use="p"))^6
powers = c(c(4:8), seq(from = 10, to=15, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
softPower <- sft$powerEstimate
line = 'softPower for ' + as.character(prefix) + ': ' + as.character(softPower) + ', ' + as.character(Sys.time())
write(line,file="softPowers.txt",append=TRUE)
adjacency <-  adjacency(datExpr, power = softPower)

# Obtaining the edge list from the correlation matrix
adjacency<-as.matrix(adjacency)
adjacency[upper.tri(adjacency)]<-NA
diag(adjacency)<-NA
edges<-melt(adjacency)
colnames(edges)<-c("source","target","score")
edges<-edges[!is.na(edges$score),]
# reduce storage complexity.
# Save to sort out edges with a smaller weight for our purpose, because these won't be in the top 5000 of weighted edges anyways
edges<-edges[(edges$score > 0.001),]
edges$score<-as.numeric(as.character(edges$score))

write.table(edges, out_path, sep='\t')


