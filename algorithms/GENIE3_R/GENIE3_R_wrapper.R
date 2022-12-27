genie3_r <- getwd()
library("GENIE3")

# read arguments
args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]

# data paths
main <- file.path(genie3_r, '..', '..')
regulator_path <- paste(main, '/temp/', prefix, '_regulators.csv', sep="")
data_path <- paste(main, '/temp/', prefix, '_expression_data.csv', sep="")
out_path <- paste(main, '/temp/', prefix, '_link_list.csv', sep="")

# read data
regulators <- read.table(regulator_path)
regulators <- regulators$V2
exprMat <- read.table(data_path, row.names=1, header = TRUE, sep='\t', as.is=TRUE)
rows = row.names(exprMat)
cols = colnames(exprMat)
exprMat <- as.matrix(exprMat)
rownames(exprMat) <- unlist(rows)
colnames(exprMat) <- unlist(cols)

# GENIE3 parameters
treeMethod <- 'RF'
K = 'sqrt'
targets <- NULL

# call GENIE3
weightMat <- GENIE3::GENIE3(exprMat, regulators=regulators, targets=targets, nTrees=1000, K=K, treeMethod=treeMethod)
linkList <- GENIE3::getLinkList(weightMat)
write.table(linkList, out_path, sep='\t')

