# Wrapper script that calls the GRN inference tool GENIE3. Please mind the README.md in this directory for
# installation notes.

genie3_r <- getwd()
library("GENIE3")

main <- file.path(genie3_r, '..', '..')
args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]

regulator_path <- paste(main, '/temp/', prefix, '_regulators.csv', sep="")
data_path <- paste(main, '/temp/', prefix, '_expression_data.csv', sep="")
out_path <- paste(main, '/temp/', prefix, '_link_list.csv', sep="")

regulators <- read.table(regulator_path)
regulators <- regulators$V2

exprMat <- read.table(data_path, row.names=1, header = TRUE, sep='\t', as.is=TRUE)
rows = row.names(exprMat)
#regulators <- intersect(regulators, rows)
cols = colnames(exprMat)
exprMat <- as.matrix(exprMat)
rownames(exprMat) <- unlist(rows)
colnames(exprMat) <- unlist(cols)

treeMethod <- 'RF'
K = 'sqrt'
targets <- NULL

weightMat <- GENIE3::GENIE3(exprMat, regulators=regulators, targets=targets, nTrees=1000, K=K, treeMethod=treeMethod)
linkList <- GENIE3::getLinkList(weightMat)
write.table(linkList, out_path, sep='\t')

