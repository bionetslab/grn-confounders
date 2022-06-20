genie3_r <- getwd()
source(paste(genie3_r, "/GENIE3.R", sep=""))

main <- file.path(genie3_r, '..', '..')
args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]
data_path <- paste(main, '/temp/', prefix, '_expression_data.csv', sep="")
out_path <- paste(main, '/temp/', prefix, '_link_list.csv', sep="")
regulators = args[2]
targets = args[3]
treeMethod <- 'RF'
K = 'sqrt'

exprMat <- read.table(data_path, row.names=1, header = TRUE, sep='\t', as.is=TRUE)
weightMat <- GENIE3(exprMat, nTrees=50, regulators=regulators, returnMatrix = TRUE, targets=targets, K=K, treeMethod=treeMethod)
linkList <- get.link.list(weightMat)
write.table(linkList, out_path, sep='\t')

