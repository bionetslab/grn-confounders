genie3_r <- dirname(rstudioapi::getSourceEditorContext()$path)
source(paste(genie3_r, "/GENIE3.R", sep=""))

args <- commandArgs(trailingOnly=TRUE)
#data_path = args[0]
data_path <- '/home/anna/BIONETS_code/temp/genie3_expression_data.csv'
#out_path = args[1]
out_path <- '/home/anna/BIONETS_code/temp/genie3_link_list.csv'

#regulators = args[2]
#targets = args[3]
#nTrees = args[4]
treeMethod <- 'RF'
K = 'sqrt'

exprMat <- read.table(data_path, row.names=1, header = TRUE, sep='\t', as.is=TRUE)

#regulators <- rownames(exprMat)
#targets <- rownames(exprMat)

weightMat <- GENIE3(exprMat, nTrees=50, regulators=regulators, returnMatrix = TRUE, targets=targets, K=K, treeMethod=treeMethod)
linkList <- get.link.list(weightMat)

write.table(linkList, out_path, sep='\t')

